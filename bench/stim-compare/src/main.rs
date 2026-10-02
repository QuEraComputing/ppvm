//! Time one surface_d30 shot on ppvm's `GeneralizedTableau` and on Stim's
//! `TableauSimulator`, interleaved so machine drift hits both sides equally.
//!
//! Usage: `stim-compare [shots] [only] [variant]` where `only` is `ppvm`, `stim`
//! or `both`, and `variant` is one of [`VARIANTS`] (default: all of them).
//! `stim-compare file <path> [shots] [only]` times one `.stim` file instead.

use std::time::{Duration, Instant};

use bnum::types::U2048;
#[cfg(feature = "legacy")]
use ppvm_stim::backend::config::indexmap::ByteFxHashF64;
use ppvm_stim::backend::prelude::*;
use ppvm_stim::{ExtendedProgram, execute_validated_with_rng, parse_extended, validate};
use rand::SeedableRng;

#[cfg(feature = "legacy")]
type Tab = GeneralizedTableau<ByteFxHashF64<237>, U2048>;
#[cfg(not(feature = "legacy"))]
type Tab = GeneralizedTableau<U2048>;

const SRC: &str = include_str!("../../../crates/ppvm-stim/examples/surface_d30.stim");
const WARMUP: usize = 2;

const MEASURE_OPS: &[&str] = &["M", "MR", "R", "DETECTOR", "OBSERVABLE_INCLUDE"];
const NOISE_OPS: &[&str] = &["DEPOLARIZE1", "DEPOLARIZE2", "X_ERROR"];

/// (name, ops stripped from surface_d30, `REPEAT 29` count). Records go with
/// the measurements because `rec[-k]` targets would dangle.
const VARIANTS: &[(&str, &[&str], &[&str], u32)] = &[
    ("full", &[], &[], 29),
    ("full-r1", &[], &[], 1),
    ("no-noise", NOISE_OPS, &[], 29),
    ("no-measure", MEASURE_OPS, &[], 29),
    ("no-measure-r1", MEASURE_OPS, &[], 1),
    ("gates-only", MEASURE_OPS, NOISE_OPS, 29),
];

/// Drop every line whose instruction name is in `a` or `b`; set the round count.
fn strip(src: &str, a: &[&str], b: &[&str], reps: u32) -> String {
    src.lines()
        .filter(|line| {
            let name = line.trim_start().split(['(', ' ']).next().unwrap_or("");
            !a.contains(&name) && !b.contains(&name)
        })
        .map(|line| match line {
            "REPEAT 29 {" => format!("REPEAT {reps} {{\n"),
            _ => format!("{line}\n"),
        })
        .collect()
}

/// Per shot: (construct time, execute time, number of `1` outcomes).
type Shot = (Duration, Duration, usize);

fn ppvm_shot(prog: &ExtendedProgram, n_qubits: usize, seed: u64) -> Shot {
    let t0 = Instant::now();
    let mut tab = Tab::new(n_qubits, 1e-10);
    let mut rng = rand::rngs::SmallRng::seed_from_u64(seed);
    let mut rec = Vec::with_capacity(prog.measurement_count());
    let t1 = Instant::now();
    execute_validated_with_rng(&prog.instructions, &mut tab, &mut rec, &mut rng);
    let t2 = Instant::now();
    (t1 - t0, t2 - t1, rec.iter().filter(|&&b| b == Some(true)).count())
}

fn stim_shot(circuit: &stim::Circuit, n_qubits: usize, seed: u64) -> Shot {
    let t0 = Instant::now();
    let mut sim = stim::TableauSimulator::with_seed(seed);
    // Pre-size so do_circuit's ensure_large_enough_for_qubits is a no-op.
    sim.set_num_qubits(n_qubits);
    let t1 = Instant::now();
    sim.do_circuit(circuit);
    let t2 = Instant::now();
    let ones = sim.current_measurement_record().iter().filter(|&&b| b).count();
    (t1 - t0, t2 - t1, ones)
}

fn median(times: &mut [Duration]) -> Duration {
    times.sort();
    times[times.len() / 2]
}

fn stats(name: &str, shots: &[Shot], n_meas: usize) -> Duration {
    let mut build: Vec<_> = shots.iter().map(|s| s.0).collect();
    let mut exec: Vec<_> = shots.iter().map(|s| s.1).collect();
    let mean_ones = shots.iter().map(|s| s.2).sum::<usize>() as f64 / shots.len() as f64;
    let (b, e) = (median(&mut build), median(&mut exec));
    println!(
        "{name:>5}: construct {b:>9.2?}  execute median {e:>9.2?} (min {:>9.2?}, max {:>9.2?})  | mean ones/shot {mean_ones:.3} of {n_meas}",
        exec[0],
        exec[exec.len() - 1],
    );
    e
}

fn run_variant(name: &str, src: &str, shots: usize, run_ppvm: bool, run_stim: bool) {
    let prog = parse_extended(src).expect("ppvm parse");
    let circuit: stim::Circuit = src.parse().expect("stim parse");
    let n_qubits = circuit.num_qubits().max(1);
    validate(&prog).expect("ppvm validate");
    let n_meas = prog.measurement_count();
    assert_eq!(n_meas as u64, circuit.num_measurements());
    println!("\n[{name}] {n_qubits} qubits, {n_meas} measurements");

    let (mut ps, mut ss) = (vec![], vec![]);
    for i in 0..WARMUP + shots {
        let seed = i as u64;
        let p = run_ppvm.then(|| ppvm_shot(&prog, n_qubits, seed));
        let s = run_stim.then(|| stim_shot(&circuit, n_qubits, seed));
        if i >= WARMUP {
            ps.extend(p);
            ss.extend(s);
        }
    }

    let p = run_ppvm.then(|| stats("ppvm", &ps, n_meas));
    let s = run_stim.then(|| stats("stim", &ss, n_meas));
    if let (Some(p), Some(s)) = (p, s) {
        println!("execute median ppvm / stim = {:.2}x", p.as_secs_f64() / s.as_secs_f64());
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let backend = if cfg!(feature = "legacy") { "legacy" } else { "traits-2" };
    if args.get(1).map(String::as_str) == Some("file") {
        let path = args.get(2).expect("file path");
        let shots: usize = args.get(3).map_or(20, |s| s.parse().expect("shots"));
        let only = args.get(4).map_or("ppvm", String::as_str);
        let src = std::fs::read_to_string(path).expect("read circuit");
        println!("{path}: ppvm backend = {backend}, {shots} shots");
        run_variant("file", &src, shots, only != "stim", only != "ppvm");
        return;
    }
    let shots: usize = args.get(1).map_or(20, |s| s.parse().expect("shots"));
    let only = args.get(2).map_or("both", String::as_str);
    let variant = args.get(3).map(String::as_str);
    println!("surface_d30: ppvm backend = {backend}, {shots} shots");

    for &(name, a, b, reps) in VARIANTS {
        if variant.is_none_or(|v| v == name) {
            let src = strip(SRC, a, b, reps);
            run_variant(name, &src, shots, only != "stim", only != "ppvm");
        }
    }
}
