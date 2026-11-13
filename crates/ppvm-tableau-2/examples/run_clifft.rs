// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Benchmark harness for the clifft-bench workloads: runs a Stim circuit,
//! noise included, on `GeneralizedTableau` one shot at a time with a fresh
//! simulator per shot.
//!
//! ```text
//! cargo run --release -p ppvm-tableau-2 --example run_clifft -- <circuit.stim> [shots] [--records <file>]
//! ```
//!
//! Clifft's plain `R_X(α)`, `R_Y(α)`, `R_Z(α)` and `U3(θ, φ, λ)` / `U(…)`
//! lines (angles in half-turns) are rewritten into `stim-parser`'s `I[…]` tags before parsing.
//! `stim-parser` drops `rec[-k]` targets, so detectors are not evaluated here;
//! `--records` writes each shot's measurement record as a line of `0`/`1` to
//! check detector statistics outside the timed loop. Workloads are fetched by
//! `scripts/fetch_clifft_bench.sh`.

use std::fmt::Write as _;
use std::time::Instant;

use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

use ppvm_tableau_2::{ActiveDormant, DenseState, GeneralizedTableau};
use ppvm_traits_2::{
    Clifford, Depolarizing2, Measure, PauliError, Reset, RotationOne, TGate, U3Gate,
};
use stim_parser::ast::{Axis, ExtendedInstruction, PauliAxis, PauliFactor, Target};
use stim_parser::instructions::{AnnotationKind, GateName, MeasureName, NoiseName};

type Sim = GeneralizedTableau<ActiveDormant<DenseState>>;

#[derive(Clone, Copy)]
enum Basis {
    X,
    Y,
    Z,
}

/// One operation on one qubit or pair; `REPEAT` blocks are unrolled.
#[derive(Clone, Copy)]
enum Op {
    One(GateName, usize),
    Two(GateName, usize, usize),
    T(usize),
    TDag(usize),
    Rotation(Axis, f64, usize),
    U3([f64; 3], usize),
    /// `X`, `Y`, `Z` with these probabilities.
    Pauli([f64; 3], usize),
    Depolarize2(f64, usize, usize),
    /// Recorded outcome flipped with probability `flip`; `reset` returns the
    /// qubit to the basis' `+1` state afterwards.
    Measure {
        basis: Basis,
        q: usize,
        flip: f64,
        reset: bool,
    },
    Reset(Basis, usize),
    /// Classically controlled Pauli: apply it to `q` if `rec[-back]` is 1.
    Feedback {
        pauli: Basis,
        back: usize,
        q: usize,
    },
    /// Pauli product `factors[at..at + len]`.
    Mpp {
        at: usize,
        len: usize,
        flip: f64,
    },
}

/// Rewrite Clifft's plain rotation and `U3` lines into `I[…]` tags. Both write
/// angles in half-turns (`R_Z(α)` rotates by `α·π`), so the numbers carry over.
fn rewrite_dialect(src: &str) -> String {
    let mut out = String::new();
    for line in src.lines() {
        let body = line.trim_start();
        let indent = &line[..line.len() - body.len()];
        let rewritten = ["R_X", "R_Y", "R_Z", "U3", "U"]
            .into_iter()
            .find_map(|name| {
                let rest = body.strip_prefix(name)?.strip_prefix('(')?;
                let (args, targets) = rest.split_once(')')?;
                let turns: Vec<f64> = args
                    .split(',')
                    .map(|a| a.trim().parse::<f64>().expect("angle"))
                    .collect();
                let params = match (name, turns.as_slice()) {
                    ("U3" | "U", [t, p, l]) => {
                        format!("U3(theta={t}*pi, phi={p}*pi, lambda={l}*pi)")
                    }
                    (_, [t]) => format!("{name}(theta={t}*pi)"),
                    _ => panic!("bad arguments in {line:?}"),
                };
                Some(format!("{indent}I[{params}]{targets}"))
            });
        writeln!(out, "{}", rewritten.as_deref().unwrap_or(line)).unwrap();
    }
    out
}

fn flip_arg(args: &[f64]) -> f64 {
    args.first().copied().unwrap_or(0.0)
}

/// Flatten the program into ops, unrolling `REPEAT` blocks. MPP factors go to
/// `factors`.
fn lower(instructions: &[ExtendedInstruction], ops: &mut Vec<Op>, factors: &mut Vec<PauliFactor>) {
    for instruction in instructions {
        match instruction {
            ExtendedInstruction::Gate(g) if g.targets.iter().any(|t| t.as_qubit().is_none()) => ops
                .extend(
                    g.targets
                        .chunks(2)
                        .map(|pair| feedback(g.name, pair[0], pair[1])),
                ),
            ExtendedInstruction::Gate(g) => {
                let qs: Vec<usize> = g
                    .targets
                    .iter()
                    .map(|t| t.as_qubit().expect("only qubit targets are supported"))
                    .collect();
                match g.name {
                    GateName::Reset | GateName::ResetZ => {
                        ops.extend(qs.iter().map(|&q| Op::Reset(Basis::Z, q)))
                    }
                    GateName::ResetX => ops.extend(qs.iter().map(|&q| Op::Reset(Basis::X, q))),
                    GateName::ResetY => ops.extend(qs.iter().map(|&q| Op::Reset(Basis::Y, q))),
                    GateName::CX
                    | GateName::ZCX
                    | GateName::CNot
                    | GateName::CY
                    | GateName::ZCY
                    | GateName::CZ
                    | GateName::ZCZ => {
                        ops.extend(qs.chunks(2).map(|p| Op::Two(g.name, p[0], p[1])))
                    }
                    GateName::Identity => {}
                    name => ops.extend(qs.iter().map(|&q| Op::One(name, q))),
                }
            }
            ExtendedInstruction::T { targets, .. } => ops.extend(targets.iter().map(|&q| Op::T(q))),
            ExtendedInstruction::TDag { targets, .. } => {
                ops.extend(targets.iter().map(|&q| Op::TDag(q)))
            }
            ExtendedInstruction::Rotation {
                axis,
                theta,
                targets,
                ..
            } => ops.extend(targets.iter().map(|&q| Op::Rotation(*axis, *theta, q))),
            ExtendedInstruction::U3 {
                theta,
                phi,
                lambda,
                targets,
                ..
            } => ops.extend(targets.iter().map(|&q| Op::U3([*theta, *phi, *lambda], q))),
            ExtendedInstruction::Noise(n) => {
                let p = flip_arg(&n.args);
                let one = |probabilities: [f64; 3]| {
                    n.targets.iter().map(move |&q| Op::Pauli(probabilities, q))
                };
                match n.name {
                    NoiseName::XError => ops.extend(one([p, 0.0, 0.0])),
                    NoiseName::YError => ops.extend(one([0.0, p, 0.0])),
                    NoiseName::ZError => ops.extend(one([0.0, 0.0, p])),
                    NoiseName::Depolarize1 => ops.extend(one([p / 3.0; 3])),
                    NoiseName::PauliChannel1 => ops.extend(one([n.args[0], n.args[1], n.args[2]])),
                    NoiseName::Depolarize2 => ops.extend(
                        n.targets
                            .chunks(2)
                            .map(|pair| Op::Depolarize2(p, pair[0], pair[1])),
                    ),
                    other => panic!("unsupported noise {other:?}"),
                }
            }
            ExtendedInstruction::Measure(m) => {
                let (basis, reset) = match m.name {
                    MeasureName::M | MeasureName::MZ => (Basis::Z, false),
                    MeasureName::MX => (Basis::X, false),
                    MeasureName::MY => (Basis::Y, false),
                    MeasureName::MR => (Basis::Z, true),
                    MeasureName::MRX => (Basis::X, true),
                    MeasureName::MRY => (Basis::Y, true),
                    other => panic!("unsupported measurement {other:?}"),
                };
                let flip = flip_arg(&m.args);
                ops.extend(m.targets.iter().map(|&q| Op::Measure {
                    basis,
                    q,
                    flip,
                    reset,
                }));
            }
            ExtendedInstruction::Mpp(mpp) => {
                let flip = flip_arg(&mpp.args);
                for product in &mpp.products {
                    ops.push(Op::Mpp {
                        at: factors.len(),
                        len: product.len(),
                        flip,
                    });
                    factors.extend_from_slice(product);
                }
            }
            ExtendedInstruction::Repeat { count, body, .. } => {
                let mut block = Vec::new();
                lower(body, &mut block, factors);
                for _ in 0..*count {
                    ops.extend_from_slice(&block);
                }
            }
            ExtendedInstruction::Annotation(a) => match a.kind {
                AnnotationKind::Detector
                | AnnotationKind::ObservableInclude
                | AnnotationKind::QubitCoords
                | AnnotationKind::ShiftCoords
                | AnnotationKind::Tick => {}
            },
            other => panic!("unsupported instruction {other:?}"),
        }
    }
}

/// `CX`/`CY`/`CZ` with a `rec[-k]` control (either side for `CZ`).
fn feedback(name: GateName, a: Target, b: Target) -> Op {
    let pauli = match name {
        GateName::CX | GateName::ZCX | GateName::CNot => Basis::X,
        GateName::CY | GateName::ZCY => Basis::Y,
        GateName::CZ | GateName::ZCZ => Basis::Z,
        other => panic!("unsupported classically controlled gate {other:?}"),
    };
    match (a, b, pauli) {
        (Target::Rec(back), Target::Qubit(q), _)
        | (Target::Qubit(q), Target::Rec(back), Basis::Z) => Op::Feedback { pauli, back, q },
        _ => panic!("unsupported targets {a:?} {b:?} for {name:?}"),
    }
}

fn apply_one(s: &mut Sim, name: GateName, q: usize) {
    match name {
        GateName::X => s.x(q),
        GateName::Y => s.y(q),
        GateName::Z => s.z(q),
        GateName::H | GateName::HXZ => s.h(q),
        GateName::S | GateName::SqrtZ => s.s(q),
        GateName::SDag | GateName::SqrtZDag => s.s_dag(q),
        GateName::SqrtX => s.sqrt_x(q),
        GateName::SqrtXDag => s.sqrt_x_dag(q),
        GateName::SqrtY => s.sqrt_y(q),
        GateName::SqrtYDag => s.sqrt_y_dag(q),
        GateName::T => s.t(q),
        GateName::TDag => s.t_dag(q),
        other => panic!("unsupported gate {other:?}"),
    }
}

fn apply_two(s: &mut Sim, name: GateName, a: usize, b: usize) {
    match name {
        GateName::CX | GateName::ZCX | GateName::CNot => s.cnot(a, b),
        GateName::CY | GateName::ZCY => s.cy(a, b),
        GateName::CZ | GateName::ZCZ => s.cz(a, b),
        other => panic!("unsupported gate {other:?}"),
    }
}

/// Map `basis` on `q` to `Z` (`to_z`) or back.
fn change_basis(s: &mut Sim, basis: Basis, q: usize, to_z: bool) {
    match (basis, to_z) {
        (Basis::Z, _) => {}
        (Basis::X, _) => s.h(q),
        (Basis::Y, true) => {
            s.s_dag(q);
            s.h(q);
        }
        (Basis::Y, false) => {
            s.h(q);
            s.s(q);
        }
    }
}

fn pauli_basis(axis: PauliAxis) -> Basis {
    match axis {
        PauliAxis::X => Basis::X,
        PauliAxis::Y => Basis::Y,
        PauliAxis::Z => Basis::Z,
    }
}

/// Run `op`, appending measurement outcomes to `record`; `rng` flips recorded
/// outcomes.
fn run(s: &mut Sim, op: Op, factors: &[PauliFactor], rng: &mut SmallRng, record: &mut Vec<bool>) {
    let mut noisy = |m: bool, flip: f64| m ^ (flip > 0.0 && rng.random::<f64>() < flip);
    match op {
        Op::One(name, q) => apply_one(s, name, q),
        Op::Two(name, a, b) => apply_two(s, name, a, b),
        Op::T(q) => s.t(q),
        Op::TDag(q) => s.t_dag(q),
        Op::Rotation(axis, theta, q) => match axis {
            Axis::X => s.rx(q, theta),
            Axis::Y => s.ry(q, theta),
            Axis::Z => s.rz(q, theta),
        },
        Op::U3([theta, phi, lambda], q) => s.u3(q, theta, phi, lambda),
        Op::Pauli(p, q) => s.pauli_error(q, p),
        Op::Depolarize2(p, a, b) => s.depolarize2(a, b, p),
        Op::Measure {
            basis,
            q,
            flip,
            reset,
        } => {
            change_basis(s, basis, q, true);
            let m = s.measure(q).expect("no loss");
            if reset && m {
                s.x(q);
            }
            change_basis(s, basis, q, false);
            record.push(noisy(m, flip));
        }
        Op::Feedback { pauli, back, q } => {
            if record[record.len() - back] {
                match pauli {
                    Basis::X => s.x(q),
                    Basis::Y => s.y(q),
                    Basis::Z => s.z(q),
                }
            }
        }
        Op::Reset(basis, q) => match basis {
            Basis::Z => s.reset(q),
            Basis::X => s.reset_x(q),
            Basis::Y => s.reset_y(q),
        },
        // Basis changes to Z, a CX ladder onto the first qubit, measure, undo.
        Op::Mpp { at, len, flip } => {
            let product = &factors[at..at + len];
            let first = product[0].qubit;
            for f in product {
                change_basis(s, pauli_basis(f.axis), f.qubit, true);
            }
            for f in &product[1..] {
                s.cnot(f.qubit, first);
            }
            let m = s.measure(first).expect("no loss");
            for f in product[1..].iter().rev() {
                s.cnot(f.qubit, first);
            }
            for f in product {
                change_basis(s, pauli_basis(f.axis), f.qubit, false);
            }
            record.push(noisy(m, flip));
        }
    }
}

fn qubits(op: &Op, factors: &[PauliFactor]) -> Vec<usize> {
    match *op {
        Op::Two(_, a, b) | Op::Depolarize2(_, a, b) => vec![a, b],
        Op::One(_, q)
        | Op::T(q)
        | Op::TDag(q)
        | Op::Rotation(_, _, q)
        | Op::U3(_, q)
        | Op::Pauli(_, q)
        | Op::Measure { q, .. }
        | Op::Reset(_, q)
        | Op::Feedback { q, .. } => vec![q],
        Op::Mpp { at, len, .. } => factors[at..at + len].iter().map(|f| f.qubit).collect(),
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: run_clifft <circuit.stim> [shots] [--records <file>]");
    let (mut shots, mut records_path) = (1000usize, None);
    while let Some(arg) = args.next() {
        if arg == "--records" {
            records_path = Some(args.next().expect("--records needs a path"));
        } else {
            shots = arg.parse().expect("shots");
        }
    }
    let src = rewrite_dialect(&std::fs::read_to_string(&path).expect("read circuit"));
    let program =
        stim_parser::parse_extended(&src).unwrap_or_else(|e| panic!("parse error: {e:?}"));
    let (mut ops, mut factors) = (Vec::new(), Vec::new());
    lower(&program.instructions, &mut ops, &mut factors);
    let n = 1 + ops
        .iter()
        .flat_map(|op| qubits(op, &factors))
        .max()
        .expect("empty circuit");

    let mut rng = SmallRng::seed_from_u64(0x9e37_79b9_7f4a_7c15);
    let (mut record, mut records) = (Vec::new(), String::new());
    let (mut parity, mut max_k, mut k_at_measurements, mut measurements) =
        (false, 0, 0usize, 0usize);
    let start = Instant::now();
    for _ in 0..shots {
        let mut s = Sim::new_with_seed(n, rng.random());
        record.clear();
        for &op in &ops {
            run(&mut s, op, &factors, &mut rng, &mut record);
            let k = s.state().active().len();
            max_k = max_k.max(k);
            if matches!(op, Op::Measure { .. } | Op::Mpp { .. }) {
                k_at_measurements += k;
                measurements += 1;
            }
        }
        parity ^= record.iter().fold(false, |p, &m| p ^ m);
        if records_path.is_some() {
            records.extend(record.iter().map(|&m| if m { '1' } else { '0' }));
            records.push('\n');
        }
    }
    let elapsed = start.elapsed().as_secs_f64();
    if let Some(p) = records_path {
        std::fs::write(p, records).expect("write records");
    }
    println!(
        "{path}: {n} qubits, {} ops ({} measurements) per shot, {shots} shots",
        ops.len(),
        record.len()
    );
    println!(
        "{:.1} µs/shot, {:.0} shots/s, max k {max_k}, mean k at measurements {:.2} (parity {parity})",
        elapsed / shots as f64 * 1e6,
        shots as f64 / elapsed,
        k_at_measurements as f64 / measurements.max(1) as f64,
    );
}
