// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_stim::css::{CssError, CssProgram};
use ppvm_stim::{ExecError, parse_extended, sample};
use ppvm_tableau::prelude::*;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

type Tab = GeneralizedTableau<ByteFxHashF64<1>, usize>;

fn compile(src: &str) -> Result<CssProgram, CssError> {
    CssProgram::compile(&parse_extended(src).expect("parse"))
}

/// Shots as rows of outcome codes: 0 / 1, or 2 for a lost qubit.
fn css_shots(src: &str, num_shots: usize, seed: u64) -> Vec<Vec<u8>> {
    let program = compile(src).expect("CSS");
    let m = program.num_measurements();
    let flat = program.sample(num_shots, seed);
    (0..num_shots)
        .map(|s| flat[s * m..][..m].to_vec())
        .collect()
}

fn tableau_shots(src: &str, num_shots: usize, seed: u64) -> Vec<Vec<u8>> {
    let program = parse_extended(src).expect("parse");
    let n = program.num_qubits().max(1);
    let shots = sample(&program, num_shots, |i| {
        Tab::new_with_seed(n, 1e-10, seed.wrapping_add(i as u64))
    })
    .expect("sample");
    shots
        .into_iter()
        .map(|shot| shot.into_iter().map(|r| r.map_or(2, u8::from)).collect())
        .collect()
}

fn assert_always(src: &str, expected: &[u8]) {
    for shot in css_shots(src, 1000, 1) {
        assert_eq!(shot, expected, "{src}");
    }
}

fn fraction(shots: &[Vec<u8>], f: impl Fn(&[u8]) -> bool) -> f64 {
    shots.iter().filter(|s| f(s)).count() as f64 / shots.len() as f64
}

fn assert_close(observed: f64, expected: f64, n: usize) {
    let sigma = (expected * (1.0 - expected) / n as f64).sqrt().max(1e-9);
    assert!(
        (observed - expected).abs() < 5.0 * sigma,
        "observed {observed}, expected {expected}"
    );
}

#[test]
fn deterministic_cliffords() {
    assert_always(
        "X 0\nRX 1\nZ 1\nM 0\nMX 1\nRX 2\nMX 2\nR 3\nY 3\nM 3",
        &[1, 1, 0, 1],
    );
    assert_always("SQRT_Y 0\nSQRT_Y 0\nM 0", &[1]);
    assert_always("SQRT_Y 0\nMX 0\nSQRT_Y_DAG 1\nMX 1", &[0, 1]);
    assert_always("X 0\nCX 0 1\nRX 2\nX 3\nCZ 3 2\nM 0 1\nMX 2", &[1, 1, 1]);
    assert_always("RX 0\nMX 0\nMPAD 0 1\nM(1) 1", &[0, 0, 1, 1]);
    assert_always("REPEAT 3 {\n    X 0\n    M 0\n}", &[1, 0, 1]);
}

#[test]
fn entangled_outcomes_agree() {
    let shots = css_shots("H 0\nCX 0 1\nH 2\nH 3\nCZ 2 3\nH 3\nM 0 1 2 3", 4000, 2);
    assert!(shots.iter().all(|s| s[0] == s[1] && s[2] == s[3]));
    assert_close(fraction(&shots, |s| s[0] == 1), 0.5, shots.len());
    assert_close(fraction(&shots, |s| s[0] == s[2]), 0.5, shots.len());
}

#[test]
fn loss_cancels_gates_until_reload() {
    assert_always("X 0\nI_ERROR[loss](1) 0\nM 0\nM 0\nR 0\nM 0", &[2, 2, 0]);
    assert_always("I_ERROR[loss](1) 0\nMR 0\nM 0", &[2, 0]);
    assert_always("X 0\nI_ERROR[loss](1) 0\nI_ERROR[reset_loss] 0\nM 0", &[0]);
    assert_always("X 0\nI_ERROR[loss](1) 0\nCX 0 1\nM 1", &[0]);
    assert_always("X 0\nI_ERROR[loss](1) 1\nCX 0 1\nM 0 1", &[1, 2]);
    assert_always("RX 1\nX 0\nI_ERROR[loss](1) 0\nCZ 0 1\nMX 1", &[0]);
    assert_always("I_ERROR[loss](1) 0\nM(1) 0", &[2]);
}

#[test]
fn leakage_pins_the_readout() {
    assert_always("X 0\nI_ERROR[leakage](1, 0) 0\nM 0\nX 0\nMX 0", &[0, 0]);
    assert_always("I_ERROR[leakage](0, 1) 0\nM 0\nCX 0 1\nM 1", &[1, 0]);
    assert_always("I_ERROR[leakage](0, 1) 0\nM(1) 0", &[0]);
    assert_always(
        "I_ERROR[leakage](0, 1) 0\nI_ERROR[reset_leakage] 0\nM 0",
        &[0],
    );
    assert_always("I_ERROR[leakage](0, 1) 0\nI_ERROR[reset_loss] 0\nM 0", &[1]);
    assert_always(
        "I_ERROR[leakage](0, 1) 0\nI_ERROR[reset_state] 0\nM 0",
        &[1],
    );
    assert_always("I_ERROR[leakage](0, 1) 0\nR 0\nM 0", &[0]);
    assert_always("I_ERROR[leakage](0, 1) 0\nI_ERROR[loss](1) 0\nM 0", &[2]);
    assert_always("I_ERROR[loss](1) 0\nI_ERROR[leakage](0, 1) 0\nM 0", &[2]);
    assert_always("X 0\nI_ERROR[reset_state] 0\nM 0", &[0]);
}

#[test]
fn feedback_uses_the_recorded_value() {
    let shots = css_shots("H 0\nM(0.1) 0\nCX rec[-1] 1\nM 1", 2000, 3);
    assert!(shots.iter().all(|s| s[0] == s[1]));
    assert_always("X 0\nI_ERROR[loss](1) 0\nM 0\nCX rec[-1] 1\nM 1", &[2, 0]);
    assert_always(
        "X 0\nM 0\nRX 1\nCZ rec[-1] 1\nMX 1\nCY rec[-2] 2\nM 2",
        &[1, 1, 1],
    );
}

#[test]
fn correlated_loss() {
    assert_always("I_ERROR[correlated_loss](1, 0, 0) 0 1\nM 0 1", &[2, 2]);
    assert_always(
        "I_ERROR[loss](1) 0\nI_ERROR[correlated_loss](0, 0, 1) 0 1\nM 0 1",
        &[2, 2],
    );
    let shots = css_shots("I_ERROR[correlated_loss](0.2, 0.3, 1) 0 1\nM 0 1", 20000, 4);
    let n = shots.len();
    assert_close(fraction(&shots, |s| s == [2, 2]), 0.2, n);
    assert_close(fraction(&shots, |s| s == [2, 0]), 0.3, n);
    assert_close(fraction(&shots, |s| s == [0, 2]), 0.3, n);
}

#[test]
fn channel_rates() {
    let n = 20000;
    let shots = css_shots(
        "DEPOLARIZE1(0.3) 0\nPAULI_CHANNEL_1(0.1, 0.2, 0.3) 1\nI_ERROR[loss](0.25) 2\n\
         I_ERROR[leakage](0.1, 0.3) 3\nX_ERROR(0.4) 4\nM 0 1 2 3 4",
        n,
        5,
    );
    assert_close(fraction(&shots, |s| s[0] == 1), 0.2, n);
    assert_close(fraction(&shots, |s| s[1] == 1), 0.3, n);
    assert_close(fraction(&shots, |s| s[2] == 2), 0.25, n);
    assert_close(fraction(&shots, |s| s[3] == 1), 0.3, n);
    assert_close(fraction(&shots, |s| s[4] == 1), 0.4, n);
    let shots = css_shots("DEPOLARIZE2(0.3) 0 1\nM 0 1", n, 6);
    assert_close(fraction(&shots, |s| s == [1, 1]), 0.08, n);
    assert_close(fraction(&shots, |s| s == [1, 0]), 0.08, n);
}

#[test]
fn rejects_non_css_programs() {
    for (src, name) in [
        ("S 0", "S"),
        ("T 0", "T"),
        ("MY 0", "MY"),
        ("RY 0", "RY"),
        ("CY 0 1", "CY"),
        ("MPP X0*X1", "MPP"),
        ("I[R_X(theta=0.1*pi)] 0", "R_X"),
    ] {
        match compile(src) {
            Err(CssError::Exec(ExecError::Unsupported {
                name: found,
                line: 1,
            })) => {
                assert_eq!(found, name, "{src}")
            }
            other => panic!("{src}: {other:?}"),
        }
    }
    assert_eq!(
        compile("CZ 0 1\nCZ 1 2\nCZ 0 2").unwrap_err(),
        CssError::NotCss { line: 3 }
    );
    assert_eq!(
        compile("CX 0 1\nH 1\nCX 0 1").unwrap_err(),
        CssError::NotCss { line: 3 }
    );
    assert!(compile("CX 0 1\nH 1\nR 1\nCX 0 1").is_ok());
    assert!(compile("H 1\nCZ 0 1\nCX 1 2\nH 0\nCX 0 1").is_ok());
    assert!(matches!(
        compile("X_ERROR(1.5) 0"),
        Err(CssError::Exec(ExecError::InvalidProbability { .. }))
    ));
}

#[test]
fn shots_depend_only_on_the_seed() {
    let program = compile("H 0 1 2\nI_ERROR[loss](0.2) 0 1 2\nM 0 1 2").expect("CSS");
    let long = program.sample(2000, 9);
    assert_eq!(program.sample(700, 9), long[..700 * 3]);
    assert_ne!(program.sample(700, 10), long[..700 * 3]);
    assert!(program.sample(0, 9).is_empty());
    assert_eq!(compile("H 0").unwrap().sample(5, 9), Vec::<u8>::new());
}

#[cfg(feature = "rayon")]
#[test]
fn shots_do_not_depend_on_the_thread_count() {
    let program = compile("H 0 1\nCX 0 1\nI_ERROR[loss](0.1) 0 1\nM 0 1").expect("CSS");
    let serial = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| program.sample(5000, 11));
    assert_eq!(program.sample(5000, 11), serial);
}

/// A random circuit that is CSS up to Hadamards, using every supported instruction.
fn random_css_circuit(rng: &mut SmallRng, n: usize, len: usize) -> String {
    let mut hadamard = vec![false; n];
    let mut measurements = 0;
    let mut lines = Vec::new();
    let p = |rng: &mut SmallRng| (rng.random::<f64>() * 0.15 * 1000.0).round() / 1000.0;
    for _ in 0..len {
        let q = rng.random_range(0..n);
        let r = (q + rng.random_range(1..n)) % n;
        let line = match rng.random_range(0..16) {
            0 => {
                hadamard[q] = false;
                format!("{} {q}", ["R", "RX"][rng.random_range(0..2)])
            }
            1 | 2 => {
                hadamard[q] ^= true;
                format!(
                    "{} {q}",
                    ["H", "SQRT_Y", "SQRT_Y_DAG"][rng.random_range(0..3)]
                )
            }
            3 => format!("{} {q}", ["X", "Y", "Z"][rng.random_range(0..3)]),
            4..=6 if hadamard[q] == hadamard[r] => format!("CX {q} {r}"),
            4..=6 => format!("CZ {q} {r}"),
            7 => {
                let name = ["M", "MX", "MR", "MRX"][rng.random_range(0..4)];
                if name.starts_with("MR") {
                    hadamard[q] = false;
                }
                measurements += 1;
                format!("{name}({}) {q}", p(rng))
            }
            8 => match rng.random_range(0..3) {
                0 => format!("DEPOLARIZE1({}) {q}", p(rng)),
                1 => format!("Y_ERROR({}) {q}", p(rng)),
                _ => format!("PAULI_CHANNEL_1({}, {}, {}) {q}", p(rng), p(rng), p(rng)),
            },
            9 => {
                if rng.random() {
                    format!("DEPOLARIZE2({}) {q} {r}", p(rng))
                } else {
                    let args: Vec<_> = (0..15).map(|_| (p(rng) / 5.0).to_string()).collect();
                    format!("PAULI_CHANNEL_2({}) {q} {r}", args.join(", "))
                }
            }
            10 => format!("I_ERROR[loss]({}) {q}", p(rng)),
            11 => format!("I_ERROR[leakage]({}, {}) {q}", p(rng), p(rng)),
            12 => format!(
                "I_ERROR[correlated_loss]({}, {}, {}) {q} {r}",
                p(rng),
                p(rng) / 2.0,
                p(rng)
            ),
            13 => {
                let tag = ["reset_loss", "reset_leakage", "reset_state"][rng.random_range(0..3)];
                if tag == "reset_state" {
                    hadamard[q] = false;
                }
                format!("I_ERROR[{tag}] {q}")
            }
            14 if measurements > 0 => {
                let k = rng.random_range(1..=measurements.min(4));
                format!(
                    "{} rec[-{k}] {q}",
                    ["CX", "CY", "CZ"][rng.random_range(0..3)]
                )
            }
            _ => {
                measurements += 1;
                format!("MPAD({}) {}", p(rng), rng.random_range(0..2))
            }
        };
        lines.push(line);
    }
    let all: Vec<_> = (0..n).map(|q| q.to_string()).collect();
    lines.push(format!("M {}\nMX {}", all.join(" "), all.join(" ")));
    lines.join("\n")
}

/// Indicator features of a shot: outcome 1 and lost, per measurement, plus parities of random
/// subsets of the outcome bits.
fn features(shot: &[u8], subsets: &[Vec<usize>]) -> Vec<bool> {
    let mut f: Vec<bool> = shot.iter().flat_map(|&c| [c == 1, c == 2]).collect();
    f.extend(
        subsets
            .iter()
            .map(|s| s.iter().filter(|&&m| shot[m] == 1).count() % 2 == 1),
    );
    f
}

/// Compare the CSS sampler with the generalized tableau on random circuits: every feature's
/// frequency and every pairwise joint frequency must agree.
#[test]
fn random_circuits_match_the_tableau() {
    let mut rng = SmallRng::seed_from_u64(12);
    let n = 20_000;
    let mut worst: f64 = 0.0;
    for circuit in 0..12 {
        let src = random_css_circuit(&mut rng, 5, 60);
        let ours = css_shots(&src, n, circuit);
        let theirs = tableau_shots(&src, n, 1000 + circuit);
        let m = ours[0].len();
        let subsets: Vec<Vec<usize>> = (0..20)
            .map(|_| (0..m).filter(|_| rng.random::<f64>() < 0.3).collect())
            .collect();
        let fa: Vec<_> = ours.iter().map(|s| features(s, &subsets)).collect();
        let fb: Vec<_> = theirs.iter().map(|s| features(s, &subsets)).collect();
        let k = fa[0].len();
        for i in 0..k {
            for j in i..k {
                let pa = fa.iter().filter(|f| f[i] && f[j]).count() as f64 / n as f64;
                let pb = fb.iter().filter(|f| f[i] && f[j]).count() as f64 / n as f64;
                let sigma = ((pa * (1.0 - pa) + pb * (1.0 - pb)) / n as f64).sqrt();
                if pa != pb {
                    worst = worst.max((pa - pb).abs() / sigma.max(1.0 / n as f64));
                }
            }
        }
        assert!(worst < 5.5, "circuit {circuit}: |z| = {worst}\n{src}");
    }
}
