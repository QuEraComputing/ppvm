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
    // Pairs act in order, a leaked qubit counts as present.
    assert_always(
        "I_ERROR[loss](1) 0\nI_ERROR[correlated_loss](0, 0, 1) 0 1 1 2\nM 0 1 2",
        &[2, 2, 2],
    );
    assert_always(
        "I_ERROR[leakage](0, 1) 0\nI_ERROR[correlated_loss](1, 0, 0) 0 1\nM 0 1",
        &[2, 2],
    );
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
    assert!(compile("CX 0 1\nH 1\nM 1\nCX 0 1").is_ok());
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
    let mut out = vec![0; 700 * 3];
    program.sample_into(&mut out, 9);
    assert_eq!(out, long[..700 * 3]);
    assert_eq!(compile("H 0").unwrap().sample(5, 9), Vec::<u8>::new());
}

#[test]
fn shots_are_independent() {
    // 40 random bits per shot: a repeated shot means two shots share their randomness.
    let qubits: Vec<_> = (0..40).map(|q| q.to_string()).collect();
    let qubits = qubits.join(" ");
    let shots = css_shots(&format!("H {qubits}\nM {qubits}"), 3000, 14);
    let distinct: std::collections::HashSet<_> = shots.iter().collect();
    assert_eq!(distinct.len(), shots.len());
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

#[test]
fn measurements_start_new_segments() {
    for shot in css_shots("H 0\nCX 0 1\nM 0\nCZ 0 1\nM 0 1", 1000, 13) {
        assert!(shot[0] == shot[1] && shot[1] == shot[2]);
    }
}

#[test]
fn repeat_is_not_unrolled() {
    let program = compile("REPEAT 1000000000 {\n    H 0\n    CX 0 1\n    R 1\n}\nM 0 1").unwrap();
    assert_eq!(program.num_measurements(), 2);
    assert_eq!(
        compile("REPEAT 1000000000 {\n    CX 0 1\n    H 1\n}").unwrap_err(),
        CssError::NotCss { line: 2 }
    );
    // Three Hadamards act as one, four as none.
    assert_eq!(
        compile("CX 0 1\nREPEAT 3 {\n    H 0\n}\nCX 0 1").unwrap_err(),
        CssError::NotCss { line: 5 }
    );
    assert!(compile("CX 0 1\nREPEAT 4 {\n    H 0\n}\nCX 0 1").is_ok());
    assert!(compile("CX 0 1\nREPEAT 1000000001 {\n    H 0\n}\nCZ 0 1").is_ok());
    assert!(compile("CX 0 1\nREPEAT 1000000001 {\n    H 0\n}\nCX 0 1").is_err());
}

/// Random circuits using every supported instruction. With `css`, every qubit segment gets a
/// hidden random Hadamard label and two-qubit gates are CX between equal labels and CZ otherwise,
/// so the circuit is CSS up to Hadamards; without, CX and CZ are chosen at random.
struct RandomCircuit {
    rng: SmallRng,
    css: bool,
    label: Vec<bool>,
    measurements: usize,
}

impl RandomCircuit {
    fn new(seed: u64, n: usize, css: bool) -> Self {
        let mut rng = SmallRng::seed_from_u64(seed);
        let label = (0..n).map(|_| rng.random()).collect();
        Self {
            rng,
            css,
            label,
            measurements: 0,
        }
    }

    fn lines(&mut self, len: usize, nested: bool) -> Vec<String> {
        let mut lines = Vec::new();
        for _ in 0..len {
            if !nested && self.rng.random::<f64>() < 0.05 {
                lines.extend(self.repeat());
            } else {
                lines.push(self.line());
            }
        }
        lines
    }

    /// A `REPEAT` whose body restores the Hadamard labels it started with.
    fn repeat(&mut self) -> Vec<String> {
        let (start, before) = (self.label.clone(), self.measurements);
        let count = self.rng.random_range(1..4);
        let len = self.rng.random_range(2..8);
        let mut lines = vec![format!("REPEAT {count} {{")];
        lines.extend(self.lines(len, true));
        for (q, &label) in start.iter().enumerate() {
            if self.css && self.label[q] != label {
                lines.push(format!("H {q}"));
                self.label[q] = label;
            }
        }
        lines.push("}".to_string());
        self.measurements = before + count * (self.measurements - before);
        lines
    }

    fn line(&mut self) -> String {
        let n = self.label.len();
        let q = self.rng.random_range(0..n);
        let r = (q + self.rng.random_range(1..n)) % n;
        let s = (r + self.rng.random_range(1..n - 1)) % n;
        let s = if s == q { (s + 1) % n } else { s };
        let pick =
            |rng: &mut SmallRng, names: &[&'static str]| names[rng.random_range(0..names.len())];
        match self.rng.random_range(0..18) {
            0 => {
                self.new_segment(q);
                format!("{} {q}", pick(&mut self.rng, &["R", "RX"]))
            }
            1 | 2 => {
                self.label[q] ^= true;
                format!(
                    "{} {q}",
                    pick(&mut self.rng, &["H", "SQRT_Y", "SQRT_Y_DAG"])
                )
            }
            3 => format!("{} {q}", pick(&mut self.rng, &["X", "Y", "Z"])),
            4..=6 => self.two_qubit(q, r),
            7 => {
                let name = pick(&mut self.rng, &["M", "MX", "MR", "MRX"]);
                self.new_segment(q);
                self.measurements += 1;
                format!("{name}({}) {q}", self.p())
            }
            8 => match self.rng.random_range(0..3) {
                0 => format!("DEPOLARIZE1({}) {q}", self.p()),
                1 => format!("Y_ERROR({}) {q}", self.p()),
                _ => format!(
                    "PAULI_CHANNEL_1({}, {}, {}) {q}",
                    self.small(),
                    self.small(),
                    self.small()
                ),
            },
            9 => {
                if self.rng.random() {
                    format!("DEPOLARIZE2({}) {q} {r}", self.small())
                } else {
                    let args: Vec<_> = (0..15).map(|_| (self.small() / 5.0).to_string()).collect();
                    format!("PAULI_CHANNEL_2({}) {q} {r}", args.join(", "))
                }
            }
            10 => format!("I_ERROR[loss]({}) {q}", self.p()),
            11 => {
                let (p0, p1) = match self.rng.random_range(0..10) {
                    0 => (1.0, 0.0),
                    1 => (0.0, 1.0),
                    _ => (self.small(), self.small()),
                };
                format!("I_ERROR[leakage]({p0}, {p1}) {q}")
            }
            12 => {
                let targets = if self.rng.random() {
                    format!("{q} {r}")
                } else {
                    format!("{q} {r} {r} {s}")
                };
                let (p0, p1, p2) = (self.small(), self.small() / 2.0, self.p());
                format!("I_ERROR[correlated_loss]({p0}, {p1}, {p2}) {targets}")
            }
            13 => {
                let tag = pick(
                    &mut self.rng,
                    &["reset_loss", "reset_leakage", "reset_state"],
                );
                if tag == "reset_state" {
                    self.new_segment(q);
                }
                format!("I_ERROR[{tag}] {q}")
            }
            14 if self.measurements > 0 => {
                let k = self.rng.random_range(1..=self.measurements.min(4));
                format!("{} rec[-{k}] {q}", pick(&mut self.rng, &["CX", "CY", "CZ"]))
            }
            15 => {
                let name = pick(
                    &mut self.rng,
                    &["H", "M", "R", "X_ERROR(0.1)", "I_ERROR[loss](0.1)"],
                );
                for t in [q, r] {
                    match name {
                        "H" => self.label[t] ^= true,
                        "M" | "R" => self.new_segment(t),
                        _ => {}
                    }
                }
                if name == "M" {
                    self.measurements += 2;
                }
                format!("{name} {q} {r}")
            }
            16 => format!("{}\n{}", self.two_qubit(q, r), self.two_qubit(r, s)),
            _ => {
                self.measurements += 1;
                format!("MPAD({}) {}", self.p(), self.rng.random_range(0..2))
            }
        }
    }

    fn two_qubit(&mut self, a: usize, b: usize) -> String {
        let cx = if self.css {
            self.label[a] == self.label[b]
        } else {
            self.rng.random()
        };
        let name = if cx { "CX" } else { "CZ" };
        format!("{name} {a} {b}")
    }

    fn new_segment(&mut self, q: usize) {
        self.label[q] = self.rng.random();
    }

    /// A probability, sometimes exactly 0 or 1.
    fn p(&mut self) -> f64 {
        match self.rng.random_range(0..20) {
            0 => 0.0,
            1 => 1.0,
            _ => self.small(),
        }
    }

    fn small(&mut self) -> f64 {
        (self.rng.random::<f64>() * 150.0).round() / 1000.0
    }
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

/// Largest z-score between two sets of event frequencies.
fn max_z(pairs: impl Iterator<Item = (f64, f64)>, n: usize) -> f64 {
    let n = n as f64;
    pairs
        .filter(|(pa, pb)| pa != pb)
        .map(|(pa, pb)| {
            let sigma = ((pa * (1.0 - pa) + pb * (1.0 - pb)) / n).sqrt();
            (pa - pb).abs() / sigma.max(1.0 / n)
        })
        .fold(0.0, f64::max)
}

/// Compare the CSS sampler with the generalized tableau on random circuits: every feature's
/// frequency and every pairwise joint frequency must agree.
#[test]
fn random_circuits_match_the_tableau() {
    let n = 20_000;
    for circuit in 0..16 {
        let mut circuit_gen = RandomCircuit::new(circuit, 5, true);
        let mut lines = circuit_gen.lines(60, false);
        lines.push("M 0 1 2 3 4\nMX 0 1 2 3 4".to_string());
        let src = lines.join("\n");
        let ours = css_shots(&src, n, circuit);
        let theirs = tableau_shots(&src, n, 1000 + circuit);
        let m = ours[0].len();
        let subsets: Vec<Vec<usize>> = (0..20)
            .map(|_| {
                (0..m)
                    .filter(|_| circuit_gen.rng.random::<f64>() < 0.3)
                    .collect()
            })
            .collect();
        let fa: Vec<_> = ours.iter().map(|s| features(s, &subsets)).collect();
        let fb: Vec<_> = theirs.iter().map(|s| features(s, &subsets)).collect();
        let k = fa[0].len();
        let frequency = |f: &[Vec<bool>], i: usize, j: usize| {
            f.iter().filter(|f| f[i] && f[j]).count() as f64 / n as f64
        };
        let pairs = (0..k).flat_map(|i| (i..k).map(move |j| (i, j)));
        let z = max_z(
            pairs.map(|(i, j)| (frequency(&fa, i, j), frequency(&fb, i, j))),
            n,
        );
        assert!(z < 5.5, "circuit {circuit}: |z| = {z}\n{src}");
    }
}

/// Compare full outcome histograms on random circuits with at most seven measurements.
#[test]
fn random_histograms_match_the_tableau() {
    let n = 20_000;
    let mut tested = 0;
    for circuit in 100.. {
        let mut circuit_gen = RandomCircuit::new(circuit, 3, true);
        let mut lines = circuit_gen.lines(15, false);
        lines.push("M 0 1 2".to_string());
        let src = lines.join("\n");
        if compile(&src).unwrap().num_measurements() > 7 {
            continue;
        }
        let histogram = |shots: Vec<Vec<u8>>| {
            let mut counts = std::collections::HashMap::new();
            for shot in shots {
                *counts.entry(shot).or_insert(0usize) += 1;
            }
            counts
        };
        let ours = histogram(css_shots(&src, n, circuit));
        let theirs = histogram(tableau_shots(&src, n, 1000 + circuit));
        let cells = ours
            .keys()
            .chain(theirs.keys())
            .collect::<std::collections::HashSet<_>>();
        let frequency = |h: &std::collections::HashMap<Vec<u8>, usize>, c: &Vec<u8>| {
            h.get(c).copied().unwrap_or(0) as f64 / n as f64
        };
        let z = max_z(
            cells
                .iter()
                .map(|c| (frequency(&ours, c), frequency(&theirs, c))),
            n,
        );
        assert!(z < 5.5, "circuit {circuit}: |z| = {z}\n{src}");
        tested += 1;
        if tested == 25 {
            break;
        }
    }
}

/// Circuits that ignore Hadamard labels: whatever the CSS check accepts must match the tableau.
#[test]
fn accepted_circuits_match_the_tableau() {
    let n = 20_000;
    let mut accepted = 0;
    for circuit in 200..300 {
        let mut circuit_gen = RandomCircuit::new(circuit, 4, false);
        let mut lines = circuit_gen.lines(25, false);
        lines.push("M 0 1 2 3\nMX 0 1 2 3".to_string());
        let src = lines.join("\n");
        if compile(&src).is_err() {
            continue;
        }
        accepted += 1;
        let ours = css_shots(&src, n, circuit);
        let theirs = tableau_shots(&src, n, 1000 + circuit);
        let k = 2 * ours[0].len();
        let fa: Vec<_> = ours.iter().map(|s| features(s, &[])).collect();
        let fb: Vec<_> = theirs.iter().map(|s| features(s, &[])).collect();
        let frequency = |f: &[Vec<bool>], i: usize, j: usize| {
            f.iter().filter(|f| f[i] && f[j]).count() as f64 / n as f64
        };
        let pairs = (0..k).flat_map(|i| (i..k).map(move |j| (i, j)));
        let z = max_z(
            pairs.map(|(i, j)| (frequency(&fa, i, j), frequency(&fb, i, j))),
            n,
        );
        assert!(z < 5.5, "circuit {circuit}: |z| = {z}\n{src}");
    }
    assert!(accepted >= 20, "only {accepted} circuits accepted");
}
