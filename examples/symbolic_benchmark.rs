// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0
//! Exact symbolic TFIM propagation; see benchmarks/symbolic/README.md.
use ppvm_pauli_sum::prelude::*;
use ppvm_sym::Term;
use std::{hint::black_box, time::Instant};
type Symbolic = PauliSum<config::fxhash::Byte<1, Term>>;
type Numeric = PauliSum<config::fxhash::Byte<1, f64>>;
#[derive(Clone, Copy)]
struct Gate {
    q: usize,
    bond: bool,
    parameter: usize,
}
fn circuit(n: usize, depth: usize, shared: bool) -> Vec<Gate> {
    let mut gates = Vec::new();
    for _ in 0..depth {
        for q in 0..n {
            gates.push(Gate {
                q,
                bond: false,
                parameter: if shared { 0 } else { gates.len() },
            });
        }
        for q in 0..n - 1 {
            gates.push(Gate {
                q,
                bond: true,
                parameter: if shared { 1 } else { gates.len() },
            });
        }
    }
    gates
}
fn angles(count: usize, sample: usize) -> Vec<f64> {
    (0..count)
        .map(|i| 0.07 + ((i * 17 + sample * 11) % 53) as f64 / 100.0)
        .collect()
}
fn symbolic(n: usize, gates: &[Gate]) -> Symbolic {
    let mut state = Symbolic::builder().n_qubits(n).build();
    state += (format!("Z{}", "I".repeat(n - 1)).as_str(), Term::from(1.0));
    for g in gates.iter().rev() {
        let mut theta = Term::var(g.parameter as u32);
        theta.set_min_eps(0.0);
        if g.bond {
            state.rzz(g.q, g.q + 1, theta);
        } else {
            state.rx(g.q, theta);
        }
    }
    state
}
fn numeric(n: usize, gates: &[Gate], values: &[f64]) -> Numeric {
    let mut state = Numeric::builder().n_qubits(n).build();
    state += (format!("Z{}", "I".repeat(n - 1)).as_str(), 1.0);
    for g in gates.iter().rev() {
        if g.bond {
            state.rzz(g.q, g.q + 1, values[g.parameter]);
        } else {
            state.rx(g.q, values[g.parameter]);
        }
    }
    state
}
fn word(w: &impl PauliWordTrait, n: usize) -> String {
    (0..n)
        .map(|q| match (w.get_xbit(q), w.get_zbit(q)) {
            (false, false) => 'I',
            (true, false) => 'X',
            (true, true) => 'Y',
            (false, true) => 'Z',
        })
        .collect()
}
fn main() {
    let get = |key: &str, default: &str| std::env::var(key).unwrap_or(default.to_owned());
    let n: usize = get("QUBITS", "3").parse().unwrap();
    assert!((1..=8).contains(&n));
    let depth: usize = get("STEPS", "3").parse().unwrap();
    let mode = get("MODE", "shared");
    assert!(mode == "shared" || mode == "independent");
    let gates = circuit(n, depth, mode == "shared");
    let count = if mode == "shared" { 2 } else { gates.len() };
    let pattern: PauliPattern = "Z?*".into();
    let state = symbolic(n, &gates);
    let expression = state.trace(&pattern);
    let mut max_error: f64 = 0.0;
    for sample in 0..5 {
        let values = angles(count, sample);
        let reference = numeric(n, &gates, &values);
        let evaluated: std::collections::HashMap<_, _> = state
            .iter()
            .map(|(w, c)| (word(w, n), c.eval(&values).unwrap()))
            .collect();
        for (w, c) in reference.iter() {
            max_error = max_error.max((evaluated.get(&word(w, n)).unwrap_or(&0.0) - c).abs());
        }
        for (w, c) in state.iter() {
            let ref_value = reference
                .iter()
                .find(|(k, _)| word(*k, n) == word(w, n))
                .map(|(_, v)| *v)
                .unwrap_or(0.0);
            max_error = max_error.max((c.eval(&values).unwrap() - ref_value).abs());
            if std::env::var_os("VERIFY").is_some() {
                println!(
                    "C,{sample},{},{:.17e}",
                    word(w, n),
                    c.eval(&values).unwrap()
                );
            }
        }
        max_error =
            max_error.max((expression.eval(&values).unwrap() - reference.trace(&pattern)).abs());
    }
    assert!(
        max_error < 1e-10,
        "symbolic/numeric coefficient mismatch: {max_error}"
    );
    if std::env::var_os("VERIFY").is_some() {
        return;
    }
    let trials: usize = get("TRIALS", "7").parse().unwrap();
    let batch: usize = get("BATCH", "100").parse().unwrap();
    assert!(trials > 0 && batch > 0);
    let assignments: Vec<_> = (0..batch).map(|sample| angles(count, sample)).collect();
    // All compilation and correctness checks are outside the timed region.
    for trial in 0..trials {
        let start = Instant::now();
        let fresh = symbolic(n, black_box(&gates));
        let prepared = fresh.trace(&pattern);
        let build_s = start.elapsed().as_secs_f64();
        black_box(&prepared);
        let start = Instant::now();
        for values in &assignments {
            black_box(expression.eval(black_box(values)).unwrap());
        }
        let eval_s = start.elapsed().as_secs_f64() / batch as f64;
        let start = Instant::now();
        for values in &assignments {
            black_box(numeric(n, &gates, black_box(values)).trace(&pattern));
        }
        let numeric_s = start.elapsed().as_secs_f64() / batch as f64;
        println!(
            "ppvm,{n},{depth},{mode},{trial},{build_s:.12e},{eval_s:.12e},{numeric_s:.12e},{},{:.17e},{max_error:.3e}",
            state.len(),
            expression.eval(&assignments[0]).unwrap()
        );
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn single_rotation_matches_cosine() {
        let pattern: PauliPattern = "Z?*".into();
        let gates = circuit(1, 1, false);
        let s = symbolic(1, &gates);
        assert!((s.trace(&pattern).eval(&[0.37]).unwrap() - 0.37_f64.cos()).abs() < 1e-12);
    }
    #[test]
    fn shared_parameters_match_independent_substitution() {
        let pattern: PauliPattern = "Z?*".into();
        let independent = circuit(3, 3, false);
        let values: Vec<_> = independent
            .iter()
            .map(|g| if g.bond { 0.27 } else { 0.13 })
            .collect();
        let a = symbolic(3, &independent)
            .trace(&pattern)
            .eval(&values)
            .unwrap();
        let b = symbolic(3, &circuit(3, 3, true))
            .trace(&pattern)
            .eval(&[0.13, 0.27])
            .unwrap();
        assert!((a - b).abs() < 1e-12);
    }
}
