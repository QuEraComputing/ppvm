// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Phase-5 perf gate for the symbolic coefficient ring: new `ppvm-sym-2` vs old
//! `ppvm-sym`, **both engines in one binary** so every reported ratio is
//! same-build.
//!
//! # Why the gate is the integration workload, not the microbench
//!
//! A tight `b.iter(|| a.clone() * b.clone())` loop lets the allocator recycle one
//! warm page, so it cannot see cumulative per-gate costs — allocation churn,
//! buffer reuse, a lost fast-path arm. That is exactly how the `mem::take`
//! zero-capacity rebuild in old's `Sum::mul_term` (perf feature 5) stays
//! invisible: it costs one realloc+rehash *per coefficient per gate*, which a
//! one-multiply microbench barely registers and a deep circuit multiplies by
//! `|support| × n_gates`. So the HEADLINE numbers here are the end-to-end
//! `sym.*` workloads; the microbenches at the bottom exist only for attribution.
//!
//! # Fair comparison
//!
//! Both sides run the **same algebraic configuration** (`ppvm_conformance_2::sym`
//! module docs): `[u8; 8]` key storage, the `Term` coefficient of the respective
//! crate, an FxHash-class seed-free monomial hasher, and no sum-level truncation
//! strategy — truncation on this path is intrinsic to the coefficient
//! (`max_sin`/`min_eps` seeded on the initial observable), identically on both
//! sides. The gate sequence is literally shared code: `bench_trotter_old` /
//! `bench_trotter_new`
//! and `replay_old`/`replay_new` are written against one circuit description, and
//! `rzz` is decomposed as `cnot; rz; cnot` on both sides rather than calling old's
//! built-in.
//!
//! # Mined from real usage
//!
//! `ppvm-sym` ships **no benches of its own**, so the workloads are derived from
//! how it is actually driven downstream: `examples/symbolic.rs` (the parametric
//! trace), `examples/tfim.rs` scaled with the layer/step structure of
//! `ppvm-pauli-sum/benches/trotter.rs`, the sweep shape of
//! `ppvm-pauli-sum/benches/truncation-weight.rs`, and the replay shape of
//! `ppvm-pauli-sum/benches/random-circuit.rs`. A symbolic-coefficient
//! `PauliSum` propagation is the realistic driver in every case.

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use std::collections::BTreeMap;

use ppvm_conformance_2::seeded_rng;
use ppvm_conformance_2::sym::*;

use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, PauliWord as NewPauliWord, Sum as NewSum};
use ppvm_sym::{Prod as OldProd, Sum as OldSumTy, Term as OldTerm2};
use ppvm_sym_2::{GaussianInt, Prod as NewProd, Sum as NewSumTy, Term as NewTerm2};

mod sym_surface;

fn fixed_angles(n: usize) -> Vec<f64> {
    (0..n).map(|i| 0.23 + 0.17 * i as f64).collect()
}

#[track_caller]
fn assert_terms_match(old: &OldTerm2, new: &NewTerm2, angles: &[&[f64]], label: &str) {
    let old_view = old_view(old);
    let new_view = new_view(new);
    assert_eq!(old_view.form, new_view.form, "[{label}] form differs");
    assert_eq!(old_view.c0, new_view.c0, "[{label}] constant differs");
    assert_eq!(
        old_view.monomials, new_view.monomials,
        "[{label}] canonical monomial coefficients differ"
    );
    for values in angles {
        let old_value = old.eval(values).unwrap();
        let new_value = new.eval(values).unwrap();
        assert!(
            (old_value - new_value).abs() < 1e-12,
            "[{label}] evaluated output differs at {values:?}: old={old_value}, new={new_value}"
        );
    }
}

#[track_caller]
fn assert_sym_sums_match(old: &OldSymSum, new: &NewSymSum, angles: &[&[f64]], label: &str) {
    let old = old_sym_support(old);
    let new = new_sym_support(new);
    assert_eq!(
        old.iter().map(|(key, _)| key).collect::<Vec<_>>(),
        new.iter().map(|(key, _)| key).collect::<Vec<_>>(),
        "[{label}] canonical Pauli support differs"
    );
    for ((key, old_coeff), (_, new_coeff)) in old.iter().zip(&new) {
        assert_terms_match(
            old_coeff,
            new_coeff,
            angles,
            &format!("{label} coefficient at {key}"),
        );
    }
}

// ---------------------------------------------------------------------------
// 1. HEADLINE — `sym.tfim.trotter`: the deep symbolic Trotter propagation.
// ---------------------------------------------------------------------------

/// The headline spec: `n = 8`, `L = 10` layers (160 rotation gates + 140 `rzz`
/// decompositions = 580 gate applications), a FRESH symbolic variable per layer
/// per gate family so the monomial space genuinely grows with depth.
fn headline(max_sin: usize) -> TrotterSpec {
    TrotterSpec {
        n: 8,
        layers: 10,
        max_sin,
        min_eps: 1e-12,
        observable: "ZIIIIIII",
    }
}

fn bench_trotter_old(spec: &TrotterSpec) -> OldSymSum {
    let mut sum = sym_surface::old_sum(spec.n);
    sum += (spec.observable, old_seed_coeff(spec.max_sin, spec.min_eps));
    for layer in 0..spec.layers {
        for q in 0..spec.n {
            ppvm_traits::traits::RotationOne::rx(&mut sum, q, OldTerm2::var(2 * layer));
        }
        for q in 0..spec.n - 1 {
            ppvm_traits::traits::Clifford::cnot(&mut sum, q, q + 1);
            ppvm_traits::traits::RotationOne::rz(&mut sum, q + 1, OldTerm2::var(2 * layer + 1));
            ppvm_traits::traits::Clifford::cnot(&mut sum, q, q + 1);
        }
    }
    sum
}

fn bench_trotter_new(spec: &TrotterSpec) -> NewSymSum {
    let mut sum = sym_surface::new_sum(spec.n);
    sum += (
        NewSymKey::from(spec.observable),
        new_seed_coeff(spec.max_sin, spec.min_eps),
    );
    for layer in 0..spec.layers {
        for q in 0..spec.n {
            ppvm_traits_2::RotationOne::rx(&mut sum, q, NewTerm2::var(2 * layer));
        }
        for q in 0..spec.n - 1 {
            ppvm_traits_2::Clifford::cnot(&mut sum, q, q + 1);
            ppvm_traits_2::RotationOne::rz(&mut sum, q + 1, NewTerm2::var(2 * layer + 1));
            ppvm_traits_2::Clifford::cnot(&mut sum, q, q + 1);
        }
    }
    sum
}

fn bench_tfim_trotter(c: &mut Criterion) {
    for k in [3usize, 4] {
        let spec = headline(k);

        // Diagnostic counters, printed once so a future change that silently
        // shrinks (or explodes) the workload is visible in the bench output.
        let ns = bench_trotter_new(&spec);
        let os = bench_trotter_old(&spec);
        let angles = fixed_angles(spec.n_vars());
        assert_sym_sums_match(&os, &ns, &[&angles], &format!("tfim trotter k={k}"));
        let total: usize = new_sym_support(&ns)
            .iter()
            .map(|(_, t)| new_view(t).n_monomials())
            .sum();
        let peak: usize = new_sym_support(&ns)
            .iter()
            .map(|(_, t)| new_view(t).n_monomials())
            .max()
            .unwrap_or(0);
        println!(
            "[sym/tfim_trotter k={k}] support={} monomials={total} peak_monomials={peak}",
            ns.len()
        );

        let mut g = c.benchmark_group(format!("sym/tfim_trotter_k{k}"));
        // The timed body is the WHOLE propagation (seed + every gate), which is
        // where the per-coefficient-per-gate allocation costs live.
        g.bench_function("new/trotter", |b| b.iter(|| bench_trotter_new(&spec)));
        g.bench_function("old/trotter", |b| b.iter(|| bench_trotter_old(&spec)));
        g.finish();
    }
}

// ---------------------------------------------------------------------------
// 2. `sym.trace.parametric` — build → 9 gates → trace → eval.
// ---------------------------------------------------------------------------

fn bench_trace_parametric(c: &mut Criterion) {
    let old = parametric_trace_old();
    let new = parametric_trace_new();
    assert_terms_match(&old, &new, &[&[1.1, 2.1]], "parametric trace");

    let mut g = c.benchmark_group("sym/trace_parametric");
    g.bench_function("new/build_propagate_trace_eval", |b| {
        b.iter(|| {
            let t = parametric_trace_new();
            t.eval(&[1.1, 2.1]).unwrap()
        })
    });
    g.bench_function("old/build_propagate_trace_eval", |b| {
        b.iter(|| {
            let t = parametric_trace_old();
            t.eval(&[1.1, 2.1]).unwrap()
        })
    });
    g.finish();
}

// ---------------------------------------------------------------------------
// 3. `sym.truncation.sweep` — the cost CURVE in `max_sin` (and `min_eps`).
// ---------------------------------------------------------------------------
//
// The SHAPE is the gate, not any one point: new must match-or-beat old at every
// `k`, and the growth in `k` must not be steeper for new.

fn bench_truncation_sweep(c: &mut Criterion) {
    let base = TrotterSpec {
        n: 6,
        layers: 6,
        max_sin: 0,
        min_eps: 1e-12,
        observable: "ZIIIII",
    };

    let mut g = c.benchmark_group("sym/truncation_sweep");
    for k in 1..=5usize {
        let spec = TrotterSpec { max_sin: k, ..base };
        let ns = bench_trotter_new(&spec);
        let os = bench_trotter_old(&spec);
        let angles = fixed_angles(spec.n_vars());
        assert_sym_sums_match(&os, &ns, &[&angles], &format!("truncation sweep k={k}"));
        let nm: usize = new_sym_support(&ns)
            .iter()
            .map(|(_, t)| new_view(t).n_monomials())
            .sum();
        let om: usize = old_sym_support(&os)
            .iter()
            .map(|(_, t)| old_view(t).n_monomials())
            .sum();
        println!("[sym/truncation_sweep k={k}] retained monomials: new={nm} old={om}");
        assert_eq!(
            nm, om,
            "canonical comparison must retain equal counts at k={k}"
        );

        g.bench_function(format!("new/k{k}"), |b| b.iter(|| bench_trotter_new(&spec)));
        g.bench_function(format!("old/k{k}"), |b| b.iter(|| bench_trotter_old(&spec)));
    }
    g.finish();

    let mut g = c.benchmark_group("sym/truncation_sweep_eps");
    for (label, min_eps) in [("eps", f64::EPSILON), ("1e-12", 1e-12), ("1e-6", 1e-6)] {
        let spec = TrotterSpec {
            max_sin: 3,
            min_eps,
            ..base
        };
        let ns = bench_trotter_new(&spec);
        let os = bench_trotter_old(&spec);
        let angles = fixed_angles(spec.n_vars());
        assert_sym_sums_match(&os, &ns, &[&angles], &format!("truncation sweep {label}"));
        g.bench_function(format!("new/{label}"), |b| {
            b.iter(|| bench_trotter_new(&spec))
        });
        g.bench_function(format!("old/{label}"), |b| {
            b.iter(|| bench_trotter_old(&spec))
        });
    }
    g.finish();
}

// ---------------------------------------------------------------------------
// 4. `sym.expectation.grid` — the VQE-shaped read-out.
// ---------------------------------------------------------------------------
//
// Two separate numbers: (a) the one-shot propagation + trace, and (b) `eval`
// throughput over a 1000-point angle grid. (b) is where the packed-monomial
// layout should show a clear win over old's two `BTreeMap` walks per monomial;
// if it does not, the layout change bought nothing on the read-out path.

fn bench_expectation_grid(c: &mut Criterion) {
    let spec = TrotterSpec {
        n: 6,
        layers: 6,
        max_sin: 3,
        min_eps: 1e-12,
        observable: "ZIIIII",
    };
    let mut rng = seeded_rng(0x9_81D);
    let grid = angle_grid(&mut rng, spec.n_vars(), 1000);
    let validation_angles: Vec<&[f64]> = grid.iter().map(Vec::as_slice).collect();

    let propagated_new = bench_trotter_new(&spec);
    let propagated_old = bench_trotter_old(&spec);
    assert_sym_sums_match(
        &propagated_old,
        &propagated_new,
        &[validation_angles[0]],
        "expectation propagation",
    );

    // (a) one-shot propagation + trace.
    let mut g = c.benchmark_group("sym/expectation_propagate");
    g.bench_function("new/propagate_trace", |b| {
        b.iter(|| {
            let s = bench_trotter_new(&spec);
            ppvm_traits_2::Trace::trace(&s, &ppvm_pauli_sum_2::PauliPattern::zero_state())
        })
    });
    g.bench_function("old/propagate_trace", |b| {
        b.iter(|| {
            let s = bench_trotter_old(&spec);
            ppvm_traits::traits::Trace::trace(
                &s,
                &ppvm_pauli_word::pattern::PauliPattern::from("Z?*"),
            )
        })
    });
    g.finish();

    // (b) `eval` throughput over the grid, on the ALREADY propagated trace.
    let nt = ppvm_traits_2::Trace::trace(
        &propagated_new,
        &ppvm_pauli_sum_2::PauliPattern::zero_state(),
    );
    let ot = ppvm_traits::traits::Trace::trace(
        &propagated_old,
        &ppvm_pauli_word::pattern::PauliPattern::from("Z?*"),
    );
    assert_terms_match(&ot, &nt, &validation_angles, "expectation trace");
    println!(
        "[sym/expectation_eval] trace monomials: new={} old={}",
        new_view(&nt).n_monomials(),
        old_view(&ot).n_monomials()
    );

    let mut g = c.benchmark_group("sym/expectation_eval");
    g.throughput(criterion::Throughput::Elements(grid.len() as u64));
    g.bench_function("new/eval_grid_1000", |b| {
        b.iter(|| {
            let mut acc = 0.0f64;
            for v in &grid {
                acc += nt.eval(v).unwrap();
            }
            acc
        })
    });
    g.bench_function("old/eval_grid_1000", |b| {
        b.iter(|| {
            let mut acc = 0.0f64;
            for v in &grid {
                acc += ot.eval(v).unwrap();
            }
            acc
        })
    });
    g.finish();
}

// ---------------------------------------------------------------------------
// 5. `sym.random.circuit` — deep heterogeneous replay, with a per-phase split.
// ---------------------------------------------------------------------------

fn bench_random_circuit(c: &mut Criterion) {
    let n = 8usize;
    let depth = 200usize;
    let n_vars = 8u32;
    let circuit = random_sym_circuit(&mut seeded_rng(0xC0FFEE), n, depth, n_vars);

    // The Clifford-only prefix (cheap re-key, coefficient untouched) and the
    // rotation-heavy remainder, so a regression can be attributed to one phase
    // without a separate microbench.
    let clifford: Vec<SymGate> = circuit
        .iter()
        .copied()
        .filter(|g| {
            matches!(
                g,
                SymGate::H(_) | SymGate::S(_) | SymGate::Cnot(..) | SymGate::Cz(..)
            )
        })
        .collect();
    let rotations: Vec<SymGate> = circuit
        .iter()
        .copied()
        .filter(|g| !clifford.contains(g))
        .collect();
    println!(
        "[sym/random_circuit] {} gates ({} Clifford, {} rotation)",
        circuit.len(),
        clifford.len(),
        rotations.len()
    );

    let seed_old = |s: &mut OldSymSum| *s += ("ZIIIIIII", old_seed_coeff(3, 1e-12));
    let seed_new =
        |s: &mut NewSymSum| *s += (NewSymKey::from("ZIIIIIII"), new_seed_coeff(3, 1e-12));

    let mut old_probe = sym_surface::old_sum(n);
    let mut new_probe = sym_surface::new_sum(n);
    seed_old(&mut old_probe);
    seed_new(&mut new_probe);
    replay_old(&mut old_probe, &circuit);
    replay_new(&mut new_probe, &circuit);
    let angles = fixed_angles(n_vars as usize);
    assert_sym_sums_match(&old_probe, &new_probe, &[&angles], "random circuit");

    let mut old_clifford_probe = sym_surface::old_sum(n);
    let mut new_clifford_probe = sym_surface::new_sum(n);
    seed_old(&mut old_clifford_probe);
    seed_new(&mut new_clifford_probe);
    replay_old(&mut old_clifford_probe, &clifford);
    replay_new(&mut new_clifford_probe, &clifford);
    assert_sym_sums_match(
        &old_clifford_probe,
        &new_clifford_probe,
        &[&angles],
        "random circuit Clifford phase",
    );

    let mut g = c.benchmark_group("sym/random_circuit");
    g.sample_size(10);
    g.measurement_time(std::time::Duration::from_secs(10));
    g.bench_function("new/full_replay", |b| {
        b.iter_batched_ref(
            || {
                let mut s = sym_surface::new_sum(n);
                seed_new(&mut s);
                s
            },
            |s| replay_new(s, &circuit),
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/full_replay", |b| {
        b.iter_batched_ref(
            || {
                let mut s = sym_surface::old_sum(n);
                seed_old(&mut s);
                s
            },
            |s| replay_old(s, &circuit),
            BatchSize::SmallInput,
        )
    });
    g.finish();

    // The Clifford-only phase (pure re-key; the coefficient is never touched).
    let mut g = c.benchmark_group("sym/random_circuit_clifford");
    g.bench_function("new/clifford_prefix", |b| {
        b.iter_batched_ref(
            || {
                let mut s = sym_surface::new_sum(n);
                seed_new(&mut s);
                s
            },
            |s| replay_new(s, &clifford),
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/clifford_prefix", |b| {
        b.iter_batched_ref(
            || {
                let mut s = sym_surface::old_sum(n);
                seed_old(&mut s);
                s
            },
            |s| replay_old(s, &clifford),
            BatchSize::SmallInput,
        )
    });
    g.finish();
}

// ---------------------------------------------------------------------------
// 5b. `sym.random.circuit` — the TRACE READ-OUT in isolation, swept over the
//     seeded `max_sin`.
// ---------------------------------------------------------------------------
//
// The propagation benches above time the write side; this one times only
// `trace(&pattern)` on an already-propagated sum. It exists because a read side
// that clones every coefficient before the pattern filter looks at it is
// invisible in a propagation total (the clone is a few percent of a gate) and
// invisible in a single-op microbench (one coefficient, warm allocator) — yet it
// cost 7×–33× against old here, growing with `max_sin` because the coefficients
// grow with it. The sweep IS the gate: the ratio must not grow with `k`. Only
// ~255 of the ~65k keys match `Z?*`, so a cloning read side does ~99.6% waste.
fn bench_trace_readout(c: &mut Criterion) {
    let n = 8usize;
    let circuit = random_sym_circuit(&mut seeded_rng(0xC0FFEE), n, 200, 8);
    let old_pat = ppvm_pauli_word::pattern::PauliPattern::from("Z?*");
    let new_pat = ppvm_pauli_sum_2::PauliPattern::zero_state();

    for k in [2usize, 3, 4, 5] {
        let mut os = sym_surface::old_sum(n);
        os += ("ZIIIIIII", old_seed_coeff(k, 1e-12));
        replay_old(&mut os, &circuit);

        let mut ns = sym_surface::new_sum(n);
        ns += (NewSymKey::from("ZIIIIIII"), new_seed_coeff(k, 1e-12));
        replay_new(&mut ns, &circuit);

        let angles = fixed_angles(8);
        assert_sym_sums_match(&os, &ns, &[&angles], &format!("trace readout input k={k}"));
        let old_trace = ppvm_traits::traits::Trace::trace(&os, &old_pat);
        let new_trace = ppvm_traits_2::Trace::trace(&ns, &new_pat);
        assert_terms_match(
            &old_trace,
            &new_trace,
            &[&angles],
            &format!("trace readout k={k}"),
        );

        let matching = new_sym_support(&ns)
            .iter()
            .filter(|(key, _)| key.chars().all(|ch| ch == 'I' || ch == 'Z'))
            .count();
        println!(
            "[sym/trace_readout k={k}] support={} matching={matching}",
            ns.len()
        );

        let mut g = c.benchmark_group(format!("sym/trace_readout_k{k}"));
        g.bench_function("new/trace", |b| {
            b.iter(|| ppvm_traits_2::Trace::trace(&ns, &new_pat))
        });
        g.bench_function("old/trace", |b| {
            b.iter(|| ppvm_traits::traits::Trace::trace(&os, &old_pat))
        });
        g.finish();
    }
}

// ---------------------------------------------------------------------------
// 6. `sym.exact.multiply` — the L4 twisted product over the EXACT ring.
// ---------------------------------------------------------------------------
//
// There is no honest OLD baseline for this path: old's phase handling is broken
// in three independent places and its `eval` discards the phase, so its numbers
// on the operator product are not a target (numeric acceptance is Lean, in
// `tests/sym_lean.rs`). The comparison here is therefore new-vs-new: the exact
// `ℤ[i]` ring against the shipped `Complex<f64>`, on the SAME support, to bound
// the constant factor the exactness costs.

type ExactSum = NewSum<HashMapStore<NewPauliWord<[u8; 8]>, GaussianInt>, NoPolicy>;
type ComplexSum = NewSum<HashMapStore<NewPauliWord<[u8; 8]>, num::Complex<f64>>, NoPolicy>;

fn bench_exact_multiply(c: &mut Criterion) {
    use ppvm_conformance_2::random_pauli_string;

    let n = 6usize;
    let terms = 256usize;
    let mut rng = seeded_rng(0xE_AC7);
    let words: Vec<String> = (0..terms)
        .map(|_| random_pauli_string(&mut rng, n))
        .collect();

    let mut ea: ExactSum = ExactSum::new(n);
    let mut eb: ExactSum = ExactSum::new(n);
    let mut ca: ComplexSum = ComplexSum::new(n);
    let mut cb: ComplexSum = ComplexSum::new(n);
    for (i, w) in words.iter().enumerate() {
        let k = NewPauliWord::<[u8; 8]>::from(w.as_str());
        let (re, im) = ((i % 7) as i64 - 3, (i % 5) as i64 - 2);
        if i % 2 == 0 {
            ea += (k, GaussianInt::new(re, im));
            ca += (k, num::Complex::new(re as f64, im as f64));
        } else {
            eb += (k, GaussianInt::new(re, im));
            cb += (k, num::Complex::new(re as f64, im as f64));
        }
    }
    println!(
        "[sym/exact_multiply] |A|={} |B|={} → {} key pairs",
        ea.len(),
        eb.len(),
        ea.len() * eb.len()
    );

    // `multiply_into` with a REUSED destination: the L4 contract is that the
    // twisted convolution accumulates into an existing store, so the timed body
    // must not pay a fresh allocation per key pair.
    let mut eacc: ExactSum = ExactSum::new(n);
    let mut cacc: ComplexSum = ComplexSum::new(n);
    ea.multiply_into(&eb, &mut eacc);
    ca.multiply_into(&cb, &mut cacc);
    let exact: BTreeMap<_, _> = eacc
        .iter()
        .map(|(key, coeff)| (key.to_string(), (coeff.re as f64, coeff.im as f64)))
        .collect();
    let complex: BTreeMap<_, _> = cacc
        .iter()
        .map(|(key, coeff)| (key.to_string(), (coeff.re, coeff.im)))
        .collect();
    assert_eq!(
        exact, complex,
        "exact and complex multiplication outputs differ"
    );

    let mut g = c.benchmark_group("sym/exact_multiply");
    g.bench_function("exact_zi/multiply_into", |b| {
        b.iter(|| ea.multiply_into(&eb, &mut eacc))
    });
    g.bench_function("complex_f64/multiply_into", |b| {
        b.iter(|| ca.multiply_into(&cb, &mut cacc))
    });
    g.finish();
}

// ---------------------------------------------------------------------------
// 7. DIAGNOSTIC microbenches (attribution only — NOT the gate).
// ---------------------------------------------------------------------------
//
// These cannot see cumulative per-gate cost (a tight `iter` loop recycles one
// warm allocator page), so a "parity" reading here means nothing on its own.
// They exist to attribute a movement in the integration numbers to one of the
// four hot primitives: the monomial product, `Sum × One` (`mul_term`, where the
// aux double-buffer lives), `Sum × Sum` (the `O(|s1|·|s2|)` loop), and `eval`.

fn build_old_sum_terms(n: usize) -> OldTerm2 {
    let mut t = OldTerm2::from(1.0);
    for v in 0..n as u32 {
        t += OldTerm2::var(v).sin();
        t += OldTerm2::var(v).cos();
    }
    t
}

fn build_new_sum_terms(n: usize) -> NewTerm2 {
    let mut t = NewTerm2::from(1.0);
    for v in 0..n as u32 {
        t += NewTerm2::var(v).sin();
        t += NewTerm2::var(v).cos();
    }
    t
}

fn bench_micro(c: &mut Criterion) {
    // --- the monomial product (`Prod × Prod`) --------------------------------
    let mut op = OldProd::sin(0);
    let mut np = NewProd::sin(0);
    for v in 1..6u32 {
        op.mul_sin(v);
        op.mul_cos(v);
        np.mul_sin(v);
        np.mul_cos(v);
    }
    let oq = OldProd::sin(7);
    let nq = NewProd::sin(7);

    let mut g = c.benchmark_group("sym/micro_prod_mul");
    g.bench_function("new/prod_mul", |b| {
        b.iter(|| std::hint::black_box(np.clone() * nq.clone()))
    });
    g.bench_function("old/prod_mul", |b| {
        b.iter(|| std::hint::black_box(op.clone() * oq.clone()))
    });
    g.finish();

    // --- `Sum × One` (`mul_term`): the aux double-buffer path ----------------
    // The monomial table is rebuilt through the multiply; old's `mem::take`
    // leaves a zero-capacity map, new ping-pongs a persistent `aux`.
    let mut os = OldSumTy::new();
    let mut ns = NewSumTy::new();
    for v in 0..24u32 {
        os.add_term(OldProd::sin(v), 1.0 + v as f64, usize::MAX, f64::EPSILON);
        ns.add_term(NewProd::sin(v), 1.0 + v as f64, usize::MAX, f64::EPSILON);
    }
    let mut g = c.benchmark_group("sym/micro_mul_term");
    g.bench_function("new/mul_term", |b| {
        b.iter_batched_ref(
            || ns.clone(),
            |s| s.mul_term(NewProd::cos(30), 1.5, usize::MAX, f64::EPSILON),
            BatchSize::SmallInput,
        )
    });
    g.bench_function("old/mul_term", |b| {
        b.iter_batched_ref(
            || os.clone(),
            |s| s.mul_term(OldProd::cos(30), 1.5, usize::MAX, f64::EPSILON),
            BatchSize::SmallInput,
        )
    });
    g.finish();

    // --- `Term × Term` on two map-backed sums (the `O(|s1|·|s2|)` loop) ------
    let ot = build_old_sum_terms(8);
    let nt = build_new_sum_terms(8);
    let mut g = c.benchmark_group("sym/micro_term_mul");
    g.bench_function("new/sum_x_sum", |b| {
        b.iter(|| std::hint::black_box(nt.clone() * nt.clone()))
    });
    g.bench_function("old/sum_x_sum", |b| {
        b.iter(|| std::hint::black_box(ot.clone() * ot.clone()))
    });
    g.finish();

    // --- `Term + Term` -------------------------------------------------------
    let mut g = c.benchmark_group("sym/micro_term_add");
    g.bench_function("new/sum_plus_sum", |b| {
        b.iter(|| std::hint::black_box(nt.clone() + nt.clone()))
    });
    g.bench_function("old/sum_plus_sum", |b| {
        b.iter(|| std::hint::black_box(ot.clone() + ot.clone()))
    });
    g.finish();

    // --- `eval` on a single term --------------------------------------------
    let vals: Vec<f64> = (0..8).map(|i| 0.3 + 0.17 * i as f64).collect();
    let mut g = c.benchmark_group("sym/micro_eval");
    g.bench_function("new/eval", |b| {
        b.iter(|| std::hint::black_box(nt.eval(&vals).unwrap()))
    });
    g.bench_function("old/eval", |b| {
        b.iter(|| std::hint::black_box(ot.eval(&vals).unwrap()))
    });
    g.finish();
}

criterion_group!(
    benches,
    bench_tfim_trotter,
    bench_trace_parametric,
    bench_truncation_sweep,
    bench_expectation_grid,
    bench_random_circuit,
    bench_trace_readout,
    bench_exact_multiply,
    bench_micro,
    sym_surface::bench,
);
criterion_main!(benches);
