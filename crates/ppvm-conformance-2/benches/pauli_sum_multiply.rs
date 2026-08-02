// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! **L4 operator product** benchmark (integration workload 6,
//! "observable-product-and-variance"): `O *= P` for a single Pauli word, and the
//! full `|A|·|B|` twisted convolution.
//!
//! # Why this is a SEPARATE bench target rather than a group in `pauli_sum_integration`
//!
//! Adding this group to `pauli_sum_integration` was measured to move the **old**
//! engine's numbers — `old/pauli_error_sweep` 4.74 µs → 5.37 µs (+13%) and
//! `old/rx_sweep` 6.85 µs → 6.21 µs (−9%) across 4 process runs each — even though
//! nothing in the old crate changed. That is pure code-layout/alignment
//! sensitivity (Mytkowicz): the two engines share one binary, so appending
//! functions shifts both. It would have re-based the standing `integration_trotter`
//! ratio from ~1.08 to ~1.13 with no underlying cause. Keeping the product in its
//! own binary leaves the baseline bench byte-identical, so its ratios stay
//! comparable across components: adding or removing this file must not move the
//! `pauli_sum_integration` numbers at all.
//!
//! (An earlier revision of this header quoted specific baseline ratios measured
//! *across* builds. They are deliberately gone: cross-build ratios are not
//! comparable — only a same-build, interleaved old-vs-new run is — and the quoted
//! figures did not reproduce same-build. The live baseline numbers belong in
//! `docs/log.md`, measured per run, not frozen in a module header.)
//!
//! # Shapes
//!
//! * `multiply_word` — `O *= P`: a bijective re-key through the phased product on
//!   old (`ppvm-pauli-sum/src/sum/ops.rs:95`) and through `RekeyBijective` on new.
//!   Both sides live, so this is a strict engine-to-engine ratio.
//! * `multiply_sum` — the full convolution into a fresh accumulator (`new/square`)
//!   and into the store's persistent aux double-buffer (`new/multiply_in_place`).
//!   **New only, and not by choice:** old's `impl MulAssign<PauliSum<T>>`
//!   (`ops.rs:70`) requires `PhasedPauliWord: for<'a> From<&'a T::PauliWordType>`,
//!   while `PhasedPauliWord`'s only word conversion is `impl<W: PauliWordTrait>
//!   From<W>` and `PauliWordTrait` is implemented only for `PauliWord` itself,
//!   never for a reference (`ppvm-pauli-word/src/word/data.rs:100`). The bound is
//!   unsatisfiable for every shipped `Config`, so `old_sum *= other_sum` does not
//!   compile: old's sum×sum is unreachable dead code and there is no baseline to
//!   ratio against. This records the new engine's absolute cost instead.
//!
//! The product runs over `Complex<f64>` on both sides: a Pauli product emits an
//! `iᵏ`, so old bounds its `Mul` on `ComplexCoefficient` and new bounds
//! `multiply_into` on `ImaginaryUnit`; `f64` satisfies neither. Storage is pinned
//! to `[u8; 8]` on both sides so the ratio is engine-to-engine.

use criterion::{Criterion, criterion_group, criterion_main};
use num::Complex;

use ppvm_pauli_sum::config::fxhash::Byte as OldByte;
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_pauli_word::word::PauliWord as OldWordT;
use ppvm_traits::traits::NoStrategy as OldNoStrategy;

use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeffThreshold, CombinedPolicy, HashMapStore,
    MaxPauliWeight as NewMaxWeight, NoPolicy, PauliWord as NewPauliWord, Sum,
};
use ppvm_traits_2::{
    Clifford as NewClifford, PauliError as NewPauliError, RotationOne as NewRotationOne,
};

const THRESHOLD: f64 = 1e-6;

type NewKey = NewPauliWord<[u8; 8]>;
type NewPolicy = CombinedPolicy<NewCoeffThreshold, NewMaxWeight>;
/// The `f64` sum used only to GROW a realistic support (the Trotter workload).
type NewGrowSum = Sum<HashMapStore<NewKey, f64>, NewPolicy>;

type OldCplxSum = OldPauliSum<OldByte<8, Complex<f64>, OldNoStrategy>>;
type OldCplxWord = OldWordT<[u8; 8]>;
type NewCplxSum = Sum<HashMapStore<NewKey, Complex<f64>>, NoPolicy>;

/// Grow a realistic ~1e2–1e3-term support by propagating the noisy-TFIM Trotter
/// circuit on the NEW engine (only the *keys* are used afterwards, so the old side
/// does not need to replay it). Params match `pauli_sum_integration`: `n = 12`,
/// `dt = 0.1`, 10 steps, `J = 1/8`, coefficient floor `1e-6`.
fn grown_support(n: usize) -> Vec<String> {
    let h = 1.0_f64;
    let dt = 0.1 / h;
    let steps = ((1.0 / h) / dt) as usize;
    let theta_x = dt * h;
    let theta_zz = dt * (h / 8.0);
    let noise = [1e-4 / 4.0; 3];

    let policy = CombinedPolicy(
        NewCoeffThreshold {
            threshold: THRESHOLD,
        },
        NewMaxWeight(usize::MAX),
    );
    let mut state = NewGrowSum::from_terms_with_policy(
        n,
        policy,
        (0..n).map(|i| {
            let s: String = (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect();
            (NewKey::from(s.as_str()), 1.0)
        }),
    );
    for _ in 0..steps {
        for i in 0..n {
            state.pauli_error(i, noise);
            state.truncate();
            state.rx(i, theta_x);
            state.truncate();
        }
        for i in 0..n - 1 {
            state.pauli_error(i + 1, noise);
            state.truncate();
            state.pauli_error(i, noise);
            state.truncate();
            state.cnot(i, i + 1);
            state.rz(i + 1, theta_zz);
            state.cnot(i, i + 1);
            state.truncate();
        }
    }
    let mut keys: Vec<String> = state.iter().map(|(k, _)| k.to_string()).collect();
    keys.sort();
    keys
}

/// Deterministic complex coefficients over the grown keys — no RNG, so both sides
/// and successive process runs see byte-identical input.
fn complex_terms(keys: &[String], cap: usize) -> Vec<(String, Complex<f64>)> {
    keys.iter()
        .take(cap)
        .enumerate()
        .map(|(i, k)| {
            let t = (i % 17) as f64;
            (k.clone(), Complex::new(0.5 + t / 32.0, -0.25 + t / 64.0))
        })
        .collect()
}

fn build_old(n: usize, terms: &[(String, Complex<f64>)]) -> OldCplxSum {
    let mut s: OldCplxSum = OldPauliSum::builder()
        .n_qubits(n)
        .capacity(n.pow(2))
        .build();
    for (w, v) in terms {
        s += (w.as_str(), *v);
    }
    s
}

fn build_new(n: usize, terms: &[(String, Complex<f64>)]) -> NewCplxSum {
    // Same explicit `n²` capacity override `build_old` passes, so the ratio is
    // engine-to-engine rather than an artifact of differently sized maps.
    let mut s = NewCplxSum::with_capacity(n, NoPolicy, n.pow(2));
    for (w, v) in terms {
        s += (NewKey::from(w.as_str()), *v);
    }
    s
}

fn bench_multiply(c: &mut Criterion) {
    let n = 12usize;
    let keys = grown_support(n);

    // --- (i) Right-multiply by a single Pauli word, over the full support. ----
    let terms = complex_terms(&keys, usize::MAX);
    println!(
        "[pauli_sum/multiply_word] benchmarked support: {} terms",
        terms.len()
    );
    let mut old_cplx = build_old(n, &terms);
    let mut new_cplx = build_new(n, &terms);

    let word: String = (0..n).map(|i| ['X', 'Y', 'Z', 'I'][i % 4]).collect();
    // Old's `MulAssign<PauliWord>` takes the word BY VALUE, but the old word is
    // `Copy`, so this is a register move, not an allocation: both sides pay the
    // same per-call word cost.
    let old_word = OldCplxWord::from(word.as_str());
    let new_word = NewKey::from(word.as_str());

    let mut g = c.benchmark_group("pauli_sum/multiply_word");
    // `P` is an involution (`P·P = +I`, `phaseExpN_self`), so repeated application
    // cycles the state with period 2 — bit-exact, no drift — and needs no
    // per-iteration clone on either side.
    g.bench_function("new/mul_word", |b| {
        b.iter(|| new_cplx.mul_word_assign(&new_word))
    });
    g.bench_function("old/mul_word", |b| b.iter(|| old_cplx *= old_word));
    g.finish();

    // --- (ii) Full sum × sum (the ⟨O²⟩ variance shape). ----------------------
    // Quadratic in the support, so cap it: 256 terms is ~65k monomial products per
    // iteration — enough for map growth and probe behaviour to dominate without a
    // multi-second sample.
    let small = complex_terms(&keys, 256);
    let a = build_new(n, &small);
    println!(
        "[pauli_sum/multiply_sum] benchmarked support: {} terms ({} pair products)",
        a.len(),
        a.len() * a.len()
    );

    let mut g = c.benchmark_group("pauli_sum/multiply_sum");
    // Allocating form: a fresh accumulator per product, then the L3 pairing.
    g.bench_function("new/square_and_overlap", |b| {
        b.iter(|| {
            let sq = a.multiply(&a);
            sq.overlap(&a)
        })
    });
    // In-place form: the accumulator is drawn from the store's persistent aux
    // double-buffer instead of a fresh map. `iter_batched_ref` because `*=` is not
    // state-preserving; criterion does not time the setup clone.
    g.bench_function("new/multiply_in_place", |b| {
        b.iter_batched_ref(
            || a.clone(),
            |s| s.multiply_in_place(&a),
            criterion::BatchSize::LargeInput,
        )
    });
    g.finish();
}

criterion_group!(benches, bench_multiply);
criterion_main!(benches);
