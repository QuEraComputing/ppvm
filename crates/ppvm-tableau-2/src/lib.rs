// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `ppvm-tableau-2`: the stabilizer tableau rebuilt on the `-2` trait tower.
//!
//! The two top-level types are [`Tableau`] (the pure Clifford frame) and
//! [`GeneralizedTableau`] (frame + sparse amplitude vector, so non-Clifford
//! gates are representable).
//!
//! Design: `docs/design/traits-2-configuration-and-hashing.md`
//! §"Pauli algebra traits", §"A third instantiation: the generalized tableau"
//! and §"Tableau indexability"; `docs/design/traits-2-implementation-plan.md`
//! Phase 4. Machine-checked semantics: `lean/PPVM/Tableau/Frame.lean`
//! (symplectic frame, measurement dichotomy), `lean/PPVM/Tableau/Bitstring.lean`
//! (the XOR relabel is a bijection), `lean/PPVM/Pauli/Symplectic.lean` and
//! `lean/PPVM/Pauli/Conjugation.lean` (per-generator bit and sign rules),
//! `lean/PPVM/Algebra/Truncation.lean` (`cutoff_mismatch`, the tableau's strict
//! `>` keep-rule vs `CoefficientThreshold`'s `>=`).
//!
//! The *algebra* — every Clifford kernel's bit and sign rules, the `u64`-packed
//! sort-merge branch coalesce, the fused batch mask sweeps, the `cz_block`
//! family, and the reusable measurement scratch — is ported from
//! `ppvm-tableau`, so values, visit order and RNG consumption stay identical.
//! The *layout* is not: old stored one compile-time-sized `BitArray` pair per
//! generator, which made every gate a strided walk over all `2n` of them and
//! capped the qubit count at the storage width. The frame now lives in one
//! aligned contiguous column-major allocation (`storage::TableauData`), where a
//! one-qubit Clifford is two `n`-bit sweeps and `n` is a runtime value. See
//! [`Tableau`]'s "Storage" section and `docs/design/tableau-data-structure.md`.
//!
//! # Quick example
//!
//! Prepare a Bell pair and verify the two measurements are perfectly
//! correlated:
//!
//! ```
//! use ppvm_tableau_2::prelude::*;
//! use rand::SeedableRng;
//!
//! let mut tab: GeneralizedTableau = GeneralizedTableau::new(2, 1e-12);
//! let mut rng = rand::rngs::SmallRng::seed_from_u64(0);
//! tab.h(0);
//! tab.cnot(0, 1);
//!
//! let r0 = tab.measure(0, &mut rng);
//! let r1 = tab.measure(1, &mut rng);
//! assert_eq!(r0, r1);
//! ```
//!
//! # Truncation contract
//!
//! Unlike `ppvm-pauli-sum-2`, where truncation is caller-driven, the tableau's
//! **gates auto-truncate**: the magnitude cutoff is applied inline while the new
//! amplitude vector is written, and there is no `truncate()` entry point. Three
//! cutoff semantics coexist, all reproduced verbatim from the old crate:
//!
//! | site | keep-rule |
//! |:--|:--|
//! | gates (`t`, `t_dag`, `rotate_1`, the relabel) | `\|c\|² > threshold²` (absolute, **strict**) |
//! | `rotate_2` merge | `\|c\| > \|threshold\|` (absolute magnitude, **strict**) |
//! | case-a measurement | `\|c\|² > threshold² · ‖v‖²` (**relative** to the current norm) |
//! | case-b measurement | *no magnitude filter at all* |
//!
//! Only **two** of those are mathematically distinct rules: for a non-negative
//! threshold `|c|² > t²` and `|c| > |t|` are the same predicate
//! (`cutoff_abs_iff_sq` in `lean/PPVM/Algebra/Truncation.lean`), so `rotate_2`'s
//! spelling differs from the gates' only in float rounding and in paying a
//! `hypot` per element. The genuine split is absolute versus **relative** — and
//! the relative form is the one the error bounds (`l1_bound`,
//! `l2_bound_normalized`) are stated for. All four sites are reproduced verbatim
//! from old regardless.
//!
//! The strict `>` differs from `ppvm-pauli-sum-2`'s `CoefficientThreshold`, whose
//! keep-rule is `magnitude() >= threshold`; the two disagree exactly at the
//! boundary, machine-checked in `lean/PPVM/Algebra/Truncation.lean`
//! (`cutoff_mismatch`). Reusing `CoefficientThreshold` here would silently flip
//! that boundary, so the amplitude store carries the tableau's rule inline
//! instead.
//!
//! # Deferrals
//!
//! Items deliberately **not** ported in this component, with the reason:
//!
//! 1. **`trace(&PauliPattern)` — resolved in Phase 7.** The `-2` pattern now
//!    exists in `ppvm-pauli-sum-2`; [`GeneralizedTableau::trace`] enumerates its
//!    accepted words and delegates each leaf to the audited expectation kernel.
//! 2. **The `rayon` feature.** The old crate gated an opt-in parallel branch
//!    path on `RAYON_COEFF_THRESHOLD = 16384` coefficients with a
//!    `rayon::current_thread_index().is_none()` nesting guard. The default build
//!    never took it (below the threshold the sequential path always wins), so
//!    omitting it changes no default behaviour; it is recorded here rather than
//!    dropped silently, and the sequential kernels are written so the parallel
//!    map can be reinstated around them unchanged.
//! 3. **A generic real coefficient.** The old crate was generic in `T::Coeff`,
//!    but every shipped configuration instantiated it at `f64`, and
//!    `ppvm-traits-2` provides [`Coefficient`](ppvm_traits_2::Coefficient) only
//!    for `f64` / `Complex<f64>` — a generic-`R` tableau could not implement the
//!    `-2` gate traits at all. Amplitudes are therefore `Complex<f64>` and
//!    probabilities `f64`.
//! 4. **`Accumulate` on [`Amplitudes`].** The store implements the graded
//!    algebra's `Support` (L0), `Scale` (L2) and `Retain`. `Accumulate` (L1) is
//!    *not* implemented: the tableau's coalesce is the inline `u64`-packed
//!    sort-merge in `branch_with_coefficients`, and a generic
//!    `accumulate_batch` over a `Vec` store would have to fall back to the
//!    linear-scan `add_or_insert`, which is `O(m²)` on a 65 536-branch state —
//!    exactly the cost the sort-merge exists to avoid.
//! 5. **`Sum<Vec<(Bitstring, Complex<C>)>, CoefficientThreshold>`.** The Phase-4
//!    plan types the amplitudes as `ppvm-pauli-sum-2`'s `Sum`, but `Sum` bounds
//!    `S::Key: Word + Indexable` and a bitstring is neither (a point already
//!    acknowledged in `ppvm-traits-2::batch`'s friction note). [`Amplitudes`] is
//!    the same graded object with a crate-local `Vec` backend.
//! 6. **Suspected old-crate defects, reproduced pending sign-off.** Each is
//!    documented at its call site: `asymmetric_loss_channel` does not pop the
//!    record entry its internal `measure` pushed (unlike `loss_channel` and
//!    `reset`); the loss channels fire on `p >= r` while the depolarizing family
//!    fires on `p > r`; and the two-qubit Pauli channels return early on loss
//!    **without** drawing, breaking the single-qubit siblings' documented
//!    "preserve the seeded RNG stream across loss events" invariant. All three
//!    are behaviour-changing for seeded runs if corrected, so they are
//!    reproduced verbatim.
//! 7. **Lazy-digest cache representation.** The design sketches the frame's
//!    lazy hash as an `OnceLock<u64>`; [`Tableau`] uses the sentinel
//!    `AtomicU64` that `ppvm-pauli-word-2` already settled on. Same contract
//!    (lazy, interior-mutable, `Send + Sync`, same finalized digest value —
//!    which is all the contract fixes), measurably cheaper to *invalidate*,
//!    which every Clifford gate must do. See [`Tableau`]'s "Cache
//!    representation" note for the A/B.
//! 8. **`rotate_2` merge order — the one adjudicated divergence.** Old merged
//!    into a `std::collections::HashMap` seeded from process entropy, so the
//!    resulting amplitude *order* was randomized per process and no old order
//!    exists to preserve. This port merges deterministically (ascending by
//!    index, as the `T` path does). Support and per-index values are unchanged;
//!    see `rotate_2`'s doc comment.

/// Clifford row operations: the symplectic/phase primitives, the fused
/// `Clifford`/`CliffordExtensions` impls, the frame primitives, and the batched
/// mask sweeps.
pub mod clifford;
/// The tableau representation: rows, the frame, the amplitude store, and the
/// generalized tableau.
pub mod data;
/// `Display` rendering of the frame and the generalized tableau.
pub mod display;
/// Pauli-string expectation values.
pub mod expectation;
/// Non-Clifford gates: `T`, the rotations, `R` and `U3`.
pub mod gates;
/// The inverse tableau's sign algebra: one rule per Clifford, and the `O(1)`
/// decomposition phase that replaces a fold of `k` generators.
pub(crate) mod inverse;
#[cfg(test)]
mod inverse_tests;
/// Z-basis measurement, `measure_all`/`measure_many`, and the reusable scratch.
pub mod measure;
/// Classical probability mixtures over complete generalized-tableau states.
pub mod mixture;
/// Stochastic channels: Pauli errors, depolarizing, and loss.
pub mod noise;
/// The frame's physical storage: one aligned contiguous allocation holding four
/// square X/Z quadrants, two phase bit planes and the loss plane.
pub(crate) mod storage;

pub use bnum::types::{U256, U512, U1024, U2048};
pub use data::{Amplitudes, Bitstring, GeneralizedTableau, Tableau};
pub use measure::MeasureScratch;
pub use mixture::{GeneralizedTableauMixture, GeneralizedTableauSum, MixtureSampler};
pub use noise::TableauLike;

/// Convenience re-exports for downstream code.
pub mod prelude {
    pub use crate::data::{Amplitudes, Bitstring, GeneralizedTableau, Tableau};
    pub use crate::measure::MeasureScratch;
    pub use crate::mixture::{GeneralizedTableauMixture, GeneralizedTableauSum, MixtureSampler};
    pub use crate::noise::TableauLike;
    pub use bnum::types::{U256, U512, U1024, U2048};
    pub use ppvm_traits_2::prelude::*;
}
