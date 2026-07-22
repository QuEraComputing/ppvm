// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lattice translation symmetry groups for operator-space Pauli evolution.
//!
//! A [`TranslationGroup`] represents a finite abelian group `G` acting on
//! qubit positions by permutations. Given such a group, every Pauli word
//! belongs to a translation orbit, and operator dynamics that commute
//! with `G` can be tracked using **one canonical representative per
//! orbit** instead of all `|G|` orbit members — reducing per-step memory
//! and compute by a factor up to `|G|`.
//!
//! Following Teng, Chang, Rudolph, and Holmes (arXiv:2512.12094), this
//! module implements **plain (real-coefficient) merging** of Pauli sums
//! into orbit-representative form — see [`canonicalize_pauli_sum`] and
//! [`symmetry_merge_pauli_sum`]. This handles observables in the trivial
//! (`k=0`) symmetry sector, e.g. sums of single-Z operators over the
//! lattice.
//!
//! **Non-trivial momentum sectors (`k ≠ 0`)** are handled by
//! [`canonicalize_pauli_sum_complex`], which folds with the character
//! phase `χ_k(g)` of each translation. On the Python side, an operator in
//! sector `k` is carried as a *real pair* (real + imaginary components, two
//! real `PauliSum`s) and merged via `PauliSum.momentum_merge`, which reuses
//! this routine — letting gate-based Trotter evolution stay symmetry-
//! compressed in any momentum sector with real coefficients throughout.
//!
//! ## Data model
//!
//! A `TranslationGroup` is specified by a list of generator permutations
//! and their cyclic orders. The group order is the product of the orders.
//! For instance, a 2D `L × L` torus has two generators (translation in
//! x and y) each of order `L`.
//!
//! ## Canonicalization
//!
//! [`TranslationGroup::canonicalize`] returns the **lex-minimum** Pauli
//! word reachable from the input via group action. The ordering is the
//! standard `Ord` impl on `PauliWord` (compare `xbits`, then `zbits`).
//! All orbit members canonicalize to the same representative; orbits are
//! disjoint by construction, so the rep uniquely identifies the orbit.
//!
//! ## Merging
//!
//! [`canonicalize_pauli_sum`] takes parallel `Vec<Word>` / `Vec<f64>`
//! buffers (the representation used by ppvm-lindblad's adaptive
//! evolution) and replaces each Pauli by its canonical rep, summing
//! coefficients for collisions. The output is an orbit-rep basis with
//! coefficients equal to the sum of the input coefficients over each
//! orbit's members. For dynamics that commute with `G` and initial
//! states that are also `G`-invariant, this preserves the expectation
//! value of any `G`-invariant observable (Theorem 1 of arXiv:2512.12094).
//!
//! See the dedicated tests for correctness against full-basis evolution
//! on small systems with no truncation.

mod group;
mod merge;
mod momentum;

pub use group::TranslationGroup;
pub use merge::{canonicalize_pauli_sum, symmetry_merge_pauli_sum};
pub use momentum::{SectorCheckError, canonicalize_pauli_sum_complex, check_momentum_sector};

#[cfg(test)]
mod tests;
