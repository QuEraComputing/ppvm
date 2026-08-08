// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`ColumnStore`] — the **structure-of-arrays** storage backend for
//! [`Sum`](crate::Sum), a sibling of [`HashMapStore`](crate::HashMapStore).
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Backends are containers;
//! columnar is expressible from day one":
//!
//! > **`ColumnStore`** — the one backend that *must* be a new struct, because it
//! > is structure-of-arrays: coefficients in one contiguous column, keys in plane
//! > blocks. Same hash-join build as the `HashMap`, but `scale` is one vectorized
//! > `*=`, `reduce` is a prefix-sum compaction, and `probe_batch` uses coalesced
//! > gathers. It is `HashMap`'s data re-laid for SIMD / GPU, not a third
//! > collection. Requires `K: Indexable + Columnar`.
//!
//! and `word-data-structures.md` §"Key columns (structure-of-arrays batches)"
//! (the column is "two plane blocks and a shared width"), whose
//! [`PauliKeyColumn`](ppvm_pauli_word_2::PauliKeyColumn) this backend **consumes**
//! rather than reinventing (implementation-plan Phase 6: "Requires `K: Columnar`
//! (already impl'd in Phase 2)").
//!
//! # Layout
//!
//! ```text
//! keys    : K::Column          the SoA plane blocks (X plane, Z plane, width)
//! coeffs  : Vec<C>             the contiguous coefficient column
//! hashes  : Vec<u64>           the parallel finalized-digest column
//! live    : Vec<u8>            the parallel physical-row liveness column
//! index   : HashTable<(u32,u64)>  open-addressed digest → physical row
//! ```
//!
//! The four data columns are **parallel**: row `i` is the term
//! `(keys[i], coeffs[i])` with digest `hashes[i]`. The keys never live in the
//! index — it holds `u32` slot numbers only — so a probe walks a `u32` table and
//! confirms against the key planes, and a *coefficient* pass
//! ([`Scale::scale`](ppvm_traits_2::Scale), [`Accumulate::reduce`]) touches one
//! contiguous `Vec<C>` with **no key-sized stride at all**. That is the whole
//! point of the layout, and it is the property the hash map cannot have.
//!
//! Retain operations first mark rejected rows dead without moving the physical
//! columns. The index may temporarily point at such a tombstone; public probes
//! filter it, while re-insertion uses it to repoint the existing bucket to a new
//! row appended at the physical tail. Once tombstones reach one eighth of the
//! physical rows (rounded up), a stable compaction preserves live row order and
//! rebuilds the index.
//!
//! # Behaviour is identical to [`HashMapStore`](crate::HashMapStore)
//!
//! A backend swap is observationally a no-op (implementation-plan Phase 6: "a
//! backend swap must be observationally identical"). Every contract the hash-join
//! backend documents holds here verbatim and for the same reasons:
//!
//! * **no gate truncates** — the policy runs only from
//!   [`Sum::truncate`](crate::Sum::truncate);
//! * **no operation runs `reduce` implicitly**, and **no operation drops a zero
//!   coefficient**: a term scaled to exactly `0.0`, an exactly-cancelling
//!   accumulation, and an exactly-zero rotation branch all stay in the support
//!   and count in [`Support::len`] (old has no `reduce` at all and its exact-map
//!   `PartialEq` depends on zeros surviving —
//!   `ppvm-pauli-sum/tests/loss.rs::test_reset_channel`);
//! * [`RotateInPlace`] keeps the **two-pass** ordering (scale *all* diagonals,
//!   then merge *all* branches), which is not an optimization —
//!   `eagerWalk_ne_twoPass` in `lean/PPVM/Instantiations/Rotation.lean` exhibits
//!   the divergence an interleaved columnar walk would produce;
//! * [`ApplyProducer`] **replaces** rather than merges (`reset` between produce
//!   and accumulate), per `pushforward_eq_reset_accumulate` in
//!   `lean/PPVM/Algebra/GradedMap.lean`.
//!
//! Iteration order differs from the hash map's (this backend iterates in
//! *insertion* order), which is unobservable: both orders are unspecified, and
//! the only order-sensitive numerics are float summations, whose reassociation
//! the differential bars already treat as relative-tolerance.
//!
//! # The double-buffer is preserved
//!
//! Architecture feature 1 (the persistent `(primary, aux)` ping-pong) and
//! feature 2 (the reusable branch `scratch`) are carried here exactly as in
//! [`HashMapStore`](crate::HashMapStore): `aux` is a second full column set,
//! allocated once and never freed, used by [`MultiplyInPlace`]; `scratch` and
//! `batch` are the rotation/producer buffers. [`RekeyBijective`] needs *neither*
//! — see its impl, which rewrites the key planes **in place** and never moves a
//! coefficient at all.
//!
//! # Measured performance
//!
//! The final Phase-6 gate runs both backends in one binary with identical
//! `PauliWord<[u8; 8]>`, `f64`, policy and capacity. Four fresh-process medians:
//!
//! | target | column/hash ratio |
//! |---|---|
//! | steady paired `rx` | **0.94–0.97** |
//! | sustained `rx` growth | **0.93–0.96** |
//! | `from_terms` | **0.89–0.92** |
//! | active 408→390 truncation | **0.94–0.96** |
//! | realistic-support `rx` | **0.95–0.98** |
//! | realistic-support native `rzz` | **0.93–0.95** |
//! | noisy TFIM Trotter, n=12, 10 steps | **0.98–0.99** |
//!
//! The former insertion regressions were implementation gaps, not necessary SoA
//! trade-offs. The fixes are structural: a grouped row index, packed-bit
//! rotation kernels with dense and sparse monomorphs, paired closed-support
//! updates, lazy auxiliary allocation, and tombstoned retention with stable
//! compaction at 1/8 dead rows. Sparse rotation traverses lazily cached live-row
//! runs; dense coefficient passes remain contiguous and branch-free.
//!
//! Bijective Clifford re-key still rewrites key planes in place, while the hash
//! backend must drain and reinsert. Exact-zero and insertion-order semantics are
//! unchanged: liveness is independent of coefficient value, removed keys are
//! absent immediately, and a removed key re-added before compaction appends at
//! the end by repointing its stale index entry.
//!
//! # Friction: the whole-map capability traits keep SoA *expressible*, but their `&K` callbacks tax it
//!
//! The design's stated reason for making [`SignFlipByKey`], [`ScaleByKey`] and
//! [`RotateInPlace`] whole-map capabilities rather than per-slot callbacks was
//! "so a columnar backend stays expressible". That claim **holds**: every one of
//! them is implemented here with no change to its signature, and no signature
//! anywhere forced an array-of-structs layout (design rule 4). But *expressible*
//! is not *free*: each of those callbacks — and [`Support::for_each_ref`], and
//! [`Retain::retain`] — is handed a `&K`, which a structure-of-arrays backend
//! can only produce by **materializing** the key from its planes. That
//! materialization still distinguishes generic callback paths from the
//! column-native rotation and coefficient-only fast paths. A future generic
//! cursor shape could remove it without exposing an array-of-structs layout.

use ppvm_traits_2::{
    Accumulate, Coefficient, Columnar, Conjugate, ImaginaryUnit, KeyBatch, KeyColumn, KeyProduct,
    Multiply, Pair, PauliBits, Retain, Scale, Support, TermBatch, TermProducer,
};

use crate::store::{
    AddTerm, ApplyProducer, BranchInPlace, InsertTerm, MultiplyInPlace, RekeyBijective,
    RotateInPlace, ScaleByKey, SignFlipByKey, StoreAlloc, product_capacity_hint,
};

mod columns;
mod gates;
mod graded;
mod lifecycle;
mod rotations;

pub use lifecycle::ColumnStore;

#[cfg(test)]
mod tests;
