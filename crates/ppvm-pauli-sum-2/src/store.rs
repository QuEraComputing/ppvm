// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The storage-lifecycle helper the [`Sum`](crate::Sum) engine needs *beyond*
//! the graded algebra: allocate a fresh empty store at a capacity hint, and
//! reset a store to empty support **keeping its allocation**.
//!
//! # Friction: `apply` needs a clear/replace the graded traits do not expose
//!
//! The design's `apply` sketch (`traits-2-configuration-and-hashing.md`
//! §"apply") reads:
//!
//! ```text
//! for (k, c) in self.storage.iter() { producer.produce(&k, &c, &mut batch); }
//! self.storage.accumulate_batch(&batch);   // <-- merges onto the *existing* support
//! self.storage.reduce();
//! ```
//!
//! For a **bijective Clifford re-key** the producer emits `(φ(k), c)` for every
//! stored `(k, c)`; merging that batch onto the un-cleared support leaves both
//! the old `k` *and* the new `φ(k)`, double-counting. The old
//! `ppvm-pauli-sum::map_add` avoided this by writing into a **cleared** auxiliary
//! map and swapping. The producer path must reproduce that "replace, not merge"
//! semantics, so [`ApplyProducer`] resets the support between producing and
//! accumulating. (It also owns the produced-term `batch`, which is per-call state
//! `Sum` has nowhere to keep — see that trait's friction note.)
//!
//! [`Accumulate`](ppvm_traits_2::Accumulate) exposes no `clear`, so this trait
//! supplies one. It is a **local** trait `impl`'d on the foreign `Vec`/`HashMap`
//! containers (legal: the trait is local), mirroring the graded-algebra impls in
//! `ppvm-traits-2`. `reset` clears in place so the allocation is reused across
//! millions of gates — matching the old aux-map reuse and avoiding a per-gate
//! reallocation.

use std::collections::HashMap;

use ppvm_traits_2::{
    Accumulate, Coefficient, Conjugate, IdentityBuildHasher, ImaginaryUnit, Indexable, KeyBatch,
    KeyProduct, Multiply, Pair, Retain, Scale, Support, TermBatch, TermProducer,
};

/// The provided hash-join storage backend for [`Sum`](crate::Sum): a primary
/// support map plus a persistent **auxiliary** map and **scratch** buffer, all
/// keyed through the pass-through [`IdentityBuildHasher`] (so a key's finalized
/// `key_hash()` digest reaches hashbrown untouched — Design: §"The pass-through
/// storage contract").
///
/// This relocates the old `ppvm-pauli-sum` double-buffer (`PauliSum` held
/// `map: (primary, aux)` + a reusable `scratch` Vec) **into the storage backend**,
/// so `Sum` stays a pure engine over an abstract store — and the same backend
/// composes into the generalized tableau — while the hot fast paths still reuse
/// their buffers across millions of gates instead of reallocating per gate:
///
/// * [`RekeyBijective`] (Clifford re-key) drains the primary through the re-key
///   closure into the cleared `aux`, then swaps — the old `map_add`'s
///   clear→write→swap: zero per-gate allocation, zero key clones. The two map
///   allocations ping-pong across gates and are never freed.
/// * [`RotateInPlace`] (rotation branch) buffers the ≤`N` branch terms in the
///   persistent `scratch` before merging — the old `map_insert`'s reused scratch.
/// * [`ApplyProducer`] (the generic [`TermProducer`] path behind
///   [`Sum::apply`](crate::Sum::apply)) stages the produced terms in the
///   persistent `batch` — the design's "driver-owned reusable batch", owned by
///   the driver that has somewhere to keep it.
///
/// `aux`, `scratch` and `batch` are **transient**: all three are empty between
/// operations (cleared on entry; left empty by the drain/swap on exit), so
/// [`Clone`] copies only the primary support and equality/iteration observe only
/// it.
#[derive(Debug)]
pub struct HashMapStore<K, C> {
    /// The reduced-canonical support the graded traits read and write.
    primary: HashMap<K, C, IdentityBuildHasher>,
    /// Persistent double-buffer for [`RekeyBijective`] (clear → write → swap).
    /// Empty between operations.
    aux: HashMap<K, C, IdentityBuildHasher>,
    /// Persistent branch/term scratch for [`RotateInPlace`]. Empty between
    /// operations.
    scratch: Vec<(K, C)>,
    /// Persistent produced-term batch for [`ApplyProducer`]. Empty between
    /// operations.
    batch: TermBatch<K, C>,
}

impl<K: Clone, C: Clone> Clone for HashMapStore<K, C> {
    /// Clone only the primary support; `aux`/`scratch`/`batch` are transient
    /// (empty between ops), so a clone starts them fresh rather than copying
    /// workspace.
    #[inline]
    fn clone(&self) -> Self {
        Self {
            primary: self.primary.clone(),
            aux: HashMap::with_hasher(IdentityBuildHasher),
            scratch: Vec::new(),
            batch: TermBatch::new(),
        }
    }
}

/// Allocate and reset a graded-sum storage container. Not a graded-algebra
/// operation — see the module docs on why `apply` needs the reset the graded
/// traits omit.
pub trait StoreAlloc {
    /// A fresh, empty store pre-sized for roughly `cap` terms.
    fn with_capacity(cap: usize) -> Self;

    /// Reset to empty support, keeping the backing allocation for reuse.
    fn reset(&mut self);
}

/// An in-place, **move-based** bijective re-key — the fast path a Clifford
/// conjugation drives instead of the [`apply`](crate::Sum::apply) batch
/// round-trip.
///
/// # Friction: `apply`'s batch clones every key twice; a Clifford cannot afford it
///
/// The generic `apply` reads the support through [`Support::iter`], which
/// **clones** every `(k, c)`, produces the re-keyed terms into a `TermBatch`,
/// then merges the batch with [`Accumulate::accumulate_batch`], which clones
/// each key a **second** time (`entry(k.clone())`). For a Pauli-keyed sum the new
/// word dropped `Copy` (its lazy `OnceLock` hash cache), so each clone is a real
/// cost; against the old crate's single in-place aux-map swap
/// (`ppvm-pauli-sum::PauliSum::map_add`, one clone per term) the batch path
/// benchmarked ~4.7× slower on a Clifford gate.
///
/// A Clifford conjugation is a **bijection** on Pauli words, so the transformed
/// support is just the old support with each key relabelled and its coefficient
/// scaled by `±1`: no new keys collide and no coefficient is zeroed. That lets us
/// skip the batch entirely and rebuild the support by **moving** each term
/// through the re-key closure — zero key clones, zero coefficient clones — while
/// reusing the backing allocation. This is a `ppvm-pauli-sum-2`-local trait
/// `impl`'d on the foreign `Vec`/`HashMap` containers (legal: the trait is
/// local), mirroring [`StoreAlloc`] and the graded-algebra impls in
/// `ppvm-traits-2`.
///
/// Design: `traits-2-configuration-and-hashing.md` §"apply" (the "driver-owned
/// reusable batch" / in-place-swap note) and §"Pauli algebra traits" (the sum
/// "applies the one-row action pointwise and drains each term's phase delta to
/// its coefficient"). The bijectivity this relies on is machine-checked in
/// `lean/PPVM/Pauli/Symplectic.lean` (`*_bijective`).
pub trait RekeyBijective<K, C> {
    /// Re-key every term by moving it through `f` and re-inserting the result,
    /// reusing the backing allocation. `f` maps `(k, c) ↦ (φ(k), c')`.
    ///
    /// The caller guarantees `f` is **injective on keys** (a Clifford
    /// conjugation is a Pauli bijection, machine-checked as `*_bijective` in
    /// `lean/PPVM/Pauli/Symplectic.lean`), so re-keyed terms never collide.
    ///
    /// Implementations may therefore insert without aggregating. That is a real
    /// precondition, not a hint: violating it **drops** a term rather than summing
    /// it. The `HashMapStore` impl `debug_assert!`s that no insert displaces an
    /// existing entry, so a caller-side bug trips loudly in debug builds. The
    /// defensive `and_modify`/`or_insert` this replaced cost ~1.3× on the whole
    /// re-key path (it compiled out of line), which is too much to pay for
    /// guarding a case the type-level contract already excludes.
    fn rekey_bijective<F>(&mut self, f: F)
    where
        F: FnMut(K, C) -> (K, C);
}

/// An in-place, **whole-map** per-key sign flip — the fast path a *pure-sign*
/// Clifford conjugation (`X`/`Y`/`Z`) drives instead of the move-based
/// [`RekeyBijective`] re-key.
///
/// # Friction: `X`/`Y`/`Z` leave the word fixed, so a re-key is pure waste
///
/// A single-qubit `X`/`Y`/`Z` conjugation is a **pure sign**: `XPX = (−1)^z P`,
/// `YPY = (−1)^{x⊕z} P`, `ZPZ = (−1)^x P` — the Pauli **word is unchanged**, only
/// each term's coefficient picks up a `±1` read off that term's own X/Z bits. The
/// generic [`RekeyBijective`] path still rebuilds the whole `HashMap` (fresh
/// alloc, move every entry into a new bucket) even though every key maps to
/// itself; the old crate (`ppvm-pauli-sum::PauliSum::x`/`y`/`z`) instead scaled
/// each existing entry's coefficient in place (`scale` with a conditional
/// `*= -1.0`), touching **zero** allocations and never moving a key.
///
/// This trait restores that zero-allocation fast path: walk the stored entries in
/// place and multiply each coefficient by `f(&k) ∈ {−1, +1}`, computed from that
/// key's bits. It is deliberately a **whole-map** capability rather than a
/// per-slot `&mut (K, C)` / `&mut [C]` callback so a columnar backend (bits in one
/// column, coefficients in another) stays expressible — the impl owns the walk and
/// may read the key column and scale the coefficient column independently. This is
/// a `ppvm-pauli-sum-2`-local trait `impl`'d on the foreign `Vec`/`HashMap`
/// containers (legal: the trait is local), mirroring [`RekeyBijective`] and
/// [`StoreAlloc`].
///
/// Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits" (the
/// sum "applies the one-row action pointwise and drains each term's phase delta to
/// its coefficient"). The pure-sign conjugation signs are machine-checked in
/// `lean/PPVM/Pauli/Conjugation.lean` and the phased word's
/// `PhaseTrack::x_phase`/`y_phase`/`z_phase`
/// (`ppvm-phased-pauli-word-2/src/clifford.rs`); the word being fixed is the
/// identity symplectic bit map of `lean/PPVM/Pauli/Symplectic.lean`.
pub trait SignFlipByKey<K, C> {
    /// Multiply every term's coefficient in place by `f(&k) ∈ {−1, +1}`, keyed on
    /// that term's own key — no key movement, no reallocation.
    fn sign_flip_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> i8;
}

/// An in-place, **whole-map** per-key coefficient scale by an arbitrary ring
/// factor — the fast path a *diagonal* unital Pauli channel ([`PauliError`]) drives
/// instead of a batch rebuild.
///
/// # Friction: noise scales by a real λ_P, not the `±1` [`SignFlipByKey`] allows
///
/// A unital single-qubit Pauli channel acts **diagonally** in the Pauli basis:
/// each term `P` is multiplied by a real transfer eigenvalue
/// `λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q` that depends only on `P`'s Pauli at
/// the target qubit — no key moves, no branch, no rebuild (the old
/// `ppvm-pauli-sum::PauliSum::pauli_error`'s in-place `scale`). [`SignFlipByKey`]
/// is the `{−1, +1}` special case of exactly this walk, but a noise factor is a
/// general `C`, so it cannot go through `mul_sign(i8)`; this trait generalizes the
/// per-key walk to a caller-supplied in-place mutation.
///
/// The callback mutates `&mut C` directly rather than returning a factor. That is
/// deliberately the old crate's `PauliSum::scale` shape (`Fn(&W, &mut C)`), and it
/// matters for two reasons. **Behaviour:** a walk that can only *mutate* cannot
/// remove, so a term scaled to exactly zero stays in the support — which is what
/// old does (old has no `reduce` at all, and its exact-map `PartialEq` depends on
/// zero terms surviving; see `ppvm-pauli-sum/tests/loss.rs::test_reset_channel`).
/// An earlier `Option<C>`-returning shape had to run `retain` to drop zeros, which
/// diverged from old. **Performance:** `retain` drags in hashbrown's erase
/// machinery for every slot; the plain `iter_mut` walk is ~2× faster on the noise
/// sweep and matches old. A no-op slot (identity at the qubit, or a lost qubit) is
/// expressed by the callback simply not touching `c`.
///
/// Like its siblings it is deliberately a **whole-map** capability rather than a
/// per-slot `&mut (K, C)` callback so a columnar backend stays expressible, and it
/// is a `ppvm-pauli-sum-2`-local trait `impl`'d on the foreign `Vec`/`HashMap`
/// containers (legal: the trait is local).
///
/// Design: `traits-2-configuration-and-hashing.md` §"Behavioral traits"
/// (`PauliError`) and §"Pauli algebra traits". The eigenvalue this scales by is
/// machine-checked in `lean/PPVM/Algebra/Noise.lean` (`pauli_channel_eigenvalue`,
/// `pauli_channel_eigenvalue_omega`).
pub trait ScaleByKey<K, C> {
    /// Apply `f(&k, &mut c)` to every term in place, keyed on that term's own
    /// key — no key movement, no reallocation, and **no removal**: a coefficient
    /// driven to exactly zero stays in the support, as in old. A no-op slot is
    /// expressed by leaving `c` untouched.
    fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K, &mut C);
}

/// An in-place, **fused branching** re-key — the fast path a single-qubit
/// rotation drives instead of the general [`apply`](crate::Sum::apply) batch
/// round-trip.
///
/// # Friction: `apply` re-hashes the whole 2N-term fan-out; a rotation's diagonal never moves
///
/// A single-qubit rotation `exp(−i·θ/2·G)` conjugates each term `(P, c)` to
/// `(P, c·cosθ)` **plus**, when `G` anticommutes with `P`, a genuinely-new branch
/// `(iGP, c·sinθ·ε)`. The **diagonal** `P` keeps its key — only its coefficient
/// is scaled by `cosθ` — so its cached structural hash stays valid and it should
/// never leave its bucket. Routing a rotation through `apply`, though, clones
/// every key on the read side ([`Support::iter`]), `reset`s the whole map, then
/// re-accumulates **all** ~`2N` produced terms — re-hashing and re-inserting even
/// the `N` untouched diagonals. Against the old crate's fused single pass
/// (`ppvm-pauli-sum::sum::rot1` over `PauliSum::map_insert`), which mutates each
/// diagonal coefficient in place and hashes only the ≤`N` branch terms, the batch
/// path benchmarked ~2.4× slower on an `rx` over a ~1000-term sum.
///
/// This trait restores that fused single pass: walk the stored entries in place,
/// letting `f(&k, &mut c)` scale each diagonal coefficient **in place** (no key
/// clone, no re-hash) and return the optional branch term; buffer the branch
/// terms (cheap, un-hashed pushes) and merge only those back through the map's
/// hash-join, aggregating any collision (a branch may land on another term). It
/// is deliberately a **whole-map** capability rather than a per-slot callback so a
/// columnar backend stays expressible, and it is a `ppvm-pauli-sum-2`-local trait
/// `impl`'d on the foreign `Vec`/`HashMap` containers (legal: the trait is local),
/// mirroring [`SignFlipByKey`], [`ScaleByKey`], and [`RekeyBijective`].
///
/// **Every** produced branch term is merged, including one whose coefficient is
/// exactly zero (an identity rotation `θ = 0` has `sinθ = 0`, so `R₀` on a
/// `Z`-bearing term produces a `0.0`-coefficient `Y` branch). That is old's
/// behaviour and it is user-visible: old's `map_insert` unconditionally
/// `add_assign`s every term its first pass produced
/// (`ppvm-pauli-sum/src/sum/data.rs`), and `add_assign` is
/// `entry().and_modify().or_insert(v)` — it inserts the zero. Old's exact-map
/// `PartialEq` counts that key, so skipping it would change `len()`, `get()`,
/// the key set and equality. No whole-map `reduce` scan runs either: a generic
/// rotation's *collision* cancellations are left in place — like the old crate's
/// `map_insert`, which leaves any residue for a later truncation, the fast path
/// skips the retain [`apply`](crate::Sum::apply) would run. A physical near-zero
/// is dropped by the policy's truncation the caller runs afterward.
///
/// Design: `traits-2-configuration-and-hashing.md` §"apply" and §"Every gate is a
/// producer feeding `accumulate`" (the fused in-place alternative to the batch
/// round-trip), and §"Behavioral traits" (`RotationOne`). The branch
/// `c·P → cos·c·P + sin·c·(iGP)` — a genuinely-new key, norm-preserving and
/// angle-additive — is machine-checked in `lean/PPVM/Instantiations/Rotation.lean`
/// (`anticommute_new_key`, `rot_norm_sq`, `rot_rot`).
pub trait RotateInPlace<K, C> {
    /// Walk every term `(k, c)` in place: `f(&k, &mut c)` scales the diagonal
    /// coefficient in place (the `cosθ` factor) and returns `Some((k′, c′))` for
    /// the branch term to merge, or `None` when the term commutes. Colliding
    /// branch keys are aggregated; an exactly-zero branch is merged like any
    /// other (old inserts it — see the trait docs) and no whole-map `reduce` scan
    /// runs (collision cancellations are left for the policy's truncation).
    ///
    /// **Every implementer owes the two-pass ordering**: *all* diagonals must be
    /// scaled before *any* branch is merged. It is not an optimization — a branch
    /// produced from `k` can land on a different key `k′` that has not been scaled
    /// yet (`rx` on a support holding `Z` and `Y` at one site sends `Z ↦ Y` and
    /// `Y ↦ Z` at once), so an implementation that merged eagerly inside the walk,
    /// or a parallel/columnar backend that interleaved the passes, would compute a
    /// different map. Lean: `accumulate_rotBatch` in
    /// `lean/PPVM/Instantiations/Rotation.lean` (the two-pass result equals the
    /// one-pass produced batch, for every walk order), with `eagerWalk_ne_twoPass`
    /// exhibiting the interleaved variant's divergence.
    fn rotate_in_place<F>(&mut self, f: F)
    where
        F: FnMut(&K, &mut C) -> Option<(K, C)>;
}

/// An in-place **L4 operator product** `A ← A · B` that rebuilds the support
/// through the store's own double-buffer instead of allocating an accumulator.
///
/// # Friction: `Multiply::multiply_into` needs a *third* map, and a fresh one would allocate
///
/// [`Multiply::multiply_into`](ppvm_traits_2::Multiply) is deliberately
/// three-address (`acc += self · other`): the twisted convolution's outer product
/// writes each `(p, q)` contribution to key `p·q`, which may collide with a key
/// still to be read from either operand, so it *cannot* be folded back into an
/// operand in place. Old's `MulAssign<PauliSum>` tried to and got the wrong
/// algebra (it computed the product **chain** `A·b₀P₀·b₁P₁` instead of the
/// bilinear sum — `ppvm-pauli-sum/src/sum/ops.rs:70`; see `crate::multiply`).
///
/// A correct `A *= B` therefore needs a scratch map — and the store already owns
/// one. This trait routes the accumulator to the persistent `aux` of the
/// double-buffer (architecture feature 1: "any NEW support-rebuilding operation
/// — notably `Multiply`'s accumulator — must draw from the same aux rather than
/// allocating"), so the in-place product is the old `map_add` shape exactly:
/// clear `aux` → write the whole convolution into it → swap. The two map
/// allocations keep ping-ponging and are never freed.
///
/// Like its siblings this is a `ppvm-pauli-sum-2`-local trait `impl`'d on the
/// foreign `Vec`/`HashMap` containers (legal: the trait is local).
/// Accumulate **one** term into the support — the storage capability behind the
/// engine's `sum += (key, coeff)` operator.
///
/// # Friction: `Accumulate` is batch-only, and a single `+=` should not allocate
///
/// [`Accumulate::accumulate_batch`](ppvm_traits_2::Accumulate) is the graded
/// algebra's merge, and it takes a whole [`TermBatch`]. Routing a one-term `+=`
/// through it means building a batch (two `Vec`s) per call, where the old crate's
/// `sum += (word, coeff)` is a single `entry().and_modify().or_insert()` on the
/// map (`ACMapAddAssign::add_assign`, `ppvm-traits/src/map/hashmap.rs`). Seeding
/// the headline Trotter observable `Σᵢ Zᵢ` is `n` such calls, and a caller may
/// stream terms one at a time, so the allocation is not amortized away.
///
/// This is the single-term door onto the same merge: **accumulate**, never
/// replace, and — like old's `add_assign` — it inserts a zero coefficient rather
/// than dropping it.
pub trait AddTerm<K, C> {
    /// Add `coeff` onto the coefficient at `key`, inserting the term if absent.
    fn add_term(&mut self, key: K, coeff: C);
}

pub trait MultiplyInPlace<K, C> {
    /// Replace this support with the twisted convolution `self · other`,
    /// reusing the backing allocations. Accumulates colliding products; runs no
    /// `reduce` and no truncation (an exact-zero cancellation survives).
    fn multiply_in_place(&mut self, other: &Self);
}

/// The generic [`TermProducer`] step — "produce the whole support through the
/// gate, then **replace** the support with the accumulated batch" — with the
/// produced-term batch owned by the **store**.
///
/// # Friction: the batch is per-call state, and `Sum` has nowhere to keep it
///
/// The design's `apply` sketch (§"apply") stages the produced terms in a
/// `TermBatch` and notes it should be a "driver-owned reusable batch". `Sum` is
/// deliberately pure data (storage + policy + width — §"There is no `SumStorage`
/// trait, and no owned workspace"), so the driver that can own it is the storage
/// backend, which already owns the `aux` double-buffer and the rotation
/// `scratch`. Building the batch inside `Sum::apply` instead means two `Vec`
/// allocations **per gate** — precisely the per-gate allocation churn
/// architecture features 1–2 exist to remove (measured at ~+13% end-to-end for
/// the aux map, i.e. ~6× what the single-gate microbench reported).
///
/// It is latent today only because all four hot gate families bypass `apply`
/// ([`RekeyBijective`], [`RotateInPlace`], [`ScaleByKey`], [`SignFlipByKey`]);
/// the first producer-based gate that lands (the multiply producer, the lossy
/// branching channels, a tableau-keyed gate) would otherwise reintroduce it with
/// no bench to catch it. Putting the batch here means it cannot come back.
///
/// Like its siblings this is a `ppvm-pauli-sum-2`-local trait `impl`'d on the
/// foreign `Vec`/`HashMap` containers (legal: the trait is local).
pub trait ApplyProducer<K, C> {
    /// Produce every stored term through `producer` and **replace** the support
    /// with the accumulated result.
    ///
    /// "Replace, not merge" is load-bearing: for a bijective re-key the producer
    /// emits `(φ(k), c)` for every `(k, c)`, and merging onto the un-cleared
    /// support would leave both (see the module docs). The support is therefore
    /// reset between producing and accumulating.
    ///
    /// Lean: `pushforward_eq_reset_accumulate` in
    /// `lean/PPVM/Algebra/GradedMap.lean` — reset-then-accumulate *is* the
    /// pushforward `mapDomain φ ∘ mapRange g` (`pushforward_apply` gives the
    /// per-key form for an injective `φ`, i.e. a Clifford's symplectic
    /// bijection) — and `merge_without_reset_ne_pushforward`, which exhibits the
    /// double-count that dropping the reset produces. The design's `apply`
    /// sketch omits this step; that omission is the `apply`-path analogue of
    /// `eagerWalk_ne_twoPass`.
    ///
    /// Runs **no** `reduce` and **no** truncation — a produced term whose
    /// coefficient is exactly zero survives, and two sub-threshold contributions
    /// to one key merge rather than being dropped. Both are caller-driven
    /// ([`Sum::reduce`](crate::Sum::reduce), [`Sum::truncate`](crate::Sum::truncate)).
    fn apply_producer<TP>(&mut self, producer: TP)
    where
        TP: TermProducer<K, C>;
}

/// The ceiling on the eagerly-reserved product accumulator: `2²⁰` entries
/// (~32 MB at a `[u8; 8]` key + `Complex<f64>` coefficient). Past this the
/// convolution falls back to hashbrown's doubling, trading a few rehashes for a
/// bounded up-front commitment.
const PRODUCT_CAPACITY_CEILING: usize = 1 << 20;

/// Capacity hint for the accumulator of a twisted convolution `A · B` whose
/// operands have `a_len` and `b_len` terms.
///
/// The convolution emits `|A|·|B|` contributions, landing on **up to** `|A|·|B|`
/// distinct keys, so that product — not `max(|A|, |B|)` — is the right size for
/// the accumulator. `max(|A|, |B|)` is only the *lower* bound (the key product is
/// a bijection in each argument, so the support is at least as large as either
/// operand), and starting there makes the accumulator pay a chain of
/// doubling + full-rehash rounds over a growing map, which is exactly the
/// mid-run rehash stall architecture feature 6 exists to prevent.
///
/// Measured, one variable held (accumulator initial capacity), same process,
/// 256 random 12-qubit terms squared, `|A·A| = 32611`:
///
/// | initial capacity | time |
/// |---|---|
/// | `\|A\|·\|B\|` (65536) | 444 µs |
/// | 16384 (the `NoPolicy` hint a fresh accumulator used to get) | 582 µs |
/// | 256 (`max(\|A\|, \|B\|)`, the previous heuristic) | 778 µs |
///
/// i.e. the old heuristic cost **1.75×** against a correctly-sized accumulator,
/// and is what made the buffer-reusing `multiply_in_place` *slower* than the
/// allocating product. The [`PRODUCT_CAPACITY_CEILING`] keeps the eager
/// reservation bounded for very large operands, where the doubling chain is the
/// cheaper end of the trade.
///
/// # Why `7/8` of the pair count and not the pair count itself
///
/// hashbrown sizes its table to `next_pow2(⌈requested · 8/7⌉)` buckets, so
/// asking for `|A|·|B|` *slots* rounds the **table** up to twice the next power
/// of two above `|A|·|B|` — and iteration, `clear` and drop are all `O(buckets)`,
/// not `O(len)`. Asking for `7/8` of the pair count instead lands the table
/// exactly on `next_pow2(|A|·|B|)`, which still holds `|A|·|B|` entries at
/// hashbrown's `7/8` load factor. Measured on the `multiply_sum` bench (256
/// 12-qubit terms, `|A·A| = 32611`), the difference is real in both directions:
/// requesting the full pair count made the *fresh-accumulator* product 14%
/// **slower** (434 µs → 493 µs — a 131072-bucket table walked by the subsequent
/// `overlap` and drop), while the `7/8` form keeps the in-place product's win and
/// leaves the fresh one at parity. The floor keeps the lower bound
/// (`max(|A|, |B|)`) for tiny operands, where `7/8` of a handful rounds to
/// nothing.
#[inline]
pub(crate) fn product_capacity_hint(a_len: usize, b_len: usize) -> usize {
    let pairs = a_len.saturating_mul(b_len).min(PRODUCT_CAPACITY_CEILING);
    (pairs - pairs / 8).max(a_len.max(b_len))
}

impl<K, C> ScaleByKey<K, C> for Vec<(K, C)>
where
    C: Coefficient,
{
    #[inline]
    fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K, &mut C),
    {
        // A pure in-place walk — the old crate's `scale`. Nothing is removed, so a
        // term scaled to exactly zero survives (old keeps it).
        for (k, c) in self.iter_mut() {
            f(k, c);
        }
    }
}

impl<K, C> ScaleByKey<K, C> for HashMap<K, C, IdentityBuildHasher>
where
    C: Coefficient,
{
    #[inline]
    fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K, &mut C),
    {
        // Walk the existing buckets in place: read each key, let `f` mutate its
        // coefficient. No bucket is moved and the backing allocation is untouched
        // — the old crate's in-place `scale`, restored. Using `iter_mut` rather
        // than `retain` keeps hashbrown's erase machinery out of the hot walk AND
        // preserves old's semantics that a zero coefficient is never removed.
        for (k, c) in self.iter_mut() {
            f(k, c);
        }
    }
}

impl<K, C> RotateInPlace<K, C> for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    #[inline]
    fn rotate_in_place<F>(&mut self, mut f: F)
    where
        F: FnMut(&K, &mut C) -> Option<(K, C)>,
    {
        // Pass 1: scale each diagonal coefficient in place, buffering the
        // anticommuting branch terms. The diagonal keeps its slot — no clone.
        let mut branches: Vec<(K, C)> = Vec::with_capacity(self.len());
        for (k, c) in self.iter_mut() {
            if let Some(term) = f(k, c) {
                branches.push(term);
            }
        }
        // Pass 2: merge only the branch terms back (linear-scan hash-join, the
        // small-support cost model this backend targets), aggregating collisions.
        // Every produced branch is merged, *including* an exactly-zero one (the
        // identity rotation `θ = 0` has `sinθ = 0`): old's `map_insert` merges
        // unconditionally through `add_assign`, which inserts the zero, and its
        // exact-map equality counts that key. No whole-map `reduce` scan runs —
        // a generic rotation's collision cancellations are left for the policy's
        // truncation, matching the old crate's `map_insert`.
        // The merge is [`AddTerm`] — the single definition of old's `add_assign`.
        for (nk, nc) in branches {
            AddTerm::add_term(self, nk, nc);
        }
    }
}

impl<K, C> SignFlipByKey<K, C> for Vec<(K, C)>
where
    C: Coefficient,
{
    #[inline]
    fn sign_flip_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> i8,
    {
        for (k, c) in self.iter_mut() {
            *c = c.mul_sign(f(k));
        }
    }
}

impl<K, C> SignFlipByKey<K, C> for HashMap<K, C, IdentityBuildHasher>
where
    C: Coefficient,
{
    #[inline]
    fn sign_flip_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> i8,
    {
        // Walk the existing buckets in place: read each key, scale its coefficient
        // by the `±1` that key's bits demand. No bucket is moved and the backing
        // allocation is untouched — the old crate's in-place `scale`, restored.
        for (k, c) in self.iter_mut() {
            *c = c.mul_sign(f(k));
        }
    }
}

impl<K, C> RekeyBijective<K, C> for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    #[inline]
    fn rekey_bijective<F>(&mut self, mut f: F)
    where
        F: FnMut(K, C) -> (K, C),
    {
        // Coordinate-list backend (small support): move every term out with
        // `mem::take` (no clones) and rebuild by re-keying each in turn. A
        // bijection maps distinct keys to distinct keys, so a positional rebuild
        // needs no aggregation.
        let taken = std::mem::take(self);
        *self = taken.into_iter().map(|(k, c)| f(k, c)).collect();
    }
}

impl<K, C> MultiplyInPlace<K, C> for Vec<(K, C)>
where
    K: KeyProduct,
    C: ImaginaryUnit,
{
    #[inline]
    fn multiply_in_place(&mut self, other: &Self) {
        // The coordinate-list backend has no persistent aux; take the current
        // support out (no clone) and convolve it into `self`.
        //
        // `mem::take` leaves `self` as `Vec::new()` — capacity **zero**: the
        // buffer moved into `lhs` and is dropped at the end of this call, so the
        // accumulator would otherwise regrow from nothing on every product.
        // Reserve the product hint up front instead (the same estimate the
        // hash-join backend uses), which is also the larger allocation the
        // convolution actually needs.
        let lhs = std::mem::take(self);
        self.reserve(product_capacity_hint(lhs.len(), other.len()));
        Multiply::multiply_into(&lhs, other, self);
    }
}

impl<K, C> ApplyProducer<K, C> for Vec<(K, C)>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn apply_producer<TP>(&mut self, producer: TP)
    where
        TP: TermProducer<K, C>,
    {
        // The coordinate-list backend carries no workspace (it is the
        // small-support/reference backend, and a `Vec` field would be dead weight
        // in the layout every other operation walks), so the batch is local here.
        // The hash-join backend — the one every workload runs — keeps a
        // persistent one.
        let mut batch = TermBatch::with_capacity(self.len());
        for (k, c) in Support::iter(&*self) {
            producer.produce(&k, &c, &mut batch);
        }
        StoreAlloc::reset(self);
        Accumulate::accumulate_batch(self, &batch);
    }
}

impl<K, C> AddTerm<K, C> for Vec<(K, C)>
where
    K: Eq,
    C: Coefficient,
{
    #[inline]
    fn add_term(&mut self, key: K, coeff: C) {
        if let Some(slot) = self.iter_mut().find(|(ek, _)| *ek == key) {
            slot.1 += coeff;
        } else {
            self.push((key, coeff));
        }
    }
}

impl<K, C> AddTerm<K, C> for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn add_term(&mut self, key: K, coeff: C) {
        // Old's `ACMapAddAssign::add_assign` verbatim: `and_modify` on a hit,
        // `or_insert` otherwise — so an exactly-zero coefficient is *inserted*,
        // not dropped.
        self.entry(key)
            .and_modify(|v| *v += coeff.clone())
            .or_insert(coeff);
    }
}

impl<K, C> StoreAlloc for Vec<(K, C)> {
    #[inline]
    fn with_capacity(cap: usize) -> Self {
        Vec::with_capacity(cap)
    }

    #[inline]
    fn reset(&mut self) {
        self.clear();
    }
}

// ===========================================================================
// HashMapStore — the double-buffered hash-join backend (primary + aux + scratch).
//
// The graded algebra (`Support`/`Accumulate`/`Scale`/`Pair`/`Retain`) and the
// pure-diagonal capabilities (`ScaleByKey`/`SignFlipByKey`) delegate to the
// `primary` map's existing impls; only the buffer-using fast paths
// (`StoreAlloc`/`RekeyBijective`/`RotateInPlace`) are bespoke.
// ===========================================================================

impl<K, C> Support for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    type Key = K;
    type Coeff = C;

    #[inline]
    fn len(&self) -> usize {
        Support::len(&self.primary)
    }

    #[inline]
    fn get(&self, key: &K) -> Option<C> {
        Support::get(&self.primary, key)
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = (K, C)> {
        Support::iter(&self.primary)
    }

    #[inline]
    fn for_each_ref(&self, f: impl FnMut(&K, &C)) {
        Support::for_each_ref(&self.primary, f);
    }
}

impl<K, C> Accumulate for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn accumulate_batch(&mut self, terms: &TermBatch<K, C>) {
        Accumulate::accumulate_batch(&mut self.primary, terms);
    }

    #[inline]
    fn reduce(&mut self) {
        Accumulate::reduce(&mut self.primary);
    }
}

impl<K, C> Scale for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn scale(&mut self, s: &C) {
        Scale::scale(&mut self.primary, s);
    }
}

impl<K, C> Pair for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn probe_batch(&self, keys: &KeyBatch<K>, out: &mut [Option<C>]) {
        Pair::probe_batch(&self.primary, keys, out);
    }

    #[inline]
    fn overlap(&self, other: &Self) -> C {
        Pair::overlap(&self.primary, &other.primary)
    }

    #[inline]
    fn hermitian_overlap(&self, other: &Self) -> C
    where
        C: Conjugate,
    {
        Pair::hermitian_overlap(&self.primary, &other.primary)
    }
}

impl<K, C> Retain<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        Retain::retain(&mut self.primary, keep);
    }
}

impl<K, C> Multiply for HashMapStore<K, C>
where
    K: Indexable + KeyProduct,
    C: ImaginaryUnit,
{
    /// The twisted convolution, accumulated into `acc`'s primary support. Reads
    /// both operands' primaries; the aux/scratch buffers stay empty (this form
    /// allocates nothing of its own — the accumulator is the caller's).
    ///
    /// The accumulator is sized up front from [`product_capacity_hint`] so a
    /// `10³ × 10³` convolution does not walk a doubling chain of full rehashes
    /// over a growing ~10⁶-entry map. `reserve` is a no-op when `acc` already has
    /// the room (e.g. a caller who pre-sized it, or the second `A·C` accumulated
    /// onto an `A·B` for bilinearity).
    #[inline]
    fn multiply_into(&self, other: &Self, acc: &mut Self) {
        acc.primary.reserve(product_capacity_hint(
            self.primary.len(),
            other.primary.len(),
        ));
        Multiply::multiply_into(&self.primary, &other.primary, &mut acc.primary);
    }
}

impl<K, C> MultiplyInPlace<K, C> for HashMapStore<K, C>
where
    K: Indexable + KeyProduct,
    C: ImaginaryUnit,
{
    #[inline]
    fn multiply_in_place(&mut self, other: &Self) {
        // The old `map_add` shape, with the accumulator drawn from the store's
        // persistent double-buffer rather than a fresh allocation: clear `aux`,
        // write the whole outer product into it, swap. `aux` is left holding the
        // (now stale) previous support and is cleared on its next use, so the two
        // allocations ping-pong across products and are never freed.
        self.aux.clear();
        // Size the accumulator from the *product* estimate, not from
        // `max(|A|, |B|)`: the latter is only the lower bound (a Pauli product is
        // a bijection in each argument, so the support is at least as large as
        // either operand) and starting there made this buffer-reusing path 1.75×
        // slower than a correctly-sized accumulator — see
        // [`product_capacity_hint`] for the controlled A/B.
        self.aux.reserve(product_capacity_hint(
            self.primary.len(),
            other.primary.len(),
        ));
        Multiply::multiply_into(&self.primary, &other.primary, &mut self.aux);
        std::mem::swap(&mut self.primary, &mut self.aux);
        // Restore the "aux is empty between operations" invariant: after the swap
        // it holds the stale pre-product support. (`RekeyBijective` gets this for
        // free because it *drains* the primary; the convolution has to read the
        // primary repeatedly, so it cannot.)
        self.aux.clear();
    }
}

impl<K, C> StoreAlloc for HashMapStore<K, C> {
    #[inline]
    fn with_capacity(cap: usize) -> Self {
        Self {
            primary: HashMap::with_capacity_and_hasher(cap, IdentityBuildHasher),
            // Size the aux to the same hint so the first re-key never resizes.
            aux: HashMap::with_capacity_and_hasher(cap, IdentityBuildHasher),
            scratch: Vec::with_capacity(cap),
            // The `apply` batch is deliberately **not** pre-sized: every workload
            // shipped here drives the fused fast paths instead, so pre-sizing it
            // from the hint would charge each sum two more `cap`-sized columns
            // (at the gate bench's `1 << 20` that is tens of megabytes) for a
            // path it never takes. It still keeps its allocation across calls —
            // it grows once, on the first `apply`, and is only `clear`ed after.
            batch: TermBatch::new(),
        }
    }

    /// Reset to empty support, keeping the primary's allocation (the aux/scratch
    /// are already empty between ops).
    #[inline]
    fn reset(&mut self) {
        self.primary.clear();
    }
}

impl<K, C> AddTerm<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn add_term(&mut self, key: K, coeff: C) {
        AddTerm::add_term(&mut self.primary, key, coeff);
    }
}

impl<K, C> ApplyProducer<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn apply_producer<TP>(&mut self, producer: TP)
    where
        TP: TermProducer<K, C>,
    {
        // The produced terms are staged in the store's **persistent** batch — the
        // design's "driver-owned reusable batch", sitting alongside `aux` and
        // `scratch`. `clear` keeps both of its columns' allocations, so a gate
        // stream pays the batch allocation once rather than once per gate.
        self.batch.clear();
        for (k, c) in Support::iter(&self.primary) {
            producer.produce(&k, &c, &mut self.batch);
        }
        // Replace, not merge: a bijective re-key must not retain the old keys.
        // `reset` keeps the primary's allocation across the clear.
        StoreAlloc::reset(self);
        Accumulate::accumulate_batch(&mut self.primary, &self.batch);
        // Restore the "workspace is empty between operations" invariant the
        // struct docs state (and `Clone`/`PartialEq` rely on). `clear` keeps the
        // capacity, so the next gate still finds a warm buffer.
        self.batch.clear();
    }
}

/// Equality is the **primary support only**, exactly: same keys, same
/// coefficients — zero-coefficient terms included, since old's map equality
/// counts them (`ppvm-pauli-sum/src/sum/data.rs`, `PartialEq for PauliSum`). The
/// transient `aux`/`scratch` buffers are deliberately **not** compared: they hold
/// whatever the last operation left behind and are not part of the value.
impl<K, C> PartialEq for HashMapStore<K, C>
where
    K: Indexable,
    C: PartialEq,
{
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.primary == other.primary
    }
}

impl<K, C> ScaleByKey<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K, &mut C),
    {
        ScaleByKey::scale_by_key(&mut self.primary, f);
    }
}

impl<K, C> SignFlipByKey<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn sign_flip_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> i8,
    {
        SignFlipByKey::sign_flip_by_key(&mut self.primary, f);
    }
}

impl<K, C> RekeyBijective<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn rekey_bijective<F>(&mut self, mut f: F)
    where
        F: FnMut(K, C) -> (K, C),
    {
        // The old crate's `map_add`, with the aux living in the store: clear the
        // persistent aux, **move** every term through the re-key into it
        // (`drain` yields owned keys — zero key clones), then swap. `drain` empties
        // the primary but keeps its allocation, and the swap hands it to `aux` — so
        // the two allocations ping-pong across gates and are never freed (the
        // double-buffer, recovered).
        //
        // The merge is a plain `insert`, not `entry(..).and_modify(..).or_insert(..)`.
        // `f` is required to be **injective on keys** (this trait's contract; for a
        // Clifford conjugation it is the symplectic bijection machine-checked as
        // `*_bijective` in `lean/PPVM/Pauli/Symplectic.lean`), so two terms can never
        // land on the same re-keyed slot and there is nothing to accumulate — the
        // `entry` form was guarding a collision that cannot occur.
        //
        // That guard was not free. On this crate's monomorphization the `entry`
        // chain compiled **out of line** (a standalone `hashbrown::rustc_entry`
        // frame took 41% of this loop's samples, where the old crate's identical
        // source construct inlines fully into `map_add_assign`). It cost ~1.44x
        // against old *and* made the whole re-key bimodal across processes: with
        // `entry` ~40% of processes landed in a ~1.5x-slower mode, with `insert`
        // none of 34 did. See `docs/log.md` (`ps2.cnot.rekey.perf`,
        // `ps2.rekey.bimodal`).
        self.aux.clear();
        self.aux.reserve(self.primary.len());
        for (k, c) in self.primary.drain() {
            let (nk, nc) = f(k, c);
            let displaced = self.aux.insert(nk, nc);
            debug_assert!(
                displaced.is_none(),
                "RekeyBijective requires an injective re-key: two distinct keys \
                 collided, so a term was silently dropped"
            );
        }
        std::mem::swap(&mut self.primary, &mut self.aux);
    }
}

impl<K, C> RotateInPlace<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn rotate_in_place<F>(&mut self, mut f: F)
    where
        F: FnMut(&K, &mut C) -> Option<(K, C)>,
    {
        // Pass 1: scale each diagonal coefficient in place (the key is unchanged,
        // so its cached hash stays valid and the entry never moves), buffering the
        // anticommuting branch terms in the persistent `scratch` — the old crate's
        // reused `map_insert` scratch, now living in the store (no per-gate Vec
        // allocation).
        self.scratch.clear();
        for (k, c) in self.primary.iter_mut() {
            if let Some(term) = f(k, c) {
                // Pre-warm the fresh branch key's structural digest HERE, in pass 1,
                // before it is buffered. `with_bits_toggled` builds the key with an
                // empty (`HASH_UNCACHED`) cache, so without this the 3-round finalize
                // fold would fire lazily inside pass-2's `entry()` — *on the bucket-
                // index critical path*, where the mul-chain latency stalls the
                // dependent bucket load with nothing to hide it. Computing the (same)
                // digest here instead lets it overlap with the walk's other work
                // (the next term's diagonal scale + branch build), so pass-2 probes a
                // cache hit. Semantic no-op (identical digest); measured ~8% off `rx`
                // on a ~1000-term sum, matching the old crate's eager-hash `map_insert`
                // (which likewise hashes in its first pass).
                let _ = term.0.key_hash();
                self.scratch.push(term);
            }
        }
        // Pass 2: merge only the branch terms through the hash-join, aggregating
        // any collision with an existing (already-scaled) diagonal or a sibling
        // branch — one hash pass over ≤N terms, not the whole 2N fan-out.
        //
        // Every produced branch is merged, *including* an exactly-zero one. Old's
        // `map_insert` pass 2 calls `data.add_assign(k, v)` for every term pass 1
        // produced (`ppvm-pauli-sum/src/sum/data.rs`), and `add_assign` is
        // `entry().and_modify().or_insert(v)` — it inserts a `0.0`. So an identity
        // rotation `rx(q, 0.0)` on a `Z`-bearing term leaves old with
        // `{Z: 1.0, Y: 0.0}`, and old's exact-map `PartialEq` counts the `Y`;
        // skipping it here would change `len()`, `get()`, the key set and equality.
        // No `reduce` scan runs either: a generic rotation's collision
        // cancellations are — like the old crate's `map_insert`
        // (`ppvm-pauli-sum::sum::rot1`), which leaves any residue for a later
        // truncation — left for the policy's truncation (or the caller's magnitude
        // floor), not dropped here.
        //
        // `drain(..)` empties `scratch` but keeps its allocation for the next gate.
        // The merge is [`AddTerm`] — the single definition of old's `add_assign`.
        for (nk, nc) in self.scratch.drain(..) {
            AddTerm::add_term(&mut self.primary, nk, nc);
        }
    }
}
