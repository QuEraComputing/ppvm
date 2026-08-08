# Word data structures

Status: design sketch

## Purpose

This document describes concrete data structures for standalone algebraic
words. The shared trait design is specified in
[`traits-2-configuration-and-hashing.md`](traits-2-configuration-and-hashing.md).

The trait-level `Word` is representation-free and **read-only**: it defines the
common indexed *inspection* operations, not mutation. Structural mutation lives
on algebra-specific traits (`PauliBits`, `SymplecticColumns`, `PhaseTrack`; see
[`traits-2-configuration-and-hashing.md`](traits-2-configuration-and-hashing.md#pauli-algebra-traits-symplectic-structure-and-phase)):

```rust
pub trait Word {
    type Site;

    fn n_sites(&self) -> usize;
    fn get(&self, index: usize) -> Self::Site;
    fn weight(&self) -> usize;
    fn iter(&self) -> impl Iterator<Item = Self::Site>;
}
```

There is no `Word::Storage` associated type. Packed arrays, loss masks, ordered
factors, hash fields, and validity flags are private fields of concrete types.
Generic propagation and collection code uses behavioral traits and never names
the backing memory.

`Word` does not extend `Indexable`. A word used as an `ACMap` key implements
both traits, while a mutable intermediate such as `Phased<W>` can implement
`Word` without pretending to be a valid map key. This separation replaces the
current `REHASH = false` use case with an explicitly non-indexable type.

This document initially focuses on:

- packed Pauli words;
- lossy Pauli words;
- phased-word composition; and
- lazy structural hash caches.

Fermion-word storage will receive a separate design when fermionic propagation
is implemented. It will use the same `Word` interface with a different `Site`:
the word index records product order, while the fermionic site value records
the physical mode and creation or annihilation action.

`weight()` counts non-identity factors according to the selected site type.
For a representation that stores only non-identity factors, it may equal
`n_sites()`. The structural mutation boundary is `PauliBits::set_x_bit` /
`set_z_bit` (and the `SymplecticColumns` / `PhaseTrack` column primitives), each
of which preserves concrete invariants and lazily invalidates the affected hash
components. `PauliWord` and `LossyPauliWord` implement `PauliBits`, and also
`SymplecticColumns + PhaseTrack` plus the opt-in marker `BlanketClifford`, so they
pick up the shared blanket `Clifford` impl. A bare word carries no phase field, so
its `PhaseTrack` is a phase-discarding no-op and its `Clifford` is the pure
`Sp(2n, 2)` bit map (matching the old bare-word bit-only Clifford).
`PhasedPauliWord` instead carries a hand-written **fused** `impl Clifford`: it
reads each inner X/Z bit once via `PauliBits`, computes the `ℤ₄` conjugation sign
(`Y ↦ −Y`, …), applies the bit update reusing those reads, and folds the sign into
the stored phase — so it does **not** implement `BlanketClifford` (the marker that
gates the blanket), avoiding the blanket's redundant double bit read while keeping
the exact old-kernel signs. This matches the authoritative trait assignment
in [`traits-2-configuration-and-hashing.md`](traits-2-configuration-and-hashing.md#pauli-algebra-traits-symplectic-structure-and-phase),
where a Clifford gate on a sum applies the one-row `SymplecticColumns` action
pointwise and drains each term's phase delta to its coefficient.

## Logical Pauli model

An ordinary Pauli word is a fixed-width tensor product over `I`, `X`, `Y`, and
`Z`. Each site is represented by two logical bits:

| Pauli | X bit | Z bit |
| --- | --- | --- |
| `I` | 0 | 0 |
| `X` | 1 | 0 |
| `Z` | 0 | 1 |
| `Y` | 1 | 1 |

A lossy Pauli word adds a `Lost` state. Loss is exclusive with the four Pauli
states; a lost site has canonical bits `(x, z, lost) = (0, 0, 1)`.

The exact Pauli alphabet and its lossy extension are site types rather than
word subtraits:

```rust
pub enum Pauli {
    I,
    X,
    Y,
    Z,
}

pub enum LossySite<S> {
    Present(S),
    Lost,
}
```

An ordinary word implements `Word<Site = Pauli>`. `LossyPauliWord` implements
`Word<Site = LossySite<Pauli>>` and returns `Lost` for marked sites. This keeps
the word interface independent of the chosen operator alphabet.

## `PauliWord` packed representation

The initial packed implementation stores parallel fixed-size X and Z arrays:

```rust
pub struct PauliWord<A, H> {
    xbits: BitArray<A>,
    zbits: BitArray<A>,
    nqubits: usize,
    hash_cache: OnceLock<u64>,
    _hasher: PhantomData<fn() -> H>,
}
```

`A` is an implementation parameter such as `[u8; N]` or `[usize; N]`, and `H`
is the private internal digest algorithm (`fxhash`, `gxhash`, …) used to compute
`key_hash()`. `H` is no longer a `BuildHasher` associated on `Indexable` — the
direct-digest model gives the map no hasher to pick — but a representation
parameter of the concrete word, on the same footing as `A`. Neither is exposed
through `Word` or `Indexable`. The implementation validates that `nqubits` fits
the arrays and ignores or canonicalizes unused high bits.

The structural identity is:

```text
(nqubits, logical X bits, logical Z bits)
```

Equality and hashing exclude unused capacity, cache contents, cache validity,
and the `PhantomData` marker.

All fields are private. Mutations go through methods that preserve unused-bit
invariants and invalidate the hash when a logical Pauli site changes.

## Lossy Pauli word

The first prototype keeps the established lossy Pauli word as a flattened
packed representation:

```rust
pub struct LossyPauliWord<A, H> {
    xbits: BitArray<A>,
    zbits: BitArray<A>,
    lbits: BitArray<A>,
    nqubits: usize,
    xz_hash_cache: AtomicU64,
    loss_hash_cache: AtomicU64,
    combined_hash_cache: AtomicU64,
    _hasher: PhantomData<fn() -> H>,
}
```

Inlining all three planes avoids wrapper nesting and keeps the lossy hot path
and component hashes direct. A generic loss wrapper is not introduced until a
second real word representation demonstrates that it needs the same
composition. `A` and `H` remain private implementation parameters from the
perspective of `Word`.

### Canonical loss invariant

A lost site must contain identity in its X/Z planes:

```text
lost[q] = 1  =>  xbits[q] = 0 and zbits[q] = 0
```

`set_lost(q)` first clears the X/Z bits, then sets the loss bit.
`set(q, LossySite::Present(p))` clears the loss bit and then writes `p`. This
prevents multiple physical encodings from representing the same logical lossy
word.

The structural identity is:

```text
(nqubits, logical X bits, logical Z bits, logical loss bits)
```

`weight()` counts `X`, `Y`, `Z`, and `Lost`; `loss_weight()` counts only lost
sites.

A Clifford gate must leave lost qubits untouched (a lost qubit carries no operator
to conjugate). Because the blanket `Clifford` composes `SymplecticColumns`
primitives directly (there is no gate-level hook), the guard lives in each
primitive: `LossyPauliWord`'s `SymplecticColumns` ops (`swap_xz`,
`xor_z_from_x`, `xor_x_col`, `xor_z_col`, `cz_bits`) are **no-ops when any
involved qubit is lost**, reproducing the old crate's whole-gate skip and
preserving the canonical `lost ⇒ (x, z) = (0, 0)` invariant. So the lossy word's
`SymplecticColumns` is a *loss-guarded* `Sp(2n, 2)` map, not the literal pure map
the bare word uses. This is machine-checked in `lean/PPVM/Pauli/Symplectic.lean`:
each guarded primitive preserves the loss invariant (`xorXColL_preserves_loss`,
`xorZColL_preserves_loss`), the two guarded primitives compose to the atomic
whole-gate skip (`xorZColL_xorXColL_eq_cnotActL`), and on present qubits the
guarded gate is still the proven `Sp` isometry (`cnotActL_present_isometry`,
`czActL_present_isometry`). The same holds for `CY`, which the blanket
`CliffordExtensions` decomposes into `s(t); cnot(c,t); s_dag(t)` — three
primitives whose guards deliberately differ (`sActL` tests `lost t` alone,
`cnotActL` tests `lost c ∨ lost t`), so a **lost control with a present target**
skips the atomic gate while still running two `S(t)` conjugations that must
cancel exactly: `sActL_cnotActL_sActL_eq_cyActL` proves the guarded composite
equals the old crate's atomic `cy` skip (`cyActL`) on every loss configuration,
with `cyActL_preserves_loss` and `cyActL_present_isometry` the invariant and
isometry halves.

### Loss-specific behavior

Generic lossy Pauli propagation sees `LossySite::Lost` through `Word` and
preserves or skips the site according to the operation's semantics. Operations
that create, clear, or count loss use inherent `LossyPauliWord` methods:

```rust
impl<A, H> LossyPauliWord<A, H> {
    pub fn is_lost(&self, qubit: usize) -> bool;
    pub fn set_lost(&mut self, qubit: usize);
    pub fn clear_loss(&mut self, qubit: usize);
    pub fn loss_weight(&self) -> usize;
}
```

Loss channels and maximum-loss-weight truncation specialize directly on the
concrete lossy word:

```rust
impl<S, P, A, H> LossChannel for Sum<S, P>
where
    S: SumStorage<Word = LossyPauliWord<A, H>>,
    P: Policy<S::Word, S::Coeff>,
{
    // ...
}
```

The specialization is expressed through the storage's associated `Word` type
rather than a free word parameter: a `Sum` is a `LossChannel` exactly when its
`SumStorage` keys on a `LossyPauliWord`.

There are no traits named `PauliWord`, `LossyPauliWord`, or `FermionWord`.
`PauliWord` and `LossyPauliWord` remain concrete domain type names. Algorithms
select the algebra through `Word::Site`; loss-only operations remain inherent
to `LossyPauliWord`.

## Phased words

Phase is another orthogonal wrapper:

```rust
pub struct Phased<W>
where
    W: Word,
{
    word: W,
    phase: Phase,
}
```

`Phased` is the generalized wrapper implementation. Pauli-facing code retains
the established name through an alias:

```rust
pub type PhasedPauliWord<W> = Phased<W>;
```

For Pauli use, `Phase` represents `+1`, `+i`, `-1`, and `-i`. The wrapper may
wrap both ordinary and lossy Pauli words. It may also be useful for other word
algebras, but algebra-specific multiplication is implemented only under the
appropriate specialized word bound.

Loss and phase compose without a new combined representation:

```rust
PhasedPauliWord<LossyPauliWord<A, H>>
```

No phased word is a production map key in the first prototype, so `Phased<W>`
does not implement `Hash` or `Indexable` and stores no hash mode or cache.

## Hash ownership

Every hash-enabled word implements `Indexable` and privately owns its internal
digest algorithm, the fields used to cache its structural hash, and the
finalization fold that makes `key_hash()` avalanche-quality. Cache
representation, the algorithm choice, and invalidation are not exposed through
`Indexable`, which surfaces only the finalized digest value via `key_hash()`.

The shipped indexable words use relaxed-atomic sentinel caches. `key_hash()`
(and the `Hash` impl, which is just `state.write_u64(self.key_hash())`) can
populate them through `&self`; structural mutators clear affected cells through
`&mut self`. The cached value is a pure function of immutable structural fields,
so racing misses compute the same digest and require no compare-and-swap. This
preserves `Send + Sync` without imposing either bound on `Indexable`.

The finalization fold is applied per-algorithm and per-width by the private
`HashFinalize` helper: narrow storage under a weak hasher (`[u8; 8]` +
`fxhash`) folds `raw ^ (raw >> 32)` so the low bits — hashbrown's bucket —
decorrelate, while a strong hasher (`gxhash`) folds nothing. This helper lives
in `ppvm-pauli-word`, not in the algebra-agnostic trait crate.

The same component caches back the bulk `hash_into` on a word's key column (see
[Key columns](#key-columns-structure-of-arrays-batches)): filling a batch's hash
column mostly gathers already-cached values, so the group-prefetch loop in the
[batch contract](traits-2-configuration-and-hashing.md#batch-execution-and-the-hash-join-contract)
pays no hashing on its critical path.

## Component hashes

Hash composition follows the logical wrappers:

```text
packed Pauli hash = hash(X bits, Z bits)
lossy hash        = combine(Pauli hash, loss hash)
```

Width remains part of structural equality but is omitted from the digest,
matching the legacy bucket distribution. Different-width words may therefore
collide, which is valid under the hash contract; a `Sum` already enforces one
common width across its support.

`combine` must be ordered and domain-separated. It must not be an
unqualified XOR of arbitrary component digests.

Loss masks may be large, so `LossyPauliWord` caches the loss component
separately from the X/Z component. A loss-only mutation then avoids rehashing
X/Z. `Phased<W>` is absent from this composition because it is not indexable.

## Invalidation rules

| Mutation | X/Z component | Loss component |
| --- | --- | --- |
| Change ordinary Pauli site | invalidate | preserve |
| Mark identity site lost | preserve | invalidate |
| Mark nonidentity site lost | invalidate | invalidate |
| Clear loss to identity | preserve | invalidate |
| Replace loss with Pauli | invalidate if nonidentity | invalidate |

Constructors leave caches empty. Cloning may copy a valid cached value because
the clone initially has identical structural contents.

## Key columns (structure-of-arrays batches)

The bulk map contract in
[`traits-2-configuration-and-hashing.md`](traits-2-configuration-and-hashing.md#batch-execution-and-the-hash-join-contract)
consumes a structure-of-arrays column of keys through the `Columnar` trait's
associated `Column`. Each concrete word owns the planar layout of its column;
generic code sees only the `KeyColumn` capability, never the planes. `Columnar`
is separate from `Indexable`, so this adds no surface to the minimal hashing
contract.

For the packed Pauli word the column is two plane blocks and a shared width,
with the hash column stored parallel to it:

```text
x_plane:  [ word 0 X-bits | word 1 X-bits | ... | word M-1 X-bits ]
z_plane:  [ word 0 Z-bits | word 1 Z-bits | ... | word M-1 Z-bits ]
nqubits:  shared width
hashes:   [ h(0) | h(1) | ... | h(M-1) ]   (in the owning KeyBatch)
```

The flattened `LossyPauliWord` extends its column with a parallel loss-bit
plane, matching its X/Z/loss structural identity. `Phased<W>` has no column: it
is not indexable, so it is never a batch key.

The layout is chosen for the three access patterns the trait contract must
serve:

- **Vertical SIMD.** A contiguous plane lets a backend hash or compare the same
  machine-word lane across many keys in one vector instruction, instead of
  gathering scattered `[u8; N]` words.
- **Coalesced device access.** A contiguous plane maps to a coalesced load for a
  GPU warp; an array of packed words would not.
- **Bandwidth-minimal probe.** The probe streams the X/Z planes and the hash
  column and never loads coefficients, which live in the caller's `TermBatch`.

Two representation choices are therefore fixed here rather than in the trait
layer:

- planes are stored as flat plane blocks, not `Vec<[u8; N]>`, and each key's
  plane slot is padded to the alignment the backend's widest SIMD lane (or a
  device's coalescing width) requires; and
- `hash_into` computes the whole hash column directly from the planes using the
  same structural fold the scalar word caches in relaxed atomics, so the
  plane-parallel hash and scalar `Hash` in
  [Component hashes](#component-hashes) agree bit for bit.

`gather` selects or permutes columns by index for radix partitioning, truncation
compaction, and device staging; it copies plane by plane and never materializes
scalar words.

## Ordering and serialization

`Eq`, `Ord`, `Hash`, display, parsing, and serialization must agree on logical
identity:

- Pauli sites compare in a documented order.
- Loss participates after the underlying Pauli content or through an explicit
  `LossySite<Pauli>` ordering.
- Phase participates in equality and serialization for `Phased<W>`, but not in
  map-key hashing because the wrapper is not indexable.
- Unused bits and cache state never participate.

Serialization uses logical symbols and lengths, not raw native-word memory, so
it remains stable across storage widths and platforms.

**First-prototype scope.** The shipped `-2` words implement `Eq`, `Hash`,
`Display`, and parsing (the surface `Indexable`/`KeyColumn` require —
`Clone + Eq + Hash`), but **not** `Ord` or (de)serialization yet: no in-scope
consumer needs them for the prototype. When added, they must satisfy the
agreement above (a documented site order; loss and phase participating as
described). Deferred deliberately, tracked here rather than silently dropped.

## Prototype validation

The prototype should include:

- round-trip parsing tests for ordinary and lossy symbols;
- property tests comparing packed operations with a simple symbol vector;
- tests enforcing `lost => X/Z identity` after every mutator;
- equality/hash agreement tests for ordinary and lossy indexable words;
- tests showing loss-only changes preserve the X/Z hash component;
- equality and serialization tests for ordinary and lossy `Phased<W>` values;
- tests proving unused high bits do not affect identity;
- `Send`/`Sync` assertions for shipped indexable words; and
- benchmarks comparing uncached structural hashing with the private lazy
  component caches.

## Open questions

1. Should the loss plane use the same packed array width as the X/Z planes, or
   a separately selected private width?
2. What plane alignment and padding does the key column use, and must it vary by
   target (for example AVX-512 vs NEON vs a GPU coalescing width)? A per-target
   column would need the plane block parameterized by alignment.
3. Does the `LossyPauliWord` column store its loss plane inline beside X/Z (one
   allocation, one `gather`), matching the flattened scalar representation?
