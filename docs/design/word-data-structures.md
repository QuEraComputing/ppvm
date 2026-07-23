# Word data structures

Status: design sketch

## Purpose

This document describes concrete data structures for standalone algebraic
words. The shared trait design is specified in
[`traits-2-configuration-and-hashing.md`](traits-2-configuration-and-hashing.md).

The trait-level `Word` is representation-free but defines the common indexed
product operations:

```rust
pub trait Word {
    type Site;

    fn n_sites(&self) -> usize;
    fn get(&self, index: usize) -> Self::Site;
    fn set(&mut self, index: usize, site: Self::Site);
    fn weight(&self) -> usize;
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
`n_sites()`. `set()` is the shared structural mutation boundary and must
preserve concrete invariants and invalidate the affected hash components.

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
is the build hasher associated through `Indexable`. Neither is exposed through
`Word`. The implementation validates that `nqubits` fits the arrays and
ignores or canonicalizes unused high bits.

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
    xz_hash_cache: OnceLock<u64>,
    loss_hash_cache: OnceLock<u64>,
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
impl<C, A, H, S, P> LossChannel<C>
    for OperatorSum<C, LossyPauliWord<A, H>, S, P>
{
    // ...
}
```

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

Every hash-enabled word implements `Indexable`, associates its build hasher,
and privately owns the fields and algorithm used to cache its structural
hash. Cache representation and invalidation are not exposed through
`Indexable`.

The shipped indexable words use component `OnceLock<u64>` caches. `Hash::hash`
can populate them through `&self`; structural mutators clear affected cells
through `&mut self`. This preserves `Send + Sync` without imposing either
bound on the `Indexable` trait.

## Component hashes

Hash composition follows the logical wrappers:

```text
packed Pauli hash = hash(nqubits, X bits, Z bits)
lossy hash        = combine(Pauli hash, loss hash)
```

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
