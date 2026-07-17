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
both traits, while a mutable intermediate such as `Phased<_, NoHash>` can still
implement `Word` without pretending to be a valid map key. This separation
replaces the current `REHASH = false` use case with an explicit non-key mode.

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

An ordinary word implements `Word<Site = Pauli>`. `WithLoss` implements
`Word<Site = LossySite<Pauli>>` and returns `Lost` for marked sites. This keeps
the word interface independent of the chosen operator alphabet.

## `PauliWord` packed representation

The initial packed implementation stores parallel fixed-size X and Z arrays:

```rust
pub struct PauliWord<A, S, M = LazyHashU64> {
    xbits: BitArray<A>,
    zbits: BitArray<A>,
    nqubits: usize,
    hash_cache: M,
    _phantom: PhantomData<S>,
}
```

`A` is an implementation parameter such as `[u8; N]` or `[usize; N]`. It is
not exposed through `Word`. The implementation validates that `nqubits` fits
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

Loss is best modeled as an orthogonal wrapper around a Pauli word:

```rust
pub struct WithLoss<W, A, M = CompositeHash>
where
    W: Word<Site = Pauli>,
{
    word: W,
    lbits: BitArray<A>,
    loss_hash: M,
}
```

`A` is the packed bit-array backing selected by the concrete type. The public
alias retains the current `LossyPauliWord` name:

```rust
pub type LossyPauliWord<A, H> = WithLoss<PauliWord<A, H>, A>;
```

`WithLoss` is a justified new implementation name because the proposal changes
the current standalone lossy representation into a generic wrapper. The
established domain-facing alias does not need to change.

Neither `A` nor the underlying word's array type appears in `Word`.

### Canonical loss invariant

A lost site must contain identity in the underlying Pauli word:

```text
lost[q] = 1  =>  word[q] = I
```

`set_lost(q)` first sets the underlying word to `I`, then sets the loss bit.
`set(q, LossySite::Present(p))` clears the loss bit and then writes `p`. This
prevents multiple physical encodings from representing the same logical lossy
word.

The structural identity is:

```text
(underlying Pauli word, logical loss bits)
```

`weight()` counts `X`, `Y`, `Z`, and `Lost`; `loss_weight()` counts only lost
sites.

### Loss-specific behavior

Generic lossy Pauli propagation sees `LossySite::Lost` through `Word` and
preserves or skips the site according to the operation's semantics. Operations
that create, clear, or count loss use inherent `WithLoss` methods:

```rust
impl<W, A, M> WithLoss<W, A, M>
where
    W: Word<Site = Pauli>,
{
    pub fn is_lost(&self, qubit: usize) -> bool;
    pub fn set_lost(&mut self, qubit: usize);
    pub fn clear_loss(&mut self, qubit: usize);
    pub fn loss_weight(&self) -> usize;
}
```

Loss channels and maximum-loss-weight truncation specialize directly on the
generic wrapper family:

```rust
impl<C, W, A, M, S, P> LossChannel<C>
    for OperatorSum<C, WithLoss<W, A, M>, S, P>
where
    W: Word<Site = Pauli>,
{
    // ...
}
```

There are no traits named `PauliWord`, `LossyPauliWord`, or `FermionWord`.
`PauliWord` and `LossyPauliWord` remain concrete domain type names. Algorithms
select the algebra through `Word::Site`; loss-only operations remain inherent
to `WithLoss`.

## Phased words

Phase is another orthogonal wrapper:

```rust
pub struct Phased<W, M = CompositeHash>
where
    W: Word,
    M: HashMode<W>,
{
    word: W,
    phase: Phase,
    hash_state: M::State,
}
```

`Phased` is the generalized wrapper implementation. Pauli-facing code retains
the established name through an alias:

```rust
pub type PhasedPauliWord<W, M = CompositeHash> = Phased<W, M>;
```

For Pauli use, `Phase` represents `+1`, `+i`, `-1`, and `-i`. The wrapper may
also be useful for other word algebras, but algebra-specific multiplication is
implemented only under the appropriate specialized word bound.

Loss and phase compose without a new combined representation:

```rust
PhasedPauliWord<WithLoss<PauliWord<A, H>, A>, Mode>
```

The wrapper fields are private so phase mutation cannot bypass component-cache
invalidation.

## Hash ownership

Every hash-enabled word implements `Indexable` and selects:

- a build hasher;
- an `Indexable::HashCache` associated type;
- concrete cache fields; and
- the algorithm that hashes its private structural data.

The associated cache type does not require every implementation to use the
same representation. Examples include:

- a plain unsigned integer with a separate validity flag for eager or explicit
  recomputation;
- `Cell<u64>` plus `Cell<bool>` for lazy single-threaded hashing;
- atomics for a `Sync` key; and
- `()` for a mode that deliberately provides no cache.

Because `Hash::hash` receives `&self`, a cache populated lazily from that method
requires interior mutability. Such a representation will generally prevent
the containing word from being `Copy`.

## Component hashes

Hash composition follows the logical wrappers:

```text
packed Pauli hash = hash(nqubits, X bits, Z bits)
lossy hash        = combine(Pauli hash, loss hash)
phased hash       = combine(inner-word hash, phase hash)
phased lossy hash = combine(Pauli hash, loss hash, phase hash)
```

`combine` must be ordered and domain-separated. It must not be an
unqualified XOR of arbitrary component digests.

The phase has only four values, so `CompositeHash` can normally compute its
contribution from a small table or mixer without another cache. `CachedHash`
is justified only if profiling identifies a caller that benefits from caching
the combined phased value.

Loss masks may be large, so `WithLoss` can cache the loss component separately
from the Pauli component. A loss-only mutation then avoids rehashing X/Z.

## Invalidation rules

| Mutation | Pauli component | Loss component | Phase component |
| --- | --- | --- | --- |
| Change ordinary Pauli site | invalidate | preserve | preserve |
| Mark identity site lost | preserve | invalidate | preserve |
| Mark nonidentity site lost | invalidate | invalidate | preserve |
| Clear loss to identity | preserve | invalidate | preserve |
| Replace loss with Pauli | invalidate if nonidentity | invalidate | preserve |
| Change phase | preserve | preserve | recompute or invalidate |

Constructors compute caches eagerly or mark them invalid according to the
selected cache mode. Cloning copies a valid cache because the clone initially
has identical structural contents.

## Hash modes

Wrappers that are sometimes used only as mutable intermediate values should be
generic over hash behavior:

```rust
pub trait HashMode<T> {
    type State: Clone;
}

pub struct NoHash;
pub struct CompositeHash;
pub struct CachedHash<C>(PhantomData<C>);
```

- `NoHash` stores no wrapper cache and does not implement `Hash` or
  `Indexable` for the wrapper.
- `CompositeHash` implements `Hash` by combining cached inner components on
  demand, with no combined cache.
- `CachedHash` stores a lazily or eagerly maintained combined cache.

The concrete policy protocol should be finalized while implementing the first
two modes. It should not expose cache state to propagation algorithms.

## Ordering and serialization

`Eq`, `Ord`, `Hash`, display, parsing, and serialization must agree on logical
identity:

- Pauli sites compare in a documented order.
- Loss participates after the underlying Pauli content or through an explicit
  `LossySite<Pauli>` ordering.
- Phase participates only when the phased wrapper's identity includes it.
- Unused bits and cache state never participate.

Serialization uses logical symbols and lengths, not raw native-word memory, so
it remains stable across storage widths and platforms.

## Prototype validation

The prototype should include:

- round-trip parsing tests for ordinary and lossy symbols;
- property tests comparing packed operations with a simple symbol vector;
- tests enforcing `lost => underlying I` after every mutator;
- equality/hash agreement tests for all wrappers and modes;
- tests showing loss-only changes preserve the Pauli hash component;
- tests showing phase-only changes preserve Pauli and loss components;
- tests proving unused high bits do not affect identity;
- `Send`/`Sync` assertions for cache modes intended for concurrent maps; and
- benchmarks comparing eager, lazy, composite, and uncached hashing.

## Open questions

1. Should `WithLoss` permit a loss-mask storage width different from its inner
   word's packed width?
2. Does any real phased-word caller benefit from a combined `CachedHash` mode?
3. Does `Indexable::HashCache` have a generic consumer, or should the cache
   type also remain entirely private like word storage?
