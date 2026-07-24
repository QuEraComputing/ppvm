# `ppvm-traits-2`: type composition, indexable values, and cached hashing

Status: design sketch

## Motivation

The current `ppvm_traits::Config` bundles choices from several unrelated
layers:

- coefficient type;
- packed Pauli storage;
- Pauli-word representation;
- key hasher;
- truncation strategy; and
- concrete map implementation.

This makes the foundational configuration specific to `PauliSum`, even though
the gate and noise traits are shared by other algorithms. In particular, a map
and a truncation strategy are not properties of quantum data. They are choices
made by a particular algorithm.

The second trait-system experiment should separate:

1. the coefficient type, passed directly as a generic parameter;
2. concrete quantum-data representations;
3. the hashing contract of values that can be used as keys; and
4. algorithm-specific storage and policy choices.

The redesign should remain compile-time generic. It should not introduce
runtime trait objects or runtime dispatch for these choices.

This proposal was compared against the current definitions in
`ppvm-traits/src/config.rs`, `traits/map.rs`, `traits/strategy.rs`, and
`traits/word_trait.rs`; `ppvm-pauli-sum/src/sum/data.rs`; the concrete word
types in `ppvm-pauli-word`; and `ppvm-tableau-sum`'s `EntryStore`, `VecStorage`,
and `MapStorage`. Existing names are retained below unless the abstraction or
responsibility changes.

## Type-composition layers

### Coefficient type

There is no algorithm-agnostic `Config` trait. Algorithms use the coefficient
type directly:

```rust
pub struct SomeAlgorithm<C: Coefficient> {
    // ...
}
```

A trait containing only `type Coeff` adds indirection without providing useful
composition. Maps, policies, Pauli-word storage, Pauli-word
implementations, and hashers are also selected independently rather than being
collected into a replacement global configuration trait.

### Representation types

There is no separate global storage configuration, and representation storage
does not appear as an associated type. A concrete value encapsulates its own
fields; generic algorithms use behavioral methods instead of naming or
inspecting the backing memory.

`Word` is the common concept for an indexed algebraic monomial. The old
Pauli-word operations are not Pauli-specific: every supported word has a site
alphabet, an indexed extent, site access and mutation, and a weight:

```rust
pub trait Word {
    type Site;

    fn n_sites(&self) -> usize;
    fn get(&self, index: usize) -> Self::Site;
    fn set(&mut self, index: usize, site: Self::Site);
    fn weight(&self) -> usize;
}
```

`Site` selects the operator alphabet without introducing `PauliWord` or
`FermionWord` subtraits. For example, the relevant concrete types may be:

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

pub struct FermionSite {
    pub mode: usize,
    pub action: FermionAction,
}
```

Thus an ordinary packed Pauli word implements `Word<Site = Pauli>`, a
concrete packed lossy word implements `Word<Site = LossySite<Pauli>>`, and a
future ordered fermionic product can implement `Word<Site = FermionSite>`. A
fermionic word's index denotes factor order; `FermionSite` carries the physical
mode. For a dense Pauli word, the index is the qubit and `n_sites()` is its
width.
`weight()` is the number of non-identity factors according to the concrete site
alphabet; an ordered representation that stores no explicit identities may
therefore have `weight() == n_sites()`. Implementations of `set()` preserve
their representation invariants and invalidate the affected hash components.

The concrete `LossyPauliWord` stores packed X, Z, and loss planes directly and
provides loss mutation and `loss_weight()` as inherent methods. Loss channels
and loss-specific truncation specialize directly on that concrete type; there
is no one-implementation `LossyPauliWord` capability trait. A generic loss
wrapper should be reconsidered only after a second real word representation
needs the same composition.

Concrete packed Pauli, lossy, phased, and hash-cache layouts are described in
[`word-data-structures.md`](word-data-structures.md). None of those layouts is
visible through `Word`.

A tableau is an independent concrete representation. It does not contain a
public row type, implement `Word`, or select an associated word
implementation. Its X/Z matrices, phases, orientation, and contiguous backing
allocation are private implementation details described in
[`tableau-data-structure.md`](tableau-data-structure.md).

### Behavioral traits

Shared traits describe operations, not representation layout. Clifford gates
need no coefficient parameter. Operations that consume numeric parameters use
the coefficient type directly:

```rust
pub trait Clifford {
    fn h(&mut self, qubit: usize);
    fn cnot(&mut self, control: usize, target: usize);
    // ...
}

pub trait RotationOne<C: Coefficient> {
    fn rx(&mut self, qubit: usize, theta: C);
    // ...
}

pub trait PauliError<C: Coefficient> {
    fn pauli_error(&mut self, qubit: usize, probabilities: [C; 3]);
}

pub trait Measure {
    fn measure(&mut self, qubit: usize) -> Option<bool>;

    fn measure_many(&mut self, targets: &[usize]) -> Vec<Option<bool>> {
        targets.iter().map(|&q| self.measure(q)).collect()
    }
}
```

The same concrete tableau may implement a numeric trait for every supported
coefficient type without storing that coefficient type. Measurement and reset
traits likewise expose behavior and result types without exposing tableau
rows, packed matrix blocks, or matrix orientation. `Measure` is loss-aware for
both `Tableau` and `GeneralizedTableau`: `Some(false)` and `Some(true)` denote
computational-basis outcomes and `None` denotes a lost qubit. The former
`Measure -> bool` and `LossyMeasure -> Option<bool>` split is removed.
Python may continue translating this Rust representation to its
`MeasurementResult` enum.

Sharing the result type does not share the measurement algorithm. `Tableau`
uses the pure Clifford measurement procedure. `GeneralizedTableau` performs
its coefficient-aware stabilizer/destabilizer decomposition and update, which
is always \(O(n^2)\). The shared behavioral trait must not force the latter
algorithm or its scratch state into the concrete tableau.

There is no shared `TableauStorage` trait in the first design. If multiple
tableau implementations are later useful, each concrete type can implement
the same behavioral traits. A storage abstraction should only be introduced
after two implementations demonstrate a common interface.

### Algorithm and storage parameters

An algorithm should take its independent choices as direct type parameters.
An associated-type bundle is not useful merely because it replaces two type
parameters with one. In particular, there is no `PauliSumAlgorithm` trait that
bundles a term map with a policy: storage layout and policy are orthogonal
choices.

The reusable sparse-sum shape is:

```rust
pub struct OperatorSum<C, W, S, P = NoPolicy>
where
    C: Coefficient,
    W: Word + Indexable,
    S: SumStorage<W, C>,
    P: Policy<W, C>,
{
    storage: S,
    policy: P,
    n_sites: usize,
}
```

Here `C`, `W`, `S`, and `P` respectively mean coefficient domain, algebraic
word, concrete sparse-sum storage engine, and algorithm policy. Each
parameter has an independent meaning. Propagation methods select their algebra
through the site type, for example `W: Word<Site = Pauli>` or
`W: Word<Site = FermionSite>`.

`Policy` is the proposed name for the current `Strategy` concept. It retains
the current responsibilities: predicting initial capacity and truncating the
sum. Existing concrete strategies become policies without otherwise changing
their meaning; `NoStrategy` and `CombinedStrategy` become `NoPolicy` and
`CombinedPolicy`, while `MaxPauliWeight` and `CoefficientThreshold` keep their
established names:

```rust
pub trait Policy<W, C>: Default + Clone + Copy
where
    W: Word + Indexable,
    C: Coefficient,
{
    fn capacity(&self, n_sites: usize) -> usize;

    fn truncate<M>(&self, map: &mut M)
    where
        M: ACMap<W, C>;
}
```

`Policy` and its concrete implementations belong to the sparse-sum crate. This
removes the current split where the `Strategy` trait lives in `ppvm-traits` but
its concrete strategies live in `ppvm-pauli-sum`; the policy is not an
algorithm-agnostic `ppvm-traits-2` concern.

`ACMap` remains the name of the associative coefficient-map capability already
implemented by `HashMap`, `IndexMap`, and `DashMap`. Its generic signature can
be simplified after `PauliStorage` and the separately supplied build hasher are
removed, but coefficient accumulation, iteration, insertion, retention, and
consumption are the same concept. `ACMap` moves with the sparse-sum engine (the
existing `ppvm-pauli-sum` initially) rather than being renamed or kept in the
algorithm-agnostic trait crate. Its existing capability names—such as
`ACMapBase`, `ACMapIter`, `ACMapAddAssign`, `ACMapInsert`, `ACMapRetain`, and
`ACMapConsume`—should also remain unless implementation work shows that two
capabilities should actually be merged or split.

A `SumStorage` is a new abstraction extracted from the fields currently owned
directly by `PauliSum`: its maps and reusable workspace. It is an actual value,
not a marker configuration:

```rust
pub trait SumStorage<W, C>: Clone
where
    W: Word + Indexable,
    C: Coefficient,
{
    type Map: ACMap<W, C>;

    fn data(&self) -> &Self::Map;
    fn data_mut(&mut self) -> &mut Self::Map;

    fn map_insert<F>(&mut self, f: F)
    where
        F: Fn(&W, &mut C) -> Option<(W, C)>;

    fn map_insert_multiple<F>(&mut self, f: F)
    where
        F: Fn(&W, &mut C) -> Option<Vec<(W, C)>>;

    fn map_add<F>(&mut self, f: F)
    where
        F: Fn(&W, &C) -> (W, C);

    fn consume(&mut self);
}
```

The exact closure bounds and support for multiple produced terms remain an
implementation detail for the prototype. The important boundary is that the
trait preserves the current semantic operation names without exposing physical
auxiliary maps or scratch buffers.

`SumStorage` owns the semantic whole-map operations and its reusable workspace.
It delegates to the lower-level `ACMap` batch kernels without restoring the
removed map-to-map `ACMapInsert::map_insert` method:

```text
SumStorage::map_insert           -> ACMapInsert::map_insert_vec
SumStorage::map_insert_multiple  -> ACMapInsert::map_insert_multiple
SumStorage::map_add              -> ACMapAddAssign::map_add_assign
```

This boundary is compatible with
`refactor/shrink-internal-trait-surface`: the higher-level sparse-sum operation
remains, while the dead low-level primitive stays removed.

The generalized engine may be named `OperatorSum`, but the Pauli specialization
retains the existing domain-facing `PauliSum` name. This is a new internal
generalization, not a requirement to rename Pauli call sites.

A classical tableau mixture follows the same principle and takes its
`EntryStore` directly rather than introducing a one-associated-type
`TableauMixtureAlgorithm` bundle:

```rust
pub struct GeneralizedTableauSum<C, T, S>
where
    C: Coefficient,
    T: Indexable,
    S: EntryStore<T, C>,
{
    entries: S,
}
```

`GeneralizedTableauSum` and `EntryStore` retain their current names because the
proposal does not change their underlying roles. `GeneralizedTableauSum` and
`OperatorSum` are both sparse linear combinations of indexable keys, so they
may eventually share an implementation. This iteration deliberately keeps them
separate: their mutation, branching, normalization, and storage requirements
have not yet been reduced to a proven common interface. The next design
iteration should look for the smallest useful common factor and merge only that
factor, rather than assuming that the two complete algorithms are identical.

Every keyed store must use the build hasher associated with its key type.

### Compatibility with current names

The redesign is not a vocabulary reset. The following names are retained or
changed according to whether their underlying responsibility changes:

| Current implementation | Proposal | Rationale |
| --- | --- | --- |
| `Config` | removed | The bundle itself is removed; this is not a rename. |
| `PauliWordTrait` | `Word` plus `Indexable` where used as a key | Word operations are generalized through `Word::Site`; hashing becomes a separate capability. |
| `n_qubits`, `get`, `set`, `weight` | `n_sites`, `get`, `set`, `weight` | Only the Pauli-specific extent name changes; the other operation names stay. |
| concrete `PauliWord` | `PauliWord` | The packed X/Z word is the same domain concept. |
| concrete `LossyPauliWord` | `LossyPauliWord` | The packed X/Z/loss representation remains concrete and flattened. |
| `PhasedPauliWord` | `PhasedPauliWord` alias over non-indexable `Phased` | The wrapper is generic over ordinary and lossy words but is not a production map key. |
| `rehash` | private cache invalidation | Recalculation changes from eager mutation-time work to lazy demand-time work without exposing cache mechanics through `Indexable`. |
| `Strategy` | `Policy` | Intentional terminology change requested for this redesign. |
| `ACMap` | `ACMap` | The associative coefficient map has the same role. |
| `PauliSum::data`, `map_insert`, `map_add` | same method names on `SumStorage` | These semantic operations already match the proposed boundary. |
| `PauliSum` map pair and `scratch` fields | `SumStorage` | A new abstraction is extracted from currently unnamed storage state. |
| `PauliSum` | `PauliSum` over generalized `OperatorSum` machinery | Pauli-facing code keeps its established name; `OperatorSum` names the new cross-algebra engine. |
| `GeneralizedTableauSum` | `GeneralizedTableauSum` | The classical mixture algorithm remains the same concept. |
| `EntryStore`, `VecStorage`, `MapStorage` | unchanged | The proposal uses the existing storage boundary and implementations. |
| `BuildHasher` | associated with `Indexable` | Hasher ownership moves from `Config` to each indexable key type. |
| `HashFinalize` | removed from the shared contract | Concrete keys may finalize or compose hashes privately. |
| `PauliStorage` | removed | Packed backing storage becomes private to the concrete word representation. |

Names such as `Word`, `Indexable`, `SumStorage`, and `OperatorSum` are
therefore new because they denote abstractions that do not
exist in the current implementation, not because the existing API is being
renamed wholesale.

### Trait admission rule

A proposed trait belongs in the design only when generic code consumes it and
there are multiple meaningful implementation families, or when it is an
established behavioral boundary implemented by different backends. A trait is
not justified merely to name inherent methods on one generic struct.

This rule keeps:

- `Indexable`, consumed by keyed stores and implemented by hash-enabled words
  and tableaus;
- `Word`, consumed by the shared sparse-sum engine and propagation algorithms;
- gate and noise traits, implemented across propagation and tableau backends;
- `SumStorage` and `ACMap`, implemented by genuinely different storage engines
  and collections;
- `EntryStore`, already implemented by `VecStorage` and `MapStorage`; and
- `Policy`, implemented by independent capacity and truncation behaviors.

It rejects the removed global `Config`, `PauliSumAlgorithm`,
`TableauMixtureAlgorithm`, and `TableauStorage` traits, as well as word
subtraits named `PauliWord`, `FermionWord`, or `LossyPauliWord`. Their
distinctions are expressed by `Word::Site` or by concrete types instead
of one-alphabet subtraits; the concrete `PauliWord` and `LossyPauliWord` type
names remain available.

### Sparse-sum branch staging

A propagation rule can turn one term into multiple terms. For example, a
Pauli rotation may produce:

```text
c P -> c cos(theta) P + c sin(theta) P'
```

The existing entry can be updated while the active map is traversed, but a new
key cannot generally be inserted into that same map during mutable iteration.
New terms must therefore be staged and merged after the traversal. New keys
may collide with each other or with existing keys, so the merge accumulates
their coefficients.

Only one staging mechanism is required for correctness. The current engine
uses two because they serve different performance paths:

- an auxiliary map supports whole-map rewrites, combines output collisions as
  they are produced, and retains its allocation across operations; and
- a reusable `Vec<(W, C)>` scratch buffer stages additional terms when existing
  entries remain in place, avoiding an auxiliary-map insertion followed by a
  second map probe during the merge.

Conceptually:

```text
whole-map rewrite:       active map -> auxiliary map -> swap
in-place branching:      active map + scratch buffer -> merge into active map
```

Neither physical mechanism is part of the public storage contract. A default
`SumStorage` implementation may privately contain both maps and the reusable
vector, matching the current `PauliSum` layout. A simpler backend may implement
both `map_insert` and `map_add` using only an auxiliary map; another backend
may use per-thread buffers. Whether retaining both mechanisms is worthwhile is
a benchmark decision, not a type-system requirement.

## Indexable values

Hash-enabled `Word` values and `Tableau` can both be expensive, mostly-stable
map keys. Their hashing contract should be expressed independently of any map.

The common capability is intentionally minimal:

```rust
pub trait Indexable: Clone + Eq + Hash {
    type BuildHasher: std::hash::BuildHasher + Clone + Default;
}
```

The important points are:

- `BuildHasher` retains the current associated-type name and moves from
  `Config` to the key type;
- the build hasher is associated with the key type, not with a configuration
  bundle;
- cache layout and invalidation are private representation invariants of the
  concrete type; and
- equality and hashing cover only structural key identity, never cache fields
  or incidental runtime state.

No generic consumer needs to name a cache type or request invalidation.
Structural mutation already occurs through `&mut self`, so each mutator can
clear the affected private cache as part of maintaining its concrete
invariants.

### Lazy hashing and interior mutability

Rust's `Hash::hash` receives `&self`, so shipped indexable words and tableaus
use private component `OnceLock<u64>` caches. `Hash::hash` may populate a cache
through shared access, while structural mutators clear affected cells through
their exclusive `&mut self` access. This preserves `Send + Sync` for the
shipped representations. `Indexable` itself does not require either
concurrency bound.

### Key mutation invariant

An indexable value must not be structurally mutated while it is stored as a map
key. Cache invalidation makes a value correct for its next insertion or lookup;
it cannot make in-place mutation of an existing map key valid.

Structural fields should therefore be private in the new representations.
Mutation must go through operations that invalidate the affected cache, or
through mutation guards that invalidate on completion.

## Concrete word hashing

Every concrete `Word` owns its `BuildHasher`, private cache representation,
structural hash algorithm, and invalidation logic. Pauli words hash
their X/Z content, lossy words compose Pauli and loss components, and future
fermion words hash their ordered factors. Factor order is part of fermionic
identity.

The trait layer does not expose packed storage to support hashing. Concrete
implementations hash their private fields and may apply hasher-specific
finalization internally. Detailed layouts and component invalidation rules are
in [`word-data-structures.md`](word-data-structures.md).

## Tableau indexability

A tableau may itself be used as a key by a classical-mixture algorithm, so the
concrete tableau implements `Indexable` directly and owns a tableau-specific
hasher and cache representation. This does not imply that a tableau is a
`Word`; they only share the `Indexable` key capability.

The tableau's structural hash is composed from its logical X/Z matrix, phase
plane, and per-qubit loss plane. It excludes RNG, padding, cache state, and
physical matrix orientation. Separate component caches allow phase-only
changes to avoid rehashing the X/Z matrix. Physical transposition is a layout
change, not a logical mutation, and does not invalidate the structural hash.

The concrete memory layout, component invalidation table, and canonical hash
order live in [`tableau-data-structure.md`](tableau-data-structure.md) so they
do not leak into the shared trait-system design.

## Expected generic composition

The intended composition is explicit:

```rust
OperatorSum<Coeff, Word, Storage, Policy>
Tableau
GeneralizedTableauSum<Coeff, TableauType, EntryStorage>
```

Domain-specific aliases or wrappers can preserve `PauliSum` and introduce
`FermionSum` without rebuilding a monolithic configuration trait.

## Non-goals for the first prototype

- Migrating the existing crates to `ppvm-traits-2` immediately.
- Merging `GeneralizedTableauSum` and `OperatorSum` in this iteration; only a
  smaller proven common factor should be considered in the next iteration.
- Defining one collection interface shared by all algorithms.
- Requiring every sparse-sum storage backend to physically contain both an
  auxiliary map and a scratch buffer.
- Exposing cache representation or invalidation through `Indexable`.
- Preserving `Copy` at the expense of correct lazy caching.
- Adding runtime dispatch for storage, hashing, or algorithm policies.

## Open design questions

1. Do benchmarks justify retaining both the auxiliary-map and vector-staging
   fast paths in the default sparse-sum storage backend?
