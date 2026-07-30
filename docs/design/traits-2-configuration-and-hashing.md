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

### Coefficient, angle, and truncation

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

The decomposition must not stop at `Config`; the current `Coefficient` trait is
itself a bundle. Its own documentation states that it "bundles every arithmetic
operation," and it carries two responsibilities that are not value-domain
arithmetic: `sin_cos` (a rotation concern) and `cutoff` (a truncation concern —
the very policy this redesign extracts). `Coefficient` is therefore split into a
value domain, a separate angle domain, and a policy predicate.

The value domain keeps only arithmetic and a magnitude, replacing `cutoff` with
a property that a policy can threshold:

```rust
/// Value-domain arithmetic only — no rotation, no truncation.
pub trait Coefficient:
    PartialEq + Clone + num::Zero
    + Neg<Output = Self>
    + Add<Self, Output = Self> + Sub<Self, Output = Self>
    + Mul<Self, Output = Self> + Mul<f64, Output = Self>
    + AddAssign<Self> + MulAssign<Self>
    + std::iter::Sum + Send + Sync
{
    fn mul_sign(&self, sign: i8) -> Self;
    fn half(&self) -> Self;

    /// Nonnegative magnitude. Exposes a property of the value for a `Policy` to
    /// threshold; it does not itself decide any cutoff. Replaces `cutoff`.
    fn magnitude(&self) -> f64;
}
```

The rotation angle becomes a separate domain, so it is no longer welded to the
coefficient. The angle type defaults to the coefficient, recovering today's
`rx(theta: C)` for the common case while permitting an `f64`-coefficient sum
driven by a symbolic or parametric angle:

```rust
/// A rotation angle that yields `(sin, cos)` in coefficient domain `C`.
pub trait Angle<C: Coefficient> {
    fn sin_cos(&self) -> (C, C);
}

impl Angle<f64> for f64 {
    fn sin_cos(&self) -> (f64, f64) { num::traits::Float::sin_cos(*self) }
}
```

The rotation traits (see [Behavioral traits](#behavioral-traits)) take the angle
domain as a defaulted parameter. Whether the "concrete coefficient, symbolic
angle" combination is supported or deliberately forbidden (via a `where A = C`
bound) is an explicit decision recorded in the open questions, not a silent
property of one fused trait.

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

pub trait RotationOne<C: Coefficient, A: Angle<C> = C> {
    fn rx(&mut self, qubit: usize, theta: A);
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

The choices that are *intrinsic to a storage instance* — which word it keys on
and which coefficient it accumulates — are associated types of that storage
rather than free parameters, the same way `HashMap<K, V>` fixes `K` and `V` per
concrete map. `SumStorage` (defined below) therefore exposes `type Word` and
`type Coeff`, and the reusable sparse-sum shape needs only its storage and its
policy:

```rust
pub struct Sum<S, P = NoPolicy>
where
    S: SumStorage,
    P: Policy<S::Word, S::Coeff>,
{
    storage: S,
    policy: P,
    /// Invariant: every key `k` in `storage` satisfies `k.n_sites() ==
    /// n_sites`. An empty sum has no word to derive the width from, so the
    /// field is carried explicitly and checked by a `debug_assert!` on every
    /// insertion path.
    n_sites: usize,
}
```

`Word` is no longer a free parameter with no home; it is `S::Word`, so there is
no phantom axis and no way to pair a storage type with the wrong word. This also
resolves the tension with the earlier note that associated types are not
justified merely to shrink a parameter count: `Word` and `Coeff` are genuinely
intrinsic to a storage instance, so making them associated types is proper
encapsulation, not a false bundle. Propagation methods still select their
algebra through the site type by bounding `S::Word: Word<Site = Pauli>` or
`Word<Site = FermionSite>`.

#### The convenience bundle

Collapsing to `Sum<S, P>` restores ergonomics through a type-alias family, not a
foundational trait. Each alias fixes a domain's axes and defaults the common
knobs, so a one-token `PauliSum` name returns:

```rust
pub type PauliSum<C = f64, P = NoPolicy>      = Sum<HashStore<PauliWord, C>, P>;
pub type LossyPauliSum<C = f64, P = NoPolicy> = Sum<HashStore<LossyPauliWord, C>, P>;
pub type FermionSum<C = f64, P = NoPolicy>    = Sum<HashStore<FermionWord, C>, P>;
```

A user who wants to pin a reusable configuration writes one `type` line instead
of an `impl Config`:

```rust
type MyPauliSum = Sum<DashStore<PauliWord, Complex<f64>>, MaxPauliWeight>;
```

This is deliberately a bundle, and it is safe to be one because it commits
neither of `Config`'s two sins:

- **It is not load-bearing at the trait level.** The shared behavioral traits
  (`Clifford`, `RotationOne`, `PauliError`, `Measure`) are implemented on `Sum`
  and take `S::Coeff` directly. No gate or noise trait is generic over the
  alias, so an algorithm that needs only a coefficient never names storage or
  policy.
- **Its axes stay independently overridable.** Because it is a plain alias over
  defaulted parameters, changing one axis is `Sum<OtherStore, _>`, not a new
  trait impl. `Config`'s failure was fusing orthogonal axes into one impl you
  had to rewrite wholesale to vary a single choice; a defaulted alias has the
  opposite property.

The name `Sum` shares a bare spelling with the `std::iter::Sum` *trait* that
`Coefficient` bounds on. They coexist because one is a type and one is a trait,
and the bound is written fully qualified; `PauliSum` and `FermionSum` remain the
public domain-facing names in any case.

`Policy` is the proposed name for the current `Strategy` concept. It retains
the current responsibilities: predicting initial capacity and truncating the
sum. Existing concrete strategies become policies without otherwise changing
their meaning; `NoStrategy` and `CombinedStrategy` become `NoPolicy` and
`CombinedPolicy`, while `MaxPauliWeight` and `CoefficientThreshold` keep their
established names:

```rust
pub trait Policy<W, C>: Default + Clone
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

`Policy` does not require `Copy`. The bound bought nothing — policies are used
through `&self` — and would forbid a stateful policy such as a per-region
threshold vector or an adaptive budget; requiring it would also contradict this
proposal's own non-goal of preserving `Copy` at the expense of correctness.

Truncation policy also reclaims the cutoff decision from `Coefficient`. The
value type now exposes only `magnitude()`, a property of the value; the policy
owns the comparison:

```rust
#[derive(Default, Clone)]
pub struct CoefficientThreshold {
    pub threshold: f64,
}

impl<W: Word + Indexable, C: Coefficient> Policy<W, C> for CoefficientThreshold {
    fn capacity(&self, _n_sites: usize) -> usize {
        0
    }

    fn truncate<M: ACMap<W, C>>(&self, map: &mut M) {
        map.retain(|_word, coeff| coeff.magnitude() >= self.threshold);
    }
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
pub trait SumStorage: Clone {
    type Word: Word + Indexable;
    type Coeff: Coefficient;
    type Map: ACMap<Self::Word, Self::Coeff>;

    fn data(&self) -> &Self::Map;
    fn data_mut(&mut self) -> &mut Self::Map;

    fn map_insert<F>(&mut self, f: F)
    where
        F: Fn(&Self::Word, &mut Self::Coeff) -> Option<(Self::Word, Self::Coeff)>;

    fn map_insert_multiple<F>(&mut self, f: F)
    where
        F: Fn(&Self::Word, &mut Self::Coeff) -> Option<Vec<(Self::Word, Self::Coeff)>>;

    fn map_add<F>(&mut self, f: F)
    where
        F: Fn(&Self::Word, &Self::Coeff) -> (Self::Word, Self::Coeff);

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

Those `ACMap` kernels are the scalar spelling of a bulk operation. The
batch-first, structure-of-arrays layout the kernels consume on hot paths — and
that admits prefetched, SIMD, threaded, and offloaded backends — is specified in
[Batch execution and the hash-join contract](#batch-execution-and-the-hash-join-contract)
below.

The generalized engine is named `Sum`, but the Pauli specialization retains the
existing domain-facing `PauliSum` name (a defaulted type alias over `Sum`). This
is a new internal generalization, not a requirement to rename Pauli call sites.

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
`Sum` are both sparse linear combinations of indexable keys, so they
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
| `Strategy` | `Policy` | Intentional terminology change requested for this redesign; the `Copy` bound is dropped. |
| `Coefficient::cutoff` | removed | The truncation predicate moves to `Policy`; `Coefficient::magnitude` exposes only the value property the policy thresholds. |
| `Coefficient::sin_cos` | `Angle<C>` | The rotation angle becomes a separate domain, defaulting to the coefficient type. |
| `ACMap` | `ACMap` | The associative coefficient map has the same role. |
| `PauliSum::data`, `map_insert`, `map_add` | same method names on `SumStorage` | These semantic operations already match the proposed boundary. |
| `PauliSum` map pair and `scratch` fields | `SumStorage` | A new abstraction is extracted from currently unnamed storage state. |
| `PauliSum` | `PauliSum` over generalized `Sum` machinery | Pauli-facing code keeps its established name; `Sum` names the new cross-algebra engine. |
| `GeneralizedTableauSum` | `GeneralizedTableauSum` | The classical mixture algorithm remains the same concept. |
| `EntryStore`, `VecStorage`, `MapStorage` | unchanged | The proposal uses the existing storage boundary and implementations. |
| `BuildHasher` | associated with `Indexable` | Hasher ownership moves from `Config` to each indexable key type. |
| `HashFinalize` | removed from the shared contract | Concrete keys may finalize or compose hashes privately. |
| `PauliStorage` | removed | Packed backing storage becomes private to the concrete word representation. |
| (new) | `Columnar`, `KeyColumn`, `KeyBatch`, `TermBatch`, `ACMapBatch` | Structure-of-arrays batch/hash-join contract, kept as a separate trait so `Indexable` stays minimal; see [Batch execution and the hash-join contract](#batch-execution-and-the-hash-join-contract). |

Names such as `Word`, `Indexable`, `SumStorage`, and `Sum` are
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

## Batch execution and the hash-join contract

The refactor's stated goal is a map contract friendly to CPU SIMD, memory
prefetching, multi-threading, and GPUs. None of those are properties of a
particular collection; they are properties of the *shape of work the collection
is asked to do*. This section fixes that shape — bulk, and hash-join-shaped —
and leaves every implementation choice open, matching the proposal's "contract
now, benchmark the backend later" stance. It keeps `Indexable` minimal: the
structure-of-arrays capability is a separate trait, not another associated type
on the key's hashing contract.

### The merge is a hash join

The branch-staging merge described above — probe each produced term against the
current sum, accumulate its coefficient on a match, insert it on a miss — is a
*hash join with aggregation* (equivalently a group-by aggregate). The current
operator map is the build side; the terms a gate produces are the probe side;
coefficient accumulation is the aggregate. It is the hottest operation in
propagation.

Naming it explicitly matters because its performance is governed by the same
cost model as a database hash join: at the sizes of interest the build table
does not fit in cache, so each probe is a random access that misses to main
memory, and the probe phase is bound by memory latency rather than by
arithmetic.

The `SumStorage` closures (`map_insert`, `map_add`) and the `ACMap` kernels they
delegate to still describe *what* the merge computes. But a closure applied to
one `(key, coeff)` at a time exposes no batch to prefetch, no homogeneous run to
vectorize, no partition to distribute across threads, and no bulk kernel to hand
to a device. Whatever the eventual implementation, the *layout* the kernels
consume must be a batch, or those backends cannot be written against it.

### Prefetching the probe phase

The canonical treatment of the memory-latency problem is Chen, Ailamaki,
Gibbons, and Mowry, *Improving Hash Join Performance through Prefetching* (ICDE
2004; extended in ACM TODS 32(3), 2007). Its insight is that one probe cannot be
made faster — the chain `hash -> load bucket -> compare -> load payload` is
fully data-dependent — but many *independent* probes can be overlapped. The
paper gives two schedules:

- **group prefetching**: take a group of probe keys, issue the bucket prefetch
  for all of them, then do the comparisons once the lines have arrived; and
- **software-pipelined prefetching**: stagger the stages of consecutive groups
  so the prefetches of later keys are in flight while earlier keys are compared.

Both turn latency-bound probing into throughput-bound probing by keeping many
cache misses outstanding at once, up to the hardware's fill-buffer limit.
Neither is expressible through a scalar lookup; both need the whole group of
keys up front. This is the concrete reason the layout is batch-first.

The private structural-hash cache (the shipped words' `OnceLock<u64>`
components) is a prerequisite: because each key already carries its computed
hash, a bulk `hash_into` mostly gathers cached values, so the group-prefetch
loop reduces to "prefetch the bucket now, confirm later" with no hashing on the
critical path.

### The batch contract

Work is presented as a *structure-of-arrays* batch of terms, with the execution
strategy left entirely to the implementation. A batch is columns split along the
two phases of the join, so probe and aggregation touch disjoint memory:

```rust
/// Keys plus their precomputed structural hashes, in parallel columns. The
/// probe side of the join; it carries no coefficients, so a probe streams only
/// key and hash memory.
pub struct KeyBatch<W: Columnar> {
    keys: W::Column,   // structure-of-arrays key planes, owned by the word type
    hashes: Vec<u64>,  // one structural hash per key, parallel to `keys`
}

/// A `KeyBatch` with the coefficient column attached: the produced terms
/// awaiting merge. Coefficients are a separate column, touched only when a
/// probe resolves to an aggregate.
pub struct TermBatch<W: Columnar, C> {
    keys: KeyBatch<W>,
    coeffs: Vec<C>,
}
```

The key column is itself structure-of-arrays and is owned by the concrete key
type, because only that type knows its planes (a packed Pauli word splits into
an X-bit plane block and a Z-bit plane block; the flattened `LossyPauliWord`
adds a loss plane). That capability is a **separate trait**, so `Indexable`
stays minimal and no column type leaks through the hashing contract:

```rust
/// A key that can be laid out as a structure-of-arrays column. Separate from
/// `Indexable` so the minimal hashing contract is unchanged: a batched key is
/// both `Indexable` (a valid map key) and `Columnar` (has a column layout).
pub trait Columnar: Indexable {
    type Column: KeyColumn<Key = Self>;
}

pub trait KeyColumn: Default + Clone {
    type Key: Columnar;

    fn len(&self) -> usize;
    fn with_capacity(n: usize) -> Self;

    /// Append one produced key; the column keeps each plane contiguous.
    fn push(&mut self, key: Self::Key);

    /// Bulk structural hash of the whole column into a parallel hash column.
    /// Where per-plane SIMD hashing lives, and what feeds the group-prefetch
    /// loop; it must agree bit for bit with the scalar `Hash` of each key.
    fn hash_into(&self, out: &mut [u64]);

    /// Join confirm: compare element `i` against a build-side key after a hash
    /// or tag match, without materializing the whole element.
    fn key_eq(&self, i: usize, other: &Self::Key) -> bool;

    /// Select or permute elements into a new column: the primitive for radix
    /// partitioning across threads, compaction during truncation, and staging a
    /// sub-batch for a device. Operates plane by plane, never scalar.
    fn gather(&self, indices: &[u32]) -> Self;

    /// Scalar materialization of one element — a naive backend's fallback,
    /// never the hot path.
    fn get(&self, i: usize) -> Self::Key;
}
```

The bulk map operations then consume columns. They are the batch spelling of the
existing `ACMap` kernels, added alongside them rather than replacing the
semantic `SumStorage` methods:

```rust
pub trait ACMapBatch<W, C>
where
    W: Word + Columnar,
    C: Coefficient,
{
    /// Merge a batch into the sum: for each term, accumulate onto an existing
    /// key or insert a new one. The build/probe-with-group-by of a hash join.
    /// The implementation owns whether it runs scalar, group- or
    /// pipeline-prefetched, SIMD-vectorized, hash-partitioned across threads,
    /// or offloaded; the contract fixes only the result and the layout.
    fn upsert_batch(&mut self, batch: &TermBatch<W, C>);

    /// Read-only probe of a key column, for the overlap and expectation paths.
    /// Takes a `KeyBatch` so the precomputed hash column drives the prefetch
    /// with no coefficient column in the working set.
    fn probe_batch(&self, keys: &KeyBatch<W>, out: &mut [Option<C>]);
}
```

Three points make this layout the performant one rather than a cosmetic
reshuffle:

- **Column separation follows the join phases.** The probe touches keys and
  hashes; aggregation touches coefficients. Separate columns keep the probe's
  bandwidth and cache footprint minimal and let a backend place the coefficient
  column on a different device or NUMA node from the keys.
- **The hash column is precomputed and contiguous.** The group-prefetch kernel
  streams the hash column, prefetches the matching buckets, and only then
  confirms against the key column — the private hash cache above keeps that
  column cheap to fill.
- **Plane-level structure-of-arrays inside the key column** lets a backend hash
  and compare the same machine-word lane across many keys with vertical SIMD,
  and gives a GPU coalesced loads. The concrete plane layout, alignment, and
  padding live in [`word-data-structures.md`](word-data-structures.md); they are
  never visible through `KeyColumn`.

`upsert_batch` is where the branch-staging merge lands: the reusable scratch
buffer of the previous section becomes a `TermBatch` — the probe side of the
join — and the `SumStorage` merge desugars to it. The scalar `SumStorage`
closures remain the semantic surface for a naive backend, but generic hot paths
build a `TermBatch` and hand it to the batch kernel.

Gate and rotation traits change accordingly, and this is the interface change
users feel first: term *production* is separated from term *insertion*. A
rotation appends produced terms into a `TermBatch` — filling the key column and
the coefficient column as it goes — instead of mutating the map through a
callback. That decoupling, plus the columnar layout, is what lets the produced
batch be prefetched, vectorized, partitioned, or shipped to a device before it
ever touches the table, and it keeps the propagation rule (the physics)
independent of the merge strategy (the systems concern).

### What the batch contract asks of the other contracts

- `Indexable` is unchanged and stays minimal. The batch surface lives on the
  separate `Columnar` trait, which a key type implements only when it is used in
  a batch. A packed Pauli word's column is its X and Z plane blocks and the
  flattened `LossyPauliWord` adds a loss plane; `Phased<W>` is not indexable, so
  it is not `Columnar` and never appears in a batch.
- `Word`: unaffected — `set` and friends still mutate one value; columns are
  built by appending produced keys, not by mutating in place.
- `Policy`: truncation already operates on the whole map and is naturally bulk;
  it should be expressible as a batch retain (via `KeyColumn::gather`) so it
  composes with a partitioned or offloaded table.
- `SumStorage`: its staging fields become the join's probe-side buffer and its
  merge desugars to `upsert_batch`. Whether it keeps both an auxiliary map and a
  scratch buffer stays a backend decision, unchanged from above.

### Parallel and offloaded backends

Recognizing the operation as a hash join fixes the parallel and GPU stories
without further interface changes. A partitioned hash join radix-partitions both
sides by the high bits of the key hash so each partition is a disjoint
sub-table; partitions then merge independently — one per thread, or one per GPU
block — with no cross-partition synchronization. The batch contract is precisely
what a partitioner consumes: `upsert_batch` may internally `gather` a
`TermBatch` into per-partition batches and run them concurrently, and a device
backend may copy a `TermBatch` across the host/device boundary and run the probe
as a kernel. None of this is visible in the contract; all of it is precluded by
a scalar one. This is the partition-then-merge shape the `DashMap` backend
already gestures at, now stated as the contract that makes such backends
interchangeable rather than separate code paths.

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
Sum<Storage, Policy>   // e.g. PauliSum<f64> = Sum<HashStore<PauliWord, f64>>
Tableau
GeneralizedTableauSum<Coeff, TableauType, EntryStorage>
```

Domain-specific aliases or wrappers can preserve `PauliSum` and introduce
`FermionSum` without rebuilding a monolithic configuration trait.

## Non-goals for the first prototype

- Migrating the existing crates to `ppvm-traits-2` immediately.
- Merging `GeneralizedTableauSum` and `Sum` in this iteration; only a
  smaller proven common factor should be considered in the next iteration.
- Defining one collection interface shared by all algorithms.
- Requiring every sparse-sum storage backend to physically contain both an
  auxiliary map and a scratch buffer.
- Exposing cache representation or invalidation through `Indexable`.
- Preserving `Copy` at the expense of correct lazy caching.
- Adding runtime dispatch for storage, hashing, or algorithm policies.
- Committing to a single batch-execution backend. The hash-join contract admits
  scalar, prefetched, SIMD, thread-partitioned, and device backends; which one
  ships is a later benchmark decision, not a type-system decision.
- Implementing the prefetch or partition schedule itself in this iteration. The
  batch contract is defined here; the group/pipeline schedule and partitioner
  are implementation work behind it.

## Open design questions

1. Do benchmarks justify retaining both the auxiliary-map and vector-staging
   fast paths in the default sparse-sum storage backend?
2. At what granularity does the key column split into planes, and how are those
   planes padded and aligned — to the widest supported SIMD lane, to a cache
   line, or to a device's coalescing width? The contract fixes that the column
   is structure-of-arrays; the plane granularity and alignment are open, and may
   need to vary by target.
3. What batch (group) size does each backend use, and should the contract expose
   it at all, or keep it entirely internal to the implementation?
4. Should `probe_batch` yield coefficients, or opaque slot handles a later read
   resolves, so a caller can probe once and then mutate the located slot?

## References

- Shimin Chen, Anastassia Ailamaki, Phillip B. Gibbons, and Todd C. Mowry.
  "Improving Hash Join Performance through Prefetching." *Proceedings of the
  20th International Conference on Data Engineering (ICDE)*, 2004. Extended
  version: *ACM Transactions on Database Systems* 32(3), Article 17, 2007. The
  origin of group prefetching and software-pipelined prefetching for the probe
  phase, and the cost model the batch contract is built to serve.
