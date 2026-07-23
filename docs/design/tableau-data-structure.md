# Contiguous tableau data structure

Status: design sketch

## Purpose

The tableau should be a specialized bit-matrix data structure. It should not
be represented as `Vec<PhasedPauliWord<...>>`, and its internal rows should not
affect the shared PPVM trait system.

This design separates:

- the logical stabilizer-tableau model;
- its physical, contiguous memory layout;
- orientation changes used to accelerate different operations; and
- structural hashing used when a tableau is a classical-mixture key.

Gate, noise, measurement, reset, and `Indexable` traits observe tableau
behavior. They do not expose matrix blocks, row values, strides, alignment, or
orientation.

## Goals

- Store the X/Z tableau bits in contiguous, aligned, bit-packed memory.
- Make one- and two-qubit column operations efficient.
- Permit temporary transposition for row-oriented elimination and collapse.
- Keep phases and per-qubit loss state independently addressable and hashable.
- Hash the logical tableau independently of its current physical orientation.
- Avoid parameterizing the trait system by a row type or Pauli-word storage.
- Leave room for SIMD-width and padding changes without changing public
  behavioral traits.

## Non-goals

- Reusing the standalone `PauliWord` representation for tableau rows.
- Exposing borrowed tableau rows as the primary public interface.
- Standardizing a general-purpose matrix trait in `ppvm-traits-2`.
- Maintaining both row-major and column-major copies before benchmarks show
  that the memory and synchronization cost is worthwhile.
- Deciding in this document whether PPVM should store a forward or inverse
  tableau. Orientation and inversion are independent design choices.

## Logical model

For `n` qubits, a stabilizer/destabilizer tableau has `2n` generators. Each
generator contains an X bit and a Z bit for every qubit, plus a phase:

```text
                     qubit
                0   1   2   ... n-1
generator 0    x/z x/z x/z
generator 1    x/z x/z x/z
    ...
generator 2n-1 x/z x/z x/z

phase          one value per generator
```

The logical state consists of:

- an X matrix of shape `(2n, n)`;
- a Z matrix of shape `(2n, n)`;
- a phase plane of length `2n`; and
- a per-qubit loss plane of length `n`.

This model does not imply a physical row object.

## Physical storage

The first prototype should use one aligned contiguous allocation for the bit
planes, divided by computed offsets:

```rust
pub struct TableauData<Block> {
    blocks: AlignedVec<Block>,
    x_offset: usize,
    z_offset: usize,
    phase_offset: usize,
    loss_offset: usize,
    major_stride: usize,
    n_qubits: usize,
    orientation: Orientation,
}
```

`Block` is an internal implementation choice such as `u64` or a SIMD-width
block. It is not an associated type of a public tableau trait. Offsets and
strides account for alignment and padding, while logical accessors enforce the
actual `(2n, n)` dimensions.

Keeping all planes in one allocation improves cloning and locality and avoids
one allocation per generator. It also lets a mixture branch copy a single
contiguous region. If benchmarks favor separate aligned allocations for X and
Z, that change remains internal to `TableauData`.

Padding must either be kept zero or excluded from equality and hashing. Zeroed
padding is preferable because it permits bulk comparison and hashing of
canonical ranges.

## Loss ownership

The existing concrete `Tableau` always owns the per-qubit loss plane. This is a
capability of the same tableau type, not a `LossyTableau` variant or a
`Tableau<LossMode>` parameter:

```rust
pub struct Tableau {
    data: TableauData,
    lost_count: usize,
    // hash caches and RNG
}

pub struct GeneralizedTableau<C, I, S> {
    tableau: Tableau,
    coefficients: S,
    // threshold and measurement record
}
```

`GeneralizedTableau` therefore has no separate `is_lost: Vec<bool>`. Its gate,
noise, and measurement algorithms query and mutate the loss plane owned by the
inner tableau.

`lost_count` is derived metadata used to preserve a fast
`lost_count == 0` path. It is excluded from equality and hashing; debug builds
should verify that it equals the population count of the logical loss plane.
When no qubit is lost, gate kernels enter their existing lossless bulk path.
When loss is present, a one-qubit gate skips a lost target and a two-qubit gate
skips the operation if either target is lost, matching the current generalized
tableau semantics. Batch kernels should mask or skip lost targets without
allocating filtered target vectors.

This ownership enables a pure Clifford-plus-loss simulation with the same
`Tableau`. Loss events use the pure Clifford collapse procedure before marking
the affected qubit lost. Generalized loss events still use the
coefficient-aware generalized measurement procedure before marking the loss;
moving the bit does not move that algorithm. The pure path covers loss models
whose conditional trajectory remains a stabilizer state; faithful
non-Clifford survival back-action remains generalized.

## Column-major orientation

The default orientation should make the generator dimension contiguous for a
fixed qubit. In other words, the X and Z planes are column-major with respect
to the logical `(generator, qubit)` matrix:

```text
qubit 0: X bits for generators 0..2n, then padding
qubit 1: X bits for generators 0..2n, then padding
...

qubit 0: Z bits for generators 0..2n, then padding
qubit 1: Z bits for generators 0..2n, then padding
...
```

This layout makes a selected qubit column contiguous. Single-qubit gates can
load the X and Z columns, update them with bitwise operations, and update the
affected phase bits. Two-qubit gates operate on two pairs of contiguous
columns. Measurement can scan the selected anticommutation column without
stepping across separately allocated row objects.

Stim uses aligned SIMD bit tables and documents its tableau layout as
column-major, with output-observable iteration following contiguous memory. It
also provides an explicit quadrant transpose and a transposition guard for
operations that need the opposite orientation:

- [Stim `Tableau`](https://github.com/quantumlib/Stim/blob/main/src/stim/stabilizers/tableau.h)
- [Stim `TableauSimulator`](https://github.com/quantumlib/Stim/blob/main/src/stim/simulators/tableau_simulator.h)
- [Stim `simd_bit_table`](https://github.com/quantumlib/Stim/blob/main/src/stim/mem/simd_bit_table.h)

Column-major storage is only one part of Stim's performance strategy. PPVM
should benchmark its own gate, measurement, and sampling workloads instead of
assuming the same total performance from layout alone.

## Temporary transposition

Row multiplication and elimination need long generator rows and therefore
prefer the opposite orientation. The initial design should transpose the X/Z
quadrants temporarily instead of storing two permanent copies:

```rust
pub enum Orientation {
    ColumnMajor,
    RowMajor,
}

pub struct TransposedTableau<'a> {
    tableau: &'a mut Tableau,
}
```

Creating the guard transposes the bit matrices and marks the physical
orientation. Dropping it restores the canonical column-major orientation:

```rust
impl Drop for TransposedTableau<'_> {
    fn drop(&mut self) {
        self.tableau.restore_column_major();
    }
}
```

Operations that require row-major access receive the guard, making the
orientation precondition explicit. Public methods return with the tableau in
canonical orientation. Panic safety must be preserved by the guard's `Drop`
implementation.

Transposition is a physical reordering of the same logical bits. It does not
invalidate structural hashes or change equality. Hashing should occur only
through logical access or while the tableau is in its public canonical state.

Maintaining both orientations with dirty flags remains a future alternative
for workloads that switch frequently enough to amortize the doubled storage.

## Tableau API boundary

The public API exposes logical operations instead of storage slices:

```rust
impl Tableau {
    pub fn n_qubits(&self) -> usize;

    pub fn h(&mut self, qubit: usize);
    pub fn cnot(&mut self, control: usize, target: usize);
    pub fn measure_z(&mut self, qubit: usize) -> Option<bool>;

    pub fn x_bit(&self, generator: usize, qubit: usize) -> bool;
    pub fn z_bit(&self, generator: usize, qubit: usize) -> bool;
    pub fn phase(&self, generator: usize) -> Phase;
    pub fn is_lost(&self, qubit: usize) -> bool;
}
```

Logical bit accessors are useful for tests, serialization, debugging, and
interoperability. They are not intended as the hot gate implementation path.
Bulk import/export methods may use a canonical serialized representation
without exposing the in-memory layout.

There should be no `stabilizers_mut() -> &mut [Row]` or equivalent escape hatch
that bypasses hash invalidation and orientation invariants. Specialized
internal row operations operate through `TableauData` or a transposition guard.

## Measurement algorithms

The Rust behavioral boundary is one loss-aware trait:

```rust
pub trait Measure {
    fn measure(&mut self, qubit: usize) -> Option<bool>;

    fn measure_many(&mut self, targets: &[usize]) -> Vec<Option<bool>> {
        targets.iter().map(|&q| self.measure(q)).collect()
    }
}
```

`Some(false)` and `Some(true)` represent computational-basis outcomes; `None`
represents a lost qubit. This is the existing core representation used by
`GeneralizedTableau`. The Python binding may continue mapping it to
`MeasurementResult::{ZERO, ONE, LOST}`. The old public split between
`Measure -> bool` and `LossyMeasure -> Option<bool>` is removed. The bare
boolean Clifford measurement routine becomes a private helper called only
after the public implementation has established that the target is present.

The common trait and result type do not imply a common measurement algorithm:

- `Tableau::measure` checks the loss plane and then uses the pure Clifford
  stabilizer measurement procedure. It does not decompose against a sparse
  generalized-state coefficient basis.
- `GeneralizedTableau::measure` checks the inner tableau's loss plane, then
  decomposes the measured Pauli into stabilizers and destabilizers and updates
  the sparse coefficients. This path is fundamentally \(O(n^2)\).

The physical tableau must make both stabilizer and destabilizer generators
available to the generalized decomposition, but that requirement must not be
promoted into the `Measure` trait or force the generalized algorithm onto the
pure Clifford implementation. Pure and generalized measurement performance
must be benchmarked separately; Stim's inverse-tableau measurement
optimizations are not complexity promises for the generalized algorithm.

## Gate access patterns

The implementation should categorize mutations by physical access and hash
effect:

| Operation | Preferred access | X/Z changed | Phase changed |
| --- | --- | --- | --- |
| Pauli `X`, `Y`, `Z` | column | no | possibly |
| `H`, `S` | one column pair | yes | possibly |
| `CNOT`, `CZ` | two column pairs | yes | possibly |
| Find measurement pivot | column | no | no |
| Row multiplication | transposed row | yes | possibly |
| Collapse/elimination | column scan + transposed rows | yes | yes |
| Physical transpose | bulk matrix | logically no | no |

This table describes logical mutations. A gate may determine that no phase bit
changed and preserve the phase cache in that special case, but conservative
component invalidation is correct for the first implementation.

## Structural hashing

A tableau used in a classical mixture is an `Indexable` key. It owns its own
hasher and cache representations; neither is inherited from `PauliWord`.

The structural hash is composed from independent logical components:

```text
tableau hash = combine(xz hash, phase hash, loss hash)
```

The X and Z planes share an `xz_hash` cache because most Clifford mutations
update them together. The phase plane has a separate cache so Pauli
conjugations and sign changes do not force a matrix rehash. The independently
mutable loss plane uses a third cache.

```rust
pub struct Tableau {
    data: TableauData,
    lost_count: usize,
    xz_hash: OnceLock<u64>,
    phase_hash: OnceLock<u64>,
    loss_hash: OnceLock<u64>,
    rng: SmallRng,
}
```

The cache fields are private representation choices made by the tableau
author. `Indexable` exposes only the associated build hasher; it does not name
cache types or expose invalidation.

Equality and hashing include logical qubit count, generator order, all logical
X/Z bits, phases, and loss state. They exclude:

- RNG state;
- allocation capacity;
- alignment padding;
- cache values and validity flags; and
- current physical orientation.

The component invalidation rules are:

| Mutation | X/Z cache | Phase cache | Loss cache |
| --- | --- | --- | --- |
| Pauli `X`, `Y`, `Z` | preserve | invalidate if changed | preserve |
| Direct phase change | preserve | invalidate | preserve |
| `H`, `S`, `CNOT`, `CZ` | invalidate | invalidate if changed | preserve |
| Row multiplication | invalidate | invalidate if changed | preserve |
| Toggle a loss bit | preserve | preserve | invalidate |
| Physical transpose | preserve | preserve | preserve |
| RNG update | preserve | preserve | preserve |

The current `ppvm-tableau-sum` split between `word_fingerprint` and
`phase_loss_hash` is evidence that component hashing matters. The new tableau
owns X/Z, phase, and loss components directly.

## Cloning and mixture use

Classical-mixture branching can clone tableaus frequently. Contiguous backing
storage makes cloning a bulk memory copy. A clone may copy valid hash caches
because it initially has identical logical contents. Subsequent mutation of
the branch invalidates only the affected components.

Copy-on-write backing storage may be evaluated later, but it should not be part
of the initial design: most branches are mutated immediately, which may turn
reference counting and deferred copying into overhead.

An indexable tableau must not be structurally mutated while stored as a map
key. Mixture storage removes or clones a key before applying gates and inserts
the resulting tableau under its updated hash.

## Sampling implications

Column-major X/Z planes make fixed-qubit queries and bit-parallel gate updates
efficient. They are also compatible with scanning many generators during a
measurement. Temporary transposition makes elimination and row products
contiguous when required.

This layout should be evaluated separately from higher-level sampling
algorithms. Stim also uses an inverse tableau and reference-frame sampling;
those algorithmic choices are not consequences of column-major storage and do
not belong in the PPVM trait system.

## Prototype validation

The prototype should include:

- property tests comparing gates and measurements with the existing tableau;
- differential tests for the pure Clifford and generalized measurement
  algorithms, including lost targets;
- tests that gates skip lost targets and retain the `lost_count == 0` fast
  path;
- round-trip tests for column-major -> row-major -> column-major transpose;
- equality and hash tests across physical orientations;
- tests proving phase-only changes preserve the X/Z hash cache;
- tests proving padding never affects equality or hashing;
- benchmarks for one- and two-qubit Clifford gates;
- separate benchmarks for pure Clifford and generalized deterministic and
  random measurement paths;
- benchmarks for lossless gates, sparse loss, and Clifford-plus-loss sampling;
- benchmarks for clone-and-mutate mixture branching; and
- benchmarks comparing permanent row-major, permanent column-major, and
  temporary-transpose variants on representative circuits.

## Open questions

1. What block width and alignment should the first implementation use?
2. Should phases occupy one or two bits per generator in the tableau model?
3. Which operations should receive a transposition guard versus performing
   column-strided work directly?
4. Does a dual-orientation representation outperform temporary transposition
   for PPVM's measurement-heavy workloads?
5. Should PPVM ultimately store a forward or inverse tableau?
