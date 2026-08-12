# Contiguous tableau data structure

Status: **partly implemented** in `ppvm-tableau-2` as of 2026-08-11. The
contiguous column-major storage, the transposition guard, the measurement route
and the inverse tableau are shipped (`crates/ppvm-tableau-2/src/storage/`,
`src/inverse.rs`); the loss-plane cutover is not. Sections below are marked
where the shipped code decided a question this sketch left open, or contradicted
it.

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

Orientation (column- vs row-major) and inversion (forward vs inverse) are
independent choices; this document fixes both — column-major storage of the
**inverse** tableau (see [Inverse tableau](#inverse-tableau)).

**Correction (2026-08-11).** Those two choices are not independent, and the
sketch's framing of them as a *choice* is the error: forward columns are inverse
rows, so column-major storage of the forward frame **is** row-major storage of
the inverse. The frame ships as the forward tableau and the inverse is read out
of the same arena, with only its `2n` signs separately tracked.

The stated reason for the pairing — that it lets gates and measurement share one
orientation, so no measurement transposes — was also wrong, and Stim's own source
is the counterexample: `collapse_qubit_z` runs under a `TableauTransposedRaii`
([`tableau_simulator.inl:96`, `:1341`](https://github.com/quantumlib/Stim/blob/main/src/stim/simulators/tableau_simulator.inl)),
which transposes on construction and again on destruction. What the inverse
removes is the *search* (the determinism check becomes one sign bit) and the
`ℤ/4` fold in `compute_decomposition` (one sign bit plus a popcount);
elimination still multiplies whole generators, so the guard stays — see
[Temporary transposition](#temporary-transposition), where a row-oriented pass
either gathers the `k` generators it needs out of the column-major arena or takes
the guard, whichever is cheaper.

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

**Shipped (2026-08-11).** `u64` blocks in a 32-byte-aligned allocation, and the
block type is private rather than a generic parameter, so a 256-bit SIMD block
is a later drop-in that changes no signature (open question 1, resolved). Rust
has no stable portable SIMD, so widening now would cost a dependency or a
nightly feature for no measured gain: the gate kernels are already
memory-bound at the widths that matter. The four X/Z quadrants are stored
**square** (`n × n` each, not `2n × n` per plane) because an in-place blockwise
transpose needs square quadrants — the same reason Stim's
`do_transpose_quadrants` is written over quadrants rather than the whole table.

Keeping all planes in one allocation improves cloning and locality and avoids
one allocation per generator. It also lets a mixture branch copy a single
contiguous region. If benchmarks favor separate aligned allocations for X and
Z, that change remains internal to `TableauData`.

Padding must either be kept zero or excluded from equality and hashing. Zeroed
padding is preferable because it permits bulk comparison and hashing of
canonical ranges.

## Loss ownership

**Superseded (2026-08-11).** The section below says loss is a capability of the
one `Tableau` type and explicitly not a `LossyTableau` variant. That is
overruled: loss becomes a **tower**, `Tableau` (lossless) → `LossyTableau`
(`Tableau` plus a packed loss plane) → `GeneralizedTableau` generic over the
frame. This mirrors `ppvm-pauli-word-2`'s `PauliWord` / `LossyPauliWord` split,
where `PauliBits::is_lost` is a const-`false` default that the lossy word
overrides, so one loss-aware kernel serves both types and the lossless build
folds every `if is_lost(q)` branch away at monomorphization. Making
`GeneralizedTableau` frame-generic is what carries that property up: a lossless
simulation stops paying for loss checks it cannot need, which the always-owned
plane below cannot express.

Everything the section says about *representation* still holds — the loss plane
is bit-packed, never transposed, excluded from the X/Z digest, and carries a
`lost_count` fast path. Only its ownership moves. The shipped `TableauData`
already allocates the plane and its accessors, marked reserved; the live flags
are still `GeneralizedTableau::is_lost: Vec<bool>` until the tower lands.

The existing concrete `Tableau` always owns the per-qubit loss plane. This is a
capability of the same tableau type, not a `LossyTableau` variant or a
`Tableau<LossMode>` parameter:

```rust
pub struct Tableau {
    data: TableauData,   // the *inverse* tableau (see below), column-major
    lost_count: usize,
    // hash caches only — NO rng field; see "Randomness is injected"
}

pub struct GeneralizedTableau<C = f64, P = CoefficientThreshold> {
    frame: Tableau,                                    // the Clifford frame U; owns loss
    amplitudes: Sum<Vec<(Bitstring, Complex<C>)>, P>,  // C[bitstring] — the graded algebra
    // measurement record
}
```

`GeneralizedTableau` therefore has no separate `is_lost: Vec<bool>`. Its gate,
noise, and measurement algorithms query and mutate the loss plane owned by the
inner `frame` tableau. Its `amplitudes` field is the graded-algebra `Sum` over
bitstring keys (see
[A third instantiation](traits-2-configuration-and-hashing.md#a-third-instantiation-the-generalized-tableau)):
Clifford gates update `frame` only, and non-Clifford gates branch `amplitudes`.

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

With the inverse tableau, ordinary measurement is a row read in the canonical
column-major orientation, so it does **not** transpose — the per-measurement
thrash that motivated a permanent second copy is gone. What still prefers the
opposite orientation is *bulk* row multiplication and elimination (canonical
form, batched frame work). Those receive a temporary transpose of the X/Z
quadrants rather than a permanently stored second copy. Because a non-square
`2n × n` bit-matrix transpose is not a swap, the guard is expected to work over
square, padded blocks (as Stim does) and may use a scratch buffer; that padding
and scratch are budgeted here rather than assumed away:

**Shipped (2026-08-11), and the cost model this sketch was missing.** The guard
exists as described and is re-entrant, so a batch of measurements amortizes one
transpose (`GeneralizedTableau::measure_all` / `measure_many` hold it, as Stim's
`TableauSimulator` does around a run of collapses). Two facts the sketch does
not account for:

- The transpose had a **floor**, not a cost proportional to `n`: it moved whole
  `64 × 64` blocks, so an `n = 12` frame cost the same as an `n = 64` one, and
  the pair (enter plus restore) was roughly `3072 · ⌈n/64⌉²` word operations.
  Paying that per measurement is what made a `.stim` `MR` sweep — one target at
  a time, no batch to amortize over — 46× slower than the row-major engine it
  replaced. **The floor is gone as of 2026-08-12**: a block populated only in
  its top-left `e × e` corner is transposed over a span of
  `e.next_power_of_two()` rows, because a shift-mask round with `j ≥ e` is the
  identity and every smaller round stays inside the corner. The pair is now
  `8 · ⌈n/64⌉² · (span/2) · log₂ span` exchanges, which for `n = 2` is eight and
  for `n ≥ 64` is what it always was.
- A fold over `k` selected generators can dodge the guard entirely by
  **gathering** those `k` generators out of the column-major arena
  (`TableauData::gather_row`, `O(n)` strided bit reads each). `blocks::prefer_gather`
  picks the cheaper route, which bounds the worst case at the transpose cost
  either way.

Gathering is a *mitigation, not a fix*: at `k · n` bit reads it is still
asymptotically worse than the row-major engine's `k · ⌈n/64⌉` word operations.
The [inverse tableau](#inverse-tableau) removed the folds it was mitigating, so
what still reaches this trade is the *elimination* — the projection's inverse-sign
update, whose site reads are generators — and it is settled the same way
(`prefer_gather`: a sparse frame gathers, a dense one takes the guard). See
[Open questions](#open-questions) item 5 for what that leaves.

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
    fn measure<R: Rng>(&mut self, qubit: usize, rng: &mut R) -> Option<bool>;

    fn measure_many<R: Rng>(&mut self, targets: &[usize], rng: &mut R) -> Vec<Option<bool>> {
        targets.iter().map(|&q| self.measure(q, rng)).collect()
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

### Randomness is injected, not stored

`measure` takes `rng: &mut R`; the tableau owns **no** RNG field. Storing the
RNG inside the tableau is a correctness hazard: classical-mixture branching
clones a tableau into two branches, and a cloned RNG would make both branches
draw the *same* random outcomes — a silent sampling bias. With randomness
injected, `clone` is a pure data copy (no stream to duplicate) and the caller
derives an independent stream per branch (`SmallRng::from_rng(&mut *rng)`). It
also makes `Tableau` trivially `Send + Sync` and removes RNG from the
structural-hash exclusion list — there is nothing to exclude.

### Inverse tableau

**Shipped 2026-08-11** as `crates/ppvm-tableau-2/src/inverse.rs` (the algebra)
and `src/storage/inverse.rs` (the storage), with two corrections to the sketch
below, both of which made it cheaper than designed:

1. **No second bit matrix, and no dual update.** The inverse's *bits are already
   in the arena.* `x`-bit `j` of `ix_q = U†X_qU` is `ω(X_q, s_j)`, the Z bit of
   stabilizer `j` at qubit `q` — which is bit `j` of the contiguous major
   `major(Stab, Z, q)`. All four planes of a qubit's two inverse rows are forward
   majors, so the canonical column-major frame *is* the inverse tableau held
   row-contiguously (Stim's `inv_state` layout). What is not free is signs:
   `U†X_qU` is Hermitian, so `2n` bits, which the gates maintain by prepending
   `G†`'s one-qubit table before their own sweep.
2. **The measurement projection does not need the fold either.** The elimination
   is a sequence of Cliffords *appended* to the frame (`U ↦ U·V`, `V` fixing
   `|0…0⟩`), so it conjugates each inverse row by `V†` — which the existing
   `blocks` kernels do, driven over site-planes instead of qubit columns. The
   signs therefore survive a case-a projection instead of being abandoned, which
   is what keeps the `O(1)` reads available across a measurement sweep.

The `X`/`Z` sign reads are orientation-free (a sign is one bit, not a row), so
they serve a batched sweep under the row guard as well as a standalone
measurement. Only the `Y` correction, a product of the qubit's two rows, needs
the canonical orientation.

The rationale, as designed:

For forward rows `dᵢ = U Xᵢ U†`, `sᵢ = U Zᵢ U†` and a Pauli `P`,

```text
ω(P, sᵢ) = ω(U†PU, Zᵢ) = x-bit i of U†PU
ω(P, dᵢ) = ω(U†PU, Xᵢ) = z-bit i of U†PU
```

and `U†PU` is exactly one row of the inverse tableau. So
`GeneralizedTableau::compute_decomposition` — two column reads for the
anticommutation masks plus a fold of `k` whole generators for the `ℤ/4` residual
— becomes **one contiguous row read for both masks** plus an `O(1)` phase
correction. That fold was the single hottest thing in the crate: it runs on every
measurement, every `T`, every rotation, every expectation value. The correction
is the fixed stabilizers-then-destabilizers ordering convention
(`stab_destab_commute_sign`: reordering shifts the phase by
`2·⟨destab_anticomm, stab_anticomm⟩`), so it is a popcount rather than new
mathematics. The determinism check becomes a contiguous scan, as in Stim's
`is_deterministic_z` (`!inv_state.zs[t].xs.not_zero()`).

The forward tableau must stay: the amplitude-carrying algorithm is stated over
forward rows. So this is a **dual reading of one matrix** (forward columns are
inverse rows), not a dual orientation and not a second copy — which is the answer
to open question 4.

The inverse is a *private representation choice* of `Tableau`: the `Clifford`
kernels encode the sign update rule and the measurement route reads the signs,
but the behavioral traits above them are unchanged. The public
`SymplecticColumns` / `PhaseTrack` primitives move bits and signs
*independently*, so they correspond to no single Clifford and **abandon** the
signs (`invalidate_inverse`); so does `StabilizerFrame::row_multiply`, which can
multiply two anticommuting generators and leave one non-Hermitian. When the flag
is clear every reader falls back to the fold it replaced, so a missing rule costs
speed and never correctness. Reference-frame sampling composes on top for the
sampling workload.

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

### Shared algebra traits

The gate rows above are exactly the Pauli algebra primitives defined in
[`traits-2-configuration-and-hashing.md`](traits-2-configuration-and-hashing.md#pauli-algebra-traits-symplectic-structure-and-phase),
realized at tableau width. The tableau implements:

- `SymplecticColumns` — `swap_xz`, `xor_x_col`, `xor_z_col` as SIMD-block
  operations over the `2n` generator rows (the "column pair" access above);
- `PhaseTrack` — the \(\mathbb{Z}_2\) sign plane plus the Aaronson–Gottesman `g`
  rule for row products; and
- `StabilizerFrame` — the role-exclusive operations (measurement, pivot search,
  `row_multiply`, canonicalization) that read the rows as a stabilizer/
  destabilizer basis.

Implementing the first two yields the blanket `Clifford` impl for free, so the
tableau shares one audited copy of the symplectic sign logic with
`PhasedPauliWord`; only the phase algebra and the frame operations are
tableau-specific. `swap_xz` at tableau width is a block swap and `xor_x_col` is a
whole-column `⊕`, which is why the column-major layout matters: the shared gate
*sequence* is width-agnostic, but the tableau's implementation of each primitive
is the SIMD path the word never needs.

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
    // One sentinel `AtomicU64` in the shipped crate, not three `OnceLock`s: the
    // component split measured slower (three cells to copy and invalidate per
    // clone-and-mutate) than one flattened digest. See
    // `word-data-structures.md` "Hash ownership".
    hash_cache: AtomicU64,
    // no rng — randomness is injected at `measure` (see above)
}
```

The cache fields are private representation choices made by the tableau
author. `Indexable` exposes only the finalized digest via `key_hash()` — which
here composes the component caches (`combine(xz_hash, phase_hash, loss_hash)`)
and applies the tableau's own finalization fold; it does not name cache types or
expose invalidation.

Because equality and hashing compare *canonical ranges* in bulk (zeroed
padding), they require the tableau to be in its **canonical column-major
orientation**. This is not left to discipline: the transposition guard holds
`&mut self` for its whole lifetime, so the borrow checker forbids any shared
`&self` hash or comparison while the tableau is transposed. Hashing therefore
only ever observes the canonical orientation.

Equality and hashing include logical qubit count, generator order, all logical
X/Z bits, phases, and loss state. They exclude:

- allocation capacity;
- alignment padding;
- cache values and validity flags; and
- current physical orientation.
  (There is no RNG state to exclude — it is not stored.)

The component invalidation rules are:

| Mutation | X/Z cache | Phase cache | Loss cache |
| --- | --- | --- | --- |
| Pauli `X`, `Y`, `Z` | preserve | invalidate if changed | preserve |
| Direct phase change | preserve | invalidate | preserve |
| `H`, `S`, `CNOT`, `CZ` | invalidate | invalidate if changed | preserve |
| Row multiplication | invalidate | invalidate if changed | preserve |
| Toggle a loss bit | preserve | preserve | invalidate |
| Physical transpose | preserve | preserve | preserve |

The current `ppvm-tableau-sum` split between `word_fingerprint` and
`phase_loss_hash` is evidence that component hashing matters. The new tableau
owns X/Z, phase, and loss components directly.

## Cloning and mixture use

Classical-mixture branching can clone tableaus frequently. Contiguous backing
storage makes cloning a bulk memory copy. Because the RNG is not stored (it is
injected at `measure`), the clone is *pure data*: there is no random stream to
duplicate, so branches cannot end up statistically correlated — the caller
derives an independent stream per branch. A clone may copy valid hash caches
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
algorithms. PPVM follows Stim in reading the inverse tableau (see
[Inverse tableau](#inverse-tableau) — shipped), and reference-frame sampling
composes on top; both are internal choices of the concrete `Tableau` and do not
surface in the PPVM trait system.

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

1. ~~What block width and alignment should the first implementation use?~~
   **Resolved:** `u64` blocks, 32-byte alignment, block type private rather than
   generic. See [Physical storage](#physical-storage).
2. ~~Should phases occupy one or two bits per generator in the tableau model?~~
   **Resolved:** two bit planes, `(low, high)`, so a phase is a `ℤ/4` value and
   the `g`-rule's carry-save accumulation is two XOR planes. One bit cannot hold
   the imaginary residual that row multiplication produces.
3. ~~Which operations should receive a transposition guard versus performing
   column-strided work directly?~~ **Resolved by measurement, not by taste:**
   the projection runs column-strided (it multiplies one pivot into many
   generators, which transposes cleanly into plane expressions); folds of many
   generators into one Pauli choose per call between gathering and the guard;
   batched measurement holds one guard for the whole batch.
4. ~~Does a dual-orientation representation still add anything now that the
   inverse tableau keeps measurement in the gate orientation?~~ **Resolved:**
   neither. Forward columns *are* inverse rows, so one matrix serves both
   readings and only `2n` signs are duplicated. See
   [Inverse tableau](#inverse-tableau).
5. ~~Row-oriented folds are the layout's remaining cost~~ — **resolved for the
   folds** by the [inverse tableau](#inverse-tableau), which deleted them:
   `surface_d30` went 29 s/shot → 590 ms/shot (the row-major engine it replaces
   is 630 ms), and a `measure_all` sweep is 2.9–3.8× faster than before the
   inverse at `n = 85…1889`. **What is left is the elimination, in one specific
   shape:** an *unbatched* case-a measurement on a *dense* frame, where the
   inverse-sign update's site reads take the guard and the transpose pair is ~69%
   of the measurement (`tableau-micro/msd_sweep_loop`, 9× the row-major engine;
   the batched form of the same sweep is 2.5×).

   **Partly resolved (2026-08-12), for small `n`.** The 69% was two separate
   mistakes, both now fixed and both cheap:

   - The transpose was floored at the block granularity, so a 2-qubit guard cost
     the same as a 64-qubit one — a flat ~1.3 µs per unbatched case-a
     measurement from `n = 2` to `n = 16`, against 130–200 ns for the row-major
     engine. `transpose_square` now truncates the shift-mask ladder to the
     populated corner (see [Transposition guard](#transposition-guard)), which
     removed `transpose_block` from the profile entirely: at `n = 2` the
     measurement went 1.37 µs → 0.52 µs, and `ppvm-stim`'s small circuits
     went `bell_pair` 1.37 → 0.53 µs, `feedback_cx` 1.31 → 0.47 µs,
     `ghz` 1.57 → 0.78 µs, `repetition_code_d3_r3` 2.21 → 1.42 µs.
   - `prefer_gather`'s cost model was `n²/16`, which is neither the floored
     sub-64 cost nor the `⌈n/64⌉²` step above it. It now compares `k·n` bit reads
     against the transpose pair's actual exchange count. That is worth 15% on the
     85-qubit unbatched MSD sweep (2772 → 2361 ns/measurement) and 5–8% on
     `cultivation_d5`; it does not change the decision for the large sparse case
     at all (`surface_d30` is 38.2 → 38.4 ms/shot, i.e. unmoved).

   **What remains at small `n`** is no longer the transpose: profiling a 2-qubit
   case-a measurement puts **46.7% of it inside `libsystem_malloc`**, plus ~6% in
   `RawVec` glue. `update_tableau_according_to_outcome`'s column-major branch
   allocates about ten `Vec`s per projection (two selectors, two
   `PhaseAccumulator`s of three planes each, two `delta_lo.to_vec()` copies)
   where `project_inverse` next to it already carves the same shapes out of the
   reusable `TableauData::take_inv_scratch` arena for exactly this reason. Doing
   the same there is the next step, and would need `PhaseAccumulator` to borrow
   its planes rather than own them. It would not close the gap on its own — 517 ns
   minus ~270 ns of allocator is still ~2× the row-major engine's 130 ns at
   `n = 2` — but it is the largest single item left.

   **What remains at large `n`** is the dense unbatched sweep, unchanged in
   shape: `msd_attrib new` is 2361 ns/measurement against the row-major engine's
   324 ns. Two candidate fixes for it, still neither attempted:
   - **A per-row scan.** In the canonical orientation an inverse *row* is
     contiguous over sites, and the selectors are contiguous masks over sites, so
     the whole sign update can run one row at a time with a prefix-parity scan
     (the running `z_p` parity is a carry) — `O(n · ⌈n/64⌉)` word operations
     instead of `k` gathers or a transpose. It replaces the site-plane kernels
     with a scan per row, so the shipped implementation becomes its oracle.
   - **Lazy orientation.** Leave the frame row-major after a projection and
     restore the canonical orientation only when a gate or a `&self` observation
     (hash, equality) needs it. Cheaper to write, but it gives up the invariant
     that a public method always returns canonically oriented — which is what
     currently makes "hashing only ever observes the canonical orientation" a
     borrow-checker fact rather than a convention.

   Secondary, from the same profile: the projection's inverse update copies each
   site's four planes (`gather_row`) only so the appends can mutate them, but the
   site's own updated bits are needed for exactly one term (`CZ`'s `x'_i`); a
   fused two-append kernel would read the arena directly and drop the copies.
6. Do the batched Clifford sweeps (`h_many`, `cz_block`, 1.7–4.3× slower than the
   row-major engine — one column pass per qubit where it masked all qubits in one
   pass) want a fused whole-plane kernel or a different batching primitive?
