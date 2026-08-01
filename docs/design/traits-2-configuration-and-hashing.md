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

One more bound comes out for the same reason. The old `Coefficient` also required
`Mul<f64, Output = Self>` (scale a coefficient by a real). Once `sin_cos` moves to
`Angle<C>` — whose `sin_cos(&self) -> (C, C)` returns the rotation amplitudes
*already in the coefficient domain* — nothing in the trait set scales a
coefficient by a bare `f64` anymore, so `Mul<f64>` is vestigial. It is also the
**one bound that forecloses exact coefficient rings**: `Complex<f64>` satisfies
`Mul<f64>`, but `GaussianInt` (`ℤ[i]`), `Complex<Rational>`, and cyclotomic
integers cannot (`0.5·(1+i)` leaves the ring). Dropping it lets those exact rings
be `Coefficient`s — which matters because `ppvm-sym` exists precisely for
exact/symbolic coefficients, and formalizing L4 (below) showed the operator
product needs no floats at all. `magnitude() -> f64` stays: an exact ring can
still report a real norm for a `Policy` to threshold.

The value domain keeps only ring arithmetic and a magnitude, replacing `cutoff`
with a property that a policy can threshold:

```rust
/// Value-domain ring arithmetic only — no rotation, no truncation, no bare-`f64`
/// scaling (so exact rings like `GaussianInt` qualify; see above).
pub trait Coefficient:
    PartialEq + Clone + num::Zero
    + Neg<Output = Self>
    + Add<Self, Output = Self> + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + AddAssign<Self> + MulAssign<Self>
    + std::iter::Sum + Send + Sync
{
    fn mul_sign(&self, sign: i8) -> Self;

    /// Nonnegative magnitude. Exposes a property of the value for a `Policy` to
    /// threshold; it does not itself decide any cutoff. Replaces `cutoff`.
    fn magnitude(&self) -> f64;
}
```

Halving (`0.5·x`) is deliberately **not** on `Coefficient` for the same reason
`Mul<f64>` was dropped: `0.5·(1+i)` leaves `ℤ[i]`, so an exact ring could satisfy
a `half()` bound only with a lossy integer `/2` for which `half(x)+half(x) != x`
— re-foreclosing the exact rings this redesign exists to admit. Halving is only
needed by the projective computational-basis measurement kernel (the `(I ± Z)/2`
projectors), which runs over `f64`/`Complex<f64>`, so it becomes a separate
capability, exactly like `sin_cos`→`Angle` and `i`→`ImaginaryUnit`:

```rust
/// Total, exact halving — the capability the `(I ± Z)/2` measurement projector
/// needs. Split from `Coefficient` so exact rings (`GaussianInt`, …) that cannot
/// halve still qualify as `Coefficient`s. Impls must satisfy
/// `x.half() + x.half() == x`.
pub trait Halvable: Coefficient {
    fn half(&self) -> Self;
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

`Word` is the common **read-only** concept for an indexed algebraic monomial. It
is an *inspection* interface — extent, per-site read, weight, and iteration —
consumed by display, serialization, tests, and the sparse-sum plumbing. It is
deliberately **not** the propagation interface: the real gate kernels operate at
*sub-site* (individual X/Z bit) granularity and are algebra-specific, so they
live on the Pauli-specific traits in
[Pauli algebra traits](#pauli-algebra-traits-symplectic-structure-and-phase),
not on `Word`. Keeping `Word` free of mutation also lets an ordered algebra (a
normal-ordered fermionic product) implement it honestly: an in-place
`set(index, site)` at a factor position would violate that algebra's canonical
form, so mutation is relocated to each algebra's own traits rather than asserted
universally here.

```rust
pub trait Word {
    type Site;

    fn n_sites(&self) -> usize;
    fn get(&self, index: usize) -> Self::Site;
    fn weight(&self) -> usize;
    fn iter(&self) -> impl Iterator<Item = Self::Site>;
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
width. Because `Word` is read-only, these differing index meanings (qubit vs
factor order) never have to agree on a shared `set` contract.
`weight()` is the number of non-identity factors according to the concrete site
alphabet; an ordered representation that stores no explicit identities may
therefore have `weight() == n_sites()`. It remains a Pauli-motivated read (the
`MaxPauliWeight` policy needs it) that other algebras may define as is natural.
Structural mutation, and the hash-component invalidation it triggers, belong to
the algebra-specific mutation traits below, not to `Word`.

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
```

A unital Pauli channel acts diagonally in the Pauli basis, `P ↦ λ_P·P`, and its
transfer eigenvalue collapses (using `Σ_Q p_Q = 1`) to
`λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q`, where anticommutation is the
symplectic form `ω(P,Q) = 1`. That algebraic form is machine-checked in
`lean/PPVM/Algebra/Noise.lean` (`pauli_channel_eigenvalue`, and
`pauli_channel_eigenvalue_omega` tying `anti` to `PPVM.Symplectic.omega`); the
zero-state read-out `⟨0ⁿ|ρ|0ⁿ⟩ = Σ_{P X-free} c_P` is `overlap_with_zero_xfree`.

```rust
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

`Clifford` is not implemented by hand on each standard type. It is a *derived*
behavioral trait, blanket-implemented once over the Pauli algebra primitives
described next for the types that opt in via the `BlanketClifford` marker (the
phaseless words and `Tableau`); the sole hand-written exception is
`PhasedPauliWord`'s read-once fused `impl Clifford` (see the marker note below).
`RotationOne` is implemented on `PauliSum` in terms of `PauliBits`, also below.

### Pauli algebra traits: symplectic structure and phase

`Word` is read-only, so the mutation that gate propagation performs lives on
Pauli-specific traits. Their shape is dictated by the algebra rather than
invented: a Pauli operator modulo phase is a vector in the symplectic space
\(\mathrm{GF}(2)^{2n}\) (the X and Z bit planes), Pauli multiplication is
vector addition (`⊕`) *up to phase* (the phase is the \(\mathbb{Z}_4\) cocycle of
the extension below), commutation is the symplectic form
\(\omega(P,Q) = x_P\cdot z_Q \oplus z_P\cdot x_Q\), and conjugation by a Clifford
factors into a symplectic map on the bits together with a phase. Concretely, the
Pauli group is a **non-split central extension**
\(1 \to \mathbb{Z}_4 \to \mathcal{P}_n \to \mathrm{GF}(2)^{2n} \to 1\): the
phases are central, and the extension does *not* split — \(\mathcal{P}_n\) is
nonabelian, so it is a genuine central extension and **not** a semidirect
product \(\ltimes\) (an earlier draft wrote \(\mathrm{Sp}(2n,2)\ltimes\text{phases}\),
which is inaccurate on both counts). The symplectic group enters one level up:
conjugation by a Clifford realizes \(\mathrm{Sp}(2n,2)\) acting on the quotient,
\(\mathcal{C}_n/\mathcal{P}_n \cong \mathrm{Sp}(2n,2)\). The non-split extension
is machine-checked **at \(n = 1\)**: `lean/PPVM/Pauli/Phase.lean` builds
\(\mathcal{P}_1\) as a `Group` with the quotient homomorphism
`PhasedPauli.toSymplectic` and the non-split witness `not_commutative`. For
general \(n\), the two facts that make the *packed* multi-qubit product
well defined are lifted — the phase 2-cocycle `phaseExpN_cocycle` (associativity)
and the commutation law `phaseExpN_sub_comm` (`lean/PPVM/Pauli/Word.lean`) — and
the packed \(n\)-qubit phase exponent `phaseExpN` (the summed
`(2*sign_count + imag_count) mod 4` that `key_mul` accumulates) is grounded
against genuine \(n\)-fold tensor-product matrices: `tensorPauli_mul`
(`lean/PPVM/Pauli/Matrix.lean`) proves it is the base-\(i\) exponent of the real
\(2^n\times 2^n\) operator product \(g(p)\cdot g(q)\), lifting the single-qubit
`pauliMat_mul` through the tensor-product phase-multiplicativity
\(\prod_i i^{k_i} = i^{\sum_i k_i}\) (`prod_iuPow`). The n-qubit
`Group`/quotient/non-split objects are not reconstructed, as the
single-qubit case already exhibits the extension's structure. The
\(\mathrm{Sp}\) action is checked in one direction — each Clifford generator is
an \(\mathrm{Sp}\)-isometry, so conjugation *lands in* \(\mathrm{Sp}(2n,2)\)
(`lean/PPVM/Pauli/Symplectic.lean`). Each generator's phase-stripped bit map is
moreover an **involution**, hence a **bijection** of the `Sp n` word space
(`hAct_involutive`/`sAct_involutive`/`cnotAct_involutive`/`czAct_involutive` and
their `*_bijective` corollaries — note `S`'s bit map has order 2 even though the
phased `conjS` has order 4): this is the no-collision guarantee the `PauliSum`
Clifford re-key relies on (`crates/ppvm-pauli-sum-2/src/{producer,clifford}.rs`,
"A Clifford re-key is a bijection, so colliding re-keyed terms never occur").
This per-generator conjugation is also
exhibited as a signed symplectic automorphism, phase included — at \(n = 1\) for
the single-qubit generators (`conjHHom`/`conjSHom` in
`lean/PPVM/Pauli/Conjugation.lean`, e.g. `conjH_Y` gives \(HYH = -Y\)). `S` is the
sole generator whose conjugation sign is convention-sensitive, so the **backward**
direction the phased-word simulator actually runs (\(S^\dagger P S\)) is pinned
separately as `conjSdag_sign` (sign \(x\wedge\lnot z\), i.e. \(S^\dagger X S = -Y\)),
matching the fused `Clifford::s` sign
(`crates/ppvm-phased-pauli-word-2/src/clifford.rs`); `conjS_conjSdag` fixes it as
the genuine inverse of the forward `conjS`. And at
\(n = 2\) for the independent two-qubit generators `CNOT`/`CZ` (`conjCNOTHom`/
`conjCZHom` on the group \(\mathcal{P}_2\) in the same file). The latter pin the
\(\mathbb{Z}_4\) conjugation-*phase* delta that `Symplectic.lean`'s bit-only
`cnotAct_isometry`/`czAct_isometry` do not cover: `conjCNOT_sign`/`conjCZ_sign`
match the fused `Clifford::cnot`/`Clifford::cz` signs
(`crates/ppvm-phased-pauli-word-2/src/clifford.rs`, ported verbatim from the old
`crates/ppvm-pauli-word/src/phase/clifford.rs`), and `conjCNOT_Xc`…/`conjCZ_Xc`…
check the maps against the standard tableau tables so the phase is forced, not
chosen. Every one of these conjugation deltas is an even \(\mathbb{Z}_4\) value
(\(\in\{0,2\}\), a real \(\pm1\), never \(\pm i\)), so a phase that starts real
stays real: `conjH_isRealPhase`/`conjS_isRealPhase`/`conjSdag_isRealPhase` and
`conjCNOT_isRealPhase`/`conjCZ_isRealPhase` (over `IsRealPhase` +
`isRealPhase_zero`) prove exactly the reachability invariant that keeps
`crates/ppvm-pauli-sum-2/src/clifford.rs`'s `clifford_sign` \(\pm1\) drain total —
its `PosI`/`NegI` "bug" branch is unreachable, not merely `debug_assert!`ed away.
The surjectivity that upgrades this containment to the full
isomorphism \(\mathcal{C}_n/\mathcal{P}_n \cong \mathrm{Sp}(2n,2)\) is stated here
but not formalized. The **loss-guarded** variant of this action — each generator
is a no-op when any operand qubit is lost, as in
`crates/ppvm-lossy-pauli-word-2/src/clifford.rs` — is machine-checked separately
in the same file (`lean/PPVM/Pauli/Symplectic.lean`, the `…ActL` definitions):
the guard preserves the canonical loss invariant `lost[q] ⇒ x[q]=0 ∧ z[q]=0`
(`hActL_preserves_loss`/`sActL_preserves_loss`/`cnotActL_preserves_loss`/
`czActL_preserves_loss`, with the critical present-control/lost-target `CNOT`
case as `cnotActL_lost_target_stays_identity`), and on the present-qubit
sub-block it coincides with the `Sp(2n,2)` isometry above
(`hActL_present_isometry`/`sActL_present_isometry`/`cnotActL_present_isometry`/
`czActL_present_isometry`). Because the blanket `Clifford` does not skip a
`CNOT` atomically but composes two independently guarded column primitives
(`xor_x_col(c,t)` then `xor_z_col(t,c)`, both testing `lost c ∨ lost t`), the
per-primitive form of the guard is checked too: each guarded column preserves the
loss invariant on its own (`xorXColL_preserves_loss`/`xorZColL_preserves_loss`)
and their composition equals the atomic whole-gate skip
(`xorZColL_xorXColL_eq_cnotActL`) — the machine-checked form of the crate's
"reproduces the old whole-gate skip" claim, and the link that would break if
either guard were weakened in isolation. (`CZ` emits a single primitive
`cz_bits`, which is exactly `czAct`, so `czActL` already models it.) That
factorization sorts every gate operation into two buckets, and the traits follow
the buckets.

**Role-independent (the \(\mathrm{Sp}\) part).** Conjugating a Pauli by a
Clifford gate does the same bit-plane algebra whether the operator is a lone
word (one row) or one of a tableau's \(2n\) stabilizer/destabilizer generators.
`H` swaps the X and Z columns; `S` does `z ⊕= x`; `CNOT` does
`x_t ⊕= x_c, z_c ⊕= z_t`. This logic is written **once**.

**Role-dependent (the phase extension).** A standalone phased Pauli lives in
\(\mathbb{Z}_4\) (`Y = iXZ` needs the `i`); a stabilizer tableau stores a
\(\mathbb{Z}_2\) sign and recovers the `i`'s through the Aaronson–Gottesman `g`
rule during row multiplication. Same gate, different phase algebra — so phase
bookkeeping is written **per type**.

Two primitive traits capture the two buckets, kept separate so they mirror
\(\mathrm{Sp}\) and its extension (and so a future phase-free classical register
can implement only the first):

```rust
/// Sp-part: bit-plane column algebra. `PhasedPauliWord` uses 1-bit columns;
/// `Tableau` uses SIMD blocks over its 2n rows. Same meaning, different width.
pub trait SymplecticColumns {
    fn n_qubits(&self) -> usize;
    fn swap_xz(&mut self, q: usize);
    fn xor_x_col(&mut self, ctrl: usize, tgt: usize); // x_tgt ⊕= x_ctrl
    fn xor_z_col(&mut self, tgt: usize, ctrl: usize); // z_ctrl ⊕= z_tgt
    // Completed in ppvm-traits-2 (the `// ...` below is not exhaustive):
    // `xor_z_from_x(q)` for S (z_q ⊕= x_q) and `cz_bits(a, b)` for CZ.
}

/// Extension-part: the phase algebra. Z4 for a phased word, Z2 + g for a tableau.
pub trait PhaseTrack {
    fn flip_phase_where_xz(&mut self, q: usize);
    fn cnot_phase(&mut self, ctrl: usize, tgt: usize);
    // ...one phase delta per gate; the tableau's g-rule lives behind these.
    // Completed in ppvm-traits-2: `s_phase`, `cz_phase`, and the pure-sign
    // `x_phase`/`y_phase`/`z_phase` (X/Y/Z are phase-only, no bit change).
}
```

`Clifford` is then the role-independent structure, blanket-implemented once for
every type that **opts in** via the empty marker `BlanketClifford`. The *sequence*
of primitives per gate is identical across roles even though the phase primitive
it calls is not — the single audited copy of the symplectic sign logic that would
otherwise be duplicated and drift:

```rust
/// Opt-in marker selecting the shared blanket Clifford (below). A phaseless word
/// implements it; `Phased<W>` deliberately does not (see the fused-impl note).
pub trait BlanketClifford {}

impl<T: SymplecticColumns + PhaseTrack + BlanketClifford> Clifford for T {
    fn h(&mut self, q: usize) { self.flip_phase_where_xz(q); self.swap_xz(q); }
    fn cnot(&mut self, c: usize, t: usize) {
        self.cnot_phase(c, t);
        self.xor_x_col(c, t);  // x_t ⊕= x_c
        self.xor_z_col(t, c);  // z_c ⊕= z_t  (arg order: xor_z_col(tgt, ctrl))
    }
    // ...
}
```

**The `BlanketClifford` marker and the fused phased override.** The blanket runs
the phase primitive and the column primitives as *separate* steps, so it reads
each symplectic bit twice on a type whose bits and phase live apart — once in
`cnot_phase`/`cz_phase` to compute the sign, again in `xor_x_col`/`xor_z_col` to
apply the bit map. For the phaseless words (`PauliWord`, `LossyPauliWord`) and the
future `Tableau` that is free (the phase primitive is a no-op or SIMD-wide), so
they opt in and share the audited copy. But for `Phased<W>` the double read
benchmarked ~1.6–1.8× slower than the old *fused* `PhasedPauliWord::cnot`, which
reads each bit once. On stable Rust a hand-written `impl Clifford for Phased<W>`
may not coexist with an *unconditional* blanket the type would also satisfy
(E0119), so the blanket is gated on `BlanketClifford`: `Phased<W>` stays out of
the marker and instead provides its own **fused** `impl Clifford` in
`crates/ppvm-phased-pauli-word-2/src/clifford.rs`, which reads each inner X/Z bit
once via `PauliBits`, computes the `ℤ₄` sign from those reads, applies the bit
update reusing them, and folds the sign into the stored phase. The signs are
byte-for-byte the old kernel (`crates/ppvm-pauli-word/src/phase/clifford.rs`), so
correctness is unchanged; only the redundant reads are gone (phased `cnot` back at
parity, new/old ≈ 0.86×). The marker keeps the blanket the single audited copy for
the standard types while letting `Phased` win the read-once fusion.

The role-*exclusive* operations — those that interpret the rows as a symplectic
basis rather than as independent operators — do not belong in the shared tower.
They are a tableau-only trait a word never implements. `StabilizerFrame` holds
the frame **primitives**, not `measure` itself: `measure` is the public
`Measure` trait, and its two algorithms (`Tableau` pure-Clifford,
`GeneralizedTableau` coefficient-aware `O(n^2)`) are built *on* these primitives:

```rust
pub trait StabilizerFrame {
    /// Find a generator that anticommutes with the measured Pauli (the pivot).
    fn anticommuting_pivot(&self, qubit: usize) -> Option<usize>;
    /// Multiply generator `src` into `dst` (uses the Aaronson–Gottesman g-rule).
    fn row_multiply(&mut self, src: usize, dst: usize);
    /// Restore canonical form after elimination.
    fn canonicalize(&mut self);
}
```

The frame these primitives operate on is a genuine symplectic basis, and this is
machine-checked in `lean/PPVM/Tableau/Frame.lean`: the `2n` generators satisfy
the symplectic-basis relations `ω(dᵢ,sⱼ) = δᵢⱼ` (`IsSymplecticFrame`), are
linearly independent (`frame_linearIndependent`), start as one
(`isSymplecticFrame_identity`), and stay one under every Clifford generator
(`isSymplecticFrame_hAct`/`sAct`/`cnotAct`/`czAct` via `IsSymplecticFrame.map`).
The `anticommuting_pivot` search rests on the measurement dichotomy
(`measurement_dichotomy`): the outcome is deterministic exactly when the measured
Pauli commutes with every stabilizer (`measure_deterministic_iff_xfree`).

The blanket's implementers (the `BlanketClifford` opt-ins) are the phaseless
words `PauliWord` and `LossyPauliWord` (one row) and `Tableau` (all three).
`PhasedPauliWord` (`Phased<PauliWord>`) is **not** a blanket implementer: it
supplies its own read-once fused `impl Clifford` (see the marker note above), so
it needs neither `SymplecticColumns` nor `PhaseTrack`. Note that `PauliSum` is
deliberately *not* an
implementer: a Clifford gate on a sum re-keys every term (each Pauli maps to a
different Pauli, so the map is rebuilt, not updated in place), so the sum applies
the one-row action pointwise and drains each term's phase delta to its
coefficient. Non-Clifford rotations, which branch one term into several, stay on
a separate mutation primitive:

```rust
/// Mutable single-vector X/Z access — a point of GF(2)^{2n}. Hosts the
/// rotation/branching kernels, which flip individual bits and ship the sign to
/// the coefficient. Implemented by `PauliWord` and `LossyPauliWord`.
///
/// The supertrait is `Word`, **not** `Word<Site = Pauli>`: `LossyPauliWord` also
/// implements `PauliBits` but its `Word::Site` is `LossySite<Pauli>` (a lost site
/// is not a bare `Pauli`), so a `Site = Pauli` bound would exclude it. These
/// methods are alphabet-agnostic (raw X/Z bits + a loss flag); code that needs
/// Pauli *propagation* re-adds `Word<Site = Pauli>` on its own methods, the same
/// pattern the graded algebra uses (bound on `Indexable`, propagation re-adds
/// `Word`/`PauliBits`). `PauliWord` still implements `Word<Site = Pauli>`.
pub trait PauliBits: Word {
    fn x_bit(&self, i: usize) -> bool;
    fn z_bit(&self, i: usize) -> bool;
    fn set_x_bit(&mut self, i: usize, v: bool); // invalidates the hash lazily
    fn set_z_bit(&mut self, i: usize, v: bool);
    fn is_lost(&self, i: usize) -> bool { false } // LossyPauliWord overrides
}
```

`PauliBits` is the narrow bit-level slice of the retired `PauliWordTrait`,
separated from key identity (`Indexable`), inspection (`Word`), and phase — and
it passes the trait admission rule because generic rotation kernels consume it
across two implementers. Loss *reads* are a defaulted predicate; loss *writes*
and the binary Pauli product (`Mul`) stay inherent on the concrete words, since
each has a single implementer.

An algorithm should take its independent choices as direct type parameters.
An associated-type bundle is not useful merely because it replaces two type
parameters with one. In particular, there is no `PauliSumAlgorithm` trait that
bundles a term map with a policy: storage layout and policy are orthogonal
choices.

The choices that are *intrinsic to a storage instance* — which key it maps and
which coefficient it accumulates — are associated types of that storage rather
than free parameters, the same way `HashMap<K, V>` fixes `K` and `V` per concrete
map. `Accumulate` (the graded map trait, below) therefore exposes `type Key` and
`type Coeff`, and the reusable sparse-sum shape needs only its storage container
and its policy:

The storage parameter is simply an `Accumulate` container (the graded algebra of
[The map is a graded algebra over `C[K]`](#the-map-is-a-graded-algebra-over-ck))
— there is no separate `SumStorage` trait and, crucially, **no owned workspace**:

```rust
pub struct Sum<S, P = NoPolicy>
where
    S: Accumulate,
    P: Policy<S::Key, S::Coeff>,
{
    storage: S,     // Vec<(K, C)>, HashMap<K, C, IdentityBuildHasher>, or ColumnStore
    policy: P,
    /// Invariant: every key `k` in `storage` satisfies `k.n_sites() ==
    /// n_sites`. An empty sum has no key to derive the width from, so the
    /// field is carried explicitly and checked by a `debug_assert!` on every
    /// insertion path.
    n_sites: usize,
}
```

`Sum` owns **only** its storage, policy, and width — no auxiliary map, no scratch
buffer. Clone is therefore pure data, which matters because a stabilizer mixture
clones frequently. The key and coefficient are `S::Key` / `S::Coeff` (associated
types of `Accumulate`), so there is no phantom axis and no way to pair a storage
type with the wrong key. Pauli propagation re-adds the `Word` / `PauliBits` bound
on *its* methods (`S::Key: Word<Site = Pauli>`), not on the engine.

#### The convenience bundle

Collapsing to `Sum<S, P>` restores ergonomics through a type-alias family, not a
foundational trait. Each alias fixes a domain's axes and defaults the common
knobs, so a one-token `PauliSum` name returns:

```rust
// HashMapStore<K, C> = HashMap<K, C, IdentityBuildHasher>; the storage is the
// bare container, the alias only names it and bakes in the pass-through hasher.
pub type PauliSum<C = f64, P = NoPolicy>      = Sum<HashMapStore<PauliWord, C>, P>;
pub type LossyPauliSum<C = f64, P = NoPolicy> = Sum<HashMapStore<LossyPauliWord, C>, P>;
pub type FermionSum<C = f64, P = NoPolicy>    = Sum<HashMapStore<FermionWord, C>, P>;
```

A user who wants to pin a reusable configuration writes one `type` line instead
of an `impl Config`:

```rust
type MyPauliSum = Sum<DashMap<PauliWord, Complex<f64>, IdentityBuildHasher>, MaxPauliWeight>;
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
        M: Retain<W, C>;
}
```

`truncate` bounds on `Retain`, not the map algebra: dropping terms still in the
support breaks module exactness, so it is the one non-algebraic map operation
and lives on its own capability that `Policy` — not the algebra — consumes.
`Retain` is `fn retain(&mut self, keep: impl Fn(&W, &C) -> bool)`, implemented by
both storage backends.

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

    fn truncate<M: Retain<W, C>>(&self, map: &mut M) {
        map.retain(|_word, coeff| coeff.magnitude() >= self.threshold);
    }
}
```

The keep-rule is `magnitude() >= threshold`, so a term whose magnitude is
*exactly* the threshold is **kept**. The stabilizer-tableau path instead keeps on
a strict `magnitude() > threshold`, so the two backends disagree at
`|c| == threshold` (kept by `CoefficientThreshold`, dropped by the tableau). This
boundary mismatch — at every threshold — is machine-checked in
`lean/PPVM/Algebra/Truncation.lean` (`cutoff_mismatch`). The error a truncation
incurs is bounded in the same file: an `ℓ¹` bound `|error| ≤ Σ_{dropped} |c_P|`
for the `PauliSum` path (`l1_bound`, under `|⟨P⟩| ≤ 1`) and, for the tableau path,
the *unconditional* Cauchy–Schwarz `ℓ²` bound `l2_bound`, sharpened to
`error² ≤ (Σ_{dropped} c_P²)·|D|` under `|⟨P⟩| ≤ 1` in `l2_bound_normalized`.

`Policy` and its concrete implementations belong to the sparse-sum crate. This
removes the current split where the `Strategy` trait lives in `ppvm-traits` but
its concrete strategies live in `ppvm-pauli-sum`; the policy is not an
algorithm-agnostic `ppvm-traits-2` concern.

#### The map is a graded algebra over `C[K]`

The associative coefficient map — implemented today by `HashMap`, `IndexMap`,
and `DashMap` — is, algebraically, the **free `C`-module on a set of keys `K`**:
a finitely-supported function `K ⇀ C`, an element of `C[K]`. The keys are
whatever is `Indexable` — `PauliWord` for `PauliSum`, `Tableau` for a stabilizer
mixture — so the *same* algebra serves both. Every operation it must support is a
module (or, at the top, ring) operation, so its capabilities are **graded by
algebraic strength** rather than named ad hoc. This replaces the flat
`ACMapBase` / `ACMapIter` / `ACMapAddAssign` / `ACMapInsert` / `ACMapRetain` /
`ACMapConsume` split — which the design's own admission rule would reject as six
grandfathered names — with layers, each justified by a distinct algebraic
property *and* a distinct consumer:

| Layer | Algebra | Trait | Justifying consumer |
| --- | --- | --- | --- |
| L0 | finite partial function `K ⇀ C` | `Support` | everything; read-out / export |
| L1 | abelian-monoid formation + canonical support | `Accumulate` (`accumulate_batch` + `reduce`) | forming any linear combination |
| L2 | the `C`-module action | `Scale` | normalization, global factors |
| L3 | the trace pairing (symmetric bilinear `overlap`; sesquilinear `hermitian_overlap`) | `Pair` (`probe_batch`, `overlap`, `hermitian_overlap`) | expectation read-out; complex state overlap |
| L4 | the ring product (a *twisted* group algebra for Paulis: the key product emits a phase) | `Multiply` (needs a key product) | operator composition, squaring an observable |

The **minimum** is L0 + L1: a finitely-supported `K ⇀ C` you can accumulate
into and reduce. That *is* the algebraic essence of a sparse sum. Truncation is
deliberately **absent** from this table: dropping terms that are still in the
support is an approximation that breaks module exactness, so it is not an
algebraic operation — it belongs to `Policy` (a `Retain` capability), sitting
outside the algebra. That is why truncation was always awkward to place.

The key type is **`Indexable`, not `Word`**: `C[K]` is the free module over any
index set, so the algebra needs only a valid map key. Requiring `Word` here would
leak Pauli-specificity into the general algebra and block the tableau mixture,
whose key is a `Tableau` (an `Indexable` that is not a `Word`). Pauli-specific
propagation re-adds the `Word` / `PauliBits` bound on *its* methods, not on the
algebra.

```rust
/// L0 — the container. Note: no `&mut (K, C)` and no `&mut [C]` slot access is
/// exposed, so a columnar (structure-of-arrays) backend is expressible. `iter`
/// is read-only export; a SoA backend synthesizes the pairs from its columns.
pub trait Support {
    type Key: Eq + Clone;           // minimal; hash backends add `Indexable`
    type Coeff: Coefficient;
    fn len(&self) -> usize;
    fn get(&self, key: &Self::Key) -> Option<Self::Coeff>;
    fn iter(&self) -> impl Iterator<Item = (Self::Key, Self::Coeff)>;
}

/// L1 — the module core: form linear combinations, then canonicalize.
pub trait Accumulate: Support {
    /// Build side of the hash join: merge a produced batch, accumulating onto an
    /// existing key or inserting a new one. Columnar in; the scalar
    /// `accumulate(k, c)` is provided sugar over a batch of one.
    fn accumulate_batch(&mut self, terms: &TermBatch<Self::Key, Self::Coeff>);

    /// Canonicalize to reduced finite-support form: drop every key whose
    /// coefficient `is_zero()`. First-class and run **only here** — see below.
    fn reduce(&mut self);
}

/// L2 — the C-module action: a pure elementwise map over the coefficients.
pub trait Scale: Support {
    fn scale(&mut self, s: &Self::Coeff);   // ∀ k. c_k *= s
}

/// L3 — the read side of the hash join. Two pairings live here, differing only in
/// whether the first operand is conjugated.
///
/// `overlap` is the **symmetric bilinear** Hilbert–Schmidt trace pairing
/// `⟨A, B⟩ = ∑_k a_k b_k = Tr(A B)/2ⁿ` on the Hermitian Pauli basis — bilinear,
/// *not* conjugated. Correct for expectation values of Hermitian observables (the
/// coefficients carry the physical structure). Its full `C`-bilinearity —
/// biadditivity (`overlap_add_left`/`overlap_add_right`), homogeneity in each slot
/// (`overlap_smul_left`/`overlap_smul_right`), and symmetry over a commutative
/// coefficient ring (`overlap_comm`) — is machine-checked in
/// `lean/PPVM/Algebra/GradedMap.lean`; and orthonormality of the basis monomials
/// under this abstract pairing (`PPVM.Noise.overlap_single_single`) is checked in
/// `lean/PPVM/Algebra/Noise.lean` — the model pairing `∑_k a_k b_k` stands in for
/// the normalized trace `Tr(A B)/2ⁿ` (which is *not* itself constructed in Lean),
/// so `overlap_single_single` is the model form of `Tr(P Q)/2ⁿ = δ`, not that
/// matrix identity verbatim. Finally, the semantic link to the Clifford path is
/// closed: a Heisenberg re-key `P ↦ φ_G(P)` by the `Sp(2n,2)` bijection with the
/// drained `±1` sign folded into the coefficient **preserves this pairing**
/// (`overlap(conj_G A, conj_G B) = overlap(A, B)`), machine-checked as
/// `clifford_conjugation_preserves_overlap` in `lean/PPVM/Algebra/GradedMap.lean`
/// (over `overlap_eq_fintype_sum`, composing `Symplectic.*_bijective` with the
/// sign reality `s_P² = 1` of `Conjugation.*_isRealPhase`) — the guarantee tying
/// `Sum::overlap` to `Sum::{h,s,cnot,cz}` (`ppvm-pauli-sum-2/src/clifford.rs`).
///
/// `hermitian_overlap` is the **sesquilinear** inner product
/// `⟨φ | ψ⟩ = ∑_k conj(a_k)·b_k`, conjugate-linear in the first argument. This is
/// the physical state/amplitude overlap — `GeneralizedTableau`'s `C[Bitstring]`
/// amplitude vector needs it to compute `⟨φ|ψ⟩` correctly over a complex ring. It
/// requires the coefficient ring to carry a conjugation, so it is bounded on
/// `Conjugate`. Over a real ring `conj` is the identity and the two pairings
/// coincide; the machine-checked properties (conjugate symmetry, sesquilinearity,
/// `⟨f,f⟩ ≥ 0`) are in `lean/PPVM/Algebra/GradedMap.lean`.
pub trait Pair: Support {
    fn probe_batch(&self, keys: &KeyBatch<Self::Key>, out: &mut [Option<Self::Coeff>]);

    fn overlap(&self, other: &Self) -> Self::Coeff;

    fn hermitian_overlap(&self, other: &Self) -> Self::Coeff
    where
        Self::Coeff: Conjugate;
}

/// A coefficient ring carrying a ring involution (a commutative `*`-ring):
/// complex conjugation on `Complex<f64>` / `GaussianInt` (`ℤ[i]`) / cyclotomic
/// integers, and the identity on real rings. It supplies exactly the conjugation
/// the sesquilinear `Pair::hermitian_overlap` needs; nothing in propagation
/// requires it, so — like `ImaginaryUnit` for L4 — it is a separate capability,
/// not a `Coefficient` bound (keeping the base trait minimal and exact-ring
/// friendly).
///
/// Laws (commutative `*`-ring): `conj(conj(a)) == a`, `conj(a + b) ==
/// conj(a) + conj(b)`, `conj(a · b) == conj(a) · conj(b)`; and, when the ring is
/// also `ImaginaryUnit`, `conj(i) == −i`. The exact ring `GaussianInt` (`ℤ[i]`)
/// realizes this capability with its genuine `StarRing` conjugation, and the
/// `conj(i) == −i` law on it is machine-checked in `lean/PPVM/Pauli/Matrix.lean`
/// (`star_iU`).
pub trait Conjugate: Coefficient {
    fn conj(&self) -> Self;
}

/// A key whose set carries a product — the (projective) group structure that
/// lifts `C[K]` from a module to an algebra. `PauliWord` implements it: the
/// phased Pauli product is `v·w = i^k (v⊕w)` with `k = phaseExp(v,w) ∈ ℤ/4ℤ`
/// (no separate `±`: `i^k` already spans `{1, i, −1, −i}`). Note the keys form a
/// group only *up to phase* — the product is **not closed on keys**, it emits an
/// `i^k` — which is exactly why `key_mul` returns `(Self, Phase)` rather than
/// `Self`, and why `C[PauliWord]` is a **2-cocycle-twisted** group algebra, not
/// the plain group algebra of `(𝔽₂^{2n}, ⊕)`. That `phaseExp` is a genuine
/// 2-cocycle (hence the twisted product is associative) is machine-checked in
/// `lean/PPVM/Pauli/Phase.lean` and `lean/PPVM/Algebra/Twisted.lean`. A bare
/// `Tableau` mixture key carries no such product.
pub trait KeyProduct: Eq + Clone {
    /// Product of two keys, with the phase it produces (folded onto the coeff).
    fn key_mul(&self, other: &Self) -> (Self, Phase);
}

/// The phase capability L4 needs, over a **commutative** coefficient ring: a
/// distinguished **primitive fourth root of unity** `i` (`i·i == −one()`, hence
/// `i⁴ = 1`). `key_mul` folds `i^k` (`k ∈ ℤ/4`) onto the coefficient, and over any
/// commutative coefficient ring the twisted product is associative —
/// machine-checked in `lean/PPVM/Algebra/Twisted.lean` (`tmul_assoc`, under
/// `[CommRing C]` and `i⁴ = 1`), with the exact ring `GaussianInt` given as a
/// worked instance in `lean/PPVM/Pauli/Matrix.lean`. Commutativity is
/// load-bearing, not incidental: the associativity proof commutes the scalar
/// `i^k` factors past coefficients. It holds for every ring L4 targets — `ℂ`,
/// `ℤ[i]`, `Complex<Rational>`, cyclotomic integers are all commutative — so L4
/// assumes commutative coefficient multiplication (the general `Coefficient`
/// trait does not itself require it). The `i⁴ = 1` requirement is strictly weaker
/// than the earlier `ComplexCoefficient` (`Complex<f64>`) bound: `GaussianInt`
/// (`ℤ[i]`), `Complex<Rational>`, and cyclotomic integers all satisfy it, so L4
/// does **not** foreclose exact or symbolic Pauli multiplication.
pub trait ImaginaryUnit: Coefficient + num::One {
    /// The imaginary unit `i`; impls must satisfy
    /// `Self::imaginary_unit() * Self::imaginary_unit() == -Self::one()`.
    fn imaginary_unit() -> Self;
}

/// L4 — the ring product. The only layer that needs the *key* to carry a
/// product; it stays optional and is not implemented for a key type that has
/// none. The Pauli product injects powers of `i`, so the coefficient must absorb
/// phase — bounded on `ImaginaryUnit`, the minimal requirement (a primitive
/// fourth root of unity), **not** the stronger `ComplexCoefficient`.
pub trait Multiply: Accumulate
where
    Self::Key: KeyProduct,
    Self::Coeff: ImaginaryUnit,
{
    fn multiply_into(&self, other: &Self, acc: &mut Self);
}
```

#### `reduce()` is first-class, and runs only at finalize

Reduce (drop zero coefficients) is **not** an inline check during accumulation.
A coefficient `c1 + c2 + c3` is one value; the fact that a partial sum `c1 + c2`
transiently hits zero during a merge is not the element reaching zero. Pruning
mid-accumulation would delete a key that a later contribution re-creates — churn
at best, and unsafe if it mutates the map during the traversal that is still
producing those contributions. So canonicalization runs **once, after all
coefficients for every word are accumulated**. Making `reduce` a named
operation guarantees that by construction, and it earns its place three ways:

- **Correctness** — there is no inline drop to get wrong (the caution above).
- **It is a bulk primitive with a backend-specific implementation** — a scalar
  `retain(|_, c| !c.is_zero())` on a `HashMap`, but a prefix-sum **stream
  compaction** kernel on a `ColumnStore`. Folded into `consume`/swap it could
  not have a bulk form.
- **It amortizes** — several `accumulate_batch` calls (or a `multiply_into`
  outer product) can precede a single `reduce`, instead of canonicalizing after
  every step.

#### Backends are containers; columnar is expressible from day one

The graded traits are `impl`'d **directly on the container** — no wrapper types.
The lookup strategy is a *collection* choice by support size, and the memory
layout is a separate *layout* choice by execution target:

- **`Vec<(K, C)>`** — an unsorted coordinate list, linear-scan `accumulate`. Best
  for small support (the `GeneralizedTableau` amplitude vector); this is today's
  `SparseVector`. Requires only `K: Eq + Clone` — it never hashes.
- **`HashMap<K, C, IdentityBuildHasher>`** — hash-join `accumulate`. Best for
  large support (`PauliSum`); this is today's `ACMap` path. (AoS, scalar CPU.)
  Requires `K: Indexable` (the direct digest).
- **`ColumnStore`** — the one backend that *must* be a new struct, because it is
  structure-of-arrays: coefficients in one contiguous column, keys in plane
  blocks. Same hash-join build as the `HashMap`, but `scale` is one vectorized
  `*=`, `reduce` is a prefix-sum compaction, and `probe_batch` uses coalesced
  gathers. It is `HashMap`'s data re-laid for SIMD / GPU, not a third collection.
  Requires `K: Indexable + Columnar`.

So `Indexable` is a *per-backend* requirement, not a universal key bound: the
`Vec` path needs only `Eq + Clone`, and only the hash-indexed backends demand the
digest.

The whole point of the refactor is to unblock SIMD / GPU, so the traits are
designed so `ColumnStore` is expressible from day one — **only because no trait
signature leaks the array-of-structs layout**. The rules the whole storage
contract holds to:

1. coefficients are a contiguous column and keys are plane blocks (SoA) — so L2
   Scale, L1 Reduce, and L3 Pair vectorize and GPU loads coalesce;
2. term I/O is a columnar `TermBatch` / `KeyBatch` (see
   [The batch contract](#the-batch-contract)), never a scalar per-element
   callback — so the `HashMap` accepts a batch and processes it scalar internally
   while `ColumnStore` vectorizes it, from the *same* call;
3. the mutating operations are whole-map (`scale`, `reduce`) or batch
   (`accumulate_batch`), never per-slot; and
4. **no signature exposes `&mut (W, C)` or `&mut [C]`** — the one rule that keeps
   the AoS layout from becoming load-bearing and foreclosing `ColumnStore`.

The columnar term types (`TermBatch`, `KeyBatch`, `KeyColumn`) that these layers
consume are the same ones specified in
[Batch execution and the hash-join contract](#batch-execution-and-the-hash-join-contract);
that section is the columnar spelling of L1 and L3, and the graded layers here
are where it plugs into the algebra.

#### Every gate is a producer feeding `accumulate`

`map_add` and `map_insert` are not two different map operations; they are one
(`accumulate`) fed by two different **producers**. Clifford, rotation, and
multiply all reduce to "produce `(w, c)` terms, `accumulate_batch`, `reduce`":

| Algebra op | Producer | Term shape |
| --- | --- | --- |
| Clifford (H/S/CNOT) | pushforward along a Pauli bijection `w ↦ φ(w)` | one term per input (injective, no collision) |
| Rotation / noise | extend-by-linearity `w ↦ cos·w + sin·w'` | a small fan-out per input (branch) |
| Multiply (L4) | outer product over two operands' support | one term per `(v, w)` pair |

The producer difference lives entirely on the **term-production side** — the
`TermSink` of [The batch contract](#the-batch-contract) — not in the map, which
only ever accumulates produced batches. That is what dissolves the old
`map_add` / `map_insert` split and its `Vec<(W, C)>` staging leak.

A producer is a **monomorphized, inlinable** type, never `dyn` — this is a hot
loop, so the abstraction must compile to nothing:

```rust
pub trait TermProducer<K, C> {
    /// Push the produced terms for one existing (key, coeff) into the sink.
    fn produce<S: TermSink<K, C>>(&self, key: &K, coeff: &C, sink: &mut S);
}

/// Bijective re-key (Clifford): one produced term per input.
pub struct RekeyProducer<F> { f: F }         // ZST when the closure captures by copy

impl<K, C, F: Fn(&K, &C) -> (K, C)> TermProducer<K, C> for RekeyProducer<F> {
    #[inline(always)]
    fn produce<S: TermSink<K, C>>(&self, key: &K, coeff: &C, sink: &mut S) {
        let (k, c) = (self.f)(key, coeff);
        sink.push(k, c);
    }
}
```

Because `apply<P>` and `produce<S>` are generic, each gate call site
monomorphizes and the `#[inline]` body folds into the accumulate loop; the sink's
column is pre-sized from `Policy::capacity`, so there is no per-term allocation.
The only inherent cost is the `key.clone()` inside the closure — a stack copy for
`PauliWord`, a memcpy for `Tableau` — which a hand-written loop pays too.

#### One engine, two key types: Pauli propagation and the tableau mixture

Because the algebra is over `Key: Indexable` and a gate is a `TermProducer` that
transforms each key, **Pauli propagation and the classical stabilizer mixture are
the same `Sum` engine** — they differ only in the key type and the key's own
conjugation:

```rust
// The storage is a bare container with the algebra traits impl'd on it; the
// alias just fixes the container and bakes in the pass-through hasher.
pub type HashMapStore<K, C> = HashMap<K, C, IdentityBuildHasher>;

pub type PauliSum<C = f64, P = NoPolicy>       = Sum<HashMapStore<PauliWord, C>, P>;
pub type TableauMixture<C = f64, P = NoPolicy> = Sum<HashMapStore<Tableau,  C>, P>;
// TableauMixture replaces the former GeneralizedTableauSum<C, T, S: EntryStore>.

// A Clifford gate on ANY sum: re-key each key via its own `Clifford` impl.
impl<S, P> Clifford for Sum<S, P>
where
    S: Accumulate,
    S::Key: Clifford + Clone,      // PauliWord's is symplectic bits; Tableau's is the inverse-tableau update
    P: Policy<S::Key, S::Coeff>,
{
    fn h(&mut self, q: usize) {
        self.apply(RekeyProducer::new(|k: &S::Key, c: &S::Coeff| {
            let mut next = k.clone();
            next.h(q);             // dispatched to the key type's Clifford impl
            (next, c.clone())      // (+ any phase the conjugation puts on the coeff)
        }));
    }
    // s, cnot, cz likewise
}
```

**Correction (the sign drain the sketch glosses).** The snippet above is correct
only for a key type whose *own* `Clifford` is **phase-complete** — e.g. `Tableau`,
which tracks its `ℤ₂` sign internally, so `next.h(q)` produces the right key. For
`S::Key = PauliWord` it is **unsound as written**: the bare word's `Clifford` is
the *bit-only* blanket (its `PhaseTrack` is a no-op), so `next.h(q)` computes the
conjugated Pauli's bits but **drops the `±1` sign** — `HYH` would come back `+Y`.
So `ppvm-pauli-sum-2`'s `Clifford for Sum` does **not** dispatch to the bare key's
`Clifford`; for each key it wraps the word in a `Phased<PauliWord>` at `+1`,
conjugates via the *audited fused* phased `Clifford` (which tracks the sign),
extracts the resulting `±1` (a Clifford never emits `i` — machine-checked, the
`*_isRealPhase` lemmas), and multiplies the coefficient by that sign. The `±1`
drain is total and the re-key is a bijection (`Symplectic.*_bijective`), so no two
terms collide. Pure-sign gates `X`/`Y`/`Z` (word unchanged) take an in-place
`SignFlipByKey` pass, not a map rebuild.

For `S::Key = Tableau`, `next.h(q)` is the tableau's own phase-complete `Clifford`.
Same engine, same producer, same `accumulate` / `reduce` — only the key's
conjugation (and, for `PauliWord`, the sign-drain wrapper) varies. This is the
"smallest useful common factor" the earlier iteration deferred: the key-agnostic
graded algebra unlocked by relaxing the key bound from `Word` to `Indexable`. The
mixture's coefficient-aware `O(n^2)` measurement remains a `Tableau`-specific
algorithm (via `StabilizerFrame`), layered on top rather than merged into the
engine.

#### A third instantiation: the generalized tableau

The graded algebra also captures the *single* coefficient-aware tableau, and
finding that it does is strong validation: the current `GeneralizedTableau`
already carries a **second hand-rolled copy of the algebra** in its
`SparseVector` trait, whose `add_or_insert` / `retain` / `iter` are exactly
`Accumulate` / `Retain` / `Support`.

A `GeneralizedTableau` represents `U|c⟩`: a Clifford **frame** `U` (a `Tableau`)
times a sparse superposition `|c⟩ = ∑_b c_b |b⟩` over bitstrings. That amplitude
vector is `C[Bitstring]` — the same graded algebra, a third key type:

```rust
pub struct GeneralizedTableau<C = f64, P = CoefficientThreshold> {
    frame: Tableau,                                      // owns loss (Step 7)
    amplitudes: Sum<Vec<(Bitstring, Complex<C>)>, P>,    // C[bitstring], Vec backend
    measurement_record: Vec<Option<bool>>,
}
```

Its gate semantics fall out of the same machinery:

- **Clifford gate → the frame only.** `G·U|c⟩ = (GU)|c⟩`, so `self.frame.h(q)`
  updates the tableau and the amplitude vector is untouched — free on
  `C[Bitstring]`.
- **Non-Clifford gate → branch the amplitude `Sum`.** Read the rotation axis from
  the frame (the `O(n^2)` decomposition), then `self.amplitudes.apply(producer)`
  branches `c → cos·c + i·sin·c'` and accumulates — the identical branch-producer
  → `accumulate` → drop-below-threshold pattern `PauliSum` rotations use, down to
  the current code spilling the merge into a temporary map.
- **`coefficient_threshold`** is `CoefficientThreshold` `Policy`.

Both gate semantics are validated in
`lean/PPVM/Instantiations/Bitstring.lean`, at two different strengths. The
non-Clifford branch's key relabel `b ↦ b ⊕ s` is *derived* to be a bijection
(`xorRelabel_bijective`) that lifts to a linear isomorphism of the amplitude
module (`relabelAmp_bijective`) — an emergent theorem, so branching only *moves*
amplitude weight, never loses a term. The "Clifford touches the frame only" split
(`cliffordStep_amplitudes`, the `G·U|c⟩ = (GU)|c⟩` factorization) is instead a
*modeling* statement: it holds `rfl`-by-construction of the state as a
frame×amplitude pair, pinning down that the two factors are independent, not
proving an emergent fact.

The `Vec` backend is the right one here: the support is T-count-bounded (small),
so linear-scan `accumulate` beats hashing — the concrete motivation for the `Vec`
container. As with the mixture, the frame↔amplitude coupling and the `O(n^2)`
measurement stay `GeneralizedTableau`-specific; only the *storage and accumulate*
unify — exactly the level at which `PauliSum` fits.

There is **no `SumStorage` trait, and no owned workspace.** An earlier draft had
a `SumStorage` value bundling the map with a reusable *auxiliary map* and scratch
buffer — but that aux map was an artifact of the *old* mutation model, where a
gate iterated the live map and inserted new keys, and so needed a second map to
stage into and swap. The producer / `TermBatch` model removes that aliasing: a
`TermProducer` *reads* the map through `&` and *writes* produced terms into a
`TermSink` (the batch), then `accumulate_batch` iterates the **batch** — a
separate buffer — and merges into the map. Nothing is mutated while it is
iterated, so there is no aux map to keep. Bundling one into the storage would
double every sum's allocation and clone it on every mixture branch, for a
workspace the design no longer needs.

`apply` is therefore a thin method on `Sum` itself, over any `Accumulate`
container:

```rust
impl<S: Accumulate, P: Policy<S::Key, S::Coeff>> Sum<S, P> {
    /// Produce terms into a batch, merge, canonicalize, truncate. `TP` is a
    /// type parameter, never `dyn`, so this monomorphizes and the producer's
    /// `#[inline]` `produce` folds into the loop — no per-term allocation.
    fn apply<TP: TermProducer<S::Key, S::Coeff>>(&mut self, producer: TP) {
        let mut batch = /* transient, or a driver-owned reusable buffer */;
        for (k, c) in self.storage.iter() { producer.produce(&k, &c, &mut batch); }
        self.storage.accumulate_batch(&batch);
        self.storage.reduce();
        self.policy.truncate(&mut self.storage);
    }
}
```

The only staging that remains is the `TermBatch` the producer fills. It is
**transient by default** (allocated per `apply`), which is exactly right for a
frequently-cloned, small `GeneralizedTableau`. A long-lived driver that applies
millions of gates — `PauliSum` — may instead own **one** reusable batch (and, for
a whole-map Clifford rewrite, a second map to swap into) to amortize allocation.
That reuse is an **opt-in optimization of the driver**, never a field of `Sum`
and never imposed on the tableau side.

Because the batch a producer stages is exactly the structure-of-arrays layout
prefetched, SIMD, threaded, and offloaded backends consume, the execution
schedule for that batch — group prefetching, radix partitioning, device
offload — is specified in
[Batch execution and the hash-join contract](#batch-execution-and-the-hash-join-contract)
below.

#### The pass-through storage contract

Because a key's [`key_hash()`](#indexable-values) is already the finalized,
avalanche-quality digest, the provided storage aliases **bake in** an identity
pass-through hasher so the map consumes the digest directly:

```rust
#[derive(Default, Clone)]
pub struct IdentityHasher(u64);
impl std::hash::Hasher for IdentityHasher {
    fn write_u64(&mut self, n: u64) { self.0 = n; }     // store the digest
    fn write(&mut self, _: &[u8]) { unreachable!() }     // keys only write a u64
    fn finish(&self) -> u64 { self.0 }                   // hand it back verbatim
}
```

The `HashMapStore` / `DashMapStore` aliases therefore fix their `HashMap` /
`DashMap` to `IdentityBuildHasher`, so `finish() == key.key_hash()` and the
digest reaches hashbrown untouched. This is part of the storage contract, not a
user responsibility: the user selects the key's *internal* digest algorithm (the
`H` in `PauliWord<A, H>`), which governs distribution quality, and never selects
a map hasher — the direct-digest model leaves none to choose. Re-hashing the
digest with a general hasher would be wasted work and could re-correlate the
low bits the key's finalization fold just decorrelated.

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

Every keyed store must consume its key's `key_hash()` digest directly, through
the identity pass-through described in
[The pass-through storage contract](#the-pass-through-storage-contract).

### Compatibility with current names

The redesign is not a vocabulary reset. The following names are retained or
changed according to whether their underlying responsibility changes:

| Current implementation | Proposal | Rationale |
| --- | --- | --- |
| `Config` | removed | The bundle itself is removed; this is not a rename. |
| `PauliWordTrait` | split into `Word` (read-only inspection), `Indexable` (key), `PauliBits` (mutation) | The old bundle is decomposed by concern; bit-level mutation is a narrow Pauli-specific trait, hashing is separate, inspection is algebra-agnostic. |
| `n_qubits`, `get`, `set`, `weight` | `n_sites`, `get`, `weight` on `Word`; `set` removed | `Word` is read-only; positional mutation (`set`) moves to `PauliBits`, since it is ill-defined for ordered algebras. |
| bit accessors `get_xbit`/`set_xbit`/… | `PauliBits::x_bit`/`set_x_bit`/… | The rotation hot path is sub-site; it keeps a dedicated Pauli trait rather than the generic `Word`. |
| word-level Clifford (blanket over `PauliWordTrait`) | `Clifford` blanket over `SymplecticColumns` + `PhaseTrack` + `BlanketClifford` (opt-in marker) | The symplectic sign logic is written once and shared by the phaseless words and `Tableau`; `PhasedPauliWord` opts *out* and supplies a read-once fused `impl Clifford` instead (avoids the blanket's double bit read). |
| (new) | `SymplecticColumns`, `PhaseTrack`, `BlanketClifford`, `StabilizerFrame` | The symplectic-bits + phase-extension decomposition (role-independent column algebra, role-dependent phase, role-exclusive frame ops), plus the empty `BlanketClifford` marker that selects the shared blanket so a fused override stays coherence-legal. |
| concrete `PauliWord` | `PauliWord` | The packed X/Z word is the same domain concept. |
| concrete `LossyPauliWord` | `LossyPauliWord` | The packed X/Z/loss representation remains concrete and flattened. |
| `PhasedPauliWord` | `PhasedPauliWord` alias over non-indexable `Phased` | The wrapper is generic over ordinary and lossy words but is not a production map key. |
| `rehash` | private cache invalidation | Recalculation changes from eager mutation-time work to lazy demand-time work without exposing cache mechanics through `Indexable`. |
| `Strategy` | `Policy` | Intentional terminology change requested for this redesign; the `Copy` bound is dropped. |
| `Coefficient::cutoff` | removed | The truncation predicate moves to `Policy`; `Coefficient::magnitude` exposes only the value property the policy thresholds. |
| `Coefficient::sin_cos` | `Angle<C>` | The rotation angle becomes a separate domain, defaulting to the coefficient type. |
| `Coefficient: Mul<f64>` | removed | Vestigial once `sin_cos` returns coefficient-domain amplitudes; it was the sole bound excluding exact rings, so dropping it lets `GaussianInt` / cyclotomic coefficients be `Coefficient`s. |
| `Coefficient::half` | `Halvable<C>` capability | `0.5·x` is partial on exact rings (`0.5·(1+i)` leaves `ℤ[i]`) — the same escape `Mul<f64>` had — and is needed only by the `(I ± Z)/2` measurement projector, so it splits into a capability rather than re-foreclosing exact `Coefficient`s. |
| L4 `Multiply` bound `ComplexCoefficient` | `ImaginaryUnit` (primitive fourth root of unity) | The Lean proof shows the twisted product needs only `i⁴ = 1`, so the bound is loosened from `Complex<f64>` to admit exact rings (`lean/PPVM/Algebra/Twisted.lean`, `lean/PPVM/Pauli/Matrix.lean`). |
| `ACMap` (`ACMapBase`/`Iter`/`AddAssign`/`Insert`/`Retain`/`Consume`) *and* `SparseVector` | graded `Support` / `Accumulate` / `Scale` / `Pair` / `Multiply`, `impl`'d directly on `Vec`/`HashMap` | Two hand-rolled copies of the same abstraction collapse into one over `C[K]`; `Retain` (truncation) leaves the algebra for `Policy`. |
| `PauliSum::map_insert`, `map_add` | one `Sum::apply(producer)` → `accumulate_batch` + `reduce` | They differed only in their producer, not the map op; the `Vec<(W,C)>` staging leak is removed. |
| `PauliSum` aux-map + `scratch` fields | removed | The aux map was an artifact of mutate-while-iterate; the producer/`TermBatch` model needs no owned workspace, so `Sum` owns none. |
| `PauliSum` | `PauliSum` over generalized `Sum` machinery | Pauli-facing code keeps its established name; `Sum` names the new cross-algebra engine. |
| `GeneralizedTableauSum<C, T, S: EntryStore>` | `TableauMixture = Sum<HashMapStore<Tableau, C>, P>` | The mixture is `C[Tableau]` — the same graded algebra as `PauliSum`, keyed on `Tableau`; storage unifies while the `O(n^2)` measurement stays tableau-specific. |
| `GeneralizedTableau.coefficients: SparseVector` | `Sum<Vec<(Bitstring, Complex<C>)>, CoefficientThreshold>` | The amplitude vector is `C[Bitstring]` — the same algebra, `Vec` backend; see [A third instantiation](#a-third-instantiation-the-generalized-tableau). |
| `EntryStore`, `VecStorage`, `MapStorage`, `SparseVector` backings | the container itself (`Vec`, `HashMap`, `ColumnStore`) | No storage-wrapper trait; the algebra traits are `impl`'d on the collection directly. |
| `Config::BuildHasher` | removed; `Indexable::key_hash() -> u64` | The direct-digest model leaves the map no hasher to choose; the key's internal algorithm is a private representation parameter, and the finalized digest is exposed as a value. |
| `HashFinalize` | retained, but private to `ppvm-pauli-word` | The per-algorithm/per-width finalization fold still runs inside `key_hash()`; it leaves the algebra-agnostic contract. |
| (new) | `IdentityBuildHasher` + pass-through storage aliases | The provided `HashMapStore` alias consumes `key_hash()` directly instead of re-hashing it. |
| `PauliStorage` | removed | Packed backing storage becomes private to the concrete word representation. |
| (new) | `Conjugate` + `Pair::hermitian_overlap` | The sesquilinear state/amplitude inner product `∑ conj(a_k) b_k` (needed by `GeneralizedTableau`'s complex `C[Bitstring]`), split from the symmetric bilinear `overlap`; the conjugation is a separate `*`-ring capability, not a `Coefficient` bound. Machine-checked in `lean/PPVM/Algebra/GradedMap.lean`. |
| (new) | `Columnar`, `KeyColumn`, `KeyBatch`, `TermBatch`, `TermSink`/`TermProducer` | Structure-of-arrays term types; the columnar spelling of `Accumulate`/`Pair`. See [Batch execution and the hash-join contract](#batch-execution-and-the-hash-join-contract). |
| (new) | `ColumnStore` (SoA) backend | The one storage that must be a new struct (SoA planes); `Vec`/`HashMap` need none. Expressible from day one because no signature leaks AoS. |

Names such as `Word`, `Indexable`, `Accumulate`, and `Sum` are
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
- `Word`, consumed by display, serialization, tests, and the sparse-sum
  plumbing as a read-only inspection interface (not the propagation interface);
- `PauliBits`, consumed by the rotation kernels and implemented by `PauliWord`
  and `LossyPauliWord`;
- `SymplecticColumns` and `PhaseTrack`, consumed by the blanket `Clifford` and
  implemented by the phaseless words (`PauliWord`/`LossyPauliWord`) and `Tableau`,
  which opt into the blanket via the `BlanketClifford` marker; `PhasedPauliWord`
  instead carries its own read-once fused `impl Clifford`;
- `StabilizerFrame`, the role-exclusive tableau operations;
- gate and noise traits, implemented across propagation and tableau backends;
- the graded map layers `Support` / `Accumulate` / `Scale` / `Pair` / `Multiply`,
  each admitted by a distinct algebraic property *and* a distinct consumer, and
  `impl`'d directly on `Vec<(K,C)>` (small support), `HashMap` (large), and
  `ColumnStore` (SIMD/GPU);
- `Conjugate`, the `*`-ring involution consumed by the sesquilinear
  `Pair::hermitian_overlap` and implemented by every complex/exact coefficient
  ring (`Complex<f64>`, `GaussianInt`, cyclotomic integers) — a separate
  capability like `ImaginaryUnit`, not a `Coefficient` bound;
- `Policy`, implemented by independent capacity and truncation behaviors, and
  consuming the non-algebraic `Retain` capability rather than the map algebra.

It rejects the removed global `Config`, `PauliSumAlgorithm`,
`TableauMixtureAlgorithm`, and `TableauStorage` traits, as well as word
subtraits named `PauliWord`, `FermionWord`, or `LossyPauliWord`. Their
alphabet distinctions are expressed by `Word::Site` or by concrete types instead
of one-alphabet subtraits; the concrete `PauliWord` and `LossyPauliWord` type
names remain available. Note the contrast with `PauliBits`: it is admitted not
as an alphabet subtrait but as a narrow *bit-mutation* capability that genuinely
has multiple implementers and generic consumers.

### Sparse-sum branch staging

A propagation rule can turn one term into multiple terms. For example, a
Pauli rotation may produce:

```text
c P -> c cos(theta) P + c sin(theta) P'
```

This branch is validated in `lean/PPVM/Instantiations/Rotation.lean`: when `G`
and `P` anticommute the new key `P' = iGP` has symplectic bits `G ⊕ P` and is
genuinely distinct from both operands (`anticommute_new_key`), so exactly one
fresh term is produced; and the `(c_P, c_{P'})` coefficient update is a 2-D
rotation that is norm-preserving (`rot_norm_sq`), reversible (`rot_neg_rot`), and
angle-additive (`rot_rot`) — the last being the identity a rotation-merging
Trotter optimization relies on.

A `TermProducer` *reads* the live map through `&` and *writes* the produced
terms into a `TermSink` — a separate buffer, never the map itself — so there is
no mutate-while-iterate hazard and **no auxiliary map is needed.** After the
producer finishes, `accumulate_batch` merges that buffer into the map, combining
colliding coefficients. This is a deliberate simplification over the old model,
which iterated the live map and inserted into it, and so needed a second map to
stage into and swap.

The only staging that remains is the single `TermSink` buffer:

```text
old model:  iterate active map, insert new keys -> needs aux map + swap
new model:  read map -> produce into TermSink buffer -> accumulate_batch into map
```

The buffer is **transient by default**, allocated per operation — right for a
small, frequently-cloned `GeneralizedTableau`. A long-lived driver (`PauliSum`)
may own **one** reusable buffer (and, for a whole-map Clifford rewrite, a second
map to swap into) purely to amortize allocation across millions of gates. That
reuse is an opt-in optimization of the driver, not a field of `Sum` and not part
of the storage contract. Whether it is worthwhile is a benchmark decision.

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

`Sum::apply` and the `Accumulate`/`Pair` layers describe *what* the merge
computes. But a producer that emitted one `(key, coeff)` at a time would expose
no batch to prefetch, no homogeneous run to vectorize, no partition to distribute
across threads, and no bulk kernel to hand to a device. That is why a
`TermProducer` stages into a `TermSink` and the map consumes an
`accumulate_batch`: the *layout* the merge consumes is a batch, so those backends
can be written against it.

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
    /// loop. Its obligation is now an equation, not a phrase: `out[i]` must
    /// equal the `i`-th key's `Indexable::key_hash()`.
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

The bulk map operations that consume these columns are not a new trait: they are
exactly the columnar methods of the graded algebra layers —
[`Accumulate::accumulate_batch`](#the-map-is-a-graded-algebra-over-cw) is the
build side and [`Pair::probe_batch`](#the-map-is-a-graded-algebra-over-cw) is the
read side. Restated here at the column level for reference:

```rust
// Accumulate::accumulate_batch — the build/probe-with-group-by of a hash join.
// The implementation owns whether it runs scalar, group- or pipeline-prefetched,
// SIMD-vectorized, hash-partitioned across threads, or offloaded; the contract
// fixes only the result and the columnar layout.
fn accumulate_batch(&mut self, batch: &TermBatch<W, C>);

// Pair::probe_batch — read-only probe of a key column, for the overlap and
// expectation paths. Takes a `KeyBatch` so the precomputed hash column drives
// the prefetch with no coefficient column in the working set.
fn probe_batch(&self, keys: &KeyBatch<W>, out: &mut [Option<C>]);
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

`accumulate_batch` is where the branch-staging merge lands: the reusable scratch
buffer of the previous section becomes a `TermBatch` — the probe side of the
join. A naive backend may still have its `TermSink` collect into a scalar `Vec`
and loop, but generic hot paths build a columnar `TermBatch` and hand it to
`accumulate_batch`.

Gate and rotation traits change accordingly, and this is the interface change
users feel first: term *production* is separated from term *insertion*. A
`TermProducer` (a rotation, a Clifford re-key, a multiply operand) appends
produced terms into a `TermSink` — filling the key column and the coefficient
column as it goes — instead of mutating the map through a callback. That decoupling, plus the columnar layout, is what lets the produced
batch be prefetched, vectorized, partitioned, or shipped to a device before it
ever touches the table, and it keeps the propagation rule (the physics)
independent of the merge strategy (the systems concern).

### What the batch contract asks of the other contracts

- `Indexable` is unchanged and stays minimal. The batch surface lives on the
  separate `Columnar` trait, which a key type implements only when it is used in
  a batch. A packed Pauli word's column is its X and Z plane blocks and the
  flattened `LossyPauliWord` adds a loss plane; `Phased<W>` is not indexable, so
  it is not `Columnar` and never appears in a batch.
- `Word`: unaffected — it is read-only inspection either way. Per-value mutation
  happens through `PauliBits` on a single word; columns are built by appending
  produced keys, not by mutating in place.
- `Policy`: truncation already operates on the whole map and is naturally bulk;
  it should be expressible as a batch retain (via `KeyColumn::gather`) so it
  composes with a partitioned or offloaded table.
- `Sum::apply`: the producer's `TermSink` *is* the join's probe-side buffer, and
  `apply` desugars to `accumulate_batch` + `reduce`. The buffer is transient by
  default; a driver may reuse one to amortize — a benchmark decision, unchanged
  from above.

### Parallel and offloaded backends

Recognizing the operation as a hash join fixes the parallel and GPU stories
without further interface changes. A partitioned hash join radix-partitions both
sides by the high bits of the key hash so each partition is a disjoint
sub-table; partitions then merge independently — one per thread, or one per GPU
block — with no cross-partition synchronization. The batch contract is precisely
what a partitioner consumes: `accumulate_batch` may internally `gather` a
`TermBatch` into per-partition batches and run them concurrently, and a device
backend may copy a `TermBatch` across the host/device boundary and run the probe
as a kernel. None of this is visible in the contract; all of it is precluded by
a scalar one. This is the partition-then-merge shape the `DashMap` backend
already gestures at, now stated as the contract that makes such backends
interchangeable rather than separate code paths.

## Indexable values

Hash-enabled `Word` values and `Tableau` can both be expensive, mostly-stable
map keys. Their hashing contract should be expressed independently of any map.

`Indexable` is **not** the universal key bound, though. The universal
requirement — what the graded algebra's `type Key` bounds on, and all any
container needs to *find* and *store* a key — is just `Eq + Clone`: a linear-scan
`Vec` backend compares by equality and never hashes. Hashing, and therefore the
digest, is a property only a *hash* backend needs. So `Indexable` is required
solely on the `HashMap` / `ColumnStore` `Accumulate` impls (and on `Columnar`),
never on `type Key` itself. A `Bitstring` used only in a `Vec`-backed
`GeneralizedTableau` is `Eq + Clone` and provides no `key_hash()` it would never
use; `PauliWord` / `LossyPauliWord` / `Tableau`, used in hash backends, are
`Indexable`.

With that scoping, `Indexable` makes the one load-bearing value — the finalized
structural digest — first class:

```rust
pub trait Indexable: Clone + Eq + Hash {
    /// The finalized structural digest of this key: avalanche-quality in both
    /// the low bits (the hashbrown bucket) and the top 7 (the control tag),
    /// so it can be consumed *directly* as the map hash. Contracts:
    ///
    ///   * `Hash for Self` is exactly `state.write_u64(self.key_hash())`;
    ///   * structurally equal keys return equal digests; and
    ///   * `KeyColumn::hash_into` reproduces this value bit for bit.
    ///
    /// This exposes the digest *value*, not the cache mechanics — there is no
    /// cache type or invalidation hook in the contract.
    fn key_hash(&self) -> u64;
}
```

The important points are:

- the digest is used **directly** as the map hash — the map does not re-hash it
  (see [the pass-through storage contract](#the-pass-through-storage-contract)),
  so the value `key_hash()` returns must already be well distributed;
- the *choice of internal hash algorithm* (`fxhash`, `gxhash`, …) is a private
  representation parameter of the concrete key (the `H` in `PauliWord<A, H>`),
  not an associated type on `Indexable`. There is no per-key `BuildHasher` on
  the contract, because the direct-digest model leaves the map no hasher to
  pick;
- avalanche quality is the key's responsibility: a weak-but-fast algorithm on a
  short key must apply a finalization fold inside `key_hash()` (see
  [Concrete word hashing](#concrete-word-hashing));
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
cache their digest through interior mutability. `Hash::hash` may populate a cache
through shared access, while structural mutators clear affected cells through
their exclusive `&mut self` access. This preserves `Send + Sync` for the shipped
representations. `Indexable` itself does not require either concurrency bound.

The *mechanism* is a private representation choice. `PauliWord` realizes it with a
sentinel `AtomicU64` (a distinguished `HASH_UNCACHED` value = "not yet computed";
`key_hash()` does a relaxed load, computes-and-relaxed-stores on a miss, and
mutators reset it) rather than the `OnceLock<u64>` this section originally
sketched: the digest is a pure function of the immutable content, so a racing
miss just recomputes *the same* value — a relaxed atomic is correct and avoids
`Once`'s CAS init path, which measurably dominated the `PauliSum` Clifford re-key
hot loop (every freshly built key hits the cold init once). A key whose true
digest equals the sentinel is simply recomputed each time — a `1`-in-`2⁶⁴`
perf non-event with no correctness effect. Both realizations satisfy the same
contract; the choice is invisible through `Indexable`.

### Key mutation invariant

An indexable value must not be structurally mutated while it is stored as a map
key. Cache invalidation makes a value correct for its next insertion or lookup;
it cannot make in-place mutation of an existing map key valid.

This invariant is not left entirely to discipline: the sparse-sum engine
enforces its structural half through the producer interface. A `TermProducer`
reads live keys through a *shared* `&Word` borrow and *emits* new `(Word, Coeff)`
terms into the `TermSink`; it is never handed `&mut Word` for a stored key. So
the propagation kernels physically cannot mutate a key in place — they clone,
mutate the clone through `PauliBits`, and the engine accumulates it under the
updated digest. That the producer only ever *reads* keys and *produces* terms is
a deliberate safety choice, not an incidental one.

Structural fields should nonetheless be private in the new representations, and
any remaining mutation must go through operations that invalidate the affected
cache, or through mutation guards that invalidate on completion.

## Concrete word hashing

Every concrete `Word` owns its internal hash algorithm, private cache
representation, and invalidation logic, and produces the finalized digest that
`Indexable::key_hash` returns. Pauli words hash their X/Z content, lossy words
compose Pauli and loss components, and future fermion words hash their ordered
factors. Factor order is part of fermionic identity.

The trait layer does not expose packed storage to support hashing. Concrete
implementations hash their private fields and apply a **finalization fold**
internally before caching, so `key_hash()` is avalanche-quality even for a short
key consumed directly as the hashbrown bucket-and-tag. That fold is per-algorithm
and per-width — a couple of `fxhash` rounds on an `[u8; 8]` word leave the low
bits correlated and need `raw ^ (raw >> 32)`, whereas AES-based `gxhash`
avalanches an 8-byte key already and folds nothing. The retained `HashFinalize`
helper that encodes those rules is a **private utility inside `ppvm-pauli-word`**,
not part of the algebra-agnostic `Indexable` contract; the public contract is
only "`key_hash()` is well distributed," checked by a distribution property test
rather than the type system. Detailed layouts and component invalidation rules
are in [`word-data-structures.md`](word-data-structures.md).

## Tableau indexability

A tableau may itself be used as a key by a classical-mixture algorithm, so the
concrete tableau implements `Indexable` directly, owning its internal hash
algorithm and cache representation and returning its own finalized
`key_hash()`. This does not imply that a tableau is a `Word`; they only share
the `Indexable` key capability.

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
Sum<Container, Policy>
// PauliSum<f64>          = Sum<HashMapStore<PauliWord, f64>, P>        // C[PauliWord]
// TableauMixture<f64>    = Sum<HashMapStore<Tableau,   f64>, P>        // C[Tableau]
// GeneralizedTableau amps = Sum<Vec<(Bitstring, Complex<f64>)>, CoefficientThreshold>  // C[Bitstring]
Tableau
```

All three are `Sum` over the same graded algebra, differing only in the key type
and the container backend. Domain aliases preserve `PauliSum` and introduce
`FermionSum` without a monolithic configuration trait.

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
