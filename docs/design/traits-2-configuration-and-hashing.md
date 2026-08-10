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
    ///
    /// Law — `magnitude` must be an **absolute value** on the coefficient ring:
    /// `N x ≥ 0`, `N x == 0 ⟺ x == 0`, `N (x + y) ≤ N x + N y`, and
    /// `N (x · y) == N x · N y`. This is not decoration: it is exactly what the
    /// `ℓ¹` truncation bound consumes.
    fn magnitude(&self) -> f64;
}
```

`magnitude` carries a **law**, and it is the one law the whole truncation
guarantee rests on: it must be an absolute value on the coefficient ring
(nonnegative, zero only at `0`, subadditive, multiplicative). Nonnegativity alone
is not enough — `N x = x²` is nonnegative, zero only at `0`, and even
multiplicative, yet the `ℓ¹` bound fails for it
(`lean/PPVM/Algebra/Truncation.lean`, `l1_bound_needs_subadditive`). The bound
itself is proved for an arbitrary coefficient ring carrying such an `N`
(`l1_bound_abv`), with the shipped `Complex<f64>`/`magnitude() = norm()`
configuration covered by the normed-field instance (`l1_bound_norm`,
`l1_bound_complex`).

The law is stated at its *sufficient* strength, not its necessary one, and one
shipped coefficient ring cannot meet it: **no** absolute value `N : ℝ[sᵢ, cᵢ] → ℝ`
exists on `ppvm-sym-2`'s symbolic `Term`, because the natural `ℓ¹` coefficient
norm is only *sub*-multiplicative (`(1+x)(1−x) = 1−x²` gives `2·2 = 4` against
`2`). `lean/PPVM/Algebra/Truncation.lean` adjudicates the resulting choice:

* `l1_bound_seminorm` — the bound **survives** weakening `AbsoluteValue` to a
  nonnegative, `0`-vanishing, subadditive, *sub*-multiplicative seminorm
  (`seminorm_weaker_than_abv` witnesses that this is strictly weaker). So an `ℓ¹`
  `magnitude` on the symbolic ring would carry the full error guarantee; only
  behaviour parity with old (which never truncated a symbolic coefficient)
  weighs against it.
* `l1_bound_seminorm_needs_zero` — the surviving clause is `N 0 = 0`, and it is
  exactly the one `ppvm-sym-2`'s parity-preserving `magnitude() = f64::INFINITY`
  breaks. That choice is therefore a documented **law exemption** carrying no
  `ℓ¹` guarantee (`CoefficientThreshold` is inert on symbolic coefficients, which
  is old's behaviour), not an approximation of one.

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

#### The symbolic coefficient ring is *free*, so the rotation laws need `eval`

`ppvm-sym-2` instantiates both domains with the same type: `impl Angle<Term> for
Term`. Its coefficient ring is the **free** polynomial ring `ℝ[sᵢ, cᵢ]`, modelled
in `lean/PPVM/Instantiations/Symbolic.lean` as `AddMonoidAlgebra ℝ (ℕ →₀ ℕ × ℕ)`
— the same `CMap` object as the sum itself, one level down. Two consequences are
machine-checked there and neither is optional reading for that instantiation:

* **`sin² + cos² = 1` does not hold** (`pythagorean_ne_one`): `sin(x).square() +
  cos(x).square()` is a genuine two-monomial `Sum` and no `Term`-level operation
  reduces it. The rotation guarantees of
  `lean/PPVM/Instantiations/Rotation.lean` (`rot_norm_sq`, `rot_rot`) consume
  exactly that relation, so they transfer to the symbolic domain **only after
  evaluation**, pointwise in `θ`: `evalHom_symRot` is the commuting square
  (`ev ∘ symRot = rot θ ∘ ev`), `symRot_norm_sq_after_eval` the transferred norm
  preservation, and `symRot_norm_sq_ne_symbolically` the witness that the
  unqualified claim is false.
* **The sine degree is a grading** (`sinDeg_add`), so the `max_sin` cutoff spans a
  monomial *ideal* (`truncIdeal_mul_right`). That is what licenses the two
  constructs the symbolic coefficient is built around: drop-at-accumulate inside
  `Sum::add_term` is *exact*, equal to truncating the finished product
  (`mulMono_drop_at_insert_eq_drop_at_end`, via
  `GradedMap.batchMap_filter_key` — a key-only keep-rule is additive, unlike
  `retain` in general), and `Sum::mul_term`'s whole-table `clear()` shortcut is
  sound (`mulMono_clear_sound`). The companion negative result
  `eps_drop_at_insert_ne_drop_at_end` shows the *coefficient*-magnitude axis
  (`min_eps`) is **not** interchangeable with a post-pass, so unlike `max_sin` it
  may not be relocated out of the accumulation loop.
* **The `min_eps` arm of the same `clear()` shortcut is *not* the degree arm.**
  `Sum::mul_term` also clears the whole table when `|coeff| < min_eps`, and that
  arm is an **over-truncation**, not an equality: `epsClear_ne_retain_pointwise`
  exhibits a table (one entry of magnitude `10⁶`, multiplier `10⁻¹³`,
  `min_eps = 10⁻¹²`) whose product monomial the per-monomial rule in
  `Sum::add_term` keeps and the shortcut discards. What licenses keeping the
  shortcut anyway is the `ℓ¹` bound `epsClear_l1_eq` / `epsClear_l1_lt` and its
  read-out corollary `epsClear_error_lt` (stated against
  `PPVM.Truncation.l1_bound`): the discarded mass is exactly `|c|·ℓ¹(A)`, hence
  strictly under `min_eps·ℓ¹(A)`. Citing only `mulMono_clear_sound` for the whole
  shortcut would over-claim.
* **The implemented ring is `ℤ/4`-graded, and its complex evaluation is not
  injective.** `Prod` stores a phase byte that is part of its `Hash`/`Eq`, so the
  ring is `PhasedSymRing = AddMonoidAlgebra ℝ (Mono × ZMod 4)`, not `SymRing`.
  `Term::eval_complex` is the `ℝ`-algebra hom `evalC` into `ℂ` (`evalC_mul`,
  built from `Twisted.iPow_add` on the grading and `monoValue_add` on the
  exponents), and `evalC_not_injective` shows its kernel is non-trivial. That is
  the machine-checked content of the `ImaginaryUnit` law exemption
  `ppvm-sym-2/src/coeff.rs` documents: `i·i` is the key `phase 2`
  (`iSym_sq_ne_neg_one`) while `−one()` is `−1` on `phase 0`, denotationally
  equal (`evalC_iSym_sq_eq_neg_one`) but distinct hash keys. The same
  non-injectivity is why symbolic truncation is *coarser* than truncation on the
  values: `phaseTwo_cancel_ne_zero` gives two summands that cancel in `ℂ` yet
  occupy different keys, so `min_eps` thresholds them independently.
  `Conjugate for Term` is the phase-negating ring involution `conjSym`
  (`conjSym_conjSym`), correct because `evalC ∘ conj = star ∘ evalC`
  (`evalC_conjSym`, whence `conj i = −i`) — the ring-level form of
  `Pauli/Matrix.lean`'s `star_iU`, which on its own is a fact about a 2×2 matrix,
  not about this ring.
* **`mul_phase` is a key relabelling, and that relabelling *is* multiplication by
  `iᵏ`.** `Term::mul_phase k` touches no coefficient: it adds `k` to every
  monomial's phase byte. `phaseFold_eq_iSym_pow_mul` proves that the relabelling
  equals the ring product `iᵏ · x` in `PhasedSymRing`, with the read-out corollary
  `evalC_phaseFold` (`evalC θ (mul_phase k x) = iᵏ · evalC θ x`). This is what
  turns "phase *every* summand, the constant one included" into an identity rather
  than a plausible choice — `phaseFold_const` is the `(0, 0) ↦ (0, k)` arm the new
  `Sum::add_term` keeps out of its `c₀` short-circuit, and
  `phaseFold_drop_const_ne` shows old's behaviour (leave the constant on key
  `(0, 0)`) computes a *different function*, not another representation of the same
  one. It is the machine-checked justification for the one deliberate behaviour
  divergence on this path (`oldSuspectedBugs` #3); the previously cited
  `Twisted.twistedConv_add_left`/`_right` are statements about the twisted
  convolution on *Pauli* keys and do not by themselves cover the symbolic fold.
* **`max_sin` is a property of the representation, not of the value.** The
  truncation theorems above describe the *map-backed* accumulation only —
  `Sum::add_term`/`Sum::mul_term` are the only sites that consult the bound. The
  shipped `Term` has three non-map fast arms (`One × One`, `Const × One`,
  `Const × Const`; perf feature 1) that never consult it, so the implemented
  product does **not** factor through the denotation: `mulImpl_not_wellDefined`
  exhibits `den a₁ = den a₂` with `den (mulImpl 2 a₁ b) ≠ den (mulImpl 2 a₂ b)`
  (`a₁ = One(sin(x₀)², 1)`, `a₂` the map-backed `Sum` denoting the same
  polynomial, `b = One(sin(x₁), 1)`). The positive half is
  `mulImpl_one_one_untruncated`: the fast arm computes the *untruncated* ring
  product for every `k`, hence `fastArm_escapes_bound` — a coefficient that stays a
  single monomial escapes `set_max_sin` without limit. Consequences for citation
  hygiene: `mulMono_drop_at_insert_eq_drop_at_end` must **not** be read as an
  end-to-end guarantee that the propagated coefficient equals the truncated ring
  product, and `set_max_sin` is not a hard degree bound on the result. Consequence
  for maintenance: unifying the four `Inner` arms onto one map-backed
  representation would make the product well-defined and therefore change numbers
  — a spec violation, not a tidy-up.

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
    fn rotate_1(&mut self, axis: Pauli, qubit: usize, theta: A);

    fn rx(&mut self, qubit: usize, theta: A) {
        self.rotate_1(Pauli::X, qubit, theta)
    }
    // ry/rz likewise; rx_many/ry_many/rz_many loop over them (`where A: Clone`)
}

// Stochastic operations take the randomness as a parameter; see
// "Where the randomness lives" below.
pub trait PauliError<C: Coefficient> {
    fn pauli_error<R: rand::Rng + ?Sized>(
        &mut self,
        qubit: usize,
        probabilities: [C; 3],
        rng: &mut R,
    );
}
```

The sketch above names only the traits whose *shape* changes. The rest of the
current behavioral surface crosses over unchanged and is listed in
[Compatibility with current names](#compatibility-with-current-names):
`CliffordExtensions` and the batched `CliffordBatch` /
`CliffordExtensionsBatch`; `Reset` with its `reset_x` = `reset` + `h` /
`reset_y` = `reset` + `h` + `s` defaults; the rest of `RotationOne` (the required
axis-generic `rotate_1` the `rx`/`ry`/`rz` defaults and the `*_many` batch loops
sit on top of) together with the other gate traits `RotationTwo`, `RotXY`,
`CRx`, `U3Gate`, `TGate`, and `Projection`; `Trace`; the stim
`x_error`/`y_error`/`z_error` aliases and the `*_many` forms on `PauliError`;
and the channel family
(`PauliErrorAll`, `TwoQubitPauliError`, `Depolarizing`, `Depolarizing2`,
`AmplitudeDamping`, `LossChannel`, `CorrelatedLossChannel`, `ResetLossChannel`,
`AsymmetricLossChannel`). Two global edits apply to them, and nothing else:

1. The `Config` bundle parameter becomes the coefficient type itself, and an
   operation with no numeric parameter (like `ResetLossChannel`) carries none.
2. Every **stochastic** method takes the generator as a parameter,
   `rng: &mut R` (`R: rand::Rng + ?Sized`), rather than drawing from state the
   simulator owns — see [Where the randomness lives](#where-the-randomness-lives).
   The full inventory across this trait surface: `Measure` (both methods), all
   eight `Reset` methods, all eight `PauliError` methods, and every method of
   `PauliErrorAll`, `TwoQubitPauliError`, `Depolarizing`, `Depolarizing2`,
   `LossChannel`, `CorrelatedLossChannel` and `AsymmetricLossChannel`.

   The edit is deliberately **not** uniform across the channel family: two of its
   members draw no randomness and therefore take no `rng` —
   `AmplitudeDamping::amplitude_damping` (a deterministic Kraus rescaling of the
   coefficients) and `ResetLossChannel::reset_loss_channel` (an unconditional
   clear of the loss bit). "Is it a channel?" is the wrong test for this
   parameter; "does it draw?" is the right one.

Default bodies are otherwise reproduced call-for-call — `reset_x` is still
`reset` then `h` — they simply thread `rng` through to the required method.

A unital Pauli channel acts diagonally in the Pauli basis, `P ↦ λ_P·P`, and its
transfer eigenvalue collapses (using `Σ_Q p_Q = 1`) to
`λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q`, where anticommutation is the
symplectic form `ω(P,Q) = 1`. That algebraic form is machine-checked in
`lean/PPVM/Algebra/Noise.lean` (`pauli_channel_eigenvalue`, and
`pauli_channel_eigenvalue_omega` tying `anti` to `PPVM.Symplectic.omega`); the
zero-state read-out `⟨0ⁿ|ρ|0ⁿ⟩ = Σ_{P X-free} c_P` is `overlap_with_zero_xfree`.

The eigenvalue is also **contractive**, which is what licenses the diagonal
channel's `ScaleByKey` fast path skipping *both* a truncation pass and a
Pauli-weight re-check: for a sub-stochastic `[p_X, p_Y, p_Z]` (`p ≥ 0`,
`Σ p ≤ 1`), `|λ_P| ≤ 1` (`pauli_channel_eigenvalue_abs_le_one`, and
`pauli_channel_eigenvalue_omega_abs_le_one` on the symplectic form), so the
channel is an `ℓ¹` contraction (`l1_contractive`) — the hypothesis that makes
`PPVM.Truncation.l1_bound` compose across a *noisy* propagation rather than only
bound one truncation — and it never introduces or moves a key
(`scaleByKey_support_subset`), so no surviving term's weight can have grown.
Sub-stochasticity is a real precondition on the channel constructors, not a
formality: `eigenvalue_abs_le_one_needs_substochastic` exhibits an
over-normalized probability vector that breaks the bound.

That eigenvalue formula says *which* `p_Q` a factor must sum, and
`two_qubit_pauli_error` is where the port could get it wrong silently: the old
crate hard-codes, for each of the 16 observed pairs, a hand-written list of eight
indices into `p: [Coeff; 15]` (`crates/ppvm-pauli-sum/src/sum/noise.rs:50-104`)
with no derivation in the source, and the shipped tests only probe one-hot
probability vectors, on which a transposed index is invisible.
`twoQubitPauliError_indices_anticommuting` checks all fifteen lists against
genuine two-qubit anticommutation in the crate's documented probability order
`{IX, IY, IZ, XI, …, ZZ}` (`crates/ppvm-traits/src/traits/noise.rs:73`), with
`…_length`/`…_nodup` making `contains`-agreement genuine set equality. All fifteen
are correct as shipped; the port copies them.

**Loss channels are channels, and the trace statement is about `I + L`, not `I`.**
A site of `LossyPauliWord` carries `I, X, Y, Z, L`, and `I` is the identity on the
*qubit subspace only* (zero on the loss level) while `L` is the loss projector, so
the identity of the full space is `𝟙 = I + L` per site and trace preservation of
the Heisenberg transfer is `Λ*(𝟙) = 𝟙`. `lossChannel_trace_preserving`,
`resetLossChannel_trace_preserving` and `correlatedLossChannel_trace_preserving`
(same file) discharge that for `loss_channel(p)`, `reset_loss_channel` and
`correlated_loss_channel(p₀,p₁,p₂)` — the last for *arbitrary* `(p₀,p₁,p₂)`, with
no normalization hypothesis. That is the oracle the loss port (workload 6) needs
before it lands, and it **vindicates** the old arithmetic the baseline flagged: in
the one-already-lost arms the `1 − p₂` survivor scale and the `p₁` branch weight
belong to different columns of `Λ*(𝟙)` (loss of the second qubit vs. gain from the
both-in-subspace sector), and the `(L,L)` arm's unscaled survivor is *required*,
because loss is irreversible. `lossChannel_nonLossy_scales_by_one_minus_p` records
the complementary convention: on the plain `PauliWord`, where `L` is
unrepresentable and the branch arm is dead code, `loss_channel` scales the trace by
`1 − p` — the documented behaviour, with `ResetLossChannel` as the
trace-preserving variant.

```rust
pub trait Measure {
    fn measure<R: rand::Rng + ?Sized>(&mut self, qubit: usize, rng: &mut R) -> Option<bool>;

    fn measure_many<R: rand::Rng + ?Sized>(
        &mut self,
        targets: &[usize],
        rng: &mut R,
    ) -> Vec<Option<bool>> {
        targets.iter().map(|&q| self.measure(q, rng)).collect()
    }
}
```

The `rng` parameter is the point of [Where the randomness lives](#where-the-randomness-lives):
a measurement is the canonical stochastic operation, and the state it measures
owns no generator. This is the same signature
[`tableau-data-structure.md`](tableau-data-structure.md) specifies for the
concrete tableau.

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

#### Where the randomness lives

**Randomness is injected, never owned.** A generator is not part of quantum
data, so no state type on this trait surface holds one. Every stochastic method
takes `rng: &mut R` (`R: rand::Rng + ?Sized`); the caller decides where the
stream comes from. This is the same kind of separation the rest of this document
performs on `Config`: a per-instance RNG is an *algorithm-driver* choice that the
old design had welded onto the data.

The old crates embed a `SmallRng` in the tableau. That made construction
entropy-consuming, forced `new_with_seed` / `fork(seed)` to exist purely to
manage the field, and left every clone silently deciding whether the copy
replays the stream or diverges from it. It also meant a value advertised as pure
state was neither reproducible from its own contents nor safely `Clone`.

Three consequences make the inversion load-bearing rather than stylistic:

- **Clone is honest.** A classical-mixture branch clones a tableau into two
  branches. With an owned generator both branches inherit the same stream and
  draw identical outcomes — a correctness bug that hides as a statistical one.
  With injection there is no stream to duplicate, so `fork()` is a plain
  `clone()` and the caller derives independent streams explicitly.
- **The states are `Send`-friendly and checkpointable.** Structural identity
  (`PartialEq`, `Debug`, `key_hash`) needs no RNG carve-out, because there is no
  RNG field to exclude. A state is fully determined by its data, so it can be
  serialized, replayed, or shipped to another thread.
- **Seeded reproducibility becomes a caller contract.** Per-shot stream
  derivation is specified where the shots are scheduled, not inherited from
  whatever generator a state happened to carry.

The rule, one layer at a time:

| Layer | Owns a generator? | Surface |
| --- | --- | --- |
| `ppvm-traits-2` (this document) | no | every stochastic method takes `rng: &mut R` |
| the sparse sum and the tableau | no | pure state; `new` is deterministic, `fork()` is `clone()` |
| a stabilizer **mixture** | yes, for **sampling only** | it owns the stream its `sampler()` draws branches from; its channels still take an injected `rng` and branch analytically without drawing |
| frontends (Python, circuit executors) | **yes** | this is the point of the split |

The frontends are where a "global generator for convenience" belongs: a Python
user writing `tab.measure(0)` should not thread a generator, so each binding
wrapper owns one, seeded from its constructor, and passes it into every call that
needs it. Pushing that ownership to the boundary is what lets the core stay pure
without changing the user-facing semantics.

Two members of the channel family take no `rng` because they draw nothing —
`AmplitudeDamping` and `ResetLossChannel`; see the enumeration in
[Behavioral traits](#behavioral-traits) above.

[`tableau-data-structure.md`](tableau-data-structure.md) applies this rule to the
concrete tableau (no RNG field, no RNG in the structural hash).
[`traits-2-implementation-plan.md`](traits-2-implementation-plan.md) records the
migration specifics the rule implies — per-layer ownership during cutover, and
how a given seed is kept reproducing the same outcomes across the old and new
backends.

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
Conjugation by a **Pauli generator** `G ∈ {X,Y,Z}` is the degenerate case that
moves no word at all — a pure `(−1)^{ω(G,P)}` sign — and is the one the `PauliSum`
takes as a hand-written fast path (`Sum::x`/`y`/`z` via `flip_sign_by_key`,
`crates/ppvm-pauli-sum-2/src/clifford.rs`) *instead of* the audited phased-word
`Clifford` impl. Those duplicated signs are pinned directly against the group
product in `lean/PPVM/Pauli/Conjugation.lean`: `conjX`/`conjY`/`conjZ` prove
\(XPX = (-1)^{z}P\), \(YPY = (-1)^{x\oplus z}P\), \(ZPZ = (-1)^{x}P\) as
`G·P·G⁻¹` in \(\mathcal{P}_1\). The surjectivity that upgrades this containment
to the full isomorphism \(\mathcal{C}_n/\mathcal{P}_n \cong \mathrm{Sp}(2n,2)\) is
stated here but not formalized.

The six **derived extension gates** (`CliffordExtensions`) are implemented once,
for every `BlanketClifford` opt-in, as *products of these audited generators*
rather than as fresh hand-written bit-and-sign rules
(`crates/ppvm-traits-2/src/pauli.rs`) — which is only sound if the products
really do reproduce the standard conjugation table, phase included, for a
phase-*carrying* opt-in as well. That composition is discharged in
`lean/PPVM/Pauli/Conjugation.lean`. Because the crate conjugates **backward**
(\(P \mapsto U^\dagger P U\)), its `s` gate is `conjSdag` and its self-adjoint
`h`/`z` gates are `conjH`/`conjPauliZ` (the latter tied to the group product by
`conjPauliZ_eq_conj`); the five single-qubit rows are then the composites
`extSdag` (= `s;z`), `extSqrtX` (= `h;s;h`), `extSqrtXdag`, `extSqrtY` (= `h;z`),
`extSqrtYdag`. Each is a `MonoidHom` **by construction** (`extSdagHom`,
`extSqrtXHom`, … are literally `MonoidHom.comp`s of the audited generator homs,
`*_apply` by `rfl`), and its table is checked against the old crate's:
`extSdag_eq_conjS` (the crate's `s_dag` *is* the forward \(S\)-conjugation),
`extSqrtX_X`/`_Y`/`_Z` (\(X\mapsto X\), \(Y\mapsto -Z\), \(Z\mapsto Y\)) and the
`extSqrtXdag_*`/`extSqrtY_*`/`extSqrtYdag_*` rows, with the dagger pairs proven
mutually inverse (`extSqrtXdag_extSqrtX`, `extSqrtYdag_extSqrtY`,
`extSdag_conjSdag`) and `extSqrtX_sq`/`extSqrtY_sq` recovering
\(\sqrt X^2 = X\), \(\sqrt Y^2 = Y\). Those tables fix each composite on the three
basis Paulis; what a **port** has to reproduce is the closed form, because the old
crate ships these gates as hand-written two-bit `map_add` kernels with no
derivation in the source (`crates/ppvm-pauli-sum/src/sum/clifford.rs:107-221`), and
a sign transposed between a gate and its dagger is invisible on one-hot test
vectors. So each closed form is derived from the generator product too:
`extSqrtX_bits`/`extSqrtX_sign` (\(x \leftarrow x\oplus z\), \(-1\) on \(x\wedge z\)),
`extSqrtXdag_bits`/`_sign` (same bit map, \(-1\) on \(\neg x\wedge z\)),
`extSqrtY_bits`/`_sign` (the swap, \(-1\) on \(\neg x\wedge z\)) and
`extSqrtYdag_bits`/`_sign` (the swap, \(-1\) on \(x\wedge\neg z\)) — the
\(\sqrt Y\) pair sharing a bit map and differing *only* in that predicate is
exactly the transposition the theorems rule out. (`s_dag` needs no entry:
`extSdag_eq_conjS` collapses it to `conjS_bits`/`conjS_sign`.) Two-qubit `CY` is the same story one level
up: `conjCY` is *defined* as the crate's call sequence
\((I\otimes S)\cdot\mathrm{CNOT}\cdot(I\otimes S^\dagger)\) on \(\mathcal{P}_2\)
(`conjSdagT` / `conjCNOT` / `conjST`, with `conjCY_calls` showing the literal
four-primitive sequence `s(t); cnot; s(t); z(t)` collapses to it and
`conjZT_conjSdagT_eq_conjST` the inlined `s_dag`), it is a `MonoidHom` for free
(`conjCYHom`), and `conjCY_bits` + `conjCY_sign` together *are* the old crate's
16-entry table — the bit rule \(z_c \mathbin{\oplus}= x_t\oplus z_t\),
\(x_t \mathbin{\oplus}= x_c\), \(z_t \mathbin{\oplus}= x_c\) and a \(-1\) on
exactly \(X\!\otimes\!X\) and \(Y\!\otimes\!Z\) (`conjCY_Xc`…, `conjCY_XcXt`,
`conjCY_YcZt` name the generator and signed entries). Corollary, and the reason
the \(\pm1\) drain stays total on the extension gates too: every composite delta
is still real — `extSdag_isRealPhase`/`extSqrtX_isRealPhase`/…/
`conjCY_isRealPhase`.

The **loss-guarded** variant of this action — each generator
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
`cz_bits`, which is exactly `czAct`, so `czActL` already models it.) **`CY` is
the strictly harder case of the same claim** and gets its own theorem: the
blanket has no `CY` primitive at all, it emits `s(t); cnot(c,t); s_dag(t)`, and
those guards do *not* agree — `sActL` tests `lost t` alone while `cnotActL` tests
`lost c ∨ lost t`. So with a **lost control and a present target** the atomic gate
is skipped yet two of the three primitives still run, and correctness rests on an
exact cancellation rather than a uniform skip.
`sActL_cnotActL_sActL_eq_cyActL` (`lean/PPVM/Pauli/Symplectic.lean`) proves the
guarded composite equals the old crate's atomic whole-gate skip `cyActL` on every
loss configuration (`sAct_cnotAct_sAct_eq_cyAct` is the unguarded half, `cyAct`
being the old `fn cy` bit rule verbatim; `zActL` records that the crate's `z` is a
bit-level no-op), with `cyActL_preserves_loss` and `cyActL_present_isometry`/
`cyAct_isometry` the loss-invariant and `Sp` halves. The *phase* half of that
cancellation — a phase-carrying opt-in, where the two guarded `S` conjugations
must cancel sign-and-all — is `conjS_conjSdag`/`conjSdag_conjS` in
`Conjugation.lean`. That
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
linearly independent (`frame_linearIndependent`), **span** (`frame_surjective`,
by the `|𝔽₂^{2n}| = |Sp n| = 4ⁿ` counting argument on the injective coordinate
map `frameCombine`), start as one (`isSymplecticFrame_identity`), and stay one
under every Clifford generator
(`isSymplecticFrame_hAct`/`sAct`/`cnotAct`/`czAct` via `IsSymplecticFrame.map`).
The `anticommuting_pivot` search rests on the measurement dichotomy
(`measurement_dichotomy`): the outcome is deterministic exactly when the measured
Pauli commutes with every stabilizer (`measure_deterministic_iff_xfree`).

`compute_decomposition` (Yoder-2012 Lemma 5) is likewise a theorem, not a cited
claim: `frame_coordinate_expansion` proves

```text
v = Σᵢ ω(v, sᵢ)·dᵢ + Σᵢ ω(v, dᵢ)·sᵢ
```

for every Pauli `v`, i.e. the anticommutation bitmasks the routine accumulates
(`stab_anticomm_bits` = the `ω(v, sᵢ)` coordinates, `destab_anticomm_bits` = the
`ω(v, dᵢ)` coordinates) really *are* `v`'s coordinates in the frame basis. The
Rust returns `p_word.phase` without ever asserting that the residual word has
collapsed to the identity; this theorem is what guarantees it has. Everything
downstream inherits it — the branch relabel `idx ^ stab_anticomm_bits`,
`get_deterministic_outcome`, `expectation` and `compute_decomposition_word`.

`canonicalize` is a **no-op** on the tableau backend, and that is now a theorem
rather than a claim. `IsSymplecticFrame.map` covers only the *unitary*
generators — it needs an `ω`-isometry, which the Aaronson–Gottesman measurement
projection is not (it overwrites two rows and is not injective on `GF(2)^{2n}`).
Since that projection is the only non-unitary frame mutation in the crate
(`update_tableau_according_to_outcome`), `Frame.lean` proves it separately:
`isSymplecticFrame_projectFrame` shows that multiplying the pivot stabilizer into
every other row with the measured `x`-bit set, then setting `dₚ := sₚ` and
`sₚ := (−1)^b Z_q`, again satisfies `IsSymplecticFrame` (`projectFrame`,
`rowUpdate_eq_ite` ties the `𝔽₂`-scalar form to the crate's conditional multiply;
the outcome sign is a phase and so invisible on `Sp n`). So the `2n` rows are a
symplectic basis after **every** public operation, unitary or not, and each
subsequent `compute_decomposition` (Yoder-2012 Lemma 5) runs on a valid basis.

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
#[derive(Clone)]
pub struct CoefficientThreshold {
    pub threshold: f64,
}

// The default threshold is a USER-FACING value, carried over unchanged from
// `ppvm-pauli-sum::strategy::CoefficientThreshold` (1e-12). A derived `Default`
// would give 0.0 and silently keep terms the current crate drops.
impl Default for CoefficientThreshold {
    fn default() -> Self {
        Self { threshold: 1e-12 }
    }
}

impl<W: Word + Indexable, C: Coefficient> Policy<W, C> for CoefficientThreshold {
    // Also carried over: the capacity hint is the current strategy's `n * 10`,
    // not 0. Both maps of the double-buffer are sized from it at construction,
    // and the caller can override it (`Sum::with_capacity`, the port of the
    // current builder's `.capacity(..)`, which every real workload passes).
    fn capacity(&self, n_sites: usize) -> usize {
        n_sites * 10
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
for the `PauliSum` path (`l1_bound`, under `|⟨P⟩| ≤ 1`; and `l1_bound_abv` /
`l1_bound_norm` / `l1_bound_complex` for a general coefficient ring whose
`magnitude` is an absolute value, which is what makes this keep-rule meaningful
for `PauliSum<Complex<f64>>` and not only for real coefficients) and, for the
tableau path,
the *unconditional* Cauchy–Schwarz `ℓ²` bound `l2_bound`, sharpened to
`error² ≤ (Σ_{dropped} c_P²)·|D|` under `|⟨P⟩| ≤ 1` in `l2_bound_normalized`.

The tableau spells its own cutoff three ways — `|c|² > t²` in the gates,
`|c| > |t|` in `rotate_2`, `|c|² > t²·‖v‖²` in the case-a measurement — and
`cutoff_abs_iff_sq` (same file) settles that only **two** of them are different
rules: for `t ≥ 0` the first two are the same predicate, so `rotate_2`'s spelling
differs only in float rounding and in paying a `hypot` per element. The real
split is absolute versus **relative**, and it is the relative form that the error
bounds above are stated for. All three are reproduced verbatim from old under the
behaviour-preservation directive; unifying them is a sign-off item.

Two properties of `retain` itself — not of any one policy — are what license the
shipped `Policy` implementations, and both are machine-checked in
`lean/PPVM/Algebra/GradedMap.lean` (§"`Retain`"):

* `CombinedPolicy::truncate` runs its two members as **two sequential `retain`
  passes** (a structural difference from PauliPropagation.jl's single combined
  walk, carried over from `CombinedStrategy`). `retain_seq_eq_retain_and` proves
  the sequential form computes the **conjunction** of the two keep-rules, and
  `retain_comm` that the surviving key set does not depend on the pass order — so
  the surviving set is a property of the policy pair, and a future backend
  (parallel, columnar) may fuse or reorder the passes without re-litigating.
* `MaxPauliWeight::truncate` (and `MaxLossWeight::truncate`) **early-returns** on
  the `usize::MAX` disable sentinel rather than running an always-true walk —
  the headline `CombinedPolicy(CoefficientThreshold(1e-6), MaxPauliWeight(MAX))`
  configuration is why. `retain_of_all_true` / `retain_le_top_eq_self` are its
  soundness: retain-all is the identity *pointwise*, so skipping the pass is
  observationally exact, zero-coefficient terms included.

`Sum::truncate`'s **preserved-key post-filter** (the port of the builder's
`preserve_strings`) is adjudicated in `lean/PPVM/Algebra/Truncation.lean`
(§"The preserved-key post-filter is a widened keep-rule"). Its three steps —
snapshot the preserved keys' *pre-truncate* coefficients, run the policy
verbatim, re-insert what the policy dropped through the accumulating
`AddTerm::add_term` under a membership guard — collapse to a single pass with the
widened keep-rule `keep k c ∨ k ∈ P` (`truncate_preserve_eq_widened_retain`).
That equivalence is what pins the two otherwise invisible guards (without the
membership test a survivor's coefficient would be **doubled**; with the snapshot
taken after the policy it would restore a post-truncate residue), and it gives
the two facts a caller relies on: a preserved key keeps exactly its pre-truncate
coefficient (`truncatePreserve_apply_of_mem` — old's `Σᵢ Zᵢ` conservation test),
the empty keep-set is the policy verbatim so the hot-path short-circuit is exact
(`truncatePreserve_empty`), and the dropped set is `D \ P`, which is the set the
`ℓ¹` bound above then applies to.

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

    /// Borrowing scan — the read side for a reader that REJECTS most terms
    /// (`trace` against a pattern, a predicate count, a weight scan). `iter`
    /// hands out owned pairs so that a columnar backend can synthesize them,
    /// which forces a `Coeff::clone()` on every key before the reader looks at
    /// it: free for `f64`, ruinous for a coefficient owning a heap table (the
    /// symbolic `Term` measured 7×–33× slower than the old crate's borrowing
    /// fold on a 65k-key support with 255 matches, and the gap grew with
    /// coefficient size). A callback rather than `impl Iterator<Item = (&K,&C)>`
    /// so no backend needs a lending iterator; a SoA backend passes
    /// `(&keys[i], &coeffs[i])`, so this does not cost the columnar option.
    /// Defaulted through `iter`, so it is optional for a backend to specialize.
    fn for_each_ref(&self, f: impl FnMut(&Self::Key, &Self::Coeff)) { … }
}

/// L1 — the module core: form linear combinations, then canonicalize.
pub trait Accumulate: Support {
    /// Build side of the hash join: merge a produced batch, accumulating onto an
    /// existing key or inserting a new one. Columnar in; the scalar
    /// `accumulate(k, c)` is provided sugar over a batch of one.
    ///
    /// The batch is algebraically a **multiset** of terms: folding it in is a
    /// homomorphism from the free commutative monoid, so the result is
    /// independent of term order (`accumulateTerms_perm`) and of how the batch
    /// is split across partitions/threads (`accumulateTerms_add`), and the
    /// scalar sugar is definitionally a batch of one
    /// (`accumulateTerms_singleton`) — machine-checked in
    /// `lean/PPVM/Algebra/GradedMap.lean`. That is the licence for a backend to
    /// `gather` a batch into per-partition sub-batches and run them
    /// concurrently.
    ///
    /// A collision must **sum**, never overwrite: old's `AddAssign<PauliSum>`
    /// routes through `HashMap::extend`, which replaces the shared key's
    /// coefficient, and that is adjudicated wrong by the same file
    /// (`accumulate_apply`, witness `accumulate_ne_overwrite`). This engine
    /// diverges from old there by design.
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
/// `lean/PPVM/Algebra/Noise.lean`. The `= Tr(A B)/2ⁿ` half of the label is no
/// longer an assertion: `overlap_eq_trace_div` in `lean/PPVM/Pauli/Matrix.lean`
/// proves `Tr(Â B̂) = 2ⁿ · ⟨A, B⟩` for the *genuine* `2ⁿ×2ⁿ` operators
/// `Â = ∑_p a_p g(p)` over the exact ring `ℤ[i]` (via `trace_tensorPauli_mul`,
/// the real-matrix Pauli orthonormality `Tr(g(p) g(q)) = 2ⁿ δ_pq` that
/// `overlap_single_single` was only the model form of). One level up, inside
/// `C[K]` itself, the same fact reads: the pairing **is** the identity
/// coefficient of the L4 twisted product, `⟨A, B⟩ = (A · B)_I`
/// (`twistedConv_apply_id` in `lean/PPVM/Algebra/Twisted.lean`) — which is what
/// ties L3 to `Multiply` instead of leaving them unrelated layers. Finally, the
/// semantic link to the Clifford path is
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
///
/// Law (the obligation *every* impl owes, not just `PauliWord`): writing
/// `key_mul(u, v) = (u·v, i^{β(u,v)})`, the key product `·` must be associative
/// and `β : K × K → ℤ/4` must be a 2-cocycle,
/// `β(u,v) + β(u·v, w) == β(v,w) + β(u, v·w)`. Under exactly those two
/// hypotheses the twisted product on `C × K` is associative for any commutative
/// `C` with `i⁴ = 1` — machine-checked key-agnostically in
/// `lean/PPVM/Algebra/Twisted.lean` (`gtmul_assoc`, over an abstract `kmul` and
/// `IsCocycle`), with `PauliWord` recovered as the instance
/// (`phaseExp_isCocycle`, `tmul_assoc_of_gtmul`). A future ordered fermionic-word
/// key must discharge the same two hypotheses; it does not inherit associativity
/// from the Pauli proof.
///
/// Third, *independent* obligation — **right-cancellativity**: `p ↦ key_mul(p,q)`
/// must be injective for every fixed `q`. `Sum::mul_word_assign` re-keys the
/// whole support through that map and merges with the plain-`insert`
/// `RekeyBijective` path, where a collision silently **drops** a term instead of
/// summing it. This does *not* follow from associativity + the cocycle law —
/// `lean/PPVM/Algebra/Twisted.lean` (`isRightCancellative_independent`) exhibits
/// a constant key product satisfying both while collapsing two keys onto one.
/// `PauliWord` discharges it (`mulWord_right_injective`,
/// `mulWord_isRightCancellative`); a key that cannot must not take that path.
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
    ///
    /// A ring that carries `i^k` as *data* satisfies this only denotationally,
    /// and the exemption is machine-checked rather than asserted: for
    /// `ppvm-sym-2`'s `Term` the implemented ring is the ℤ/4-graded
    /// `PhasedSymRing`, on which `iSym_sq_ne_neg_one` shows `i·i` and `−one()`
    /// are different keys, `evalC_iSym_sq_eq_neg_one` shows they have the same
    /// complex value, and `evalC_not_injective` shows why both can hold at once
    /// (`lean/PPVM/Instantiations/Symbolic.lean`). The law is discharged
    /// literally by the exact `GaussianInt` witness instead.
    fn imaginary_unit() -> Self;

    /// Multiply by `i` — semantically `self * imaginary_unit()`, but an OVERRIDE
    /// POINT: on `Complex<f64>` the ring product spells the rotation as
    /// `(re·0 − im·1, re·1 + im·0)`, which is `NaN`-contaminating and loses the
    /// sign of zero, where old's `ComplexCoefficient::mul_phase` swapped the
    /// components by hand. The impl restores the swap.
    fn mul_i(&self) -> Self { self.clone() * Self::imaginary_unit() }

    /// Multiply by `i^k` (mod 4) — what `Phase::apply` delegates to, and the
    /// second override point. The default is the four-arm
    /// `{clone, mul_i, neg, neg∘mul_i}` fold, right for any ring whose values
    /// are numbers. A ring that carries `i^k` AS DATA (the symbolic `Term`, whose
    /// monomials hold a ℤ/4 phase byte) must override it: old's `mul_phase`
    /// promoted `Const(c)` to `One(i⁰, c)` unconditionally, including at `k = 0`,
    /// and `Term`'s `PartialEq`/`Display` are representational, so folding
    /// through `clone()`/`neg()` instead would be a user-visible divergence.
    fn mul_i_pow(&self, k: u8) -> Self { … }
}

/// L4 — the ring product. The only layer that needs the *key* to carry a
/// product; it stays optional and is not implemented for a key type that has
/// none. The Pauli product injects powers of `i`, so the coefficient must absorb
/// phase — bounded on `ImaginaryUnit`, the minimal requirement (a primitive
/// fourth root of unity), **not** the stronger `ComplexCoefficient`.
///
/// L4 and L3 are the same structure read two ways: the twisted convolution on
/// `C[PauliWord]` (`twistedConv`, the outer product `multiply_into` computes)
/// has, as its **identity-key coefficient**, exactly `Pair::overlap` —
/// `(A · B)_I = ⟨A, B⟩`, machine-checked as `twistedConv_apply_id` in
/// `lean/PPVM/Algebra/Twisted.lean` (the outer product collapses to the diagonal
/// by `mulWord_eq_id_iff`, where the `i^k` twist vanishes by `phaseExpN_self`).
/// That is the correctness spec a container `Multiply` impl is checked against,
/// and it is why `overlap` deserves the name Hilbert–Schmidt pairing.
///
/// `multiply_into` is **bilinear**: `(A + B)·D == A·D + B·D` and
/// `A·(B + D) == A·B + A·D`, machine-checked as `twistedConv_add_left` /
/// `twistedConv_add_right` in `lean/PPVM/Algebra/Twisted.lean`. So each monomial
/// product must be accumulated into a *fresh* accumulator; folding rhs terms
/// back into an operand computes the product **chain** `A·b₀P₀·b₁P₁` instead of
/// the sum `A·b₀P₀ + A·b₁P₁` (which is precisely the latent bug in old's
/// `MulAssign<PauliSum>` — Lean adjudicates old wrong here, see
/// `crates/ppvm-pauli-sum-2/src/multiply.rs`). Biadditivity is also what lifts
/// the *monomial* `tmul_assoc` to the whole-map law `twistedConv_assoc`
/// (`(A·B)·D == A·(B·D)`, same file): monomial associativity alone does **not**
/// imply the convolution's — bilinearity reduces the general case to monomials
/// and the 2-cocycle settles them.
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

The same licence — and only that licence — permits the `RekeyBijective` plain
`insert` (no accumulation probe, `debug_assert!` only) on any *other* rewrite,
because a non-injective map would silently **drop** a term there rather than sum
it. L4's right-multiplication by a single Pauli word, `p ↦ p·q`
(`Sum::mul_word_assign`, old's `MulAssign<PauliWord>`), qualifies: injectivity is
machine-checked as `mulWord_right_injective` in `lean/PPVM/Pauli/Word.lean`, the
word-product counterpart of the Clifford `Symplectic.*_bijective` /
`Conjugation.conj*_injective`. Any future re-key wanting the fast path owes the
same one-line theorem.

Injectivity alone is not the whole obligation, though: `mul_word_assign` reaches
the *product* by a different code path from `multiply_into`, so it also owes that
the re-key **is** the L4 product against a one-term map.
`lean/PPVM/Algebra/Twisted.lean` closes that: `twistedConv_single_right` shows
`A · (b·q)` collapses to the pushforward of `A` along `p ↦ p·q` (coefficients
scaled by `b` and twisted by `i^{phaseExpN p q}`), and
`twistedConv_single_right_apply` shows the coefficient landing on `p·q` is exactly
the one source term's contribution — i.e. the re-key needs **no** aggregation,
which is precisely what makes `insert` (rather than `entry().and_modify()`)
sound. Without injectivity those two differ and the plain `insert` drops a term
instead of summing it (`isRightCancellative_independent`).

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

Neither of those says what the **case-a measurement** computes, and that is the
component's most load-bearing arithmetic. `lean/PPVM/Tableau/Projection.lean`
covers it. Modeling the frame-conjugated `Z_q` as `M|k⟩ = i^{φ k}|k ⊕ s⟩` (an
XOR shift with a `ℤ/4` phase, `s = stab_anticomm_bits`), it proves:
`rustTerm_eq` — the overlap merge's four-way `ℤ/4` dispatch
(`0 ⇒ +re_w, 1 ⇒ +im_w, 2 ⇒ −re_w, 3 ⇒ −im_w`) is exactly
`Re(conj(iᵠ·a)·b)`, *not* `Re(iᵠ·conj a·b)` — the odd branches carry the
conjugated phase, which is precisely where a sign slip would silently change a
probability; `shiftOp_involutive` / `shiftOp_selfAdjoint` — `M² = I` and
`M† = M`; `overlap_eq_inner` — the crate's `z_overlap_re` is `Re⟨c, M c⟩`;
`proj_add` / `proj_idem` — `P₀ + P₁ = I` and `P_b² = P_b` for
`P_b = (I + (−1)^b M)/2`, so the case-a step really is a projective measurement;
`probOne_eq` — `prob_1 = 0.5 − 0.5·z_overlap_re` **is** the Born probability
`⟨c, P₁ c⟩` for a normalized amplitude vector; and
`projectRaw_eq_two_proj` — the keep-`A`/transform-`B`/merge map is `2·P_b` on the
surviving half, the factor `2` being what the subsequent unconditional
`normalize()` removes. One link is deliberately left open and is documented in
that file's scope note: the projection's phase `ψ` (`alpha + 2·⟨idx, destab⟩`)
omits the odd-phase-destabilizer term the overlap folds in through
`compute_phase_with_mask_static`, because the two are read in different frames
(pre- and post-projection); relating them needs a Hilbert-space model of the
frame that this development does not have. Old and `-2` agree there verbatim, so
it is a specification gap, not a port divergence.

Every one of those theorems is stated under the hypothesis `SelfInverse s φ`, and
`Projection.lean` keeps `φ` *abstract* — so on its own it says nothing about the
`φ` the crate computes. `lean/PPVM/Tableau/BranchPhase.lean` closes that link. It
models the amplitude basis as `|j⟩ = ∏_l d_l^{j_l}|ψ₀⟩` and *derives* the crate's
formula: `frameOp_eq_shiftOp` proves that `i^{pd}·D_L·S_G` acting on amplitudes is
exactly `phase_decomp + 2·⟨destab_anticomm, j⟩ + 2·popcount(j ∧ stab_anticomm ∧
odd_phase_mask)` — i.e. `compute_phase_with_mask_static` (`data.rs`) is the
composite of the generator actions, not a convention. On top of that:
`frameOp_sq` / `frameOp_involutive_iff` give `M² = i^{2·pd + 2(⟨G,L⟩+⟨L,m⟩)}·I`,
so `M² = I` **iff** the `ℤ/2` *frame identity*
`phase_decomp + ⟨destab, stab⟩ + popcount(stab ∧ mask) ≡ 0`;
`selfInverse_branchPhase_iff` proves that identity is *equivalent* to
`SelfInverse`, and `selfInverse_branchPhase` discharges it — with
`shiftOp_involutive_crate`, `proj_idem_crate` and `probOne_eq_crate` restating
`M² = I`, `P_b² = P_b` and the Born rule for the concrete phase. Specialized to
case b (`stab_anticomm = 0`), `frameInvolution_zero_iff` says the identity *is*
`phase_decomp ∈ {0, 2}` — the crate's `debug_assert!` "Measurement result cannot
be imaginary!" is now a theorem. `destabAction_sq` (`d_L² = (−1)^{popcount(L∧m)}`)
plus `add_phase_eight_sub` (`8 − 2·ph = 2·ph` in `ℤ/4`) are the algebraic content
of `compute_decomposition`'s "multiply the generator in and divide its phase
squared out" step, and `stab_destab_commute_sign` shows the *opposite* visit order
would shift `φ` by `2·⟨G, L⟩` — so the two-loop order ("all stabilizers first") is
a real convention, and the crate's use of the **original** index in the
`⟨destab, ·⟩` term is the one that matches it.
`crates/ppvm-conformance-2/tests/tableau_lean.rs` checks the frame identity on
every decomposition of a random Clifford+`T` sweep, so the discharge is pinned to
the code, not just the model. (That sweep also confirms `odd_phase_mask` is always
zero on a valid frame — rows are Hermitian, hence even-phased — so the mask term
is a defensive no-op kept for behaviour parity; the theorems hold for arbitrary
masks regardless.)

`rotate_2` is covered by the same closed forms. It never builds a two-site
decomposition: it applies `compute_coefficients_after_pauli_apply` at `b` and then
at `a`, two independent single-site relabels. `shiftOp_comp` proves the composite
is again an operator of the same shape — shift `L_a ⊕ L_b`, weight `w_a + w_b`,
phase `pd_a + pd_b + 2⟨w_a, L_b⟩` — so the sequential relabel **is** the
frame-conjugated `P_a ⊗ P_b`, phase included, and `rotate_2` therefore really is
`cos(θ/2)·I − i·sin(θ/2)·(P_a ⊗ P_b)`. `dot_crateWeight_order` /
`shiftOp_comp_order_iff` isolate the entire order dependence in
`⟨G_a, L_b⟩ + ⟨G_b, L_a⟩` (the mask term is symmetric), which
`omega_eq_frame_coords` (`Frame.lean`, new) identifies as `ω(P_a, P_b)` read in
frame coordinates; since `rotate_2`'s two Paulis sit on distinct qubits they
commute (`omega_disjoint_support`), so `rot2_order_irrelevant` /
`rot2_order_irrelevant_of_commuting` prove the `b`-before-`a` order does **not**
affect the accumulated `ℤ/4` phase. The port keeps that order anyway: what it
pins is the float summation order, which is a separate concern.

Multi-qubit `expectation(&word)` goes through `compute_decomposition_word`, which
conjugates each site separately and folds the single-site results in the
canonical form `iᵠ Xˣ Zᶻ`, correcting with `(−1)^{popcount(z_running ∧ x_new)}`.
`lean/PPVM/Pauli/Word.lean` proves that fold is the genuine operator product:
`phaseExpN_eq_canon` reconciles the canonical `Xˣ Zᶻ` normalization with the
`g(x,z) = i^{x·z} Xˣ Zᶻ` one that `phaseExpN` (hence
`PauliMatrix.tensorPauli_mul`, hence real `2ⁿ×2ⁿ` matrices) is stated for, so the
cross term is operator multiplication rather than a convention;
`Canon.toG_mul` states that step-wise (`toG (a·b) = toG a + toG b + phaseExpN`);
and `crossPhase_cocycle` / `Canon.mul_assoc` / `Canon.foldl_eq_prod` show the
left-fold with a *running* Z-mask equals the ordered product of the per-site
conjugates. This is the only piece of the expectation path with no single-qubit
oracle — the cross term vanishes identically at weight 1.

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

The sketch above is missing one step, and the omission is load-bearing: a
re-keying producer emits `(φ(k), c')` while the map still holds `(k, c)`, so
`accumulate_batch` **onto the live support leaves both** — a silent double-count
of the entire support on every Clifford gate. The support must therefore be
**reset between producing and accumulating** (`ApplyProducer::apply_producer` in
`ppvm-pauli-sum-2/src/store.rs` does this). Machine-checked in
`lean/PPVM/Algebra/GradedMap.lean`: `pushforward_eq_reset_accumulate` shows the
reset-then-accumulate composite equals the pushforward `mapDomain φ ∘ mapRange g`
(with `pushforward_apply` giving the injective, per-key form for a Clifford's
symplectic bijection), and `merge_without_reset_ne_pushforward` exhibits the
double-count that merging without the reset produces. This is the `apply`-path
analogue of `eagerWalk_ne_twoPass`, and it matters most where it would be hardest
to spot: the Phase-4 tableau-keyed reuse of `apply`.

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
| `CliffordExtensions` (`s_dag`, `sqrt_x`, `sqrt_x_dag`, `sqrt_y`, `sqrt_y_dag`, `cy`, `zcy`) | `CliffordExtensions`, unchanged shape; blanket-implemented over the same `BlanketClifford` opt-ins | Same user-facing gate set and the same required/defaulted split. Only the blanket's *derivation* changes: each gate is a product of audited `Clifford` generators (`S† = S·Z`, `√X ≃ H·S·H`, `√Y ≃ H·Z`, `CY = (I⊗S)·CNOT·(I⊗S†)`) rather than a fresh hand-written bit rule, so the blanket is correct for phase-carrying opt-ins (`Tableau`) without six new unproved `PhaseTrack` deltas. Diffed against the old tables in `ppvm-conformance-2`; `PhasedPauliWord` keeps a read-once fused override, as for `Clifford`. |
| `CliffordBatch`, `CliffordExtensionsBatch` | retained verbatim | Loop defaults over the single-qubit methods, overridable by the tableau's fused sweeps; the "empty `impl` opts in to the defaults" convention is preserved. |
| `Reset` (`reset`, `reset_z`, `reset_x`, `reset_y`, `*_many`) | retained; each method gains `rng: &mut R` | The basis variants are *behaviour*, not sugar (`reset_x` = `reset` then `h`, `reset_y` = `reset` then `h` then `s`), so the default bodies and the `Clifford + CliffordExtensions` supertrait bound are reproduced call-for-call — they just thread the injected generator through. Reset is stochastic (it measures, then corrects), hence the parameter; see [Where the randomness lives](#where-the-randomness-lives). |
| `RotationOne<T: Config>` (`rotate_1`, `rx`/`ry`/`rz`, `rx_many`/`ry_many`/`rz_many`) | `RotationOne<C, A = C>`, same required/defaulted split | The whole surface crosses over: `rotate_1(axis, qubit, theta)` stays the *required* axis-generic entry point (it is what `rotate_2` and the tableau backends compose with), `rx`/`ry`/`rz` stay one-line defaults over it (a backend with per-axis fast paths overrides them, as `PauliSum` does), and the three `*_many` batch loops stay (the Python bindings call them). Two edits follow from the shape change: the `Pauli::L` axis panic is unrepresentable because `L` is not a `Pauli` any more (loss is a `LossySite`), and `theta` is the angle domain `A` rather than `impl Into<T::Coeff>` — with `A` a free trait parameter an `Into` conversion would be uninferable at the call site, so the one instantiation callers used (`sum.rx(0, 0.1)` on a complex sum) is preserved by `impl Angle<Complex<f64>> for f64` instead. The batch defaults clone the angle, so they carry an explicit `where A: Clone` (the old crate got it free from `Coefficient: Clone`). |
| `RotationTwo` (`rotate_2` + the `rxx`…`rzz` family and their `*_many`) | `RotationTwo<C, A = C>`, unchanged shape | User-facing gate surface implemented today by `ppvm-pauli-sum` and `ppvm-tableau-sum`; the `[x, z]` axis encoding, the nine named gates and the batch loops are ported verbatim. Only the angle domain and the `where A: Clone` on the batch defaults change. |
| `TGate` (`t`, `t_dag`, `*_many`), `Projection` (`p0`, `p1`) | retained; both unparameterized | Neither takes a numeric argument — old `TGate<T: Config>` never used its parameter — so the rule that leaves `Clifford` and `ResetLossChannel` unparameterized applies. Default bodies verbatim, with one exception: `Projection`'s halving. `lean/PPVM/Instantiations/Projector.lean` (`projLin_add`/`projLin_smul`/`projLin_idem` vs `oldStep_not_additive`/`oldProj_not_idem`, `oldStep_eq_half_iff`) adjudicates old's `let half = v.half(); *v *= half` — which computes `c ↦ c²/2` — as a genuine defect; the Lean-correct halving is the ring constant `½`. |
| `U3Gate` (`u3`), `RotXY` (`r`), `CRx` (`crx`) | `U3Gate<C, A = C>`, `RotXY<C, A = C>`, `CRx<C, A = C>` | Angle-carrying single/two-qubit gates; retained with the angle-domain edit only. `RotXY::r` keeps its documented `RZ(φ)·RX(θ)·RZ(−φ)` decomposition as *behaviour* of the implementing backend — the Heisenberg (backward) sub-rotation order, machine-checked in `lean/PPVM/Instantiations/Rotation.lean` (`rotXY_heisenberg_order`, `rotXY_halfPi_eq_ry`). |
| `Trace<'a, RHS>` | retained (in `graded.rs`, next to `Pair`) | **Not** subsumed by `Pair::overlap`: `overlap` pairs a map with another map of the same type, while `Trace` is the heterogeneous `tr(self·value)` the old crate used against a `PauliPattern`. Kept with its free right-hand type and `Output`; implementers land with the pattern port. |
| `PauliWordTrait::anticommutes_at` | provided method on `PauliBits` | Derivable from the two bit reads (`ω(P,Q) = x_P·z_Q ⊕ z_P·x_Q`), and consumed by the tableau's measurement pivot search, so it stays — as a default body, adding no required method. |
| `PauliWordTrait::get_multiple`/`get_slice`/`set_multiple`/`set_new`/`set_new_2` | removed from the trait; inherent on the concrete word | Every one of them *constructs* a word (`Self::new(Q)`, clone-and-edit), and construction is deliberately concrete in this design — `Word` is read-only inspection and has no constructor, `PauliBits` is bit mutation. They remain available as inherent methods on `PauliWord`/`LossyPauliWord`, where the width and backing storage are known; no generic consumer needs them. |
| `PauliError<T: Config>` stim aliases (`x_error`/`y_error`/`z_error`, `*_many`) | retained on `PauliError<C: Coefficient>`; each method gains `rng: &mut R` | One-hot defaults over `pauli_error`, unchanged in body. Two edits to the signatures: the `Config` parameter becomes the coefficient type itself, and the generator is injected ([Where the randomness lives](#where-the-randomness-lives)). |
| noise channels: `PauliErrorAll`, `TwoQubitPauliError`, `Depolarizing`, `Depolarizing2`, `AmplitudeDamping`, `LossChannel`, `CorrelatedLossChannel`, `ResetLossChannel`, `AsymmetricLossChannel` | retained, each `<T: Config>` → `<C: Coefficient>`; the seven **drawing** channels gain `rng: &mut R` | The channel family is user-facing behaviour, so the defaults and probability orderings move across unchanged. Two signature edits, applied by different rules. Coefficient: `ResetLossChannel` consumes none, so it drops the parameter entirely — the same rule that leaves `Clifford` unparameterized. Generator: injected into every channel that *draws* ([Where the randomness lives](#where-the-randomness-lives)), which excludes `AmplitudeDamping` (a deterministic Kraus rescaling) and `ResetLossChannel` (an unconditional bit clear). The two carve-outs are not the same set, so neither parameter can be inferred from the other. |
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
Trotter optimization relies on. The branch's per-axis `±1` sign `ε` is *not* a
free convention: `iGP = i^{1 + phaseExp(G,P)} · g(G⊕P)`, so `ε` is the real phase
`i^{1+phaseExp}` of the anticommuting single-qubit product `G·P` (the leading `i`
of `iGP` cancels the product's `i`, `branchExp_isRealPhase`). `rx_eps_from_product`
/`ry_eps_from_product`/`rz_eps_from_product` check that this phase-derived `±1`
equals the hand-ported table in `crates/ppvm-pauli-sum-2/src/producer.rs:141-143`
(`RotationProducer::produce`, `rx: ε=−1 iff x`; `ry: ε=−1 iff z`; `rz: ε=+1 iff z`),
grounding the one propagation sign the abstract `rot` model does not derive.

The **two-qubit** branch (`RotationTwo`: the headline workload's native `rzz`, plus
`rxx`/`ryy`) is not a trivial special case of that, and it is the one gate family
whose fast path the old crate validates only by a Rust diff against its own
generic `comm_2`/`rotate_2` (`rot2.rs::rxx_matches_generic`) — agreement between
two implementations, not an oracle. The branch key toggles bits at **two** sites
at once and \(\varepsilon\) comes from a two-factor product, so the sign
convention has two independent places to be wrong. The same file now supplies the
oracle. `anti2` is the two-site branch key
\(w_2(P,G) = \omega(P_a,G_a)\oplus\omega(P_b,G_b)\) (`antiBit_eq_omega` ties the
kernels' two-bit expression to \(\omega\)); `branchExp2` is the two-site
\(i^{1+\text{phaseExp}_a+\text{phaseExp}_b}\), real whenever the branch fires
(`branchExp2_isRealPhase`), and `anticommute_new_key2` again gives exactly one
fresh term. The collapse onto a *single* site's \(\varepsilon\) column is
`phaseExp_of_commute` — a commuting single-site pair carries **no** phase at all —
after which each native kernel is pinned expression-for-expression:
`rzz_anti`/`rzz_branch_key`/`rzz_eps_from_product`, and the `rxx_*` / `ryy_*` rows
(`sum/rot2.rs:62-184`: the commutation test, the bits toggled, and
`eps = if z_anti {1} else {-1}` etc.). `accumulate_rotBatch_two` records that the
fused two-pass shape lifts to the two-site branch unchanged — `accumulate_rotBatch`
is stated over an abstract key type and predicate, so the two-qubit kernel is an
instance of it, not a new obligation.

The **generic** `rotate_2` — the shipped public method that accepts an arbitrary
`[x, z]` axis pair — does not use those three kernels at all: it runs the
branch-free `comm_2` (`crates/ppvm-pauli-sum-2/src/rotation.rs:109-140`) over a
hand-rolled 16-entry orientation mask `SIGN_NEG = 0x2840`, and then applies the
*opposite* sign convention to the fast paths (`sin.mul_sign(-eps)` against their
`sin.mul_sign(eps)`). The only checks on any of that — old's `*_matches_generic`
tests and their ports — compare the two paths at exactly the three **diagonal**
axes; every off-diagonal pair (`XZ`, `YX`, `ZY`, …) was untested on both sides and
unproven, and a `+eps`/`−eps` asymmetry is precisely what can be correct by
coincidence on the diagonal. `lean/PPVM/Instantiations/Rotation.lean` (§"The
generic `rotate_2` kernel") now transcribes `comm_2` line for line (`signNegMask`,
`signNegIdx`, `comm2Coeff`, `comm2Key`) and closes it over all `2⁸` (axis, key)
bit patterns: `comm2Coeff_eq_zero_iff` — the `if eps == 0 { return None }`
early-out is exactly the `anti2` branch predicate the native kernels test;
`comm2_generic_sign_eq_branchExp2` — the sign the generic path actually applies
(`−ε`) **is** the real branch prefactor `i^{1+phaseExp_a+phaseExp_b}` of
`−i·[G_a ⊗ G_b, P]/2`, so both the mask and the sign flip are right off the
diagonal too; and `comm2_key_eq_mulBits2` — its `present`-masked output bits are
the two-site product key, so `anticommute_new_key2` applies to the generic path
verbatim.

`RotXY::r` is the one rotation-family contract whose **order** (not its per-axis
arithmetic) is user-visible: the crate emits `rz(φ); rx(θ); rz(−φ)` — Heisenberg
(backward) — so the composite is `M_z(−φ) ∘ M_x(θ) ∘ M_z(φ)`, and a
forward-ordered implementation yields `ry(q, −θ)` at `φ = π/2` while passing every
other rotation test. Same file, §"`RotXY::r`": the per-axis `ε` columns are
assembled into the `3 × 3` action on the coefficient triple `(c_X, c_Y, c_Z)`
(`mz`/`mx`/`my`, whose off-diagonal entries are checked against `mulBits` and
`branchExp` by `mz_from_kernel`/`mx_from_kernel`/`my_from_kernel`, i.e. read off
the kernel rather than modeled), and `rotXY_heisenberg_order` proves the composite
is Rodrigues rotation about the in-plane axis `cos φ·X + sin φ·Y` by `θ`
(`rotAxis`). The two behavioural corollaries drop out: `rotXY_zero_eq_rx`
(`r(q,0,θ) = rx(q,θ)`) and `rotXY_halfPi_eq_ry` (`r(q,π/2,θ) = ry(q,θ)`, **not**
`ry(q,−θ)`) — contract 10's order detector, previously pinned only by two ported
example values.

Those are all *per-term* facts. The fused `RotateInPlace` fast path
(`crates/ppvm-pauli-sum-2/src/store.rs`) that every rotation actually takes is a
**two-pass** whole-map walk — pass 1 scales every diagonal in place and buffers
the branch terms in `scratch`, pass 2 merges only the buffered branches — and its
risk is *cross-term*: a branch key produced from `k` can collide with a different
key `k'` that pass 1 has not scaled yet (`rx` on a support holding both `Z` and
`Y` at one site swaps them). `anticommute_new_key` is a near miss here; it only
rules out a branch colliding with its own source. The whole-map licence is
`accumulate_rotBatch` in `lean/PPVM/Instantiations/Rotation.lean`: the two-pass
map equals folding the whole `≤ 2N`-term produced batch into an empty map (the
generic producer → `accumulate_batch` semantics), for **every** walk order, since
the batch is a `Multiset` and `GradedMap.accumulateTerms_perm`/`_add` apply. The
two-pass structure is load-bearing rather than an optimization:
`eagerWalk_ne_twoPass` exhibits a support on which merging each branch eagerly
*inside* the walk — the tidier single-pass refactor, and what a backend
interleaving the passes computes — gives a different map.

`Projection` (`p0`/`p1`) was the remaining shipped gate on `Sum` with **no**
oracle at all: `lean/PPVM/Tableau/Projection.lean` is about the generalized
tableau's amplitude vector, not the Heisenberg action of a computational-basis
projector on `C[K]`, and no old test or bench exercises `p0`/`p1`.
`lean/PPVM/Instantiations/Projector.lean` supplies it, and it adjudicates the old
kernel **wrong**. The intended map — halve by the ring constant `½` and add the
`Z`-toggled partner with sign `ε` — is linear (`projLin_add`, `projLin_smul`) and
idempotent (`projLin_idem`), i.e. a genuine projector on `C[K]`. Old instead reads
`let half = v.half()` (a *value*, `c/2`) and then does `*v *= half`, so both the
survivor and the branch come out at `c²/2`: `oldStep_not_additive`,
`oldStep_not_homogeneous` and `oldProj_not_idem` show that map is neither linear
nor idempotent (on `2·I` the "projector" *grows* the state), and
`oldStep_eq_half_iff` pins the blind spot exactly — `c²/2 = c/2 ↔ c ∈ {0,1}`, and
unit-coefficient stabilizer sums were old's only usage. The Lean-correct value is
`c/2`, so the implementation builds the ring's `½` once outside the walk. A
*second*, independent defect is corrected by the same Lean-governed exception:
over honest `ℤ[i]`
matrices with `2Π = I + Z`, `twoProj_conj_I`/`twoProj_conj_Z` agree with `projLin`
on the `I`/`Z` block but `twoProj_conj_X`/`twoProj_conj_Y` give `Π X Π = Π Y Π = 0`,
where the crate's `_ => None` leaves `X`/`Y` untouched — `projLin_p0_add_p1` shows
the consequence, `p0 + p1` is the identity on `I`/`Z` but *doubles* `X`/`Y` where
completeness `Π₀ + Π₁ = 1` forces the dephasing channel. The new kernel zeros
those coefficients in place and leaves their keys present until caller-driven
`reduce()`, preserving the engine's explicit-reduction contract.

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
