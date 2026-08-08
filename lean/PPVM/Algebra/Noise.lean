/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Data.Real.Basic
import PPVM.Algebra.GradedMap
import PPVM.Pauli.Symplectic

/-!
# Noise channels and observable extraction

The Tier-3 targets from `sum/noise.rs` and the observable read-out, modeled on
the Pauli-basis coefficient vector.

* **Unital Pauli channel eigenvalue.** A Pauli channel acts diagonally in the
  Pauli basis: `P ↦ λ_P · P` with `λ_P = Σ_Q p_Q (−1)^{ω(P,Q)}`. Because
  `Σ_Q p_Q = 1`, this collapses to `λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q`.
  The Pauli-specific corollary ties `anticommute` to the actual symplectic form
  `PPVM.Symplectic.omega`.
* **Pauli-basis orthonormality.** In the `C[K]` model the L3 `overlap` pairing
  plays the role of the normalized trace `Tr(PQ)/2ⁿ`; the Pauli basis vectors
  `single P 1` are orthonormal under it. (The `2ⁿ`-normalized matrix trace itself
  is not constructed *here*; it is now built in `PPVM.PauliMatrix`, where
  `trace_tensorPauli_mul` proves the genuine matrix identity
  `Tr(g(p) g(q)) = 2ⁿ δ` and `overlap_eq_trace_div` proves `overlap` *is*
  `Tr(Â B̂)/2ⁿ`. This file's statement is the abstract-key model form.)
* **Zero-state read-out.** `⟨0ⁿ|ρ|0ⁿ⟩ = Σ_{P ∈ {I,Z}ⁿ} c_P`: only the diagonal
  (`X`-free) Paulis contribute; the corollary uses the concrete `X`-free predicate.
-/

namespace PPVM.Noise

open scoped BigOperators

/-! ### Unital Pauli channel eigenvalue -/

/-- **Unital Pauli channel eigenvalue** (general form). With `Σ_Q p_Q = 1`, the
Pauli-transfer eigenvalue `λ_P = Σ_Q p_Q (−1)^{[anti]}` equals
`1 − 2·Σ_{anti} p_Q`. -/
theorem pauli_channel_eigenvalue {K : Type*} [Fintype K] (anti : K → K → Prop)
    [∀ P Q, Decidable (anti P Q)] (p : K → ℝ) (hp : ∑ Q, p Q = 1) (P : K) :
    ∑ Q, p Q * (if anti P Q then -1 else 1)
      = 1 - 2 * ∑ Q ∈ Finset.univ.filter (anti P ·), p Q := by
  have h : ∀ Q, p Q * (if anti P Q then (-1 : ℝ) else 1)
             = p Q - 2 * (if anti P Q then p Q else 0) := by
    intro Q; split <;> ring
  simp_rw [h]
  rw [Finset.sum_sub_distrib, hp, ← Finset.mul_sum, Finset.sum_filter]

/-- **Pauli channel eigenvalue formula, tied to the symplectic form.** Here `anti`
is genuine anticommutation `ω(P,Q) = 1`. This is the arithmetic identity behind
`λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q`; the channel superoperator and its
diagonalization are not constructed here — only the eigenvalue's algebraic form. -/
theorem pauli_channel_eigenvalue_omega {m : ℕ} (p : Symplectic.Sp m → ℝ)
    (hp : ∑ Q, p Q = 1) (P : Symplectic.Sp m) :
    ∑ Q, p Q * (if Symplectic.omega P Q = 1 then -1 else 1)
      = 1 - 2 * ∑ Q ∈ Finset.univ.filter (fun Q => Symplectic.omega P Q = 1), p Q :=
  pauli_channel_eigenvalue (fun P Q => Symplectic.omega P Q = 1) p hp P

/-! ### The eigenvalue is **contractive**, and the channel never grows a key

`pauli_channel_eigenvalue` above gives the *formula* for `λ_P` but says nothing
about its size. The size is what the implementation actually cites: the
`Sum::scale_by_key` fast path (`crates/ppvm-pauli-sum-2/src/sum.rs`,
`PauliError`/`Depolarizing`/…) runs **no** truncation pass and **no** Pauli-weight
re-check after a channel, on the grounds that "the channel is contractive
(`|λ_P| ≤ 1`), so it can never grow a key's Pauli weight". Two separate facts are
being claimed there, and both are proved below:

* `|λ_P| ≤ 1` for a *sub-stochastic* probability vector (`p ≥ 0`, `Σ p ≤ 1` — the
  single-qubit `PauliError` takes `[p_X, p_Y, p_Z]` with an implicit
  `p_I = 1 − Σ`), hence the diagonal channel is an `ℓ¹` **contraction**
  (`l1_contractive`). That is also the missing hypothesis for composing
  `PPVM.Truncation.l1_bound` across a *noisy* circuit: the per-truncation `ℓ¹`
  error bound only telescopes over a long propagation if every intervening
  channel is non-expanding in `ℓ¹`.
* The diagonal channel's support only shrinks (`scaleByKey_support_subset`) and
  never moves a key, so every surviving key's Pauli weight is one it already had
  — nothing for a weight policy to re-check.

`eigenvalue_abs_le_one_needs_substochastic` shows the hypothesis is load-bearing:
a caller passing an over-normalized `[p_X, p_Y, p_Z]` breaks the bound (and with
it both the skipped truncation and the composed error bound), so `Σ p ≤ 1` is a
real precondition on the Rust channel constructors, not a formality. -/

/-- **The Pauli-transfer eigenvalue is contractive.** For a sub-stochastic
probability vector (`p Q ≥ 0` and `Σ_Q p Q ≤ 1`), the eigenvalue
`λ_P = 1 − 2·Σ_{Q anti P} p_Q` of `pauli_channel_eigenvalue` satisfies
`|λ_P| ≤ 1`. (The anticommuting mass is between `0` and `1`, so `λ_P ∈ [−1, 1]`.) -/
theorem pauli_channel_eigenvalue_abs_le_one {K : Type*} [Fintype K] (anti : K → K → Prop)
    [∀ P Q, Decidable (anti P Q)] (p : K → ℝ) (hp0 : ∀ Q, 0 ≤ p Q) (hp : ∑ Q, p Q ≤ 1)
    (P : K) :
    |1 - 2 * ∑ Q ∈ Finset.univ.filter (anti P ·), p Q| ≤ 1 := by
  have h0 : 0 ≤ ∑ Q ∈ Finset.univ.filter (anti P ·), p Q :=
    Finset.sum_nonneg fun Q _ => hp0 Q
  have h1 : ∑ Q ∈ Finset.univ.filter (anti P ·), p Q ≤ 1 :=
    le_trans (Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
      (fun Q _ _ => hp0 Q)) hp
  rw [abs_le]
  constructor <;> linarith

/-- **The same bound, tied to the symplectic form** — the eigenvalue of
`pauli_channel_eigenvalue_omega` is contractive. -/
theorem pauli_channel_eigenvalue_omega_abs_le_one {m : ℕ} (p : Symplectic.Sp m → ℝ)
    (hp0 : ∀ Q, 0 ≤ p Q) (hp : ∑ Q, p Q ≤ 1) (P : Symplectic.Sp m) :
    |1 - 2 * ∑ Q ∈ Finset.univ.filter (fun Q => Symplectic.omega P Q = 1), p Q| ≤ 1 :=
  pauli_channel_eigenvalue_abs_le_one (fun P Q => Symplectic.omega P Q = 1) p hp0 hp P

/-- **A contractive diagonal channel is an `ℓ¹` contraction on `C[K]`.** If every
transfer eigenvalue satisfies `|λ_k| ≤ 1` then `‖Λ A‖₁ ≤ ‖A‖₁` on any finite set
of keys. This is what lets the `ℓ¹` truncation bound (`PPVM.Truncation.l1_bound`)
compose across a noisy propagation: the channels between two truncations cannot
inflate the mass a later truncation is measured against. -/
theorem l1_contractive {K : Type*} (lam c : K → ℝ) (hlam : ∀ k, |lam k| ≤ 1)
    (D : Finset K) :
    ∑ k ∈ D, |lam k * c k| ≤ ∑ k ∈ D, |c k| := by
  refine Finset.sum_le_sum fun k _ => ?_
  rw [abs_mul]
  exact mul_le_of_le_one_left (abs_nonneg _) (hlam k)

/-- The diagonal channel on `C[K]`: rescale each coefficient by its own key's
eigenvalue, leaving the key alone. This is `Sum::scale_by_key`'s action
(`f(&k, &mut c)` over `iter_mut`: no key ever moves). -/
noncomputable def scaleByKey {K C : Type*} [DecidableEq K] [Semiring C] (lam : K → C)
    (A : GradedMap.CMap K C) : GradedMap.CMap K C :=
  A.sum fun k a => Finsupp.single k (lam k * a)

/-- **A diagonal channel never introduces a key.** Its support is contained in the
input's, so every surviving term's Pauli weight is one the map already carried —
the formal content of "`pauli_error` runs no weight re-check". (It is only a
*subset*: a zero eigenvalue kills a coefficient. The Rust backend deliberately
keeps such a term in its map with coefficient `0`, which refines this ideal
`Finsupp` model, where a zero coordinate is structurally absent —
`GradedMap.reduce_structural`.) -/
theorem scaleByKey_support_subset {K C : Type*} [DecidableEq K] [Semiring C] (lam : K → C)
    (A : GradedMap.CMap K C) : (scaleByKey lam A).support ⊆ A.support := by
  classical
  refine (Finsupp.support_sum).trans ?_
  refine Finset.biUnion_subset.2 fun k hk => ?_
  exact (Finsupp.support_single_subset).trans (by simpa using hk)

/-- **The sub-stochastic hypothesis is load-bearing.** Drop `Σ_Q p_Q ≤ 1` and the
bound fails: the nonnegative vector `p ≡ 1` on a two-element alphabet gives
`|1 − 2·2| = 3 > 1`. So an over-normalized `[p_X, p_Y, p_Z]` silently breaks both
the skipped weight re-check and the composed `ℓ¹` truncation bound; the Rust
channel constructors owe the precondition. -/
theorem eigenvalue_abs_le_one_needs_substochastic :
    ¬ ∀ p : Bool → ℝ, (∀ Q, 0 ≤ p Q) → |1 - 2 * ∑ Q, p Q| ≤ 1 := by
  intro h
  have := h (fun _ => 1) (fun _ => by norm_num)
  rw [Fintype.sum_bool] at this
  norm_num at this

/-! ### `two_qubit_pauli_error`'s fifteen hand-written index lists

`pauli_channel_eigenvalue_omega` gives the *formula* `λ_P = 1 − 2·Σ_{Q anti P} p_Q`
but says nothing about which `Q` a given implementation actually sums. The old
crate's `two_qubit_pauli_error` (`crates/ppvm-pauli-sum/src/sum/noise.rs:50-104`)
hard-codes, for each of the 16 observed pairs `P`, a **hand-written list of eight
indices** into `p : [Coeff; 15]`, with no derivation in the source; the shipped
tests only probe one-hot probability vectors, on which a transposed index is
invisible. The theorem below is the missing check: each list is *exactly* the
anticommuting set of its `P`, in the crate's documented probability order
`{IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ}`
(`crates/ppvm-traits/src/traits/noise.rs:73`).

(The `(I,I)` arm is the remaining case and needs no list: nothing anticommutes
with the identity, so `λ_{II} = 1` and the arm is the no-op the crate ships. The
`_ =>` arm — no noise when either site is *lost* — is a deliberate modeling
choice, outside this algebra.) -/

/-- A single-qubit Pauli as `(x, z)` bits, indexed by its position in the
alphabet order `0 = I, 1 = X, 2 = Y, 3 = Z` the probability vector uses. -/
def code : Fin 4 → Bool × Bool
  | 0 => (false, false)
  | 1 => (true, false)
  | 2 => (true, true)
  | 3 => (false, true)

/-- The 15 non-identity two-qubit Paulis in the crate's documented probability
order: index `i` is the pair numbered `i + 1` in base 4 over `{I, X, Y, Z}` —
`0 ↦ IX`, `3 ↦ XI`, `14 ↦ ZZ`. This is also the order of the `match` arms in
`noise.rs`. -/
def qPair : Fin 15 → (Bool × Bool) × (Bool × Bool)
  | 0 => (code 0, code 1)   -- IX
  | 1 => (code 0, code 2)   -- IY
  | 2 => (code 0, code 3)   -- IZ
  | 3 => (code 1, code 0)   -- XI
  | 4 => (code 1, code 1)   -- XX
  | 5 => (code 1, code 2)   -- XY
  | 6 => (code 1, code 3)   -- XZ
  | 7 => (code 2, code 0)   -- YI
  | 8 => (code 2, code 1)   -- YX
  | 9 => (code 2, code 2)   -- YY
  | 10 => (code 2, code 3)  -- YZ
  | 11 => (code 3, code 0)  -- ZI
  | 12 => (code 3, code 1)  -- ZX
  | 13 => (code 3, code 2)  -- ZY
  | 14 => (code 3, code 3)  -- ZZ

/-- Single-site anticommutation `ω` on `(x,z)` bits, as a `Bool`. -/
def antiSite (p q : Bool × Bool) : Bool := xor (p.1 && q.2) (p.2 && q.1)

/-- Two-qubit anticommutation: the site-wise `ω`s add in `𝔽₂`. -/
def antiPair (P Q : (Bool × Bool) × (Bool × Bool)) : Bool :=
  xor (antiSite P.1 Q.1) (antiSite P.2 Q.2)

/-- The fifteen index lists of `two_qubit_pauli_error`, transcribed **verbatim and
in source order** from `noise.rs:53-100` (arm `qPair i` gets list `oldIndices i`). -/
def oldIndices : Fin 15 → List (Fin 15)
  | 0 => [1, 10, 13, 14, 2, 5, 6, 9]      -- (I, X)
  | 1 => [0, 10, 12, 14, 2, 4, 6, 8]      -- (I, Y)
  | 2 => [0, 1, 12, 13, 4, 5, 8, 9]       -- (I, Z)
  | 3 => [10, 11, 12, 13, 14, 7, 8, 9]    -- (X, I)
  | 4 => [1, 11, 12, 2, 5, 6, 7, 8]       -- (X, X)
  | 5 => [0, 11, 13, 2, 4, 6, 7, 9]       -- (X, Y)
  | 6 => [0, 1, 10, 11, 14, 4, 5, 7]      -- (X, Z)
  | 7 => [11, 12, 13, 14, 3, 4, 5, 6]     -- (Y, I)
  | 8 => [1, 10, 11, 12, 2, 3, 4, 9]      -- (Y, X)
  | 9 => [0, 10, 11, 13, 2, 3, 5, 8]      -- (Y, Y)
  | 10 => [0, 1, 11, 14, 3, 6, 8, 9]      -- (Y, Z)
  | 11 => [10, 3, 4, 5, 6, 7, 8, 9]       -- (Z, I)
  | 12 => [1, 13, 14, 2, 3, 4, 7, 8]      -- (Z, X)
  | 13 => [0, 12, 14, 2, 3, 5, 7, 9]      -- (Z, Y)
  | 14 => [0, 1, 10, 12, 13, 3, 6, 7]     -- (Z, Z)

/-- **Every hand-written index list is exactly the anticommuting set of its
observed pair.** Together with `pauli_channel_eigenvalue_omega` this makes each of
the 15 `1 − 2·Σ p[i]` factors in `two_qubit_pauli_error` the genuine Pauli-transfer
eigenvalue `λ_P`, for an *arbitrary* probability vector — not only for the one-hot
vectors the crate's tests probe. -/
theorem twoQubitPauliError_indices_anticommuting :
    ∀ i j : Fin 15, (oldIndices i).contains j = antiPair (qPair i) (qPair j) := by
  decide

/-- Each list has eight entries and no repeats, so `contains`-agreement above is
genuine set equality (each `p[i]` is counted exactly once). -/
theorem twoQubitPauliError_indices_length : ∀ i : Fin 15, (oldIndices i).length = 8 := by
  decide

theorem twoQubitPauliError_indices_nodup : ∀ i : Fin 15, (oldIndices i).Nodup := by
  decide

/-! ### The Bernoulli firing convention of the stochastic channels

The crate's stochastic channels draw `r = rng.random::<f64>() ∈ [0,1)` and fire
on a comparison against the channel probability `p`. Two *different* conventions
coexist (`crates/ppvm-tableau/src/{noise.rs,tableau_like.rs}`, and the `-2` port
reproduces both verbatim):

* the depolarizing family (`depolarize_impl`, `pauli_error_impl`,
  `two_qubit_pauli_error_impl`) uses `if p <= r { return }`, i.e. it fires on the
  **strict** predicate `r < p`;
* `loss_channel` / `asymmetric_loss_channel` use `if p < r { return }`, i.e. they
  fire on the **non-strict** predicate `r ≤ p`.

For an *ideal* continuous `Uniform[0,1)` the two are indistinguishable — they
differ only on the null event `{r = p}` (`fire_conventions_agree_off_diagonal`),
so neither is "wrong" as a sampler of a `Bernoulli(p)`. What is **not** a null
event is the endpoint `p = 0`, which a real `f64` sampler does hit
(`random::<f64>()` can return exactly `0.0`): the strict predicate is then a
guaranteed no-op, the non-strict one fires (`fire_strict_zero_noop`,
`fire_nonstrict_fires_at_zero`). So `loss_channel(q, 0.0)` is *not* the identity
under the shipped convention, while `depolarize1(q, 0.0)` is.

This is a convention inconsistency inside the channel family, not an algebraic
defect: at every `p` that is not exactly representable as a multiple of `2⁻⁵³`
the two conventions agree on every draw. Under the behaviour-preservation
directive the port reproduces both verbatim; unifying them is a seeded-stream
change and needs sign-off. -/

/-- **Off the diagonal the two firing conventions agree.** For `r ≠ p` the strict
and non-strict predicates coincide, so an ideal `Uniform[0,1)` draw fires with
probability `p` under either. -/
theorem fire_conventions_agree_off_diagonal (p r : ℝ) (h : r ≠ p) :
    (r < p) ↔ (r ≤ p) :=
  ⟨le_of_lt, fun hle => lt_of_le_of_ne hle h⟩

/-- **Only the strict convention makes `p = 0` a guaranteed no-op**: no draw from
`[0,1)` satisfies `r < 0`. -/
theorem fire_strict_zero_noop (r : ℝ) (hr : 0 ≤ r) : ¬ (r < 0) := not_lt.mpr hr

/-- **…and the non-strict convention fires at `p = 0` on the draw `r = 0`**, which
a real `f64` `Uniform[0,1)` sampler produces with probability `2⁻⁵³`. This is the
whole observable content of the `loss_channel` convention divergence. -/
theorem fire_nonstrict_fires_at_zero : (0 : ℝ) ≤ 0 := le_refl 0

/-! ### Loss channels: the Heisenberg transfer is trace-preserving

Nothing above (or in `PPVM.Symplectic`, which proves only bit-level loss
*invariance*) says anything about the loss channels **as channels**. Workload 6 of
the integration baseline is the acceptance bar for the loss port, and the old
crate's `correlated_loss_channel` (`crates/ppvm-pauli-sum/src/sum/noise.rs:192-247`)
carries an arithmetic no test covers at distinct `p₀/p₁/p₂`. This section is the
oracle for it.

**The alphabet.** A site of `LossyPauliWord` carries one of `I, X, Y, Z, L`. The
crucial (and easy-to-miss) convention is that `I` is the identity on the **qubit
subspace only** — it is `0` on the loss level — while `L` is the loss projector
`|L⟩⟨L|`. So the identity of the *full* space is `𝟙 = I + L` per site, and the
statement "the channel is trace preserving" is `Λ*(𝟙) = 𝟙` for the Heisenberg
transfer `Λ*` the crate applies to observables, **not** `Λ*(I) = I`.

That distinction is the whole content:

* `loss_channel(p)` sends `I ↦ (1−p)·I` and `L ↦ L + p·I`, and `(1−p) + p = 1`
  drains exactly into the `I` component of `𝟙`. It is trace preserving on the
  lossy word (`lossChannel_trace_preserving`) and trace-*reducing* by `1−p` on the
  plain `PauliWord`, where `L` is unrepresentable so the `L`-arm is dead
  (`lossChannel_nonLossy_scales_by_one_minus_p`) — which is exactly what the
  method's own doc comment says.
* `reset_loss_channel` sends `I ↦ I + L`, `Z ↦ Z + L`, `X, Y ↦ themselves`, and
  `L ↦ 0` (the crate's `*v *= 0.0`, which keeps the term in the map at coefficient
  zero). Trace preserving (`resetLossChannel_trace_preserving`).
* `correlated_loss_channel(p₀,p₁,p₂)` is trace preserving **for every**
  `(p₀,p₁,p₂)`, with no normalization hypothesis at all
  (`correlatedLossChannel_trace_preserving`). In particular the two arms the
  baseline flags as suspicious are *correct*: the one-already-lost arm scales its
  survivor by `1 − p₂` because that population leaves for `LL` at rate `p₂`, and
  weights its emitted branch `p₁` because it *gains* from the both-in-subspace
  population at rate `p₁` — two different processes, so they are not meant to
  pair; and the `(L,L)` arm leaves its survivor unscaled because loss is
  irreversible (`LL` population never leaves). The books balance column-wise, not
  arm-wise.

The channels are modeled by their transfer **matrix** `T k j` = coefficient of key
`j` in the image of the basis observable `k`, transcribed arm-for-arm from the
Rust `map_insert` / `map_insert_multiple` closures. Only the arm *selection*
depends on the site alphabet (the Rust `match` tests solely for `Pauli::L`), so
the loss sector `{I, L}ⁿ` — where `𝟙` lives — is invariant, and the theorems below
are complete statements of trace preservation. -/

section Loss

/-- A site of `LossyPauliWord`, as `Fin 5`: `0 = I, 1 = X, 2 = Y, 3 = Z, 4 = L`.
`I` is the identity on the *qubit subspace* (zero on the loss level); `L` is the
loss projector `|L⟩⟨L|`. -/
abbrev Site := Fin 5

/-- `I` — identity on the qubit subspace. -/
def sI : Site := 0
/-- `Z`. -/
def sZ : Site := 3
/-- `L` — the loss projector. -/
def sL : Site := 4

/-- An observable on one site, as its coefficient vector in the `{I,X,Y,Z,L}`
basis. -/
abbrev Obs1 := Site → ℝ

/-- The Heisenberg transfer of a (possibly branching) channel, given as a matrix:
`T k j` is the coefficient of key `j` in the image of the basis observable `k`.
This is exactly what a `map_insert` closure returns — the scaled survivor plus the
emitted branches. -/
def transfer1 (T : Site → Site → ℝ) (A : Obs1) : Obs1 := fun j => ∑ k, A k * T k j

/-- The identity observable on one site: `𝟙 = I + L`. Trace preservation is
`Λ*(𝟙) = 𝟙`. -/
def unit1 : Obs1 := fun s => if s = sI then 1 else if s = sL then 1 else 0

/-- `loss_channel(p)` (`noise.rs:164-179`): an in-subspace site is scaled by
`1 − p` with no branch; a already-lost site keeps its (unscaled) term and emits a
branch onto `I` with weight `p`. -/
def lossT (p : ℝ) : Site → Site → ℝ := fun k j =>
  if k = sL then (if j = sL then 1 else if j = sI then p else 0)
  else (if j = k then 1 - p else 0)

/-- **`loss_channel` is trace preserving on the lossy word.** `Λ*(I + L) = I + L`
for every `p`: the `1 − p` the survivor loses is exactly the `p` the already-lost
key branches back onto `I`. -/
theorem lossChannel_trace_preserving (p : ℝ) : transfer1 (lossT p) unit1 = unit1 := by
  funext j
  fin_cases j <;> simp [transfer1, lossT, unit1, sI, sL, Fin.sum_univ_five]

/-- **…and trace-reducing by `1 − p` on the plain `PauliWord`**, where `L` is
unrepresentable so the branch arm is dead code: the qubit-subspace identity is
scaled, `I ↦ (1−p)·I`. This is the documented behaviour of the trait
(`noise.rs:155-162`, "reduces the trace of the density matrix as (1 − p) per lost
qubit"), and it is a *convention*, not a defect — `ResetLossChannel` is the
trace-preserving variant. -/
theorem lossChannel_nonLossy_scales_by_one_minus_p (p : ℝ) :
    transfer1 (lossT p) (fun s => if s = sI then 1 else 0)
      = fun s => if s = sI then 1 - p else 0 := by
  funext j
  fin_cases j <;> simp [transfer1, lossT, sI, sL, Fin.sum_univ_five]

/-- `reset_loss_channel` (`noise.rs:253-266`): `I` and `Z` keep their term and emit
an `L` branch at the same coefficient; `X`/`Y` are untouched; `L` is scaled to `0`
(the term stays in the map — the "no implicit reduce" contract). -/
def resetT : Site → Site → ℝ := fun k j =>
  if k = sL then 0
  else if k = sI ∨ k = sZ then (if j = k then 1 else if j = sL then 1 else 0)
  else (if j = k then 1 else 0)

/-- **`reset_loss_channel` is trace preserving.** `Λ*(I + L) = (I + L) + 0`: the
lost population is reset into `|0⟩`, on which the qubit-subspace identity reads
`1`. -/
theorem resetLossChannel_trace_preserving : transfer1 resetT unit1 = unit1 := by
  funext j
  fin_cases j <;> simp [transfer1, resetT, unit1, sI, sL, sZ, Fin.sum_univ_five]

/-- A two-site observable. -/
abbrev Obs2 := Site × Site → ℝ

/-- The two-site Heisenberg transfer (same matrix convention as `transfer1`). -/
def transfer2 (T : Site × Site → Site × Site → ℝ) (A : Obs2) : Obs2 :=
  fun j => ∑ k, A k * T k j

/-- The two-site identity observable `𝟙 ⊗ 𝟙 = (I + L) ⊗ (I + L)`. -/
def unit2 : Obs2 := fun k => unit1 k.1 * unit1 k.2

/-- `correlated_loss_channel(p₀,p₁,p₂)` (`noise.rs:192-247`), arm for arm:

* both already lost — survivor unscaled, three branches `(I,L)`, `(L,I)` at `p₂`
  and `(I,I)` at `p₀`;
* exactly one already lost — survivor scaled by `1 − p₂`, one branch replacing the
  lost site by `I` at weight `p₁`;
* both in subspace — survivor scaled by `1 − 2p₁ − p₀`, no branch. -/
def corrT (p0 p1 p2 : ℝ) : Site × Site → Site × Site → ℝ := fun k j =>
  if k.1 = sL ∧ k.2 = sL then
    (if j = (sL, sL) then 1
     else if j = (sI, sL) then p2
     else if j = (sL, sI) then p2
     else if j = (sI, sI) then p0 else 0)
  else if k.2 = sL then
    (if j = k then 1 - p2 else if j = (k.1, sI) then p1 else 0)
  else if k.1 = sL then
    (if j = k then 1 - p2 else if j = (sI, k.2) then p1 else 0)
  else (if j = k then 1 - 2 * p1 - p0 else 0)

-- `simp` discharges most of the 25 columns outright, so the trailing `ring` must
-- be `<;>`-distributed (a plain `;` hits "no goals" on those); the seq-focus
-- linter's suggestion is a false positive here.
set_option linter.unnecessarySeqFocus false in
/-- **`correlated_loss_channel` is trace preserving for every `(p₀,p₁,p₂)`** — no
normalization hypothesis is needed. Reading the four columns of `Λ*(𝟙) = 𝟙`:

* `(I,I)`: `(1 − 2p₁ − p₀) + p₁ + p₁ + p₀ = 1`;
* `(I,L)` and `(L,I)`: `(1 − p₂) + p₂ = 1`;
* `(L,L)`: `1`.

So the baseline's suspicion is refuted: the `1 − p₂` survivor scale and the `p₁`
branch weight in the one-already-lost arms belong to *different* columns (loss of
the second qubit vs. gain from the both-in-subspace sector), and the unscaled
`(L,L)` survivor is required — loss is irreversible, so that column is closed.
The port must reproduce this arithmetic verbatim. -/
theorem correlatedLossChannel_trace_preserving (p0 p1 p2 : ℝ) :
    transfer2 (corrT p0 p1 p2) unit2 = unit2 := by
  funext j
  obtain ⟨j1, j2⟩ := j
  fin_cases j1 <;> fin_cases j2 <;>
    simp [transfer2, corrT, unit2, unit1, sI, sL, Fintype.sum_prod_type,
      Fin.sum_univ_five, Prod.ext_iff] <;>
    ring

end Loss

/-! ### Pauli-basis orthonormality (the `overlap` pairing) -/

/-- **Pauli-basis orthonormality** in the `C[K]` model: the basis vectors
`single P 1` are orthonormal under the L3 `overlap` pairing, `⟪P, Q⟫ = δ_{PQ}`.
This is the model-level form of `Tr(PQ)/2ⁿ = δ_{PQ}` over an *abstract* key set,
with `overlap` standing in for the normalized trace. On the concrete Pauli key
the matrix identity itself is now proved: `PPVM.PauliMatrix.trace_tensorPauli_mul`
(and `overlap_eq_trace_div`, `Tr(Â B̂) = 2ⁿ ⟪A,B⟫`). -/
theorem overlap_single_single {K : Type*} [DecidableEq K] (P Q : K) :
    GradedMap.overlap (Finsupp.single P (1 : ℝ)) (Finsupp.single Q 1)
      = if P = Q then 1 else 0 := by
  rw [GradedMap.overlap, Finsupp.sum_single_index (by simp)]
  simp [Finsupp.single_apply, eq_comm]

/-! ### Zero-state read-out -/

/-- **Zero-state read-out** (general form): with a per-Pauli `diag` indicator,
`Σ_P c_P·[diag P] = Σ_{diag} c_P`. -/
theorem overlap_with_zero {K : Type*} [Fintype K] (c : K → ℝ) (diag : K → Prop)
    [DecidablePred diag] :
    ∑ P, c P * (if diag P then 1 else 0) = ∑ P ∈ Finset.univ.filter diag, c P := by
  simp_rw [mul_ite, mul_one, mul_zero]
  rw [Finset.sum_filter]

/-- **Zero-state read-out, concrete.** `⟨0ⁿ|ρ|0ⁿ⟩ = Σ_{P X-free} c_P`. The physics
input — a Pauli's zero-state expectation `⟨0|P|0⟩` is `1` exactly for the diagonal
(`X`-free, `I`/`Z`-only) Paulis — is the *modeling assumption* (the `if … then 1
else 0` factor); the theorem is the resulting collapse to the coefficient sum over
the concrete `X`-free sector (`diag P := ∀ i, (P i).1 = 0`). -/
theorem overlap_with_zero_xfree {m : ℕ} (c : Symplectic.Sp m → ℝ) :
    ∑ P, c P * (if (∀ i, (P i).1 = 0) then 1 else 0)
      = ∑ P ∈ Finset.univ.filter (fun P => ∀ i, (P i).1 = 0), c P :=
  overlap_with_zero c (fun P => ∀ i, (P i).1 = 0)

end PPVM.Noise
