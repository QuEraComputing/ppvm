/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Algebra.MonoidAlgebra.Basic
import Mathlib.Analysis.Complex.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import PPVM.Algebra.GradedMap
import PPVM.Algebra.Truncation
import PPVM.Algebra.Twisted
import PPVM.Instantiations.Rotation

/-!
# The symbolic coefficient ring `ℝ[sᵢ, cᵢ]`

`ppvm-sym-2`'s `Term` instantiates `ppvm_traits_2::Coefficient` (and
`Angle<Term>`) with a **free polynomial ring** in the formal symbols
`sᵢ = sin(xᵢ)` and `cᵢ = cos(xᵢ)`: a `Prod` is a monomial (per-variable sine and
cosine exponents, with cached degree totals `sin_pow`/`cos_pow`), and a `Sum` is
an `FxHashMap<Prod, f64>` plus a constant `c₀`. That is exactly
`AddMonoidAlgebra ℝ (ℕ →₀ ℕ × ℕ)`, whose underlying module is the
`GradedMap.CMap` this repo already models.

Two facts about that ring are load-bearing in the Rust and had no oracle.

## 1. The sine degree is additive, so `max_sin` truncation is a monomial ideal

`Sum::add_term` drops a produced monomial the instant `p.sin_pow() > max`
(drop-at-insert, never materialized), and `Sum::mul_term` short-circuits the
*whole* table to empty in `O(1)` when the multiplier `p` already exceeds the
bound. Both rest on `sinDeg (p·q) = sinDeg p + sinDeg q` (`sinDeg_add`), which
makes `Dₖ = {p | sinDeg p > k}` closed under multiplication — a monomial ideal
(`truncIdeal_mul_left`/`_right`). Consequences, both proved below:

* **drop-at-insert = drop-at-end** (`mulMono_drop_at_insert_eq_drop_at_end`), via
  `GradedMap.batchMap_filter_key`: the degree bound reads only the *key*, so
  filtering during the product loop is exact, not an approximation;
* **the `clear()` shortcut is sound** (`mulMono_clear_sound`): if `sinDeg p > k`
  then every monomial of `A·p` is in `Dₖ`, so the retained map is `0`.

And the contrast that keeps the *other* truncation axis honest:
`eps_drop_at_insert_ne_drop_at_end` shows a coefficient-magnitude rule (`min_eps`)
is **not** interchangeable with a post-pass — so `min_eps` must stay where old put
it, inside the accumulation loop.

## 2. Evaluation is a ring hom, and Pythagoras holds only after it

`impl Angle<Term> for Term` drives the engine's rotation kernel over this ring.
The rotation guarantees of `PPVM.Rotation` (`rot_norm_sq`, `rot_rot`) are stated
over `ℝ` and rest on `sin²θ + cos²θ = 1`. In the **free** ring that relation is
false: `sᵢ² + cᵢ² ≠ 1` (`pythagorean_ne_one`), and no `Term`-level operation ever
reduces it. It becomes true only in the quotient by `(sᵢ² + cᵢ² − 1)`, equivalently
only under evaluation at an angle vector (`evalHom_pythagorean`).

The bridge is therefore a *commuting square*: the symbolic rotation on a pair of
coefficients evaluates to the real 2-D rotation, pointwise in `θ`
(`evalHom_symRot`), from which norm preservation transfers
(`symRot_norm_sq_after_eval`) — while `symRot_norm_sq_ne_symbolically` witnesses
that it genuinely fails before evaluation.

## 3. The ring is `ℤ/4`-**graded**, and the complex evaluation is not injective

`Prod` carries a `phase` byte `k ∈ ℤ/4` which is part of its `Hash`/`Eq`, so the
ring `ppvm-sym-2` implements is not `SymRing` but the phase-graded
`PhasedSymRing = AddMonoidAlgebra ℝ (Mono × ZMod 4)` (§"The `ℤ/4` phase grading"
below). `Term::eval_complex` is the algebra hom `evalC θ` into `ℂ` that sends the
grading to `iᵏ` (`evalC_mul`, built on `Twisted.iPow_add` and `monoValue_add`),
and it is **not injective** (`evalC_not_injective`): representational equality on
this ring is strictly finer than denotational equality in `ℂ`. That is exactly
the documented `ImaginaryUnit` law exemption in `ppvm-sym-2/src/coeff.rs` —
`i·i` is the key `(1, 2)` while `−one()` is the key `(1, 0)` scaled by `−1`, two
distinct monomials with the same value (`iSym_sq_ne_neg_one`,
`evalC_iSym_sq_eq_neg_one`) — and it is also why `i²·p` never cancels against
`−p` in the table, so the `min_eps` rule thresholds the two summands
independently. `conjSym` (`Conjugate for Term`) is the phase-negating ring
involution, with `evalC ∘ conj = star ∘ evalC` (`evalC_conjSym`) — the ring-level
form of `Pauli/Matrix.lean`'s `star_iU`.

## 4. The `min_eps` arm of the `clear()` shortcut is coarser, but `ℓ¹`-bounded

`Sum::mul_term` clears the whole table when `|coeff| < min_eps`, standing in for
the per-monomial rule in `Sum::add_term`. Unlike the degree arm
(`mulMono_clear_sound`, an *equality*) this one genuinely discards monomials the
per-monomial rule keeps (`epsClear_ne_retain_pointwise`), so it is an
over-truncation, not an exact shortcut. What licenses keeping it is the `ℓ¹`
bound `epsClear_l1_eq` / `epsClear_l1_lt` and its observable-error corollary
`epsClear_error_lt`, stated against `PPVM.Truncation.l1_bound`.
-/

namespace PPVM.Symbolic

open PPVM.GradedMap

/-! ### Monomials and their two degrees -/

/-- A monomial of `ℝ[sᵢ, cᵢ]`: for each variable `i`, the pair
`(sine exponent, cosine exponent)`. This is `ppvm-sym-2`'s `Prod` (a canonical
ascending list of `Factor { var, sin, cos }`), written additively — monomial
*multiplication* is `+` here. -/
abbrev Mono := ℕ →₀ ℕ × ℕ

/-- `Prod::sin_pow` — the total sine degree, the quantity `max_sin` bounds. The
Rust caches it incrementally in a `u32` field so the truncation test is `O(1)`;
`sinDeg_add` is the invariant that cache must maintain. -/
def sinDeg (m : Mono) : ℕ := m.sum fun _ e => e.1

/-- `Prod::cos_pow` — the total cosine degree. -/
def cosDeg (m : Mono) : ℕ := m.sum fun _ e => e.2

@[simp] theorem sinDeg_zero : sinDeg 0 = 0 := rfl

@[simp] theorem cosDeg_zero : cosDeg 0 = 0 := rfl

/-- **The sine degree is additive under monomial multiplication.**
`sinDeg (p · q) = sinDeg p + sinDeg q` — the grading that makes `max_sin` mean
anything, and the invariant the cached `Prod::sin_pow` counter must preserve
across `mul_sin` / `mul_cos` / `MulAssign<Prod>`. -/
theorem sinDeg_add (m₁ m₂ : Mono) : sinDeg (m₁ + m₂) = sinDeg m₁ + sinDeg m₂ :=
  Finsupp.sum_add_index' (fun _ => rfl) (fun _ _ _ => rfl)

/-- **The cosine degree is additive** — the same statement for `Prod::cos_pow`. -/
theorem cosDeg_add (m₁ m₂ : Mono) : cosDeg (m₁ + m₂) = cosDeg m₁ + cosDeg m₂ :=
  Finsupp.sum_add_index' (fun _ => rfl) (fun _ _ _ => rfl)

/-- Multiplying can only raise the sine degree (there are no inverses). -/
theorem sinDeg_le_mul (m₁ m₂ : Mono) : sinDeg m₁ ≤ sinDeg (m₁ + m₂) := by
  rw [sinDeg_add]; exact Nat.le_add_right _ _

/-- **`Dₖ = {p | sinDeg p > k}` is closed under multiplication on the right** —
so it spans a monomial *ideal*, not merely a set. This is the algebraic fact
`Sum::mul_term`'s whole-table `clear()` fast path is built on. -/
theorem truncIdeal_mul_right {k : ℕ} {p : Mono} (hp : k < sinDeg p) (q : Mono) :
    k < sinDeg (p + q) :=
  lt_of_lt_of_le hp (sinDeg_le_mul p q)

/-- …and on the left (`Mono` is commutative, but both spellings are used). -/
theorem truncIdeal_mul_left {k : ℕ} {p : Mono} (hp : k < sinDeg p) (q : Mono) :
    k < sinDeg (q + p) := by
  rw [add_comm]; exact truncIdeal_mul_right hp q

/-! ### The ring, and its two truncation axes -/

/-- The symbolic coefficient ring `ℝ[sᵢ, cᵢ]` — `ppvm-sym-2`'s `Term`. Its
underlying module is the `GradedMap.CMap Mono ℝ` this repo already models, i.e.
the `FxHashMap<Prod, f64>` (plus `c₀`, the coefficient of the empty monomial). -/
abbrev SymRing := AddMonoidAlgebra ℝ Mono

/-- **The ring product multiplies coefficients and adds monomial exponents** —
`GradedMap.multiply_single` at `K = Mono`. This is what ties the `+` used in the
degree theorems above to the multiplication `Sum::mul_term` actually performs. -/
theorem single_mul_single (m₁ m₂ : Mono) (c d : ℝ) :
    (AddMonoidAlgebra.single m₁ c : SymRing) * AddMonoidAlgebra.single m₂ d
      = AddMonoidAlgebra.single (m₁ + m₂) (c * d) :=
  AddMonoidAlgebra.single_mul_single _ _ _ _

/-- `Sum::add_term`'s `max_sin` keep-rule, as a `Retain` predicate: keep the
monomial `m` iff `sinDeg m ≤ k`. It reads only the **key**. -/
def keepSinDeg (k : ℕ) : Mono → ℝ → Bool := fun m _ => decide (sinDeg m ≤ k)

/-- The batch `Sum::mul_term` produces when multiplying the table by a single
monomial `p`: one term `(m · p, A m)` per stored monomial. -/
noncomputable def mulMonoBatch (p : Mono) (A : CMap Mono ℝ) : Multiset (Mono × ℝ) :=
  A.support.val.map fun m => (m + p, A m)

/-- **Drop-at-insert equals drop-at-end for `max_sin`.** Rejecting each produced
monomial inside the product loop (what `Sum::add_term` does — the rejected
monomial is never materialized, so the table never grows to hold it) yields
exactly the map obtained by forming the full product and truncating afterwards.

An instance of `GradedMap.batchMap_filter_key`; it applies precisely because the
sine-degree bound is a function of the monomial alone. -/
theorem mulMono_drop_at_insert_eq_drop_at_end (k : ℕ) (p : Mono) (A : CMap Mono ℝ) :
    batchMap ((mulMonoBatch p A).filter fun t => sinDeg t.1 ≤ k)
      = retain (keepSinDeg k) (batchMap (mulMonoBatch p A)) :=
  batchMap_filter_key (fun m => sinDeg m ≤ k) _

/-- **The whole-sum-to-zero `clear()` shortcut is sound.** If the multiplier `p`
already exceeds the bound, then *every* monomial of `A · p` lies in the truncation
ideal `Dₖ` (`truncIdeal_mul_left`), so the truncated product is the zero map — and
`Sum::mul_term` may clear the entire table in `O(1)` instead of walking it. -/
theorem mulMono_clear_sound {k : ℕ} {p : Mono} (hp : k < sinDeg p) (A : CMap Mono ℝ) :
    batchMap ((mulMonoBatch p A).filter fun t => sinDeg t.1 ≤ k) = 0 := by
  have hnil : (mulMonoBatch p A).filter (fun t => sinDeg t.1 ≤ k) = 0 := by
    rw [Multiset.filter_eq_nil]
    rintro ⟨m, c⟩ hm
    rw [mulMonoBatch, Multiset.mem_map] at hm
    obtain ⟨j, _, hj⟩ := hm
    have : m = j + p := congrArg Prod.fst hj.symm
    subst this
    exact Nat.not_le_of_lt (truncIdeal_mul_left hp j)
  rw [hnil]
  simp [batchMap]

/-- **…and the same statement through the `Retain` keep-rule**: the truncated
product of a sum by an over-degree monomial is the zero map. -/
theorem mulMono_retain_clear {k : ℕ} {p : Mono} (hp : k < sinDeg p) (A : CMap Mono ℝ) :
    retain (keepSinDeg k) (batchMap (mulMonoBatch p A)) = 0 := by
  rw [← mulMono_drop_at_insert_eq_drop_at_end, mulMono_clear_sound hp]

/-- **The `min_eps` axis is *not* interchangeable with a post-pass.** A
coefficient-magnitude rule looks at the value, not the key, so
`batchMap_filter_key` does not apply — and it genuinely fails: two produced terms
on the same key, each below the threshold, are both dropped at insert (result
`0`) but survive as their sum when the rule runs after accumulation (result `2`).

So old's drop-at-insert `min_eps` is a *different function*, not an optimization
of a post-pass: relocating it to a `truncate()` would change results. The
`max_sin` axis, by contrast, may be moved freely
(`mulMono_drop_at_insert_eq_drop_at_end`). -/
theorem eps_drop_at_insert_ne_drop_at_end :
    batchMap ((({(true, (1 : ℤ)), (true, 1)} : Multiset (Bool × ℤ))).filter
        fun t => 2 ≤ |t.2|)
      ≠ retain (fun _ c => decide (2 ≤ |c|))
        (batchMap ({(true, (1 : ℤ)), (true, 1)} : Multiset (Bool × ℤ))) := by
  classical
  have hfilter :
      (({(true, (1 : ℤ)), (true, 1)} : Multiset (Bool × ℤ)).filter fun t => 2 ≤ |t.2|) = 0 := by
    norm_num [Multiset.filter_cons]
  rw [hfilter]
  intro h
  have hval := DFunLike.congr_fun h true
  norm_num [batchMap, Finsupp.single_apply] at hval

/-! ### Evaluation is a ring homomorphism -/

/-- The value a monomial takes under the substitution `sᵢ ↦ a i`, `cᵢ ↦ b i`. -/
noncomputable def monoValue (a b : ℕ → ℝ) (m : Mono) : ℝ :=
  m.prod fun i e => a i ^ e.1 * b i ^ e.2

@[simp] theorem monoValue_zero (a b : ℕ → ℝ) : monoValue a b 0 = 1 := rfl

/-- **Monomial evaluation is multiplicative** — the exponent vectors add, so the
values multiply. -/
theorem monoValue_add (a b : ℕ → ℝ) (m₁ m₂ : Mono) :
    monoValue a b (m₁ + m₂) = monoValue a b m₁ * monoValue a b m₂ :=
  Finsupp.prod_add_index' (fun _ => by simp) (fun i e₁ e₂ => by
    simp only [Prod.fst_add, Prod.snd_add, pow_add]; ring)

/-- Monomial evaluation as a monoid hom, the datum `AddMonoidAlgebra.lift`
consumes. -/
noncomputable def monoValueHom (a b : ℕ → ℝ) : Multiplicative Mono →* ℝ where
  toFun m := monoValue a b m.toAdd
  map_one' := monoValue_zero a b
  map_mul' m₁ m₂ := monoValue_add a b m₁.toAdd m₂.toAdd

/-- **Substitution is an `ℝ`-algebra homomorphism** `ℝ[sᵢ, cᵢ] → ℝ`, for *any*
choice of values for the symbols. Universality of the free ring: nothing
constrains `a i` and `b i` to be a sine/cosine pair. -/
noncomputable def substHom (a b : ℕ → ℝ) : SymRing →ₐ[ℝ] ℝ :=
  AddMonoidAlgebra.lift ℝ ℝ Mono (monoValueHom a b)

@[simp] theorem substHom_single (a b : ℕ → ℝ) (m : Mono) (c : ℝ) :
    substHom a b (AddMonoidAlgebra.single m c) = c * monoValue a b m := by
  rw [substHom, AddMonoidAlgebra.lift_single, smul_eq_mul]
  rfl

/-- `Term::eval` — evaluation at an angle vector `θ`, i.e. the substitution
`sᵢ ↦ sin (θ i)`, `cᵢ ↦ cos (θ i)`. -/
noncomputable def evalHom (θ : ℕ → ℝ) : SymRing →ₐ[ℝ] ℝ :=
  substHom (fun i => Real.sin (θ i)) (fun i => Real.cos (θ i))

/-- **`eval` is additive** — `ev (p + q) = ev p + ev q`. -/
theorem evalHom_add (θ : ℕ → ℝ) (p q : SymRing) :
    evalHom θ (p + q) = evalHom θ p + evalHom θ q := map_add _ _ _

/-- **`eval` is multiplicative** — `ev (p · q) = ev p · ev q`. Together with
`evalHom_add` (and `map_one`) this is the ring-homomorphism property every
"evaluate the propagated symbolic coefficient" read-out silently assumes. -/
theorem evalHom_mul (θ : ℕ → ℝ) (p q : SymRing) :
    evalHom θ (p * q) = evalHom θ p * evalHom θ q := map_mul _ _ _

/-! ### The Pythagorean relation holds only after evaluation -/

/-- The formal symbol `sᵢ = sin(xᵢ)` — `Term::var i |>.sin()`. -/
noncomputable def sinVar (i : ℕ) : SymRing :=
  AddMonoidAlgebra.single (Finsupp.single i (1, 0)) 1

/-- The formal symbol `cᵢ = cos(xᵢ)` — `Term::var i |>.cos()`. -/
noncomputable def cosVar (i : ℕ) : SymRing :=
  AddMonoidAlgebra.single (Finsupp.single i (0, 1)) 1

@[simp] theorem substHom_sinVar (a b : ℕ → ℝ) (i : ℕ) : substHom a b (sinVar i) = a i := by
  rw [sinVar, substHom_single, monoValue, Finsupp.prod_single_index (by simp)]
  simp

@[simp] theorem substHom_cosVar (a b : ℕ → ℝ) (i : ℕ) : substHom a b (cosVar i) = b i := by
  rw [cosVar, substHom_single, monoValue, Finsupp.prod_single_index (by simp)]
  simp

@[simp] theorem evalHom_sinVar (θ : ℕ → ℝ) (i : ℕ) :
    evalHom θ (sinVar i) = Real.sin (θ i) := substHom_sinVar _ _ i

@[simp] theorem evalHom_cosVar (θ : ℕ → ℝ) (i : ℕ) :
    evalHom θ (cosVar i) = Real.cos (θ i) := substHom_cosVar _ _ i

/-- **`sᵢ² + cᵢ² = 1` is FALSE in the ring `ppvm-sym-2` implements.** The
coefficient ring is the *free* polynomial ring: `sin(x).square() + cos(x).square()`
is a genuine two-monomial `Sum`, and no `Term`-level operation reduces it. The
witness is the substitution `sᵢ ↦ 0`, `cᵢ ↦ 0`, which is a ring hom out of the
free ring (`substHom`) and sends the left side to `0` and the right to `1`.

Equivalently: the Pythagorean relation lives in the quotient by the ideal
`(sᵢ² + cᵢ² − 1)`, not in `ℝ[sᵢ, cᵢ]`. -/
theorem pythagorean_ne_one (i : ℕ) : sinVar i ^ 2 + cosVar i ^ 2 ≠ (1 : SymRing) := by
  intro h
  have := congrArg (substHom (fun _ => 0) (fun _ => 0)) h
  simp only [map_add, map_pow, substHom_sinVar, substHom_cosVar, map_one] at this
  norm_num at this

/-- **…but it holds after evaluation, pointwise in `θ`.** This is the only sense
in which `Term` "knows" `sin² + cos² = 1`: as an identity in the image of
`evalHom`, i.e. in the quotient the evaluation factors through. -/
theorem evalHom_pythagorean (θ : ℕ → ℝ) (i : ℕ) :
    evalHom θ (sinVar i ^ 2 + cosVar i ^ 2) = 1 := by
  simp only [map_add, map_pow, evalHom_sinVar, evalHom_cosVar]
  exact Real.sin_sq_add_cos_sq (θ i)

/-! ### Transferring the rotation guarantees through `eval`

`impl Angle<Term> for Term` instantiates the engine's rotation kernel over this
ring, so the coefficient pair `(c_P, c_{P'})` a rotation updates is a pair of
*polynomials*, not reals. `PPVM.Rotation.rot_norm_sq` / `rot_rot` are stated over
`ℝ` and consume `sin²θ + cos²θ = 1`, which `pythagorean_ne_one` has just shown is
unavailable symbolically. The commuting square below is the honest bridge. -/

/-- The symbolic rotation on a coefficient pair: exactly `PPVM.Rotation.rot`, but
with the formal symbols `sᵢ`, `cᵢ` in place of `sin θ`, `cos θ` — i.e. what the
engine computes when the coefficient domain is `Term` and the angle is
`Term::var i`. -/
noncomputable def symRot (i : ℕ) (v : SymRing × SymRing) : SymRing × SymRing :=
  (cosVar i * v.1 - sinVar i * v.2, sinVar i * v.1 + cosVar i * v.2)

/-- **The bridge: evaluation intertwines the symbolic and the real rotation.**
`ev_θ ∘ symRot i = rot (θ i) ∘ ev_θ`. Every guarantee `PPVM.Rotation` proves about
`rot` therefore transfers to the symbolic coefficient domain **after** evaluation,
pointwise in `θ` — and only then. -/
theorem evalHom_symRot (θ : ℕ → ℝ) (i : ℕ) (v : SymRing × SymRing) :
    (evalHom θ (symRot i v).1, evalHom θ (symRot i v).2)
      = PPVM.Rotation.rot (θ i) (evalHom θ v.1, evalHom θ v.2) := by
  simp [symRot, PPVM.Rotation.rot, map_sub, map_add, map_mul]

/-- **Norm preservation, transferred.** The symbolic rotation preserves the `ℓ²`
norm of the evaluated coefficient pair, for every angle vector. This is the
precise form of the claim `impl Angle<Term> for Term` makes when it cites
`rot_norm_sq`: it is a statement about `eval`s, not about `Term`s. -/
theorem symRot_norm_sq_after_eval (θ : ℕ → ℝ) (i : ℕ) (v : SymRing × SymRing) :
    evalHom θ (symRot i v).1 ^ 2 + evalHom θ (symRot i v).2 ^ 2
      = evalHom θ v.1 ^ 2 + evalHom θ v.2 ^ 2 := by
  have h := PPVM.Rotation.rot_norm_sq (θ i) (evalHom θ v.1, evalHom θ v.2)
  have hpair := evalHom_symRot θ i v
  rw [show evalHom θ (symRot i v).1 = (PPVM.Rotation.rot (θ i)
        (evalHom θ v.1, evalHom θ v.2)).1 from congrArg Prod.fst hpair,
    show evalHom θ (symRot i v).2 = (PPVM.Rotation.rot (θ i)
        (evalHom θ v.1, evalHom θ v.2)).2 from congrArg Prod.snd hpair]
  simpa using h

/-- **…and it fails before evaluation.** Rotating the pure basis pair `(1, 0)`
gives symbolic coefficients whose squared norm is the polynomial `cᵢ² + sᵢ²`,
which is *not* the polynomial `1`. So a `Term`-level "norm" is not conserved by
the rotation kernel, and citing `rot_norm_sq` for the symbolic domain without the
`eval` qualifier is unsound. -/
theorem symRot_norm_sq_ne_symbolically (i : ℕ) :
    (symRot i (1, 0)).1 ^ 2 + (symRot i (1, 0)).2 ^ 2 ≠ (1 : SymRing) := by
  intro h
  refine pythagorean_ne_one i ?_
  rw [← h, symRot]
  simp
  ring

/-! ### The `ℤ/4` phase grading — the ring `Term` *actually* implements

`SymRing` above drops a field `ppvm-sym-2` really stores: `Prod` carries a
`phase : u8` holding `k ∈ ℤ/4` in `iᵏ`, and that byte is part of the derived
`Hash`/`Eq` keyed into the `FxHashMap<Prod, f64>`. So the implemented coefficient
ring is the **`ℤ/4`-graded** algebra `ℝ[sᵢ, cᵢ][ℤ/4]`, and its evaluation lands
in `ℂ`, not `ℝ`.

The distinction is not cosmetic: it is the whole content of the `ImaginaryUnit`
law exemption `ppvm-sym-2/src/coeff.rs` documents in prose. `ImaginaryUnit`
requires `i·i == −one()`; here `i·i` is the monomial `(1, phase 2)` while
`−one()` is `−1` on the monomial `(1, phase 0)` — two *distinct hash keys* with
the same complex value. That is possible exactly because the evaluation
`PhasedSymRing → ℂ` is a surjective-but-not-injective algebra hom, which is what
this section proves. Three consequences ride on it:

* the law exemption is sound (the two sides agree denotationally,
  `evalC_iSym_sq_eq_neg_one`, and differ representationally,
  `iSym_sq_ne_neg_one`);
* `eval_complex` is multiplicative (`evalC_mul`), which is what makes the Rust
  `l4_operator_product_runs_over_the_symbolic_ring` test a valid oracle for
  `X·Y = iZ`;
* truncation on this ring is *coarser* than on `ℂ`: `phaseTwo_cancel_ne_zero`
  exhibits two summands that cancel in `ℂ` but occupy different keys here, so
  `min_eps` thresholds them independently and never sees the cancellation.
-/

section Phased
open PPVM.Twisted

/-- The monomial key `ppvm-sym-2`'s `Prod` really is: the sine/cosine exponent
vector **together with** the phase byte `k ∈ ℤ/4` of `iᵏ`. Both components are
part of `Prod`'s `Hash`/`Eq`. -/
abbrev PMono := Mono × ZMod 4

/-- **The coefficient ring `ppvm-sym-2` implements**: `ℝ[sᵢ, cᵢ]` graded by
`ℤ/4`, i.e. `FxHashMap<Prod, f64>` with `Prod`'s phase byte included in the key.
`SymRing` is its degree-`0` part. -/
abbrev PhasedSymRing := AddMonoidAlgebra ℝ PMono

theorem I_pow_four : Complex.I ^ 4 = 1 := by
  rw [show (4 : ℕ) = 2 * 2 from rfl, pow_mul, Complex.I_sq]
  norm_num

@[simp] theorem iPowI_zero : iPow Complex.I 0 = 1 := by
  rw [iPow, show (0 : ZMod 4).val = 0 from rfl, pow_zero]

@[simp] theorem iPowI_one : iPow Complex.I 1 = Complex.I := by
  rw [iPow, show (1 : ZMod 4).val = 1 from rfl, pow_one]

@[simp] theorem iPowI_two : iPow Complex.I 2 = -1 := by
  rw [iPow, show (2 : ZMod 4).val = 2 from rfl, Complex.I_sq]

@[simp] theorem iPowI_three : iPow Complex.I 3 = -Complex.I := by
  rw [iPow, show (3 : ZMod 4).val = 3 from rfl, show (3 : ℕ) = 2 + 1 from rfl, pow_succ,
    Complex.I_sq, neg_one_mul]

/-- The value a *phased* monomial takes at the angle vector `θ`:
`iᵏ · sin(θ)^… · cos(θ)^…`. This is `Prod::eval_complex`. -/
noncomputable def phasedValue (θ : ℕ → ℝ) (mk : PMono) : ℂ :=
  iPow Complex.I mk.2 *
    ((monoValue (fun i => Real.sin (θ i)) (fun i => Real.cos (θ i)) mk.1 : ℝ) : ℂ)

/-- **Phased-monomial evaluation is multiplicative** — `iPow_add` on the `ℤ/4`
grading (which is where `i⁴ = 1` is spent) times `monoValue_add` on the
exponents. -/
noncomputable def phasedValueHom (θ : ℕ → ℝ) : Multiplicative PMono →* ℂ where
  toFun mk := phasedValue θ mk.toAdd
  map_one' := by
    change phasedValue θ (0 : PMono) = 1
    simp [phasedValue]
  map_mul' a b := by
    change phasedValue θ (a.toAdd + b.toAdd) = phasedValue θ a.toAdd * phasedValue θ b.toAdd
    simp only [phasedValue, Prod.fst_add, Prod.snd_add, iPow_add Complex.I I_pow_four,
      monoValue_add, Complex.ofReal_mul]
    ring

/-- **`Term::eval_complex` — an `ℝ`-algebra homomorphism `PhasedSymRing → ℂ`.**
`single (m, k) c ↦ c · iᵏ · monoValue θ m`, extended linearly (`eval.rs:77, 100,
127`). Being a hom is what licenses reading a propagated symbolic coefficient off
as a complex number at the end of a circuit. -/
noncomputable def evalC (θ : ℕ → ℝ) : PhasedSymRing →ₐ[ℝ] ℂ :=
  AddMonoidAlgebra.lift ℝ ℂ PMono (phasedValueHom θ)

@[simp] theorem evalC_single (θ : ℕ → ℝ) (mk : PMono) (c : ℝ) :
    evalC θ (AddMonoidAlgebra.single mk c) = (c : ℂ) * phasedValue θ mk := by
  rw [evalC, AddMonoidAlgebra.lift_single, Complex.real_smul]
  rfl

/-- **`eval_complex` is additive.** -/
theorem evalC_add (θ : ℕ → ℝ) (x y : PhasedSymRing) :
    evalC θ (x + y) = evalC θ x + evalC θ y := map_add _ _ _

/-- **`eval_complex` is multiplicative.** This is the specification of
`Term::eval_complex` on a product, and the reason the Rust exact-ring operator
test can use it as an oracle for `X·Y = iZ`. -/
theorem evalC_mul (θ : ℕ → ℝ) (x y : PhasedSymRing) :
    evalC θ (x * y) = evalC θ x * evalC θ y := map_mul _ _ _

/-- Coefficient read-out on the graded ring — the `FxHashMap<Prod, f64>` lookup. -/
theorem phased_add_apply (x y : PhasedSymRing) (k : PMono) : (x + y) k = x k + y k := rfl

/-- `i` as `ppvm-sym-2` represents it: the phase-only monomial `i¹`
(`Term::imaginary_unit()` is `Const(1.0).mul_phase(1)`, i.e. `One(Prod{phase: 1},
1.0)` — *not* a numeric constant). -/
noncomputable def iSym : PhasedSymRing := AddMonoidAlgebra.single ((0 : Mono), (1 : ZMod 4)) 1

@[simp] theorem evalC_iSym (θ : ℕ → ℝ) : evalC θ iSym = Complex.I := by
  simp [iSym, phasedValue]

/-- `i·i` folds the two phase bytes: the *key* becomes `phase 2`, and the value
`1` stays put. (`MulAssign<Prod> for Prod` composes the phase — the correction of
`oldSuspectedBugs` #2, which left it untouched.) -/
theorem iSym_sq : iSym * iSym = AddMonoidAlgebra.single ((0 : Mono), (2 : ZMod 4)) 1 := by
  rw [iSym, AddMonoidAlgebra.single_mul_single, mul_one,
    show ((0 : Mono), (1 : ZMod 4)) + ((0 : Mono), (1 : ZMod 4)) = ((0 : Mono), (2 : ZMod 4)) from
      Prod.ext (by simp) (by decide)]

/-- **The `ImaginaryUnit` law fails *representationally*.** `i·i` is the monomial
`(1, phase 2)` while `−one()` is `−1` on `(1, phase 0)`; the phase byte is part of
the key, so the two are different elements of the ring — which is precisely why
`ppvm-sym-2`'s `PartialEq` (representational, behavioural contract 5) rejects
`i * i == -Term::one()`. -/
theorem iSym_sq_ne_neg_one : iSym * iSym ≠ (-1 : PhasedSymRing) := by
  classical
  intro h
  have hval : (iSym * iSym) ((0 : Mono), (0 : ZMod 4))
      = (-1 : PhasedSymRing) ((0 : Mono), (0 : ZMod 4)) := by rw [h]
  rw [iSym_sq, Finsupp.single_apply,
    if_neg (show ¬(((0 : Mono), (2 : ZMod 4)) = ((0 : Mono), (0 : ZMod 4))) from fun hc =>
      absurd (congrArg Prod.snd hc) (show ¬((2 : ZMod 4) = (0 : ZMod 4)) by decide)),
    AddMonoidAlgebra.neg_apply, AddMonoidAlgebra.one_def, Finsupp.single_apply,
    if_pos (show (0 : PMono) = ((0 : Mono), (0 : ZMod 4)) from rfl)] at hval
  norm_num at hval

/-- **…but it holds *denotationally*.** Both sides evaluate to `−1` in `ℂ`, so
the law exemption documented in `ppvm-sym-2/src/coeff.rs` is sound: the
implemented ring surjects onto the complex symbolic ring, it is just not
isomorphic to it. -/
theorem evalC_iSym_sq_eq_neg_one (θ : ℕ → ℝ) :
    evalC θ (iSym * iSym) = evalC θ (-1 : PhasedSymRing) := by
  rw [map_mul, map_neg, map_one, evalC_iSym, Complex.I_mul_I]

/-- **Two summands that cancel in `ℂ` are distinct keys here.** `i²·p + p` is a
genuine two-monomial table — the phase byte separates them — so the accumulator
never sees the cancellation and `min_eps` thresholds each summand on its own
`|c|`. Truncation on this ring is therefore *coarser* than truncation on the
complex values it denotes. -/
theorem phaseTwo_cancel_ne_zero (m : Mono) {c : ℝ} (hc : c ≠ 0) :
    AddMonoidAlgebra.single (m, (2 : ZMod 4)) c
        + AddMonoidAlgebra.single (m, (0 : ZMod 4)) c ≠ (0 : PhasedSymRing) := by
  classical
  intro h
  have hval : ((AddMonoidAlgebra.single (m, (2 : ZMod 4)) c
      + AddMonoidAlgebra.single (m, (0 : ZMod 4)) c : PhasedSymRing)) (m, (0 : ZMod 4))
      = 0 := by rw [h]; rfl
  rw [phased_add_apply, Finsupp.single_apply, Finsupp.single_apply, if_neg (by
      intro hc'
      exact absurd (congrArg Prod.snd hc')
        (show ¬((2 : ZMod 4) = (0 : ZMod 4)) by decide))] at hval
  exact hc (by simpa using hval)

/-- …while their common image in `ℂ` is `0`. -/
theorem evalC_phaseTwo_cancel (θ : ℕ → ℝ) (m : Mono) (c : ℝ) :
    evalC θ (AddMonoidAlgebra.single (m, (2 : ZMod 4)) c
        + AddMonoidAlgebra.single (m, (0 : ZMod 4)) c) = 0 := by
  rw [map_add, evalC_single, evalC_single]
  simp only [phasedValue, iPowI_zero, iPowI_two]
  ring

/-- **`eval_complex` is not injective.** The kernel is non-trivial, so
representational equality on the implemented ring is *strictly finer* than
denotational equality in `ℂ`. This is the machine-checked form of the argument
`ppvm-sym-2/src/coeff.rs` gives in prose for the `ImaginaryUnit` law exemption:
the ring surjects onto but is not isomorphic to the complex symbolic ring, so a
law stated up to `PartialEq` cannot be expected to hold. -/
theorem evalC_not_injective (θ : ℕ → ℝ) : ¬ Function.Injective (evalC θ) := by
  intro hinj
  refine phaseTwo_cancel_ne_zero (0 : Mono) (c := (1 : ℝ)) one_ne_zero (hinj ?_)
  rw [evalC_phaseTwo_cancel, map_zero]

/-! #### `mul_phase` is multiplication by `iᵏ` — including on the constant summand

`Term::mul_phase k` (`ppvm-sym-2/src/coeff.rs`) does not multiply by anything: it
**relabels keys**, adding `k` to every monomial's phase byte. That is only the
right implementation of "multiply by `iᵏ`" if the relabelling really *is* the ring
multiplication by `iᵏ`, and the two differ exactly at the summand old dropped.

`phaseFold_eq_iSym_pow_mul` is that identity, and it is the sole justification for
the one deliberate behaviour divergence from old on this path (`oldSuspectedBugs`
#3): old's `Sum::add_term` short-circuited on `pow() == 0` alone, so `mul_phase`'s
`Sum` arm fed it a phase-only monomial carrying the constant part and the phase was
folded away — `phaseFold_drop_const_ne` is the machine-checked statement that
leaving the constant summand unphased is *not* multiplication by `iᵏ`. The new
`add_term` additionally requires `p.phase() == 0` before folding into `c₀`, which
is precisely `phaseFold_const`. -/

/-- **`Term::mul_phase k` at the key level**: the relabelling `(m, j) ↦ (m, j + k)`
applied to *every* monomial, the constant summand `(0, 0)` included — a
`Finsupp.mapDomain`, so it touches no coefficient. -/
noncomputable def phaseFold (k : ZMod 4) (x : PhasedSymRing) : PhasedSymRing :=
  Finsupp.mapDomain (fun mk => (mk.1, mk.2 + k)) x

@[simp] theorem phaseFold_single (k : ZMod 4) (mk : PMono) (c : ℝ) :
    phaseFold k (AddMonoidAlgebra.single mk c)
      = AddMonoidAlgebra.single (mk.1, mk.2 + k) c :=
  Finsupp.mapDomain_single

theorem phaseFold_zero (k : ZMod 4) : phaseFold k 0 = 0 := Finsupp.mapDomain_zero

theorem phaseFold_add (k : ZMod 4) (x y : PhasedSymRing) :
    phaseFold k (x + y) = phaseFold k x + phaseFold k y := Finsupp.mapDomain_add

/-- `iⁿ` as a key: the phase-only monomial `(0, n mod 4)`. -/
theorem iSym_pow (n : ℕ) :
    iSym ^ n = AddMonoidAlgebra.single ((0 : Mono), ((n : ZMod 4))) 1 := by
  induction n with
  | zero =>
    rw [pow_zero, AddMonoidAlgebra.one_def]
    exact congrArg (fun j => AddMonoidAlgebra.single ((0 : Mono), j) (1 : ℝ)) Nat.cast_zero.symm
  | succ n ih =>
    rw [pow_succ, ih, iSym, AddMonoidAlgebra.single_mul_single, mul_one]
    congr 1
    exact Prod.ext (add_zero _) (by simp only [Prod.snd_add]; push_cast; ring)

theorem iSym_pow_val (k : ZMod 4) :
    iSym ^ k.val = AddMonoidAlgebra.single ((0 : Mono), k) 1 := by
  rw [iSym_pow, show ((k.val : ℕ) : ZMod 4) = k by revert k; decide]

/-- **The phase relabelling *is* multiplication by `iᵏ`.** `mul_phase k x` — a
pure key rewrite touching no coefficient — equals the ring product `iᵏ · x` in
`PhasedSymRing`. This is what makes "phase every summand, including the constant"
an *identity* rather than a plausible-looking choice. -/
theorem phaseFold_eq_iSym_pow_mul (k : ZMod 4) (x : PhasedSymRing) :
    phaseFold k x = iSym ^ k.val * x := by
  induction x using AddMonoidAlgebra.induction_linear with
  | zero => rw [phaseFold_zero, mul_zero]
  | add f g hf hg => rw [phaseFold_add, mul_add, hf, hg]
  | single mk c =>
    rw [phaseFold_single, iSym_pow_val, AddMonoidAlgebra.single_mul_single, one_mul]
    congr 1
    exact Prod.ext (zero_add mk.1).symm (add_comm _ _)

/-- **…hence it scales the complex value by `iᵏ`.** `evalC θ (mul_phase k x) =
iᵏ · evalC θ x`, the read-out form of the identity. -/
theorem evalC_phaseFold (θ : ℕ → ℝ) (k : ZMod 4) (x : PhasedSymRing) :
    evalC θ (phaseFold k x) = iPow Complex.I k * evalC θ x := by
  rw [phaseFold_eq_iSym_pow_mul, map_mul, map_pow, evalC_iSym, iPow]

/-- **The constant summand is phased like every other**: the key `(0, 0)` becomes
`(0, k)`. This is the arm the new `Sum::add_term` keeps out of the `c₀`
short-circuit (`p.pow() == 0 && p.phase() == 0`). -/
theorem phaseFold_const (k : ZMod 4) (c : ℝ) :
    phaseFold k (AddMonoidAlgebra.single ((0 : Mono), (0 : ZMod 4)) c)
      = AddMonoidAlgebra.single ((0 : Mono), k) c := by
  rw [phaseFold_single]
  exact congrArg (fun j => AddMonoidAlgebra.single ((0 : Mono), j) c) (zero_add k)

/-- **Old's behaviour is not multiplication by `iᵏ`.** Leaving the constant
summand on the key `(0, 0)` — what old's `pow() == 0` short-circuit did — differs
from `phaseFold k` for every non-zero phase and non-zero constant. So the
divergence recorded for `oldSuspectedBugs` #3 is forced: old computed a different
function, not a different representation of the same one. -/
theorem phaseFold_drop_const_ne {k : ZMod 4} (hk : k ≠ 0) {c : ℝ} (hc : c ≠ 0) :
    phaseFold k (AddMonoidAlgebra.single ((0 : Mono), (0 : ZMod 4)) c)
      ≠ AddMonoidAlgebra.single ((0 : Mono), (0 : ZMod 4)) c := by
  classical
  rw [phaseFold_const]
  intro h
  have hval : (AddMonoidAlgebra.single ((0 : Mono), k) c : PhasedSymRing) ((0 : Mono), k)
      = (AddMonoidAlgebra.single ((0 : Mono), (0 : ZMod 4)) c : PhasedSymRing)
        ((0 : Mono), k) := by rw [h]
  rw [Finsupp.single_apply, if_pos rfl, Finsupp.single_apply,
    if_neg (fun hcon => hk (congrArg Prod.snd hcon).symm)] at hval
  exact hc hval

/-! #### Conjugation is the phase-negating ring involution -/

/-- `Conjugate for Term`'s action on a key: `q.phase = (4 − q.phase) % 4`, i.e.
`iᵏ ↦ i^{−k}`. The `f64` coefficients are real and untouched. -/
def negPhase (mk : PMono) : PMono := (mk.1, -mk.2)

noncomputable def conjSymHom : Multiplicative PMono →* PhasedSymRing where
  toFun mk := AddMonoidAlgebra.single (negPhase mk.toAdd) 1
  map_one' := by
    change AddMonoidAlgebra.single (negPhase (0 : PMono)) (1 : ℝ) = 1
    rw [show negPhase (0 : PMono) = (0 : PMono) from Prod.ext rfl neg_zero,
      AddMonoidAlgebra.one_def]
  map_mul' a b := by
    change AddMonoidAlgebra.single (negPhase (a.toAdd + b.toAdd)) (1 : ℝ)
      = AddMonoidAlgebra.single (negPhase a.toAdd) 1
        * AddMonoidAlgebra.single (negPhase b.toAdd) 1
    rw [AddMonoidAlgebra.single_mul_single, mul_one,
      show negPhase (a.toAdd + b.toAdd) = negPhase a.toAdd + negPhase b.toAdd from
        Prod.ext rfl (neg_add _ _)]

/-- **`Conjugate for Term`** (`coeff.rs:257-291`) as a ring map: negate every
monomial's phase exponent. -/
noncomputable def conjSym : PhasedSymRing →ₐ[ℝ] PhasedSymRing :=
  AddMonoidAlgebra.lift ℝ PhasedSymRing PMono conjSymHom

@[simp] theorem conjSym_single (mk : PMono) (c : ℝ) :
    conjSym (AddMonoidAlgebra.single mk c) = AddMonoidAlgebra.single (negPhase mk) c := by
  rw [conjSym, AddMonoidAlgebra.lift_single]
  change c • AddMonoidAlgebra.single (negPhase mk) (1 : ℝ) = _
  rw [AddMonoidAlgebra.smul_single', mul_one]

/-- **Conjugation is an involution** — `conj ∘ conj = id`, so it is a genuine
ring involution on the implemented ring (being an `AlgHom` already gives
additivity and multiplicativity). -/
theorem conjSym_conjSym (x : PhasedSymRing) : conjSym (conjSym x) = x := by
  induction x using AddMonoidAlgebra.induction_linear with
  | zero => simp
  | add f g hf hg => rw [map_add, map_add, hf, hg]
  | single mk c => simp [negPhase]

/-- `star (iᵏ) = i^{−k}` on `ℂ` — the scalar fact behind `evalC_conjSym`, and the
`ℤ/4` form of `Pauli/Matrix.lean`'s `star_iU`. -/
theorem star_iPowI (k : ZMod 4) : star (iPow Complex.I k) = iPow Complex.I (-k) := by
  have h : k = 0 ∨ k = 1 ∨ k = 2 ∨ k = 3 := by revert k; decide
  rcases h with rfl | rfl | rfl | rfl
  · simp
  · rw [show -(1 : ZMod 4) = 3 from by decide, iPowI_one, iPowI_three, Complex.star_def,
      Complex.conj_I]
  · rw [show -(2 : ZMod 4) = 2 from by decide, iPowI_two]
    simp
  · rw [show -(3 : ZMod 4) = 1 from by decide, iPowI_three, iPowI_one, star_neg,
      Complex.star_def, Complex.conj_I, neg_neg]

/-- **`evalC ∘ conj = star ∘ evalC`.** Conjugation on the implemented ring is
compatible with complex conjugation on its image — the property a sesquilinear
pairing needs, and the only thing that makes `Conjugate for Term` (a brand-new
impl with no old counterpart) the *right* map rather than an arbitrary one.
Specializing at `iSym` gives `conj i = −i`. -/
theorem evalC_conjSym (θ : ℕ → ℝ) (x : PhasedSymRing) :
    evalC θ (conjSym x) = star (evalC θ x) := by
  induction x using AddMonoidAlgebra.induction_linear with
  | zero => simp
  | add f g hf hg => simp only [map_add, star_add, hf, hg]
  | single mk c =>
    rw [conjSym_single, evalC_single, evalC_single]
    simp only [phasedValue, negPhase, star_mul', star_iPowI, Complex.star_def,
      Complex.conj_ofReal]

/-- `conj i = −i`, in the value domain — `Pauli/Matrix.lean`'s `star_iU`
transported to the coefficient ring. -/
theorem evalC_conjSym_iSym (θ : ℕ → ℝ) : evalC θ (conjSym iSym) = -Complex.I := by
  rw [evalC_conjSym, evalC_iSym, Complex.star_def, Complex.conj_I]

/-- …and, as with `i·i`, only denotationally: `conj i` is the key `phase 3` with
coefficient `+1`, whereas `−i` is the key `phase 1` with coefficient `−1`. -/
theorem conjSym_iSym_ne_neg_iSym : conjSym iSym ≠ -iSym := by
  classical
  intro h
  rw [iSym, conjSym_single, ← AddMonoidAlgebra.single_neg] at h
  have hval : (AddMonoidAlgebra.single (negPhase ((0 : Mono), (1 : ZMod 4))) (1 : ℝ)
      : PhasedSymRing) ((0 : Mono), (1 : ZMod 4))
      = (AddMonoidAlgebra.single ((0 : Mono), (1 : ZMod 4)) (-1 : ℝ) : PhasedSymRing)
      ((0 : Mono), (1 : ZMod 4)) := by rw [h]
  rw [Finsupp.single_apply, Finsupp.single_apply, if_neg (by
      intro hc
      exact absurd (congrArg Prod.snd hc)
        (show ¬((3 : ZMod 4) = (1 : ZMod 4)) by decide))] at hval
  norm_num at hval
end Phased

/-! ### The `min_eps` arm of `mul_term`'s `clear()` shortcut

`mul.rs:60-73` cites `mulMono_clear_sound` / `mulMono_retain_clear` for the
*whole* `clear()` shortcut, but that pair proves only the **degree** arm, and the
two arms are not the same kind of statement. The degree arm is an equality: every
product monomial provably lands in the truncation ideal, so clearing loses
nothing. The `min_eps` arm is not — it discards monomials the per-monomial rule
in `Sum::add_term` (`|A m · c| < min_eps`) would keep, because a large stored
coefficient can rescue a small multiplier.

`epsClear_ne_retain_pointwise` is that counterexample, and it is what stops a
future reader from "simplifying" the shortcut into a per-monomial loop and
silently changing results. `epsClear_l1_eq` / `epsClear_l1_lt` /
`epsClear_error_lt` are what license keeping it anyway: the shortcut is an
`ℓ¹`-controlled *over*-truncation, not an unsound one. -/

section EpsClear

open PPVM.GradedMap

/-- The batch `Sum::mul_term` produces when multiplying the table by `coeff · p`:
one term `(m·p, A m · c)` per stored monomial. (`mulMonoBatch` is the `c = 1`
case; the scalar is what the `min_eps` arm tests.) -/
noncomputable def mulScaledBatch (p : Mono) (c : ℝ) (A : CMap Mono ℝ) : Multiset (Mono × ℝ) :=
  A.support.val.map fun m => (m + p, A m * c)

/-- `Sum::add_term`'s per-monomial `min_eps` keep-rule applied to that batch —
what the shortcut stands in for. -/
noncomputable def epsRetained (eps : ℝ) (p : Mono) (c : ℝ) (A : CMap Mono ℝ) : CMap Mono ℝ :=
  batchMap ((mulScaledBatch p c A).filter fun t => eps ≤ |t.2|)

/-- **The whole-sum `clear()` is strictly coarser than the per-monomial `min_eps`
rule.** With one stored entry of magnitude `10⁶`, multiplier `c = 10⁻¹³` and
`min_eps = 10⁻¹²`, the shortcut fires (`|c| < eps`) and yields the zero map, while
the per-monomial rule retains the product monomial (`|10⁶ · 10⁻¹³| = 10⁻⁷ ≥ eps`).

So — unlike the degree arm (`mulMono_clear_sound`, an equality) — this arm is an
**over-truncation**. It is sound only in the `ℓ¹` sense proved below, and
replacing it by a per-monomial loop would change results. -/
theorem epsClear_ne_retain_pointwise :
    ∃ (A : CMap Mono ℝ) (p : Mono) (c eps : ℝ), |c| < eps ∧ epsRetained eps p c A ≠ 0 := by
  classical
  refine ⟨Finsupp.single 0 (10 ^ 6 : ℝ), 0, 1 / 10 ^ 13, 1 / 10 ^ 12, by norm_num [abs_of_pos], ?_⟩
  have hs : (Finsupp.single (0 : Mono) (10 ^ 6 : ℝ)).support = {0} :=
    Finsupp.support_single _ (by norm_num)
  have hb : mulScaledBatch 0 (1 / 10 ^ 13) (Finsupp.single (0 : Mono) (10 ^ 6 : ℝ))
      = {((0 : Mono), (10 ^ 6 : ℝ) * (1 / 10 ^ 13))} := by
    rw [mulScaledBatch, hs]
    simp
  rw [epsRetained, hb, Multiset.filter_singleton, if_pos (by norm_num [abs_of_pos])]
  intro h
  have hval : batchMap ({((0 : Mono), (10 ^ 6 : ℝ) * (1 / 10 ^ 13))} : Multiset (Mono × ℝ))
      (0 : Mono) = 0 := by rw [h]; rfl
  simp [batchMap] at hval

/-- The `ℓ¹` mass of a coefficient table. -/
noncomputable def l1 (A : CMap Mono ℝ) : ℝ := ∑ m ∈ A.support, |A m|

/-- **The discarded mass is exactly `|c| · ℓ¹(A)`.** Multiplying by the single
monomial `c·p` scales every coefficient by `c`, and the destination keys `m·p`
are distinct, so no cancellation can occur. -/
theorem epsClear_l1_eq (A : CMap Mono ℝ) (c : ℝ) :
    ∑ m ∈ A.support, |A m * c| = |c| * l1 A := by
  simp_rw [abs_mul, l1, Finset.mul_sum]
  exact Finset.sum_congr rfl fun _ _ => mul_comm _ _

/-- **…so the shortcut is an `ℓ¹`-controlled over-truncation.** When it fires
(`|c| < min_eps`) the mass it throws away is strictly under `min_eps · ℓ¹(A)` —
the sense in which discarding monomials the per-monomial rule would have kept is
still bounded, not arbitrary. -/
theorem epsClear_l1_lt (A : CMap Mono ℝ) (c eps : ℝ) (hc : |c| < eps) (hA : 0 < l1 A) :
    ∑ m ∈ A.support, |A m * c| < eps * l1 A := by
  rw [epsClear_l1_eq, mul_comm _ (l1 A), mul_comm eps]
  exact mul_lt_mul_of_pos_left hc hA

/-- **The observable error of the shortcut.** With per-key expectations bounded by
`1` (`PPVM.Truncation.l1_bound`), clearing the whole table instead of keeping the
monomials the per-monomial rule would keep costs strictly less than
`min_eps · ℓ¹(A)` in the read-out. This is the `ℓ¹` statement
`PPVM/Algebra/Truncation.lean` was built to express, applied to the one
truncation site that is *not* exact. -/
theorem epsClear_error_lt (e : Mono → ℝ) (he : ∀ k, |e k| ≤ 1) (A : CMap Mono ℝ)
    (p : Mono) (c eps : ℝ) (hc : |c| < eps) (hA : 0 < l1 A) :
    |∑ m ∈ A.support, (A m * c) * e (m + p)| < eps * l1 A :=
  lt_of_le_of_lt
    (PPVM.Truncation.l1_bound (fun m => e (m + p)) (fun _ => he _) (fun m => A m * c) A.support)
    (epsClear_l1_lt A c eps hc hA)

end EpsClear

/-! ### `max_sin` is a property of the *representation*, not of the value

Everything above about `max_sin` (`sinDeg_add`, `truncIdeal_mul_right`,
`mulMono_drop_at_insert_eq_drop_at_end`, `mulMono_clear_sound`) is a statement
about the **map-backed** accumulation: `Sum::add_term` and `Sum::mul_term` are the
only two places the bound is ever consulted. The shipped `Term` is a four-way
`Inner`, and three of its arms — `One × One`, `Const × One`, `Const × Const`
(`ppvm-sym-2/src/mul.rs`, the mandatory non-allocating fast forms of the
integration baseline's perf feature 1) — never look at `max_sin` at all.

So the truncated product is **not** a function of the value a `Term` denotes: two
representations of the same polynomial can multiply to different polynomials. That
is `mulImpl_not_wellDefined` below, and it is why
`mulMono_drop_at_insert_eq_drop_at_end` must **not** be cited as an end-to-end
guarantee that the propagated coefficient equals the truncated ring product.
`set_max_sin` is a bound on what the *hash-map* accumulation retains, not a hard
degree bound on the result: a coefficient that stays a single monomial for a whole
circuit escapes it without limit (`mulImpl_one_one_untruncated`,
`fastArm_escapes_bound` — the `sin_pow = 7` at `max_sin = 2` escapee both crates
leave on the Trotter replay, `ppvm-sym-2/src/lib.rs` §"Preserved old quirks").

Stating it as a **negative** result is the point: the "cleanup" that unifies the
four `Inner` arms onto one map-backed representation — exactly the regression perf
feature 1 warns about — would make `mulImpl` well-defined, and would therefore
change numbers. With this theorem in place that shows up as a spec violation
rather than a tidy-up. -/

section FastArms

open PPVM.GradedMap

/-- `ppvm-sym-2`'s `Inner`, cut down to the two arms the fast-path question turns
on: `One` (a single weighted monomial, no allocation) and `Sum` (the map-backed
general case). `Const c` is `One 0 c` and is not modelled separately. -/
inductive Repr where
  /-- `Inner::One(Prod, f64)`. -/
  | one : Mono → ℝ → Repr
  /-- `Inner::Sum(Sum)` — the `FxHashMap<Prod, f64>`. -/
  | sum : CMap Mono ℝ → Repr

/-- The abstraction function: the polynomial a representation denotes. -/
noncomputable def den : Repr → SymRing
  | .one m c => AddMonoidAlgebra.single m c
  | .sum A => A

/-- `Sum::mul_term` — multiply the whole table by the single monomial `c·p`,
dropping every produced monomial over the sine-degree bound at insert. -/
noncomputable def truncMulMono (k : ℕ) (p : Mono) (c : ℝ) (A : CMap Mono ℝ) : CMap Mono ℝ :=
  batchMap ((mulScaledBatch p c A).filter fun t => sinDeg t.1 ≤ k)

/-- **The product `ppvm-sym-2` actually implements**, at bound `k = max_sin`.
Note the asymmetry, which is the shipped code's asymmetry: the `One × One` arm
(`mul.rs`, `Inner::One(p.clone() * p2, c * c2)`) consults nothing, while every arm
that touches the map routes through `Sum::add_term`/`Sum::mul_term` and therefore
through the bound. -/
noncomputable def mulImpl (k : ℕ) : Repr → Repr → Repr
  | .one m c, .one m' c' => .one (m + m') (c * c')
  | .sum A, .one m c => .sum (truncMulMono k m c A)
  | .one m c, .sum A => .sum (truncMulMono k m c A)
  | .sum A, .sum B => .sum (retain (keepSinDeg k) (multiply A B))

/-- `sinDeg` of a one-variable monomial `sᵢ^a cᵢ^b`. -/
theorem sinDeg_single (i a b : ℕ) : sinDeg (Finsupp.single i (a, b)) = a :=
  Finsupp.sum_single_index rfl

/-- `sin(x₀)²` — a monomial that already sits at the bound `k = 2`. -/
noncomputable def sinSq : Mono := Finsupp.single 0 (2, 0)

/-- `sin(x₁)` — the multiplier that pushes it over. -/
noncomputable def sinOne : Mono := Finsupp.single 1 (1, 0)

theorem sinDeg_sinSq_add_sinOne : sinDeg (sinSq + sinOne) = 3 := by
  rw [sinDeg_add, sinSq, sinOne, sinDeg_single, sinDeg_single]

/-- **The fast arm computes the *untruncated* ring product, for every `k`.** It is
exact — and therefore unbounded: `max_sin` has no effect on it whatsoever. -/
theorem mulImpl_one_one_untruncated (k : ℕ) (m m' : Mono) (c c' : ℝ) :
    den (mulImpl k (.one m c) (.one m' c')) = den (.one m c) * den (.one m' c') :=
  (single_mul_single m m' c c').symm

/-- **…so a single-monomial coefficient escapes `set_max_sin` without limit.** At
`k = 2`, `sin(x₀)² · sin(x₁)` comes out of the fast arm as a live monomial of sine
degree `3`. Users reading `set_max_sin` as a hard degree bound on the propagated
coefficient are wrong; it bounds only what the map-backed accumulation retains. -/
theorem fastArm_escapes_bound :
    2 < sinDeg (sinSq + sinOne) ∧ den (mulImpl 2 (.one sinSq 1) (.one sinOne 1)) ≠ 0 := by
  classical
  refine ⟨by rw [sinDeg_sinSq_add_sinOne]; norm_num, ?_⟩
  change AddMonoidAlgebra.single (sinSq + sinOne) ((1 : ℝ) * 1) ≠ 0
  simp

/-- **The truncated product does not factor through the denotation.** With
`a₁ = One(sin(x₀)², 1)` and `a₂` the map-backed `Sum` denoting the same
polynomial, `b = One(sin(x₁), 1)` and `k = 2`:

* `den a₁ = den a₂` — the two are the same element of `SymRing`;
* `mulImpl 2 a₁ b` takes the fast arm and keeps the degree-`3` monomial;
* `mulImpl 2 a₂ b` routes through `Sum::mul_term`, which drops it, leaving `0`.

Hence `mulImpl k` is not well defined on the quotient by `den` for finite `k`, and
no theorem of the form "the propagated coefficient equals the truncated ring
product" can hold for the shipped type. (`mulImpl_one_one_untruncated` is the
positive half: the offending arm is exact, just unbounded.) -/
theorem mulImpl_not_wellDefined :
    ∃ (a₁ a₂ b : Repr) (k : ℕ),
      den a₁ = den a₂ ∧ den (mulImpl k a₁ b) ≠ den (mulImpl k a₂ b) := by
  classical
  refine ⟨.one sinSq 1, .sum (Finsupp.single sinSq 1), .one sinOne 1, 2, rfl, ?_⟩
  have hs : (Finsupp.single sinSq (1 : ℝ)).support = {sinSq} :=
    Finsupp.support_single _ one_ne_zero
  have hb : mulScaledBatch sinOne 1 (Finsupp.single sinSq (1 : ℝ))
      = {(sinSq + sinOne, (1 : ℝ) * 1)} := by
    rw [mulScaledBatch, hs]
    simp
  have hzero : den (mulImpl 2 (.sum (Finsupp.single sinSq (1 : ℝ))) (.one sinOne 1)) = 0 := by
    change truncMulMono 2 sinOne 1 (Finsupp.single sinSq (1 : ℝ)) = 0
    rw [truncMulMono, hb, Multiset.filter_singleton,
      if_neg (by simp [sinDeg_sinSq_add_sinOne])]
    simp [batchMap]
  rw [hzero]
  exact fastArm_escapes_bound.2

end FastArms

end PPVM.Symbolic
