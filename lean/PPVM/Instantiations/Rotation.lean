/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Pauli.Phase
import PPVM.Algebra.GradedMap
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-!
# Non-Clifford rotation as a branch on `C[K]`

The design's non-Clifford producer
(`traits-2-configuration-and-hashing.md`, "Every gate is a producer"): a
rotation `extend[s]-by-linearity` a single term into a small fan-out,

  `e^{-iθG/2} · (c · P) · e^{iθG/2} = c·cos θ · P + c·sin θ · (iGP)`

when `G` and `P` anticommute (`{G,P}=0`), and leaves `c·P` unchanged when they
commute. This file validates two things about that branch — **at the level of
the bits and the coefficient plane**, not by deriving the operator conjugation:

* **The new key `P' = iGP` is a distinct Pauli** whose symplectic bits are
  `G ⊕ P` — the same `⊕` key product the L4 algebra uses (`PPVM.PauliPhase`) —
  so the branch produces exactly one extra term, distinct from both `G` and `P`.
  This part *is* derived (from `mulBits` and the anticommutation `ω = 1`).
* **The coefficient update `(c_P, c_{P'})` is a 2-D rotation.** We formalize the
  standard rotation on `ℝ²` and prove it is norm-preserving, reversible, and
  angle-additive (all from `sin²θ + cos²θ = 1`). This is the abstract shape the
  branch applies to the coefficient pair; the operator identity
  `e^{-iθG/2}(cP)e^{iθG/2} = …` is *modeled* by this 2-D rotation, **not** derived
  from operator algebra here. It gives the algebraic guarantees the branch-free
  `eps` / `comm_2` / `SIGN_NEG` kernel (`sum/rot1.rs`, `sum/rot2.rs`) inherits.
-/

namespace PPVM.Rotation

open PPVM.PauliPhase

/-! ### The new key `iGP` is a distinct Pauli with bits `G ⊕ P` -/

/-- **The anticommuting branch produces a genuinely new key, distinct from both
operands.** If `G = (a,b)` and `P = (c,d)` anticommute (`ω = 1`), then `iGP`,
whose symplectic bits are `G ⊕ P` (`mulBits`), differs from `P` *and* from `G`.
So a rotation adds exactly one fresh term per input; it never silently merges
back into `P` or into `G`. -/
theorem anticommute_new_key :
    ∀ a b c d, omega a b c d = 1 →
      mulBits a b c d ≠ (c, d) ∧ mulBits a b c d ≠ (a, b) := by decide

/-- The product bits are `(a⊕c, b⊕d)` — a definitional restatement of `mulBits`.
(In the *commuting* case the rotation is inert not because of these bits but
because the `sin` branch carries coefficient `0`; `rot 0` below is the identity.) -/
theorem commute_bits (a b c d : Bool) :
    mulBits a b c d = (xor a c, xor b d) := rfl

/-! ### The branch sign `ε` is the real phase of `iGP`, derived from `phaseExp`

`Rotation.lean`'s 2-D `rot` *models* the coefficient plane but does not derive the
per-axis `±1` sign the branch coefficient carries. That sign is a genuine Pauli
phase fact: the stored branch key is the **real** word `g(G⊕P)`, and the physical
branch term is `sinθ · (iGP)`, with

  `iGP = i · (G·P) = i · i^{phaseExp(G,P)} · g(G⊕P) = i^{1 + phaseExp(G,P)} · g(G⊕P)`.

So the coefficient sign is `ε = i^{1 + phaseExp(G,P)}` as a base-`i` exponent
(`branchExp` below). When `G` and `P` anticommute (`ω = 1`) the single-qubit
product `G·P` carries an odd power of `i`, so `1 + phaseExp` is *even* and `ε` is a
real `±1` — this is why the leading `i` "cancels" and the branch stays real. The
three axis theorems then check that this derived `±1` equals the hand-ported table
in `producer.rs:141-143` (`RotationProducer::produce`) case-by-case. -/

/-- The base-`i` exponent of the prefactor sitting on the *real* branch word
`g(G⊕P)`: `iGP = i^{1 + phaseExp(G,P)} · g(G⊕P)`. -/
def branchExp (gx gz x z : Bool) : ZMod 4 := 1 + phaseExp gx gz x z

/-- **The `i` in `iGP` cancels the anticommuting product's `i`.** When `G = (gx,gz)`
and `P = (x,z)` anticommute (`ω = 1`), `branchExp` is even (`∈ {0,2}`), i.e. the
branch prefactor is a real `±1`, never `±i` — so `producer.rs` may store a real
coefficient sign `ε` on the real word `g(G⊕P)`. -/
theorem branchExp_isRealPhase (gx gz x z : Bool) (h : omega gx gz x z = 1) :
    IsRealPhase (branchExp gx gz x z) := by
  change branchExp gx gz x z = 0 ∨ branchExp gx gz x z = 2
  revert h; revert gx gz x z; decide

/-- **`rx` (`G = X = (1,0)`): the table `ε = −1 iff x`.** For an anticommuting `P`
(`ω = 1`, i.e. `z = 1`), the phase-derived branch sign `i^{1+phaseExp(X,P)}` equals
`−1` (exponent `2`) exactly when `x`, matching `producer.rs:141`
(`RotAxis::X => … eps = if x { -1 } else { 1 }`). -/
theorem rx_eps_from_product :
    ∀ x z, omega true false x z = 1 →
      branchExp true false x z = (if x then 2 else 0) := by decide

/-- **`ry` (`G = Y = (1,1)`): the table `ε = −1 iff z`.** For an anticommuting `P`
(`ω = 1`, i.e. `x ≠ z`), `i^{1+phaseExp(Y,P)}` equals `−1` exactly when `z`,
matching `producer.rs:142` (`RotAxis::Y => … eps = if z { -1 } else { 1 }`). -/
theorem ry_eps_from_product :
    ∀ x z, omega true true x z = 1 →
      branchExp true true x z = (if z then 2 else 0) := by decide

/-- **`rz` (`G = Z = (0,1)`): the table `ε = +1 iff z`.** For an anticommuting `P`
(`ω = 1`, i.e. `x = 1`), `i^{1+phaseExp(Z,P)}` equals `+1` (exponent `0`) exactly
when `z`, matching `producer.rs:143`
(`RotAxis::Z => … eps = if z { 1 } else { -1 }`). -/
theorem rz_eps_from_product :
    ∀ x z, omega false true x z = 1 →
      branchExp false true x z = (if z then 0 else 2) := by decide

/-! ### The two-site branch (`rzz` / `rxx` / `ryy`, and the generic `rotate_2`)

Everything above is stated at **one** site. `RotationTwo` rotates about a
two-site generator `G = G_a ⊗ G_b` acting on sites `(a,b)`, and its three native
kernels (`rzz`/`rxx`/`ryy`, `crates/ppvm-pauli-sum/src/sum/rot2.rs:62-184`, ported
onto `Sum`'s in-place rotate path) read only the **two bit pairs** at `(a,b)` and
consult *no* Levi-Civita table: each toggles bits at both sites at once and gets
its `ε` from a *two-factor* product, so the sign convention has two independent
places to be wrong. The old crate validates these kernels only by a Rust diff
against its own generic `comm_2`/`rotate_2` (`rot2.rs::rxx_matches_generic`) —
agreement between two implementations, not an oracle. This section is the oracle.

The two-site facts are exactly the one-site ones plus one extra ingredient:
`phaseExp_of_commute` — a *commuting* single-site pair contributes **no** phase at
all. So the two-site branch exponent collapses onto the unique anticommuting
site, and the `ε` column of each kernel is the corresponding single-site column
(`rx`/`ry`/`rz` above) evaluated there. -/

/-- The single-site anticommutation bit `ω(G,P)` as a `Bool` — the two-bit
predicate the kernels literally evaluate (`x_G ∧ z_P` xor `z_G ∧ x_P`). -/
def antiBit (gx gz x z : Bool) : Bool := xor (gx && z) (gz && x)

/-- `antiBit` is the symplectic form `ω` of `PPVM.PauliPhase`. -/
theorem antiBit_eq_omega : ∀ gx gz x z,
    (antiBit gx gz x z = true ↔ omega gx gz x z = 1) := by decide

/-- **The two-site branch key `w₂(P,G) = ω(P_a,G_a) ⊕ ω(P_b,G_b)`.** The branch
fires iff *exactly one* site anticommutes; both-commute and both-anticommute are
inert. This is the single boolean each kernel's early-`return None` tests. -/
def anti2 (gxa gza gxb gzb xa za xb zb : Bool) : Bool :=
  xor (antiBit gxa gza xa za) (antiBit gxb gzb xb zb)

/-- **A commuting single-site pair carries no phase at all.** `G·P = g(G⊕P)`
exactly (not merely up to a real sign) whenever `ω(G,P) = 0` — the reason the
two-site branch exponent below collapses onto the one anticommuting site. -/
theorem phaseExp_of_commute : ∀ gx gz x z,
    antiBit gx gz x z = false → phaseExp gx gz x z = 0 := by decide

/-- The base-`i` exponent of the prefactor sitting on the **real** two-site branch
word `g(G ⊕ P)`: `iGP = i^{1 + phaseExp_a + phaseExp_b} · g(G⊕P)`, because the two
per-site cocycles simply add (`PPVM.PauliPhase.TwoPauli.mul`). -/
def branchExp2 (gxa gza gxb gzb xa za xb zb : Bool) : ZMod 4 :=
  1 + phaseExp gxa gza xa za + phaseExp gxb gzb xb zb

/-- **The leading `i` still cancels at two sites.** When exactly one site
anticommutes (`w₂ = 1`) the two-site product `G·P` carries a single odd power of
`i`, so `branchExp2` is even and the branch prefactor is a real `±1`, never `±i`
— the two-site form of `branchExp_isRealPhase`, and the licence for the kernels
to store a real `ε` on the real word `g(G⊕P)`. -/
theorem branchExp2_isRealPhase (gxa gza gxb gzb xa za xb zb : Bool)
    (h : anti2 gxa gza gxb gzb xa za xb zb = true) :
    IsRealPhase (branchExp2 gxa gza gxb gzb xa za xb zb) := by
  change branchExp2 gxa gza gxb gzb xa za xb zb = 0
    ∨ branchExp2 gxa gza gxb gzb xa za xb zb = 2
  revert h; revert gxa gza gxb gzb xa za xb zb; decide

/-- The two-site branch key: `mulBits` at **both** sites (`G ⊕ P` site-wise). -/
def mulBits2 (gxa gza gxb gzb xa za xb zb : Bool) : (Bool × Bool) × (Bool × Bool) :=
  (mulBits gxa gza xa za, mulBits gxb gzb xb zb)

/-- **The two-site branch also produces a genuinely new key.** When `w₂ = 1` the
branch word differs from `P` *and* from `G`, so a two-qubit rotation adds exactly
one fresh term per input and never silently merges back into its own source. (The
two-pass structure of `RotateInPlace` is still what handles collisions with
*other* keys — see `accumulate_rotBatch` / `eagerWalk_ne_twoPass` below.) -/
theorem anticommute_new_key2 : ∀ gxa gza gxb gzb xa za xb zb,
    anti2 gxa gza gxb gzb xa za xb zb = true →
      mulBits2 gxa gza gxb gzb xa za xb zb ≠ ((xa, za), (xb, zb))
        ∧ mulBits2 gxa gza gxb gzb xa za xb zb ≠ ((gxa, gza), (gxb, gzb)) := by decide

/-! #### The three native kernels

`rzz`, `rxx`, `ryy` are the `G = Z⊗Z`, `X⊗X`, `Y⊗Y` instances. For each we pin
(i) the branch predicate, (ii) the bits the kernel toggles, and (iii) the `ε`
column — each against the exact two-bit expression in `sum/rot2.rs`. -/

/-- **`rzz` branch predicate**: `ZZ` commutes iff the two sites agree on carrying
an `X` component (`rot2.rs:139`, `if xa == xb { return None }`). -/
theorem rzz_anti : ∀ xa za xb zb,
    anti2 false true false true xa za xb zb = xor xa xb := by decide

/-- **`rzz` branch key**: the `x` bits are untouched and both `z` bits flip
(`set_zbit(a, !za); set_zbit(b, !zb)`). -/
theorem rzz_branch_key : ∀ xa za xb zb,
    mulBits2 false true false true xa za xb zb = ((xa, !za), (xb, !zb)) := by decide

/-- **`rzz`'s `ε` column**: `ε = +1` iff the *anticommuting* site (the one with
`x` set) carries `z` too, i.e. is `Y` rather than `X` — exactly
`let z_anti = if xa { za } else { zb }; eps = if z_anti { 1 } else { -1 }`
(`rot2.rs:145-148`). Derived, not asserted: it is `i^{1+phaseExp_a+phaseExp_b}`. -/
theorem rzz_eps_from_product : ∀ xa za xb zb, xor xa xb = true →
    branchExp2 false true false true xa za xb zb
      = (if (if xa then za else zb) then 0 else 2) := by decide

/-- **`rxx` branch predicate**: `XX` commutes iff the two sites agree on carrying
a `Z` component (`rot2.rs:177`, `if za == zb { return None }`). -/
theorem rxx_anti : ∀ xa za xb zb,
    anti2 true false true false xa za xb zb = xor za zb := by decide

/-- **`rxx` branch key**: the `z` bits are untouched and both `x` bits flip. -/
theorem rxx_branch_key : ∀ xa za xb zb,
    mulBits2 true false true false xa za xb zb = ((!xa, za), (!xb, zb)) := by decide

/-- **`rxx`'s `ε` column**: `ε = −1` iff the anticommuting site (the one with `z`
set) carries `x` too, i.e. is `Y` rather than `Z` — exactly
`let x_anti = if za { xa } else { xb }; eps = if x_anti { -1 } else { 1 }`. -/
theorem rxx_eps_from_product : ∀ xa za xb zb, xor za zb = true →
    branchExp2 true false true false xa za xb zb
      = (if (if za then xa else xb) then 2 else 0) := by decide

/-- **`ryy` branch predicate**: `YY` commutes iff the two sites agree on `x ⊕ z`
(`rot2.rs:222`, `pa = xa ^ za`, `pb = xb ^ zb`, `if pa == pb { return None }`). -/
theorem ryy_anti : ∀ xa za xb zb,
    anti2 true true true true xa za xb zb = xor (xor xa za) (xor xb zb) := by decide

/-- **`ryy` branch key**: both bits flip at both sites. -/
theorem ryy_branch_key : ∀ xa za xb zb,
    mulBits2 true true true true xa za xb zb = ((!xa, !za), (!xb, !zb)) := by decide

/-- **`ryy`'s `ε` column**: `ε = +1` iff the anticommuting site (the one with
`x ≠ z`) is `X` rather than `Z` — exactly
`let x_anti = if pa { xa } else { xb }; eps = if x_anti { 1 } else { -1 }`. -/
theorem ryy_eps_from_product : ∀ xa za xb zb, xor (xor xa za) (xor xb zb) = true →
    branchExp2 true true true true xa za xb zb
      = (if (if xor xa za then xa else xb) then 0 else 2) := by decide

/-! ### The branch is a norm-preserving 2-D rotation

Model the `(coefficient of P, coefficient of P')` pair as a point of `ℝ²`. The
rotation sends `P ↦ cos θ · P + sin θ · P'` and `P' ↦ −sin θ · P + cos θ · P'`,
i.e. the standard rotation matrix. -/

/-- The rotation acting on the `(P, P')` coefficient plane. Starting from a pure
`P` (`v = (1,0)`) gives `(cos θ, sin θ)` — the design's `cos·w + sin·w'`. -/
noncomputable def rot (θ : ℝ) (v : ℝ × ℝ) : ℝ × ℝ :=
  (Real.cos θ * v.1 - Real.sin θ * v.2, Real.sin θ * v.1 + Real.cos θ * v.2)

@[simp] theorem rot_fst (θ : ℝ) (v : ℝ × ℝ) :
    (rot θ v).1 = Real.cos θ * v.1 - Real.sin θ * v.2 := rfl

@[simp] theorem rot_snd (θ : ℝ) (v : ℝ × ℝ) :
    (rot θ v).2 = Real.sin θ * v.1 + Real.cos θ * v.2 := rfl

/-- A pure `P` rotates to `cos θ · P + sin θ · P'` — the exact branch shape of the
design's rotation producer. -/
theorem rot_basis (θ : ℝ) : rot θ (1, 0) = (Real.cos θ, Real.sin θ) := by
  simp [rot]

/-- `θ = 0` is the identity (the commuting / no-op branch). -/
theorem rot_zero (v : ℝ × ℝ) : rot 0 v = v := by
  simp [rot]

/-- **Reversibility.** `R_{-θ} ∘ R_θ = id`: rotating by `θ` and then `−θ` returns
the original coefficients. This is the algebraic core of the gate being
invertible, and it uses exactly `sin²θ + cos²θ = 1`. -/
theorem rot_neg_rot (θ : ℝ) (v : ℝ × ℝ) : rot (-θ) (rot θ v) = v := by
  have h := Real.sin_sq_add_cos_sq θ
  refine Prod.ext ?_ ?_ <;>
    simp only [rot_fst, rot_snd, Real.cos_neg, Real.sin_neg]
  · linear_combination v.1 * h
  · linear_combination v.2 * h

/-- **Norm preservation.** The 2-D rotation preserves the `ℓ²` norm of the
coefficient pair, `‖cos·P + sin·P'‖² = ‖P‖²`, from `sin²θ + cos²θ = 1`. This is
the algebraic reason the branch is norm-preserving; the operator-level unitarity
of `e^{-iθG/2}` is *modeled* by this, not derived from operator algebra. -/
theorem rot_norm_sq (θ : ℝ) (v : ℝ × ℝ) :
    (rot θ v).1 ^ 2 + (rot θ v).2 ^ 2 = v.1 ^ 2 + v.2 ^ 2 := by
  have h := Real.sin_sq_add_cos_sq θ
  simp only [rot_fst, rot_snd]
  linear_combination (v.1 ^ 2 + v.2 ^ 2) * h

/-- **Angle addition.** Composing two rotations adds the angles, so a Trotter
step of many small rotations about the same axis is one rotation — the identity
a rotation-merging optimization relies on. -/
theorem rot_rot (θ φ : ℝ) (v : ℝ × ℝ) : rot θ (rot φ v) = rot (θ + φ) v := by
  refine Prod.ext ?_ ?_ <;>
    simp only [rot_fst, rot_snd, Real.cos_add, Real.sin_add] <;> ring

/-! ### The whole-map two-pass decomposition (`RotateInPlace`)

Everything above is a *per-term* fact: one input key, one branch key, one
coefficient pair. The thing every rotation in every workload actually runs is
`RotateInPlace` (`crates/ppvm-pauli-sum-2/src/store.rs`), a **two-pass** fused
walk over the whole map:

* **pass 1** — `iter_mut` over the live support: scale each diagonal coefficient
  in place (`cos θ`, key untouched, cached hash intact) and *buffer* the branch
  term `(iGP, c·sinθ·ε)` in `scratch`;
* **pass 2** — merge the buffered branch terms into the map with `add_term`,
  aggregating collisions.

No per-term theorem licenses this, because the risk is **cross-term**: a branch
key produced from `k` can collide with a *different* key `k'` that pass 1 has not
scaled yet (`rx` on a support containing both `Z` and `Y` at one site sends
`Z ↦ Y` and `Y ↦ Z` simultaneously). `anticommute_new_key` is the near miss — it
rules out a branch colliding with *its own* source, which is a strictly weaker
statement. The correctness of the split rests on pass 1 scaling **all** diagonals
before pass 2 merges **any** branch.

`accumulate_rotBatch` below is the licence: the two-pass map equals the one-pass
"produce all `2N` terms, then accumulate the batch" semantics of the design's
generic producer/`Accumulate` path — and, since the batch is a `Multiset`, for
*every* order the walk visits keys in (so a hash-partitioned or columnar backend
may reorder freely, `GradedMap.accumulateTerms_perm`/`_add`).

`eagerWalk_ne_twoPass` shows the theorem has content: an implementation that
merged each branch *eagerly inside* the walk — the natural "tidier" single-pass
refactor, and what a backend interleaving the two passes would compute — is
observably different, because a later diagonal then gets scaled *after* an
earlier branch has already been added to it. -/

section TwoPass

open PPVM.GradedMap

variable {K C : Type*} [CommRing C]

/-- The terms one input `(k, a)` produces: on an anticommuting key the diagonal
`(k, d·a)` (`d = cos θ`) **and** the branch `(br k, s k · a)`
(`s k = sinθ·ε(G,q,k)`, `br k = iGP`'s key); on a commuting key just `(k, a)`.
This is `RotationProducer::produce`'s fan-out, as an unordered batch. -/
def rotTerms (anti : K → Prop) [DecidablePred anti] (br : K → K) (d : C) (s : K → C)
    (k : K) (a : C) : Multiset (K × C) :=
  if anti k then {(k, d * a), (br k, s k * a)} else {(k, a)}

/-- The whole rotation as **one** batch of `≤ 2N` produced terms — the design's
generic `TermProducer` → `TermSink` → `accumulate_batch` path, with no fast path.
The `Multiset` records that the walk order is not part of the semantics. -/
noncomputable def rotBatch (anti : K → Prop) [DecidablePred anti] (br : K → K) (d : C)
    (s : K → C) (A : CMap K C) : Multiset (K × C) :=
  A.support.val.bind fun k => rotTerms anti br d s k (A k)

/-- **Pass 1** — scale every diagonal coefficient in place. No key moves, so this
is a pointwise rescale of the *whole* map, done before any branch is merged. -/
noncomputable def diagPass (anti : K → Prop) [DecidablePred anti] (d : C)
    (A : CMap K C) : CMap K C :=
  A.sum fun k a => Finsupp.single k (if anti k then d * a else a)

/-- **Pass 2** — merge every buffered branch term, aggregating collisions. -/
noncomputable def branchPass (anti : K → Prop) [DecidablePred anti] (br : K → K)
    (s : K → C) (A : CMap K C) : CMap K C :=
  A.sum fun k a => if anti k then Finsupp.single (br k) (s k * a) else 0

/-- The `RotateInPlace` fast path: pass 1 then pass 2. -/
noncomputable def twoPass (anti : K → Prop) [DecidablePred anti] (br : K → K) (d : C)
    (s : K → C) (A : CMap K C) : CMap K C :=
  diagPass anti d A + branchPass anti br s A

/-- Pass 1 really is the pointwise rescale `k ↦ if anticommutes then cosθ·A k
else A k` — in particular it touches *every* key of the support, including the
ones that are also branch destinations. -/
theorem diagPass_apply (anti : K → Prop) [DecidablePred anti] (d : C) (A : CMap K C)
    (j : K) : diagPass anti d A j = if anti j then d * A j else A j := by
  classical
  rw [diagPass, Finsupp.sum, Finset.sum_apply']
  have hterm : ∀ k ∈ A.support,
      (Finsupp.single k (if anti k then d * A k else A k) : CMap K C) j
        = if k = j then (if anti j then d * A j else A j) else 0 := by
    intro k _
    rw [Finsupp.single_apply]
    by_cases hk : k = j
    · subst hk; rfl
    · simp only [if_neg hk]
  rw [Finset.sum_congr rfl hterm, Finset.sum_ite_eq' A.support j]
  by_cases hj : j ∈ A.support
  · rw [if_pos hj]
  · rw [if_neg hj, Finsupp.notMem_support_iff.mp hj, mul_zero, ite_self]

/-- **The two-pass fast path computes the one-pass batch semantics.** Folding the
whole `≤ 2N`-term batch into the empty map (the generic
producer/`accumulate_batch` path) equals scaling every diagonal in place *first*
and merging all branch terms *second*.

This is the correctness licence for `RotateInPlace`, and it holds for **every**
order the walk visits keys in: the batch is a `Multiset`, so
`GradedMap.accumulateTerms_perm` (reordering) and `accumulateTerms_add`
(partitioning across a parallel/columnar backend) apply to the left-hand side
unchanged. -/
theorem accumulate_rotBatch (anti : K → Prop) [DecidablePred anti] (br : K → K) (d : C)
    (s : K → C) (A : CMap K C) :
    accumulateTerms (rotBatch anti br d s A) 0 = twoPass anti br d s A := by
  classical
  rw [accumulateTerms_eq, accumulate, zero_add]
  have hbatch : ∀ k : K, batchMap (rotTerms anti br d s k (A k))
      = Finsupp.single k (if anti k then d * A k else A k)
        + (if anti k then Finsupp.single (br k) (s k * A k) else 0) := by
    intro k
    by_cases hk : anti k <;> simp [batchMap, rotTerms, hk]
  simp only [rotBatch, batchMap, Multiset.map_bind, Multiset.sum_bind, twoPass, diagPass,
    branchPass, Finsupp.sum, Finset.sum, ← Multiset.sum_map_add]
  exact congrArg Multiset.sum (Multiset.map_congr rfl fun k _ => hbatch k)

/-! #### Interleaving the two passes is observably wrong

An `eagerStep` visits one key, scales the coefficient the map *currently* holds
there, and merges that key's branch immediately — the single-pass shape a
"tidier" refactor (or a backend that fuses the produce and the merge) would
compute. -/

/-- One eager step: rescale the coefficient at `k` in place (`m k ↦ d · m k`, via
the `+ (d−1)·m k` correction) and merge `k`'s branch, both against the *current*
map. -/
noncomputable def eagerStep (anti : K → Prop) [DecidablePred anti] (br : K → K) (d : C)
    (s : K → C) (m : CMap K C) (k : K) : CMap K C :=
  if anti k then
    m + Finsupp.single k ((d - 1) * m k) + Finsupp.single (br k) (s k * m k)
  else m

/-- The eager walk in a given visit order. -/
noncomputable def eagerWalk (anti : K → Prop) [DecidablePred anti] (br : K → K) (d : C)
    (s : K → C) (l : List K) (m : CMap K C) : CMap K C :=
  l.foldl (eagerStep anti br d s) m

/-- **Merging branches eagerly inside the walk is wrong.** Two mutually-branching
keys (the `rx`-on-`{Z, Y}` shape: `br` swaps them, both anticommute) already
separate the two implementations: with `d = 1`, `s ≡ 1` and coefficients
`(1, 1)`, the two-pass map carries `2` at the first key while the eager walk
carries `3`, because visiting the second key rescales *and* re-branches a
coefficient that already absorbed the first key's branch.

So the two-pass structure of `RotateInPlace` is load-bearing, not a mere
optimization: `anticommute_new_key` (a branch never collides with its own source)
still holds here and does not save the eager version. -/
theorem eagerWalk_ne_twoPass :
    eagerWalk (fun _ : Bool => True) Bool.not (1 : ℤ) (fun _ => 1) [true, false]
        (Finsupp.single true 1 + Finsupp.single false 1)
      ≠ twoPass (fun _ : Bool => True) Bool.not (1 : ℤ) (fun _ => 1)
        (Finsupp.single true 1 + Finsupp.single false 1) := by
  classical
  intro h
  have hval := DFunLike.congr_fun h true
  have hbranch : branchPass (fun _ : Bool => True) Bool.not (fun _ => (1 : ℤ))
      (Finsupp.single true 1 + Finsupp.single false 1)
      = Finsupp.single false 1 + Finsupp.single true 1 := by
    rw [branchPass, Finsupp.sum_add_index' (fun _ => by simp)
      (fun _ _ _ => by simp [mul_add, Finsupp.single_add]),
      Finsupp.sum_single_index (by simp), Finsupp.sum_single_index (by simp)]
    simp
  rw [eagerWalk, twoPass, hbranch] at hval
  simp only [List.foldl, eagerStep, Bool.not_true, Bool.not_false, if_pos trivial,
    diagPass_apply, Finsupp.add_apply, Finsupp.single_apply] at hval
  norm_num at hval

/-- **The fused two-pass `RotateInPlace` shape is sound for the two-qubit branch
too.** `accumulate_rotBatch` is stated over an *abstract* key type with an
abstract anticommutation predicate and branch map, so the two-site kernels are an
instance of it: read the four bits at `(a,b)` off the key (`bits`), branch on
`anti2`, and re-key by the two-site product. The one-pass produced-batch
semantics and the two-pass walk agree, in every visit order — nothing about the
lift is single-site. -/
theorem accumulate_rotBatch_two (bits : K → Bool × Bool × Bool × Bool)
    (gxa gza gxb gzb : Bool) (br : K → K) (d : C) (s : K → C) (A : CMap K C) :
    accumulateTerms
        (rotBatch (fun k => anti2 gxa gza gxb gzb (bits k).1 (bits k).2.1
          (bits k).2.2.1 (bits k).2.2.2 = true) br d s A) 0
      = twoPass (fun k => anti2 gxa gza gxb gzb (bits k).1 (bits k).2.1
          (bits k).2.2.1 (bits k).2.2.2 = true) br d s A :=
  accumulate_rotBatch _ br d s A

end TwoPass

/-! ### The **generic** `rotate_2` kernel: `comm_2`'s `SIGN_NEG` table

Everything above pins the three *native* two-site kernels (`rzz`/`rxx`/`ryy`),
which is where the workloads spend their time. But `rotate_2` is itself a shipped
public method accepting an **arbitrary** `[x, z]` axis pair, and it does not use
those kernels: it runs the branch-free `comm_2`
(`crates/ppvm-pauli-sum-2/src/rotation.rs:109-140`, ported verbatim from
`ppvm-pauli-sum/src/sum/rot2.rs`) over a hand-rolled 16-entry orientation mask
`SIGN_NEG = 0x2840`, and then applies the *opposite* sign convention to the fast
paths: `sin.mul_sign(-eps)` where `rzz`/`rxx`/`ryy` use `sin.mul_sign(eps)`.

The only checks on any of that — old's `*_matches_generic` tests and this crate's
ports of them — compare the generic path against the fast paths at exactly the
three **diagonal** axes `ZZ`, `XX`, `YY`. Every off-diagonal pair (`XZ`, `YX`,
`ZY`, …) is untested on both sides, and a `+eps`/`−eps` asymmetry is exactly the
kind of thing that is correct by coincidence on the diagonal. This section closes
it: over all `2⁸` (axis, key) bit patterns, the generic path's sign is the real
branch prefactor `i^{1+phaseExp_a+phaseExp_b}` of `branchExp2`, and its masked
output bits are `mulBits2`. -/

section GenericTwo

/-- The `SIGN_NEG` constant of `comm_2`, verbatim (`0x2840`). Its set bits are
`{6, 11, 13}`, i.e. the negatively-oriented single-site pairs `(X, Z)`, `(Z, Y)`,
`(Y, X)`. -/
def signNegMask : ℕ := 0x2840

/-- `comm_2`'s per-site table index `(z_G ≪ 3) | (x_G ≪ 2) | (z_P ≪ 1) | x_P`. -/
def signNegIdx (gx gz x z : Bool) : ℕ :=
  (if gz then 8 else 0) + (if gx then 4 else 0) + (if z then 2 else 0) + (if x then 1 else 0)

/-- The mask lookup `((SIGN_NEG >> idx) as u8) & 1`. -/
def signNegBit (gx gz x z : Bool) : Bool :=
  (signNegMask >>> signNegIdx gx gz x z) % 2 == 1

/-- **`comm_2`'s returned `ε`**, transcribed line for line:

```text
a0 = (x_a & z_c) ^ (z_a & x_c);   a1 = (x_b & z_d) ^ (z_b & x_d);
present = a0 ^ a1;
neg0 = ((SIGN_NEG >> idx0) & 1) & a0;   neg1 = ((SIGN_NEG >> idx1) & 1) & a1;
coeff = ((1 - (neg0 << 1)) * a0 + (1 - (neg1 << 1)) * a1) * present;
```
-/
def comm2Coeff (gxa gza gxb gzb xa za xb zb : Bool) : ℤ :=
  let a0 := antiBit gxa gza xa za
  let a1 := antiBit gxb gzb xb zb
  let neg0 := signNegBit gxa gza xa za && a0
  let neg1 := signNegBit gxb gzb xb zb && a1
  ((1 - 2 * (if neg0 then 1 else 0)) * (if a0 then 1 else 0)
      + (1 - 2 * (if neg1 then 1 else 0)) * (if a1 then 1 else 0))
    * (if anti2 gxa gza gxb gzb xa za xb zb then 1 else 0)

/-- `comm_2`'s four returned bit flags `(x_out0, z_out0, x_out1, z_out1)`, each
masked by `present`. -/
def comm2Key (gxa gza gxb gzb xa za xb zb : Bool) : (Bool × Bool) × (Bool × Bool) :=
  let present := anti2 gxa gza gxb gzb xa za xb zb
  ((xor gxa xa && present, xor gza za && present),
   (xor gxb xb && present, xor gzb zb && present))

/-- A real (`±1`) `ℤ/4` phase as an integer sign. -/
def realSign (φ : ZMod 4) : ℤ := if φ = 2 then -1 else 1

/-- **`comm_2`'s `ε` is zero exactly on the commuting pairs**, so the generic
kernel's `if eps == 0 { return None }` early-out is precisely the `anti2` branch
predicate the native kernels test — the two paths agree on *which* terms branch,
at every axis pair, not only the diagonal ones. -/
theorem comm2Coeff_eq_zero_iff : ∀ gxa gza gxb gzb xa za xb zb,
    (comm2Coeff gxa gza gxb gzb xa za xb zb = 0
      ↔ anti2 gxa gza gxb gzb xa za xb zb = false) := by decide

/-- **The generic `rotate_2` sign table is the Pauli product phase.** For *every*
axis pair `(G_a, G_b)` and every two-site key `P`, the value the generic path
actually multiplies `sin` by — `−ε` with `ε` from `comm_2`, i.e. old's
`sin.mul_sign(-eps)` — equals the real branch prefactor
`i^{1 + phaseExp_a + phaseExp_b}` of `−i·[G_a ⊗ G_b, P]/2` (`branchExp2`, real by
`branchExp2_isRealPhase`).

So the `SIGN_NEG = 0x2840` mask *and* the `+eps`/`−eps` asymmetry between the
generic path and the `rzz`/`rxx`/`ryy` fast paths are both correct, and correct
off the diagonal too — the regime neither old's `*_matches_generic` tests nor
their ports ever exercise. -/
theorem comm2_generic_sign_eq_branchExp2 : ∀ gxa gza gxb gzb xa za xb zb,
    anti2 gxa gza gxb gzb xa za xb zb = true →
      -comm2Coeff gxa gza gxb gzb xa za xb zb
        = realSign (branchExp2 gxa gza gxb gzb xa za xb zb) := by decide

/-- **The generic path's branch key is the two-site product key.** `comm_2`'s
`present`-masked output bits are `mulBits2 = G ⊕ P` site-wise whenever the branch
fires, so the generic kernel re-keys exactly as the native ones do (and
`anticommute_new_key2` then applies to it verbatim). -/
theorem comm2_key_eq_mulBits2 : ∀ gxa gza gxb gzb xa za xb zb,
    anti2 gxa gza gxb gzb xa za xb zb = true →
      comm2Key gxa gza gxb gzb xa za xb zb = mulBits2 gxa gza gxb gzb xa za xb zb := by
  decide

end GenericTwo

/-! ### `RotXY::r` — the sub-rotation **order** is Heisenberg (backward)

`RotXY::r(q, φ, θ)` emits `rz(q, φ); rx(q, θ); rz(q, −φ)` in that order
(`crates/ppvm-pauli-sum-2/src/rotation.rs`, `impl RotXY`), the *reverse* of the
tableau's forward order. Behavioural contract 10 singles this out because a
forward-ordered implementation yields `ry(q, −θ)` at `φ = π/2` and passes every
*other* rotation test; up to now the order was pinned only by two ported example
values (`rot1.rs::test_r`, `tests/gate_surface.rs::r_is_heisenberg_ordered`) and
`grep -rl 'RotXY|rotXY' lean/PPVM/` was empty.

The per-axis `ε` columns above already fix each single-qubit rotation's action on
the Pauli basis; assembling them into a `3 × 3` real matrix on the coefficient
triple `(c_X, c_Y, c_Z)` — the site-restricted `C[K]` — makes the order claim a
short computation. The matrix entries are *read off the kernel*, not modeled:
`mz_from_kernel`/`mx_from_kernel`/`my_from_kernel` check each off-diagonal entry
against `mulBits` (which key the branch lands on) and `branchExp` (its `±1`),
i.e. against `rz_eps_from_product` / `rx_eps_from_product` / `ry_eps_from_product`. -/

section RotXY

open Real

/-- The coefficient triple `(c_X, c_Y, c_Z)` at one site — `C[K]` restricted to
the rotated qubit, in the basis the `ε` columns above are stated over. -/
abbrev Vec3 := ℝ × ℝ × ℝ

/-- **`rz(θ)`'s action.** From the `rz` column (`rz_eps_from_product`): `Z` is
inert, `X ↦ cosθ·X − sinθ·Y`, `Y ↦ cosθ·Y + sinθ·X`. -/
noncomputable def mz (θ : ℝ) (v : Vec3) : Vec3 :=
  (cos θ * v.1 + sin θ * v.2.1, -(sin θ) * v.1 + cos θ * v.2.1, v.2.2)

/-- **`rx(θ)`'s action.** From the `rx` column (`rx_eps_from_product`): `X` is
inert, `Z ↦ cosθ·Z + sinθ·Y`, `Y ↦ cosθ·Y − sinθ·Z`. -/
noncomputable def mx (θ : ℝ) (v : Vec3) : Vec3 :=
  (v.1, cos θ * v.2.1 + sin θ * v.2.2, -(sin θ) * v.2.1 + cos θ * v.2.2)

/-- **`ry(θ)`'s action.** From the `ry` column (`ry_eps_from_product`): `Y` is
inert, `X ↦ cosθ·X + sinθ·Z`, `Z ↦ cosθ·Z − sinθ·X`. -/
noncomputable def my (θ : ℝ) (v : Vec3) : Vec3 :=
  (cos θ * v.1 - sin θ * v.2.2, v.2.1, sin θ * v.1 + cos θ * v.2.2)

/-- **`mz`'s off-diagonal entries are the kernel's.** `rz` on `X = (1,0)` branches
to `mulBits = (1,1) = Y` with `branchExp = 2` (`ε = −1`), and on `Y = (1,1)`
branches to `(1,0) = X` with `branchExp = 0` (`ε = +1`) — the `−sinθ` and `+sinθ`
of `mz`. -/
theorem mz_from_kernel :
    mulBits false true true false = (true, true) ∧ branchExp false true true false = 2
      ∧ mulBits false true true true = (true, false)
      ∧ branchExp false true true true = 0 := by decide

/-- **`mx`'s off-diagonal entries are the kernel's.** `rx` on `Z = (0,1)` branches
to `(1,1) = Y` with `branchExp = 0` (`ε = +1`), and on `Y` branches to
`(0,1) = Z` with `branchExp = 2` (`ε = −1`). -/
theorem mx_from_kernel :
    mulBits true false false true = (true, true) ∧ branchExp true false false true = 0
      ∧ mulBits true false true true = (false, true)
      ∧ branchExp true false true true = 2 := by decide

/-- **`my`'s off-diagonal entries are the kernel's.** `ry` on `X = (1,0)` branches
to `(0,1) = Z` with `branchExp = 0` (`ε = +1`), and on `Z` branches to
`(1,0) = X` with `branchExp = 2` (`ε = −1`). -/
theorem my_from_kernel :
    mulBits true true true false = (false, true) ∧ branchExp true true true false = 0
      ∧ mulBits true true false true = (true, false)
      ∧ branchExp true true false true = 2 := by decide

/-- Rotation of the coefficient triple about a unit axis `n` by `θ`, in the same
(Heisenberg / backward) orientation as `mx`/`my`/`mz` — Rodrigues' formula with
angle `−θ`:  `v·cosθ − (n × v)·sinθ + n (n·v)(1 − cosθ)`. -/
noncomputable def rotAxis (n : Vec3) (θ : ℝ) (v : Vec3) : Vec3 :=
  ( cos θ * v.1 - sin θ * (n.2.1 * v.2.2 - n.2.2 * v.2.1)
      + n.1 * (n.1 * v.1 + n.2.1 * v.2.1 + n.2.2 * v.2.2) * (1 - cos θ),
    cos θ * v.2.1 - sin θ * (n.2.2 * v.1 - n.1 * v.2.2)
      + n.2.1 * (n.1 * v.1 + n.2.1 * v.2.1 + n.2.2 * v.2.2) * (1 - cos θ),
    cos θ * v.2.2 - sin θ * (n.1 * v.2.1 - n.2.1 * v.1)
      + n.2.2 * (n.1 * v.1 + n.2.1 * v.2.1 + n.2.2 * v.2.2) * (1 - cos θ) )

/-- **`RotXY::r` is rotation about the in-plane axis `cos φ·X + sin φ·Y`.**
Applying the crate's three sub-rotations in the order it emits them —
`rz(φ)` first, then `rx(θ)`, then `rz(−φ)`, so the composite map is
`M_z(−φ) ∘ M_x(θ) ∘ M_z(φ)` — is exactly `rotAxis (cos φ, sin φ, 0) θ`.

The forward (Schrödinger) order composes to the *inverse* rotation, which is why
contract 10 calls this out: only the backward order gives `r(q, π/2, θ) = ry(q, θ)`
rather than `ry(q, −θ)`. -/
theorem rotXY_heisenberg_order (φ θ : ℝ) (v : Vec3) :
    mz (-φ) (mx θ (mz φ v)) = rotAxis (cos φ, sin φ, 0) θ v := by
  have h : sin φ ^ 2 + cos φ ^ 2 = 1 := sin_sq_add_cos_sq φ
  refine Prod.ext ?_ (Prod.ext ?_ ?_) <;>
    simp only [mz, mx, rotAxis, cos_neg, sin_neg]
  · linear_combination (cos θ * v.1) * h
  · linear_combination (cos θ * v.2.1) * h
  · ring

/-- **`r(q, 0, θ) = rx(q, θ)`** — the `φ = 0` end of the family. -/
theorem rotXY_zero_eq_rx (θ : ℝ) (v : Vec3) : mz (-0) (mx θ (mz 0 v)) = mx θ v := by
  rw [rotXY_heisenberg_order]
  refine Prod.ext ?_ (Prod.ext ?_ ?_) <;> simp only [rotAxis, mx, cos_zero, sin_zero] <;> ring

/-- **`r(q, π/2, θ) = ry(q, θ)`** — the identity behavioural contract 10 names as
the order detector: a forward-ordered implementation returns `ry(q, −θ)` here and
is otherwise indistinguishable. -/
theorem rotXY_halfPi_eq_ry (θ : ℝ) (v : Vec3) :
    mz (-(π / 2)) (mx θ (mz (π / 2) v)) = my θ v := by
  rw [rotXY_heisenberg_order]
  refine Prod.ext ?_ (Prod.ext ?_ ?_) <;>
    simp only [rotAxis, my, cos_pi_div_two, sin_pi_div_two] <;> ring

end RotXY

end PPVM.Rotation
