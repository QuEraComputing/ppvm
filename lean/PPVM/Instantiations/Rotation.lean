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

end TwoPass

end PPVM.Rotation
