/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Pauli.Conjugation
import PPVM.Pauli.Word

/-!
# Fused batch gates: a batched sign is a site-parity

Every fused batch path in the tableau crate — `x_many` / `y_many` / `z_many` /
`s_many` / `h_many` / `s_dag_many` / `sqrt_{x,y}[_dag]_many`
(`crates/ppvm-tableau-2/src/clifford.rs`) and the `cz_block_pairs` /
`cz_block_pairs_cross_word` / `cz_many` / `cy_many` family
(`crates/ppvm-tableau-2/src/data.rs`) — replaces a *sequence* of `ℤ/4` phase
updates with a **single** `count_ones() & 1` parity, read off the row *before*
any of the gates in the batch has run. That substitution is an algebraic claim,
not a SIMD detail, and it has two parts:

1. the per-site sign contribution is `ℤ/2`-valued (a `+2` in `ℤ/4`), so `2+2 ≡ 0`
   and XOR-accumulating the parity is sound (`two_mul_natCast`); and
2. the sites are *independent*: a gate at site `q` must not change what a later
   gate at site `q' ≠ q` reads. This is what makes summation order irrelevant and
   licenses reading every sign off the original row.

Part 2 is a genuine side condition, not decoration. `seqApply_eq_batchApply`
carries `L.Nodup`; `czSeq_phase` carries `List.Pairwise Disjoint2`; and
`czSeq_phase_needs_disjoint` exhibits two `CZ` pairs sharing an endpoint on which
the sequential loop and the batched parity **disagree**, because a pair's sign
reads a `z`-bit an earlier pair already rewrote. That is exactly the failure mode
a `cz_block_pairs(base, offset, count)` call with `offset < count` would hit, so
the precondition belongs in the record rather than in an unwritten assumption of
a shipped public API.

The theorems below say nothing about words, masks or storage layout — only that a
`ℤ/2` character summed over independent sites equals its per-site parity.
-/

namespace PPVM.Tableau.Batch

open PPVM.PauliWord PPVM.PauliPhase PPVM.PauliPhase.PhasedPauli

variable {n : ℕ}

/-- A tableau row: an `n`-qubit Pauli word together with its `ℤ/4` phase
(`PhasedPauliWordNoHash` in the crate). -/
abbrev Row (n : ℕ) := ZMod 4 × Word n

/-! ### A `ℤ/2`-valued phase accumulates as a parity -/

/-- **`2 · m` in `ℤ/4` only sees `m mod 2`.** This is the `count_ones() & 1`
step: a batch of `m` sign contributions of `+2` each collapses to a single `+2`
when `m` is odd and to nothing when `m` is even. -/
theorem two_mul_natCast (m : ℕ) : 2 * (m : ZMod 4) = if m % 2 = 1 then 2 else 0 := by
  have hc : ((2 * m : ℕ) : ZMod 4) = 2 * (m : ZMod 4) := by push_cast; ring
  have h : ((2 * m : ℕ) : ZMod 4) = ((2 * m % 4 : ℕ) : ZMod 4) := (ZMod.natCast_mod _ _).symm
  have h2 : 2 * m % 4 = 2 * (m % 2) := by omega
  rcases Nat.mod_two_eq_zero_or_one m with hm | hm <;>
    rw [← hc, h, h2, hm] <;> simp

/-- The `ℤ/4` delta of one sign bit. -/
theorem two_mul_ite (b : Bool) :
    2 * ((if b then 1 else 0 : ℕ) : ZMod 4) = if b then 2 else 0 := by
  cases b <;> simp

/-! ### Single-qubit batches: `x_many` … `sqrt_y_dag_many` -/

/-- Apply a sitewise gate (`bits`, with sign predicate `sign`) at one site. -/
def siteApply (bits : Bool × Bool → Bool × Bool) (sign : Bool × Bool → Bool)
    (q : Fin n) (r : Row n) : Row n :=
  (r.1 + (if sign (r.2 q) then 2 else 0), Function.update r.2 q (bits (r.2 q)))

/-- The unfused path: apply the gate at each site of `L`, in order. -/
def seqApply (bits : Bool × Bool → Bool × Bool) (sign : Bool × Bool → Bool)
    (L : List (Fin n)) (r : Row n) : Row n :=
  L.foldl (fun r q => siteApply bits sign q r) r

/-- The fused mask sweep: **one** phase flip for the parity of the sign predicate
over the whole index set, and the sitewise bit map applied to every masked site —
everything read off the *original* row. -/
def batchApply (bits : Bool × Bool → Bool × Bool) (sign : Bool × Bool → Bool)
    (L : List (Fin n)) (r : Row n) : Row n :=
  (r.1 + 2 * ((L.countP fun q => sign (r.2 q) : ℕ) : ZMod 4),
   fun k => if k ∈ L then bits (r.2 k) else r.2 k)

/-- **The fused sweep equals the per-site loop, on DISTINCT sites.** The batched
`count_ones() & 1` phase and the sitewise bit map reproduce the sequential `ℤ/4`
updates exactly; combined with `two_mul_natCast` this is the whole justification
of `x_many`/`y_many`/`z_many`/`s_many`/`h_many`/`sqrt_*_many`. -/
theorem seqApply_eq_batchApply (bits : Bool × Bool → Bool × Bool)
    (sign : Bool × Bool → Bool) :
    ∀ (L : List (Fin n)), L.Nodup → ∀ r : Row n,
      seqApply bits sign L r = batchApply bits sign L r := by
  intro L
  induction L with
  | nil =>
    intro _ r
    simp only [seqApply, batchApply, List.foldl_nil, List.countP_nil, Nat.cast_zero,
      mul_zero, add_zero, List.not_mem_nil, if_false]
  | cons q qs ih =>
    intro hnd r
    rw [List.nodup_cons] at hnd
    obtain ⟨hq, hqs⟩ := hnd
    have hupd : ∀ k, k ∈ qs → (Function.update r.2 q (bits (r.2 q))) k = r.2 k := by
      intro k hk
      refine Function.update_of_ne (fun h => ?_) _ _
      subst h
      exact hq hk
    have hstep : seqApply bits sign (q :: qs) r
        = batchApply bits sign qs (siteApply bits sign q r) := by
      rw [seqApply, List.foldl_cons]
      exact ih hqs _
    rw [hstep]
    refine Prod.ext ?_ ?_
    · change (siteApply bits sign q r).1 + 2 * ((qs.countP fun k =>
        sign ((siteApply bits sign q r).2 k) : ℕ) : ZMod 4) = _
      have hcount : (qs.countP fun k => sign ((siteApply bits sign q r).2 k))
          = qs.countP fun k => sign (r.2 k) :=
        List.countP_congr fun k hk => by simp only [siteApply, hupd k hk]
      rw [hcount]
      change r.1 + (if sign (r.2 q) then 2 else 0) + _ = _
      simp only [batchApply, List.countP_cons, Nat.cast_add, mul_add, two_mul_ite]
      ring
    · funext k
      by_cases hk : k ∈ qs
      · simp only [batchApply, siteApply, hk, if_true, List.mem_cons, or_true, hupd k hk]
      · by_cases hkq : k = q
        · subst hkq
          simp only [batchApply, siteApply, hk, if_false, List.mem_cons, true_or, if_true,
            Function.update_self]
        · simp only [batchApply, siteApply, hk, if_false, List.mem_cons, hkq, or_self,
            if_false, Function.update_of_ne hkq]

/-! ### The concrete single-qubit sign predicates

`siteApply` is parametric in `(bits, sign)`; the crate's batched gates use the
sign predicates below. Each is read off — and checked against — the *audited*
`PhasedPauli` conjugation maps of `Pauli/Conjugation.lean`, so the batch paths
inherit those sign tables rather than re-deriving them. Each `decide` also
certifies that the map really is of the sitewise "add `2·sign` to the phase,
rewrite the bits" shape that `siteApply` assumes.

**Direction note.** The `ext*` maps are the *forward* conjugation `G·P·G†`; the
tableau rows are conjugated in the **backward** (Heisenberg) direction, so the
crate's gate names pair with the adjoint Lean map. Concretely, against
`crates/ppvm-tableau-2/src/clifford.rs`:

| crate batch gate | phase mask         | Lean witness              |
|------------------|--------------------|---------------------------|
| `h_many`         | `x & z`            | `isSitewise_conjH`        |
| `s_many`         | `x & z`            | `isSitewise_conjS`        |
| `x_many`         | `z`                | `isSitewise_conjX`        |
| `y_many`         | `x ^ z`            | `isSitewise_conjY`        |
| `z_many`         | `x`                | `isSitewise_conjZ`        |
| `sqrt_x_many`    | `z & !x`           | `isSitewise_extSqrtXdag`  |
| `sqrt_x_dag_many`| `x & z`            | `isSitewise_extSqrtX`     |
| `sqrt_y_many`    | `x & !z`           | `isSitewise_extSqrtYdag`  |
| `sqrt_y_dag_many`| `z & !x`           | `isSitewise_extSqrtY`     |

Every entry is a `ℤ/2`-valued per-site predicate, which is precisely hypothesis 1
above. -/

/-- A single-qubit conjugation is *sitewise*: it adds `2·sign(x,z)` to the phase
and rewrites `(x,z)`, uniformly in the incoming phase. -/
def IsSitewise (g : PhasedPauli → PhasedPauli) (bits : Bool × Bool → Bool × Bool)
    (sign : Bool × Bool → Bool) : Prop :=
  ∀ p : PhasedPauli,
    g p = ⟨p.phase + (if sign (p.x, p.z) then 2 else 0), (bits (p.x, p.z)).1,
           (bits (p.x, p.z)).2⟩

/-- `H`: swap the bits, sign `x ∧ z`. -/
theorem isSitewise_conjH :
    IsSitewise conjH (fun b => (b.2, b.1)) (fun b => b.1 && b.2) := by
  intro p; revert p; decide

/-- `S`: `z ⊕= x`, sign `x ∧ z`. -/
theorem isSitewise_conjS :
    IsSitewise conjS (fun b => (b.1, xor b.1 b.2)) (fun b => b.1 && b.2) := by
  intro p; revert p; decide

/-- `S†` (the backward direction the simulator runs): `z ⊕= x`, sign `x ∧ ¬z`. -/
theorem isSitewise_conjSdag :
    IsSitewise conjSdag (fun b => (b.1, xor b.1 b.2)) (fun b => b.1 && !b.2) := by
  intro p; revert p; decide

/-- `X`: bits unchanged, sign `z`. -/
theorem isSitewise_conjX :
    IsSitewise (fun p => ⟨p.phase + (if p.z then 2 else 0), p.x, p.z⟩) id
      (fun b => b.2) := by
  intro p; revert p; decide

/-- `Y`: bits unchanged, sign `x ⊕ z`. -/
theorem isSitewise_conjY :
    IsSitewise (fun p => ⟨p.phase + (if xor p.x p.z then 2 else 0), p.x, p.z⟩) id
      (fun b => xor b.1 b.2) := by
  intro p; revert p; decide

/-- `Z`: bits unchanged, sign `x`. -/
theorem isSitewise_conjZ :
    IsSitewise (fun p => ⟨p.phase + (if p.x then 2 else 0), p.x, p.z⟩) id
      (fun b => b.1) := by
  intro p; revert p; decide

/-- `√X` (`extSqrtX`, table `X ↦ X`, `Y ↦ −Z`, `Z ↦ Y`): sign `x ∧ z`. -/
theorem isSitewise_extSqrtX :
    IsSitewise extSqrtX (fun b => (xor b.1 b.2, b.2)) (fun b => b.1 && b.2) := by
  intro p; revert p; decide

/-- `√X†` (`X ↦ X`, `Y ↦ Z`, `Z ↦ −Y`): sign `z ∧ ¬x`. -/
theorem isSitewise_extSqrtXdag :
    IsSitewise extSqrtXdag (fun b => (xor b.1 b.2, b.2)) (fun b => b.2 && !b.1) := by
  intro p; revert p; decide

/-- `√Y` (`X ↦ Z`, `Y ↦ Y`, `Z ↦ −X`): swap the bits, sign `z ∧ ¬x`. -/
theorem isSitewise_extSqrtY :
    IsSitewise extSqrtY (fun b => (b.2, b.1)) (fun b => b.2 && !b.1) := by
  intro p; revert p; decide

/-- `√Y†` (`X ↦ −Z`, `Y ↦ Y`, `Z ↦ X`): swap the bits, sign `x ∧ ¬z`. -/
theorem isSitewise_extSqrtYdag :
    IsSitewise extSqrtYdag (fun b => (b.2, b.1)) (fun b => b.1 && !b.2) := by
  intro p; revert p; decide

/-! ### Two-qubit batches: `cz_block_pairs`, `cz_many`, `cy_many`

For a family of pairs the batched form is the same claim one level up: the total
sign is the XOR of the per-pair signs, all read off the original row. Here the
independence hypothesis has real teeth, because a `CZ` writes `z`-bits and the
sign predicate *reads* `z`-bits. -/

/-- The `CZ` sign on a pair — the crate's `phase_bits = xc & xt & (zc ^ zt)`. -/
def czSign (c t : Fin n) (w : Word n) : Bool :=
  (w c).1 && (w t).1 && xor (w c).2 (w t).2

/-- The `CZ` bit rule on a pair: `z_c ⊕= x_t`, `z_t ⊕= x_c`. -/
def czBits (c t : Fin n) (w : Word n) : Word n := fun k =>
  if k = c then ((w c).1, xor (w c).2 (w t).1)
  else if k = t then ((w t).1, xor (w t).2 (w c).1)
  else w k

theorem czBits_of_ne (c t : Fin n) (w : Word n) {k : Fin n} (hc : k ≠ c) (ht : k ≠ t) :
    czBits c t w k = w k := by simp only [czBits, hc, ht, if_false]

/-- One `CZ`. -/
def czApply (c t : Fin n) (r : Row n) : Row n :=
  (r.1 + (if czSign c t r.2 then 2 else 0), czBits c t r.2)

/-- The unfused path: `CZ` on each pair in order. -/
def czSeq (P : List (Fin n × Fin n)) (r : Row n) : Row n :=
  P.foldl (fun r p => czApply p.1 p.2 r) r

/-- Two pairs have disjoint supports. -/
def Disjoint2 (p q : Fin n × Fin n) : Prop :=
  p.1 ≠ q.1 ∧ p.1 ≠ q.2 ∧ p.2 ≠ q.1 ∧ p.2 ≠ q.2

/-- **The batched two-qubit sign is the XOR of the per-pair signs**, provided the
pairs have pairwise-disjoint supports — the `cz_block_pairs` / `cz_many` /
`cy_many` phase accumulator. Every sign is read off the *original* row. -/
theorem czSeq_phase :
    ∀ (P : List (Fin n × Fin n)), P.Pairwise Disjoint2 → ∀ r : Row n,
      (czSeq P r).1 = r.1 + 2 * ((P.countP fun p => czSign p.1 p.2 r.2 : ℕ) : ZMod 4) := by
  intro P
  induction P with
  | nil => intro _ r; simp [czSeq]
  | cons p Q ih =>
    intro hpw r
    rw [List.pairwise_cons] at hpw
    obtain ⟨hdis, hQ⟩ := hpw
    have hsame : ∀ q ∈ Q, czSign q.1 q.2 (czBits p.1 p.2 r.2) = czSign q.1 q.2 r.2 := by
      intro q hq
      obtain ⟨h11, h12, h21, h22⟩ := hdis q hq
      simp only [czSign, czBits_of_ne _ _ _ (Ne.symm h11) (Ne.symm h21),
        czBits_of_ne _ _ _ (Ne.symm h12) (Ne.symm h22)]
    have hstep : czSeq (p :: Q) r = czSeq Q (czApply p.1 p.2 r) := by
      rw [czSeq, List.foldl_cons]; rfl
    rw [hstep, ih hQ]
    have hcount : (Q.countP fun q => czSign q.1 q.2 (czApply p.1 p.2 r).2)
        = Q.countP fun q => czSign q.1 q.2 r.2 :=
      List.countP_congr fun q hq => by simp only [czApply, hsame q hq]
    rw [hcount]
    change r.1 + (if czSign p.1 p.2 r.2 then 2 else 0) + _ = _
    simp only [List.countP_cons, Nat.cast_add, mul_add, two_mul_ite]
    ring

/-- **Disjointness is necessary.** On `x = 111`, `z = 000` the two overlapping
pairs `(0,1)` and `(1,2)` make the sequential loop flip the phase (the second
pair's sign reads the `z`-bit the first pair just wrote) while the batched parity
counts **zero** signs. So a `cz_block_pairs`-style fused sweep is *only* valid
under the disjoint-support precondition. -/
theorem czSeq_phase_needs_disjoint :
    (czSeq [((0 : Fin 3), (1 : Fin 3)), (1, 2)] (0, fun _ => (true, false))).1 = 2 ∧
      ([((0 : Fin 3), (1 : Fin 3)), (1, 2)].countP
        fun p => czSign p.1 p.2 (fun _ => (true, false))) = 0 := by
  constructor <;> decide

end PPVM.Tableau.Batch
