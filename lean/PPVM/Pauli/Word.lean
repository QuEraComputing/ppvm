/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Pauli.Phase
import Mathlib.Algebra.BigOperators.Fin

/-!
# The n-qubit Pauli word: lifting the phase cocycle

`PPVM.PauliPhase` proved the single-qubit facts. The traits-2 design's key type
`PauliWord` is an n-qubit tensor product, and its packed kernel
(`phase/mul.rs`) accumulates a per-qubit phase and reduces the total mod 4:

```rust
sign_count += sign.count_ones();
imag_count += imag.count_ones();
// ...
self.add_phase(((2 * sign_count + imag_count) % 4) as u8);
```

So the n-qubit phase exponent is the **sum over qubits** of the single-qubit
`phaseExp`, and the n-qubit symplectic form is the sum of the single-qubit
`omega`. This file lifts the group cocycle and the commutation law to that sum,
which is exactly what makes the packed multi-qubit `MulAssign` well defined.

The product bits are pointwise (`x = a^c`, `z = b^d` per qubit).

Scope note: the single-qubit `phaseExp` is grounded against genuine `ℤ[i]`
matrices in `PPVM.PauliMatrix`. The n-qubit phase `phaseExpN` here is defined as
the *sum* of those single-qubit exponents (which is exactly what the Rust kernel
computes); it is *also* tied back to genuine `n`-fold tensor-product matrices in
`PPVM.PauliMatrix.tensorPauli_mul`, which proves `phaseExpN` is the base-`i`
exponent of the real `2ⁿ×2ⁿ` operator product `g(p)·g(q)` (via the
phase-multiplicativity `prod_iuPow`, `∏ᵢ iᵏⁱ = i^{Σ kᵢ}`).
-/

namespace PPVM.PauliWord

open PPVM.PauliPhase

variable {n : ℕ}

/-- An n-qubit Pauli word up to global phase: an `(x, z)` bit pair per qubit.
`i`-th slot is qubit `i`, matching the packed `xbits`/`zbits` planes. -/
abbrev Word (n : ℕ) := Fin n → Bool × Bool

/-- Pointwise product bits: `⊕` on each qubit's X and Z plane. -/
def mulWord (p q : Word n) : Word n :=
  fun i => mulBits (p i).1 (p i).2 (q i).1 (q i).2

/-- The n-qubit phase exponent: the sum over qubits of the per-qubit
`phaseExp`, reduced in `ℤ/4ℤ` — the spec of `(2*sign_count + imag_count) % 4`. -/
def phaseExpN (p q : Word n) : ZMod 4 :=
  ∑ i, phaseExp (p i).1 (p i).2 (q i).1 (q i).2

/-- The n-qubit symplectic form: sum over qubits of the single-qubit `omega`. -/
def omegaN (p q : Word n) : ZMod 4 :=
  ∑ i, omega (p i).1 (p i).2 (q i).1 (q i).2

/-- **The n-qubit 2-cocycle law.** The summed phase is still a 2-cocycle, so the
packed multi-qubit product is associative. Lifts `phaseExp_cocycle` termwise. -/
theorem phaseExpN_cocycle (p q r : Word n) :
    phaseExpN p q + phaseExpN (mulWord p q) r
      = phaseExpN q r + phaseExpN p (mulWord q r) := by
  simp only [phaseExpN, mulWord, mulBits]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ => phaseExp_cocycle _ _ _ _ _ _

/-- **The n-qubit commutation law.** `phaseExpN p q − phaseExpN q p = 2·ω(p,q)`,
i.e. `P·Q = (−1)^{ω(P,Q)} Q·P` for n-qubit words. Lifts `phaseExp_sub_comm`. -/
theorem phaseExpN_sub_comm (p q : Word n) :
    phaseExpN p q - phaseExpN q p = 2 * omegaN p q := by
  simp only [phaseExpN, omegaN, Finset.mul_sum]
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun i _ => phaseExp_sub_comm _ _ _ _

/-- Product bits are involutive per qubit, so `P · P` has identity bits. -/
theorem mulWord_self (p : Word n) : mulWord p p = fun _ => (false, false) := by
  funext i
  simp only [mulWord, mulBits, Bool.xor_self]

/-- **The product bits are the identity word exactly on the diagonal.**
`mulWord p q = I ↔ q = p` — the converse of `mulWord_self`. Each qubit's
`(x, z)` pair is `⊕`-cancelling, and `⊕` is injective, so the only word whose
product with `p` has trivial bits is `p` itself. This is what collapses the L4
outer product to its diagonal when one reads off the identity coefficient (see
`PPVM.Twisted.twistedConv_apply_id`). -/
theorem mulWord_eq_id_iff (p q : Word n) :
    mulWord p q = (fun _ => (false, false)) ↔ q = p := by
  have xor_cancel : ∀ a b : Bool, xor a b = false → b = a := by decide
  constructor
  · intro h
    funext i
    have hi := congrFun h i
    simp only [mulWord, mulBits, Prod.mk.injEq] at hi
    exact Prod.ext (xor_cancel _ _ hi.1) (xor_cancel _ _ hi.2)
  · rintro rfl
    exact mulWord_self _

/-- **The product bits are associative.** Per qubit both planes are `⊕`, and
`⊕` is associative, so `(P·Q)·R` and `P·(Q·R)` carry the same bits. Together
with `phaseExpN_cocycle` (the phase half) this is what makes the packed n-qubit
product associative; it is the `kmul`-associativity hypothesis of
`PPVM.Twisted.gtmul_assoc` for the Pauli key. -/
theorem mulWord_assoc (p q r : Word n) :
    mulWord (mulWord p q) r = mulWord p (mulWord q r) := by
  funext i
  simp only [mulWord, mulBits, Prod.mk.injEq]
  exact ⟨Bool.xor_assoc _ _ _, Bool.xor_assoc _ _ _⟩

/-- **Right multiplication by a fixed word is injective.**
`P₁ · Q = P₂ · Q → P₁ = P₂`, because each plane is `p ⊕ q` and `⊕ q` is a
bijection on `Bool`.

This is the machine-checked licence for the `RekeyBijective` fast path in the L4
`Multiply` component: `Sum::mul_word_assign` re-keys every term through
`k ↦ mulWord k q` and merges with a *plain* `insert` plus a `debug_assert!`,
because injectivity means two distinct source keys can never land on the same
destination key. (Injectivity is a real precondition there, not a hint: a
collision would silently drop a term in release rather than sum it.) It is the
word-product analogue of the Clifford re-key's `PPVM.Symplectic.*_bijective` /
`PPVM.Conjugation.conj*_injective`. -/
theorem mulWord_right_injective (q : Word n) :
    Function.Injective (fun p : Word n => mulWord p q) := by
  have xor_cancel : ∀ a b c : Bool, xor a c = xor b c → a = b := by decide
  intro p₁ p₂ h
  funext i
  have hi := congrFun h i
  simp only [mulWord, mulBits, Prod.mk.injEq] at hi
  exact Prod.ext (xor_cancel _ _ _ hi.1) (xor_cancel _ _ _ hi.2)

/-- `P · P` is phase-free: the summed self-phase is `0`. Combined with
`mulWord_self`, `P² = +I`. -/
theorem phaseExpN_self (p : Word n) : phaseExpN p p = 0 := by
  simp only [phaseExpN]
  rw [Finset.sum_eq_zero]
  exact fun i _ => phaseExp_self _ _

/-! ### Canonical form `iᵠ Xˣ Zᶻ` and the multi-site conjugation fold

`compute_decomposition_word` (`crates/ppvm-tableau-2/src/data.rs`) is the sole
entry point for multi-qubit `expectation(&word)`, and hence for every golden
expectation value in the acceptance bar (Bell `⟨YY⟩ = −1`, GHZ `⟨XXX⟩ = 1`, the
`trace` sums). It conjugates a Pauli **word** through the frame by conjugating
each site separately with `compute_decomposition`, then folding the single-site
results left-to-right in the **canonical form** `iᵠ Xˣ Zᶻ`:

```rust
let cross = 2 * (symplectic_inner(destab_anticomm, q_stab) as u8 % 2);
phase = (phase + q_phase + cross) % 4;
stab_anticomm  ^= q_stab;    // the X plane
destab_anticomm ^= q_destab; // the Z plane
```

The `cross` term is the commutation correction of
`Z^{z_a} X^{x_b} = (−1)^{z_a·x_b} X^{x_b} Z^{z_a}`, taken between the **running**
Z-mask and the **new** site's X-mask. The claim that folding this way reproduces
the true conjugate of the product is a statement neither
`Conjugation.lean` (per-generator signs) nor `phaseExpN` (the two-word product in
the `g(x,z) = i^{x·z} Xˣ Zᶻ` normalization) makes, and it has no single-qubit
oracle to fall back on — the cross term vanishes identically at weight 1. Two
theorems settle it:

* `phaseExpN_eq_canon` — the canonical-form cross-phase rule **is** the genuine
  Pauli product phase, once the two normalizations are reconciled. Since
  `phaseExpN` is grounded in real `2ⁿ×2ⁿ` matrices
  (`PPVM.PauliMatrix.tensorPauli_mul`), this ties the `cross` term to actual
  operator multiplication.
* `Canon.foldl_eq_prod` — the running left-fold the Rust performs equals the
  ordered product of the per-site factors (the fold is associative and unital),
  so accumulating site by site is the same as multiplying the whole word at once.
-/

/-- `∑ᵢ xᵢ·zᵢ` in `ℤ/4ℤ`: the exponent of the `iˣᶻ` normalization relating the
crate's canonical `Xˣ Zᶻ` form to the `g(x,z)` form `phaseExpN` is stated for. -/
def xzOverlap (w : Word n) : ZMod 4 := ∑ i, b4 (w i).1 * b4 (w i).2

/-- The Rust `cross` term: `2 · popcount(z_a ∧ x_b)`, i.e. the sign of
`Z^{z_a} X^{x_b} = (−1)^{z_a·x_b} X^{x_b} Z^{z_a}`. -/
def crossPhase (p q : Word n) : ZMod 4 := ∑ i, 2 * (b4 (p i).2 * b4 (q i).1)

/-- **The canonical cross-phase rule is the genuine Pauli product phase.**
Writing `A = i^{φ_a} X^{x_a} Z^{z_a}` and `B = i^{φ_b} X^{x_b} Z^{z_b}`, moving
`Z^{z_a}` past `X^{x_b}` gives `A·B = i^{φ_a + φ_b + cross} X^{x_a⊕x_b} Z^{z_a⊕z_b}`.
Converting both sides into the `g(x,z) = i^{x·z} Xˣ Zᶻ` normalization turns that
into the identity below, whose per-site content is exactly `phaseRef` — hence
`phaseExp`, hence (via `PPVM.PauliMatrix.tensorPauli_mul`) real matrices. So the
`(−1)^{popcount(z_running ∧ x_new)}` correction in `compute_decomposition_word`
is not a convention but the operator product. -/
theorem phaseExpN_eq_canon (p q : Word n) :
    phaseExpN p q
      = xzOverlap p + xzOverlap q + crossPhase p q - xzOverlap (mulWord p q) := by
  simp only [phaseExpN, xzOverlap, crossPhase, mulWord, mulBits, phaseExp_eq_ref,
    phaseRef]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]

/-- A Pauli word in the crate's **canonical form** `iᵠ Xˣ Zᶻ`: the triple
`compute_decomposition_word` accumulates (`phase`, and the `stab_anticomm` /
`destab_anticomm` masks that are its X and Z planes). -/
@[ext]
structure Canon (n : ℕ) where
  /-- Base-`i` phase exponent. -/
  phase : ZMod 4
  /-- The `(x, z)` bit planes. -/
  word : Word n

/-- The fold step: phases add, the cross term corrects the `Z`-past-`X`
commutation, and the planes `⊕`. Verbatim the Rust loop body. -/
def Canon.mul (a b : Canon n) : Canon n where
  phase := a.phase + b.phase + crossPhase a.word b.word
  word := mulWord a.word b.word

/-- The fold's seed: `phase = 0`, both planes empty. -/
def Canon.one : Canon n where
  phase := 0
  word := fun _ => (false, false)

/-- The identity word contributes no cross phase on the left (its Z plane is
empty). -/
theorem crossPhase_zero_left (q : Word n) :
    crossPhase (fun _ => (false, false)) q = 0 := by
  refine Finset.sum_eq_zero fun i _ => ?_
  simp [b4]

/-- …nor on the right (its X plane is empty). -/
theorem crossPhase_zero_right (p : Word n) :
    crossPhase p (fun _ => (false, false)) = 0 := by
  refine Finset.sum_eq_zero fun i _ => ?_
  simp [b4]

theorem Canon.one_mul (a : Canon n) : Canon.mul Canon.one a = a := by
  ext i <;>
    simp only [Canon.mul, Canon.one, crossPhase_zero_left, zero_add, add_zero, mulWord,
      mulBits, Bool.false_xor]

theorem Canon.mul_one (a : Canon n) : Canon.mul a Canon.one = a := by
  ext i <;>
    simp only [Canon.mul, Canon.one, crossPhase_zero_right, add_zero, mulWord, mulBits,
      Bool.xor_false]

private theorem crossPhase_cocycle_site (zp xq zq xr : Bool) :
    2 * (b4 zp * b4 xq) + 2 * (b4 (xor zp zq) * b4 xr)
      = 2 * (b4 zq * b4 xr) + 2 * (b4 zp * b4 (xor xq xr)) := by
  cases zp <;> cases xq <;> cases zq <;> cases xr <;> decide

/-- **The cross-phase is a 2-cocycle**, so the canonical-form product is
associative — the fold may accumulate sites in any bracketing. Note this is a
*different* cocycle from `phaseExp` (no `iˣᶻ` normalization); over `ℤ/4ℤ` the
`2 ·` factor is what lets `⊕` be replaced by `+` inside it. -/
theorem crossPhase_cocycle (p q r : Word n) :
    crossPhase p q + crossPhase (mulWord p q) r
      = crossPhase q r + crossPhase p (mulWord q r) := by
  simp only [crossPhase, mulWord, mulBits]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ => crossPhase_cocycle_site _ _ _ _

theorem Canon.mul_assoc (a b c : Canon n) :
    Canon.mul (Canon.mul a b) c = Canon.mul a (Canon.mul b c) := by
  have h := crossPhase_cocycle a.word b.word c.word
  ext i
  · simp only [Canon.mul]; linear_combination h
  · simp only [Canon.mul, mulWord_assoc]
  · simp only [Canon.mul, mulWord_assoc]

/-- The ordered product of a list of canonical factors. -/
def Canon.prod (l : List (Canon n)) : Canon n := l.foldr Canon.mul Canon.one

theorem Canon.foldl_mul (a : Canon n) :
    ∀ l : List (Canon n), l.foldl Canon.mul a = Canon.mul a (Canon.prod l)
  | [] => (Canon.mul_one a).symm
  | x :: xs => by
    rw [List.foldl_cons, Canon.foldl_mul (Canon.mul a x) xs, Canon.mul_assoc]
    rfl

/-- **The running fold is the ordered product.** `compute_decomposition_word`
accumulates `(phase, x, z)` site by site, taking the cross term against the
*running* Z-mask; that left-fold equals the product of the per-site canonical
conjugates. Together with `phaseExpN_eq_canon` this is the machine-checked
statement that the folded phase is the phase of the genuine conjugated word. -/
theorem Canon.foldl_eq_prod (l : List (Canon n)) :
    l.foldl Canon.mul Canon.one = Canon.prod l := by
  rw [Canon.foldl_mul, Canon.one_mul]

/-- The canonical phase re-expressed in the `g(x,z)` normalization
`phaseExpN` uses. -/
def Canon.toG (a : Canon n) : ZMod 4 := a.phase - xzOverlap a.word

/-- **The fold step, read in the `g` normalization, is the genuine Pauli
product.** `toG (a·b) = toG a + toG b + phaseExpN a b` — i.e. the Rust's
`phase + q_phase + cross` accumulation carries exactly the phase that the
matrix-grounded n-qubit Pauli product carries. -/
theorem Canon.toG_mul (a b : Canon n) :
    (Canon.mul a b).toG = a.toG + b.toG + phaseExpN a.word b.word := by
  simp only [Canon.toG, Canon.mul, phaseExpN_eq_canon]
  ring

end PPVM.PauliWord
