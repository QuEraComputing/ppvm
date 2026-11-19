/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
/-!
The ground-truth model for `ppvm-pauli-kernels`: exact matrices over the
Gaussian integers, the Pauli matrices, and every supported Clifford gate as a
scaled unitary `G = M / √c` with Gaussian-integer entries.

Everything here is computable, so all facts about it are checked by `decide`
in the kernel. No Mathlib.
-/

namespace PPVM.Kernels

/-! ## Gaussian integers -/

/-- `re + im·i` with integer parts. -/
structure GI where
  re : Int
  im : Int
  deriving DecidableEq, Repr

namespace GI

def zero : GI := ⟨0, 0⟩
def one : GI := ⟨1, 0⟩
def I : GI := ⟨0, 1⟩
def ofInt (n : Int) : GI := ⟨n, 0⟩

def add (a b : GI) : GI := ⟨a.re + b.re, a.im + b.im⟩
def neg (a : GI) : GI := ⟨-a.re, -a.im⟩
def mul (a b : GI) : GI := ⟨a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re⟩
def conj (a : GI) : GI := ⟨a.re, -a.im⟩

instance : Add GI := ⟨add⟩
instance : Neg GI := ⟨neg⟩
instance : Mul GI := ⟨mul⟩
instance : OfNat GI 0 := ⟨zero⟩
instance : OfNat GI 1 := ⟨one⟩

/-- `iᵏ` for `k` taken modulo four. -/
def iPow (k : Nat) : GI :=
  match k % 4 with
  | 0 => ⟨1, 0⟩
  | 1 => ⟨0, 1⟩
  | 2 => ⟨-1, 0⟩
  | _ => ⟨0, -1⟩

end GI

/-! ## Matrices -/

/-- An `n × n` matrix with Gaussian-integer entries. -/
abbrev Mat (n : Nat) := Fin n → Fin n → GI

namespace Mat

/-- Row-major entries, for decidable equality and printing. -/
def toList {n : Nat} (A : Mat n) : List (List GI) :=
  (List.finRange n).map fun i => (List.finRange n).map fun j => A i j

def mul {n : Nat} (A B : Mat n) : Mat n :=
  fun i j => (List.finRange n).foldl (fun acc k => acc + A i k * B k j) 0

/-- Conjugate transpose. -/
def adj {n : Nat} (A : Mat n) : Mat n := fun i j => (A j i).conj

def smul {n : Nat} (c : GI) (A : Mat n) : Mat n := fun i j => c * A i j

def id (n : Nat) : Mat n := fun i j => if i = j then 1 else 0

/-- Kronecker product; the first factor is the more significant index bit. -/
def kron (A B : Mat 2) : Mat 4 :=
  fun i j =>
    A ⟨i.val / 2, by omega⟩ ⟨j.val / 2, by omega⟩ * B ⟨i.val % 2, by omega⟩ ⟨j.val % 2, by omega⟩

def of2 (a b c d : GI) : Mat 2 :=
  fun i j =>
    match i.val, j.val with
    | 0, 0 => a | 0, _ => b | _, 0 => c | _, _ => d

def of4 (rows : List (List GI)) : Mat 4 :=
  fun i j => (rows.getD i.val []).getD j.val 0

end Mat

/-- Matrix equality, decided on entries. -/
def Mat.Eq {n : Nat} (A B : Mat n) : Prop := A.toList = B.toList

instance {n : Nat} (A B : Mat n) : Decidable (Mat.Eq A B) :=
  inferInstanceAs (Decidable (A.toList = B.toList))

/-! ## Paulis -/

/-- `(−1)^f` as a Gaussian integer. -/
def sgn (f : Bool) : GI := if f then ⟨-1, 0⟩ else ⟨1, 0⟩

/-- The Hermitian Pauli with bits `(x, z)`: I = (0,0), X = (1,0), Z = (0,1),
and the Hermitian Y = (1,1) (not `XZ`). -/
def pauli (x z : Bool) : Mat 2 :=
  match x, z with
  | false, false => Mat.of2 1 0 0 1
  | true, false => Mat.of2 0 1 1 0
  | false, true => Mat.of2 1 0 0 ⟨-1, 0⟩
  | true, true => Mat.of2 0 ⟨0, -1⟩ ⟨0, 1⟩ 0

/-- A two-qubit Pauli; `(xa, za)` is the more significant qubit. -/
def pauli2 (xa za xb zb : Bool) : Mat 4 := Mat.kron (pauli xa za) (pauli xb zb)

/-! ## Gates

Each gate is `G = M / √c` for an integer matrix `M` with `M† M = c · I`
(proved in `Rules.lean`). The Heisenberg action is then
`G† P G = (M† P M) / c`. -/

/-- A gate as a scaled unitary: `G = mat / √scale`. -/
structure Gate (n : Nat) where
  mat : Mat n
  scale : Int

def gX : Gate 2 := ⟨pauli true false, 1⟩
def gY : Gate 2 := ⟨pauli true true, 1⟩
def gZ : Gate 2 := ⟨pauli false true, 1⟩
def gH : Gate 2 := ⟨Mat.of2 1 1 1 ⟨-1, 0⟩, 2⟩
def gS : Gate 2 := ⟨Mat.of2 1 0 0 ⟨0, 1⟩, 1⟩
def gSdag : Gate 2 := ⟨Mat.of2 1 0 0 ⟨0, -1⟩, 1⟩
/-- `√X = ½ [[1+i, 1−i], [1−i, 1+i]]`. -/
def gSqrtX : Gate 2 := ⟨Mat.of2 ⟨1, 1⟩ ⟨1, -1⟩ ⟨1, -1⟩ ⟨1, 1⟩, 4⟩
def gSqrtXdag : Gate 2 := ⟨Mat.of2 ⟨1, -1⟩ ⟨1, 1⟩ ⟨1, 1⟩ ⟨1, -1⟩, 4⟩
/-- `√Y = ½ [[1+i, −1−i], [1+i, 1+i]]`. -/
def gSqrtY : Gate 2 := ⟨Mat.of2 ⟨1, 1⟩ ⟨-1, -1⟩ ⟨1, 1⟩ ⟨1, 1⟩, 4⟩
def gSqrtYdag : Gate 2 := ⟨Mat.of2 ⟨1, -1⟩ ⟨1, -1⟩ ⟨-1, 1⟩ ⟨1, -1⟩, 4⟩

private def o : GI := 1
private def z : GI := 0
private def m : GI := ⟨-1, 0⟩

/-- `CX` with the control as the more significant qubit. -/
def gCX : Gate 4 := ⟨Mat.of4 [[o, z, z, z], [z, o, z, z], [z, z, z, o], [z, z, o, z]], 1⟩
def gCZ : Gate 4 := ⟨Mat.of4 [[o, z, z, z], [z, o, z, z], [z, z, o, z], [z, z, z, m]], 1⟩
/-- `CY = |0⟩⟨0| ⊗ I + |1⟩⟨1| ⊗ Y`. -/
def gCY : Gate 4 :=
  ⟨Mat.of4 [[o, z, z, z], [z, o, z, z], [z, z, z, ⟨0, -1⟩], [z, z, ⟨0, 1⟩, z]], 1⟩

/-- `M† M = c · I`: the gate is unitary after scaling. -/
def Gate.Unitary {n : Nat} (g : Gate n) : Prop :=
  Mat.Eq (Mat.mul (Mat.adj g.mat) g.mat) (Mat.smul (GI.ofInt g.scale) (Mat.id n))

instance {n : Nat} (g : Gate n) : Decidable g.Unitary :=
  inferInstanceAs (Decidable (Mat.Eq _ _))

/-- `M† P M`, which is `c · G† P G`. -/
def Gate.heis {n : Nat} (g : Gate n) (P : Mat n) : Mat n :=
  Mat.mul (Mat.adj g.mat) (Mat.mul P g.mat)

end PPVM.Kernels
