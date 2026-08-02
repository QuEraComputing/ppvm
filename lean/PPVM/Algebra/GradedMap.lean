/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Algebra.MonoidAlgebra.Basic
import Mathlib.Algebra.Star.BigOperators
import Mathlib.Algebra.Order.Star.Basic

/-!
# The coefficient map is the free module `C[K]`

The traits-2 design
(`docs/design/traits-2-configuration-and-hashing.md`, "The map is a graded
algebra over `C[K]`") states that the associative coefficient map — the thing
`HashMap`, `Vec`, and `ColumnStore` all implement — is *algebraically* the

> free `C`-module on a set of keys `K`: a finitely-supported function `K ⇀ C`,
> an element of `C[K]`.

In Mathlib that object is `K →₀ C` (`Finsupp`), and when `K` carries a product
it is `MonoidAlgebra C K` / `AddMonoidAlgebra C K`. This file identifies each
graded layer the design names (L0–L4) with the corresponding Mathlib structure,
so the module and algebra axioms the design relies on are inherited rather than
re-proved, and pins down the two places the design deliberately steps *outside*
the algebra:

* `reduce` (drop zero coefficients) is **structural** in `Finsupp` — the ideal
  model has no zeros in support, which is exactly why the design says
  canonicalization "runs only at finalize."
* `truncate` (drop small coefficients) is **not additive**, so it cannot be an
  algebra operation and must live on `Policy` / `Retain`.

## The layer ↔ Mathlib dictionary

| Design layer | Trait | Mathlib |
| --- | --- | --- |
| L0 finite partial function `K ⇀ C` | `Support` | `Finsupp` (`K →₀ C`) |
| L1 abelian-monoid + canonical support | `Accumulate` | `+` and the `Finsupp` support invariant |
| L2 `C`-module action | `Scale` | `Module C (K →₀ C)` (`•`) |
| L3 symmetric bilinear trace pairing | `Pair` | the `Finsupp` pairing |
| L4 group-algebra product | `Multiply` | `AddMonoidAlgebra`/`MonoidAlgebra` |
-/

namespace PPVM.GradedMap

open Finsupp

/-- `C[K]` — the free `C`-module on the key set `K`, i.e. finitely-supported
`K ⇀ C`. This is the mathematical object every `Sum` backend refines. -/
abbrev CMap (K C : Type*) [Zero C] := K →₀ C

variable {K C : Type*}

/-! ### L0 `Support` — the finite partial function -/

section Support
variable [Zero C]

/-- `Support::get` — read-out of a coordinate. -/
def get (f : CMap K C) (k : K) : C := f k

/-- `Support::len` — the size of the canonical (zero-free) support. -/
def len (f : CMap K C) : ℕ := f.support.card

/-- **`get` reads the coordinate function**: the map really is `K ⇀ C`. -/
theorem get_eq (f : CMap K C) (k : K) : get f k = f k := rfl

end Support

/-! ### L1 `Accumulate` — form linear combinations, then canonicalize

`accumulate_batch` merges produced terms into the map; algebraically that is
addition in `K →₀ C`, whose `AddCommMonoid` instance is Mathlib's. The design's
separate `reduce()` step drops keys whose coefficient is zero — but a `Finsupp`
*is* its reduced form by construction. -/

section Accumulate
variable [AddCommMonoid C]

/-- `Accumulate::accumulate_batch` (merge two contributions) is `+` on `C[K]`.
(`noncomputable` only because `Finsupp`'s bundled `+` is; the design's backends
implement it concretely.) -/
noncomputable def accumulate (f g : CMap K C) : CMap K C := f + g

@[simp] theorem accumulate_apply (f g : CMap K C) (k : K) :
    accumulate f g k = f k + g k := by
  simp [accumulate]

/-- Accumulation is commutative and associative — the abelian-monoid formation
the design's L1 requires, inherited from `Finsupp`'s `AddCommMonoid`. -/
theorem accumulate_comm (f g : CMap K C) : accumulate f g = accumulate g f :=
  add_comm f g

theorem accumulate_assoc (f g h : CMap K C) :
    accumulate (accumulate f g) h = accumulate f (accumulate g h) :=
  add_assoc f g h

/-! #### The *batch* is a multiset: order- and partition-invariance

`accumulate_comm`/`accumulate_assoc` are about whole maps, but the operation the
contract exposes is `Accumulate::accumulate_batch(&TermBatch<K, C>)` — a batch of
`(key, coeff)` **terms**. The design deliberately licenses a backend to `gather`
that batch into per-partition sub-batches and run them concurrently (radix
hash-partitioning, `ColumnStore`), i.e. to reorder and cut the batch up. Below is
the spec every such backend must meet, and it is exactly the free-commutative-
monoid structure of accumulation: folding terms in is a monoid hom from
`Multiset (K × C)`, hence (i) order-invariant and (ii) additive under any split
of the batch. The scalar `accumulate(k, c)` sugar is the singleton case. -/

/-- Merge one produced term `(k, c)` into the map — one `accumulate_batch` step
(and the whole of the scalar `accumulate(k, c)` sugar). -/
noncomputable def accumulateTerm (t : K × C) (m : CMap K C) : CMap K C :=
  accumulate m (Finsupp.single t.1 t.2)

/-- Merging terms is **left-commutative**, which is what makes folding a batch
well defined on a `Multiset` (an *unordered* batch) at all. -/
instance : LeftCommutative (accumulateTerm (K := K) (C := C)) where
  left_comm t₁ t₂ m := by
    simp only [accumulateTerm, accumulate]
    exact add_right_comm m _ _

/-- `Accumulate::accumulate_batch` — fold a whole `TermBatch` into the map. The
batch is modelled as a `Multiset`, so the *type* already records that the call
may not depend on the order terms arrive in; the instance above is the proof
obligation that discharges it. -/
noncomputable def accumulateTerms (B : Multiset (K × C)) (m : CMap K C) : CMap K C :=
  Multiset.foldr accumulateTerm m B

/-- The batch viewed as a map in its own right: `∑ (k,c) ∈ B, c·|k⟩`. -/
noncomputable def batchMap (B : Multiset (K × C)) : CMap K C :=
  (B.map fun t => Finsupp.single t.1 t.2).sum

/-- **Folding a batch in is adding the batch's map** — `accumulate_batch` is the
L1 `+` against `batchMap B`, so everything already proved about `accumulate`
(commutativity, associativity) applies at the batch level. -/
theorem accumulateTerms_eq (B : Multiset (K × C)) (m : CMap K C) :
    accumulateTerms B m = accumulate m (batchMap B) := by
  induction B using Multiset.induction with
  | empty => simp [accumulateTerms, batchMap, accumulate]
  | cons a s ih =>
    simp only [accumulateTerms, Multiset.foldr_cons] at *
    simp only [accumulateTerm, accumulate, batchMap, Multiset.map_cons, Multiset.sum_cons, ih]
    exact (add_assoc m _ _).trans (congrArg (m + ·) (add_comm _ _))

/-- **(i) Order-invariance.** Any two orderings of the same batch — e.g. the
`Vec`/slice a backend actually walks, before and after a `gather` — accumulate to
the same map. -/
theorem accumulateTerms_perm {l₁ l₂ : List (K × C)} (h : l₁.Perm l₂) (m : CMap K C) :
    l₁.foldr accumulateTerm m = l₂.foldr accumulateTerm m := by
  have hcoe : (l₁ : Multiset (K × C)) = (l₂ : Multiset (K × C)) := Quot.sound h
  simpa [accumulateTerms, Multiset.coe_foldr] using congrArg (accumulateTerms · m) hcoe

/-- **(ii) Partition-invariance.** Splitting the batch `B = B₁ + B₂` and running
the halves in sequence (or, since `accumulate` is commutative, concurrently into
disjoint partitions and merged) gives the same map as one call on `B`. This is
the algebraic precondition for the design's hash-partitioned / columnar
`accumulate_batch` backends. -/
theorem accumulateTerms_add (B₁ B₂ : Multiset (K × C)) (m : CMap K C) :
    accumulateTerms (B₁ + B₂) m = accumulateTerms B₂ (accumulateTerms B₁ m) := by
  simp only [accumulateTerms_eq, batchMap, accumulate, Multiset.map_add, Multiset.sum_add]
  exact (add_assoc m _ _).symm

/-- **The scalar `accumulate(k, c)` really is a batch of one.** -/
theorem accumulateTerms_singleton (t : K × C) (m : CMap K C) :
    accumulateTerms {t} m = accumulate m (Finsupp.single t.1 t.2) := rfl

/-! #### Producing a batch must **reset** the support before accumulating it

`ApplyProducer::apply_producer` (`ppvm-pauli-sum-2/src/store.rs`) stages every
stored term through a `TermProducer` and then merges the batch back. The design's
`apply` sketch (§"There is no `SumStorage` trait…") writes that merge as
`self.storage.accumulate_batch(&batch)` — straight onto the **live** support. For
a re-keying producer (every Clifford gate, and right-multiplication by a Pauli
word) that is wrong: the batch carries `(φ k, g c)` while the map still carries
`(k, c)`, so the merge leaves *both*. The implementation therefore `reset`s the
support between producing and accumulating; the two theorems below are the
licence for that extra step and the witness that it is not optional.

`produceBatch` models a producer that emits exactly one term `(φ k, g c)` per
stored `(k, c)` — the shape of `RekeyProducer` — and `pushforward` is what the
re-key is *supposed* to compute: Mathlib's `mapDomain ∘ mapRange`. -/

/-- The batch a one-in/one-out re-keying `TermProducer` emits: `(φ k, g (A k))`
for each key in the support. A `Multiset`, since the walk order is not part of
the semantics (`accumulateTerms_perm`). -/
noncomputable def produceBatch (φ : K → K) (g : C → C) (A : CMap K C) : Multiset (K × C) :=
  A.support.val.map fun k => (φ k, g (A k))

/-- The re-keyed map a Clifford conjugation denotes: push the coefficients
forward along `φ` and through the coefficient map `g` (`c ↦ ±c`, so `g 0 = 0`). -/
noncomputable def pushforward (φ : K → K) (g : C → C) (hg : g 0 = 0) (A : CMap K C) :
    CMap K C :=
  Finsupp.mapDomain φ (Finsupp.mapRange g hg A)

/-- **Reset-then-accumulate is the pushforward.** Accumulating the produced batch
into the **empty** map — `StoreAlloc::reset` followed by `accumulate_batch`, which
is what `apply_producer` does — yields exactly the re-keyed map. (No injectivity
is needed here: colliding images accumulate, matching `accumulate_batch`.) -/
theorem pushforward_eq_reset_accumulate (φ : K → K) (g : C → C) (hg : g 0 = 0)
    (A : CMap K C) :
    accumulateTerms (produceBatch φ g A) 0 = pushforward φ g hg A := by
  classical
  rw [accumulateTerms_eq, accumulate, zero_add]
  rw [batchMap, produceBatch, Multiset.map_map]
  rw [pushforward, Finsupp.mapDomain,
    Finsupp.sum_of_support_subset _ (Finsupp.support_mapRange) _
      (fun _ _ => Finsupp.single_zero _)]
  simp only [Finset.sum_eq_multiset_sum, Finsupp.mapRange_apply, Function.comp_def]

/-- With `φ` **injective** — the Clifford case, `Symplectic.…_bijective` — the
pushforward is a genuine re-key: the coefficient of `k` moves to `φ k` untouched
(up to `g`), nothing is doubled and nothing is lost. -/
theorem pushforward_apply (φ : K → K) (hφ : Function.Injective φ) (g : C → C)
    (hg : g 0 = 0) (A : CMap K C) (k : K) :
    pushforward φ g hg A (φ k) = g (A k) := by
  rw [pushforward, Finsupp.mapDomain_apply hφ, Finsupp.mapRange_apply]

/-- **Merging the produced batch onto the un-reset support is not the
pushforward.** Take `A = 1·|true⟩` and the re-key `φ = not` (no coefficient
change): the batch is `{(false, 1)}`, so merging it into `A` leaves `true` *and*
`false`, while the pushforward carries only `false`. Dropping the `reset` would
therefore double-count the entire support on every Clifford gate — silently, and
with no per-term theorem to catch it (this is the `apply`-path analogue of
`Rotation.eagerWalk_ne_twoPass`). -/
theorem merge_without_reset_ne_pushforward :
    accumulateTerms (produceBatch Bool.not (id : ℤ → ℤ) (Finsupp.single true 1))
        (Finsupp.single true 1)
      ≠ pushforward Bool.not (id : ℤ → ℤ) rfl (Finsupp.single true 1) := by
  classical
  intro h
  have hval := DFunLike.congr_fun h true
  rw [accumulateTerms_eq, ← pushforward_eq_reset_accumulate (hg := rfl),
    accumulateTerms_eq, accumulate, accumulate, zero_add] at hval
  simp at hval

/-! #### Sum-plus-sum is coefficient **addition**, not a right-biased overwrite

`ppvm-pauli-sum`'s `impl AddAssign<PauliSum> for PauliSum` routes through
`HashMap::extend`, whose duplicate-key behaviour is `insert` — it **replaces** the
left operand's coefficient with the right one's. That is not the free-module
addition this layer specifies: `accumulate_apply` says the shared coordinate must
be `f k + g k`. The witness below pins the discrepancy on the smallest case, so
the adjudication is machine-checked rather than argued: `1·|k⟩ += 2·|k⟩` is
`3·|k⟩`, not `2·|k⟩`. The new engine's `Accumulate::accumulate_batch` genuinely
accumulates, so it *diverges from old here by design* — old is wrong. -/

/-- **A right-biased overwrite is not module addition.** At a key both operands
carry, `extend`-style replacement returns the right operand's coefficient while
`accumulate` returns their sum; with `1` and `2` those are `2` and `3`. -/
theorem accumulate_ne_overwrite {K : Type*} (k : K) :
    accumulate (Finsupp.single k (1 : ℤ)) (Finsupp.single k 2) k
      ≠ (Finsupp.single k (2 : ℤ)) k := by
  classical simp

/-- **`reduce()` is structural.** A key is in the canonical support *iff* its
coefficient is nonzero: the ideal `C[K]` model carries no zero coordinates to
drop, so canonicalization is meaningful only for a backend (`HashMap`) that
fails to maintain this invariant. This is the formal content of the design's
"`reduce` runs only at finalize." -/
theorem reduce_structural (f : CMap K C) (k : K) :
    f k ≠ 0 ↔ k ∈ f.support :=
  (mem_support_iff).symm

/-- Consequently a coordinate outside the support is exactly zero: there is
nothing for a correct `reduce` to remove from the ideal model. -/
theorem zero_off_support (f : CMap K C) (k : K) (hk : k ∉ f.support) : f k = 0 := by
  by_contra h
  exact hk (mem_support_iff.2 h)

end Accumulate

/-! ### L2 `Scale` — the `C`-module action -/

section Scale
variable [Semiring C]

/-- `Scale::scale` — multiply every coefficient by `s`; the `C`-module action,
Mathlib's `•` on `K →₀ C`. -/
noncomputable def scale (s : C) (f : CMap K C) : CMap K C := s • f

@[simp] theorem scale_apply (s : C) (f : CMap K C) (k : K) :
    scale s f k = s * f k := by
  simp [scale, Finsupp.smul_apply, smul_eq_mul]

/-- Scaling distributes over accumulation — the module distributive law, so L2
is compatible with L1 exactly as a module action must be. -/
theorem scale_accumulate (s : C) (f g : CMap K C) :
    scale s (f + g) = scale s f + scale s g :=
  smul_add s f g

/-- Scaling is associative with `C`-multiplication (`s·(t·f) = (s·t)·f`). -/
theorem scale_scale (s t : C) (f : CMap K C) :
    scale s (scale t f) = scale (s * t) f := by
  simp [scale, mul_smul]

end Scale

/-! ### L3 `Pair` — the bilinear form (observable read-out)

`Pair::overlap` is the `∑ₖ fₖ · gₖ` read side of the hash join, used for
expectation values. We prove it is **additive in each argument** (biadditive) —
the property the batched probe relies on to split work across produced terms —
*and* **`C`-homogeneous in each argument** (`overlap_smul_left`/`_right`), which
together upgrade biadditivity to full `C`-bilinearity. Left-homogeneity holds
over any semiring; right-homogeneity needs commutativity (it moves the scalar
past a coefficient), matching the L2 `Scale` action being over a `CommRing`.

Note this is the **symmetric bilinear** trace pairing `∑ₖ fₖ gₖ = Tr(f g)/2ⁿ`,
with *no* complex conjugation — correct for expectation values of Hermitian
observables (`PPVM.Noise.overlap_single_single` gives the Pauli orthonormality
`Tr(P Q)/2ⁿ = δ`). The `= Tr(f g)/2ⁿ` half of that label is *not* an assumption:
on the Pauli key it is proved against genuine `2ⁿ×2ⁿ` matrices over `ℤ[i]` by
`PPVM.PauliMatrix.overlap_eq_trace_div`, and inside `C[K]` itself `overlap` is
the identity coefficient of the L4 twisted product
(`PPVM.Twisted.twistedConv_apply_id`). A complex state overlap `⟨φ|ψ⟩ = ∑ₖ conj(fₖ) gₖ` is
*sesquilinear* and is a separate operation, not this one. -/

section Pair
variable [Semiring C]

/-- `Pair::overlap` — `⟪f, g⟫ = ∑ₖ fₖ gₖ`. Written as a `Finsupp.sum` over `f`'s
support, since terms with `fₖ = 0` contribute nothing. -/
def overlap (f g : CMap K C) : C := f.sum fun k a => a * g k

/-- **Left-additive** — `⟪f₁ + f₂, g⟫ = ⟪f₁, g⟫ + ⟪f₂, g⟫`. -/
theorem overlap_add_left (f₁ f₂ g : CMap K C) :
    overlap (f₁ + f₂) g = overlap f₁ g + overlap f₂ g :=
  Finsupp.sum_add_index' (fun a => zero_mul (g a)) (fun a b₁ b₂ => add_mul b₁ b₂ (g a))

/-- **Right-additive** — `⟪f, g₁ + g₂⟫ = ⟪f, g₁⟫ + ⟪f, g₂⟫`. Together with
`overlap_add_left` this is biadditivity of the pairing. -/
theorem overlap_add_right (f g₁ g₂ : CMap K C) :
    overlap f (g₁ + g₂) = overlap f g₁ + overlap f g₂ := by
  simp only [overlap, Finsupp.add_apply, mul_add]
  exact Finsupp.sum_add

/-- **Left-homogeneous** — `⟪s · f, g⟫ = s · ⟪f, g⟫`. The scalar factors out of
the first slot; with `overlap_add_left` this is `C`-linearity in the first
argument. Holds over any semiring (no commutativity needed: `s` never has to
move past a coefficient). This is the homogeneity half of the design's
"**bilinear** trace pairing" label, pairing `Scale` (L2) with `Pair` (L3). -/
theorem overlap_smul_left (s : C) (f g : CMap K C) :
    overlap (scale s f) g = s * overlap f g := by
  simp only [overlap, scale]
  rw [Finsupp.sum_smul_index' (fun k => zero_mul (g k)), Finsupp.mul_sum]
  exact Finsupp.sum_congr fun k _ => by rw [smul_eq_mul, mul_assoc]

end Pair

section PairComm
variable [CommSemiring C]

/-- `overlap` rewritten as a plain sum over the union of supports (terms outside
`f`'s support vanish). The symmetric form used to prove commutativity. -/
theorem overlap_eq_union_sum [DecidableEq K] (f g : CMap K C) :
    overlap f g = ∑ k ∈ f.support ∪ g.support, f k * g k :=
  Finsupp.sum_of_support_subset f Finset.subset_union_left _ (fun _ _ => zero_mul _)

/-- **`overlap` is symmetric** — `⟪f, g⟫ = ⟪g, f⟫` over a commutative coefficient
ring. Combined with `overlap_add_left`/`overlap_add_right`, this makes the
design's "**symmetric** bilinear trace pairing" label exact (symmetry is the
commutative half; biadditivity holds over any semiring). -/
theorem overlap_comm (f g : CMap K C) : overlap f g = overlap g f := by
  classical
  rw [overlap_eq_union_sum, overlap_eq_union_sum g f, Finset.union_comm g.support f.support]
  exact Finset.sum_congr rfl fun k _ => mul_comm _ _

/-- **Right-homogeneous** — `⟪f, s · g⟫ = s · ⟪f, g⟫`. Unlike the left slot this
needs commutativity, since factoring `s` out moves it past a coefficient `fₖ`.
Together with `overlap_smul_left`, `overlap_add_left`, and `overlap_add_right`
this is the **full `C`-bilinearity** the design attributes to the pairing (the
homogeneity half that `hermitianOverlap_smul_left`/`_right` already supply for
the sesquilinear twin). -/
theorem overlap_smul_right (s : C) (f g : CMap K C) :
    overlap f (scale s g) = s * overlap f g := by
  simp only [overlap]
  rw [Finsupp.mul_sum]
  exact Finsupp.sum_congr fun k _ => by rw [scale_apply]; ring

end PairComm

/-! ### Clifford conjugation preserves the trace pairing

`ppvm-pauli-sum-2/src/clifford.rs` propagates a Clifford gate `G` over a
`PauliSum` in the Heisenberg picture: it re-keys every term `P ↦ φ_G(P)` by the
phase-stripped symplectic bijection `φ_G` (`Symplectic.hAct_bijective`… ) and
**drains** the conjugation sign `s_P ∈ {+1,−1}` (`Conjugation.conjH_isRealPhase`…)
into the coefficient. Physical correctness of the user-facing `Sum::overlap`
read-out requires this evolution to preserve the Hilbert–Schmidt trace pairing
`⟪A,B⟫ = ∑ₖ aₖ bₖ`.

The two ingredients live in this repo separately — `φ_G` is an `Sp(2n,2)`
bijection (`Symplectic.lean`, `*_bijective`) and every drained sign is a real
`±1` (`Conjugation.lean`, `*_isRealPhase`) — but nothing composed them into the
invariance of the *abstract pairing* under the sign-carrying re-key. That is the
one uncovered link between `Sum::overlap` and the `Sum::{h,s,cnot,cz}` Clifford
path, and it is pure algebra: over a finite key type the Finsupp `overlap` is the
plain sum `∑ₖ aₖ bₖ` (`overlap_eq_fintype_sum`), and the re-key preserves it
because `(s aₖ)(s bₖ) = s² aₖ bₖ = aₖ bₖ` (the sign squares out, `s_P² = 1`) while
the bijection only permutes the summands. -/

section RekeyInvariance
variable [Fintype K] [CommRing C]

/-- Over a finite key type, the Finsupp bilinear `overlap` is the plain sum
`∑ₖ fₖ gₖ` over **all** keys (off-support terms vanish). This identifies the
named `Pair::overlap` with the fintype sum the re-key invariance below is stated
on. -/
theorem overlap_eq_fintype_sum (f g : CMap K C) :
    overlap f g = ∑ k, f k * g k := by
  simp only [overlap, Finsupp.sum]
  refine Finset.sum_subset (Finset.subset_univ _) ?_
  intro k _ hk
  rw [Finsupp.notMem_support_iff.mp hk, zero_mul]

/-- The Clifford Heisenberg re-key on the coefficient function: reindex the keys
by the symplectic bijection `φ` (`= φ_G`) and multiply each coordinate by the
conjugation sign `sgn` (`= s_P`), so `(cliffordRekey φ sgn f) (φ k) = sgn k · fₖ`.
This is exactly the `RekeyProducer` action `k ↦ φ(k)`, `c ↦ c · s_P` of
`clifford.rs`. -/
def cliffordRekey (φ : K ≃ K) (sgn : K → C) (f : K → C) : K → C :=
  fun k => sgn (φ.symm k) * f (φ.symm k)

/-- **Clifford conjugation preserves the trace pairing.** For any key bijection
`φ` and any per-key sign `sgn` with `sgn k · sgn k = 1`, the re-keyed maps have
the same overlap: `⟪conj A, conj B⟫ = ⟪A, B⟫` (both sides written as the fintype
sum of `overlap_eq_fintype_sum`).

This composes the `Sp(2n,2)` **bijectivity** of the bit map
(`Symplectic.lean`, `*_bijective`) with the **reality** `s_P² = 1` of the drained
conjugation sign (`Conjugation.lean`, `*_isRealPhase`) into the semantic
guarantee that ties `Sum::overlap` to the `Sum::{h,s,cnot,cz}` Clifford path in
`ppvm-pauli-sum-2/src/clifford.rs`: Heisenberg evolution is Hilbert–Schmidt
inner-product-preserving. -/
theorem clifford_conjugation_preserves_overlap
    (φ : K ≃ K) (sgn : K → C) (hsgn : ∀ k, sgn k * sgn k = 1) (f g : K → C) :
    (∑ k, cliffordRekey φ sgn f k * cliffordRekey φ sgn g k) = ∑ k, f k * g k := by
  have hterm : ∀ k, cliffordRekey φ sgn f k * cliffordRekey φ sgn g k
      = (sgn (φ.symm k) * sgn (φ.symm k)) * (f (φ.symm k) * g (φ.symm k)) := by
    intro k; simp only [cliffordRekey]; ring
  simp_rw [hterm, hsgn, one_mul]
  exact Equiv.sum_comp φ.symm (fun k => f k * g k)

end RekeyInvariance

/-! ### L3 (sesquilinear) `hermitian_overlap` — the complex state inner product

`overlap` above is the *symmetric bilinear* trace pairing `∑ₖ fₖ gₖ`, correct for
expectation values of Hermitian observables. The physical state/amplitude overlap
`⟨φ|ψ⟩` that a complex `C[Bitstring]` (`GeneralizedTableau`'s amplitudes) needs is
instead **sesquilinear**: `∑ₖ conj(fₖ) gₖ`, conjugate-linear in the first slot. It
requires the coefficient ring to carry a conjugation — the design's `Conjugate`
capability, i.e. a commutative `*`-ring (`StarRing`). Here we model the design's
`Pair::hermitian_overlap` and prove its defining properties: conjugate symmetry,
sesquilinearity (additive, with conjugate- and linear-homogeneity in the two slots), and
(over a `StarOrderedRing`, e.g. `ℂ`) positive semidefiniteness `⟨f,f⟩ ≥ 0`. Over a
real ring `star` is the identity and it collapses back to `overlap`. -/

section HermitianOverlap
variable [Fintype K] [CommRing C] [StarRing C]

/-- `Pair::hermitian_overlap` — the sesquilinear inner product
`⟨f, g⟩ = ∑ₖ conj(fₖ) gₖ`. The `conj` is the `Conjugate` capability (`StarRing`
`star`); over a real ring it is the identity and this collapses to `overlap`. -/
def hermitianOverlap (f g : K → C) : C := ∑ k, star (f k) * g k

/-- **Conjugate symmetry** — `conj ⟨f, g⟩ = ⟨g, f⟩`. -/
theorem hermitianOverlap_conj_symm (f g : K → C) :
    star (hermitianOverlap f g) = hermitianOverlap g f := by
  simp only [hermitianOverlap, star_sum, star_mul', star_star]
  exact Finset.sum_congr rfl fun k _ => mul_comm _ _

/-- **Additive in the second argument.** -/
theorem hermitianOverlap_add_right (f g₁ g₂ : K → C) :
    hermitianOverlap f (g₁ + g₂) = hermitianOverlap f g₁ + hermitianOverlap f g₂ := by
  simp only [hermitianOverlap, Pi.add_apply, mul_add]
  exact Finset.sum_add_distrib

/-- **Additive in the first argument** (with `hermitianOverlap_smul_left` below,
conjugate-linear). -/
theorem hermitianOverlap_add_left (f₁ f₂ g : K → C) :
    hermitianOverlap (f₁ + f₂) g = hermitianOverlap f₁ g + hermitianOverlap f₂ g := by
  simp only [hermitianOverlap, Pi.add_apply, star_add, add_mul]
  exact Finset.sum_add_distrib

/-- **Linear (homogeneous) in the second argument.** -/
theorem hermitianOverlap_smul_right (c : C) (f g : K → C) :
    hermitianOverlap f (c • g) = c * hermitianOverlap f g := by
  simp only [hermitianOverlap, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- **Conjugate-homogeneous in the first argument** — the "sesqui": scaling the
first operand pulls out `conj c`, not `c`. -/
theorem hermitianOverlap_smul_left (c : C) (f g : K → C) :
    hermitianOverlap (c • f) g = star c * hermitianOverlap f g := by
  simp only [hermitianOverlap, Pi.smul_apply, smul_eq_mul, star_mul', Finset.mul_sum]
  exact Finset.sum_congr rfl fun k _ => by ring

end HermitianOverlap

/-- **Positive semidefinite** — over a `StarOrderedRing` (e.g. `ℂ`),
`⟨f, f⟩ = ∑ₖ conj(fₖ) fₖ ≥ 0`, because each term `star (fₖ) * fₖ` is nonnegative.
This is the inner-product positivity a genuine state overlap needs, and it is
exactly why `hermitian_overlap` (not the bilinear `overlap`) is the right pairing
for a complex amplitude vector. -/
theorem hermitianOverlap_self_nonneg {K : Type*} [Fintype K] {C : Type*}
    [CommRing C] [PartialOrder C] [StarRing C] [StarOrderedRing C] (f : K → C) :
    0 ≤ hermitianOverlap f f :=
  Finset.sum_nonneg fun k _ => star_mul_self_nonneg (f k)

/-! ### L4 `Multiply` — the group-algebra product

The design's L4 needs a *`KeyProduct`* on the keys and computes an outer product
over the two operands' supports, accumulating collisions. When the key set is an
additive monoid (the Pauli symplectic bits under `⊕`), that is exactly
`AddMonoidAlgebra` convolution. -/

section Multiply
variable [Semiring C] [AddMonoid K]

/-- `Multiply::multiply_into` — the ring product on `C[K]`, i.e.
`AddMonoidAlgebra` convolution. -/
noncomputable def multiply (f g : AddMonoidAlgebra C K) : AddMonoidAlgebra C K := f * g

/-- **The product of two basis monomials is the design's `key_mul`.**
`(c·|v⟩)·(d·|w⟩) = (c·d)·|v ⊕ w⟩` — one key product `v ⊕ w`, coefficients
multiplied. This is the phase-free core of `KeyProduct::key_mul`; the Pauli case
multiplies the coefficient by the extra `i^{phaseExp}` factor from
`PPVM.PauliPhase` (a 2-cocycle *twist* of this untwisted convolution). -/
theorem multiply_single (v w : K) (c d : C) :
    multiply (Finsupp.single v c) (Finsupp.single w d)
      = Finsupp.single (v + w) (c * d) := by
  simp [multiply, AddMonoidAlgebra.single_mul_single]

end Multiply

/-! ### Truncation is *not* an algebra operation

The design keeps `truncate` off the algebra table and on `Policy`, because
dropping terms still in the support breaks module exactness. Here is the witness:
coefficient-threshold truncation, acting on a single coordinate, is not additive.
`truncMag t` keeps a coefficient only when `|c| ≥ t` — the `CoefficientThreshold`
policy's per-term decision. -/

/-- Coefficient-threshold truncation of one coordinate: keep `c` iff `t ≤ |c|`,
else drop to `0`. This is `CoefficientThreshold::truncate`'s pointwise action. -/
def truncMag (t : ℤ) (c : ℤ) : ℤ := if t ≤ |c| then c else 0

/-- **Truncation breaks additivity.** With threshold `2`, two coefficients that
each round to `0` sum to one that survives: `trunc 1 + trunc 1 = 0 ≠ 2 =
trunc (1 + 1)`. So `truncate` is not a module homomorphism and correctly lives
outside the graded algebra (on `Policy`/`Retain`), not inside it. -/
theorem truncMag_not_additive :
    truncMag 2 1 + truncMag 2 1 ≠ truncMag 2 (1 + 1) := by decide

/-! ### `Retain` — the keep-filter every `Policy` truncates through

`Policy::truncate` is written entirely as `map.retain(|k, c| …)`
(`docs/design/…`, §"`Policy`"), the one non-algebraic map capability. Two shipped
decisions in `ppvm-pauli-sum-2/src/policy.rs` are claims *about* `retain` rather
than about any particular policy, and neither had an oracle:

* `CombinedPolicy::truncate` runs its two members as **two sequential passes**
  (old documents this as a structural difference from PauliPropagation.jl's
  single combined walk). `retain_seq_eq_retain_and` says the sequential form
  computes the **conjunction** of the two keep-rules, and `retain_comm` that the
  surviving key set is a property of the policy pair, not of the pass order — so
  a future backend may fuse or reorder the passes.
* `MaxPauliWeight::truncate` early-returns on the `usize::MAX` sentinel instead
  of running an always-true walk. `retain_of_all_true` / `retain_le_top_eq_self`
  are its soundness: retain-all *is* the identity, so skipping the pass entirely
  is observationally exact.

`retain` is modelled with a `Bool`-valued keep-predicate over `(key, coeff)`,
matching `fn retain(&mut self, keep: impl Fn(&W, &C) -> bool)`. -/

section Retain
variable [Zero C]

/-- `Retain::retain` — keep the coordinate `k` when `keep k (f k)` holds, else
zero it. Coefficients are never modified, only dropped. -/
noncomputable def retain (keep : K → C → Bool) (f : CMap K C) : CMap K C :=
  Finsupp.onFinset f.support (fun k => if keep k (f k) then f k else 0)
    (by
      intro k hk
      by_cases h : keep k (f k) = true
      · rw [if_pos h] at hk; exact Finsupp.mem_support_iff.2 hk
      · rw [if_neg h] at hk; exact absurd rfl hk)

@[simp] theorem retain_apply (keep : K → C → Bool) (f : CMap K C) (k : K) :
    retain keep f k = if keep k (f k) then f k else 0 := rfl

/-- **Two sequential `retain` passes compute the conjunction.** This is
`CombinedPolicy::truncate`'s contract: running `P₁` then `P₂` keeps exactly the
terms both keep-rules accept. (It needs no hypothesis beyond the predicates
depending only on `(k, c)` and `retain` not modifying coefficients — both
enforced by the definition: on a coordinate the first pass zeroed, the second
pass returns `0` whichever way its predicate falls.) -/
theorem retain_seq_eq_retain_and (p q : K → C → Bool) (f : CMap K C) :
    retain p (retain q f) = retain (fun k c => p k c && q k c) f := by
  classical
  ext k
  simp only [retain_apply]
  cases hq : q k (f k) with
  | false => simp [ite_self]
  | true => simp

/-- **…and the two passes commute.** The surviving key set is a property of the
policy *pair*, not of the order `CombinedPolicy` happens to run them in, so a
backend is free to fuse them into one walk or run them the other way round. -/
theorem retain_comm (p q : K → C → Bool) (f : CMap K C) :
    retain p (retain q f) = retain q (retain p f) := by
  rw [retain_seq_eq_retain_and, retain_seq_eq_retain_and]
  congr 1
  funext k c
  exact Bool.and_comm _ _

/-- **Retain-all is the identity.** The equality is at *every* coordinate, not
just on the support, so a backend that stores an explicit zero (which the design
permits between `reduce`s — `Sum` never drops zeros on its own) keeps it too. -/
theorem retain_of_all_true (keep : K → C → Bool) (f : CMap K C)
    (h : ∀ k, keep k (f k) = true) : retain keep f = f := by
  ext k
  rw [retain_apply, if_pos (h k)]

/-- **The `usize::MAX` disable-sentinel is exactly the identity.** A weight cap at
the top of the weight order accepts every key, so `MaxPauliWeight::truncate`'s
`if self.0 == usize::MAX { return; }` early return — skipping the bucket walk
altogether — is observationally equal to running it (`MaxLossWeight` likewise).
Since `retain_of_all_true` is pointwise, this holds including on zero-coefficient
terms. -/
theorem retain_le_top_eq_self {W : Type*} [LinearOrder W] [OrderTop W] (w : K → W)
    (f : CMap K C) : retain (fun k _ => decide (w k ≤ ⊤)) f = f :=
  retain_of_all_true _ _ fun k => by simp

/-- The same statement for any bound that dominates every key's weight, which is
the hypothesis `usize::MAX` discharges concretely. -/
theorem retain_weight_le_eq_self (w : K → ℕ) (n : ℕ) (hn : ∀ k, w k ≤ n)
    (f : CMap K C) : retain (fun k _ => decide (w k ≤ n)) f = f :=
  retain_of_all_true _ _ fun k => by simp [hn k]

end Retain

end PPVM.GradedMap
