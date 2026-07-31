# ppvm — Lean formalization

Machine-checked proofs of the mathematics behind the `ppvm` simulator. The Rust
workspace under `crates/` is the source of truth; this Lean development states
the **mathematical spec** and proves that ppvm's bit-level implementation is a
faithful refinement of it. That "spec ⇄ implementation" framing is where the
value is — it turns the packed-bit tricks in `ppvm-pauli-word` / `ppvm-tableau`
into checked theorems and surfaces where code and intended math diverge.

## Why this exists right now: validating the `ppvm-traits-2` algebra

This development is being written **before** the `ppvm-traits-2` crate, to
examine the unifying algebra the redesign is built on. The design docs are

- [`docs/design/traits-2-configuration-and-hashing.md`](../docs/design/traits-2-configuration-and-hashing.md)
  — the graded algebra `C[K]` and the `Sp ⋉ phase` Pauli factorization;
- [`docs/design/word-data-structures.md`](../docs/design/word-data-structures.md);
  and
- [`docs/design/tableau-data-structure.md`](../docs/design/tableau-data-structure.md).

The central design claim is that **Pauli propagation, the stabilizer tableau
mixture, and the generalized tableau are the *same* engine** — the free
`C`-module `C[K]` (`K →₀ C`) — differing only in the key type (`PauliWord`,
`Tableau`, `Bitstring`) and the key's own product. The formalization checks that
this abstraction really is the standard mathematical object it claims to be, and
that the concrete Pauli / bitstring / rotation rules fall out of it.

## Setup

The Lean toolchain is managed with [mise](https://mise.jdx.dev) + `elan`. From
the worktree root:

```bash
mise run setup          # install elan, fetch the Mathlib cache, build PPVM
```

Thereafter:

```bash
mise run lean-build     # lake build (from anywhere in the worktree)
mise run lean-cache     # re-download the prebuilt Mathlib olean cache
cd lean && lake build   # equivalent to lean-build, from the lean/ dir
```

The exact Lean/Mathlib version is pinned by `lean/lean-toolchain` and
`lean/lakefile.toml` (Lean `v4.31.0`, Mathlib `v4.31.0`).

## Layout

| File | Contents | Status |
| :--- | :------- | :----- |
| `PPVM.lean` | Library root; imports every module. | — |
| `PPVM/Basic.lean` | Namespace + Mathlib-resolution smoke test. | — |
| `PPVM/Pauli.lean` | Phase-free single-qubit symplectic group law + symplectic form (`decide`). | ✅ |
| `PPVM/Pauli/Phase.lean` | The **`ℤ/4ℤ` phase cocycle**: `phaseExp` from the Rust `sign`/`imag` booleans, matrix-model reference `phaseRef`, `phaseExp = phaseRef`, the 2-cocycle (associativity) law, the single-qubit **Pauli `Group`** (a real Mathlib `Group` instance) with the **non-split central extension** `toSymplectic : 𝒫₁ →* (ℤ/2)²`, and commutation ↔ symplectic form. | ✅ |
| `PPVM/Pauli/Word.lean` | Lifts the cocycle and commutation law to the **n-qubit** `PauliWord` (per-qubit phase summed mod 4, per `phase/mul.rs`). | ✅ |
| `PPVM/Pauli/Symplectic.lean` | The **symplectic space `GF(2)^{2n}`**: Pauli mult = module `+`; `ω` as a Mathlib `LinearMap.BilinForm (ℤ/2)` proven **alternating** (`IsAlt`); `H`/`S`/`CNOT`/`CZ` proven to be `Sp(2n,2)` **isometries**. | ✅ |
| `PPVM/Pauli/Conjugation.lean` | **Clifford conjugation as signed symplectic automorphisms**: `conjH`/`conjS` are group homs of `𝒫₁` (the `Sp` bit-map + sign `β`), with `HXH=Z`, `SXS†=Y`, `H²=I`. | ✅ |
| `PPVM/Tableau/Frame.lean` | The **stabilizer tableau** (key type #2): the `2n` generators as a **symplectic basis** (`IsSymplecticFrame`); the initial `X`/`Z` frame is one; **every Clifford (`H`/`S`/`CNOT`/`CZ`) preserves it**; measurement coordinate-readout + the deterministic/random **dichotomy**. | ✅ |
| `PPVM/Algebra/GradedMap.lean` | The **abstract graded map** `C[K]`: identifies L0–L4 (`Support`/`Accumulate`/`Scale`/`Pair`/`Multiply`) with `Finsupp` + `AddMonoidAlgebra`; proves `reduce` is structural and `truncate` is not additive. | ✅ |
| `PPVM/Algebra/Twisted.lean` | **`key_mul` is associative**: the twisted convolution `(c,v)·(d,w) = (cd·i^{phaseExp}, v⊕w)` is associative over any ring with `i⁴=1`, from the 2-cocycle. | ✅ |
| `PPVM/Algebra/Truncation.lean` | **Truncation error bounds**: L1 `|error| ≤ Σ|dropped|` (triangle), L2 Cauchy–Schwarz, and the `<` vs `≥` backend cutoff-mismatch witness. | ✅ |
| `PPVM/Algebra/Noise.lean` | **Noise & observables**: unital Pauli-channel eigenvalue `λ_P = 1 − 2Σ_{anti}p_Q`; Pauli-basis orthonormality (`overlap` = δ); zero-state read-out = diagonal-sector sum. | ✅ |
| `PPVM/Instantiations/Bitstring.lean` | The generalized-tableau **`C[Bitstring]`** amplitude algebra: XOR relabel is a bijection; Clifford gates leave amplitudes fixed. | ✅ |
| `PPVM/Instantiations/Rotation.lean` | The non-Clifford **rotation branch**: `iGP` is a distinct Pauli with bits `G ⊕ P`; the `(P,P')` update is a norm-preserving 2-D rotation (`sin²+cos²=1`). | ✅ |

Everything is proved against the pinned Mathlib; `mise run lean-build` is green.

## What is proved (mapped to the design)

**The unifying algebra `C[K]` (`GradedMap.lean`).** The design's "the map is a
graded algebra over `C[K]`" is exactly Mathlib's `Finsupp` (free `C`-module) and
`AddMonoidAlgebra` (the L4 convolution product `single v c * single w d =
single (v+w) (c*d)` — the phase-free `KeyProduct::key_mul`). Two design decisions
are given formal witnesses: `reduce()` is *structural* (a `Finsupp` carries no
zero coordinate to drop), and coefficient-threshold `truncate` is *not additive*
(`trunc 1 + trunc 1 ≠ trunc 2`), so it correctly lives on `Policy`, outside the
algebra. For the *phased* Pauli key, L4 is realized directly: because
`PhasedPauli` is a `Group`, `MonoidAlgebra C 𝒫₁` is an associative unital
`C`-algebra whose basis product is the Pauli group law
(`monoidAlgebra_single_mul`) — Mathlib's group algebra, gotten for free from the
group instance.

**The Pauli phase extension (`Phase.lean`, `Word.lean`).** The headline Tier-1
target. The packed `2·sign + imag` boolean formula in
`crates/ppvm-pauli-word/src/phase/mul.rs` computes exactly the matrix-model
phase exponent, it is a genuine 2-cocycle (so the packed multiplication is
associative — a single-qubit *and* n-qubit theorem), and the phase asymmetry is
`2·ω`, i.e. `P·Q = (−1)^{ω(P,Q)} Q·P`. The single-qubit Pauli group is a genuine
Mathlib `Group`, with the forget-phase map a group homomorphism onto the
symplectic bits — the **non-split** central extension `1 → ℤ₄ → 𝒫₁ → 𝔽₂² → 1`.

**The symplectic structure (`Symplectic.lean`, `Conjugation.lean`).** The
design's "Pauli mod phase is a vector in the symplectic space `GF(2)^{2n}`" is
made literal: the space is the `ℤ/2`-module `(ℤ/2)²ⁿ`, Pauli multiplication is
its addition, and `ω` is a genuine `LinearMap.BilinForm (ℤ/2)` proven alternating
(`IsAlt`) — a Mathlib symplectic space. All four Clifford generators
`H`/`S`/`CNOT`/`CZ` are proven `ω`-isometries (elements of `Sp(2n,2)`) — the
two-qubit ones by a cross-index sum argument. Conjugation itself is realized as a
**signed symplectic automorphism** of the Pauli group `𝒫₁`: `conjH`/`conjS` are
group homomorphisms `G·P·G† = (−1)^{β(G,P)} φ_G(P)`, with `φ_G` the symplectic
bit-map and `β` the explicit sign.

**The other key types (`Tableau/Frame.lean`, `Bitstring.lean`).** The design
unifies three key types under one `Sum` engine; all three are now covered. The
stabilizer **tableau** is proven to satisfy its Aaronson–Gottesman invariant —
the `2n` generators form a **symplectic basis** of `(Sp n, ω)`, the initial
`X`/`Z` frame is one, and every `ω`-isometry (hence every Clifford, via the
`H`/`S` isometries above) preserves it. The generalized tableau's
**`C[Bitstring]`** amplitudes reuse the graded algebra: the XOR branch relabel is
a bijection, and Clifford gates leave the amplitudes fixed. The **rotation**
producer is a reversible, norm-preserving 2-D rotation on the `(P, iGP)`
coefficient plane. (As keys, `PauliWord`, `Frame`, and `Bitstring` are all just
`K` in the same `K →₀ C` — the design's "one engine, N key types.")

## Corrections made to the design doc (from deriving the theorems)

Formalizing forced one substantive fix to
`docs/design/traits-2-configuration-and-hashing.md`:

* The doc described the gate structure as *"the Clifford group is the central
  extension `Sp(2n,2) ⋉ phases`."* Deriving the `Group` instance shows this is
  inaccurate twice over: (1) `𝒫₁` is nonabelian (`not_commutative`) while its
  kernel `ℤ₄` and quotient `(ℤ/2)²` are abelian, so the phase extension is
  **non-split** — a central extension, *not* a semidirect product `⋉`; and
  (2) `Sp(2n,2)` is not a factor of the Pauli group at all — it acts one level
  up, as `Clifford / Pauli`. The doc now states the two exact sequences
  correctly and cites the Lean witnesses.

**Non-Clifford `multiply` / `key_mul` (`Twisted.lean`).** The design's L4
`key_mul` on mod-phase keys — `(c,v)·(d,w) = (c·d·i^{phaseExp(v,w)}, v⊕w)` — is
proven **associative** over any coefficient ring with a fourth root of unity
`i⁴=1` (the `ComplexCoefficient` bound), directly from the phase 2-cocycle
transported through `iᵏ`. So `C[PauliWord]` with `key_mul` is an associative
algebra on the mod-phase key itself (with unit `1·I`).

**Measurement, truncation, noise, observables.** The stabilizer measurement case
split is captured by the **coordinate read-out** (`ω(M, ·)` recovers `M`'s bits
in a symplectic frame) and the **deterministic/random dichotomy** (some
stabilizer anticommutes ⇔ the pivot search fires); measuring `Zq` is
deterministic and `Xq` random on `|0⟩`. Truncation carries an **L1** bound
`|error| ≤ Σ|dropped|` and an **L2** Cauchy–Schwarz bound, and the `<`-vs-`≥`
backend cutoff mismatch has an explicit disagreeing witness. Unital Pauli
channels have eigenvalue `λ_P = 1 − 2·Σ_{Q anticommutes} p_Q`; the Pauli basis is
orthonormal under the `overlap` pairing; and the zero-state read-out is the
diagonal (`X`-free) coefficient sum.

## Divergences (found while scoping, now **fixed**)

The formal spec pinned three genuine bugs in the current Rust; each is fixed on
this branch (`cargo test -p ppvm-sym -p ppvm-pauli-word` green):

* `crates/ppvm-sym/src/mul.rs` — the `Sum × Sum` cross-terms were gated on the
  **signed** `s2.c0 > min_eps` (and `s1.c0`), so a significant *negative* scalar
  dropped the cross terms. Now `.abs() > min_eps`, matching `add_term`.
* `crates/ppvm-sym/src/mul.rs` — the `One × Sum` branch built `new_sum` but never
  wrote it back, silently discarding the product. Now assigns
  `self.inner = Inner::Sum(new_sum)`.
* `crates/ppvm-pauli-word/src/pattern/contains.rs` — `AnyPauliOrIdentity`
  (`[XYZ]?`) matched only `I`; it now matches every Pauli, as `[XYZ]?` intends.

## Status

No open targets remain: every design-doc algebra claim and every scoped
divergence above is either machine-checked in Lean or fixed in Rust, and both
`mise run lean-build` and the affected `cargo test`s are green.
