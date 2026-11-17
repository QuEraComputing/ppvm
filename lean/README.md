# PPVM Lean formalization

Proofs for `crates/ppvm-pauli-kernels`. Core Lean only (no Mathlib); every
fact is decided in the kernel.

Run these from `lean/`; its `mise.toml` installs elan and defines the tasks
(run `mise trust` there once).

```sh
mise run lean-build    # lake build
mise run lean-tables   # regenerate crates/ppvm-pauli-kernels/tests/generated/tables.rs
```

## What is proved

`PPVM/Kernels/Model.lean` is the ground truth: exact 2×2 and 4×4 matrices over
the Gaussian integers, the Hermitian Pauli matrices (Y = `[[0, −i], [i, 0]]`),
and each gate as an integer matrix `M` with `G = M / √c`.

`PPVM/Kernels/Rules.lean`:

- `g*_unitary`: `M† M = c · I` for all 13 gates, so each `G` is unitary.
- `rule*_correct`: for every input Pauli, the rule's output bits and sign flip
  satisfy `G† P G = (−1)^flip · P'`, computed by matrix multiplication. This
  covers X, Y, Z, H, S, S†, √X, √X†, √Y, √Y†, CX, CZ, CY.
- Identity checks on the matrices: `S² = Z`, `S† S = I`, `√X² = X`,
  `√X† √X = I`, `√Y² = Y`, `√Y† √Y = I`, `H² = I`, and
  `CY = (I ⊗ S) · CX · (I ⊗ S†)`.

`PPVM/Kernels/RowProduct.lean`:

- `site_product`: at one site, `p · q = i^{2·sign + imag} · (p ⊕ q)`, with the
  kernel's `sign` and `imag` masks.
- `two_site_product`: on two sites the exponents add (all 256 cases).
- `popcount_exponent`: for any number of sites,
  `2·popcount(sign) + popcount(imag)` equals the sum of the site exponents.

All theorems depend only on `propext` and `Quot.sound`; none uses `sorry` or
`native_decide`.

## How the Rust code is tied to the proofs

The Lean rules are transcriptions of `src/heisenberg.rs`. The Rust tests do not rely
on that transcription: `scripts/GenTables.lean` builds the conjugation and
product tables from the matrix model by searching over all candidate images,
and `tests/rules.rs` and `tests/product.rs` check the Rust code against those
tables on every `bool` input.

## Gaps

1. **Rust ↔ Lean transcription.** The rule bodies in `Rules.lean` are copied by
   hand from `heisenberg.rs`; nothing links them mechanically. The exhaustive Rust
   tests against the generated tables cover the `bool` case independently.
2. **Generated tables.** `GenTables.lean` runs as compiled code (`lean --run`),
   so the tables rely on the Lean compiler, not the kernel. The file is checked
   in; regenerate it after changing the model.
3. **Integer widths.** That a rule on `u8`…`u128` equals the rule applied to each
   bit position follows from Rust's bitwise operators. It is not formalized;
   `tests/rules.rs` checks `u64` against per-bit `bool` and narrower widths
   against `u64` on random inputs.
4. **Locality.** The proofs are on the one- and two-qubit spaces. That a gate on
   sites `a`, `b` of an `n`-qubit word acts as the rule on those sites and as the
   identity elsewhere (`G ⊗ I`) is not formalized.
5. **Row product for `n` sites.** One site is proved against matrices, two sites
   against the 4×4 tensor product, and the popcount aggregation for any length.
   The general `n`-fold tensor-product statement is not formalized.
6. **Word packing and arithmetic.** That the Rust loop's per-word `count_ones`
   sums to the per-site count, and that wrapping `u32` arithmetic preserves the
   exponent mod 4 (4 divides 2³²), is argued in comments and tested, not proved.
7. **Gate identity.** The matrices are checked for unitarity and the identities
   above. That they are the gates PPVM and Stim name `SQRT_Y`, `CY`, … is a
   convention, fixed by those identities and the conjugation tables.

## Relation to #204

PR #204 has a separate, Mathlib-based development (`lean/PPVM/Pauli/*`). Its
two-qubit conjugations are defined by their bit and sign formulas and pinned
via the group-homomorphism property and textbook images of the generators.
The proofs here instead compute every conjugation from gate matrices, so no
conjugation table is assumed.
