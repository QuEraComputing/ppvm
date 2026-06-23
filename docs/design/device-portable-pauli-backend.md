<!--
SPDX-FileCopyrightText: 2026 The PPVM Authors
SPDX-License-Identifier: Apache-2.0
-->

# Device-portable `PauliBackend` abstraction

**Status:** draft / RFC — iterate on this branch (`refactor/device-portable-pauli-map`).
**Goal:** make `CuPauliSum = PauliSum<CudaConfig>` possible — a CUDA-backed map as a
drop-in `Config`, reusing the existing `PauliSum<T>` struct, with **no second Rust
type**.

## 1. Why the current abstraction blocks a GPU backend

`PauliSum<T: Config>` is already generic over its storage map (`T::Map`), and one
could imagine a `Config` whose `Map` is a cuco device hashmap. It does not work,
because the `ACMap` trait family (`crates/ppvm-traits/src/traits/map.rs`) is built
around **host (CPU) semantics** in two ways:

1. **Host closures.** Every transformation takes a Rust closure:
   - `ACMapAddAssign::map_add_assign(&self, dest, F: Fn(&W,&V)->(W,V))`
   - `ACMapInsert::map_insert{,_vec,_multiple}(&mut self, dest, F: Fn(&W,&mut V)->Option<…>)`
   - `ACMapScale::scale(F: Fn(&W,&mut V))`, `ACMapRetain::retain(F: FnMut(&W,&V)->bool)`,
     `ACMapContains::contains_with(F: Fn(&V)->bool)`

   The actual gate math lives in these closures in `ppvm-pauli-sum`
   (`sum/clifford.rs`, `sum/rot1.rs`, `sum/rot2.rs`). **You cannot ship an arbitrary
   Rust closure to a CUDA kernel.** This is the deep blocker — deeper than iteration.

2. **Host-reference iteration.** `ACMapIter::Item = (&'a W, &'a V)` and `Trace` folds
   over those host references (used by `overlap`, `trace`, `Display`). A
   device-resident map cannot hand out `&V` into GPU memory, and per-element host
   iteration over device data would erase any GPU advantage (one D2H copy per term).

The bias is therefore not "the trait happens to use a `HashMap`" — it's that the
trait's *operations are expressed as host code*.

## 2. Design principle — reify the operations

Replace open-ended host closures with a **fixed, declarative vocabulary** of
operations that each backend executes however it likes (a CPU loop or a GPU kernel
launch). This is viable because the operation set of Pauli propagation is **small,
closed, and stable**:

- single-qubit Cliffords: `X Y Z H S S† √X √X† √Y √Y†`
- two-qubit Cliffords: `CNOT CZ CY`
- single- and two-qubit Pauli **rotations** (the only *branching* ops: one term →
  up to two)
- scalar `scale_all`, accumulate-`merge`, `truncate`, and `overlap`/`trace`.

The per-term math (`levi_civita`, `comm_2`, the per-gate bit transforms) is **pure
arithmetic on `(x, z, coeff)`** — identical on CPU and GPU. It becomes the single
source of truth: the CPU backend calls the Rust functions; the CUDA backend mirrors
them in C. (Recommended: hoist those `pub fn`s out of `ppvm-pauli-sum` into a shared,
`no_std`-friendly module so both backends share one definition.)

### Key consequence: gate ops don't mention the word type

The reified ops take only **gate identity + qubit indices + scalar coeffs** — no `W`,
no closure, no host references:

```rust
fn apply_clifford_1q(&mut self, gate: Clifford1q, q: usize);
fn apply_rotation_1q(&mut self, axis: Pauli, q: usize, sin: Coeff, cos: Coeff);
```

The backend owns its key representation: the CPU backend stores
`HashMap<PauliWord, Coeff>`; the CUDA backend stores device arrays of packed `u64`
keys. Host interop happens only at the explicit boundary (`export` / import of a
canonical packed form), never in the hot path.

## 3. The trait (this PR)

See `crates/ppvm-traits/src/traits/backend.rs`. Sketch:

```rust
pub trait PauliBackend: Sized {
    type Coeff: Coefficient;

    // container primitives (closure-free)
    fn with_capacity(n_qubits: usize, capacity: usize) -> Self;
    fn n_qubits(&self) -> usize;
    fn len(&self) -> usize;
    fn clear(&mut self);
    fn scale_all(&mut self, factor: Self::Coeff);
    fn merge(&mut self, other: &mut Self);          // accumulate on key collision

    // reified gates (backend supplies a CPU loop or a GPU kernel)
    fn apply_clifford_1q(&mut self, gate: Clifford1q, q: usize);
    fn apply_clifford_2q(&mut self, gate: Clifford2q, a: usize, b: usize);
    fn apply_rotation_1q(&mut self, axis: Pauli, q: usize, sin: Self::Coeff, cos: Self::Coeff);
    fn apply_rotation_2q(&mut self, axis_a: [u8;2], axis_b: [u8;2],
                         a: usize, b: usize, sin: Self::Coeff, cos: Self::Coeff);

    // reductions / truncation / IO (no host-ref iteration)
    fn truncate_abs(&mut self, eps: f64);
    fn overlap(&self, other: &Self) -> Self::Coeff;
    fn export(&self) -> Vec<(PackedWord, Self::Coeff)>;
}
```

## 4. How each backend implements it

| op | CPU (`HashMap<PauliWord, C>`) | CUDA (cuco device map) |
|----|------------------------------|------------------------|
| `apply_clifford_1q` (X/Y/Z) | iterate, flip sign per key — keys unchanged | kernel over coeff array, flip sign by bit test |
| `apply_clifford_1q` (H/S/…) | drain → bit-transform key → re-insert (merge collisions) | kernel rewrites packed keys, then `insert_or_apply` |
| `apply_rotation_1q/2q` | per term: emit ≤2 terms, accumulate | kernel expands ≤2×, then `insert_or_apply` merge |
| `merge` | drain + add (today's `consume`) | `cuco::static_map::insert_or_apply` (already built) |
| `truncate_abs` | `retain(|_,v| v.abs() >= eps)` | stream compaction (drop `\|c\| < eps`) |
| `overlap` | dot over matching keys (today's `trace` path) | the overlap kernel (already built) |

The CUDA `merge` and `overlap` already exist in `ppvm-pauli-sum-cuda` — this PR is the
trait that lets them slot in as a `Config::Map`.

## 5. Migration plan (the async work — pick a row)

1. **Hoist pure bit-math** (`levi_civita`, `comm_2`, per-gate bit transforms) from
   `ppvm-pauli-sum` into a shared module so CPU + CUDA share one definition.
2. **Land `PauliBackend` + gate enums** (this PR) — additive, nothing removed yet.
3. **Implement `PauliBackend` for the CPU maps** (`HashMap`/`IndexMap`/`DashMap`),
   porting the gate closures from `ppvm-pauli-sum` into the impl, selected by enum.
4. **Rewire `PauliSum<T>`** so `Clifford`/`Rotation`/`Trace`/`truncate` dispatch to
   `self.map.<op>()`. The double-buffer (`aux`) folds into the backend's own `merge`,
   simplifying `PauliSum`.
5. **Retire the closure methods** from `ACMap` once unused; keep `ACMapBase` + `export`.
6. **Add `CudaConfig` + the cuco `PauliBackend` impl** in `ppvm-pauli-sum-cuda`, then
   `pub type CuPauliSum = PauliSum<CudaConfig>;`.
7. **Reify truncation strategies** (abs-coeff, weight) so the GPU can run them.

Regression bar: existing `ppvm-pauli-sum` tests must pass **bit-identical** — this is a
structural refactor, not a semantics change.

## 6. Open questions for reviewers

- **Trait name:** new `PauliBackend` vs. reshaping `ACMap` in place? (proposal: new
  trait, deprecate the closure methods.)
- **Coefficient genericity on GPU:** cuco's atomic fast path is ≤ 8-byte payload, so
  `f64` is fine but `Complex<f64>` (16 B) needs the slow path or split re/im. How do we
  keep `Config::Coeff` general while a GPU backend constrains it? (Associated `Coeff`
  with a `GpuCoeff` bound on the CUDA config only?)
- **Where do the shared pure bit-math fns live** — `ppvm-traits` or `ppvm-pauli-word`?
- **Truncation:** full reification vs. a hybrid (reified common strategies + a CPU-only
  `retain` escape hatch behind a capability trait).
- **Keys > 64 qubits:** `PackedWord` here is `(u64 x, u64 z)` (≤ 64 qubits). Wider words
  need a multi-word packed key; the trait should stay agnostic (associated `Key`?).
- **Batched ops** (`CliffordBatch`): one kernel launch for the same gate on many qubits.

## 7. Non-goals

Not changing math/semantics, not optimizing the CUDA kernels here, not removing the old
`ACMap` in this PR (additive first, migrate in follow-ups).
