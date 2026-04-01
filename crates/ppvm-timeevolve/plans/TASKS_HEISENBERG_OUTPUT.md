# Tasks: Heisenberg-Picture Output & Initial-State Overlap

Before starting, read `PLAN_HEISENBERG_OUTPUT.md` and `GUIDELINES.md` in full.

**Workflow:** Developer implements, hands off to reviewer. Reviewer is the only one who
can mark a task done. Developer commits and moves on only after explicit approval.

---

## Task 1 — Doc and naming corrections (no logic changes)

**Goal:** Fix the one wrong doc comment and rename the propagated-object parameter
throughout the Rust public API to make the Heisenberg-picture framing explicit.

**Steps:**
1. In `crates/ppvm-timeevolve/src/lindblad.rs:539`, replace:
   ```
   /// Computes `dP/dt = i[ham, P] + L(P)` and returns the result.
   ```
   with:
   ```
   /// Computes the Heisenberg-picture RHS: `dP/dt = i[H, P] + L†(P)`.
   ///
   /// `L†` is the adjoint Lindblad superoperator (observable picture). The `LindbladOp`
   /// stores pre-conjugated left operators so `apply` implements the adjoint form directly.
   ```
2. In `crates/ppvm-timeevolve/src/solve.rs`, rename all four public function parameters
   named `state` or `initial` to `observable`. Update doc comments on each to say:
   "observable O propagated in the Heisenberg picture under the adjoint master equation."
3. In `crates/ppvm-timeevolve/examples/superradiance.rs`, rename the `initial_state()`
   function to `initial_observable()` and update its inline comment to:
   "Initial observable O(0). Propagated under dO/dt = L†(O) — NOT the density matrix."

**Verification:**
- `cargo build -p ppvm-timeevolve` passes.
- `cargo clippy -p ppvm-timeevolve -- -D warnings` passes.
- `cargo build --example superradiance` passes.

**Review checklist:**
- [ ] `lindblad.rs` comment accurately describes the adjoint form.
- [ ] All four `solve` variants updated consistently.
- [ ] Doc comments correctly describe the Heisenberg picture.
- [ ] No logic changes anywhere.

---

## Task 2 — `ProductState` struct and constructors

**Goal:** Introduce the `ProductState` type that encodes a separable initial state ρ₀ as
per-qubit Bloch vectors, with all constructors but no expectation logic yet.

**Steps:**
1. Create `crates/ppvm-timeevolve/src/product_state.rs` with:

   ```rust
   pub struct ProductState {
       bloch: Vec<[f64; 3]>,  // bloch[i] = [bx, by, bz]
   }

   impl ProductState {
       pub fn all_zero(n_qubits: usize) -> Self { ... }
       pub fn all_one(n_qubits: usize) -> Self { ... }
       pub fn bitstring(bits: &[u8]) -> Self { ... }  // bits[i] ∈ {0,1}
       pub fn bloch_vectors(vectors: Vec<[f64; 3]>) -> Self { ... }
       pub fn n_qubits(&self) -> usize { self.bloch.len() }
       pub(crate) fn from_flat(flat: &[f64]) -> Self { ... }
   }
   ```

   - `all_zero`: every qubit has `[0.0, 0.0, 1.0]`.
   - `all_one`: every qubit has `[0.0, 0.0, -1.0]`.
   - `bitstring`: bits[i]=0 → `[0,0,1]`; bits[i]=1 → `[0,0,-1]`. Panic via `expect` on
     values other than 0 or 1.
   - `bloch_vectors`: accepts `Vec<[f64;3]>` directly. Print a warning via `eprintln!` if
     any vector has norm² > 1 + 1e-9 (not a valid density matrix, but proceed).
   - `from_flat`: reconstruct from flat `[bx₀,by₀,bz₀, bx₁,…]` using `chunks_exact(3)`.
     Assert (with message) that length is divisible by 3.

2. Declare `mod product_state;` in `crates/ppvm-timeevolve/src/lib.rs`. Do not export yet.

3. Write unit tests in `product_state.rs`:
   - `test_all_zero_n_qubits`: `ProductState::all_zero(3).n_qubits() == 3`.
   - `test_bitstring_encoding`: bits=[0,1] → bloch[0][2]=+1.0, bloch[1][2]=-1.0, all
     X/Y components zero.
   - `test_from_flat_roundtrip`: construct via `bloch_vectors`, convert to flat manually,
     reconstruct via `from_flat`, assert bloch arrays equal.
   - `test_bitstring_invalid`: `bitstring(&[2])` panics.

**Verification:**
- `cargo test -p ppvm-timeevolve` passes (all new tests green).
- `cargo clippy -p ppvm-timeevolve -- -D warnings` passes.

**Review checklist:**
- [ ] All constructors implement the correct Bloch-vector conventions.
- [ ] `from_flat` and `bloch_vectors` are consistent (roundtrip test).
- [ ] `bloch_vectors` warning threshold is `norm² > 1 + 1e-9`.
- [ ] `from_flat` is `pub(crate)`, all others are `pub`.
- [ ] No `unwrap()` in production paths.

---

## Task 3 — `ProductState::expectation` method

**Goal:** Implement the core computation Tr(ρ₀ O) = Σ_α c_α · Π_i weight(α_i) and
export `ProductState` from the crate.

**Steps:**
1. Add `expectation` to `crates/ppvm-timeevolve/src/product_state.rs`:

   ```rust
   pub fn expectation<T>(&self, observable: &PauliSum<T>) -> f64
   where
       T: Config,
       for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
       T::Coeff: Into<f64> + Copy,
       T::PauliWordType: PauliIter,
   {
       observable.data().iter().map(|(word, coeff)| {
           let weight: f64 = word.iter().enumerate().map(|(i, pauli)| {
               match pauli {
                   Pauli::I => 1.0,
                   Pauli::X => self.bloch[i][0],
                   Pauli::Y => self.bloch[i][1],
                   Pauli::Z => self.bloch[i][2],
                   _        => 0.0,  // Pauli::L (loss) does not contribute
               }
           }).product();
           (*coeff).into() * weight
       }).sum()
   }
   ```

2. Export `ProductState` from `crates/ppvm-timeevolve/src/lib.rs`:
   ```rust
   pub use product_state::ProductState;
   ```

3. Add unit tests in `product_state.rs`:
   - `test_all_zero_expectation`: `all_zero(2)` on sum of `ZI(1) + IZ(1) + ZZ(1) + II(1)`
     → 4.0.
   - `test_bitstring_10_expectation`: bits=[1,0], bz=[-1,+1]; check ZI→-1, IZ→+1,
     ZZ→-1, II→+1, total = 0.0.
   - `test_bloch_x_plus_expectation`: bx=1, by=bz=0; {I,X}^⊗2 strings each have weight 1,
     all others 0. Verify on a sum containing XI, IX, XX, II (each coeff 1) → 4.0, and
     that a ZI term contributes 0.
   - `test_xy_zero_for_bitstring`: a PauliSum with only X and Y terms has expectation 0
     for any bitstring state.

**Verification:**
- `cargo test -p ppvm-timeevolve` passes (all tests green).
- `cargo clippy -p ppvm-timeevolve -- -D warnings` passes.
- `cargo doc -p ppvm-timeevolve --no-deps` renders `ProductState` in the public docs.

**Review checklist:**
- [ ] `expectation` matches the formula in `PLAN_HEISENBERG_OUTPUT.md`.
- [ ] `Pauli::L` arm returns 0.0.
- [ ] `enumerate()` indices align with `self.bloch` indices (same qubit order as `PauliIter`).
- [ ] All four test cases pass and cover the correct numerical values.

---

## Task 4 — Native Python bridge: new functions and cleanup

**Goal:** Expose `ProductState::expectation` to Python via two new `#[pyfunction]`s,
remove the old `solve_timeevolve_observables`, and rename the `state` parameter in
`solve_timeevolve_states` for consistency.

**Steps:**
1. In `crates/ppvm-python-native/src/interface_timeevolve.rs`:

   a. Add `solve_timeevolve_expectation(observable, bloch_vectors, ...)`:
      - Signature mirrors `solve_timeevolve_states` but with an added `bloch_vectors: Vec<f64>`
        and return type `PyResult<(Vec<f64>, Vec<f64>)>`.
      - Before the `try_arm!` dispatch, construct:
        ```rust
        let ps = ProductState::from_flat(&bloch_vectors);
        ```
      - Inside `try_arm!`, the callback is `|_, p| ps.expectation(p)`.

   b. Add `product_state_expectation(observable, bloch_vectors: Vec<f64>)`:
      - This is a **static expectation only** — no ODE solve. Construct
        `let ps = ProductState::from_flat(&bloch_vectors);` before the dispatch.
        Use `try_arm!` solely to downcast `observable` to its concrete type, then
        call `ps.expectation(&concrete_observable)` directly and return the `f64`.
        Returns `PyResult<f64>`.

   c. Remove `solve_timeevolve_observables` entirely.

   d. Rename the `state` parameter to `observable` in `solve_timeevolve_states`.

2. In `crates/ppvm-python-native/src/lib.rs`, update the `#[pymodule]` block:
   - Remove `solve_timeevolve_observables`.
   - Register `solve_timeevolve_expectation` and `product_state_expectation`.

**Verification:**
- `cargo build -p ppvm-python-native` passes.
- `cargo clippy -p ppvm-python-native -- -D warnings` passes.

**Review checklist:**
- [ ] `ps` is constructed once before `try_arm!`, not inside the callback.
- [ ] `solve_timeevolve_observables` is fully removed (no dead code, no stray registration).
- [ ] `solve_timeevolve_states` parameter renamed to `observable`.
- [ ] Both new functions registered in the `#[pymodule]` block.
- [ ] `#[allow(clippy::too_many_arguments)]` retained where needed.

---

## Task 5 — Python `ProductState` class

**Goal:** Expose `ProductState` to Python users with the four constructors and the
`expectation` method, plus validation and warnings.

**Steps:**
1. Create `ppvm-python/src/ppvm/product_state.py` with the `ProductState` class as
   specified in `PLAN_HEISENBERG_OUTPUT.md` (Task E).
   - `__init__` validates length divisible by 3.
   - `all_zero`, `all_one`, `bitstring`, `bloch_vectors` static constructors.
   - `bloch_vectors` calls `warnings.warn` if any qubit has norm² > 1 + 1e-9.
   - `expectation` delegates to `ppvm_python_native.product_state_expectation`.
   - `n_qubits` property.

2. Add `from .product_state import ProductState` to `ppvm-python/src/ppvm/__init__.py`.

3. Write tests in `ppvm-python/test/test_product_state.py`:
   - `test_all_zero_n_qubits`: `ProductState.all_zero(4).n_qubits == 4`.
   - `test_bitstring_encoding`: `ProductState.bitstring("01")._bloch == [0,0,1, 0,0,-1]`.
   - `test_invalid_bitstring`: `ProductState.bitstring("2")` raises `ValueError`.
   - `test_bloch_warning`: `ProductState.bloch_vectors([(2,0,0)])` emits a `UserWarning`.
   - `test_expectation_all_zero`: construct a small PauliSum with known coefficients,
     call `rho0.expectation(obs)`, assert correct scalar result.

**Verification:**
- `pytest ppvm-python/test/test_product_state.py` passes (all tests green).

**Review checklist:**
- [ ] `warnings.warn` uses `stacklevel=2` so the warning points to the call site.
- [ ] `bitstring` accepts both `str` and `Sequence[int]`.
- [ ] `expectation` passes `self._bloch` (a flat `list[float]`) to the native function.
- [ ] `ProductState` exported from the package `__init__.py`.

---

## Task 6 — Update Python `solve` API

**Goal:** Rename `state` → `observable`, add `initial_state: ProductState | None`,
remove the old `observable="trace:<pat>"` dispatch.

**Steps:**
1. In `ppvm-python/src/ppvm/timeevolve.py`:
   - Rename first positional parameter `state` → `observable`.
   - Add `initial_state: ProductState | None = None` keyword argument.
   - Remove the old `observable` keyword argument (the one accepting `"trace:<pat>"` strings)
     and its entire dispatch block (the `if observable is None: ... else: ...` branching).
   - Replace dispatch with:
     ```python
     if initial_state is not None:
         times, results = ppvm_python_native.solve_timeevolve_expectation(
             observable=native_observable, bloch_vectors=initial_state._bloch, **kwargs
         )
         return times, results
     times, raw = ppvm_python_native.solve_timeevolve_states(
         observable=native_observable, **kwargs
     )
     return times, [_wrap_native(s) for s in raw]
     ```
   - Update the docstring to describe the Heisenberg picture and the two return modes.
   - Update validation code: rename internal `native_state` → `native_observable`.

**Verification:**
- `pytest ppvm-python/test/` passes, except `test_superradiance.py` which is expected to
  fail here (it calls `solve` with the old API and is rewritten in Task 8). Any other test
  calling `solve(state=..., ...)` with `state` as a keyword will also fail — treat those
  as expected breakage to be fixed in Task 8.

**Review checklist:**
- [ ] Old `observable="trace:<pat>"` keyword arg and its dispatch block are gone, including the `patterns` parsing.
- [ ] `native_observable = observable._interface` (variable renamed throughout).
- [ ] Docstring return section lists both modes correctly.
- [ ] Validation (t_span, save_at, lindblad directions) is unchanged.

---

## Task 7 — Update `test_superradiance.py` and `superradiance.rs` comments

**Goal:** Bring the example and the test in line with the Heisenberg-picture framing,
providing a concrete verified use of `ProductState` and `initial_state=`.

**Steps:**
1. In `ppvm-python/test/test_superradiance.py`:
   - Import `ProductState` from `ppvm`.
   - Construct `rho0 = ProductState.all_zero(N)`.
   - Pass `initial_state=rho0` to `solve`.
   - Include `0.0` as the first element of `save_at` so the t=0 value is recorded.
   - Assert results is `list[float]`.
   - Assert `results[0] == pytest.approx(N)` (known value: Tr(|0⟩^N · Σᵢ Zᵢ) = N at t=0,
     since each Zᵢ has expectation +1 in the all-zero state).

2. In `crates/ppvm-timeevolve/examples/superradiance.rs`, update the top-level comment
   block to state clearly that the PauliSum is an observable (not a state), and that the
   adjoint master equation is being solved.

**Verification:**
- `pytest ppvm-python/test/test_superradiance.py` passes.
- `cargo build --example superradiance` passes.

**Review checklist:**
- [ ] Test asserts the t=0 value `N` to confirm correct initial state overlap.
- [ ] Test result type is `list[float]`, not `list[PauliSum]`.
- [ ] `superradiance.rs` comment does not say "density matrix" anywhere.
