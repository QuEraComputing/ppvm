# Plan: Heisenberg-Picture Output & Initial-State Overlap

## Physical Context

**This is not a density-matrix propagator.** The `PauliSum` passed to `solve` is an
**observable** O, and the solver advances it under the **adjoint (Heisenberg-picture)
master equation**:

$$\frac{dO}{dt} = i[H, O] + \mathcal{L}^\dagger(O)$$

where $\mathcal{L}^\dagger$ is the adjoint Lindblad superoperator. This is already
correctly implemented in `lindblad.rs::apply()`: the construction at lines 189–195
conjugates the left operator (`phi_k_dag = (4 - sigma_k.phase) % 4`), giving `left = L_i†`
and `right = L_j`, so the sandwich `left · P · right = L_i† P L_j` is the adjoint form.
**The code is correct; only the doc comment at `lindblad.rs:539` is wrong** — it says
`L(P)` but should say `L†(P)`.

The observable at time $t$ is $O(t)$, represented as a (truncated) Pauli sum. To obtain
a physically meaningful number, the user computes the **expectation value**:

$$\langle O(t) \rangle = \operatorname{Tr}(\rho_0 \, O(t))$$

where $\rho_0$ is the **initial density matrix** — a fixed state that is never propagated.

A **time correlation function** follows the same pattern:

$$C(t) = \langle A(0) \, B(t) \rangle = \operatorname{Tr}(\rho_0 \, A \cdot B(t))$$

(Correlation functions require `PauliSum` multiplication, which is not yet implemented.
Defer to a future phase.)

---

## How $\operatorname{Tr}(\rho_0 O(t))$ works in the Pauli basis

For a product state $\rho_0 = \bigotimes_i \rho_0^{(i)}$, the expectation value is:

$$\langle O(t) \rangle = \sum_\alpha c_\alpha(t) \prod_{i=1}^{n} \operatorname{Tr}(\rho_0^{(i)} P_{\alpha_i})$$

where each factor is given by the **Bloch vector** $(b^x_i, b^y_i, b^z_i)$ of qubit $i$:

| $P_{\alpha_i}$ | weight |
|---|---|
| I | 1 |
| X | $b^x_i$ |
| Y | $b^y_i$ |
| Z | $b^z_i$ |

**Bitstring states** are the special case where $b^x_i = b^y_i = 0$ and $b^z_i = (-1)^{b_i}$
(+1 for $|0\rangle$, −1 for $|1\rangle$). Only $\{I,Z\}^{\otimes n}$ Pauli strings
survive, each weighted $\prod_{i:\,\alpha_i=Z}(-1)^{b_i}$.

**All-zero state** ($|0\rangle^{\otimes n}$, all $b^z_i = +1$): every $\{I,Z\}^{\otimes n}$
term contributes with coefficient $+1$. The corresponding pattern `"Z?*"` (confirmed
correct — see below) sums these terms directly.

The computation is $O(|O(t)| \times n)$ regardless of whether Bloch vectors have nonzero
X/Y components, so general product states are essentially free to support once bitstrings
are implemented.

### Pattern `"Z?*"` — confirmed correct, not a bug

`"Z?"` parses to `SingleOrIdentity(Z)` (matches I or Z). `*` makes it
`Star(SingleOrIdentity(Z))`. `match_star` in `contains.rs` iterates non-identity positions
and advances while they are Z; any X or Y stops it and the post-loop identity check returns
false. Result: `"Z?*"` matches exactly $\{I,Z\}^{\otimes n}$. Confirmed by
`contains.rs:199` (`"XYY"` does not match `"Z?*"`).

---

## Current State vs. Required State

### What the code does today

| Layer | What `save_at` returns |
|-------|----------------------|
| Rust `solve` callback | Generic `R = F(t, &PauliSum)` — already fully flexible |
| Python `solve(..., observable=None)` | Full `PauliSum` snapshots |
| Python `solve(..., observable="trace:<pat>")` | Coefficient sum over matching Pauli words |

The `"trace:<pattern>"` mode sums coefficients of all matching Pauli strings — this is
**not** the same as $\operatorname{Tr}(\rho_0 O(t))$ for general $\rho_0$.

The parameter `state` in the Python API and doc comments describes the propagated object
as "Initial density-matrix state" — this is **wrong**; it is the observable.

### What we need

1. Fix the doc comment `L(P)` → `L†(P)` in `lindblad.rs:539`.
2. Rename `state` → `observable` everywhere (Python and Rust doc comments).
3. A `ProductState` type (Rust-native) that encodes $\rho_0$ and computes
   $\operatorname{Tr}(\rho_0 O)$ for any `PauliSum` O.
4. Built-in constructors for bitstring states and general product states (Bloch vectors).
5. A `ProductState` exposed to Python, with an `expectation(O_t)` method.
6. Update the Python `solve` API: replace the `observable` kwarg with an
   `initial_state: ProductState | None` parameter; add an `output: Callable | None`
   escape hatch for custom extraction logic.
7. Update examples and test files to use the correct Heisenberg framing.

---

## Changes Required

### 1. Fix doc comment — `lindblad.rs:539`

Change:
```
/// Computes `dP/dt = i[ham, P] + L(P)` and returns the result.
```
to:
```
/// Computes `dP/dt = i[H, P] + L†(P)` (Heisenberg / adjoint picture) and returns the result.
```

**File:** `crates/ppvm-timeevolve/src/lindblad.rs:539`

### 2. Rename `state` → `observable` / update doc comments

The propagated object is an **observable**, not the density matrix.

**Rust files:**
- `crates/ppvm-timeevolve/src/solve.rs` — all four `solve` variants: rename `state`/`initial`
  parameters and update doc comments to say "observable O propagated in the Heisenberg
  picture under the adjoint master equation"
- `crates/ppvm-timeevolve/examples/superradiance.rs` — rename `initial_state()` function,
  update comments

**Python files:**
- `ppvm-python/src/ppvm/timeevolve.py:48` — rename `state` parameter → `observable`
- `ppvm-python/src/ppvm/timeevolve.py:60` — rewrite docstring

### 3. New `ProductState` type — Rust native in `ppvm-runtime`

Location: `crates/ppvm-runtime/src/product_state.rs` (new file)

```rust
/// A product initial state ρ₀ = ⊗ᵢ ρ₀⁽ⁱ⁾, encoded as per-qubit Bloch vectors.
///
/// Used to compute expectation values ⟨O(t)⟩ = Tr(ρ₀ O(t)) for a Heisenberg-picture
/// observable O(t) represented as a PauliSum.
pub struct ProductState {
    /// Per-qubit Bloch vectors (bx, by, bz).  For a pure computational-basis state |b⟩,
    /// bx = by = 0 and bz = (-1)^b (i.e. +1 for |0⟩, -1 for |1⟩).
    bloch: Vec<[f64; 3]>,  // bloch[i] = [bx_i, by_i, bz_i]
}

impl ProductState {
    /// All-zero state |0⟩^⊗n: every qubit has bz = +1.
    pub fn all_zero(n_qubits: usize) -> Self { ... }

    /// All-one state |1⟩^⊗n: every qubit has bz = -1.
    pub fn all_one(n_qubits: usize) -> Self { ... }

    /// Arbitrary computational-basis (bitstring) state.
    /// `bits[i] = 0` → |0⟩ (bz=+1); `bits[i] = 1` → |1⟩ (bz=-1).
    pub fn bitstring(bits: &[u8]) -> Self { ... }

    /// General product state via explicit Bloch vectors.
    /// `vectors[i] = [bx, by, bz]`.  |bx|²+|by|²+|bz|² ≤ 1 for mixed, = 1 for pure.
    pub fn bloch_vectors(vectors: Vec<[f64; 3]>) -> Self { ... }

    /// Computes ⟨O⟩ = Tr(ρ₀ O) = Σ_α c_α · Πᵢ weight(α_i, i).
    pub fn expectation<T: Config>(&self, observable: &PauliSum<T>) -> f64
    where
        for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
        T::Coeff: Into<f64> + Copy,
        T::PauliWordType: PauliWordTrait,
    { ... }
}
```

**`expectation` implementation sketch:**

```rust
pub fn expectation<T: Config>(&self, observable: &PauliSum<T>) -> f64 {
    observable.data().iter().map(|(word, coeff)| {
        let weight: f64 = word.iter().enumerate().map(|(i, pauli)| {
            match pauli {
                Pauli::I => 1.0,
                Pauli::X => self.bloch[i][0],
                Pauli::Y => self.bloch[i][1],
                Pauli::Z => self.bloch[i][2],
            }
        }).product();
        (*coeff).into() * weight
    }).sum()
}
```

This is $O(|O(t)| \times n)$ and handles bitstrings and general product states uniformly.

### 4. Python `ProductState` binding

Location: `ppvm-python/src/ppvm/product_state.py` (new file) + PyO3 native wrapper
in `crates/ppvm-python-native/src/`

```python
class ProductState:
    @staticmethod
    def all_zero(n_qubits: int) -> "ProductState": ...
    @staticmethod
    def all_one(n_qubits: int) -> "ProductState": ...
    @staticmethod
    def bitstring(bits: str | list[int]) -> "ProductState":
        """E.g. bitstring("0101") or bitstring([0,1,0,1]).""" ...
    @staticmethod
    def bloch_vectors(vectors: list[tuple[float, float, float]]) -> "ProductState": ...

    def expectation(self, observable: PauliSum) -> float: ...
```

The `expectation` call delegates to the Rust-native implementation.

### 5. Update Python `solve` API

```python
def solve(
    observable: PauliSum,                        # renamed from `state`
    lindblad: LindbladOp,
    t_span: tuple[float, float],
    save_at: Sequence[float],
    *,
    hamiltonian: PauliSum | None = None,
    initial_state: ProductState | None = None,   # NEW: ρ₀ for ⟨O(t)⟩
    output: Callable[[float, PauliSum], Any] | None = None,  # NEW: custom extractor
    config: SolverConfig | None = None,
) -> tuple[list[float], list]:
```

**Output dispatch (precedence order):**

| `initial_state` | `output` | Return type per save point |
|---|---|---|
| given | `None` | `float`: `initial_state.expectation(O_t)` |
| `None` | given | `Any`: `output(t, O_t)` |
| `None` | `None` | `PauliSum`: raw snapshot |

The old `observable="trace:<pat>"` dispatch is **removed** (no deprecated alias).

Native Rust solve functions in `ppvm-python-native`:
- `solve_timeevolve_states` — unchanged (used for raw snapshots)
- `solve_timeevolve_observables` — **remove** (replaced by `initial_state` path)
- Add: `solve_timeevolve_expectation(state, initial_state_weights, ...)` → `list[float]`
  where `initial_state_weights` encodes the per-qubit Bloch vectors as a flat array

  Alternatively, since `expectation` is cheap (called only at save points), compute it in
  Python by calling the native `ProductState.expectation` on each returned snapshot. This
  avoids a new native solve variant. **Prefer this simpler approach.**

### 6. Update examples and tests

**`crates/ppvm-timeevolve/examples/superradiance.rs`:**
- Rename `initial_state()` → `initial_observable()`
- Add comment explaining Heisenberg picture
- Show how to compute `⟨O(t)⟩` in the callback

**`ppvm-python/test/test_superradiance.py`:**
- Construct `ProductState.all_zero(N)` (or `bitstring("0"*N)`)
- Pass as `initial_state` to `solve`
- Verify returned scalars are expectation values, not raw coefficients

---

## Files to Change (Checklist)

| File | Change |
|------|--------|
| `crates/ppvm-timeevolve/src/lindblad.rs:539` | Fix doc comment `L(P)` → `L†(P)` |
| `crates/ppvm-timeevolve/src/solve.rs` | Rename `state`/`initial` → `observable`; update doc comments |
| `crates/ppvm-timeevolve/examples/superradiance.rs` | Rename, recomment, show expectation value |
| `crates/ppvm-runtime/src/product_state.rs` (new) | `ProductState` struct + `expectation` method |
| `crates/ppvm-runtime/src/lib.rs` | Export `ProductState` |
| `crates/ppvm-python-native/src/` | PyO3 wrapper for `ProductState` and `expectation` |
| `ppvm-python/src/ppvm/product_state.py` (new) | Python `ProductState` class |
| `ppvm-python/src/ppvm/__init__.py` | Export `ProductState` |
| `ppvm-python/src/ppvm/timeevolve.py` | Rename `state`→`observable`; redesign output dispatch |
| `ppvm-python/test/test_superradiance.py` | Rewrite with correct Heisenberg framing |

---

## Open Questions

| # | Question |
|---|----------|
| Q1 | Should the Rust parameter names in `solve.rs` change from `state`/`initial` to `observable`? This is a public API break for Rust callers. Alternatively, keep the Rust name and only update docs + Python. |
| Q2 | For `solve_timeevolve_expectation`: compute in Python (call native `expectation` on each snapshot) or add a dedicated native Rust variant? Python path is simpler; native is faster if save points are very frequent or snapshots are large. |
| Q3 | Should `ProductState::bloch_vectors` validate that each vector satisfies $\|b\| \le 1$? |
| Q4 | Correlation functions $C(t) = \operatorname{Tr}(\rho_0 A \cdot B(t))$ require `PauliSum` multiplication (operator product with phase tracking). Not yet implemented. Defer to a future phase. |
