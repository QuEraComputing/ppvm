# Plan: Heisenberg-Picture Output & Initial-State Overlap

## Physical Context

**The `PauliSum` passed to `solve` is an observable O, not a density matrix.**
It evolves under the **adjoint (Heisenberg-picture) master equation**:

$$\frac{dO}{dt} = i[H, O] + \mathcal{L}^\dagger(O)$$

This is already correctly implemented. In `lindblad.rs::LindbladOp::new` (lines 189–195)
the left operator is conjugated (`phi_k_dag = (4 - sigma_k.phase) % 4`), yielding
`left = L_i†` and `right = L_j`. The `apply` function then computes the sandwich
`L_i† · P · L_j`, which is the adjoint Lindblad form. **The doc comment at `lindblad.rs:539`
is the only thing wrong** — it says `L(P)` but should say `L†(P)`.

To obtain a physically meaningful number from the evolved observable O(t), the user computes
the **expectation value** with respect to a fixed initial state ρ₀:

$$\langle O(t) \rangle = \operatorname{Tr}(\rho_0 \, O(t))$$

ρ₀ is never propagated — it is a static object used only at output time.

A **time correlation function** has the same structure:

$$C(t) = \langle A \cdot B(t) \rangle = \operatorname{Tr}(\rho_0 \, A \cdot B(t))$$

**Defer to a future phase** — requires `PauliSum` operator multiplication (see Future Work).

---

## How Tr(ρ₀ O(t)) Works in the Pauli Basis

For a product state $\rho_0 = \bigotimes_i \rho_0^{(i)}$:

$$\langle O(t) \rangle = \sum_\alpha c_\alpha(t) \prod_{i=1}^{n} \underbrace{\operatorname{Tr}(\rho_0^{(i)} P_{\alpha_i})}_{\text{per-qubit weight}}$$

The per-qubit weight is given directly by the **Bloch vector** $(b^x_i, b^y_i, b^z_i)$
of qubit $i$:

| $P_{\alpha_i}$ | weight |
|---|---|
| I | 1 |
| X | $b^x_i$ |
| Y | $b^y_i$ |
| Z | $b^z_i$ |

**Bitstring states** are the special case $b^x_i = b^y_i = 0$, $b^z_i = (-1)^{b_i}$
(+1 for $|0\rangle$, −1 for $|1\rangle$). Only $\{I,Z\}^{\otimes n}$ strings survive.

**All-zero state** ($b^z_i = +1$ for all $i$): every $\{I,Z\}^{\otimes n}$ term contributes
with weight $+1$, which is exactly what the pattern `"Z?*"` (`Star(SingleOrIdentity(Z))`)
sums — confirming that formula and implementation agree.

The computation is $O(|O(t)| \times n)$ regardless of whether Bloch vectors have nonzero
X/Y components, so supporting general product states costs nothing beyond bitstrings.

---

## Coding Guidelines

This plan is implemented under the rules in `GUIDELINES.md` and `agents/developer.md`.
The most relevant constraints for this work:

- **No changes to `ppvm-runtime`** (GUIDELINES §3). `ProductState` lives in
  `ppvm-timeevolve`, not `ppvm-runtime`.
- **Reuse over new infrastructure** (GUIDELINES §1). The expectation loop uses
  `PauliIter::iter()` from `ppvm-runtime` traits — no new iteration abstractions.
- **No `unwrap()` in production** (developer.md). Bloch-vector length mismatch panics with
  `expect`; all Python-facing functions return `PyResult`.
- **`pub(crate)` for non-public internals** (developer.md). Only `ProductState`,
  `solve_timeevolve_expectation`, and `product_state_expectation` are public.
- **Test before review** (GUIDELINES §5). Each task includes unit tests covering the happy
  path and at least one edge case.
- **Short, focused functions** (developer.md). `ProductState::expectation` is a single
  iterator chain. The Python dispatch in `timeevolve.py` is a short match on keyword args.

---

## Changes Required

### Task A — Fix the misleading doc comment

**File:** `crates/ppvm-timeevolve/src/lindblad.rs:539`

Change `/// Computes \`dP/dt = i[ham, P] + L(P)\`` to:

```
/// Computes the Heisenberg-picture RHS: `dP/dt = i[H, P] + L†(P)`.
///
/// `L†` is the adjoint Lindblad superoperator (observable picture). The `LindbladOp`
/// stores pre-conjugated left operators so `apply` implements the adjoint form directly.
```

---

### Task B — Rename `state`/`initial` → `observable` in the Rust API

All four `solve` variants in `crates/ppvm-timeevolve/src/solve.rs` currently name the
propagated `PauliSum` parameter `state` or `initial`. Rename to `observable` and update
doc comments to say "observable O propagated in the Heisenberg picture."

Also rename `initial_state()` in `crates/ppvm-timeevolve/examples/superradiance.rs` and
update inline comments.

---

### Task C — `ProductState` in `ppvm-timeevolve`

**New file:** `crates/ppvm-timeevolve/src/product_state.rs`

#### Struct and constructors

```rust
/// A separable initial state ρ₀ = ⊗ᵢ ρ₀⁽ⁱ⁾ encoded as per-qubit Bloch vectors.
///
/// Used to compute expectation values ⟨O(t)⟩ = Tr(ρ₀ O(t)) for a Heisenberg-picture
/// observable O(t). ρ₀ is never propagated — it is evaluated only at output checkpoints.
pub struct ProductState {
    /// `bloch[i] = [bx, by, bz]`.
    /// Convention: bz = +1 for |0⟩, bz = -1 for |1⟩.
    bloch: Vec<[f64; 3]>,
}

impl ProductState {
    /// All qubits in |0⟩: bz = +1.
    pub fn all_zero(n_qubits: usize) -> Self { ... }

    /// All qubits in |1⟩: bz = -1.
    pub fn all_one(n_qubits: usize) -> Self { ... }

    /// Computational basis state from a bit slice.
    /// `bits[i] = 0` → |0⟩ (bz=+1);  `bits[i] = 1` → |1⟩ (bz=-1).
    pub fn bitstring(bits: &[u8]) -> Self { ... }

    /// Arbitrary product state. `vectors[i] = [bx, by, bz]`.
    /// Pure states satisfy |b|² = 1; mixed states |b|² < 1.
    /// A warning is printed (via `eprintln!`) if any |bᵢ| > 1 + 1e-9.
    pub fn bloch_vectors(vectors: Vec<[f64; 3]>) -> Self { ... }

    /// Returns the number of qubits this state is defined for.
    pub fn n_qubits(&self) -> usize { self.bloch.len() }

    /// Constructs a `ProductState` from a flat array `[bx₀,by₀,bz₀, bx₁,by₁,bz₁, …]`.
    /// Used by the native Python bridge, which passes Bloch vectors as a flat `Vec<f64>`.
    /// Panics if `flat.len()` is not divisible by 3.
    pub(crate) fn from_flat(flat: &[f64]) -> Self {
        assert!(flat.len() % 3 == 0, "from_flat: length must be divisible by 3");
        let bloch = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
        ProductState { bloch }
    }
}
```

#### `expectation` method

```rust
impl ProductState {
    /// Computes ⟨O⟩ = Tr(ρ₀ O) = Σ_α cα · Πᵢ weight(α_i).
    ///
    /// Runs in O(|O| × n). Called only at save checkpoints — not in the RK hot loop.
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
            let c: f64 = (*coeff).into();
            c * weight
        }).sum()
    }
}
```

`PauliIter::iter()` yields all n qubits (including I positions) as `Pauli` values
from index 0 to n−1. The `enumerate()` recovers the qubit index for the Bloch lookup.

Export `ProductState` from `crates/ppvm-timeevolve/src/lib.rs`.

#### Unit tests (in `product_state.rs`)

- `test_all_zero`: `ProductState::all_zero(2).expectation(&sum([("ZI",1),("IZ",1),("ZZ",1),("II",1)]))` = 4.
- `test_bitstring_10`: bits=[1,0], bz=[-1,+1]; only {I,Z}^⊗2 strings survive, check
  signs (e.g. ZI → weight −1, IZ → weight +1, ZZ → weight −1, II → weight +1).
- `test_bloch_x_plus`: bx=1, by=bz=0; only {I,X}^⊗2 strings survive (IX, XI, XX, II
  all have nonzero weight; Y and Z terms contribute 0).
- `test_expectation_ignores_xy_for_bitstring`: X and Y terms contribute 0 for any bitstring state.
- `test_bloch_warning`: `bloch_vectors` with |b| = 1.5 constructs successfully and returns
  a `ProductState` (eprintln warning is informational only; do not test stderr capture).

---

### Task D — Native Python bridge

**File:** `crates/ppvm-python-native/src/interface_timeevolve.rs`

#### New function: `solve_timeevolve_expectation`

Replaces `solve_timeevolve_observables`. Accepts a flat Bloch-vector array and computes
`ProductState::expectation` inside the `try_arm!` callback — no state cloning, no
heap allocation beyond the single `f64` result per save point.

```rust
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn solve_timeevolve_expectation(
    observable: &Bound<PyAny>,
    bloch_vectors: Vec<f64>,     // flat: [bx₀,by₀,bz₀, bx₁,by₁,bz₁, …], len = 3·n
    lindblad_ops: Vec<(usize, String)>,
    rates: &Bound<PyAny>,
    t_span_start: f64,
    t_span_end: f64,
    save_at: Vec<f64>,
    hamiltonian: Option<&Bound<PyAny>>,
    rtol: f64,
    atol: f64,
    h0: Option<f64>,
    hmin: f64,
    hmax: f64,
) -> PyResult<(Vec<f64>, Vec<f64>)>
```

Construct `ProductState` once before the `try_arm!` dispatch, then capture by reference
in the callback — avoids re-allocating the `Vec<[f64;3]>` at every save point:

```rust
let ps = ProductState::from_flat(&bloch_vectors);
// inside try_arm!, the callback is simply:
|_, p| ps.expectation(p)
```

#### New function: `product_state_expectation`

Standalone expectation for post-solve use on a raw `PauliSum` snapshot:

```rust
#[pyfunction]
pub fn product_state_expectation(
    observable: &Bound<PyAny>,   // a PauliSumIndexMapFxHash$N object
    bloch_vectors: Vec<f64>,     // flat: [bx₀,by₀,bz₀, …], len = 3·n
) -> PyResult<f64>
```

Uses the same `try_arm!` pattern but returns immediately after a single `expectation` call.

Remove `solve_timeevolve_observables` from `lib.rs` registration (and from the module
`#[pymodule]` block).

---

### Task E — Python `ProductState` class

**New file:** `ppvm-python/src/ppvm/product_state.py`

```python
from __future__ import annotations

import warnings
from typing import Sequence

import ppvm_python_native

from .paulisum import PauliSum


class ProductState:
    """A separable initial state ρ₀ = ⊗ᵢ ρ₀⁽ⁱ⁾, encoded as per-qubit Bloch vectors.

    Used to compute expectation values ⟨O(t)⟩ = Tr(ρ₀ O(t)) after Heisenberg-picture
    time evolution.  ρ₀ is never propagated — it is evaluated at output checkpoints only.

    The Bloch vector (bx, by, bz) for qubit i gives:
        Tr(ρ₀⁽ⁱ⁾ I) = 1,  Tr(ρ₀⁽ⁱ⁾ X) = bx,  Tr(ρ₀⁽ⁱ⁾ Y) = by,  Tr(ρ₀⁽ⁱ⁾ Z) = bz.
    Pure states satisfy |b|² = 1; mixed states |b|² < 1.
    """

    def __init__(self, bloch: list[float]) -> None:
        """Low-level constructor.  `bloch` is a flat array [bx₀,by₀,bz₀, bx₁,…]."""
        if len(bloch) % 3 != 0:
            raise ValueError(f"bloch must have length divisible by 3, got {len(bloch)}")
        self._bloch = bloch
        self._n_qubits = len(bloch) // 3

    # ------------------------------------------------------------------ constructors

    @staticmethod
    def all_zero(n_qubits: int) -> "ProductState":
        """All qubits in |0⟩: bz = +1."""
        return ProductState([v for _ in range(n_qubits) for v in (0.0, 0.0, 1.0)])

    @staticmethod
    def all_one(n_qubits: int) -> "ProductState":
        """All qubits in |1⟩: bz = -1."""
        return ProductState([v for _ in range(n_qubits) for v in (0.0, 0.0, -1.0)])

    @staticmethod
    def bitstring(bits: str | Sequence[int]) -> "ProductState":
        """Computational-basis state.

        Args:
            bits: String of '0'/'1' characters or a sequence of 0/1 integers.
                  bits[i] = 0 → |0⟩ (bz=+1);  bits[i] = 1 → |1⟩ (bz=-1).

        Example:
            ProductState.bitstring("0101")   # 4-qubit state |0101⟩
        """
        bloch = []
        for b in bits:
            bit = int(b)
            if bit not in (0, 1):
                raise ValueError(f"bitstring contains invalid character {b!r}")
            bz = 1.0 if bit == 0 else -1.0
            bloch.extend([0.0, 0.0, bz])
        return ProductState(bloch)

    @staticmethod
    def bloch_vectors(
        vectors: Sequence[tuple[float, float, float]],
    ) -> "ProductState":
        """Arbitrary product state via explicit per-qubit Bloch vectors.

        Args:
            vectors: List of (bx, by, bz) tuples, one per qubit.

        Example:
            ProductState.bloch_vectors([(0,0,1), (1,0,0)])  # |0⟩ ⊗ |+⟩
        """
        bloch = []
        for i, (bx, by, bz) in enumerate(vectors):
            norm_sq = bx**2 + by**2 + bz**2
            if norm_sq > 1.0 + 1e-9:
                warnings.warn(
                    f"Bloch vector for qubit {i} has |b|² = {norm_sq:.6g} > 1 "
                    f"(not a valid density matrix). Proceeding anyway.",
                    stacklevel=2,
                )
            bloch.extend([bx, by, bz])
        return ProductState(bloch)

    # ------------------------------------------------------------------ properties

    @property
    def n_qubits(self) -> int:
        return self._n_qubits

    # ------------------------------------------------------------------ expectation

    def expectation(self, observable: PauliSum) -> float:
        """Compute ⟨O⟩ = Tr(ρ₀ O) for a Heisenberg-picture observable O.

        Args:
            observable: A PauliSum representing the (possibly evolved) observable.

        Returns:
            Real-valued expectation value.
        """
        return ppvm_python_native.product_state_expectation(
            observable._interface, self._bloch
        )
```

Export `ProductState` from `ppvm-python/src/ppvm/__init__.py`.

---

### Task F — Update Python `solve` API

**File:** `ppvm-python/src/ppvm/timeevolve.py`

```python
def solve(
    observable: PauliSum,
    lindblad: LindbladOp,
    t_span: tuple[float, float],
    save_at: Sequence[float],
    *,
    hamiltonian: PauliSum | None = None,
    initial_state: ProductState | None = None,
    config: SolverConfig | None = None,
) -> tuple[list[float], list]:
    """Solve the Heisenberg-picture adjoint master equation.

    Propagates an observable O under dO/dt = i[H, O] + L†(O).
    To obtain expectation values, supply an `initial_state` ρ₀; the solver then
    returns ⟨O(t)⟩ = Tr(ρ₀ O(t)) at each save point without cloning the full state.

    Args:
        observable:    Initial observable O (PauliSum) in the Heisenberg picture.
        lindblad:      Dissipation operator (jump ops + rate matrix).
        t_span:        (t_start, t_end) integration interval.
        save_at:       Times at which to record results.
        hamiltonian:   Optional coherent Hamiltonian.
        initial_state: ρ₀ for computing ⟨O(t)⟩ = Tr(ρ₀ O(t)).
        config:        ODE solver parameters.

    Returns:
        (times, results) where:
            - `initial_state` given  → results is list[float]
            - neither given          → results is list[PauliSum] (raw snapshots)
    """
```

**Dispatch logic:**

```python
if initial_state is not None:
    # Fast path: no state clone, single f64 per save point computed in Rust.
    times, results = ppvm_python_native.solve_timeevolve_expectation(
        observable=native_state, bloch_vectors=initial_state._bloch, **kwargs
    )
    return times, results

# Default: raw PauliSum snapshots.
times, raw_states = ppvm_python_native.solve_timeevolve_states(
    observable=native_state, **kwargs
)
return times, [_wrap_native(s) for s in raw_states]
```

Remove the old `observable="trace:<pat>"` dispatch entirely.

---

### Task G — Update examples and tests

**`crates/ppvm-timeevolve/examples/superradiance.rs`:**
- Rename `initial_state()` → `initial_observable()`
- Update inline comments to state: "this is the initial observable O(0); it is propagated
  under the adjoint master equation dO/dt = L†(O)"

**`ppvm-python/test/test_superradiance.py`:**
- Construct `ProductState.all_zero(N)`
- Pass as `initial_state` to `solve`
- Assert returned values are `list[float]` (expectation values), not `PauliSum` objects
- Verify a known value at t=0: `⟨O(0)⟩ = Tr(|0⟩⟨0|^N · Σᵢ Zᵢ) = N`

---

## File Checklist

| File | Action |
|------|--------|
| `crates/ppvm-timeevolve/src/lindblad.rs:539` | Fix doc comment `L(P)` → `L†(P)` |
| `crates/ppvm-timeevolve/src/solve.rs` | Rename `state`/`initial` → `observable`; update docs |
| `crates/ppvm-timeevolve/src/product_state.rs` | **New.** `ProductState` + `expectation` |
| `crates/ppvm-timeevolve/src/lib.rs` | Export `ProductState` |
| `crates/ppvm-timeevolve/examples/superradiance.rs` | Rename, recomment |
| `crates/ppvm-python-native/src/interface_timeevolve.rs` | Add `solve_timeevolve_expectation`, `product_state_expectation`; remove `solve_timeevolve_observables`; rename `state` → `observable` in `solve_timeevolve_states` |
| `crates/ppvm-python-native/src/lib.rs` | Update `#[pymodule]` registration |
| `ppvm-python/src/ppvm/product_state.py` | **New.** Python `ProductState` |
| `ppvm-python/src/ppvm/timeevolve.py` | Rename param, redesign dispatch |
| `ppvm-python/src/ppvm/__init__.py` | Export `ProductState` |
| `ppvm-python/test/test_superradiance.py` | Rewrite with Heisenberg framing |

---

## Future Work (Out of Scope Here)

- **Correlation functions** — fully expressible once PauliSum operator multiplication
  is implemented:
  1. Call `solve` without `initial_state` to obtain `list[PauliSum]` snapshots of B(t).
  2. For each snapshot, compute `AB_t = A * B_t` (PauliSum operator product).
  3. Compute `rho0.expectation(AB_t)` → `Tr(ρ₀ · A · B(t))`.
  The only missing piece is step 2. `ProductState.expectation` already handles step 3
  for any `PauliSum`. This path costs memory (one full clone per save point) but is
  otherwise complete.
- **Custom output functions** (`output: Callable`): if the raw-snapshot + post-processing
  pattern above proves too memory-intensive, a future plan can add this escape hatch.
