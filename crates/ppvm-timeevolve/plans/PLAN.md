# ppvm-timeevolve: Implementation Plan

## Equation

$$\frac{dP}{dt} = i[H, P] + \mathcal{L}(P)$$

$$\mathcal{L}(P) = \sum_{ij} \frac{\gamma_{ij}}{2} \left( 2 c_i^\dagger P c_j - c_i^\dagger c_j P - P c_i^\dagger c_j \right)$$

$P$ is the observable (`PauliSum`), $H$ the Hamiltonian (`PauliSum`, optional), $c_i$
collapse operators (`PauliSum`), and $\gamma_{ij}$ a real positive-semidefinite rate matrix.

---

## Crate Structure

```
crates/ppvm-timeevolve/
  Cargo.toml
  GUIDELINES.md
  PLAN.md
  src/
    lib.rs       -- public re-exports
    lindblad.rs  -- LindbladOp, commutator_real, rhs
    dopri5.rs    -- DOPRI5 adaptive stepper
    solve.rs     -- solve() / solve_mut()
```

---

## Public Types

### `RateMatrix`

```rust
pub enum RateMatrix {
    Vector(Vec<f64>),      // one rate per collapse operator (diagonal)
    Dense(Vec<Vec<f64>>),  // full NxN real PSD matrix
}

impl From<Vec<f64>> for RateMatrix { ... }  // wraps as Vector
```

### `SolverConfig`

```rust
pub struct SolverConfig {
    pub rtol: f64,          // default: 1e-6
    pub atol: f64,          // default: 1e-9
    pub h0:   Option<f64>,  // None = auto-estimated
    pub hmin: f64,          // default: 1e-12
    pub hmax: f64,          // default: t_span length
}

impl Default for SolverConfig { ... }
```

### `CollapseOp<T: Config>`

Input type for collapse operators. Standard `PauliSum<T>` with `T::Coeff = f64` uses
`PauliWord` as map keys and cannot represent imaginary coefficients (e.g. the `iY` term in
`X + iY`). `CollapseOp` stores terms as `(PhasedPauliWord, f64)` pairs, where any imaginary
unit is carried in the `PhasedPauliWord`'s phase field and the f64 is a real magnitude.

```rust
pub struct CollapseOp<T: Config> {
    terms:    Vec<(PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>, f64)>,
    n_qubits: usize,
}

impl<T: Config> CollapseOp<T> {
    pub fn new(n_qubits: usize) -> Self;
    pub fn push(&mut self, word: PhasedPauliWord<...>, coeff: f64);
}
```

### `LindbladOp<T: Config>`

Preprocessed Lindblad superoperator. Constructed once, reused across solves.

```rust
pub struct LindbladOp<T: Config> { ... }

impl<T: Config> LindbladOp<T> {
    pub fn new(ops: Vec<CollapseOp<T>>, rates: RateMatrix) -> Self;
}
```

---

## Preprocessing (`LindbladOp::new`)

Expand each `CollapseOp` into its `(PhasedPauliWord { word: W_k, phase: φ_k }, r_ik)` terms.

For each pair of terms from `c_i†` and `c_j`:

- `sigma_k†` from `c_i†`: word `W_k`, phase `φ_k† = (4 − φ_k) % 4`, coefficient `r_ik`
- `sigma_l` from `c_j`: word `W_l`, phase `φ_l`, coefficient `r_jl`
- Rate phase: `p = (φ_k† + φ_l) % 4`  (= phase of `conj(φ_k) · φ_l`)
- **Absorb** `p` into `left`: `left = PhasedPauliWord { word: W_k, phase: (φ_k† + p) % 4 }`
- `right = PhasedPauliWord { word: W_l, phase: φ_l }`
- `weight = γ_ij / 2 · r_ik · r_jl`  (real, can be negative)
- Precompute `a_kl = left * right` (single `PhasedPauliWord` multiply)

Absorbing `p` into `left` ensures that `weight · Re(left * W_a * right)` equals the correct
real contribution of this term to `L(P)` — no complex arithmetic needed in `apply`.

Drop terms where `weight == 0`.

**Internal representation:**

```rust
struct LindbladTerm<S, H, W> {
    left:   PhasedPauliWord<S, H, W>,  // sigma_k†, rate phase absorbed
    right:  PhasedPauliWord<S, H, W>,  // sigma_l
    a_kl:   PhasedPauliWord<S, H, W>,  // left * right (precomputed)
    weight: f64,
}
```

---

## RHS Evaluation (`lindblad.rs`)

The RHS function returns a fresh `PauliSum` representing `dP/dt`:

```rust
fn rhs<T: Config>(h: Option<&PauliSum<T>>, l: &LindbladOp<T>, p: &PauliSum<T>) -> PauliSum<T>
```

It accumulates into a fresh zero-initialised `PauliSum` by calling `commutator_real` and
`apply`, then calls `truncate()`.

### `commutator_real`

Adds `i[H, P]` to `result`. `MulAssign<PauliSum>` requires `ComplexCoefficient` and does
not work for f64 configs; this function replaces it.

For each `(W_a, h_a)` in `H` and `(W_b, p_b)` in `P`:
- Multiply: `tmp = PhasedPauliWord::from(W_a.clone()); tmp *= PhasedPauliWord::from(W_b.clone())`
- `tmp.phase ∈ {0, 2}`: commuting — skip (cancels in `[H,P]`)
- `tmp.phase = 1` (`+i`): add `−2 · h_a · p_b` to `tmp.word` in result
- `tmp.phase = 3` (`−i`): add `+2 · h_a · p_b` to `tmp.word` in result

### `LindbladOp::apply` (private)

Adds `L(P)` to `result`.

For each `LindbladTerm { left, right, a_kl, weight }` and `(W_a, coeff_a)` in `P`:

**Sandwich** `2 · left · W_a · right`:
```
let mut tmp = left.clone();
tmp *= PhasedPauliWord::from(W_a.clone());
tmp *= right.clone();
// Re(phase): +1 if tmp.phase==0, -1 if tmp.phase==2, else 0
result += (tmp.word, 2.0 * weight * re_phase(tmp.phase) * coeff_a)  // skip if zero
```

**Anticommutator** `{a_kl, W_a}`:
```
let mut t1 = a_kl.clone(); t1 *= PhasedPauliWord::from(W_a.clone());
result += (t1.word, -weight * re_phase(t1.phase) * coeff_a)

let mut t2 = PhasedPauliWord::from(W_a.clone()); t2 *= a_kl.clone();
result += (t2.word, -weight * re_phase(t2.phase) * coeff_a)
```

---

## DOPRI5 (`dopri5.rs`)

Standard Dormand-Prince 4(5) with FSAL. 6 RHS evaluations per accepted step (k7 = k1 of
next step).

**Error norm:**
$$\text{err} = \frac{\sqrt{\text{overlap}(e,\, e)}}{\text{atol} + \text{rtol} \cdot \sqrt{\text{overlap}(y,\, y)}}$$

where `e = h · Σ_i (b_i − b_i*) · k_i` and `y` is the pre-step state. `PauliSum::overlap`
is the existing Pauli-basis dot product — keys absent from either operand contribute 0,
which is correct.

**Step update:** accept if `err < 1`; then
$$h_\text{new} = h \cdot \text{clamp}(\,0.9 \cdot \text{err}^{-1/5},\; 0.2,\; 10.0\,)$$

**h0 auto-estimate** (Hairer et al.):
1. `d0 = sqrt(overlap(y0, y0))`, `f0 = rhs(y0)`, `d1 = sqrt(overlap(f0, f0))`
2. `h0 = 0.01 · d0 / d1` (or `1e-6` if `d0` or `d1` too small)
3. `f1 = rhs(y0 + h0 · f0)`, `d2 = sqrt(overlap(f1−f0, f1−f0)) / h0`
4. `h1 = (0.01 / max(d1, d2))^(1/5)`
5. `h_init = min(100 · h0, h1, hmax)`

This costs one extra RHS evaluation but is the standard approach.

---

## Solve API (`solve.rs`)

```rust
/// Advances `state` in-place from t_span.0 to t_span.1.
pub fn solve_mut<T, R, F>(
    hamiltonian: Option<&PauliSum<T>>,
    lindblad:    &LindbladOp<T>,
    state:       &mut PauliSum<T>,
    t_span:      (f64, f64),
    save_at:     &[f64],            // sorted, within t_span
    callback:    F,
    config:      SolverConfig,
) -> (Vec<f64>, Vec<R>)
where T: Config, F: Fn(f64, &PauliSum<T>) -> R;

/// Clones `initial` and calls `solve_mut`.
pub fn solve<T, R, F>(
    hamiltonian: Option<&PauliSum<T>>,
    lindblad:    &LindbladOp<T>,
    initial:     &PauliSum<T>,
    t_span:      (f64, f64),
    save_at:     &[f64],
    callback:    F,
    config:      SolverConfig,
) -> (Vec<f64>, Vec<R>)
where T: Config, F: Fn(f64, &PauliSum<T>) -> R;
```

Save points are hit exactly by capping `h` when the next save time would be overshot.

---

## Notes

- **f64-only:** The real-arithmetic approach requires `T::Coeff = f64`. Complex coefficient
  configs are not supported and not needed.
- **Trait bounds:** `commutator_real` and `apply` require
  `PhasedPauliWord<...>: From<T::PauliWordType> + MulAssign + Clone` — the same bounds
  already present on `MulAssign<PauliSum>`. Standard configs satisfy these.
- **No changes to `ppvm-runtime`.**
