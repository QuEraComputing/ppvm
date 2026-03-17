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

---

## Performance Work (Tasks 11–14)

*(Tasks 11–14 are the performance-focused phase. Task 11 establishes a baseline; Tasks
12–14 each introduce one optimisation and must demonstrate measurable improvement.)*

### Motivation

For systems with many collapse operators or large Pauli sums, the solver is bottlenecked
by three hot paths:

1. **`LindbladOp::apply`** — called 6× per accepted DOPRI5 step. Dominates for large
   `|terms|` and `|p|`. For n=5 qubits with 5 lowering operators and a dense rate matrix,
   `|terms| = 100` and `|p|` can reach several hundred after truncation.
2. **`commutator_real`** — same call count; same structural problems.
3. **HashMap allocation inside `rhs`** — a fresh `PauliSum` (with two `HashMap`s) is
   constructed on every call, discarding capacity that was just earned through growth.

### Benchmark fixture (Task 11)

A Criterion benchmark in `benches/rhs.rs` using the `ByteF64<1, CoefficientThreshold>`
config. The fixture builds:

- **n = 5 qubits**, lowering operators `c_i = X_i + iY_i` for i = 0..4.
- **Dense 5×5 rate matrix** with `γ_ij = 1 / (1 + |i − j|)`.
- **Initial state** `P = Σ_i Z_i` (sum of single-qubit Z operators), threshold `1e-6`.

Two benchmarks:
- `bench_rhs`: a single call to `rhs(None, &lindblad, &p)` where `p` is a snapshot of the
  state after one warm-up solve step (so it is representative, not trivially sparse).
- `bench_solve`: `solve(None, &lindblad, &initial, (0.0, 1.0), &[0.1, 0.2, …, 1.0], …)`
  with default `SolverConfig`.

The fixture and helper to build the Lindblad operator are extracted into a shared
`benches/fixture.rs` module so subsequent benchmark files can reuse them.

Add `criterion` to `[dev-dependencies]` and a `[[bench]]` entry in `Cargo.toml`.

### Task 12 — Loop restructuring

**`commutator_real`**: `left = PhasedPauliWord::from(w_a.clone())` depends only on `w_a`
(fixed for the inner loop). Move it above the inner loop.

**`apply`**: swap loop order so `p.data().iter()` is the outer loop and `&self.terms` is
the inner. Hoist `wa_phased = PhasedPauliWord::from(w_a.clone())` to the outer scope.

*Why it is faster:* in the current order, `p`'s HashMap is traversed `|terms|` times.
After the swap it is traversed once, reducing cache pressure by a factor of `|terms|`.
`self.terms` is a contiguous `Vec` and remains cache-friendly regardless of loop order.

### Task 13 — Collapse anticommutator into one multiplication

**Mathematical basis.** For any two `PhasedPauliWord` values A and B:

```
(A * B).word == (B * A).word       // XOR is commutative on the word bits
```

The phases differ by exactly `2 × comm_parity(A, B) mod 4`, where `comm_parity` is the
parity of the number of single-qubit anti-commuting pairs:

```
comm_parity(A, B) = popcount((A.xbits & B.zbits) XOR (A.zbits & B.xbits)) mod 2
```

From the four possible combinations of `(a_kl.phase & 1, parity)`:

| `a_kl.phase & 1` | parity | `re_phase(t1) + re_phase(t2)` |
|------------------|--------|-------------------------------|
| 0                | 0      | `2 × re_phase(t1.phase)`      |
| 1                | 1      | `2 × re_phase(t1.phase)`      |
| 0                | 1      | 0                             |
| 1                | 0      | 0                             |

So: combined = `2 × re_phase(t1.phase)` when `(a_kl.phase & 1) == parity`, else 0.

**Implementation.** Replace the two-multiplication anticommutator block in `apply` with:
1. Compute `t1 = a_kl.clone() * wa_phased.clone()` (one `MulAssign`).
2. Compute `parity = comm_parity(&term.a_kl.word, &wa_phased.word)` (bitwise, O(N_bytes)).
3. If `(term.a_kl.phase & 1) == parity` and `re_phase(t1.phase) != 0.0`: accumulate
   `(-2 × term.weight × re_phase(t1.phase)) × coeff_a` into `t1.word`.

Add `#[inline] fn comm_parity<A, S>(a: &PauliWord<A, S>, b: &PauliWord<A, S>) -> u8`
using `a.xbits`, `a.zbits` (both `pub`) over the raw byte storage. Keep it in
`lindblad.rs` as a `pub(crate)` helper.

*Why it is faster:* the anticommutator drops from 2 multiplications + up to 2 HashMap
inserts to 1 multiplication + 1 bitwise check + at most 1 insert. For n=5, |terms|=100,
|p|=100: saves ~10 000 multiplications per `rhs` call.

### Task 14 — `SolverCache`: solve-level buffer pre-allocation

**Allocation budget.** Counting `PauliSum` allocations per step in the current code:

| Item                               | Count per step |
|------------------------------------|----------------|
| `y.clone()` for yi stages + y_new  | 6 (stages 2–6 + y_new) |
| `rhs()` internal alloc for k2..k7  | 6              |
| `err_vec` fresh build              | 1              |
| `k1.clone()` in `solve_mut`        | 1              |
| **Total**                          | **14**         |

After Task 14 all of these drop to zero per step. Nine `PauliSum`s are allocated once per
`solve` call (or once per user-managed `SolverCache`). The `estimate_h0` helper still
allocates via `rhs()` but is called only once per solve and is not on the hot path.

**`rhs_into` (prerequisite, added in this task).** Add
`pub(crate) fn rhs_into<T: Config>(ham, lindblad, p, result: &mut PauliSum<T>)` that:
1. Clears the output map: `result.data_mut().clear()` (retains allocated capacity).
2. Calls `commutator_real` and `lindblad.apply` accumulating into `result`.
3. Calls `result.truncate()`.

Add `T::Map: ACMapBase` to the where-clause. Keep the existing `rhs` as a one-line
wrapper (allocates a fresh buffer, calls `rhs_into`) so `estimate_h0` and tests continue
to compile unchanged.

**`SolverCache` layout.** A flat `Vec<PauliSum<T>>` of length 7 holds all k-vectors:
index 0 is the FSAL carry-over (k1), indices 1–6 are stage buffers for k2–k7. One
additional `PauliSum` holds stage states and the 5th-order solution; another holds the
error estimate.

```rust
pub struct SolverCache<T: Config> {
    pub(crate) k:         Vec<PauliSum<T>>,  // length 7; k[0] = FSAL (k1), k[1..=6] = k2..k7
    pub(crate) y_scratch: PauliSum<T>,        // reused for all yi and y_new
    pub(crate) err:       PauliSum<T>,        // error estimate vector
}

impl<T: Config> SolverCache<T> {
    pub fn new(template: &PauliSum<T>) -> Self;  // allocates all 9 PauliSums
}
```

`SolverCache::new` uses `template` only to read `n_qubits` and `strategy`; it does not
clone `template`'s data.

**`HashMap::clone_from` for stage states.** `T::Map: Clone` is satisfied by all standard
configs. `clone_from(&other)` clears the map and re-inserts from `other` without a new
`malloc`, provided the table has sufficient capacity. Replace every `let mut yi = y.clone()`
in `step` with:
```rust
cache.y_scratch.data_mut().clone_from(y.data());
// then add_scaled calls modify cache.y_scratch in place
```

**`std::mem::swap` for zero-copy state update.** After accepting a step, instead of
`*state = y_new` (drops old map, copies new map header):
```rust
std::mem::swap(state, &mut cache.y_scratch);
```
The old state lands in `cache.y_scratch` and is overwritten by `clone_from` at the start
of the next stage — no deallocation needed.

**FSAL without clone.** `cache.k[0]` holds the current k1 and is borrowed immutably
during stage computations. After k7 is written into `cache.k[6]`:
```rust
cache.k.swap(0, 6);  // O(1): k[0] ← k7, k[6] ← stale k1 (overwritten next step)
```

**Seeding `cache.k[0]` at the start of `solve`.** Before entering the step loop, call:
```rust
rhs_into(ham, lindblad, state, &mut cache.k[0]);
```
This replaces the current `let mut k1 = rhs(ham, lindblad, state)`.

**Modified `step` signature** (lives in `dopri5.rs`):
```rust
pub(crate) fn step<T: Config>(
    ham:      Option<&PauliSum<T>>,
    lindblad: &LindbladOp<T>,
    y:        &PauliSum<T>,
    dt:       f64,
    config:   &SolverConfig,
    cache:    &mut SolverCache<T>,
) -> StepResult    // StepResult loses T: simplified to Accept { h_new } | Reject { h_new }
```

`StepResult` in `dopri5.rs` is simplified (the generic parameter and the `y_new`/`k_next`
fields are removed; both Accept and Reject carry only `h_new: f64`).

**User-facing cache.** Expose `SolverCache<T>` publicly. Users calling `solve` many
times with the same system (parameter sweeps, ensemble averages) can allocate the cache
once and pass it to `solve_cached` / `solve_mut_cached`. The existing `solve` / `solve_mut`
remain as backward-compatible wrappers that create a temporary cache internally.
