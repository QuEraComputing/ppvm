# ppvm-timeevolve: Task Breakdown

Before starting any task, read `PLAN.md` and `GUIDELINES.md` in full.

**Workflow:** The developer implements the task and hands it to the reviewer. The reviewer
is the only one who can mark a task as complete. Once the reviewer explicitly approves,
the developer creates a commit and then moves on to the next task. No commit is created
and no next task is started until approval has been given.

---

## Task 0 — Crate Scaffold

**Goal:** Establish a compiling, empty crate wired into the workspace.

**Steps:**
1. Add `crates/ppvm-timeevolve` to the workspace `members` list in the root `Cargo.toml`.
2. Create `crates/ppvm-timeevolve/Cargo.toml` with a path dependency on `ppvm-runtime`
   (default features).
3. Create empty `src/lib.rs`, `src/lindblad.rs`, `src/dopri5.rs`, `src/solve.rs` and
   declare all four modules in `lib.rs`.

**Verification:**
- `cargo build -p ppvm-timeevolve` passes with zero errors and zero warnings.
- `cargo test -p ppvm-timeevolve` compiles and runs (no tests yet).

**Review checklist:**
- [ ] Workspace `Cargo.toml` updated correctly.
- [ ] Module structure matches `PLAN.md`.
- [ ] No logic written yet.

---

## Task 1 — `RateMatrix` and `SolverConfig`

**Goal:** Define the two public configuration types with no logic beyond construction and
default values.

**Steps:**
1. In `src/lindblad.rs`, define and export `RateMatrix`:
   ```rust
   pub enum RateMatrix {
       Vector(Vec<f64>),
       Dense(Vec<Vec<f64>>),
   }
   impl From<Vec<f64>> for RateMatrix { ... }
   ```
2. In `src/solve.rs`, define and export `SolverConfig`:
   ```rust
   pub struct SolverConfig {
       pub rtol: f64,        // default 1e-6
       pub atol: f64,        // default 1e-9
       pub h0:   Option<f64>,// default None (auto-estimated)
       pub hmin: f64,        // default 1e-12
       pub hmax: f64,        // default f64::INFINITY
   }
   impl Default for SolverConfig { ... }
   ```
3. Re-export both from `lib.rs`.

**Unit tests:**
- `RateMatrix::from(vec![1.0, 2.0])` produces `RateMatrix::Vector(v)` where `v == [1.0, 2.0]`.
- `SolverConfig::default()` has the exact default values listed above.
- A `Dense` rate matrix can be constructed manually and its rows accessed.

**Review checklist:**
- [ ] Types match `PLAN.md` exactly.
- [ ] `From<Vec<f64>>` is implemented and tested.
- [ ] No solver logic present.

---

## Task 2 — `CollapseOp`

**Goal:** Define the input type for collapse operators.

**Why not `PauliSum`?** Standard `PauliSum<T>` with `T::Coeff = f64` stores `(PauliWord,
f64)` pairs and cannot represent imaginary coefficients. The `iY` in `X + iY` must be
encoded as a `PhasedPauliWord` with `phase = 1` (+i) and real coefficient `1.0`. This
requires a dedicated type as described in `PLAN.md`.

**Steps:**
1. In `src/lindblad.rs`, define and export:
   ```rust
   pub struct CollapseOp<T: Config> {
       pub(crate) terms: Vec<(PhasedPauliWord<...>, f64)>,
       pub(crate) n_qubits: usize,
   }
   impl<T: Config> CollapseOp<T> {
       pub fn new(n_qubits: usize) -> Self;
       pub fn push(&mut self, word: PhasedPauliWord<...>, coeff: f64);
   }
   ```
2. Add a stub `LindbladOp` struct (empty, just the type) and a stub `LindbladOp::new`
   signature accepting `Vec<CollapseOp<T>>`. No logic yet.

**Unit tests:**
- Construct `X + iY` as a `CollapseOp` with two `push` calls:
  `(PhasedPauliWord(X, phase=0), 1.0)` and `(PhasedPauliWord(Y, phase=1), 1.0)`.
  Assert `terms.len() == 2`.
- Construct `Z` (real, phase=0) as a `CollapseOp` with one push. Assert `terms.len() == 1`.
- `n_qubits` is stored correctly.

**Review checklist:**
- [ ] Phase field encodes the imaginary unit, not the f64 coefficient.
- [ ] Stub `LindbladOp::new` accepts `Vec<CollapseOp<T>>` (not `Vec<PauliSum<T>>`).
- [ ] Type is consistent with `PLAN.md`.

---

## Task 3 — `LindbladOp` Preprocessing

**Goal:** Implement `LindbladOp::new` — expand `CollapseOp` inputs into `LindbladTerm`
entries with absorbed rate phase, as specified in `PLAN.md`.

**Steps:**
1. Define the private `LindbladTerm` struct (see `PLAN.md`).
2. Implement `LindbladOp::new`:
   - For each pair of terms from `c_i†` and `c_j` across all (i, j):
     - Compute `φ_k† = (4 − φ_k) % 4`.
     - Compute combined rate phase `p = (φ_k† + φ_l) % 4`.
     - Set `left.phase = (φ_k† + p) % 4`; set `right = sigma_l`.
     - `weight = γ_ij / 2 · r_ik · r_jl`.
     - Precompute `a_kl = left * right`.
   - Drop terms where `weight == 0.0`.

**Unit tests** (use a 1-qubit config throughout):

- **Single real op, vector rate:**
  `c = PhasedPauliWord(Z, phase=0)` with coeff `1.0`, `rates = [1.0]`.
  Expect 1 term: `left = (Z, phase=0)`, `right = (Z, phase=0)`, `a_kl = (I, phase=0)`,
  `weight = 0.5`.

- **Single imaginary op:**
  `c = PhasedPauliWord(Y, phase=1)` (`= iY`) with coeff `1.0`, `rates = [1.0]`.
  Compute expected phases manually and assert they match.
  (φ_k=1, φ_k†=3, φ_l=1, p=(3+1)%4=0, left.phase=(3+0)%4=3, weight=0.5.)

- **Two-term op:**
  `c = X + iY` (two terms). Expect 4 `LindbladTerm`s (2×2 pairs). Verify the total
  number of terms and that none have `weight == 0`.

- **Dense rate matrix with off-diagonal entries:**
  Two ops `c1 = X`, `c2 = Y`, `γ = [[1.0, 0.5], [0.5, 1.0]]`.
  Expect `2×2×2×2 = 4` cross-terms total (only single-term ops here: 1×1 per (i,j) pair,
  so 4 terms for the 4 (i,j) pairs). Verify the off-diagonal (i=0,j=1) term has
  `weight = 0.5/2 * 1.0 * 1.0 = 0.25`.

**Review checklist:**
- [ ] Phase conjugation formula is `(4 − φ_k) % 4`.
- [ ] Absorbed phase in `left` is `(φ_k† + p) % 4`.
- [ ] `a_kl` is `left * right` (with absorbed phase), not `sigma_k† * sigma_l` bare.
- [ ] Zero-weight terms are dropped.
- [ ] All four test cases pass and match manual derivations.

---

## Task 4 — `commutator_real`

**Goal:** Implement `i[H, P]` using real arithmetic, accumulating into an existing
`PauliSum`. This replaces `MulAssign<PauliSum>`, which requires `ComplexCoefficient` and
does not compile for f64 configs.

**Steps:**
1. In `src/lindblad.rs`, implement:
   ```rust
   pub(crate) fn commutator_real<T: Config>(
       ham:    &PauliSum<T>,
       p:      &PauliSum<T>,
       result: &mut PauliSum<T>,
   )
   ```
   For each `(W_a, h_a)` in `ham` and `(W_b, p_b)` in `p`: compute the `PhasedPauliWord`
   product `tmp = W_a * W_b`. If `tmp.phase == 1`, add `−2·h_a·p_b` to `tmp.word`; if
   `tmp.phase == 3`, add `+2·h_a·p_b`; otherwise skip.

**Unit tests** (1-qubit, single-term H and P):
- `i[X, X] = 0` → result is empty.
- `i[Z, X] = −2Y`: ZX has phase 1 (+i), so contribute `−2·1·1 = −2` to Y.
- `i[X, Z] = +2Y`: XZ has phase 3 (−i), so contribute `+2` to Y.
- `i[Z, Y] = +2X`: ZY has phase 3 (−i), so contribute `+2` to X.
- Multi-term: `H = 0.5·Z`, `P = X + Y`. Verify linearity holds against the above.
- Accumulation: call twice on same `result`; verify values double.

**Review checklist:**
- [ ] Phase 0 and 2 (commuting) are skipped.
- [ ] Sign for phase 1 is `−2`, for phase 3 is `+2`.
- [ ] Function accumulates into `result` without clearing it first.
- [ ] Trait bounds compile for the standard `ByteF64` config.

---

## Task 5 — `LindbladOp::apply`

**Goal:** Implement the Lindblad superoperator `L(P)`, accumulating into an existing
`PauliSum`.

**Steps:**
1. Define the helper `fn re_phase(phase: u8) -> f64` returning `1.0` (phase=0), `−1.0`
   (phase=2), or `0.0` (phase=1,3).
2. Add `pub(crate) fn apply(&self, p: &PauliSum<T>, result: &mut PauliSum<T>)` to
   `LindbladOp`. Implement the sandwich and anticommutator logic from `PLAN.md`. Skip any
   contribution where `re_phase(phase) == 0.0`.

**Unit tests** (`c = X`, `γ = [1.0]`, 1-qubit):
- `L(X) = 0`: X is invariant under X-dephasing (`2XXX − 2X = 2X − 2X`).
- `L(Z) = −4·Z`: `2·X·Z·X = 2·(−Z)`, so `−2Z − 2Z = −4Z`.
  Derivation: XZ = −iY (phase=3), then (−iY)·X = −i·YX = −i·(−iZ) = −Z. So XZX = −Z. ✓
- `L(Y) = −4·Y`: analogous (XYX = −Y).
- Accumulation: call twice; verify values double.

**Unit test — lowering operator** (`c = X + iY`, `γ = [1.0]`, 1-qubit):
- Construct `c` as a `CollapseOp` with terms `(PhasedPauliWord(X, phase=0), 1.0)` and
  `(PhasedPauliWord(Y, phase=1), 1.0)`, build `LindbladOp`, and apply to `P = Z`.
- Derivation (write this as a comment in the test):
  - `c† = X − iY`, `c†c = (X−iY)(X+iY) = XX + XiY − iYX − iY·iY = I + iXY − iYX + YY`
    `= I + i(iZ) − i(−iZ) + I = I − Z − Z + I = 2I − 2Z`.
  - Sandwich: `(X−iY)·Z·(X+iY)`.
    Step 1: `ZX = iY` (phase 1), `ZY = −iX` (phase 3).
    So `Z·(X+iY) = iY + i(−iX) = iY + X`.
    Step 2: `(X−iY)·(X + iY) = XX − iYX + iXY − i²YY = I + Z + Z + I = 2I + 2Z`.
    Wait — recompute with left acting on the right:
    `(X−iY)(iY+X)`: expand:
      `X·iY = i·XY = i·iZ = −Z`
      `X·X   = I`
      `(−iY)·iY = −i²·YY = I`
      `(−iY)·X  = −i·YX = −i·(−iZ) = −Z`
    Sum: `−Z + I + I − Z = 2I − 2Z`.
  - `L(Z) = 2·(2I−2Z) − (2I−2Z)·Z − Z·(2I−2Z)`
    `= 4I − 4Z − (2Z − 2I) − (2Z − 2I)`
    `= 4I − 4Z − 2Z + 2I − 2Z + 2I`
    `= 8I − 8Z`.
- Assert `apply` produces `{I: 8.0, Z: −8.0}`.

**Review checklist:**
- [ ] Sandwich uses `left`, then `W_a`, then `right` in that order.
- [ ] Anticommutator is subtracted (negative sign).
- [ ] `re_phase` is correct for all four phase values.
- [ ] Function accumulates without clearing `result`.

---

## Task 6 — `rhs` Function

**Goal:** Compose the full right-hand side `dP/dt = i[H,P] + L(P)` with truncation.

**Steps:**
1. In `src/lindblad.rs`, implement:
   ```rust
   pub(crate) fn rhs<T: Config>(
       ham:     Option<&PauliSum<T>>,
       lindblad: &LindbladOp<T>,
       p:       &PauliSum<T>,
   ) -> PauliSum<T>
   ```
   Create a fresh empty `PauliSum` with the same `n_qubits` and `strategy` as `p`
   (consult the `ppvm-runtime` builder API). Call `commutator_real` if `ham.is_some()`,
   then `lindblad.apply`, then `result.truncate()`. Return the result.

**Unit tests:**
- **Pure Hamiltonian:** `H = 0.5·Z`, no Lindblad terms, `P = X`.
  `i[0.5Z, X] = 0.5·(−2Y) = −Y`. Assert result is `{Y: −1.0}`.
- **Pure Lindblad:** `c = X`, `γ = [1.0]`, `P = Z`.
  Assert result is `{Z: −4.0}`.
- **Both:** `H = 0.5·Z`, `c = X`, `γ = [1.0]`, `P = X`.
  Hamiltonian gives `−Y`; Lindblad gives 0. Assert result is `{Y: −1.0}`.
- **Truncation:** use a `CoefficientThreshold` strategy with a high threshold; verify
  that a small term in the result is removed.

**Review checklist:**
- [ ] Fresh `PauliSum` is correctly initialised (not sharing map state with `p`).
- [ ] `truncate()` is called on the *result*, not on `p`.
- [ ] `ham = None` is handled without panic.

---

## Task 7 — DOPRI5 Single Step

**Goal:** Implement one adaptive Dormand-Prince 4(5) step.

**Steps:**
1. In `src/dopri5.rs`, define the Butcher tableau as constants (cross-check against
   Hairer et al. "Solving ODEs I", Table 5.2).
2. Implement:
   ```rust
   pub(crate) enum StepResult<T: Config> {
       Accept { y_new: PauliSum<T>, k_next: PauliSum<T>, h_new: f64 },
       Reject { h_new: f64 },
   }

   pub(crate) fn step<T: Config>(
       ham:      Option<&PauliSum<T>>,
       lindblad: &LindbladOp<T>,
       y:        &PauliSum<T>,
       k1:       PauliSum<T>,      // FSAL: already-computed rhs(y)
       dt:       f64,              // current step size
       config:   &SolverConfig,
   ) -> StepResult<T>
   ```
   Compute k2…k6, build `y_new` (5th-order) and error estimate `e = dt · Σ (b_i−b_i*)·k_i`,
   evaluate the error norm, accept or reject, and return the new step size.

**Unit tests:**
- **Zero RHS:** `ham = None`, empty `LindbladOp`, any `P`. `y_new` must equal `y`
  and step must always be accepted.
- **Larmor precession, small step:**
  `H = 0.5·Z`, no Lindblad, `P(0) = X`, `dt = 0.01`.
  Exact solution: `P(t) = cos(t)·X − sin(t)·Y`.
  Assert the returned `y_new` matches the first-order Taylor expansion
  `X − 0.01·Y` to within `O(dt²) ≈ 1e-4`.
- **Step rejection:** force `dt` large enough to exceed tolerances; assert `Reject` is
  returned and `h_new < dt`.

**Review checklist:**
- [ ] Butcher tableau coefficients cited against a reference and verified.
- [ ] Error norm uses `PauliSum::overlap` as in `PLAN.md`.
- [ ] FSAL: `k_next` in `Accept` is `rhs(y_new)` — not recomputed by the caller.
- [ ] Step size update formula matches `PLAN.md` exactly.
- [ ] No parameter name collision (`ham` for Hamiltonian, `dt` for step size).

---

## Task 8 — h0 Auto-Estimation

**Goal:** Implement the initial step size heuristic from `PLAN.md` (Hairer et al.).

**Steps:**
1. In `src/dopri5.rs`, implement:
   ```rust
   pub(crate) fn estimate_h0<T: Config>(
       ham:      Option<&PauliSum<T>>,
       lindblad: &LindbladOp<T>,
       y0:       &PauliSum<T>,
       t_span:   (f64, f64),
       config:   &SolverConfig,
   ) -> f64
   ```
   If `config.h0.is_some()`, return it directly. Otherwise follow the 5-step procedure
   in `PLAN.md`. Clamp the result to `[config.hmin, config.hmax]`.

**Unit tests:**
- **`h0` specified:** `config.h0 = Some(0.1)` → returns `0.1` regardless of the system.
- **Zero RHS:** `ham = None`, empty Lindblad. With `d1 ≈ 0`, the fallback (`1e-6`) is
  returned.
- **Non-trivial:** `H = 0.5·Z`, `P = X`. Estimated `h0` is positive, finite, and
  no larger than `t_span.1 − t_span.0`.

**Review checklist:**
- [ ] `config.h0 = Some(x)` short-circuits the computation.
- [ ] All 5 steps of the Hairer procedure are implemented.
- [ ] Result is clamped to `[hmin, hmax]`.

---

## Task 9 — `solve_mut` and `solve`

**Goal:** Implement the full adaptive ODE solve loop.

**Steps:**
1. In `src/solve.rs`, implement `solve_mut`:
   - Call `estimate_h0` for the initial step.
   - Main loop advancing `t` from `t_span.0` to `t_span.1`:
     - Cap `dt` so the next `save_at` point is not overshot.
     - Call `step`; on `Reject`, retry with `h_new`.
     - On `Accept`, update state, carry `k_next` forward as `k1` (FSAL), check if a
       save point was reached and invoke `callback` if so.
   - Return `(Vec<f64>, Vec<R>)` of save times and callback results.
2. Implement `solve` as a one-liner clone + `solve_mut`.

**Unit tests:**
- **Empty `save_at`:** returns `([], [])` without panic.
- **Single save at `t_end`:**
  `H = 0.5·Z`, no Lindblad, `P(0) = X`, `t_span = (0.0, 1.0)`, `save_at = [1.0]`.
  Callback: `|_, p| p.overlap(&x_sum)` where `x_sum` is a `PauliSum` with `{X: 1.0}`.
  Exact answer: `cos(1.0) ≈ 0.5403`. Assert within `1e-4`.
- **Multiple save points:**
  Same setup, `save_at = [0.25, 0.5, 0.75, 1.0]`. Assert all four values match
  `cos(t)` within `1e-4`.
- **`solve` vs `solve_mut`:** same inputs produce identical results; `initial` is
  unchanged after `solve`.

**Review checklist:**
- [ ] Save points are hit exactly (not approximated by the nearest step).
- [ ] FSAL: `k1` of step `n+1` is the `k_next` from step `n` — not recomputed.
- [ ] `solve` does not mutate `initial`.
- [ ] Returned time vector matches `save_at` exactly.

---

## Task 10 — Integration Test: Spontaneous Emission

**Goal:** Validate the full stack against a physically meaningful problem with a known
analytic solution.

**Setup:** Single qubit, collapse operator `c = X + iY` (un-normalised lowering), rate
`γ = 1.0`, no Hamiltonian, initial observable `P(0) = Z`.

**Analytic solution:** Work out `L(Z)` for this `c` by hand (or from `PLAN.md`'s equation),
derive the ODE for `<Z>(t)`, and document the closed-form solution in the test before
asserting against it. (Do not assert against a magic number without a derivation.)

**Steps:**
1. Construct `CollapseOp`, `LindbladOp`, and initial state.
2. Solve over `t ∈ [0, 2.0]` with `save_at` at several evenly-spaced points.
3. Compare each saved value to the analytic solution.

**Verification:**
- Error relative to analytic solution is within `1e-4` at all save points.

**Review checklist:**
- [ ] Analytic solution is derived and written out in a comment before the assertion.
- [ ] Default `SolverConfig` is used (no hand-tuned tolerances).
- [ ] Test is self-contained and does not access private fields.

---

## General Review Rules

Before approving any task:

1. **Plan consistency:** re-read the relevant section of `PLAN.md` and confirm the
   implementation matches.
2. **Guidelines compliance:** no other crates modified, commit made, all tests pass,
   `cargo clippy -p ppvm-timeevolve -- -D warnings` is clean.
3. **No dead code:** no unused functions, types, or imports.
4. **Sign-off:** explicit approval required before the next task starts.
