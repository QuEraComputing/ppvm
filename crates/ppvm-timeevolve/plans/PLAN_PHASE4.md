# Phase 4: Consistent Per-Stage Truncation

## Motivation

### The DOPRI5 fighting problem

In the current `dopri5.rs::step()`, truncation happens at exactly one point: after the
5th-order update `y5` is formed (line 274), and before the FSAL stage 7 evaluates
`k[6] = rhs(T(y5))`.

The Butcher tableau error estimate is:

```
e = dt * (E1·k[0] + E3·k[2] + E4·k[3] + E5·k[4] + E6·k[5] + E7·k[6])
```

with `E7 = -1/40 ≠ 0`.  The problem: `k[0]…k[5]` are derivatives of *untruncated*
stage states, while `k[6]` is the derivative of the *truncated* final state.  This
mismatch introduces a spurious "jump" into the error estimate — proportional to how much
the truncation changes the state.  DOPRI5 interprets this jump as genuine local error and
repeatedly halves the step size until `rtol` is large enough to absorb it.

The workaround used in examples (`rtol ≈ 10·min_threshold`) is a band-aid: it makes
DOPRI5 tolerant of the artificial noise, but it also prevents it from tightening below
that floor when the solution genuinely warrants it.

### The fix: per-stage truncation

Truncate `y_scratch` at every stage, not just at `y5`.  Concretely, after each
`add_scaled` accumulation that builds stage state `y_i`, call `y_scratch.truncate()`
before passing it to `rhs_into`.  This means:

- `k[i] = rhs(T(y_i))` for all i — every derivative is computed on the truncated manifold.
- The error estimate `Σ Ei·ki` is self-consistent within the truncated ODE.
- DOPRI5 sees a continuous error signal and takes large steps without fighting.
- The FSAL property is preserved: after acceptance, `k[6] = rhs(T(y5))` is exactly the
  `k[0]` the next step needs.

### Design decision: do not truncate k-buffers

The k-buffers hold derivatives (`dρ/dt`), not states.  Each `k[i]` is used only as an
additive contribution scaled by `dt * A_{s,i}` or `dt * E_i`.  After accumulation,
`y_scratch` is immediately truncated.  Truncating `k[i]` before this accumulation would
apply a further approximation — dropping terms from the *rate of change* rather than
from the *state itself* — with no clear physical justification and a different error
character.

The current design: truncate stage states; leave k-buffers un-truncated.  The strategy
embedded in `PauliSum<T>` (via `T::Strategy`) governs all state truncations uniformly.
Separate configurability for derivative truncation is not needed for Phase 4.

---

## Task 26 — Per-stage truncation in `dopri5.rs::step()`

**Goal:** Add `y_scratch.truncate()` before each `rhs_into` call in stages 2–6, making
all k_i consistent with the truncated ODE.  Remove the existing single `truncate()` at
line 274 (before stage 7) — it is subsumed by the new per-stage calls.

**Steps:**

1. In `src/dopri5.rs`, inside `step()`, after each block that accumulates `y_scratch`
   for stages 2–6 and before the `rhs_into` call, insert `y_scratch.truncate();`:

   ```rust
   // Stage 2
   {
       let (lo, hi) = k.split_at_mut(1);
       y_scratch.data_mut().clone_from(y.data());
       add_scaled(y_scratch, &lo[0], dt * A21);
       y_scratch.truncate();                      // ← new
       rhs_into(ham, lindblad, y_scratch, &mut hi[0]);
   }
   // (same pattern for stages 3, 4, 5, 6)
   ```

2. Retain the existing `y_scratch.truncate()` before stage 7 — it is still correct and
   semantically clearer to keep it explicit before the FSAL evaluation.

   After the change, `truncate()` is called 6 times per step: once at the end of each
   of stages 2–6 (before `rhs_into`) and once before stage 7.  The existing call before
   stage 7 is now redundant with the per-stage pattern, but keeping it makes the intent
   explicit: "this is the accepted state; truncate it".  If profiling shows it is a
   measurable cost, it can be removed in a follow-up.

3. Update the doc-comment on `step()` to describe the truncation behaviour.

4. Remove the `rtol` matching heuristic from the `SolverConfig` doc-comment (the
   `rtol ≈ 10·min_threshold` guidance is no longer needed).

**Cost note (known and accepted):** Truncation is now O(|P|·log|P|) per stage for
`Budget` (binary-search pass) vs O(|P|) for `CoefficientThreshold` and `MaxPauliWeight`.
For strategies with cheap `truncate`, the added 5 calls are negligible.  For `Budget`,
there is a measurable per-step overhead.  This is treated as a known cost; Task 28
quantifies it.

**Unit tests:**

- `per_stage_truncation_error_signal_is_smooth`: drive `step()` with a `Budget`
  strategy (target=5) on a system where the natural state has ~10 terms.  Without
  per-stage truncation, the step would be rejected; with it, the step should be accepted
  at the default `rtol=1e-6`.  Assert `StepResult::Accept`.
- `per_stage_truncation_fsal_consistent`: after a successful step, verify that
  `cache.k[6]` equals `rhs(y_scratch)` computed independently (existing
  `step_fsal_k_next_equals_rhs_y_new` test pattern, adapted for Budget).
- Regression guard: all existing `step_*` and `truncate_*` tests must still pass.

**Review checklist:**
- [ ] Six `truncate()` calls total: one per stage 2–6, one before stage 7.
- [ ] k-buffers are not truncated.
- [ ] `step_fsal_k_next_equals_rhs_y_new` still passes (FSAL consistency).
- [ ] `truncate_state_does_not_accumulate` still passes.
- [ ] No rtol heuristic guidance remains in doc-comments.
- [ ] `cargo test -p ppvm-timeevolve` and `cargo clippy -p ppvm-timeevolve -- -D warnings` clean.

---

## Task 27 — Remove `Budget::min_threshold`; use `CombinedStrategy`

**Goal:** Simplify `Budget` to a pure count-cap strategy.  The `min_threshold` field is
redundant with `CoefficientThreshold`; users who want both behaviours can compose them
via `CombinedStrategy<Budget, CoefficientThreshold>` from `ppvm_runtime`.

**Steps:**

1. In `src/strategy.rs`, change `Budget` to:
   ```rust
   pub struct Budget {
       pub target: usize,
   }
   ```
   Remove `min_threshold` from `Default`, `truncate`, the doc-comment, and everywhere
   else in the file.

2. In `Budget::truncate`, remove Step 1 (the `map.retain(|_, v| !v.cutoff(…))` call).
   Keep Step 2 onward (collect→binary-search→retain).  Update the exponential search
   starting point: start from `f64::MIN_POSITIVE` (there is no `min_threshold` anchor).

3. Update `Budget::capacity`: remove the `usize::MAX / 2` check; always return
   `self.target`.

4. Update all call sites inside `ppvm-timeevolve` (tests, examples, benches) to use
   `Budget { target: N }` instead of `Budget { target: N, min_threshold: T }`.

5. Update the doc-comment on `Budget` to describe it as a *pure count cap* and explain
   that combining with coefficient pruning requires
   `ByteF64<N, CombinedStrategy<Budget, CoefficientThreshold>>`.

6. Update `examples/superradiance.rs` and `examples/scaling.rs` to use the new
   `Budget { target }` form.  (The `min_threshold`-derived `rtol` heuristic is also
   removed; per-stage truncation makes it unnecessary.)

**Unit tests:**

- Update existing Budget tests to use `Budget { target: N }`.
- `budget_no_threshold_keeps_largest`: build a `PauliSum` with 10 terms of known
  magnitudes; apply `Budget { target: 5 }`; assert the 5 largest survive.
- `combined_strategy_matches_separate`: verify that
  `CombinedStrategy<Budget, CoefficientThreshold>` applied to the same map yields the
  same result as applying `CoefficientThreshold` first and then `Budget`.

**Review checklist:**
- [ ] `Budget` has only `target: usize`.
- [ ] `Budget::truncate` no longer calls `cutoff(min_threshold)`.
- [ ] All `Budget { target, min_threshold }` literals are gone from the crate.
- [ ] Doc-comment explains composition with `CombinedStrategy`.
- [ ] All existing tests pass with updated construction syntax.

---

## Task 28 — Benchmark truncation overhead and document the decision on `Budget`

**Goal:** Measure the wall-time cost of per-stage truncation for each strategy and record
the findings.  Based on the results, document whether `Budget` is recommended or not.

**What to measure:**

Run `cargo bench` before and after Tasks 26 + 27 for:

1. `bench_rhs` (existing fixture, no truncation in hot path) — should be unchanged.
2. New bench `bench_step_ct` — one full `step()` call with `CoefficientThreshold(1e-6)`,
   n=5, ~25 live terms.
3. New bench `bench_step_mpw` — one full `step()` call with `MaxPauliWeight(2)`, n=5.
4. New bench `bench_step_budget` — one full `step()` call with `Budget { target: 300 }`,
   n=5.

**Steps:**

1. Add the three new benches in `benches/` (extend `benches/solve.rs` or add a new file).
2. Run all four benches before (at the Task 26 commit) and after (at the Task 28 commit);
   record the numbers.
3. In a comment in `src/strategy.rs` above the `Budget` struct, write the decision:
   if `bench_step_budget` is ≥ 2× slower than `bench_step_ct`, state that `Budget` is
   *not recommended as a default strategy* due to its O(n log n) truncation cost per
   stage, and that `CoefficientThreshold` with an appropriate threshold is preferred.
   Do not deprecate `Budget`; leave it available for users who genuinely need a hard cap.

**Review checklist:**
- [ ] Three new benches present and runnable.
- [ ] Before/after numbers recorded in the commit message.
- [ ] Budget recommendation comment added to `src/strategy.rs`.
- [ ] No regressions in `bench_rhs`.

---

## Task 29 — Update examples and docs; restore core message

**Goal:** Update the two examples and all doc-comments to reflect the new behaviour.  The
superradiance example should revert to its original core message: *you gain performance
with acceptable accuracy by choosing any truncation strategy appropriately*.  The
`rtol`-matching heuristic is no longer needed and should be removed.

**Steps:**

1. **`examples/superradiance.rs`**: Remove `BUDGET_RTOL`, `BUD_MIN_THRESH` and the
   `SolverConfig { rtol: BUDGET_RTOL, … }` overrides.  All three variants use
   `SolverConfig::default()`.  Rewrite the doc-comment at the top to say: comparing
   Baseline (generous budget, acts like no truncation), Budget (tight cap, small
   accuracy loss), showing that truncation can accelerate integration.  The example
   should demonstrate that per-stage truncation now eliminates step rejections without
   any manual rtol tuning.

2. **`examples/scaling.rs`**: Remove `CT_RTOL`, `BUD_RTOL`, `BUD_MIN_THRESH`.  All
   `run_ct` calls use `SolverConfig::default()`.  Update the constants section and the
   trailing `Notes:` block accordingly.

3. **`src/solve.rs` (or wherever `SolverConfig` is documented)**: Remove the
   `rtol ≈ 10·min_threshold` guidance from the doc-comment.

4. **`src/strategy.rs`**: Update the `Budget` doc-comment (already touched in Task 27)
   to remove any reference to rtol coupling.

5. Confirm that both examples run and produce sensible output:
   - `cargo run --example superradiance --release`
   - `cargo run --example scaling --release`

**Review checklist:**
- [ ] No `CT_RTOL`, `BUD_RTOL`, `BUDGET_RTOL`, `BUD_MIN_THRESH` in examples.
- [ ] `SolverConfig::default()` used throughout both examples.
- [ ] Superradiance doc-comment no longer mentions "rtol coupling".
- [ ] Scaling doc-comment updated to match new strategy list.
- [ ] Both examples compile and produce output without step-rejection warnings.
- [ ] `cargo test -p ppvm-timeevolve` and `cargo clippy -p ppvm-timeevolve -- -D warnings` clean.

---

## Summary of changes across Phase 4

| Change | Task |
|--------|------|
| Per-stage `truncate()` in `step()` (stages 2–6 + y5) | 26 |
| Remove `Budget::min_threshold` | 27 |
| Benchmark and document Budget overhead | 28 |
| Remove rtol heuristics from examples and docs | 29 |

After Phase 4, `CoefficientThreshold` and `MaxPauliWeight` are the recommended strategies
for typical use.  `Budget` remains available for applications that need a strict memory
cap.  The `rtol` heuristic is eliminated: users set `rtol` based on physical accuracy
requirements, not truncation noise.
