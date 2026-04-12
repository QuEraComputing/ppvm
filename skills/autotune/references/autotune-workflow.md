# Autotune Workflow

## 1. Validate the target metric and profile

- Ensure the target metric is numeric and available from a focused benchmark.
- If no good target exists, inspect benchmarks or profile and propose one.
- Identify the exact benchmark command (e.g., `cargo bench -p crate --bench micro -- "gates/x"`).
- **Profile the benchmark to get a time breakdown.** Before the first optimization iteration, measure how time is distributed across major code sections. Use an ad-hoc timing binary with section markers (e.g., `Instant::now()` around each phase), not just the overall benchmark number. This prevents wasting iterations on non-bottlenecks. Record the breakdown in `log.md` under "Architecture Notes". Re-profile after major gains to update the breakdown — the bottleneck shifts as you optimize.

## 2. Prepare the experiment

1. Generate a memorable `<task>` name (e.g., `tableau-x-gate`).
2. Run `scripts/init_experiment.py <task>` to create `docs/autotune/<task>/metric.toml` and `log.md`.
3. Create the working branch: `git checkout -b autotune/<task>`.
4. Run the baseline benchmark, record it with `scripts/record_result.py`, and commit.

## 3. Run one iteration

### 3a. Plan the approach

1. Read `docs/autotune/<task>/metric.toml` and `log.md`.
2. Choose one hypothesis to test. Give it a short name (e.g., `simd-phase-update`).
3. Create `docs/autotune/<task>/<approach>/prompts.md` describing the approach, target metrics, and relevant prior findings.

### 3b. Implement via subagent in a worktree

Dispatch an implementation subagent using the **Agent tool with `isolation: "worktree"`**. This creates an isolated git worktree automatically. The subagent:

- Reads the prompts.md you wrote.
- Implements the approach (edits only crate `src/` directories).
- Does NOT run benchmarks — the host agent does that.
- Commits all changes before finishing.

Example Agent call:

```
Agent({
  description: "Implement <approach>",
  prompt: "<self-contained prompt with hypothesis, file scope, prior findings, and commit instruction>",
  isolation: "worktree"
})
```

The Agent tool returns the worktree path and branch name when the subagent makes changes.

### 3c. Benchmark

After the subagent returns, run the targeted microbenchmark **in the worktree directory** (use the returned path). Only run the benchmarks related to the target metric — not the full suite.

### 3d. Integrate the worktree branch

**Always merge the worktree branch to the working branch**, regardless of outcome. This preserves a record of every attempt in git history.

```bash
# Switch back to the working branch
git checkout autotune/<task>
# Merge the iteration branch
git merge <worktree-branch> --no-ff -m "merge: autotune iteration <approach>"
```

### 3e. Decide: keep or discard

- **If the metric improved (keep):** The merge stands. The code change is preserved on the working branch.
- **If the metric regressed or is unchanged (discard):** Revert the merge commit on the working branch.
  ```bash
  git revert -m 1 HEAD --no-edit
  ```
- **If the subagent crashed:** Same as discard — merge, then revert.

### 3f. Record the result

In a **separate commit** on the working branch (not mixed with code changes):

1. Append the benchmark data to `metric.toml` via `scripts/record_result.py`.
2. Append any durable finding to `log.md` via `scripts/add_log_entry.py`.
3. Commit these ledger updates.

This separation ensures that discarded experiments still have their metrics recorded even after code reverts.

### 3g. Clean up the worktree

Remove the worktree directory (git does this automatically if you used the Agent tool's `isolation: "worktree"`, but verify it's cleaned up).

## 4. Crash policy

- Fix obvious trivial failures (typo, missing import) and rerun.
- Mark fundamentally broken ideas as `crash` and continue to the next iteration.

## 5. When to escalate

When 3+ consecutive micro-optimizations show <1% improvement:

1. **Harvest gains**: Create a PR from the current working branch to lock in micro wins.
2. **Shift strategy**: Try structural optimizations — data layout changes (e.g., columnar storage), batched operations (processing multiple gates in one pass), SIMD vectorization, or algorithmic rewrites.
3. **Larger iterations are OK**: Structural changes may touch more files and take longer to implement. Use worktree isolation to keep them safe. Break large changes into compile-first-then-benchmark steps.

## 6. Role separation summary

| Responsibility | Host agent | Subagent (worktree) |
|---|---|---|
| Read ledger + choose hypothesis | Yes | No |
| Write prompts.md | Yes | No |
| Implement code changes | No | Yes |
| Run benchmarks | Yes | No |
| Merge worktree to working branch | Yes | No |
| Decide keep/discard/revert | Yes | No |
| Record results (metric.toml, log.md) | Yes | No |
| Commit ledger updates | Yes | No |
