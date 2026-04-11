# Autotune Workflow

## 1. Validate the target metric

- Ensure the target metric is numeric and available from a focused benchmark.
- If no good target exists, inspect benchmarks or profile and propose one.

## 2. Prepare the experiment

- Generate a memorable `<task>` name.
- Create `docs/autotune/<task>/`.
- Initialize `metric.toml` and `log.md`.

## 3. Run one iteration

1. Read the existing ledger and findings log.
2. Choose one hypothesis.
3. Create a worktree branch named `auto-tune-<task>-<approach>`.
4. Create `docs/autotune/<task>/<approach>/prompts.md`.
5. Dispatch a subagent to implement only the approach.
6. Run the relevant microbenchmarks.
7. Commit the code changes on the worktree branch first.
8. If the code change wins, switch back to the canonical branch, integrate the code commit from the worktree branch there, and keep that code commit separate from the later ledger update.
9. Append the measured winning result to the ledger on the canonical branch in a separate commit.
10. If the code change loses or crashes, do not integrate the code commit into the canonical branch; still append the measured result to the ledger on the canonical branch in a separate commit so the attempt is preserved.
11. Append any durable finding to `log.md`.

## 4. Crash policy

- Fix obvious trivial failures and rerun.
- Mark fundamentally broken ideas as `crash` and continue.
