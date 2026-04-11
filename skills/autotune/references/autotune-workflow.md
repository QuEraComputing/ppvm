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
8. Return to the main branch and append the measured result to the ledger in a separate commit.
9. Keep both commits if metrics improve; otherwise revert only the code commit and leave the ledger commit intact.
10. Append any durable finding to `log.md`.

## 4. Crash policy

- Fix obvious trivial failures and rerun.
- Mark fundamentally broken ideas as `crash` and continue.
