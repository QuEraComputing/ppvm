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
7. Append the measured result on the main branch ledger.
8. Keep the code change if metrics improve; otherwise revert only the code commit.
9. Append any durable finding to `log.md`.

## 4. Crash policy

- Fix obvious trivial failures and rerun.
- Mark fundamentally broken ideas as `crash` and continue.
