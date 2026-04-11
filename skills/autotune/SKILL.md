---
name: autotune
description: Run autonomous benchmark-driven performance tuning loops for this repository. Use when Codex needs to improve a specific numeric performance metric, run iterative optimization experiments, create worktree-based approaches, benchmark only targeted microbenchmarks, keep measured wins, discard regressions, and continue operating unattended while preserving experiment history in docs/autotune.
---

# Autotune

Use this skill only for a narrow, measurable performance target. If the request is vague or broad, narrow it to one or two numeric metrics before starting.

## Core Rules

- Prefer targeted microbenchmarks and profiling over full benchmark suites.
- Keep each iteration short enough to fit roughly within 10 minutes including measurement.
- Once the loop starts, continue autonomously until interrupted unless a hard blocker requires human input.
- Restrict implementation subagents to crate source directories and the active `docs/autotune/` experiment area.
- Record every attempt with `keep`, `discard`, or `crash`.
- Keep ledger updates in a separate commit from code changes so discarded results survive code reverts.

## Branching Model

The autotune loop uses a two-level branch structure to keep the user's starting branch clean:

1. **Working branch** (`autotune/<task>`): Created at the start of the experiment from the user's current HEAD. All experiment commits (code changes, metric records, log entries, reverts) happen here. The user's starting branch is never touched.
2. **Iteration branches** (created automatically by the Agent tool's `isolation: "worktree"` parameter): Each iteration's implementation work happens in an isolated worktree on its own branch. After benchmarking, the iteration branch is always merged back to the working branch — even if the result is discarded — so there is a complete record of every attempt.

### Setup sequence

```
git checkout -b autotune/<task>   # create working branch from current HEAD
```

### After each iteration

- **Always** merge the iteration's worktree branch back to the working branch. This preserves the implementation history regardless of outcome.
- **If keep**: the merge stands as-is. Record result in a separate commit.
- **If discard**: after merging, revert the merge commit on the working branch, then record the result. This way the attempt is visible in history but the code is undone.
- **If crash**: same as discard — merge, revert, record.

## Experiment Setup

Use `scripts/init_experiment.py` to create `docs/autotune/<task>/metric.toml` and `docs/autotune/<task>/log.md`.

## The Iteration Loop

See `references/autotune-workflow.md` for the detailed step-by-step. The high-level flow for each iteration:

1. **Plan**: Read the ledger and log, choose one hypothesis.
2. **Implement**: Dispatch a subagent with `isolation: "worktree"` to implement the approach. The subagent only writes code — it does not run benchmarks.
3. **Benchmark**: The host agent runs the relevant microbenchmarks on the worktree code.
4. **Integrate**: Merge the worktree branch to the working branch (always).
5. **Decide**: If the metric improved, keep. Otherwise, revert the merge on the working branch.
6. **Record**: Append the result to `metric.toml` and any findings to `log.md` in a separate commit.

### Subagent Prompts

When dispatching the implementation subagent, the prompt must be self-contained because the subagent has no memory of this conversation. Include:

- What to optimize and why (the hypothesis).
- Which files are in scope (crate `src/` directories only).
- What NOT to touch (benchmarks, tests, docs beyond `docs/autotune/`).
- Relevant findings from previous iterations (copy from `log.md`).
- The specific benchmark command that will be used to measure (so the subagent understands the target).
- Instruction to commit all changes before finishing.

### Crash Policy

- Fix obvious trivial failures (typo, missing import) and rerun.
- Mark fundamentally broken ideas as `crash` and continue.

## Recording Results

Use `scripts/record_result.py` to append benchmark data to `metric.toml`.
Use `scripts/add_log_entry.py` to append durable findings to `log.md`.

## Autonomy

Once the experiment loop has begun, do NOT pause to ask the human if you should continue. The human may be away and expects you to work autonomously until manually interrupted. If you run out of ideas, think harder — re-read the code for new angles, try combining previous near-misses, try more radical approaches. The loop runs until the human interrupts.
