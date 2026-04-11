---
name: autotune
description: Run autonomous benchmark-driven performance tuning loops for this repository. Use when Codex needs to improve a specific numeric performance metric, run iterative optimization experiments, create worktree-based approaches, benchmark only targeted microbenchmarks, keep measured wins, discard regressions, and continue operating unattended while preserving experiment history in docs/autotune.
---

# Autotune

Use this skill only for a narrow, measurable performance target. If the request is vague or broad, narrow it to one or two numeric metrics before starting.

## Core Rules

- Prefer targeted microbenchmarks and profiling over full benchmark suites.
- Once the loop starts, continue autonomously until interrupted unless a hard blocker requires human input.
- Restrict implementation subagents to crate source directories and the active `docs/autotune/` experiment area.
- Record every attempt with `keep`, `discard`, or `crash`.
- Keep ledger updates in a separate commit from code changes so discarded results survive code reverts.

## Escalation Strategy

The autotune loop has three phases. Micro-optimizations are where you start, not where you stop.

### Phase 1: Micro-optimizations
Quick, surgical changes — precomputing values, replacing algorithms for hot-path operations, avoiding allocations, changing data structures (e.g. HashMap hasher). Each iteration is small and fast to benchmark. This is the default starting phase.

### Phase 2: Harvest and escalate
When 3+ consecutive micro-optimizations show <1% improvement, the micro well is dry. Do NOT conclude "diminishing returns" and stop. Instead:
1. **Harvest**: Create a PR on a new branch from the current gains so the user can review and merge the micro wins independently.
2. **Escalate**: Return to the working branch and shift to structural optimizations — data layout changes, batched operations, algorithmic rewrites, SIMD. These changes are larger, may touch more files, and take longer per iteration. That is expected and acceptable.

### Phase 3: Architectural changes
Major refactors — columnar data layouts, new data structures, parallelism. These may require multiple subagent iterations to get right. Use worktree isolation aggressively. If a change is large, break it into sub-steps: first make it compile and pass tests, then benchmark. A failed architectural attempt that compiles and tests correctly is valuable — it narrows the search space.

The key principle: **consecutive failures at one level of abstraction are a signal to move up, not to stop.** "I ran out of micro-optimizations" means "time for structural changes", not "time to give up."

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

Once the experiment loop has begun, do NOT pause to ask the human if you should continue. The human may be away and expects you to work autonomously until manually interrupted. If you run out of ideas, think harder — re-read the code for new angles, try combining previous near-misses, try more radical approaches.

If micro-optimizations plateau (3+ consecutive <1% results):
1. Create a PR to harvest current gains (so they're not lost if structural changes break things).
2. Escalate to Phase 2/3 — try structural changes like data layout rewrites, batched operations, SIMD, or algorithmic improvements.
3. If a structural change is risky, the worktree isolation protects the working branch. Failed experiments get reverted as usual.

Never conclude "diminishing returns" as a reason to stop. Instead, escalate the approach. The loop runs until the human interrupts.
