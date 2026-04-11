---
name: autotune
description: Run autonomous benchmark-driven performance tuning loops for this repository. Use when Codex needs to improve a specific numeric performance metric, run iterative optimization experiments, create worktree-based approaches, benchmark only targeted microbenchmarks, keep measured wins, discard regressions, and continue operating unattended while preserving experiment history in docs/autotune.
---

# Autotune

Use this skill only for a narrow, measurable performance target. If the request is vague or broad, narrow it to one or two numeric metrics before starting.

Read `references/autotune-workflow.md` for the full loop, workflow, and helper script usage.

## Core Rules

- Prefer targeted microbenchmarks and profiling over full benchmark suites.
- Keep each iteration short enough to fit roughly within 10 minutes including measurement.
- Once the loop starts, continue autonomously until interrupted unless a hard blocker requires human input.
- Restrict implementation subagents to crate source directories and the active `docs/autotune/` experiment area.
- Record every attempt with `keep`, `discard`, or `crash`.
- Keep ledger updates in a separate commit from code changes so discarded results survive code reverts.

## Experiment Setup

Use `scripts/init_experiment.py` to create `docs/autotune/<task>/metric.toml` and `docs/autotune/<task>/log.md`.

## Recording Results

Use `scripts/record_result.py` to append benchmark data to `metric.toml`.
Use `scripts/add_log_entry.py` to append durable findings to `log.md`.
