---
title: Autotune Skill Design
date: 2026-04-10
status: approved-in-chat
---

# Autotune Skill Design

## Goal

Create a repo-local skill named `autotune` at `skills/autotune/` that lets Codex run autonomous, benchmark-driven performance tuning loops for this repository. The skill is intended for unattended runs, including overnight operation, so it should continue iterating until manually interrupted.

## Scope

The skill covers:

- validating that the requested target metric is numeric and cheaply measurable
- preparing an experiment area under `docs/autotune/<task>/`
- running an autonomous iteration loop based on worktrees and short experiments
- benchmarking only targeted microbenchmarks relevant to the current hypothesis
- recording `keep`, `discard`, and `crash` outcomes durably
- preserving experiment history even when a failed code change is reverted

The skill does not cover:

- full benchmark suite execution by default
- broad performance campaigns without a focused target metric
- unrestricted subagent edits across tests, docs, or unrelated files

## Design Summary

The skill will be markdown-first, with a concise `SKILL.md` for triggering and guardrails, a reference document for the detailed loop, and small helper scripts for repetitive filesystem and ledger updates.

Planned layout:

- `skills/autotune/SKILL.md`
- `skills/autotune/agents/openai.yaml`
- `skills/autotune/references/autotune-workflow.md`
- `skills/autotune/scripts/init_experiment.py`
- `skills/autotune/scripts/record_result.py`
- `skills/autotune/scripts/add_log_entry.py`

## Triggering And Intent

The skill should trigger when the user wants Codex to:

- improve runtime or benchmark performance in this repository
- run autonomous tuning or optimization loops
- explore multiple implementation approaches to improve a measurable metric
- operate unattended while continuing to benchmark, keep wins, and discard regressions

The frontmatter description should explicitly mention autonomous benchmark-driven tuning, worktree-based experimentation, and durable experiment logging.

## Workflow

### 1. Validate the target metric

Require the user to specify one or two numeric metrics that can be obtained from a focused benchmark or microbenchmark. Reject vague goals such as “make it faster” when no measurable target is attached.

If the requested metric is not clearly available:

- inspect the existing benchmarks
- optionally profile to identify bottlenecks
- suggest a smaller, measurable metric

### 2. Create the experiment area

Generate a memorable task name and create:

- `docs/autotune/<task>/metric.toml`
- `docs/autotune/<task>/log.md`

`metric.toml` is an append-only ledger of benchmark results. Each entry records:

- commit
- status (`keep`, `discard`, `crash`)
- description
- numeric metric keys matching benchmark output

`log.md` records only durable findings worth preserving for later study.

### 3. Run the autonomous loop

For each iteration:

1. read the existing ledger and log
2. choose one concrete hypothesis
3. create a worktree on `auto-tune-<task>-<approach>`
4. create `docs/autotune/<task>/<approach>/prompts.md`
5. spawn a subagent to implement only the approach
6. run only the relevant microbenchmarks
7. record the measured result in the main experiment ledger
8. keep the code if the target metric improves, otherwise revert the code change
9. append any durable finding to `log.md`
10. continue without asking for permission unless blocked

### 4. Preserve metrics across reverts

Discarded results must remain in history even when the code change is reverted. To achieve that:

- keep experiment records on the user branch as the canonical ledger
- separate code changes from ledger updates into different commits
- if an approach regresses, revert only the code change commit
- keep the ledger update commit so `discard` and `crash` results remain visible

This separation is mandatory for the skill design.

## Guardrails

The skill should explicitly instruct Codex to:

- avoid full benchmark suites unless the user specifically requests them
- prefer experiments that fit within roughly 10 minutes including benchmarking
- use subagents for implementation work when appropriate
- restrict implementation edits to crate source directories and the relevant `docs/autotune/` area
- avoid large documentation output beyond the experiment ledger and concise findings log
- treat crashes pragmatically: fix trivial mistakes and rerun, but log fundamentally broken ideas as `crash`

## Helper Scripts

The scripts support the workflow but do not replace agent judgment.

### `init_experiment.py`

Responsibilities:

- normalize the task name
- create `docs/autotune/<task>/`
- initialize `metric.toml`
- initialize `log.md`

### `record_result.py`

Responsibilities:

- append a metric entry to `metric.toml`
- support arbitrary numeric metric keys
- record status, commit, and description in a consistent format

### `add_log_entry.py`

Responsibilities:

- append dated findings to `log.md`
- keep formatting consistent and concise

## Agent Metadata

Generate `agents/openai.yaml` for the skill with:

- a user-facing display name
- a short description focused on autonomous performance tuning
- a default prompt that asks Codex to improve a specific benchmark metric using iterative experiments

## Validation

Validate the finished skill with Ion and the repository’s normal skill validation flow. The implementation should be considered complete only after:

- the scaffold exists under `skills/autotune/`
- the frontmatter and metadata are valid
- the helper scripts run successfully for representative inputs
- the skill structure reflects this design

## Recommended Implementation Approach

Use `ion skill new --path skills/autotune` to scaffold the skill, then replace the template with the repo-specific `autotune` instructions and helper scripts.
