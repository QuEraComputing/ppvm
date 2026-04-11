---
title: Autotune Binary Skill Design
date: 2026-04-10
status: approved-in-chat
---

# Autotune Binary Skill Design

## Goal

Create a standalone Rust binary skill at `/Users/roger/Code/rust/autotune` that provides a deterministic, validated command surface for experiment-driven performance tuning. The agent remains responsible for the research loop and hypothesis generation, while the CLI owns the deterministic operations needed to run that loop safely and repeatably.

The binary must fit Ion's binary-skill model so agents invoke it via `ion run autotune ...` and it provides the standard `self skill`, `self info`, `self check`, and `self update` commands through `ionem`.

## Scope

The first version should:

- be generic rather than `ppvm`-specific
- read one checked-in TOML config file from the tuned repository
- manage deterministic experiment state, worktrees, benchmark execution, metric parsing, result recording, and reporting
- provide terminal-native reporting with ASCII output and JSON output for agents
- avoid requiring ad hoc Python, grep, or similar external tooling during normal agent operation

The first version should not:

- own the autonomous research loop end-to-end
- choose hypotheses or decide when to stop
- depend on `ppvm` benchmark conventions in its code
- require browser-based or image-based reporting

## High-Level Model

The system is split into two layers:

### Agent responsibilities

The agent remains the researcher. It:

- decides which metric or metrics to optimize
- proposes the next approach or hypothesis
- decides when enough iterations have run
- invokes deterministic `autotune` commands in the right order

### CLI responsibilities

The `autotune` binary is the deterministic execution engine. It:

- validates repository config and command inputs
- creates and updates experiment state
- manages worktrees and branch naming
- materializes prompt and experiment files
- runs configured benchmark commands
- parses named numeric metrics from benchmark output
- computes scores and guardrail checks
- integrates or discards code changes according to deterministic policy
- renders progress reports

This boundary is intentional. The agent keeps the flexible reasoning; the binary owns the reproducible state machine transitions.

## Packaging And Ion Integration

The binary should be a standalone Cargo crate at `/Users/roger/Code/rust/autotune`.

It should follow Ion's binary-skill conventions:

- expose a `self skill` subcommand that prints a valid binary SKILL.md
- expose `self info`, `self check`, and `self update`
- use `ionem` for the self-management plumbing
- be invokable by agents through `ion run autotune ...`

The generated SKILL.md should position the binary as a deterministic tuning orchestrator used by agents, not as an autonomous research agent itself.

## Repository Config Model

Each tuned repository provides a single checked-in TOML file, for example `.autotune.toml`.

That config file defines:

- experiment storage paths
- prompt materialization paths
- canonical branch policy
- worktree and approach naming rules
- benchmark commands that are allowed to run
- named metric extractors for benchmark output
- primary metrics
- guardrail metrics
- optional report defaults

The binary code should not encode repository-specific benchmark semantics. Repository behavior is driven by config.

## Metric Model

### Primary metrics

The config supports multiple primary metrics.

Each primary metric declares:

- metric name
- direction: `maximize` or `minimize`

For ranking, each metric is compared against the current baseline and normalized according to its direction. The total rank score is the sum of the normalized improvements across all primary metrics.

This lets the first version support multi-metric optimization without introducing arbitrary weights.

### Guardrail metrics

Guardrail metrics are separate from rank scoring.

Each guardrail metric declares:

- metric name
- direction: `maximize` or `minimize`
- regression threshold or limit

Guardrails must not regress beyond their configured tolerance. A candidate that fails any guardrail cannot be kept, regardless of its primary-metric score.

## Benchmark Interface

The binary should treat benchmark execution generically:

- run a configured command
- capture stdout and stderr
- parse named numeric metrics from structured output according to config

This avoids framework-specific code in the first version and lets `ppvm` act as the first target profile without polluting the core abstraction.

## Deterministic Command Surface

The binary should expose explicit subcommands that the agent composes into a loop.

Illustrative first-version commands:

- `autotune init`
  - validate config
  - initialize experiment state
  - record canonical branch and baseline context

- `autotune prepare-approach`
  - create deterministic approach metadata
  - create worktree and branch
  - materialize prompt or approach files

- `autotune benchmark`
  - run the configured benchmark command
  - parse named metrics
  - emit human-readable and JSON results

- `autotune record`
  - validate parsed metrics
  - compute score versus baseline
  - evaluate guardrails
  - append the iteration result to the structured ledger

- `autotune apply`
  - apply deterministic keep/discard rules
  - integrate winning code changes onto the canonical branch
  - preserve losing/crashing attempts in the ledger without integrating code

- `autotune report`
  - print terminal-native progress summaries and charts
  - optionally emit JSON

The exact names may change during implementation, but the interface should stay narrow and deterministic.

## Git And Worktree Model

The CLI should manage git and worktree operations directly so the agent does not need to orchestrate them with fragile shell glue.

### Canonical branch

The canonical branch is read from config or inferred during initialization and stored in experiment state.

### Approach branches and worktrees

For each approach, the CLI should:

- create a deterministic worktree branch name
- create a deterministic worktree path
- associate that branch and path with the approach in experiment state

### Winning iterations

For a winning attempt:

1. the implementation commit is created on the approach branch
2. the CLI integrates that code commit onto the canonical branch deterministically
3. the CLI records the winning benchmark result in a separate ledger update step

This preserves the invariant that code integration and experiment-record updates are separate actions.

### Losing or crashing iterations

For a losing or crashing attempt:

1. the implementation branch remains an isolated approach artifact
2. the code commit is not integrated onto the canonical branch
3. the measured result is still recorded on the canonical branch ledger so the attempt remains visible

This preserves experiment history without contaminating the canonical branch with losing code.

## State Model

The binary should keep structured experiment state in a repo-local experiment directory configured by the repo TOML.

State should include:

- experiment identity
- canonical branch
- iteration order
- approach metadata
- worktree and branch metadata
- benchmark command metadata
- parsed metric values
- computed score
- guardrail pass/fail state
- keep/discard/crash outcome

Commands should reconstruct their behavior from this persisted state rather than from hidden process memory.

## Reporting

`autotune report` should be terminal-native.

First-version report output should include:

- summary header with experiment name, repo, canonical branch, and baseline
- iteration table with approach, status, score, primary metrics, and guardrail result
- ASCII chart of score or a selected metric over iteration number
- `--metric <name>` to focus the chart on one metric
- `--json` for structured agent consumption

This keeps reporting self-contained and avoids browser or graphics dependencies in the MVP.

## Validation And Safety

The CLI should validate inputs before every mutating action.

Examples:

- config file schema is valid
- benchmark command is one of the configured allowed commands
- parsed metrics are finite numeric values
- primary and guardrail metric names are consistent with config
- branch and worktree state is consistent with recorded experiment state
- the requested operation is valid for the current iteration status

If a deterministic command cannot safely proceed, it should fail with a clear, typed CLI error rather than silently guessing.

## MVP Orientation

The binary is generic, but the MVP should be developed against `ppvm` as the first real consumer.

That means:

- design the core abstractions to be generic
- test early against a real `ppvm` config
- avoid embedding `ppvm`-specific definitions into the crate

`ppvm` should validate the design, not define it.

## Testing Strategy

The implementation plan should cover:

- unit tests for config parsing, metric parsing, scoring, guardrails, and report rendering
- integration tests for git/worktree flows in temp repositories
- binary-skill tests for `self skill`, `self info`, `self check`, and `self update` wiring
- fixture-driven tests for benchmark command parsing
- report tests for ASCII output and JSON output

The final binary should be trustworthy enough that agents can rely on it as their primary deterministic orchestration surface.

## Recommended Implementation Direction

Use Ion's binary-skill scaffolding for the new crate and build the CLI around a small number of deterministic state transitions. The first implementation should optimize for explicitness and auditability over automation magic.
