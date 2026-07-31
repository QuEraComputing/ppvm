# `ppvm-*-2` refactor log

Status: temp / machine-maintained scratchpad.

**Ownership:** this file is written by the refactor **orchestrator** (the workflow
script), not by individual agents. Agents *read* it for context and *return*
structured findings; the orchestrator appends/updates the gap tables below. This
keeps writes race-free when review and test/bench agents run in parallel.

**Lifecycle:** resolved gaps stay (struck through) as the migration record. At each
crate/phase boundary, surviving lessons are distilled into
[`traits-2-configuration-and-hashing.md`](design/traits-2-configuration-and-hashing.md)
or the implementation plan, and this scratch is trimmed.

**Gap schema:** each gap is one row.

| field | meaning |
| --- | --- |
| `id` | `<crate>.<component>.<n>` |
| `type` | `impl-friction` \| `correctness` \| `perf-drift` \| `missing-proof` \| `missing-test` |
| `sev` | high \| medium \| low |
| `routed-to` | impl \| proof \| test \| design \| human |
| `status` | open \| closed |
| `commit` | resolving commit sha (when closed) |

Routing rule: `correctness`/`impl-friction`→impl (or design); `perf-drift`→impl or
human sign-off (hard gate + allowlist); `missing-proof`→proof; `missing-test`→test.

Convergence rule: a component is done when a full cycle yields **0 new gaps** and
`cargo build` + differential/property tests + `lake build` are all green. Max 4
iterations, then escalate to human.

---

## Phase 0 — scaffolding

Status: **complete.** Deliverables:
- `crates/ppvm-traits-2` — stub lib (trait modules land in Phase 1).
- `crates/ppvm-conformance-2` — test-only harness (`publish = false`): seeded RNG,
  `random_pauli_string`, replayable `GateOp`/`random_circuit`, `old_word_from_str`,
  `assert_close`; 4 smoke tests; `harness_smoke` Criterion bench.
- Both wired into `[workspace].members`.

Gates green: `cargo build --all-targets`, `cargo test -p ppvm-conformance-2` (4/4),
`cargo clippy -- -D warnings`, `cargo fmt --check`, `cargo machete`, and the
Criterion bench (`random_circuit/5q/1024` ≈ 9.2 µs). No gaps.

_(no gap entries)_

## Phase 1 — ppvm-traits-2

_(no entries yet)_

## Phase 2 — ppvm-pauli-word-2

_(no entries yet)_

## Phase 3 — ppvm-pauli-sum-2

_(no entries yet)_

## Phase 4 — ppvm-tableau-2

_(no entries yet)_

## Phase 5 — ppvm-sym-2

_(no entries yet)_

## Phase 6 — ColumnStore backend

_(no entries yet)_

## Phase 7 — downstream + cutover

_(no entries yet)_

---

## Perf-drift allowlist (human-approved regressions)

Regressions listed here are accepted as designed-in trade-offs and do **not** block
a crate. Anything not listed is a hard gate.

| id | component | metric | accepted ratio | justification | approved-by |
| --- | --- | --- | --- | --- | --- |
| _(none yet)_ | | | | | |
