---
title: ppvm-stim Test Corpus Expansion
date: 2026-04-28
status: approved-in-chat
---

# ppvm-stim Test Corpus Expansion

## Goal

Expand the `ppvm-stim` test suite from its current shape (8 hand-written fixtures, ~120 unit tests) into a substantive corpus covering the full set of phase-1-supported instructions, edge cases, and noise channels — plus a structurally-identical "phase-1-unsupported" corpus that flips to positive correctness tests as phase-2 lifts restrictions.

The expansion uses Stim's own sampler as a cross-check oracle for every fixture that exercises only Stim-supported features. Cross-checking happens at dev-time (regenerator script run by humans with `stim` installed); test-time runs are fully self-contained, deterministic, and flake-free.

## Non-Goals

- Reaching 100% line coverage. Coverage emerges as a side effect, not a target.
- Cross-checking ppvm-specific dialect instructions (`I[R_X(theta=…)]`, `S[T]`) against Stim — Stim doesn't simulate them. Those use deterministic-mode reference outputs derived by hand or by running ppvm itself.
- **Loss in the corpus.** Stim does not simulate loss (`I_ERROR[loss]` / `I_ERROR[correlated_loss]`), so loss circuits have no Stim oracle and don't fit the corpus's cross-check architecture. Loss is exercised by standalone Rust tests in `crates/ppvm-stim/tests/executor.rs` (already covering the main paths) and by Python tests in `ppvm-python/test/`. The corpus's measurement assertions assume every measurement returns `Some(bool)`.
- Property-based testing with `proptest` / `quickcheck`. Considered and deferred — useful but a sizable infra add and orthogonal to corpus expansion. Track separately.
- CI tiering. We aren't splitting fast vs slow tests yet. If runtime becomes a problem we can revisit.
- Performance regression benchmarks. Already partially covered by `tableau-msd-stim`; threshold work is bench-suite scope.

## Approach in One Picture

```
                 ┌─────────────────────┐
                 │ regen-stim CLI      │  Run by humans (uv-managed venv, `stim` dep).
                 │ (dev-time only)     │  Not invoked in CI.
                 └────────┬────────────┘
                          │ produces:
                          ▼
       ┌──────────────────────────────────────────────────┐
       │  crates/ppvm-stim/tests/data/                    │
       │    <category>/                                   │
       │      <name>.stim          (committed)            │
       │      <name>.expected.json (committed)            │
       └────────────┬─────────────────────────────────────┘
                    │ consumed every cargo test by:
                    ▼
            ┌───────────────────┐
            │ stim_corpus.rs    │  No external deps, fully deterministic.
            │ (test harness)    │  Bit-exact assertions.
            └───────────────────┘
```

## Two-Tier Verification

A central design decision: **the cross-check against Stim happens at regen time, not at test time.** Tests assert ppvm matches its own committed reference output bit-for-bit; the "do we agree with Stim?" signal lives in the regen script and is enforced before a fixture is committed.

Why: ppvm and Stim use different RNGs internally, so even with identical numerical seeds they produce different bit streams. Bit-exact comparison between the two simulators is impossible. Statistical comparison is, but introduces flake risk if done at test time.

The two-tier scheme keeps both signals:

- **Regen time** — script runs Stim with `stim_num_shots` to compute a high-confidence reference distribution, then tries ppvm seeds until it finds one whose empirical distribution at `num_shots` is within tolerance of Stim's reference. If no seed in `[0..32]` passes, regen errors and refuses to commit — that's a real correctness divergence worth investigating.
- **Test time** — harness runs ppvm at the committed seed for `num_shots` shots, computes per-bit means, asserts bit-exact f64 equality with the committed `ppvm_bit_means`. Same code + same seed + same shot count → identical f64 means every run, by construction.

Failure mode at test time can only happen if ppvm's behavior changes — which is exactly what we want to detect. The Stim cross-check is preserved as the "did our committed reference data drift from Stim?" signal at regen time.

## Fixture Format

Every fixture is two committed files: `<name>.stim` (the source circuit) and `<name>.expected.json` (the test's expectation). The JSON uses one of three modes:

### Mode 1: deterministic

For circuits whose measurement outcomes are uniquely determined regardless of RNG seed — i.e. no noise instructions **and** no quantum measurement randomness. Examples: `X 0; M 0` (always 1), `H 0; H 0; M 0` (always 0), `CX 0 1; M 0 1` after preparing both qubits. One shot suffices.

This mode is narrower than "no noise" — `H 0; M 0` produces a 50/50 random outcome despite having no noise, and goes in distribution mode, not here.

```jsonc
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true, false, false, true, true, false]
  // length matches the circuit's measurement count
  // every entry is a concrete bool — see Non-Goals re: loss
}
```

Test action: parse + normalize + execute once with `ppvm_seed`; for each measurement result, expect `Some(bool)` matching the corresponding `bitstring` entry. Encountering `None` (loss) is a hard test failure — the corpus excludes loss.

### Mode 2: distribution

For circuits with any source of randomness — either noise instructions or quantum measurement randomness from non-trivial superpositions. This is the **default mode** for most fixtures; deterministic is the narrow special case.

```jsonc
{
  "mode": "distribution",
  "num_shots": 256,
  "ppvm_seed": 7,
  "ppvm_bit_means": [0.4961, 0.5078, 0.0, 0.0078, ...],
  // Documentation + regen reproducibility — not consumed by tests:
  "stim_seed": 0,
  "stim_num_shots": 10000,
  "stim_bit_means": [0.5, 0.5, 0.0, 0.01, ...],
  "tolerance_sigma_at_regen": 5.0,
  "stim_version": "1.15.0"
}
```

Test action: parse + normalize; call `sample(prog, num_shots, || GeneralizedTableau::new_with_seed(n, 1e-10, ppvm_seed))`; compute per-measurement empirical means; assert exact f64 equality with `ppvm_bit_means`.

### Mode 3: unsupported

For fixtures using a phase-1-unsupported instruction. The Stim reference data is pre-recorded (Stim *can* simulate these instructions today), so flipping a fixture to "supported" in phase-2 requires only adding a `ppvm_seed` + `ppvm_bit_means` to the JSON — no `.stim` change.

```jsonc
{
  "mode": "unsupported",
  "awaiting_phase2_instruction": "SWAP",
  // Pre-recorded for the day phase-2 enables this instruction:
  "stim_seed": 0,
  "stim_num_shots": 10000,
  "stim_bit_means": [0.5, 0.5, 0.0, ...],
  "tolerance_sigma_at_regen": 5.0,
  "stim_version": "1.15.0"
}
```

Test action: parse + normalize; assert `Err(NormalizeError::Unsupported { name: awaiting_phase2_instruction, .. })`.

When phase-2 lands, regen on this fixture:
1. Picks a `ppvm_seed`, runs ppvm at `num_shots`, verifies vs the pre-recorded `stim_bit_means` within tolerance.
2. Adds `num_shots`, `ppvm_seed`, `ppvm_bit_means`.
3. Changes `mode` from `"unsupported"` to `"distribution"` (or `"deterministic"` if the fixture has no noise), removes `awaiting_phase2_instruction`.

Test starts checking real correctness automatically; the `.stim` file never changes.

## Categories

Fixtures live in subdirectories under `crates/ppvm-stim/tests/data/`. Each subdir has its own README documenting provenance and regen invocation.

| Subdir | Source | Approx count | Purpose |
|---|---|---|---|
| `generated/codes/` | `stim gen` with sweeps over `(code, distance, rounds, basis, noise)` | 80–100 | Real-world fault-tolerance circuits; noisy variants exercise hundreds of `DEPOLARIZE1/2`, `X_ERROR`, reset-flip noise instructions per circuit |
| `generated/noise_sweeps/` | Programmatic, per-channel sweeps | 30–50 | Per-`NoiseKind` parameter sweeps producing dozens of generated fixtures across probability values — covers `PauliChannel1/2` and skewed parameter spaces that `stim gen` doesn't naturally emit |
| `generated/random/` | Our own random-walk Python generator | 30–40 | Random sequences of supported instructions; parameterized to vary noise density (low/medium/high) so a meaningful fraction is noise-heavy |
| `generated/dialect/` | Programmatic via Python + Stim's API | 10–15 | ppvm-specific dialect (`I[R_X/Y/Z(...)]`, `I[U3(...)]`, `S[T]`, `S_DAG[T]`); mostly deterministic mode (no Stim oracle for these) |
| `unsupported/` | `stim gen` natural emission + programmatic templates | 20–25 | One fixture per phase-1-unsupported instruction; flips in phase-2 |
| `edge_cases/` | Hand-written | 20–25 | Empty programs, deeply nested REPEAT, every Pi-expression form, every tag form, dense/sparse measurement patterns, comment/whitespace stress, boundary noise probabilities (p=0.0, p=1.0) |
| `noise_channels/` | Hand-written | 8–10 | Hand-written corner cases the sweeps don't naturally cover (e.g. `M(0.0)` ≡ `M`; `MR(1.0)` always-flips-recorded-bit; readout-noise-on-already-flipped bits) |

**Total: ~200–250 committed fixtures** (`.stim` + `.expected.json` pairs).

Noise gets multi-source coverage: every realistic-noise surface code in `generated/codes/` (~30–50 of them) exercises DEPOLARIZE1/2 and X_ERROR heavily; `generated/noise_sweeps/` adds explicit per-channel breadth; `generated/random/` adds noise-heavy random programs; and `noise_channels/` + `edge_cases/` cover boundary values. Cumulatively this is the bulk of the fixture corpus.

### `generated/codes/` — `stim gen` sweeps

The regen script invokes `stim gen` for combinations of:

- `--code` ∈ {`surface_code`, `repetition_code`, `color_code`}
- `--task` ∈ all task variants the code supports (e.g. `surface_code` → `unrotated_memory_x`, `unrotated_memory_z`, `rotated_memory_x`, `rotated_memory_z`)
- `--distance` ∈ {3, 5, 7}
- `--rounds` ∈ {1, 3, 5}
- noise flags: either none, or `--after_clifford_depolarization=p --before_round_data_depolarization=p --before_measure_flip_probability=p --after_reset_flip_probability=p` for `p` ∈ {0.001, 0.01}

Some `(code, task)` combinations naturally emit phase-1-unsupported instructions (e.g. `unrotated_memory_z` uses `MR`; certain transversal tasks use `SWAP`; stabilizer-measurement tasks use `MPP`). Those land in `unsupported/` instead of `generated/codes/` — the regen script auto-routes by inspecting the emitted source.

### `generated/dialect/` — ppvm-specific instructions

Stim itself cannot simulate `I[R_X(theta=…)]`, `I[R_Y(theta=…)]`, `I[R_Z(theta=…)]`, `I[U3(theta=…, phi=…, lambda=…)]`, `S[T]`, or `S_DAG[T]`, so Stim cross-check is impossible. Each fixture combines one or more dialect instructions in a small circuit whose expected outcome is hand-derivable (e.g. `I[R_X(theta=1.0*pi)] 0; M 0` → always 1; `S[T] 0; S_DAG[T] 0; M 0` after `X 0` → always 1).

Where the resulting circuit is uniquely determined (no superposition at measurement time), the fixture is **deterministic mode**. Where it has measurement randomness (e.g. `H 0; I[R_X(theta=0.5*pi)] 0; M 0`), it's **distribution mode** with `ppvm_bit_means` recorded by running ppvm itself — no oracle, but the test still locks down the output against future regression. The dialect-only restriction means we can't catch ppvm bugs by comparing against another simulator, only by detecting drift from a previously-recorded ppvm output.

~10–15 fixtures total, programmatically generated by a per-instruction template that emits a known-outcome circuit.

### `generated/noise_sweeps/` — per-channel parameter sweeps

A Python script emits one `.stim` per `(NoiseKind, qubit_count, probability)` triple. The circuit shape per fixture is small and uniform: prepare a known state → apply the noise channel under test N times → measure all qubits. Sweeping the probability axis gives the regen-time Stim cross-check enough breadth to catch slope/sign bugs on each channel.

Sweep axes:

| Channel | Probabilities | Qubit counts | Approx fixtures |
|---|---|---|---|
| `DEPOLARIZE1` | 0.001, 0.01, 0.1, 0.5 | 1, 4 | 8 |
| `DEPOLARIZE2` | 0.001, 0.01, 0.1, 0.5 | 2, 4 | 8 |
| `PAULI_CHANNEL_1` | 3 hand-picked skewed parameter sets | 1, 4 | 6 |
| `PAULI_CHANNEL_2` | 3 hand-picked skewed parameter sets (15-arg) | 2, 4 | 6 |
| `X_ERROR` / `Y_ERROR` / `Z_ERROR` | 0.001, 0.01, 0.5 | 1 | 9 |
| readout-noise: `M(p)` / `MR(p)` | 0.001, 0.01, 0.5 | 1 | 6 |

~30–50 fixtures total (the upper bound includes optional sweeps over different state-prep contexts: `M(p)` after `X 0` vs after `H 0` etc., to validate Stim agreement on both rare-1 and 50/50 outcomes).

Stim cross-check works on every entry except `Loss`/`CorrelatedLoss` (which Stim doesn't support — those are out of corpus scope per Non-Goals). Distribution mode for everything except possibly the `M(0.0)`/`MR(0.0)` cases which collapse to deterministic.

### `generated/random/` — random-walk programs

A Python script emits random sequences of supported instructions. Parameterized by:
- Number of qubits ∈ {2, 4, 8, 16}
- Number of instructions ∈ {10, 50, 200}
- Mix ratios — three regimes: `clifford-only`, `clifford+noise (high noise density, ~30% of instructions)`, `clifford+noise+measurement-readout`
- RNG seed (8 distinct seeds per combo for diversity)

The high-noise-density regime is intentional: it produces fixtures with many noise instructions per circuit, complementing `noise_sweeps/`'s "single channel per circuit" coverage with "many mixed channels per circuit".

Cross-check via Stim works for any sequence using only Stim-supported instructions; falls back to deterministic mode for sequences mixing in dialect-specific ones (though the random generator can be configured to omit dialect instructions for full Stim coverage).

### `unsupported/` — flips in phase-2

One fixture per phase-1-unsupported instruction (~20 instructions: `SWAP`, `ISWAP`, `ISWAP_DAG`, `SQRT_XX`, `SQRT_YY`, `SQRT_ZZ`, `CXSWAP`, `SWAPCX`, `XCX`, `XCY`, `XCZ`, `YCX`, `YCY`, `YCZ`, `C_XYZ`, `C_ZYX`, `H_XY`, `H_YZ`, `MX`, `MY`, `MRX`, `MRY`, `MXX`, `MYY`, `MZZ`, `MPP`, `HERALDED_ERASE`, `HERALDED_PAULI_CHANNEL_1`, `CORRELATED_ERROR`, `ELSE_CORRELATED_ERROR`).

Each fixture is a small circuit that exercises the instruction in context: prepare a state → apply the unsupported gate → measure. Sourced either from `stim gen`'s natural emission (when applicable) or from a programmatic per-instruction template. Pre-recorded Stim reference is committed in the JSON so phase-2 flipping doesn't require Stim being reinstalled.

### `edge_cases/` — hand-written corner cases

Specific cases covered:
- Empty program (`""`).
- Whitespace-only program.
- Comments and blank lines interspersed; line-number assertions.
- Every `pi_expr` form: bare `pi`, `<coeff>*pi`, plain `f64`, very small / very large coefficients.
- Every tag shape: bare ident, single positional, single named, multiple named (with order preservation), multiple positional, mixed positional + named.
- Single-line REPEAT, multi-line REPEAT, nested REPEAT (depth 3), REPEAT with single-instruction body.
- Annotation `rec[-k]` targets (which the parser tolerates and discards for annotations only).
- Dense measurement: `M 0 1 2 3 ... 63` on a 64-qubit tableau.
- Sparse measurement: many `M`s with one target each interleaved with gates.
- Maximum-target-count circuits (limited by tableau size).
- Programs that mix every supported gate family in a single circuit.

### `noise_channels/` — hand-written corner cases

This category is the smallest of the noise-exercising subdirs — the per-channel breadth lives in `generated/noise_sweeps/`. What stays here are small hand-written corner cases that the parameter sweeps don't naturally produce:
- Boundary probabilities: `M(0.0) 0` (must equal noiseless `M`), `MR(1.0) 0` (always-flips-recorded-bit but resets correctly), `X_ERROR(1.0) 0` followed by readout (recorded bit is forced).
- Ordering invariants: noise applied between two measurements of the same qubit; verifies the recorded bit reflects state at that moment.
- Composition: `DEPOLARIZE1` followed by `M(p)` on the same qubit; verifies the two noise sources compose multiplicatively as expected.

Higher `num_shots` (4096) than other categories because the circuits are tiny (single shot is fast) and statistical drift in noise is the exact failure mode these locked-down corner cases exist to detect. ~8–10 fixtures.

**Loss tests are deliberately absent.** Per Non-Goals, loss has no Stim oracle and is exercised by Rust unit tests in `crates/ppvm-stim/tests/executor.rs` (which already cover loss-marking, correlated-loss, and the `Option::None` measurement path).

## Per-Category Shot Defaults

Tests bit-exact-compare ppvm's empirical means against committed `ppvm_bit_means`, so `num_shots` only affects test runtime — not test precision. The defaults below balance regen-time cross-check tightness against test-time speed:

Per-fixture, the JSON's `num_shots` field is the source of truth — the regen script picks defaults at fixture-creation time per the table below. The category-level defaults are:

| Category | num_shots | per-shot cost (release) | per-fixture wall time |
|---|---|---|---|
| `edge_cases/` deterministic | 1 | <0.1 ms | <1 ms |
| `edge_cases/` distribution | 256 | ~0.1–0.5 ms | ~50 ms |
| `generated/codes/` distance ≤ 5 | 256 | ~0.5 ms | ~125 ms |
| `generated/codes/` distance = 7 | 64 | ~5 ms | ~325 ms |
| `generated/noise_sweeps/` | 4096 | ~0.05 ms | ~200 ms |
| `generated/dialect/` deterministic | 1 | <0.1 ms | <1 ms |
| `generated/dialect/` distribution | 256 | ~0.5 ms | ~125 ms |
| `generated/random/` | 128 | ~1 ms | ~130 ms |
| `noise_channels/` | 4096 | ~0.05 ms | ~200 ms |
| `unsupported/` | n/a | n/a | <1 ms (parse + normalize only) |

Note that "deterministic" is selected based on whether the circuit's outcome is uniquely determined by the program (no measurement randomness, no noise) — not just on the absence of noise instructions. The regen script detects this by inspecting the AST and falls back to distribution mode whenever there's any source of randomness.

**Estimated total fixture-test runtime: ~30–40 seconds**, on top of the existing ~2 seconds of unit tests. Larger than the previous estimate because the `noise_sweeps/` category adds 30–50 statistical fixtures at 4096 shots each.

If runtime later becomes a problem, the simplest knobs are: dropping the `d=7` codes fixtures, reducing `noise_sweeps/`'s shot count to 2048, or splitting the suite into core/extended tiers. We aren't tiering yet; revisit if the budget is exceeded in practice.

## Regen Tool

Lives at `crates/ppvm-stim/tests/regen-stim/` — its own uv-managed Python project with `pyproject.toml`, `uv.lock`, and `.python-version` for reproducibility. Standalone from `ppvm-python/` deliberately: `regen-stim` only needs `stim`; `ppvm-python` doesn't need `stim` as a runtime dep. Keeping them separate avoids contaminating the ppvm wheel's transitive dep graph and lets `regen-stim` pin its own Python version (currently 3.12).

Direct dep: `stim>=1.15.0`. Indirectly the regen tool also needs `ppvm` installed to drive ppvm's sampler — added at implementation time as a path source pointing back at `../../../../../ppvm-python` (or via PYTHONPATH; the implementation plan picks the cleanest option).

Entry point: a CLI (`regen-stim` script in `pyproject.toml`'s `[project.scripts]`) with subcommands:

```bash
regen-stim codes      [--distance 3,5,7] [--rounds 1,3,5] [--noise 0,0.001,0.01]
regen-stim noise-sweeps  [--channel <name>]
regen-stim dialect    [--name <name>]
regen-stim random     [--seed 0] [--count 25]
regen-stim unsupported [--name <name>]
regen-stim refresh    <fixture-path>      # re-record one fixture (e.g. after bumping Stim)
regen-stim verify     <fixture-path>      # re-run cross-check, error on mismatch, write nothing
regen-stim all                              # regenerate everything
```

Run from the regen-stim folder via `uv run regen-stim …` (uv handles venv + dep resolution).

Output: writes `.stim` and `.expected.json` files into `crates/ppvm-stim/tests/data/<category>/` (the corpus directory). The path from regen-stim to corpus is `../data/`, hard-coded in the script (we're committing both, they live next to each other in the repo).

Loop for distribution-mode fixtures:

```
1. Generate or load .stim source.
2. Run Stim with stim_num_shots (default 10_000) at stim_seed=0
   → reference distribution (per-bit means).
3. For ppvm_seed in [0..32]:
       Run ppvm with num_shots (per-tier default).
       Compute empirical per-bit means.
       For each bit:
           sigma = sqrt(stim_mean*(1-stim_mean) / num_shots)  if stim_mean ∈ (0,1)
                 = sqrt(1.0 / num_shots)                       if stim_mean ∈ {0, 1} (worst-case)
           if |ppvm_mean - stim_mean| > tolerance_sigma * sigma:
               break to next seed
       else:
           commit this seed; record ppvm_bit_means.
   If no seed passes, error loudly. Don't commit.
4. Write .expected.json with all metadata.
```

Loop for unsupported-mode fixtures: same as distribution-mode through step 2, then write JSON without ppvm fields, with `mode: "unsupported"` and `awaiting_phase2_instruction` set.

Loop for deterministic-mode fixtures: just run ppvm once at `ppvm_seed=0`, record the bitstring.

Stim version is recorded in every JSON so a future Stim release can be detected as the cause of any mismatch.

## Test Harness Changes

`tests/stim_corpus.rs` (today: 95 lines, 2 tests) grows to support the new modes. The two existing tests (`corpus_table_covers_every_file`, `corpus_obeys_expectations`) generalize:

- `corpus_table_covers_every_file` walks the entire `data/` tree (recursively into subdirs), asserts every `.stim` has a matching `.expected.json` and vice versa.
- `corpus_obeys_expectations` walks the same tree, parses the JSON, dispatches to one of three test paths based on `mode`.

The harness uses `serde_json` for JSON parsing — a new dev-dep on `ppvm-stim`. (We considered TOML to avoid the dep; JSON wins because the regen script is Python and `json` is the path of least friction there.)

The legacy `Expect` enum (`Ok` / `NormalizeUnsupported` / `ParseFails`) becomes obsolete and is removed; the JSON's `mode` field is the source of truth.

## Snapshot Tests

Separate test file: `crates/ppvm-stim/tests/parser_snapshots.rs`. Uses `insta` (already a dev-dep transitively but needs to be added to ppvm-stim's `Cargo.toml`).

About 15 representative programs, one per AST shape:
- Bare gate (`H 0`).
- Tagged gate (`S[T] 0`).
- Tagged gate with named params (`I[R_X(theta=0.5*pi)] 0`).
- Multi-tag (`S[T,debug] 0`).
- Args-only noise (`DEPOLARIZE1(0.5) 0`).
- Tag + args (`I_ERROR[correlated_loss](0.1, 0.2, 0.3) 0 1`).
- Single-target measurement (`M 0`).
- Multi-target measurement with noise (`M(0.001) 0 1 2`).
- Annotation with rec target (`DETECTOR rec[-1]`).
- Annotation with args (`OBSERVABLE_INCLUDE(0)`).
- Empty REPEAT body.
- Multi-instruction REPEAT.
- Nested REPEAT.
- Comment-heavy program.
- Whitespace-stress program.

Each snapshot captures the parsed `Program` debug-formatted via `insta::assert_debug_snapshot!`. Re-running flags any unexpected AST drift cheaply.

## Implementation Order

The corpus expansion is a substantial chunk of work; phasing it lets the test improvements land incrementally. The `regen-stim/` Python project is already scaffolded (uv-managed, `stim>=1.15.0`); steps 5+ build on top of it.

1. **Harness rewrite** — extend `stim_corpus.rs` to walk subdirectories recursively, dispatch on JSON `mode`, support all three modes; add `serde_json` dev-dep on `ppvm-stim`; convert the existing 8 fixtures to the new schema (most will be deterministic mode).
2. **Edge cases corpus** — fully hand-written, no Stim dependency, ~20-25 fixtures. Ships value immediately.
3. **Noise channels corpus** — hand-written corner cases (p=0.0/1.0 boundaries, ordering invariants); Stim cross-check the non-trivial ones at regen time. ~8-10 fixtures.
4. **Snapshot tests** — `tests/parser_snapshots.rs` with `insta`, ~15 representative programs. Small, independent, useful immediately.
5. **regen-stim CLI scaffolding** — fill in the empty `regen-stim/` project: shared library code (Stim invocation, ppvm invocation, JSON read/write, seed search loop, tolerance check), `[project.scripts]` entry, common subcommand framework. No fixtures generated yet — just the plumbing.
6. **Generated codes corpus** — `regen-stim codes` with full sweeps. Bulk of the new fixtures (~80-100).
7. **Generated noise_sweeps corpus** — `regen-stim noise-sweeps` with per-channel sweeps (~30-50 fixtures).
8. **Generated dialect corpus** — Python templates for ppvm-specific instructions; deterministic-mode predominantly. (~10-15 fixtures)
9. **Generated random corpus** — random-walk generator with three regimes (clifford-only / +noise / +readout). (~30-40 fixtures)
10. **Unsupported corpus** — auto-routed from `regen-stim codes` natural emission + per-instruction templates. ~20-25 fixtures.
11. **Documentation** — README per category subdir, top-level `tests/data/README.md` describing the regen workflow, regen-stim's own README.

Each step ends with a passing `cargo test`. Steps 5–10 require `stim` installed in the regen-stim venv (which `uv sync` from `regen-stim/` handles automatically).

## Out of Scope

- **Property-based testing.** Considered for the parser (round-trip random AST → string → AST), considered deferred. If we revisit, the natural place is a sibling `tests/parser_proptest.rs` with `proptest` as a dev-dep.
- **Performance regression assertions.** `criterion` benches in `crates/ppvm-tableau/benches/tableau-msd-stim.rs` already cover the perf side; thresholding is bench-suite scope.
- **Stim version drift handling.** When Stim updates and changes its RNG semantics, regen on every distribution-mode fixture breaks loudly. The fix is "rerun regen, commit refreshed JSONs". We don't try to be RNG-version-agnostic.
- **Cross-language tests.** Python wrapper tests (`ppvm-python/test/`) stay focused on Python API surface; the corpus drives the Rust side only. The Python tests will pick up correctness improvements transitively via the wheel build.

## Open Questions for Phase-2 Test Work

- When phase-2 implements an unsupported instruction (e.g. SWAP), `uv run regen-stim refresh ../data/unsupported/swap.stim` flips it to a passing fixture. We need to commit to running this on every unsupported instruction phase-2 implements — easy to forget. Worth adding a CI check that asserts every `unsupported/` JSON's `awaiting_phase2_instruction` actually maps to an unsupported variant (would catch a stale mode after phase-2 lifted the restriction).
- A future "fuzz" tier that generates malformed input and asserts the parser never panics — independent of the corpus, lives next to `parser_proptest.rs` if we adopt proptest.
- Cross-check fidelity at higher `stim_num_shots` (e.g. 100k) is cheap but produces slightly more accurate `stim_bit_means`. We default to 10k; can revisit if cross-check tolerance starts producing false positives.
