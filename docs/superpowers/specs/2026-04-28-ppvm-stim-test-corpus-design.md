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
- Cross-checking ppvm-specific dialect instructions (`I[R_X(theta=…)]`, `S[T]`, `I_ERROR[loss]`) against Stim — Stim doesn't simulate them. Those use deterministic-mode reference outputs derived by hand or by running ppvm itself.
- Property-based testing with `proptest` / `quickcheck`. Considered and deferred — useful but a sizable infra add and orthogonal to corpus expansion. Track separately.
- CI tiering. We aren't splitting fast vs slow tests yet. If runtime becomes a problem we can revisit.
- Performance regression benchmarks. Already partially covered by `tableau-msd-stim`; threshold work is bench-suite scope.

## Approach in One Picture

```
                 ┌─────────────────────┐
                 │ regen.py            │  Run by humans with `pip install stim`.
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

For circuits with no noise instructions. The bitstring is unique given the circuit; one shot suffices.

```jsonc
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true, false, false, true, null, true]
  // length matches the circuit's measurement count
  // null entries represent measurements on lost qubits
}
```

Test action: parse + normalize + execute once with `ppvm_seed`; assert returned `Vec<Option<bool>>` equals `bitstring` (mapping `null` → `None`).

### Mode 2: distribution

For circuits with noise. Stim cross-check is meaningful here.

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
  "stim_version": "1.13.0"
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
  "stim_version": "1.13.0"
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
| `generated/codes/` | `stim gen` with sweeps over `(code, distance, rounds, basis, noise)` | 80–100 | Real-world fault-tolerance circuits with realistic noise |
| `generated/dialect/` | Programmatic via Python + Stim's circuit-construction API | 10–15 | ppvm-specific dialect (`I[R_*]`, `S[T]`, `I_ERROR[loss]`); deterministic-mode only |
| `generated/random/` | Our own random-walk Python generator | 20–30 | Random sequences of supported instructions; broad shape coverage |
| `unsupported/` | `stim gen` natural emission + programmatic templates | 20–25 | One fixture per phase-1-unsupported instruction; flips in phase-2 |
| `edge_cases/` | Hand-written | 20–25 | Empty programs, deeply nested REPEAT, every Pi-expression form, every tag form, dense/sparse measurement patterns, comment/whitespace stress |
| `noise_channels/` | Hand-written | 10–12 | One fixture per supported `NoiseKind`, plus a few channel combinations; statistical-power-sensitive |

**Total: ~160–200 committed fixtures** (`.stim` + `.expected.json` pairs).

### `generated/codes/` — `stim gen` sweeps

The regen script invokes `stim gen` for combinations of:

- `--code` ∈ {`surface_code`, `repetition_code`, `color_code`}
- `--task` ∈ all task variants the code supports (e.g. `surface_code` → `unrotated_memory_x`, `unrotated_memory_z`, `rotated_memory_x`, `rotated_memory_z`)
- `--distance` ∈ {3, 5, 7}
- `--rounds` ∈ {1, 3, 5}
- noise flags: either none, or `--after_clifford_depolarization=p --before_round_data_depolarization=p --before_measure_flip_probability=p --after_reset_flip_probability=p` for `p` ∈ {0.001, 0.01}

Some `(code, task)` combinations naturally emit phase-1-unsupported instructions (e.g. `unrotated_memory_z` uses `MR`; certain transversal tasks use `SWAP`; stabilizer-measurement tasks use `MPP`). Those land in `unsupported/` instead of `generated/codes/` — the regen script auto-routes by inspecting the emitted source.

### `generated/dialect/` — ppvm-specific instructions

Stim itself cannot simulate `I[R_X(theta=…)]`, `S[T]`, or `I_ERROR[loss]`, so cross-check via Stim is impossible. These fixtures use **deterministic mode**: the regen script runs ppvm itself and records the bitstring. Verification is structural (the circuit produces the bitstring derivable from inspection) rather than oracle-driven. ~15 hand-templated programmatic fixtures: each combines a few rotation/noise instructions in a small circuit with a known expected outcome.

### `generated/random/` — random-walk programs

A Python script that emits random sequences of supported instructions. Parameterized by:
- Number of qubits ∈ {2, 4, 8, 16}
- Number of instructions ∈ {10, 50, 200}
- Mix ratios (Clifford : non-Clifford : noise : measurement : annotation)
- RNG seed

For each combination the script emits a `.stim` file. Cross-check via Stim works for any sequence using only Stim-supported instructions; falls back to deterministic mode for sequences mixing in dialect-specific ones.

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

### `noise_channels/` — statistical-power-sensitive

One small (1–3 qubit) fixture per supported `NoiseKind`:
- `Depolarize1` at p=0.5 over a single qubit.
- `Depolarize2` at p=0.5 over a single pair.
- `PauliChannel1` with skewed probabilities `[0.1, 0.2, 0.3]`.
- `PauliChannel2` with all 15 args set.
- `XError` / `YError` / `ZError` at extreme probabilities (`p ∈ {0.0, 0.5, 1.0}` to lock down the boundary cases).
- `I_ERROR[loss]` at p=0.5 → measurement must produce mixed `Some(_)`/`None`.
- `I_ERROR[correlated_loss]` at `(p_x=0.3, p_y=0.2, p_z=0.1)`.
- `M(p)` / `MR(p)` with `p ∈ {0.0, 0.5, 1.0}` to verify the readout-noise path that the recent `MZ(p)` bug fix exercised.

Higher `num_shots` (4096) than other categories because the circuits are tiny and statistical drift in noise is the failure mode this category exists to detect.

## Per-Category Shot Defaults

Tests bit-exact-compare ppvm's empirical means against committed `ppvm_bit_means`, so `num_shots` only affects test runtime — not test precision. The defaults below balance regen-time cross-check tightness against test-time speed:

| Category | num_shots | per-shot cost (release) | per-fixture wall time |
|---|---|---|---|
| `edge_cases/` (mostly deterministic) | 1 | <0.1 ms | <1 ms |
| `generated/codes/` distance ≤ 5 | 256 | ~0.5 ms | ~125 ms |
| `generated/codes/` distance = 7 | 64 | ~5 ms | ~325 ms |
| `generated/dialect/` no noise | 1 | <0.1 ms | <1 ms |
| `generated/dialect/` with noise | 256 | ~0.5 ms | ~125 ms |
| `generated/random/` | 128 | ~1 ms | ~130 ms |
| `noise_channels/` | 4096 | ~0.05 ms | ~200 ms |
| `unsupported/` | n/a | n/a | <1 ms (parse + normalize only) |

**Estimated total fixture-test runtime: ~25–30 seconds**, on top of the existing ~2 seconds of unit tests.

If runtime later becomes a problem, the simplest knob is dropping `d=7` fixtures or reducing their `num_shots` to 32. We aren't tiering the suite into core/extended yet; revisit if the budget is exceeded in practice.

## Regen Script

Single Python entry point at `crates/ppvm-stim/tests/data/regen.py`. Subcommands:

```bash
regen.py codes      [--distance 3,5,7] [--rounds 1,3,5] [--noise 0,0.001,0.01]
regen.py dialect    [--name <name>]
regen.py random     [--seed 0] [--count 25]
regen.py unsupported [--name <name>]
regen.py refresh    <fixture-path>     # re-record one fixture (e.g. after bumping Stim)
regen.py verify     <fixture-path>     # re-run cross-check, error on mismatch, write nothing
regen.py all                            # regenerate everything
```

Inputs the script needs:
- `stim` Python package installed in the regen environment.
- A working ppvm install (the script imports `ppvm` to drive ppvm's sampler).

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

The corpus expansion is a substantial chunk of work; phasing it lets the test improvements land incrementally.

1. **Harness rewrite** — extend `stim_corpus.rs` to handle the three modes; add `serde_json` dep; keep the existing 8 fixtures working in their new schema.
2. **Edge cases corpus** — fully hand-written, no Stim dependency, ~25 fixtures. Ships value immediately.
3. **Noise channels corpus** — fully hand-written for deterministic-mode corner cases (p=0.0, p=1.0). Stim cross-check the rest at regen time.
4. **Snapshot tests** — small, independent, useful immediately.
5. **`regen.py` skeleton** — supports all three modes, but only wired to a few fixtures initially.
6. **Generated codes corpus** — `regen.py codes` with full sweeps. Bulk of the new fixtures.
7. **Generated dialect corpus** — Python templates for ppvm-specific instructions.
8. **Generated random corpus** — random-walk generator.
9. **Unsupported corpus** — auto-routed from `regen.py codes` (natural emission) plus per-instruction templates.
10. **Documentation** — README per category, top-level `tests/data/README.md` describing the regen workflow.

Each step ends with a passing `cargo test`. Steps 6–9 require Stim installed.

## Out of Scope

- **Property-based testing.** Considered for the parser (round-trip random AST → string → AST), considered deferred. If we revisit, the natural place is a sibling `tests/parser_proptest.rs` with `proptest` as a dev-dep.
- **Performance regression assertions.** `criterion` benches in `crates/ppvm-tableau/benches/tableau-msd-stim.rs` already cover the perf side; thresholding is bench-suite scope.
- **Stim version drift handling.** When Stim updates and changes its RNG semantics, regen on every distribution-mode fixture breaks loudly. The fix is "rerun regen, commit refreshed JSONs". We don't try to be RNG-version-agnostic.
- **Cross-language tests.** Python wrapper tests (`ppvm-python/test/`) stay focused on Python API surface; the corpus drives the Rust side only. The Python tests will pick up correctness improvements transitively via the wheel build.

## Open Questions for Phase-2 Test Work

- When phase-2 implements an unsupported instruction (e.g. SWAP), `regen.py refresh tests/data/unsupported/swap.stim` flips it to a passing fixture. We need to commit to running this on every unsupported instruction phase-2 implements — easy to forget. Worth adding a CI check that asserts every `unsupported/` JSON's `awaiting_phase2_instruction` actually maps to an unsupported variant (would catch a stale mode after phase-2 lifted the restriction).
- A future "fuzz" tier that generates malformed input and asserts the parser never panics — independent of the corpus, lives next to `parser_proptest.rs` if we adopt proptest.
- Cross-check fidelity at higher `stim_num_shots` (e.g. 100k) is cheap but produces slightly more accurate `stim_bit_means`. We default to 10k; can revisit if cross-check tolerance starts producing false positives.
