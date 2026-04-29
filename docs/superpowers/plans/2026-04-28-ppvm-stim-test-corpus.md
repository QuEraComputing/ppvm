# ppvm-stim Test Corpus Expansion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expand `ppvm-stim`'s test suite from 8 hand-written fixtures into a substantive corpus (~200–250 fixtures across 7 categories) cross-checked against Stim's reference simulator at regen time, plus a parser snapshot test suite.

**Architecture:** Two-tier verification. A `regen-stim` Python CLI (uv-managed; not run in CI) calls Stim and ppvm to produce committed `.stim` + `.expected.json` pairs under `crates/ppvm-stim/tests/data/<category>/`. Tests load those JSONs and bit-exact-compare ppvm's output against the committed reference; the Stim cross-check happens before each fixture is committed, not at test time. Fixtures using phase-1-unsupported instructions live in their own `mode: "unsupported"` schema and flip to "supported" automatically when phase-2 lands the gate.

**Tech Stack:** Rust (existing `ppvm-stim` crate; new dev-deps: `serde_json`, `walkdir`, `insta`). Python 3.12 in `crates/ppvm-stim/tests/regen-stim/` (`stim>=1.15.0`, `ppvm` as path-source). JSON for fixture metadata. `insta` for parser snapshots.

**Spec:** `docs/superpowers/specs/2026-04-28-ppvm-stim-test-corpus-design.md`. Read this for category details, sweep axes, and per-category counts. The plan implements the spec's "Implementation Order" section.

---

## File Structure

**Crate (Rust) changes:**
- Modify: `crates/ppvm-stim/Cargo.toml` — add `serde_json`, `walkdir`, `insta`, `serde` dev-deps.
- Modify: `crates/ppvm-stim/tests/stim_corpus.rs` — full rewrite; recursive walk, JSON-driven mode dispatch.
- Create: `crates/ppvm-stim/tests/parser_snapshots.rs` — `insta`-driven snapshots for ~15 representative programs.
- Modify: each existing `crates/ppvm-stim/tests/data/*.stim` — add a sibling `<name>.expected.json`. Move existing top-level fixtures into category subdirs (`edge_cases/` or `unsupported/`).
- Create: `crates/ppvm-stim/tests/data/<category>/` subdirectory tree (categories: `edge_cases/`, `noise_channels/`, `unsupported/`, `generated/codes/`, `generated/noise_sweeps/`, `generated/dialect/`, `generated/random/`).
- Modify: `crates/ppvm-stim/tests/data/README.md` — describe the new schema and regen workflow.

**Regen tool (Python) changes:**
- Modify: `crates/ppvm-stim/tests/regen-stim/pyproject.toml` — add `[project.scripts]`, ppvm path source, dev tooling.
- Modify: `crates/ppvm-stim/tests/regen-stim/README.md` — fill in (currently empty).
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/__init__.py`
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/cli.py` — entry point + subcommand dispatch.
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/core.py` — Stim invocation, ppvm invocation, JSON read/write, seed-search loop, tolerance check.
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/codes.py` — `stim gen` sweep generator (Task 7).
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/noise_sweeps.py` — per-channel sweep generator (Task 8).
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/dialect.py` — ppvm-specific dialect templates (Task 9).
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/random_walk.py` — random program generator (Task 10).
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/unsupported.py` — phase-1-unsupported fixtures (Task 11).
- Create: `crates/ppvm-stim/tests/regen-stim/test/test_core.py` — pytest tests for the seed-search loop and tolerance math.

**Documentation:**
- Modify: `crates/ppvm-stim/tests/data/README.md` (top-level corpus README).
- Create: `crates/ppvm-stim/tests/data/<category>/README.md` for each category (provenance, regen invocation).

**Files NOT changed:**
- `crates/ppvm-stim/src/**` — the corpus expansion is test-only; Phase-1 source is frozen for this plan.
- `ppvm-python/**` — the corpus drives the Rust side. Python wrapper tests stay in their own world.
- CI configuration — no tiering; harness runs as part of the existing `cargo test` in CI.

---

### Task 1: Harness rewrite — JSON-driven schema, three modes, recursive walk

**Files:**
- Modify: `crates/ppvm-stim/Cargo.toml` — add dev-deps.
- Modify: `crates/ppvm-stim/tests/stim_corpus.rs` — full rewrite.

The new harness walks the entire `data/` tree recursively, requires every `.stim` file to have a sibling `<name>.expected.json`, parses the JSON, and dispatches on the `mode` field. We need three JSON modes (`deterministic`, `distribution`, `unsupported`) plus a strict "every `.stim` has a JSON, every JSON has a `.stim`" pairing assertion.

- [ ] **Step 1: Add dev-deps to `Cargo.toml`**

Edit `crates/ppvm-stim/Cargo.toml`. Replace the existing `[dev-dependencies]` block with:

```toml
[dev-dependencies]
# Used by tests/data fixtures and statistical assertions.
rand = "0.10"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
walkdir = "2.5"
insta = { version = "1.43.2", features = ["yaml"] }
```

Run:

```bash
cargo build -p ppvm-stim --tests
```

Expected: build succeeds; new deps resolve.

- [ ] **Step 2: Write a failing test that exercises the new harness**

Replace the *entire* contents of `crates/ppvm-stim/tests/stim_corpus.rs` with the version below. The two existing tests (`corpus_table_covers_every_file`, `corpus_obeys_expectations`) are removed because the table-driven `CASES` constant is obsolete — the JSON files are now the source of truth.

```rust
//! Corpus harness for `crates/ppvm-stim/tests/data/`.
//!
//! Walks the entire `data/` tree recursively. Every `.stim` file must have a
//! sibling `<name>.expected.json`. The JSON's `mode` field selects one of
//! three test paths (deterministic, distribution, unsupported). Bit-exact
//! comparison against committed `ppvm_bit_means` is the test signal — the
//! Stim cross-check happens at regen time, not here.
//!
//! See `docs/superpowers/specs/2026-04-28-ppvm-stim-test-corpus-design.md`
//! for the schema reference and the rationale.

use std::path::{Path, PathBuf};

use ppvm_runtime::config::indexmap::ByteFxHashF64;
use ppvm_stim::{NormalizeError, execute, normalize, parse, sample};
use ppvm_tableau::prelude::*;
use serde::Deserialize;
use walkdir::WalkDir;

type Tab = GeneralizedTableau<ByteFxHashF64<8>, usize>;

/// 64-qubit tableau is enough for every fixture we plan to commit
/// (max distance-7 surface code is ~50 data qubits + ancillas).
const N_QUBITS: usize = 64;
const COEFF_THRESHOLD: f64 = 1e-10;

#[derive(Debug, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case")]
enum Expected {
    Deterministic {
        ppvm_seed: u64,
        bitstring: Vec<bool>,
    },
    Distribution {
        num_shots: usize,
        ppvm_seed: u64,
        ppvm_bit_means: Vec<f64>,
        // Documentation-only fields; harness does not use them.
        #[serde(default)]
        stim_seed: Option<u64>,
        #[serde(default)]
        stim_num_shots: Option<usize>,
        #[serde(default)]
        stim_bit_means: Option<Vec<f64>>,
        #[serde(default)]
        tolerance_sigma_at_regen: Option<f64>,
        #[serde(default)]
        stim_version: Option<String>,
    },
    Unsupported {
        awaiting_phase2_instruction: String,
        // Pre-recorded for phase-2 flip; harness does not use them.
        #[serde(default)]
        stim_seed: Option<u64>,
        #[serde(default)]
        stim_num_shots: Option<usize>,
        #[serde(default)]
        stim_bit_means: Option<Vec<f64>>,
        #[serde(default)]
        tolerance_sigma_at_regen: Option<f64>,
        #[serde(default)]
        stim_version: Option<String>,
    },
}

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
}

fn relpath(p: &Path) -> String {
    p.strip_prefix(data_dir())
        .unwrap_or(p)
        .to_string_lossy()
        .into_owned()
}

/// Walk the corpus, returning every `.stim` path under `data/`.
fn collect_stim_files() -> Vec<PathBuf> {
    let mut out: Vec<PathBuf> = WalkDir::new(data_dir())
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .map(|e| e.into_path())
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("stim"))
        .collect();
    out.sort();
    out
}

fn collect_json_files() -> Vec<PathBuf> {
    let mut out: Vec<PathBuf> = WalkDir::new(data_dir())
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .map(|e| e.into_path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.ends_with(".expected.json"))
        })
        .collect();
    out.sort();
    out
}

#[test]
fn corpus_pairs_every_stim_with_expected_json() {
    let stim_files = collect_stim_files();
    let json_files = collect_json_files();

    // Every .stim must have a sibling <name>.expected.json.
    let mut missing_json: Vec<String> = Vec::new();
    for s in &stim_files {
        let expected = s.with_extension("expected.json");
        if !expected.exists() {
            missing_json.push(relpath(s));
        }
    }

    // Every .expected.json must have a sibling <name>.stim.
    let mut missing_stim: Vec<String> = Vec::new();
    for j in &json_files {
        // <name>.expected.json → <name>.stim
        let stem = j
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap()
            .strip_suffix(".expected.json")
            .unwrap();
        let stim = j.with_file_name(format!("{stem}.stim"));
        if !stim.exists() {
            missing_stim.push(relpath(j));
        }
    }

    assert!(
        missing_json.is_empty(),
        "fixtures missing .expected.json: {missing_json:?}"
    );
    assert!(
        missing_stim.is_empty(),
        "expected JSONs missing .stim: {missing_stim:?}"
    );
    assert!(
        !stim_files.is_empty(),
        "no fixtures found under {}; corpus must have at least one fixture",
        data_dir().display()
    );
}

#[test]
fn corpus_obeys_expectations() {
    let stim_files = collect_stim_files();
    let mut failures: Vec<String> = Vec::new();

    for stim_path in &stim_files {
        let label = relpath(stim_path);
        let json_path = stim_path.with_extension("expected.json");

        let src = match std::fs::read_to_string(stim_path) {
            Ok(s) => s,
            Err(e) => {
                failures.push(format!("{label}: read .stim failed: {e}"));
                continue;
            }
        };
        let json_str = match std::fs::read_to_string(&json_path) {
            Ok(s) => s,
            Err(e) => {
                failures.push(format!("{label}: read .expected.json failed: {e}"));
                continue;
            }
        };
        let expected: Expected = match serde_json::from_str(&json_str) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("{label}: malformed expected.json: {e}"));
                continue;
            }
        };

        if let Err(msg) = run_one(&label, &src, &expected) {
            failures.push(msg);
        }
    }

    assert!(
        failures.is_empty(),
        "corpus failures:\n  - {}",
        failures.join("\n  - ")
    );
}

fn run_one(label: &str, src: &str, expected: &Expected) -> Result<(), String> {
    let prog = parse(src).map_err(|e| format!("{label}: parse failed: {e}"))?;

    match expected {
        Expected::Unsupported {
            awaiting_phase2_instruction,
            ..
        } => match normalize::to_tableau(&prog) {
            Err(NormalizeError::Unsupported { name, .. }) => {
                if name != *awaiting_phase2_instruction {
                    return Err(format!(
                        "{label}: expected Unsupported({awaiting_phase2_instruction}), got Unsupported({name})"
                    ));
                }
                Ok(())
            }
            Err(other) => Err(format!(
                "{label}: expected Unsupported, got {other:?}"
            )),
            Ok(_) => Err(format!(
                "{label}: expected Unsupported({awaiting_phase2_instruction}), but normalize succeeded"
            )),
        },

        Expected::Deterministic {
            ppvm_seed,
            bitstring,
        } => {
            let tprog = normalize::to_tableau(&prog)
                .map_err(|e| format!("{label}: normalize failed: {e}"))?;
            let mut tab: Tab =
                GeneralizedTableau::new_with_seed(N_QUBITS, COEFF_THRESHOLD, *ppvm_seed);
            let results =
                execute(&tprog, &mut tab).map_err(|e| format!("{label}: execute failed: {e:?}"))?;
            if results.len() != bitstring.len() {
                return Err(format!(
                    "{label}: bitstring length mismatch: got {} bits, expected {}",
                    results.len(),
                    bitstring.len()
                ));
            }
            for (i, (got, want)) in results.iter().zip(bitstring).enumerate() {
                match got {
                    Some(b) if b == want => {}
                    Some(b) => {
                        return Err(format!(
                            "{label}: bit {i}: got {b}, expected {want}"
                        ));
                    }
                    None => {
                        return Err(format!(
                            "{label}: bit {i}: got None (loss), expected {want} — \
                             deterministic-mode fixtures must not contain loss"
                        ));
                    }
                }
            }
            Ok(())
        }

        Expected::Distribution {
            num_shots,
            ppvm_seed,
            ppvm_bit_means,
            ..
        } => {
            let tprog = normalize::to_tableau(&prog)
                .map_err(|e| format!("{label}: normalize failed: {e}"))?;
            // Single seed reused across all shots — sample() builds a fresh
            // tableau per shot, but seeds them from a counter starting at
            // ppvm_seed. We rely on `new_with_seed` advancing the seed by 1
            // per call so ppvm's internal RNG stream is fully determined by
            // (program, ppvm_seed, num_shots).
            let mut next_seed = *ppvm_seed;
            let shots = sample(&tprog, *num_shots, || {
                let s = next_seed;
                next_seed = next_seed.wrapping_add(1);
                GeneralizedTableau::<ByteFxHashF64<8>, usize>::new_with_seed(
                    N_QUBITS,
                    COEFF_THRESHOLD,
                    s,
                )
            })
            .map_err(|e| format!("{label}: sample failed: {e:?}"))?;

            if shots.is_empty() {
                return Err(format!("{label}: sample returned 0 shots"));
            }
            let bits = shots[0].len();
            if bits != ppvm_bit_means.len() {
                return Err(format!(
                    "{label}: bit count mismatch: got {bits}, expected {} per ppvm_bit_means",
                    ppvm_bit_means.len()
                ));
            }
            for shot in &shots {
                if shot.len() != bits {
                    return Err(format!(
                        "{label}: shot bit-count drift: {} vs {bits}",
                        shot.len()
                    ));
                }
                if shot.iter().any(|m| m.is_none()) {
                    return Err(format!(
                        "{label}: distribution-mode fixture must not contain loss \
                         (got None measurement); use Mode::Unsupported or \
                         exclude loss per spec Non-Goals"
                    ));
                }
            }

            // Per-bit empirical mean across shots; bit-exact f64 compare.
            for bit in 0..bits {
                let sum: f64 = shots
                    .iter()
                    .map(|shot| if shot[bit].unwrap() { 1.0 } else { 0.0 })
                    .sum();
                let mean = sum / (*num_shots as f64);
                let want = ppvm_bit_means[bit];
                if mean.to_bits() != want.to_bits() {
                    return Err(format!(
                        "{label}: bit {bit}: ppvm_bit_means drift: got {mean}, \
                         expected {want} (bit-exact f64 compare)"
                    ));
                }
            }
            Ok(())
        }
    }
}
```

- [ ] **Step 3: Run the harness — expect `pairs_every_stim_with_expected_json` to fail**

Run:

```bash
cargo test -p ppvm-stim --test stim_corpus
```

Expected: `corpus_pairs_every_stim_with_expected_json` fails because none of the existing 8 `.stim` files have a `.expected.json` sibling. Output should list all 8 by name. `corpus_obeys_expectations` may also fail or be skipped depending on iteration order — both will be fixed by Task 2.

- [ ] **Step 4: Commit the harness rewrite**

```bash
git add crates/ppvm-stim/Cargo.toml crates/ppvm-stim/tests/stim_corpus.rs
git commit -m "test(stim): rewrite corpus harness for JSON-driven mode dispatch"
```

(Cargo.lock will also have changed. Stage it too if needed: `git add Cargo.lock`.)

---

### Task 2: Convert existing 8 fixtures to the new schema

**Files:**
- Move: `crates/ppvm-stim/tests/data/{x_only,bell_pair,ghz,repeat_block,repetition_code_d3_r3,depolarize_smoke}.stim` → `crates/ppvm-stim/tests/data/edge_cases/`.
- Move: `crates/ppvm-stim/tests/data/{swap_unsupported,mx_unsupported}.stim` → `crates/ppvm-stim/tests/data/unsupported/`.
- Create: a sibling `<name>.expected.json` for each.
- Modify: `crates/ppvm-stim/tests/data/README.md` — update provenance after moves.

These fixtures are already written; we only need to relocate them into category subdirs and write their JSON expectation files. The classification:
- Pure-deterministic (no randomness, no noise): `x_only`, `bell_pair`, `ghz`, `repetition_code_d3_r3`, `repeat_block` — all have predictable bitstrings.
- Distribution: `depolarize_smoke` — has noise, needs `ppvm_bit_means`.
- Unsupported: `swap_unsupported` (SWAP), `mx_unsupported` (MX).

We don't need regen-stim for this task — we generate `ppvm_bit_means` for `depolarize_smoke` by running the harness once (Step 6 below) and pasting the output.

- [ ] **Step 1: Move existing fixtures into category subdirs**

```bash
mkdir -p crates/ppvm-stim/tests/data/edge_cases crates/ppvm-stim/tests/data/unsupported
git mv crates/ppvm-stim/tests/data/x_only.stim crates/ppvm-stim/tests/data/edge_cases/x_only.stim
git mv crates/ppvm-stim/tests/data/bell_pair.stim crates/ppvm-stim/tests/data/edge_cases/bell_pair.stim
git mv crates/ppvm-stim/tests/data/ghz.stim crates/ppvm-stim/tests/data/edge_cases/ghz.stim
git mv crates/ppvm-stim/tests/data/repeat_block.stim crates/ppvm-stim/tests/data/edge_cases/repeat_block.stim
git mv crates/ppvm-stim/tests/data/repetition_code_d3_r3.stim crates/ppvm-stim/tests/data/edge_cases/repetition_code_d3_r3.stim
git mv crates/ppvm-stim/tests/data/depolarize_smoke.stim crates/ppvm-stim/tests/data/edge_cases/depolarize_smoke.stim
git mv crates/ppvm-stim/tests/data/swap_unsupported.stim crates/ppvm-stim/tests/data/unsupported/swap_unsupported.stim
git mv crates/ppvm-stim/tests/data/mx_unsupported.stim crates/ppvm-stim/tests/data/unsupported/mx_unsupported.stim
```

- [ ] **Step 2: Write `expected.json` for the deterministic fixtures**

Both `x_only.stim` and `bell_pair.stim` are pure Clifford preparations with measurements that have unique outcomes. Verify by reading each `.stim`:
- `x_only.stim`: `X 0; M 0` — measures 1.
- `bell_pair.stim`: `H 0; CX 0 1; M 0 1` — has measurement randomness, **not** deterministic. Goes in distribution mode.
- `ghz.stim`: `H 0; CX 0 1; CX 1 2; CX 2 3; M 0 1 2 3` — has measurement randomness, distribution mode.
- `repeat_block.stim`: `REPEAT 5 { X 0 }; M 0` — measures 1 (5 X gates → flipped state).
- `repetition_code_d3_r3.stim`: noisy reset behavior; the file uses `MR 3 4` which has measurement randomness on initial reset state. Distribution mode. **Wait — `R` followed by `MR` on a freshly-reset qubit is deterministic 0**. The file's measurement count is computed by walking the program; let's count them precisely below.

Let me re-read each file before writing JSON.

```bash
cat crates/ppvm-stim/tests/data/edge_cases/repetition_code_d3_r3.stim
```

The program is:
```
R 0 1 2 3 4
REPEAT 3 {
    CX 0 3
    CX 1 3
    CX 1 4
    CX 2 4
    MR 3 4
    DETECTOR rec[-1]
    DETECTOR rec[-2]
    TICK
}
M 0 1 2
OBSERVABLE_INCLUDE(0) rec[-1]
```

After `R 0 1 2 3 4`, all qubits are in |0⟩. Each round: a few CX onto qubits 3 and 4 (which start in |0⟩), then MR (measure-and-reset). All five qubits start in |0⟩, so the parity sums on qubits 3 and 4 stay deterministic 0 across all rounds. Final `M 0 1 2` measures three reset qubits → 0 0 0. Total measurements = 3 rounds × 2 (MR per round) + 3 (final M) = 9. All zeros.

`x_only.stim` — measurements: 1. Bitstring: `[true]`.
`repeat_block.stim` — measurements: 1. Bitstring: `[true]` (5 X = X = flip, then M 0).
`repetition_code_d3_r3.stim` — measurements: 9. Bitstring: `[false; 9]`.

`bell_pair.stim` and `ghz.stim` and `depolarize_smoke.stim` are distribution mode.

Write `crates/ppvm-stim/tests/data/edge_cases/x_only.expected.json`:

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

Write `crates/ppvm-stim/tests/data/edge_cases/repeat_block.expected.json`:

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

Write `crates/ppvm-stim/tests/data/edge_cases/repetition_code_d3_r3.expected.json`:

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [false, false, false, false, false, false, false, false, false]
}
```

- [ ] **Step 3: Write `expected.json` for the two unsupported fixtures**

Write `crates/ppvm-stim/tests/data/unsupported/swap_unsupported.expected.json`:

```json
{
  "mode": "unsupported",
  "awaiting_phase2_instruction": "SWAP"
}
```

Write `crates/ppvm-stim/tests/data/unsupported/mx_unsupported.expected.json`:

```json
{
  "mode": "unsupported",
  "awaiting_phase2_instruction": "MX"
}
```

(Stim cross-check metadata for these will be added by `regen-stim refresh` once the regen tool is built — Tasks 6+. Phase-2 doesn't need it; `awaiting_phase2_instruction` is sufficient for the harness.)

- [ ] **Step 4: Bootstrap distribution-mode JSONs by capturing ppvm output**

We need `ppvm_bit_means` for `bell_pair`, `ghz`, and `depolarize_smoke`. Write a tiny throwaway test, run it, paste the output. Add this temporarily to `crates/ppvm-stim/tests/stim_corpus.rs` at the bottom (we will remove it in Step 6):

```rust
#[test]
#[ignore]
fn bootstrap_distribution_means() {
    type Tab = GeneralizedTableau<ByteFxHashF64<8>, usize>;
    let cases: &[(&str, u64, usize)] = &[
        ("edge_cases/bell_pair.stim", 0, 256),
        ("edge_cases/ghz.stim", 0, 256),
        ("edge_cases/depolarize_smoke.stim", 0, 256),
    ];
    for (rel, seed, num_shots) in cases {
        let path = data_dir().join(rel);
        let src = std::fs::read_to_string(&path).unwrap();
        let prog = parse(&src).unwrap();
        let tprog = normalize::to_tableau(&prog).unwrap();
        let mut next = *seed;
        let shots = sample(&tprog, *num_shots, || {
            let s = next;
            next = next.wrapping_add(1);
            GeneralizedTableau::<ByteFxHashF64<8>, usize>::new_with_seed(N_QUBITS, COEFF_THRESHOLD, s)
        })
        .unwrap();
        let bits = shots[0].len();
        let mut means = Vec::with_capacity(bits);
        for b in 0..bits {
            let sum: f64 = shots
                .iter()
                .map(|s| if s[b].unwrap() { 1.0 } else { 0.0 })
                .sum();
            means.push(sum / (*num_shots as f64));
        }
        eprintln!("{rel} bits={bits} num_shots={num_shots} ppvm_seed={seed} means={means:?}");
    }
    panic!("bootstrap output above");
}
```

Run it (the `--ignored` flag opts in; the `--nocapture` flag prints `eprintln!` to stdout):

```bash
cargo test -p ppvm-stim --test stim_corpus -- --ignored --nocapture bootstrap_distribution_means 2>&1 | tee /tmp/bootstrap.log
```

Expected: panics with the means printed to stderr above the panic. Pull the three `means=[...]` arrays from the output.

- [ ] **Step 5: Write the three distribution-mode JSONs using captured means**

Using the values from `/tmp/bootstrap.log`, write:

`crates/ppvm-stim/tests/data/edge_cases/bell_pair.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 256,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<paste bit 0>, <paste bit 1>]
}
```

`crates/ppvm-stim/tests/data/edge_cases/ghz.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 256,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<paste bit 0>, <paste bit 1>, <paste bit 2>, <paste bit 3>]
}
```

`crates/ppvm-stim/tests/data/edge_cases/depolarize_smoke.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 256,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<paste bit 0>]
}
```

Replace `<paste bit N>` with the exact f64 from the bootstrap output (e.g. `0.4921875`). Don't round; the harness does bit-exact compare.

- [ ] **Step 6: Remove the bootstrap test and run the harness**

Delete the `bootstrap_distribution_means` test from `crates/ppvm-stim/tests/stim_corpus.rs`. Run:

```bash
cargo test -p ppvm-stim --test stim_corpus
```

Expected: `corpus_pairs_every_stim_with_expected_json` passes (every `.stim` has a sibling `.expected.json`). `corpus_obeys_expectations` passes (all 8 fixtures match their expected JSON).

- [ ] **Step 7: Update the corpus README**

Replace `crates/ppvm-stim/tests/data/README.md` with:

```markdown
# Stim corpus

Fixtures consumed by `crates/ppvm-stim/tests/stim_corpus.rs`. Each fixture is
two committed files:

- `<name>.stim` — the source circuit.
- `<name>.expected.json` — declares the test's expected behavior in one of three
  modes: `deterministic`, `distribution`, or `unsupported`. See the spec
  (`docs/superpowers/specs/2026-04-28-ppvm-stim-test-corpus-design.md`) for
  the schema.

The harness walks the directory tree recursively and asserts that every
`.stim` has a sibling `.expected.json` (and vice versa). It then dispatches
on the JSON's `mode` field.

## Categories

| Subdir | Source | Purpose |
|---|---|---|
| `edge_cases/` | hand-written | Empty programs, REPEAT, every tag/Pi-expression form, dense/sparse measurement, comments/whitespace stress. |
| `noise_channels/` | hand-written | Boundary probabilities (p=0.0, p=1.0) and ordering corner cases. |
| `unsupported/` | hand-written + `regen-stim unsupported` | One fixture per phase-1-unsupported instruction. Flips to `distribution` in phase-2. |
| `generated/codes/` | `regen-stim codes` | `stim gen` sweeps over surface/repetition/color codes. |
| `generated/noise_sweeps/` | `regen-stim noise-sweeps` | Per-channel parameter sweeps. |
| `generated/dialect/` | `regen-stim dialect` | ppvm-specific `I[R_X(...)]`, `S[T]`, etc. |
| `generated/random/` | `regen-stim random` | Random-walk programs. |

## Provenance

- Hand-written fixtures (`edge_cases/`, `noise_channels/`, hand-written
  `unsupported/`) are authored by the ppvm team.
- Generated fixtures are produced by the `regen-stim` Python CLI in
  `crates/ppvm-stim/tests/regen-stim/`. See its README for invocation.

## Regenerating

`regen-stim` is a uv-managed Python tool. Cross-check against Stim happens at
regen time, not at test time:

```bash
cd crates/ppvm-stim/tests/regen-stim
uv sync
uv run regen-stim all          # regenerate everything
uv run regen-stim codes        # subcommand-specific
uv run regen-stim refresh ../data/<category>/<name>.stim
```

The committed `ppvm_bit_means` are bit-exact-compared against ppvm's output
at `cargo test` time. Bit drift here means ppvm's behavior changed — that's
the signal we want.
```

- [ ] **Step 8: Verify and commit**

```bash
cargo test -p ppvm-stim --test stim_corpus
```

Expected: both tests pass; 8 fixtures all green.

```bash
git add crates/ppvm-stim/tests/data/
git add crates/ppvm-stim/tests/stim_corpus.rs
git commit -m "test(stim): convert seed fixtures to JSON-driven schema"
```

---

### Task 3: Edge cases corpus — hand-written, ~20 fixtures

**Files:**
- Create: `crates/ppvm-stim/tests/data/edge_cases/<name>.stim` and `.expected.json` for ~14 new fixtures (in addition to the 6 already moved in Task 2).
- Create: `crates/ppvm-stim/tests/data/edge_cases/README.md`.

The spec lists exactly which corner cases must be covered. We hand-write `.stim` sources, derive expected outputs by hand for deterministic cases, and capture ppvm output for distribution cases (using the same temporary harness pattern as Task 2 Step 4).

The fixture list (one per bullet from spec section "edge_cases/ — hand-written corner cases"):

| Filename | Purpose | Mode |
|---|---|---|
| `empty_program.stim` | empty source | deterministic (0 measurements, empty bitstring) |
| `whitespace_only.stim` | only blanks/comments | deterministic |
| `comments_and_blanks.stim` | comments + blanks interspersed | deterministic |
| `pi_expr_forms.stim` | `I[R_X(theta=pi)]`, `0.5*pi`, `2.0*pi`, plain `1.5708` | deterministic |
| `tag_bare.stim` | `S[T] 0; M 0` after `X` | deterministic |
| `tag_named_param.stim` | `I[R_X(theta=1.0*pi)] 0; M 0` | deterministic |
| `tag_multi_named.stim` | `I[U3(theta=1.0*pi, phi=0.0, lambda=0.0)] 0; M 0` | deterministic |
| `tag_multi_positional.stim` | `I_ERROR[correlated_loss](0.0, 0.0, 0.0) 0 1; M 0 1` | deterministic |
| `tag_mixed.stim` | `S[T,debug] 0; M 0` after `X` | deterministic |
| `repeat_single_line.stim` | `REPEAT 3 { X 0 } M 0` | deterministic |
| `repeat_nested_d3.stim` | three-deep nesting | deterministic |
| `repeat_single_instr.stim` | `REPEAT 1 { H 0 } REPEAT 1 { H 0 } M 0` | deterministic |
| `annotation_rec_target.stim` | `M 0; DETECTOR rec[-1]` | deterministic (already covered by `repetition_code_d3_r3` but smaller-scale fixture is useful) |
| `dense_measurement.stim` | `M 0 1 2 ... 31` after assorted prep | distribution |
| `sparse_measurement.stim` | many single-target `M`s interleaved with gates | distribution |
| `mix_all_gate_families.stim` | one of each Phase-1 gate | distribution |

- [ ] **Step 1: Write each `.stim` source**

Examples (write all of them):

`crates/ppvm-stim/tests/data/edge_cases/empty_program.stim`:

```
```

(yes, completely empty.)

`crates/ppvm-stim/tests/data/edge_cases/whitespace_only.stim`:

```


   

```

(blank lines + a line of spaces.)

`crates/ppvm-stim/tests/data/edge_cases/comments_and_blanks.stim`:

```
# Top-of-file comment.

# Another comment.
X 0
# Comment between gates.

M 0
# Trailing comment.
```

`crates/ppvm-stim/tests/data/edge_cases/pi_expr_forms.stim`:

```
# Every pi-expression shape that the parser must accept.
I[R_X(theta=pi)] 0
I[R_X(theta=0.5*pi)] 1
I[R_X(theta=2.0*pi)] 2
I[R_X(theta=1.5708)] 3
M 0 1 2 3
```

`crates/ppvm-stim/tests/data/edge_cases/tag_bare.stim`:

```
X 0
S[T] 0
S_DAG[T] 0
M 0
```

`crates/ppvm-stim/tests/data/edge_cases/tag_named_param.stim`:

```
I[R_X(theta=1.0*pi)] 0
M 0
```

`crates/ppvm-stim/tests/data/edge_cases/tag_multi_named.stim`:

```
I[U3(theta=1.0*pi, phi=0.0, lambda=0.0)] 0
M 0
```

`crates/ppvm-stim/tests/data/edge_cases/tag_multi_positional.stim`:

```
# I_ERROR[correlated_loss] takes 3 positional args (p_x, p_y, p_z).
I_ERROR[correlated_loss](0.0, 0.0, 0.0) 0 1
M 0 1
```

(p=0 means no loss; the measurements are still deterministic 0.)

`crates/ppvm-stim/tests/data/edge_cases/tag_mixed.stim`:

```
X 0
S[T,debug] 0
S_DAG[T,debug] 0
M 0
```

`crates/ppvm-stim/tests/data/edge_cases/repeat_single_line.stim`:

```
REPEAT 3 { X 0 }
M 0
```

(odd repeat count → deterministic 1.)

`crates/ppvm-stim/tests/data/edge_cases/repeat_nested_d3.stim`:

```
REPEAT 2 {
    REPEAT 2 {
        REPEAT 2 {
            X 0
        }
    }
}
M 0
```

(2*2*2 = 8 X applications → deterministic 0.)

`crates/ppvm-stim/tests/data/edge_cases/repeat_single_instr.stim`:

```
REPEAT 1 { H 0 }
REPEAT 1 { H 0 }
M 0
```

(H followed by H = identity → deterministic 0.)

`crates/ppvm-stim/tests/data/edge_cases/annotation_rec_target.stim`:

```
X 0
M 0
DETECTOR rec[-1]
OBSERVABLE_INCLUDE(0) rec[-1]
```

(annotations are no-ops; deterministic 1.)

`crates/ppvm-stim/tests/data/edge_cases/dense_measurement.stim`:

```
# 32 qubits, half flipped, dense measurement.
H 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
X 16 17 18 19 20 21 22 23
M 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
```

(qubits 16–23 are 1; 24–31 are 0; 0–15 are 50/50. Distribution mode.)

`crates/ppvm-stim/tests/data/edge_cases/sparse_measurement.stim`:

```
# Many single-target M's interleaved with gates.
H 0
M 0
X 1
M 1
H 0
M 0
H 1
M 1
S 0
M 0
```

(distribution mode.)

`crates/ppvm-stim/tests/data/edge_cases/mix_all_gate_families.stim`:

```
# At least one of every phase-1 supported gate.
R 0 1 2 3
H 0
X 1
Y 2
Z 3
S 0
S_DAG 1
SQRT_X 2
SQRT_X_DAG 3
SQRT_Y 0
SQRT_Y_DAG 1
S[T] 2
S_DAG[T] 3
I[R_X(theta=0.5*pi)] 0
I[R_Y(theta=0.5*pi)] 1
I[R_Z(theta=0.5*pi)] 2
I[U3(theta=0.0, phi=0.0, lambda=0.0)] 3
CX 0 1
CY 1 2
CZ 2 3
DEPOLARIZE1(0.001) 0
DEPOLARIZE2(0.001) 1 2
X_ERROR(0.001) 0
PAULI_CHANNEL_1(0.001, 0.001, 0.001) 0
M 0 1 2 3
```

(distribution mode.)

- [ ] **Step 2: Write deterministic JSONs by hand**

For each deterministic fixture, work out the bitstring by hand:

`empty_program.expected.json`:

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": []
}
```

`whitespace_only.expected.json`:

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": []
}
```

`comments_and_blanks.expected.json`:

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

`pi_expr_forms.expected.json`: `R_X(pi)` flips |0⟩ → |1⟩ (with a phase that doesn't affect Z-measurement). `R_X(0.5*pi)` produces 50/50 — wait, that means `pi_expr_forms` is **distribution mode**, not deterministic. Replace with deterministic combinations only:

Re-do the source:

```
# All pi-expression shapes; each rotation cancels with its inverse.
I[R_X(theta=pi)] 0
I[R_X(theta=pi)] 0
I[R_X(theta=2.0*pi)] 1
I[R_X(theta=0.5*pi)] 2
I[R_X(theta=1.5*pi)] 2
M 0 1 2
```

Now: bit 0 = 0 (R_X(π) twice = identity up to phase = stays |0⟩ → 0). Bit 1 = 0 (R_X(2π) = identity → stays |0⟩ → 0). Bit 2 = 0 (R_X(π/2) followed by R_X(3π/2) = R_X(2π) = identity → 0).

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [false, false, false]
}
```

`tag_bare.expected.json`: `X 0; S[T] 0; S_DAG[T] 0; M 0` → S[T] is T, S_DAG[T] is T_dag; T then T_dag = identity (up to phase invisible to Z-measurement). So result is just `X 0; M 0` = 1.

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

`tag_named_param.expected.json`: `R_X(π)` flips |0⟩ → |1⟩. Result = 1.

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

`tag_multi_named.expected.json`: `U3(theta=π, phi=0, lambda=0)` is a Y gate up to global phase, sending |0⟩ → e^{iφ}|1⟩. Z-measurement → 1.

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

`tag_multi_positional.expected.json`: `I_ERROR[correlated_loss](0,0,0)` with all probabilities 0 → no-op. M 0 1 with both qubits in |0⟩ → 0,0.

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [false, false]
}
```

`tag_mixed.expected.json`: same logic as `tag_bare` (the `debug` tag is unused metadata; `T` is the dialect tag).

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

`repeat_single_line.expected.json`:

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

`repeat_nested_d3.expected.json`:

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [false]
}
```

`repeat_single_instr.expected.json`:

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [false]
}
```

`annotation_rec_target.expected.json`:

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

- [ ] **Step 3: Distribution-mode JSONs — bootstrap then commit**

Write a bootstrap test (same shape as Task 2 Step 4) to capture means for `dense_measurement`, `sparse_measurement`, `mix_all_gate_families`. Per the spec, edge_cases/distribution defaults to `num_shots=256`.

Add the test to `tests/stim_corpus.rs` temporarily:

```rust
#[test]
#[ignore]
fn bootstrap_edge_cases_distribution() {
    let cases: &[(&str, u64)] = &[
        ("edge_cases/dense_measurement.stim", 0),
        ("edge_cases/sparse_measurement.stim", 0),
        ("edge_cases/mix_all_gate_families.stim", 0),
    ];
    for (rel, seed) in cases {
        let path = data_dir().join(rel);
        let src = std::fs::read_to_string(&path).unwrap();
        let prog = parse(&src).unwrap();
        let tprog = normalize::to_tableau(&prog).unwrap();
        let mut next = *seed;
        let shots = sample(&tprog, 256, || {
            let s = next; next = next.wrapping_add(1);
            GeneralizedTableau::<ByteFxHashF64<8>, usize>::new_with_seed(N_QUBITS, COEFF_THRESHOLD, s)
        }).unwrap();
        let bits = shots[0].len();
        let means: Vec<f64> = (0..bits).map(|b| {
            let sum: f64 = shots.iter().map(|s| if s[b].unwrap() { 1.0 } else { 0.0 }).sum();
            sum / 256.0
        }).collect();
        eprintln!("{rel} num_shots=256 ppvm_seed={seed} ppvm_bit_means={means:?}");
    }
    panic!("bootstrap output above");
}
```

Run:

```bash
cargo test -p ppvm-stim --test stim_corpus -- --ignored --nocapture bootstrap_edge_cases_distribution 2>&1 | tee /tmp/bootstrap_edge.log
```

Paste each printed `ppvm_bit_means` array into the corresponding JSON:

`dense_measurement.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 256,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<bit 0>, ..., <bit 31>]
}
```

`sparse_measurement.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 256,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<bit 0>, ..., <bit 9>]
}
```

`mix_all_gate_families.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 256,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<bit 0>, ..., <bit 3>]
}
```

Delete the `bootstrap_edge_cases_distribution` test from `tests/stim_corpus.rs`.

- [ ] **Step 4: Run the harness**

```bash
cargo test -p ppvm-stim --test stim_corpus
```

Expected: every fixture green; both tests pass.

- [ ] **Step 5: Write `edge_cases/README.md`**

Write `crates/ppvm-stim/tests/data/edge_cases/README.md`:

```markdown
# edge_cases/

Hand-written corpus of corner cases. No regen tooling required.

Fixture inventory:

- Empty / whitespace-only / comment-only programs.
- Every `pi_expr` shape: bare `pi`, `<coeff>*pi`, plain f64.
- Every tag shape: bare ident, single named, multiple named, multiple positional, mixed.
- REPEAT shapes: single-line, multi-line, three-deep nested, single-instruction body.
- Annotation `rec[-k]` targets (parser tolerates and discards for annotations).
- Dense measurement (`M 0 1 ... 31`).
- Sparse measurement (single-target `M`s interleaved with gates).
- One fixture exercising every Phase-1 gate family.

Distribution-mode fixtures use `num_shots=256`. Bootstrapped at fixture-creation time
by running ppvm and recording per-bit means.

Migrating from the previous flat layout, fixtures `bell_pair.stim`, `ghz.stim`,
`x_only.stim`, `repeat_block.stim`, `repetition_code_d3_r3.stim`, and
`depolarize_smoke.stim` were relocated here.
```

- [ ] **Step 6: Commit**

```bash
git add crates/ppvm-stim/tests/data/edge_cases/
git commit -m "test(stim): expand edge_cases corpus with corner-case fixtures"
```

---

### Task 4: Noise channels corpus — hand-written, ~8 fixtures

**Files:**
- Create: `crates/ppvm-stim/tests/data/noise_channels/<name>.stim` and `.expected.json` for ~8 fixtures.
- Create: `crates/ppvm-stim/tests/data/noise_channels/README.md`.

These are boundary-probability and ordering-invariant fixtures the parameter sweeps don't naturally cover. Per spec, `num_shots=4096` (circuits are tiny, so the wall time per fixture is still ~200ms).

| Filename | Purpose | Mode |
|---|---|---|
| `m_zero_noise.stim` | `M(0.0)` ≡ noiseless `M` | deterministic |
| `mr_one_flips.stim` | `MR(1.0)` always flips recorded bit, but reset still works | distribution |
| `xerror_one_then_m.stim` | `X_ERROR(1.0)` followed by `M` — recorded bit forced to 1 | deterministic |
| `m_then_noise_then_m.stim` | M; X_ERROR(p); M — second M reflects state after error | distribution |
| `depolarize1_then_m_noisy.stim` | DEPOLARIZE1 followed by M(p) on same qubit — multiplicative composition | distribution |
| `xerror_zero_no_op.stim` | `X_ERROR(0.0)` is a no-op | deterministic |
| `depolarize2_then_m.stim` | `DEPOLARIZE2(0.5)` between two qubits | distribution |
| `pauli_channel_1_uniform.stim` | `PAULI_CHANNEL_1(p, p, p)` symmetric | distribution |

- [ ] **Step 1: Write `.stim` sources**

`crates/ppvm-stim/tests/data/noise_channels/m_zero_noise.stim`:

```
X 0
M(0.0) 0
```

`crates/ppvm-stim/tests/data/noise_channels/mr_one_flips.stim`:

```
X 0
MR(1.0) 0
M 0
```

(MR(1.0): measure-with-flip-recorded-bit then reset. Recorded bit = NOT(measured value). Then the second M reads the post-reset state.)

`crates/ppvm-stim/tests/data/noise_channels/xerror_one_then_m.stim`:

```
X_ERROR(1.0) 0
M 0
```

`crates/ppvm-stim/tests/data/noise_channels/m_then_noise_then_m.stim`:

```
M 0
X_ERROR(0.3) 0
M 0
```

`crates/ppvm-stim/tests/data/noise_channels/depolarize1_then_m_noisy.stim`:

```
DEPOLARIZE1(0.1) 0
M(0.05) 0
```

`crates/ppvm-stim/tests/data/noise_channels/xerror_zero_no_op.stim`:

```
X 0
X_ERROR(0.0) 0
M 0
```

`crates/ppvm-stim/tests/data/noise_channels/depolarize2_then_m.stim`:

```
DEPOLARIZE2(0.5) 0 1
M 0 1
```

`crates/ppvm-stim/tests/data/noise_channels/pauli_channel_1_uniform.stim`:

```
PAULI_CHANNEL_1(0.1, 0.1, 0.1) 0
M 0
```

- [ ] **Step 2: Deterministic JSONs by hand**

`m_zero_noise.expected.json`: X then M → 1.

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

`xerror_one_then_m.expected.json`: X_ERROR(1.0) is a deterministic X.

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

`xerror_zero_no_op.expected.json`: X then no-op then M → 1.

```json
{
  "mode": "deterministic",
  "ppvm_seed": 0,
  "bitstring": [true]
}
```

- [ ] **Step 3: Bootstrap distribution JSONs at `num_shots=4096`**

Add to `tests/stim_corpus.rs` temporarily:

```rust
#[test]
#[ignore]
fn bootstrap_noise_channels() {
    let cases: &[(&str, u64)] = &[
        ("noise_channels/mr_one_flips.stim", 0),
        ("noise_channels/m_then_noise_then_m.stim", 0),
        ("noise_channels/depolarize1_then_m_noisy.stim", 0),
        ("noise_channels/depolarize2_then_m.stim", 0),
        ("noise_channels/pauli_channel_1_uniform.stim", 0),
    ];
    for (rel, seed) in cases {
        let path = data_dir().join(rel);
        let src = std::fs::read_to_string(&path).unwrap();
        let prog = parse(&src).unwrap();
        let tprog = normalize::to_tableau(&prog).unwrap();
        let mut next = *seed;
        let shots = sample(&tprog, 4096, || {
            let s = next; next = next.wrapping_add(1);
            GeneralizedTableau::<ByteFxHashF64<8>, usize>::new_with_seed(N_QUBITS, COEFF_THRESHOLD, s)
        }).unwrap();
        let bits = shots[0].len();
        let means: Vec<f64> = (0..bits).map(|b| {
            let sum: f64 = shots.iter().map(|s| if s[b].unwrap() { 1.0 } else { 0.0 }).sum();
            sum / 4096.0
        }).collect();
        eprintln!("{rel} num_shots=4096 ppvm_seed={seed} ppvm_bit_means={means:?}");
    }
    panic!("bootstrap output above");
}
```

Run:

```bash
cargo test -p ppvm-stim --test stim_corpus -- --ignored --nocapture bootstrap_noise_channels 2>&1 | tee /tmp/bootstrap_noise.log
```

Paste each printed array into:

`mr_one_flips.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 4096,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<bit 0 (the MR result)>, <bit 1 (the post-reset M)>]
}
```

`m_then_noise_then_m.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 4096,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<bit 0>, <bit 1>]
}
```

`depolarize1_then_m_noisy.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 4096,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<bit 0>]
}
```

`depolarize2_then_m.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 4096,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<bit 0>, <bit 1>]
}
```

`pauli_channel_1_uniform.expected.json`:

```json
{
  "mode": "distribution",
  "num_shots": 4096,
  "ppvm_seed": 0,
  "ppvm_bit_means": [<bit 0>]
}
```

Delete the `bootstrap_noise_channels` test.

- [ ] **Step 4: Verify**

```bash
cargo test -p ppvm-stim --test stim_corpus
```

Expected: green.

- [ ] **Step 5: Write `noise_channels/README.md`**

```markdown
# noise_channels/

Hand-written noise corner cases. The per-channel breadth lives in
`generated/noise_sweeps/`; this category covers boundary probabilities
(p=0.0, p=1.0), ordering invariants, and noise-channel composition that
sweeps don't naturally produce.

Distribution-mode fixtures use `num_shots=4096` because the circuits are
tiny and statistical drift is the exact signal we lock down.

Loss tests are deliberately absent (per spec Non-Goals: Stim has no oracle
for `I_ERROR[loss]`). Loss is exercised by Rust unit tests in
`crates/ppvm-stim/tests/executor.rs`.
```

- [ ] **Step 6: Commit**

```bash
git add crates/ppvm-stim/tests/data/noise_channels/
git commit -m "test(stim): add noise-channels corpus with boundary fixtures"
```

---

### Task 5: Parser snapshot tests

**Files:**
- Create: `crates/ppvm-stim/tests/parser_snapshots.rs`
- (`insta` was already added in Task 1.)

About 15 representative programs, one per AST shape. `insta` writes initial snapshots into `crates/ppvm-stim/tests/snapshots/` on first run; we review and accept them.

- [ ] **Step 1: Write `tests/parser_snapshots.rs`**

```rust
//! Snapshot tests for the Stim parser.
//!
//! Each test parses a small representative program and snapshots the resulting
//! `Program` debug representation. Re-running flags any AST drift cheaply.
//!
//! On first run insta writes pending snapshots into `tests/snapshots/`.
//! Review with `cargo insta review` (install with `cargo install cargo-insta`).

use ppvm_stim::parse;

fn snap(name: &str, src: &str) {
    let prog = parse(src).unwrap_or_else(|e| panic!("parse failed for {name}: {e}"));
    insta::assert_debug_snapshot!(name, prog);
}

#[test]
fn snapshot_bare_gate() {
    snap("bare_gate", "H 0\n");
}

#[test]
fn snapshot_tagged_gate() {
    snap("tagged_gate", "S[T] 0\n");
}

#[test]
fn snapshot_tagged_gate_with_named_params() {
    snap(
        "tagged_gate_with_named_params",
        "I[R_X(theta=0.5*pi)] 0\n",
    );
}

#[test]
fn snapshot_multi_tag() {
    snap("multi_tag", "S[T,debug] 0\n");
}

#[test]
fn snapshot_args_only_noise() {
    snap("args_only_noise", "DEPOLARIZE1(0.5) 0\n");
}

#[test]
fn snapshot_tag_plus_args() {
    snap(
        "tag_plus_args",
        "I_ERROR[correlated_loss](0.1, 0.2, 0.3) 0 1\n",
    );
}

#[test]
fn snapshot_single_target_measurement() {
    snap("single_target_measurement", "M 0\n");
}

#[test]
fn snapshot_multi_target_measurement_with_noise() {
    snap("multi_target_measurement_with_noise", "M(0.001) 0 1 2\n");
}

#[test]
fn snapshot_annotation_with_rec_target() {
    snap("annotation_with_rec_target", "DETECTOR rec[-1]\n");
}

#[test]
fn snapshot_annotation_with_args() {
    snap("annotation_with_args", "OBSERVABLE_INCLUDE(0)\n");
}

#[test]
fn snapshot_empty_repeat_body() {
    snap("empty_repeat_body", "REPEAT 5 { }\n");
}

#[test]
fn snapshot_multi_instruction_repeat() {
    snap(
        "multi_instruction_repeat",
        "REPEAT 3 {\n    H 0\n    CX 0 1\n    M 0 1\n}\n",
    );
}

#[test]
fn snapshot_nested_repeat() {
    snap(
        "nested_repeat",
        "REPEAT 2 {\n    REPEAT 3 {\n        X 0\n    }\n}\n",
    );
}

#[test]
fn snapshot_comment_heavy_program() {
    snap(
        "comment_heavy_program",
        "# Top comment.\n\nH 0  # trailing comment\n# Mid comment.\nM 0\n",
    );
}

#[test]
fn snapshot_whitespace_stress_program() {
    snap(
        "whitespace_stress_program",
        "  \t  H   0  \t \nCX  0   1\n   M  0  1  \n",
    );
}
```

- [ ] **Step 2: Run, then review and accept snapshots**

Run:

```bash
cargo test -p ppvm-stim --test parser_snapshots 2>&1 | tee /tmp/snap.log
```

Expected: 15 tests fail with `insta` reporting pending snapshots.

Install `cargo-insta` if needed and accept:

```bash
cargo install cargo-insta
cargo insta accept --workspace
```

Re-run to verify:

```bash
cargo test -p ppvm-stim --test parser_snapshots
```

Expected: all 15 pass.

- [ ] **Step 3: Commit**

```bash
git add crates/ppvm-stim/tests/parser_snapshots.rs crates/ppvm-stim/tests/snapshots/
git commit -m "test(stim): add parser snapshot tests for representative AST shapes"
```

---

### Task 6: regen-stim CLI scaffolding — shared library + CLI plumbing

**Files:**
- Modify: `crates/ppvm-stim/tests/regen-stim/pyproject.toml` — add ppvm path-source, dev tools, scripts entry point.
- Modify: `crates/ppvm-stim/tests/regen-stim/README.md`.
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/__init__.py`
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/cli.py`
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/core.py`
- Create: `crates/ppvm-stim/tests/regen-stim/test/test_core.py`

This task fills in the plumbing that every subsequent regen subcommand uses: invoke Stim, invoke ppvm, JSON read/write, the seed-search loop, the tolerance check, and a CLI argparse skeleton with subcommand dispatch. No fixtures get generated yet.

- [ ] **Step 1: Update `pyproject.toml`**

Edit `crates/ppvm-stim/tests/regen-stim/pyproject.toml`:

```toml
[project]
name = "regen-stim"
version = "0.1.0"
description = "Generate ppvm-stim test corpus fixtures by cross-checking against quantumlib/Stim."
readme = "README.md"
requires-python = ">=3.12"
dependencies = [
    "stim>=1.15.0",
    "ppvm",
]

[project.scripts]
regen-stim = "regen_stim.cli:main"

[build-system]
requires = ["uv_build>=0.9.21,<0.10.0"]
build-backend = "uv_build"

[tool.uv.sources]
ppvm = { path = "../../../../ppvm-python", editable = true }

[dependency-groups]
dev = [
    "pytest>=8.0",
]

[tool.ruff]
line-length = 100

[tool.ruff.lint]
select = ["E", "F", "B", "UP", "SIM", "I"]
```

Run:

```bash
cd crates/ppvm-stim/tests/regen-stim
uv sync
uv run python -c "import stim; import ppvm; print(stim.__version__)"
```

Expected: prints a version number ≥ 1.15.0; no import errors. (The `ppvm` import requires the wheel to be built first; if it isn't, run `uv run --project ../../../../ppvm-python --group dev maturin develop --uv` from the regen-stim directory or rely on the editable install — see implementer notes.)

- [ ] **Step 2: Write `core.py` — Stim/ppvm invocation, JSON helpers, seed-search loop**

Create `crates/ppvm-stim/tests/regen-stim/src/regen_stim/__init__.py`:

```python
"""Regenerate ppvm-stim test fixtures."""

__version__ = "0.1.0"
```

Create `crates/ppvm-stim/tests/regen-stim/src/regen_stim/core.py`:

```python
"""Shared library for the regen-stim CLI.

Layout: every subcommand emits one or more (stim_source, category, name) triples
to ``write_fixture``, which:

  1. Runs Stim on the source to compute reference per-bit means.
  2. Searches ppvm seeds in [0, max_seed) until ppvm's empirical means at
     ``num_shots`` are within tolerance_sigma * sqrt(p*(1-p)/N) of Stim's.
  3. Writes ``<name>.stim`` and ``<name>.expected.json`` into
     ``<corpus_root>/<category>/``.

Distribution mode is the default. ``write_unsupported_fixture`` and
``write_deterministic_fixture`` handle the two narrower cases.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import stim
import ppvm

DEFAULT_STIM_SHOTS = 10_000
DEFAULT_MAX_SEED = 32
DEFAULT_TOLERANCE_SIGMA = 5.0


@dataclass(frozen=True)
class CorpusPaths:
    """Where fixtures live on disk. Defaults to the corpus root next to this tool."""

    root: Path

    @classmethod
    def default(cls) -> "CorpusPaths":
        # regen-stim/src/regen_stim/core.py → ../../../data
        here = Path(__file__).resolve()
        root = here.parents[3] / "data"
        return cls(root=root)

    def category_dir(self, category: str) -> Path:
        return self.root / category


@dataclass
class StimReference:
    """Per-bit means computed by running Stim."""

    bit_means: list[float]
    num_shots: int
    seed: int
    stim_version: str


def run_stim(source: str, num_shots: int = DEFAULT_STIM_SHOTS, seed: int = 0) -> StimReference:
    """Run Stim and return per-bit empirical means."""
    circuit = stim.Circuit(source)
    sampler = circuit.compile_sampler(seed=seed)
    samples = sampler.sample(shots=num_shots)
    # samples shape: (num_shots, num_measurements). Per-bit mean over shots.
    n_meas = samples.shape[1] if samples.ndim == 2 else 0
    means = [float(samples[:, i].mean()) if n_meas > 0 else 0.0 for i in range(n_meas)]
    return StimReference(
        bit_means=means,
        num_shots=num_shots,
        seed=seed,
        stim_version=stim.__version__,
    )


@dataclass
class PpvmRun:
    """Per-bit means computed by running ppvm."""

    bit_means: list[float]
    num_shots: int
    seed: int


def run_ppvm(source: str, num_shots: int, seed: int) -> PpvmRun:
    """Run ppvm and return per-bit empirical means.

    Counts qubits by inspecting Stim's parsed circuit (Stim has the parser, ppvm
    expects an n_qubits arg) and uses the same value when constructing tableaux.
    """
    n_qubits = max(64, _max_qubit_in_source(source) + 1)
    prog = ppvm.StimProgram.parse(source)
    shots = ppvm.sample_stim(prog, n_qubits=n_qubits, num_shots=num_shots, seed=seed)
    if not shots:
        return PpvmRun(bit_means=[], num_shots=num_shots, seed=seed)
    n_meas = len(shots[0])
    means: list[float] = []
    for i in range(n_meas):
        s = 0
        for shot in shots:
            v = shot[i]
            if v is None:
                raise ValueError(
                    f"ppvm returned None (loss) for bit {i}; corpus excludes loss "
                    "(see spec Non-Goals)"
                )
            s += 1 if v else 0
        means.append(s / num_shots)
    return PpvmRun(bit_means=means, num_shots=num_shots, seed=seed)


def _max_qubit_in_source(source: str) -> int:
    """Walk Stim's parsed circuit to find the highest qubit index referenced."""
    circuit = stim.Circuit(source)
    max_q = -1
    for inst in circuit.flattened():
        for t in inst.targets_copy():
            if t.is_qubit_target:
                max_q = max(max_q, t.qubit_value)
    return max_q


def per_bit_sigma(stim_means: list[float], num_shots: int) -> list[float]:
    """Worst-case binomial sigma per bit at the test-time num_shots."""
    out: list[float] = []
    for p in stim_means:
        if 0.0 < p < 1.0:
            out.append(math.sqrt(p * (1.0 - p) / num_shots))
        else:
            # p ∈ {0, 1}: empirical mean is always p in the noiseless limit, but
            # any noise term we can't see in stim_means could raise the variance
            # by one order. Use the worst-case bound 1/N as a guardrail.
            out.append(math.sqrt(1.0 / num_shots))
    return out


def within_tolerance(
    ppvm_means: list[float],
    stim_means: list[float],
    test_num_shots: int,
    tolerance_sigma: float = DEFAULT_TOLERANCE_SIGMA,
) -> bool:
    if len(ppvm_means) != len(stim_means):
        return False
    sigmas = per_bit_sigma(stim_means, test_num_shots)
    for p, s, sig in zip(ppvm_means, stim_means, sigmas):
        if abs(p - s) > tolerance_sigma * sig:
            return False
    return True


@dataclass
class FixtureMeta:
    name: str
    category: str
    source: str
    test_num_shots: int
    stim_num_shots: int = DEFAULT_STIM_SHOTS
    stim_seed: int = 0
    max_ppvm_seed: int = DEFAULT_MAX_SEED
    tolerance_sigma: float = DEFAULT_TOLERANCE_SIGMA
    extra_metadata: dict = field(default_factory=dict)


def write_distribution_fixture(meta: FixtureMeta, paths: CorpusPaths) -> Path:
    """Generate a distribution-mode fixture; raises on irreconcilable cross-check."""
    ref = run_stim(meta.source, num_shots=meta.stim_num_shots, seed=meta.stim_seed)
    chosen_seed: int | None = None
    chosen_means: list[float] = []
    for seed in range(meta.max_ppvm_seed):
        ppvm_run = run_ppvm(meta.source, num_shots=meta.test_num_shots, seed=seed)
        if within_tolerance(
            ppvm_run.bit_means, ref.bit_means, meta.test_num_shots, meta.tolerance_sigma
        ):
            chosen_seed = seed
            chosen_means = ppvm_run.bit_means
            break
    if chosen_seed is None:
        raise RuntimeError(
            f"{meta.category}/{meta.name}: no ppvm seed in [0, {meta.max_ppvm_seed}) "
            f"agrees with Stim within {meta.tolerance_sigma} sigma at "
            f"num_shots={meta.test_num_shots}. Stim means: {ref.bit_means}. "
            "This is a real correctness divergence — do not commit."
        )

    payload = {
        "mode": "distribution",
        "num_shots": meta.test_num_shots,
        "ppvm_seed": chosen_seed,
        "ppvm_bit_means": chosen_means,
        "stim_seed": meta.stim_seed,
        "stim_num_shots": meta.stim_num_shots,
        "stim_bit_means": ref.bit_means,
        "tolerance_sigma_at_regen": meta.tolerance_sigma,
        "stim_version": ref.stim_version,
        **meta.extra_metadata,
    }
    return _emit(meta, paths, payload)


def write_deterministic_fixture(meta: FixtureMeta, paths: CorpusPaths) -> Path:
    """One-shot ppvm run; record the bitstring."""
    ppvm_run = run_ppvm(meta.source, num_shots=1, seed=0)
    if not ppvm_run.bit_means:
        bitstring: list[bool] = []
    else:
        # bit_means at num_shots=1 are 0.0 or 1.0; ppvm.sample_stim returns
        # actual bools per shot. Re-run via the underlying API to get bools.
        n_qubits = max(64, _max_qubit_in_source(meta.source) + 1)
        prog = ppvm.StimProgram.parse(meta.source)
        shots = ppvm.sample_stim(prog, n_qubits=n_qubits, num_shots=1, seed=0)
        bitstring = [bool(b) for b in shots[0]]

    payload = {
        "mode": "deterministic",
        "ppvm_seed": 0,
        "bitstring": bitstring,
        **meta.extra_metadata,
    }
    return _emit(meta, paths, payload)


def write_unsupported_fixture(
    meta: FixtureMeta,
    paths: CorpusPaths,
    awaiting_phase2_instruction: str,
) -> Path:
    """Pre-record Stim reference for the day phase-2 lifts the restriction."""
    ref = run_stim(meta.source, num_shots=meta.stim_num_shots, seed=meta.stim_seed)
    payload = {
        "mode": "unsupported",
        "awaiting_phase2_instruction": awaiting_phase2_instruction,
        "stim_seed": meta.stim_seed,
        "stim_num_shots": meta.stim_num_shots,
        "stim_bit_means": ref.bit_means,
        "tolerance_sigma_at_regen": meta.tolerance_sigma,
        "stim_version": ref.stim_version,
        **meta.extra_metadata,
    }
    return _emit(meta, paths, payload)


def _emit(meta: FixtureMeta, paths: CorpusPaths, payload: dict) -> Path:
    cat_dir = paths.category_dir(meta.category)
    cat_dir.mkdir(parents=True, exist_ok=True)
    stim_path = cat_dir / f"{meta.name}.stim"
    json_path = cat_dir / f"{meta.name}.expected.json"
    stim_path.write_text(meta.source)
    json_path.write_text(json.dumps(payload, indent=2) + "\n")
    return json_path
```

- [ ] **Step 3: Write `cli.py`**

Create `crates/ppvm-stim/tests/regen-stim/src/regen_stim/cli.py`:

```python
"""regen-stim entry point."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="regen-stim", description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("codes", help="(Task 7) Generate generated/codes/ via stim gen sweeps")
    sub.add_parser("noise-sweeps", help="(Task 8) Generate generated/noise_sweeps/ per-channel")
    sub.add_parser("dialect", help="(Task 9) Generate generated/dialect/ ppvm-specific")
    sub.add_parser("random", help="(Task 10) Generate generated/random/ random-walk")
    sub.add_parser("unsupported", help="(Task 11) Generate unsupported/ phase-1-rejected")
    p_refresh = sub.add_parser("refresh", help="Re-record one fixture by path")
    p_refresh.add_argument("path", type=Path)
    p_verify = sub.add_parser("verify", help="Re-run cross-check; error on mismatch; write nothing")
    p_verify.add_argument("path", type=Path)
    sub.add_parser("all", help="Run every regen subcommand")

    args = parser.parse_args(argv)

    if args.cmd == "codes":
        from . import codes
        return codes.run()
    if args.cmd == "noise-sweeps":
        from . import noise_sweeps
        return noise_sweeps.run()
    if args.cmd == "dialect":
        from . import dialect
        return dialect.run()
    if args.cmd == "random":
        from . import random_walk
        return random_walk.run()
    if args.cmd == "unsupported":
        from . import unsupported
        return unsupported.run()
    if args.cmd == "refresh":
        return _refresh(args.path)
    if args.cmd == "verify":
        return _verify(args.path)
    if args.cmd == "all":
        rc = 0
        for mod_name in ("codes", "noise_sweeps", "dialect", "random_walk", "unsupported"):
            mod = __import__(f"regen_stim.{mod_name}", fromlist=["run"])
            rc |= mod.run()
        return rc

    parser.print_help()
    return 2


def _refresh(path: Path) -> int:
    """Re-emit one fixture given a path to <name>.stim or <name>.expected.json."""
    from . import core

    stim_path = path if path.suffix == ".stim" else path.with_suffix("").with_suffix(".stim")
    if stim_path.suffix == ".expected" and stim_path.stem.endswith(".expected"):
        stim_path = stim_path.with_name(stim_path.stem.removesuffix(".expected") + ".stim")
    if not stim_path.exists():
        print(f"refresh: no .stim file at {stim_path}", file=sys.stderr)
        return 1

    # Read existing JSON for metadata; re-derive payload.
    json_path = stim_path.with_name(stim_path.stem + ".expected.json")
    if not json_path.exists():
        print(f"refresh: no expected.json at {json_path}", file=sys.stderr)
        return 1

    import json
    existing = json.loads(json_path.read_text())
    src = stim_path.read_text()
    paths = core.CorpusPaths.default()
    category = stim_path.parent.relative_to(paths.root).as_posix()
    name = stim_path.stem

    mode = existing.get("mode")
    if mode == "distribution":
        meta = core.FixtureMeta(
            name=name,
            category=category,
            source=src,
            test_num_shots=existing["num_shots"],
            stim_num_shots=existing.get("stim_num_shots", core.DEFAULT_STIM_SHOTS),
            stim_seed=existing.get("stim_seed", 0),
            tolerance_sigma=existing.get("tolerance_sigma_at_regen", core.DEFAULT_TOLERANCE_SIGMA),
        )
        core.write_distribution_fixture(meta, paths)
    elif mode == "deterministic":
        meta = core.FixtureMeta(name=name, category=category, source=src, test_num_shots=1)
        core.write_deterministic_fixture(meta, paths)
    elif mode == "unsupported":
        meta = core.FixtureMeta(
            name=name,
            category=category,
            source=src,
            test_num_shots=0,
            stim_num_shots=existing.get("stim_num_shots", core.DEFAULT_STIM_SHOTS),
            stim_seed=existing.get("stim_seed", 0),
        )
        core.write_unsupported_fixture(
            meta, paths, awaiting_phase2_instruction=existing["awaiting_phase2_instruction"]
        )
    else:
        print(f"refresh: unknown mode {mode!r}", file=sys.stderr)
        return 1
    print(f"refreshed {category}/{name}")
    return 0


def _verify(path: Path) -> int:
    """Re-run cross-check without writing; error on mismatch."""
    from . import core

    stim_path = path if path.suffix == ".stim" else path.with_suffix("").with_suffix(".stim")
    json_path = stim_path.with_name(stim_path.stem + ".expected.json")
    import json
    existing = json.loads(json_path.read_text())
    src = stim_path.read_text()
    if existing["mode"] != "distribution":
        print("verify: only meaningful for distribution mode", file=sys.stderr)
        return 0  # nothing to do for det/unsup
    ref = core.run_stim(
        src,
        num_shots=existing["stim_num_shots"],
        seed=existing.get("stim_seed", 0),
    )
    ppvm_run = core.run_ppvm(
        src,
        num_shots=existing["num_shots"],
        seed=existing["ppvm_seed"],
    )
    ok = core.within_tolerance(
        ppvm_run.bit_means,
        ref.bit_means,
        existing["num_shots"],
        existing.get("tolerance_sigma_at_regen", core.DEFAULT_TOLERANCE_SIGMA),
    )
    if not ok:
        print(f"verify: mismatch at {stim_path}", file=sys.stderr)
        print(f"  stim_means: {ref.bit_means}", file=sys.stderr)
        print(f"  ppvm_means: {ppvm_run.bit_means}", file=sys.stderr)
        return 1
    print(f"verify: {stim_path} OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

The CLI references modules `codes`, `noise_sweeps`, `dialect`, `random_walk`, `unsupported` that don't exist yet — those land in Tasks 7–11. The CLI imports them lazily so unused subcommands don't break the build.

- [ ] **Step 4: Write `test/test_core.py`**

Create `crates/ppvm-stim/tests/regen-stim/test/test_core.py`:

```python
"""Unit tests for the regen-stim core helpers.

Tests the math (per_bit_sigma, within_tolerance) and a smoke test for
run_stim. ppvm-side smoke (run_ppvm, write_distribution_fixture) lives
in the integration tests but exercising it requires the ppvm wheel built.
"""

import math

import pytest

from regen_stim import core


def test_per_bit_sigma_interior_probabilities():
    sigmas = core.per_bit_sigma([0.5, 0.1, 0.9], num_shots=100)
    assert math.isclose(sigmas[0], math.sqrt(0.5 * 0.5 / 100))
    assert math.isclose(sigmas[1], math.sqrt(0.1 * 0.9 / 100))
    assert math.isclose(sigmas[2], math.sqrt(0.9 * 0.1 / 100))


def test_per_bit_sigma_boundary_probabilities():
    sigmas = core.per_bit_sigma([0.0, 1.0], num_shots=100)
    # Worst-case 1/N guardrail at the boundary.
    assert math.isclose(sigmas[0], math.sqrt(1.0 / 100))
    assert math.isclose(sigmas[1], math.sqrt(1.0 / 100))


def test_within_tolerance_tight_match():
    assert core.within_tolerance([0.5], [0.5], test_num_shots=1024, tolerance_sigma=5.0)


def test_within_tolerance_drift():
    # 50% off at p=0.5 with N=1024: sigma = sqrt(0.25/1024) ≈ 0.0156.
    # 5*sigma ≈ 0.078. Drift 0.10 is outside tolerance.
    assert not core.within_tolerance(
        [0.6], [0.5], test_num_shots=1024, tolerance_sigma=5.0
    )


def test_within_tolerance_length_mismatch_is_failure():
    assert not core.within_tolerance(
        [0.5, 0.5], [0.5], test_num_shots=1024, tolerance_sigma=5.0
    )


def test_run_stim_smoke_x_then_m():
    ref = core.run_stim("X 0\nM 0\n", num_shots=64)
    assert ref.bit_means == [1.0]
    assert ref.num_shots == 64
    assert ref.stim_version  # non-empty


def test_run_stim_smoke_h_then_m_is_random():
    ref = core.run_stim("H 0\nM 0\n", num_shots=10_000, seed=42)
    # 10k shots at p=0.5 → sigma = 0.005; 0.4 < mean < 0.6 trivially.
    assert 0.4 < ref.bit_means[0] < 0.6


def test_max_qubit_in_source_indexes_correctly():
    assert core._max_qubit_in_source("M 0\n") == 0
    assert core._max_qubit_in_source("CX 0 7\nM 0 7\n") == 7
```

- [ ] **Step 5: Run the regen-stim tests**

```bash
cd crates/ppvm-stim/tests/regen-stim
uv sync
uv run pytest test/
```

Expected: 8 tests pass.

Note: this depends on `stim` being installed. If `uv sync` complains that `ppvm` can't be resolved (because the wheel isn't built), that's fine — the unit tests don't import `ppvm`. If pytest imports it transitively via `core.py` and fails, edit `core.py` to defer the `import ppvm` to the function bodies that need it. The cleanest fix is `import ppvm` inside `run_ppvm` only.

- [ ] **Step 6: Smoke-run the CLI**

```bash
uv run regen-stim --help
```

Expected: prints subcommand list. `regen-stim codes` will fail importing `codes` (not yet implemented) — that's expected and gets fixed in Task 7.

- [ ] **Step 7: Update README**

Edit `crates/ppvm-stim/tests/regen-stim/README.md`:

```markdown
# regen-stim

Generates ppvm-stim test corpus fixtures by cross-checking against
[`quantumlib/Stim`](https://github.com/quantumlib/Stim). Run by humans, not in CI.

## Setup

```bash
cd crates/ppvm-stim/tests/regen-stim
uv sync
```

`uv sync` resolves the editable `ppvm` Python package. The first sync requires
the `ppvm-python-native` extension built — run from the repo root once:

```bash
uv run --project ppvm-python --group dev maturin develop --uv
```

## Commands

```bash
uv run regen-stim codes        # generated/codes/  (stim gen sweeps)
uv run regen-stim noise-sweeps # generated/noise_sweeps/
uv run regen-stim dialect      # generated/dialect/
uv run regen-stim random       # generated/random/
uv run regen-stim unsupported  # unsupported/
uv run regen-stim refresh ../data/<category>/<name>.stim
uv run regen-stim verify  ../data/<category>/<name>.stim
uv run regen-stim all          # everything
```

## Dev

```bash
uv run pytest test/   # unit tests for the seed-search loop and tolerance math
```

## How it works

For distribution-mode fixtures, the seed-search loop in
`src/regen_stim/core.py` runs Stim at a high shot count to compute reference
per-bit means, then tries ppvm seeds in [0, 32) until ppvm's empirical means
(at the test-time `num_shots`) are within `tolerance_sigma * sigma` of Stim's.

Test-time runs use the committed seed at the committed shot count and
bit-exact-compare against the committed `ppvm_bit_means`. The cross-check
against Stim only happens here, at regen time.
```

- [ ] **Step 8: Commit**

```bash
git add crates/ppvm-stim/tests/regen-stim/
git commit -m "tools(stim): scaffold regen-stim CLI and core library"
```

---

### Task 7: regen-stim codes — generate `generated/codes/` corpus

**Files:**
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/codes.py`
- Create: `crates/ppvm-stim/tests/data/generated/codes/<many>.stim` and `<many>.expected.json`
- Create: `crates/ppvm-stim/tests/data/generated/codes/README.md`

The `codes` subcommand sweeps `stim gen` over `(code, task, distance, rounds, noise)` per the spec. Some `(code, task)` combinations emit phase-1-unsupported instructions (e.g. `unrotated_memory_z` uses `MR`, transversal tasks use `SWAP`, stabilizer-measurement tasks use `MPP`). The script auto-routes those into `unsupported/` instead.

- [ ] **Step 1: Write `codes.py`**

Create `crates/ppvm-stim/tests/regen-stim/src/regen_stim/codes.py`:

```python
"""generated/codes/: stim gen sweeps over surface, repetition, and color codes."""

from __future__ import annotations

import re

import stim

from . import core

CODES_TASKS: dict[str, list[str]] = {
    "surface_code": [
        "unrotated_memory_x",
        "unrotated_memory_z",
        "rotated_memory_x",
        "rotated_memory_z",
    ],
    "repetition_code": ["memory"],
    "color_code": ["memory_xyz"],
}

DISTANCES = [3, 5, 7]
ROUNDS = [1, 3, 5]
NOISE_VALUES: list[float | None] = [None, 0.001, 0.01]

# Phase-1 supported instruction names (parser + normalizer accepts each).
PHASE1_SUPPORTED = {
    # Resets / single-qubit Cliffords / two-qubit Cliffords.
    "R", "RZ", "I", "X", "Y", "Z", "H", "H_XZ",
    "S", "S_DAG", "SQRT_Z", "SQRT_Z_DAG", "SQRT_X", "SQRT_X_DAG",
    "SQRT_Y", "SQRT_Y_DAG",
    "CX", "ZCX", "CNOT", "CY", "ZCY", "CZ", "ZCZ",
    # Measurements.
    "M", "MZ", "MR",
    # Noise.
    "DEPOLARIZE1", "DEPOLARIZE2", "PAULI_CHANNEL_1", "PAULI_CHANNEL_2",
    "X_ERROR", "Y_ERROR", "Z_ERROR", "I_ERROR",
    # Annotations.
    "DETECTOR", "OBSERVABLE_INCLUDE", "TICK", "QUBIT_COORDS", "SHIFT_COORDS",
    "REPEAT",
}

INSTR_RE = re.compile(r"^([A-Z_][A-Z0-9_]*)(?:\[|\(|\s|$)", re.MULTILINE)


def first_unsupported_instruction(source: str) -> str | None:
    """Return the first phase-1-unsupported instruction name, or None."""
    for m in INSTR_RE.finditer(source):
        name = m.group(1)
        if name == "REPEAT":
            continue
        if name not in PHASE1_SUPPORTED:
            return name
    return None


def gen_circuit(code: str, task: str, distance: int, rounds: int, noise: float | None) -> str:
    kwargs: dict = dict(code=code, task=task, distance=distance, rounds=rounds)
    if noise is not None:
        kwargs.update(
            after_clifford_depolarization=noise,
            before_round_data_depolarization=noise,
            before_measure_flip_probability=noise,
            after_reset_flip_probability=noise,
        )
    circuit = stim.Circuit.generated(**kwargs)
    return str(circuit) + "\n"


def fixture_name(code: str, task: str, distance: int, rounds: int, noise: float | None) -> str:
    noise_str = "noiseless" if noise is None else f"p{noise:g}".replace("0.", "")
    return f"{code}_{task}_d{distance}_r{rounds}_{noise_str}"


def shot_count_for(distance: int) -> int:
    # Per spec: distance ≤ 5 → 256 shots; distance = 7 → 64 shots.
    return 64 if distance >= 7 else 256


def run() -> int:
    paths = core.CorpusPaths.default()
    failures: list[str] = []
    written: list[str] = []
    for code, tasks in CODES_TASKS.items():
        for task in tasks:
            for distance in DISTANCES:
                for rounds in ROUNDS:
                    for noise in NOISE_VALUES:
                        try:
                            src = gen_circuit(code, task, distance, rounds, noise)
                        except Exception as e:
                            failures.append(
                                f"stim gen failed for {code}/{task}/d{distance}/r{rounds}: {e}"
                            )
                            continue
                        name = fixture_name(code, task, distance, rounds, noise)
                        unsupported = first_unsupported_instruction(src)
                        if unsupported is not None:
                            meta = core.FixtureMeta(
                                name=name,
                                category="unsupported",
                                source=src,
                                test_num_shots=0,
                            )
                            try:
                                core.write_unsupported_fixture(
                                    meta, paths, awaiting_phase2_instruction=unsupported
                                )
                                written.append(f"unsupported/{name}")
                            except Exception as e:
                                failures.append(f"unsupported/{name}: {e}")
                            continue
                        # Supported → distribution mode.
                        test_shots = shot_count_for(distance)
                        meta = core.FixtureMeta(
                            name=name,
                            category="generated/codes",
                            source=src,
                            test_num_shots=test_shots,
                        )
                        try:
                            core.write_distribution_fixture(meta, paths)
                            written.append(f"generated/codes/{name}")
                        except Exception as e:
                            failures.append(f"generated/codes/{name}: {e}")
    print(f"regen-stim codes: wrote {len(written)} fixtures")
    if failures:
        print("regen-stim codes: failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0
```

- [ ] **Step 2: Run the regen subcommand**

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim codes
```

Expected: prints "wrote N fixtures" with N in the 80–100 range; some fixtures land in `unsupported/` (auto-routed by `first_unsupported_instruction`). May take 5–30 minutes depending on shot count and the seed-search loop.

If the script errors with "no ppvm seed in [0, 32) agrees with Stim within 5.0 sigma", that's a real correctness divergence between ppvm and Stim — investigate before proceeding. Common causes:
- ppvm noise-channel parameterization differs from Stim's (e.g. `DEPOLARIZE2` arg ordering).
- ppvm's reset noise (`after_reset_flip_probability`) interpretation differs.
- `MR`'s recorded-bit semantics differ.

If a single fixture fails but others pass, examine `stim.Circuit.generated` output for that fixture by hand and compare against ppvm's interpretation. If multiple fixtures fail with the same channel, fix the underlying ppvm bug.

- [ ] **Step 3: Run the cargo harness**

```bash
cd ../../../..  # back to repo root
cargo test -p ppvm-stim --test stim_corpus
```

Expected: every fixture green. Per-fixture wall time at 256 shots: ~125ms; total runtime for the codes corpus: ~10–15 seconds.

- [ ] **Step 4: Write `generated/codes/README.md`**

```markdown
# generated/codes/

`stim gen` sweep over fault-tolerance circuits. Generated by:

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim codes
```

Sweep axes:
- `--code` ∈ {`surface_code`, `repetition_code`, `color_code`}
- `--task` ∈ all task variants (some emit phase-1-unsupported instructions and land in `../unsupported/` instead).
- `--distance` ∈ {3, 5, 7}
- `--rounds` ∈ {1, 3, 5}
- noise: noiseless, p=0.001, p=0.01

Distance ≤ 5 fixtures use `num_shots=256`; distance = 7 uses `num_shots=64`.

The `regen-stim codes` script auto-routes fixtures whose `stim gen` output
contains a phase-1-unsupported instruction (e.g. `MR`, `SWAP`, `MPP`) into
`../unsupported/`.
```

- [ ] **Step 5: Commit**

```bash
git add crates/ppvm-stim/tests/regen-stim/src/regen_stim/codes.py
git add crates/ppvm-stim/tests/data/generated/codes/
git add crates/ppvm-stim/tests/data/unsupported/  # if codes auto-routed any here
git commit -m "test(stim): generate codes/ corpus from stim gen sweeps"
```

---

### Task 8: regen-stim noise-sweeps — generate `generated/noise_sweeps/` corpus

**Files:**
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/noise_sweeps.py`
- Create: `crates/ppvm-stim/tests/data/generated/noise_sweeps/<many>.stim` + `<many>.expected.json`
- Create: `crates/ppvm-stim/tests/data/generated/noise_sweeps/README.md`

Per-channel parameter sweeps with circuits of uniform shape: prepare-state → apply-channel-N-times → measure-all. `num_shots=4096` per spec.

- [ ] **Step 1: Write `noise_sweeps.py`**

Create `crates/ppvm-stim/tests/regen-stim/src/regen_stim/noise_sweeps.py`:

```python
"""generated/noise_sweeps/: per-channel parameter sweeps."""

from __future__ import annotations

from . import core

NUM_SHOTS = 4096

# (channel_name, args, n_qubits_arg, qubit_counts, prep)
# args is a callable taking p and producing the parens body.

def depolarize1_args(p: float) -> str: return f"({p})"
def depolarize2_args(p: float) -> str: return f"({p})"
def xerror_args(p: float) -> str: return f"({p})"
def yerror_args(p: float) -> str: return f"({p})"
def zerror_args(p: float) -> str: return f"({p})"

def pauli_channel_1_args(params: tuple[float, float, float]) -> str:
    return f"({params[0]}, {params[1]}, {params[2]})"

def pauli_channel_2_args(params: tuple[float, ...]) -> str:
    # 15 parameters: II is implicit; IX, IY, IZ, XI, ..., ZZ.
    assert len(params) == 15
    return "(" + ", ".join(str(p) for p in params) + ")"


SINGLE_QUBIT_PROBS = [0.001, 0.01, 0.1, 0.5]
TWO_QUBIT_PROBS = [0.001, 0.01, 0.1, 0.5]
SMALL_BIG_PROBS = [0.001, 0.01, 0.5]

PAULI1_PARAM_SETS: list[tuple[float, float, float]] = [
    (0.05, 0.0, 0.0),    # X-only
    (0.0, 0.0, 0.1),     # Z-only
    (0.05, 0.05, 0.05),  # symmetric
]

# 15-arg Pauli channel: IX IY IZ XI XX XY XZ YI YX YY YZ ZI ZX ZY ZZ
PAULI2_PARAM_SETS: list[tuple[float, ...]] = [
    (0.01,) + (0.0,) * 14,                              # IX-only
    (0.0,) * 7 + (0.01,) + (0.0,) * 7,                   # YI
    tuple([1.0 / 16.0] * 15),                            # symmetric (uniform pairwise)
]

READOUT_PROBS = [0.001, 0.01, 0.5]


def fixture_source_per_qubit_channel(
    channel: str, args: str, n_qubits: int, repeat: int = 5
) -> str:
    """Build: H 0 1 ... n-1; (channel <args> 0 1 ... n-1) × repeat; M 0 1 ... n-1."""
    qs = " ".join(str(i) for i in range(n_qubits))
    body = "\n".join(f"{channel}{args} {qs}" for _ in range(repeat))
    return f"H {qs}\n{body}\nM {qs}\n"


def fixture_source_pair_channel(
    channel: str, args: str, n_qubits: int, repeat: int = 5
) -> str:
    """Two-qubit channel; targets pair up."""
    assert n_qubits % 2 == 0
    pairs = " ".join(str(i) for i in range(n_qubits))
    body = "\n".join(f"{channel}{args} {pairs}" for _ in range(repeat))
    qs = " ".join(str(i) for i in range(n_qubits))
    return f"H {qs}\n{body}\nM {qs}\n"


def fixture_source_readout_noise(measure: str, p: float, n_qubits: int) -> str:
    """Apply readout noise via M(p) or MR(p) on a prepared state."""
    qs = " ".join(str(i) for i in range(n_qubits))
    return f"X {qs}\n{measure}({p}) {qs}\n"


def run() -> int:
    paths = core.CorpusPaths.default()
    failures: list[str] = []
    written = 0

    # DEPOLARIZE1: 4 probs × 2 qubit counts = 8.
    for n in (1, 4):
        for p in SINGLE_QUBIT_PROBS:
            src = fixture_source_per_qubit_channel("DEPOLARIZE1", f"({p})", n)
            name = f"depolarize1_n{n}_p{p:g}".replace("0.", "")
            meta = core.FixtureMeta(
                name=name, category="generated/noise_sweeps", source=src,
                test_num_shots=NUM_SHOTS,
            )
            try:
                core.write_distribution_fixture(meta, paths); written += 1
            except Exception as e:
                failures.append(f"{name}: {e}")

    # DEPOLARIZE2: 4 probs × 2 even qubit counts = 8.
    for n in (2, 4):
        for p in TWO_QUBIT_PROBS:
            src = fixture_source_pair_channel("DEPOLARIZE2", f"({p})", n)
            name = f"depolarize2_n{n}_p{p:g}".replace("0.", "")
            meta = core.FixtureMeta(
                name=name, category="generated/noise_sweeps", source=src,
                test_num_shots=NUM_SHOTS,
            )
            try:
                core.write_distribution_fixture(meta, paths); written += 1
            except Exception as e:
                failures.append(f"{name}: {e}")

    # PAULI_CHANNEL_1: 3 sets × 2 qubit counts = 6.
    for n in (1, 4):
        for i, params in enumerate(PAULI1_PARAM_SETS):
            src = fixture_source_per_qubit_channel(
                "PAULI_CHANNEL_1", pauli_channel_1_args(params), n
            )
            name = f"pauli_channel_1_set{i}_n{n}"
            meta = core.FixtureMeta(
                name=name, category="generated/noise_sweeps", source=src,
                test_num_shots=NUM_SHOTS,
            )
            try:
                core.write_distribution_fixture(meta, paths); written += 1
            except Exception as e:
                failures.append(f"{name}: {e}")

    # PAULI_CHANNEL_2: 3 sets × 2 even qubit counts = 6.
    for n in (2, 4):
        for i, params in enumerate(PAULI2_PARAM_SETS):
            src = fixture_source_pair_channel(
                "PAULI_CHANNEL_2", pauli_channel_2_args(params), n
            )
            name = f"pauli_channel_2_set{i}_n{n}"
            meta = core.FixtureMeta(
                name=name, category="generated/noise_sweeps", source=src,
                test_num_shots=NUM_SHOTS,
            )
            try:
                core.write_distribution_fixture(meta, paths); written += 1
            except Exception as e:
                failures.append(f"{name}: {e}")

    # X_ERROR / Y_ERROR / Z_ERROR: 3 channels × 3 probs × 1 qubit count = 9.
    for ch in ("X_ERROR", "Y_ERROR", "Z_ERROR"):
        for p in SMALL_BIG_PROBS:
            src = fixture_source_per_qubit_channel(ch, f"({p})", 1)
            name = f"{ch.lower()}_p{p:g}".replace("0.", "")
            meta = core.FixtureMeta(
                name=name, category="generated/noise_sweeps", source=src,
                test_num_shots=NUM_SHOTS,
            )
            try:
                core.write_distribution_fixture(meta, paths); written += 1
            except Exception as e:
                failures.append(f"{name}: {e}")

    # M(p) / MR(p): 2 measurements × 3 probs = 6.
    for measure in ("M", "MR"):
        for p in READOUT_PROBS:
            src = fixture_source_readout_noise(measure, p, n_qubits=1)
            name = f"{measure.lower()}_readout_p{p:g}".replace("0.", "")
            meta = core.FixtureMeta(
                name=name, category="generated/noise_sweeps", source=src,
                test_num_shots=NUM_SHOTS,
            )
            try:
                core.write_distribution_fixture(meta, paths); written += 1
            except Exception as e:
                failures.append(f"{name}: {e}")

    print(f"regen-stim noise-sweeps: wrote {written} fixtures")
    if failures:
        print("regen-stim noise-sweeps: failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0
```

- [ ] **Step 2: Run the subcommand**

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim noise-sweeps
```

Expected: ~43 fixtures generated. Wall time: ~5–10 minutes (4096 shots is slow per fixture but each circuit is tiny).

- [ ] **Step 3: Verify the harness passes**

```bash
cd ../../../..
cargo test -p ppvm-stim --test stim_corpus
```

Expected: green.

- [ ] **Step 4: Write `generated/noise_sweeps/README.md`**

```markdown
# generated/noise_sweeps/

Per-channel parameter sweeps. Generated by:

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim noise-sweeps
```

Each fixture has uniform shape: H-prepare → apply-channel-N-times → measure-all.
The probability axis sweep gives the regen-time Stim cross-check enough breadth
to catch slope/sign bugs on each channel.

Fixtures are distribution-mode at `num_shots=4096` (statistical signal needs
many shots; circuits are tiny so the wall time is acceptable).

Loss channels (`I_ERROR[loss]`) are absent — Stim has no oracle for them.
See spec Non-Goals.
```

- [ ] **Step 5: Commit**

```bash
git add crates/ppvm-stim/tests/regen-stim/src/regen_stim/noise_sweeps.py
git add crates/ppvm-stim/tests/data/generated/noise_sweeps/
git commit -m "test(stim): generate noise_sweeps/ corpus with per-channel sweeps"
```

---

### Task 9: regen-stim dialect — generate `generated/dialect/` corpus

**Files:**
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/dialect.py`
- Create: `crates/ppvm-stim/tests/data/generated/dialect/<many>.stim` and `<many>.expected.json`
- Create: `crates/ppvm-stim/tests/data/generated/dialect/README.md`

ppvm-specific dialect: `I[R_X(theta=…)]`, `I[R_Y(theta=…)]`, `I[R_Z(theta=…)]`, `I[U3(theta=…, phi=…, lambda=…)]`, `S[T]`, `S_DAG[T]`. Stim cannot simulate these. Each fixture combines one or more dialect instructions in a circuit whose expected outcome is hand-derivable. Mostly deterministic; distribution mode for cases with measurement randomness, recorded directly from ppvm.

- [ ] **Step 1: Write `dialect.py`**

Create `crates/ppvm-stim/tests/regen-stim/src/regen_stim/dialect.py`:

```python
"""generated/dialect/: ppvm-specific dialect (I[R_X(...)], S[T], etc.).

Stim cannot simulate these so there's no oracle. Where the outcome is uniquely
determined we use deterministic mode; where there's measurement randomness we
record ppvm's output to lock down regression behavior.
"""

from __future__ import annotations

from . import core


# (name, source, mode_hint, bitstring_for_deterministic, num_shots_for_distribution)
DIALECT_FIXTURES = [
    # Deterministic: rotations that net out to identity or a Pauli.
    (
        "rx_pi_flips",
        "I[R_X(theta=1.0*pi)] 0\nM 0\n",
        "deterministic", [True], 0,
    ),
    (
        "ry_pi_flips",
        "I[R_Y(theta=1.0*pi)] 0\nM 0\n",
        "deterministic", [True], 0,
    ),
    (
        "rz_pi_no_flip",
        "I[R_Z(theta=1.0*pi)] 0\nM 0\n",
        "deterministic", [False], 0,
    ),
    (
        "rx_2pi_no_flip",
        "I[R_X(theta=2.0*pi)] 0\nM 0\n",
        "deterministic", [False], 0,
    ),
    (
        "rx_then_inverse",
        "I[R_X(theta=0.5*pi)] 0\nI[R_X(theta=1.5*pi)] 0\nM 0\n",
        "deterministic", [False], 0,
    ),
    (
        "u3_x_equivalent",
        "I[U3(theta=1.0*pi, phi=0.0, lambda=0.0)] 0\nM 0\n",
        "deterministic", [True], 0,
    ),
    (
        "t_then_t_dag_identity",
        "X 0\nS[T] 0\nS_DAG[T] 0\nM 0\n",
        "deterministic", [True], 0,
    ),
    (
        "t_eight_times_identity",
        "X 0\nS[T] 0\nS[T] 0\nS[T] 0\nS[T] 0\nS[T] 0\nS[T] 0\nS[T] 0\nS[T] 0\nM 0\n",
        "deterministic", [True], 0,
    ),
    # Distribution: rotations into a superposition.
    (
        "rx_half_pi_random",
        "I[R_X(theta=0.5*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "ry_half_pi_random",
        "I[R_Y(theta=0.5*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "u3_h_equivalent_random",
        "I[U3(theta=0.5*pi, phi=0.0, lambda=1.0*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "h_then_t_random",
        "H 0\nS[T] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "rx_pi_quarter_random",
        "I[R_X(theta=0.25*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
    # Combo of multiple dialect ops.
    (
        "rx_ry_rz_combo",
        "I[R_X(theta=0.5*pi)] 0\nI[R_Y(theta=0.5*pi)] 0\nI[R_Z(theta=0.5*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
    (
        "u3_arbitrary",
        "I[U3(theta=0.5*pi, phi=0.5*pi, lambda=0.5*pi)] 0\nM 0\n",
        "distribution", None, 256,
    ),
]


def run() -> int:
    paths = core.CorpusPaths.default()
    failures: list[str] = []
    written = 0
    for name, src, mode, bitstring, num_shots in DIALECT_FIXTURES:
        if mode == "deterministic":
            meta = core.FixtureMeta(
                name=name, category="generated/dialect", source=src, test_num_shots=1,
            )
            try:
                # We override write_deterministic_fixture's bitstring with the
                # caller-asserted one; simpler is: write json by hand here.
                core.write_deterministic_fixture(meta, paths)
                # Verify the recorded bitstring matches the asserted one.
                import json
                json_path = paths.category_dir("generated/dialect") / f"{name}.expected.json"
                got = json.loads(json_path.read_text())["bitstring"]
                if got != bitstring:
                    failures.append(
                        f"{name}: ppvm produced {got}, asserted {bitstring} — "
                        "either ppvm has a bug or the asserted bitstring is wrong"
                    )
                written += 1
            except Exception as e:
                failures.append(f"{name}: {e}")
        else:  # distribution
            # No Stim oracle. Run ppvm at seed=0 and record the means directly,
            # bypassing the seed-search loop in core.write_distribution_fixture.
            ppvm_run = core.run_ppvm(src, num_shots=num_shots, seed=0)
            payload = {
                "mode": "distribution",
                "num_shots": num_shots,
                "ppvm_seed": 0,
                "ppvm_bit_means": ppvm_run.bit_means,
                "note": "no Stim oracle (ppvm dialect); means recorded directly from ppvm",
            }
            cat_dir = paths.category_dir("generated/dialect")
            cat_dir.mkdir(parents=True, exist_ok=True)
            (cat_dir / f"{name}.stim").write_text(src)
            import json
            (cat_dir / f"{name}.expected.json").write_text(
                json.dumps(payload, indent=2) + "\n"
            )
            written += 1

    print(f"regen-stim dialect: wrote {written} fixtures")
    if failures:
        print("regen-stim dialect: failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0
```

- [ ] **Step 2: Run, verify, document**

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim dialect
cd ../../../..
cargo test -p ppvm-stim --test stim_corpus
```

Expected: 15 fixtures green.

If a deterministic fixture's asserted bitstring doesn't match ppvm's output, investigate carefully — either the asserted outcome is wrong (e.g. you confused R_Y with R_X), or there's a real bug in ppvm's rotation gates. Re-verify the math by hand against the ppvm gate reference in `crates/ppvm-tableau/src/gates/rot1.rs`.

Write `crates/ppvm-stim/tests/data/generated/dialect/README.md`:

```markdown
# generated/dialect/

ppvm-specific dialect instructions: `I[R_X(theta=…)]`, `I[R_Y(theta=…)]`,
`I[R_Z(theta=…)]`, `I[U3(theta=…, phi=…, lambda=…)]`, `S[T]`, `S_DAG[T]`.

Stim cannot simulate these — there is no oracle here. Deterministic-mode
fixtures are hand-derivable (e.g. `R_X(π)` flips |0⟩ to |1⟩); distribution-mode
fixtures record ppvm's output directly to lock down regression behavior.

This category catches ppvm dialect bugs only by detecting drift from
previously-recorded ppvm output, not by cross-checking against another simulator.

Generated by:

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim dialect
```
```

- [ ] **Step 3: Commit**

```bash
git add crates/ppvm-stim/tests/regen-stim/src/regen_stim/dialect.py
git add crates/ppvm-stim/tests/data/generated/dialect/
git commit -m "test(stim): generate dialect/ corpus for ppvm-specific instructions"
```

---

### Task 10: regen-stim random — generate `generated/random/` corpus

**Files:**
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/random_walk.py`
- Create: `crates/ppvm-stim/tests/data/generated/random/<many>.stim` and `<many>.expected.json`
- Create: `crates/ppvm-stim/tests/data/generated/random/README.md`

Random sequences of supported instructions. Three regimes (`clifford-only`, `+noise`, `+readout`), four qubit counts {2, 4, 8, 16}, three lengths {10, 50, 200}, eight RNG seeds. Cross-checked against Stim for the regimes that don't include dialect instructions.

- [ ] **Step 1: Write `random_walk.py`**

Create `crates/ppvm-stim/tests/regen-stim/src/regen_stim/random_walk.py`:

```python
"""generated/random/: random-walk sequences of supported instructions."""

from __future__ import annotations

import random as pyrand

from . import core

CLIFFORD_GATES_1Q = ["H", "S", "S_DAG", "X", "Y", "Z", "SQRT_X", "SQRT_Y"]
CLIFFORD_GATES_2Q = ["CX", "CY", "CZ"]
NOISE_GATES_1Q = [
    ("DEPOLARIZE1", lambda: f"({pyrand.choice([0.001, 0.01, 0.1])})"),
    ("X_ERROR", lambda: f"({pyrand.choice([0.001, 0.01, 0.1])})"),
    ("Y_ERROR", lambda: f"({pyrand.choice([0.001, 0.01, 0.1])})"),
    ("Z_ERROR", lambda: f"({pyrand.choice([0.001, 0.01, 0.1])})"),
]


def gen_program(n_qubits: int, n_instructions: int, regime: str, seed: int) -> str:
    pyrand.seed(seed)
    lines = [f"R {' '.join(str(i) for i in range(n_qubits))}"]
    for _ in range(n_instructions):
        roll = pyrand.random()
        if regime == "clifford-only":
            if roll < 0.6 or n_qubits < 2:
                gate = pyrand.choice(CLIFFORD_GATES_1Q)
                q = pyrand.randrange(n_qubits)
                lines.append(f"{gate} {q}")
            else:
                gate = pyrand.choice(CLIFFORD_GATES_2Q)
                a, b = pyrand.sample(range(n_qubits), 2)
                lines.append(f"{gate} {a} {b}")
        elif regime == "+noise":
            if roll < 0.30:  # ~30% noise as spec specifies high-density
                gate, args_fn = pyrand.choice(NOISE_GATES_1Q)
                q = pyrand.randrange(n_qubits)
                lines.append(f"{gate}{args_fn()} {q}")
            elif roll < 0.7 or n_qubits < 2:
                gate = pyrand.choice(CLIFFORD_GATES_1Q)
                q = pyrand.randrange(n_qubits)
                lines.append(f"{gate} {q}")
            else:
                gate = pyrand.choice(CLIFFORD_GATES_2Q)
                a, b = pyrand.sample(range(n_qubits), 2)
                lines.append(f"{gate} {a} {b}")
        elif regime == "+readout":
            if roll < 0.10:
                p = pyrand.choice([0.001, 0.01, 0.1])
                q = pyrand.randrange(n_qubits)
                lines.append(f"M({p}) {q}")
            elif roll < 0.30:
                gate, args_fn = pyrand.choice(NOISE_GATES_1Q)
                q = pyrand.randrange(n_qubits)
                lines.append(f"{gate}{args_fn()} {q}")
            elif roll < 0.7 or n_qubits < 2:
                gate = pyrand.choice(CLIFFORD_GATES_1Q)
                q = pyrand.randrange(n_qubits)
                lines.append(f"{gate} {q}")
            else:
                gate = pyrand.choice(CLIFFORD_GATES_2Q)
                a, b = pyrand.sample(range(n_qubits), 2)
                lines.append(f"{gate} {a} {b}")
        else:
            raise ValueError(f"unknown regime {regime}")
    lines.append(f"M {' '.join(str(i) for i in range(n_qubits))}")
    return "\n".join(lines) + "\n"


REGIMES = ["clifford-only", "+noise", "+readout"]
QUBIT_COUNTS = [2, 4, 8, 16]
LENGTHS = [10, 50, 200]
SEEDS = list(range(8))


def run() -> int:
    paths = core.CorpusPaths.default()
    failures: list[str] = []
    written = 0
    for regime in REGIMES:
        for n in QUBIT_COUNTS:
            for length in LENGTHS:
                for seed in SEEDS:
                    src = gen_program(n, length, regime, seed)
                    name = f"{regime.replace('+','plus_').replace('-','_')}_n{n}_len{length}_s{seed}"
                    meta = core.FixtureMeta(
                        name=name, category="generated/random",
                        source=src, test_num_shots=128,
                    )
                    try:
                        core.write_distribution_fixture(meta, paths); written += 1
                    except Exception as e:
                        failures.append(f"{name}: {e}")
    print(f"regen-stim random: wrote {written} fixtures")
    if failures:
        print("regen-stim random: failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0
```

The full sweep is 3 × 4 × 3 × 8 = 288 — too many; the spec asks for ~30–40. To stay in budget, sub-sample:

- Edit `REGIMES`, `QUBIT_COUNTS`, `LENGTHS`, `SEEDS` to a subset that produces ~36 total. Suggested subset: `REGIMES=["clifford-only", "+noise", "+readout"]`, `QUBIT_COUNTS=[2, 4, 8]`, `LENGTHS=[10, 50]`, `SEEDS=[0, 1]`. Total: 3*3*2*2 = 36.

Make that edit now and re-run.

- [ ] **Step 2: Run, verify, document**

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim random
cd ../../../..
cargo test -p ppvm-stim --test stim_corpus
```

Expected: ~36 fixtures green. Test runtime: ~5 seconds added at 128 shots each.

Write `crates/ppvm-stim/tests/data/generated/random/README.md`:

```markdown
# generated/random/

Random-walk programs over supported phase-1 instructions. Three regimes:
- `clifford-only`: just Clifford gates.
- `+noise`: ~30% of instructions are noise channels.
- `+readout`: also includes M(p) readout-noise measurements.

Generated by:

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim random
```

Fixture count: 36 (3 regimes × 3 qubit counts × 2 lengths × 2 seeds).
Distribution mode at `num_shots=128`.

The high-noise-density regime complements `noise_sweeps/` (which has one
channel per circuit) by producing fixtures with many mixed channels per
circuit.
```

- [ ] **Step 3: Commit**

```bash
git add crates/ppvm-stim/tests/regen-stim/src/regen_stim/random_walk.py
git add crates/ppvm-stim/tests/data/generated/random/
git commit -m "test(stim): generate random/ corpus with random-walk programs"
```

---

### Task 11: regen-stim unsupported — generate `unsupported/` corpus

**Files:**
- Create: `crates/ppvm-stim/tests/regen-stim/src/regen_stim/unsupported.py`
- Create: `crates/ppvm-stim/tests/data/unsupported/<many>.stim` and `<many>.expected.json`
- Create: `crates/ppvm-stim/tests/data/unsupported/README.md`

Per-instruction templates: prepare a state → apply the unsupported gate → measure. Stim reference is pre-recorded so phase-2 flipping requires only adding ppvm-side fields.

The phase-1-unsupported instruction list (from spec section "unsupported/"): `SWAP`, `ISWAP`, `ISWAP_DAG`, `SQRT_XX`, `SQRT_YY`, `SQRT_ZZ`, `CXSWAP`, `SWAPCX`, `XCX`, `XCY`, `XCZ`, `YCX`, `YCY`, `YCZ`, `C_XYZ`, `C_ZYX`, `H_XY`, `H_YZ`, `MX`, `MY`, `MRX`, `MRY`, `MXX`, `MYY`, `MZZ`, `MPP`, `HERALDED_ERASE`, `HERALDED_PAULI_CHANNEL_1`, `CORRELATED_ERROR`, `ELSE_CORRELATED_ERROR`. (Two existing fixtures from Task 2 already cover `SWAP` and `MX`; the script should be idempotent — overwriting existing fixtures is fine, but make sure the asserted instruction name still matches.)

- [ ] **Step 1: Write `unsupported.py`**

Create `crates/ppvm-stim/tests/regen-stim/src/regen_stim/unsupported.py`:

```python
"""unsupported/: one fixture per phase-1-unsupported instruction.

Each fixture: prep → apply unsupported gate → measure. Stim reference is
pre-recorded so phase-2 lifting only needs ppvm fields added.
"""

from __future__ import annotations

from . import core

# (instruction_name, source_template, fixture_name)
UNSUPPORTED_FIXTURES: list[tuple[str, str, str]] = [
    ("SWAP", "X 0\nSWAP 0 1\nM 0 1\n", "swap_unsupported"),
    ("ISWAP", "X 0\nISWAP 0 1\nM 0 1\n", "iswap_unsupported"),
    ("ISWAP_DAG", "X 0\nISWAP_DAG 0 1\nM 0 1\n", "iswap_dag_unsupported"),
    ("SQRT_XX", "H 0 1\nSQRT_XX 0 1\nM 0 1\n", "sqrt_xx_unsupported"),
    ("SQRT_YY", "H 0 1\nSQRT_YY 0 1\nM 0 1\n", "sqrt_yy_unsupported"),
    ("SQRT_ZZ", "H 0 1\nSQRT_ZZ 0 1\nM 0 1\n", "sqrt_zz_unsupported"),
    ("CXSWAP", "X 0\nCXSWAP 0 1\nM 0 1\n", "cxswap_unsupported"),
    ("SWAPCX", "X 0\nSWAPCX 0 1\nM 0 1\n", "swapcx_unsupported"),
    ("XCX", "X 0\nXCX 0 1\nM 0 1\n", "xcx_unsupported"),
    ("XCY", "X 0\nXCY 0 1\nM 0 1\n", "xcy_unsupported"),
    ("XCZ", "X 0\nXCZ 0 1\nM 0 1\n", "xcz_unsupported"),
    ("YCX", "X 0\nYCX 0 1\nM 0 1\n", "ycx_unsupported"),
    ("YCY", "X 0\nYCY 0 1\nM 0 1\n", "ycy_unsupported"),
    ("YCZ", "X 0\nYCZ 0 1\nM 0 1\n", "ycz_unsupported"),
    ("C_XYZ", "H 0\nC_XYZ 0\nM 0\n", "c_xyz_unsupported"),
    ("C_ZYX", "H 0\nC_ZYX 0\nM 0\n", "c_zyx_unsupported"),
    ("H_XY", "H 0\nH_XY 0\nM 0\n", "h_xy_unsupported"),
    ("H_YZ", "H 0\nH_YZ 0\nM 0\n", "h_yz_unsupported"),
    ("MX", "H 0\nMX 0\n", "mx_unsupported"),
    ("MY", "H 0\nMY 0\n", "my_unsupported"),
    ("MRX", "H 0\nMRX 0\n", "mrx_unsupported"),
    ("MRY", "H 0\nMRY 0\n", "mry_unsupported"),
    ("MXX", "H 0 1\nMXX 0 1\n", "mxx_unsupported"),
    ("MYY", "H 0 1\nMYY 0 1\n", "myy_unsupported"),
    ("MZZ", "H 0 1\nMZZ 0 1\n", "mzz_unsupported"),
    ("MPP", "H 0 1\nMPP X0*X1\n", "mpp_unsupported"),
    ("HERALDED_ERASE", "X 0\nHERALDED_ERASE(0.1) 0\nM 0\n", "heralded_erase_unsupported"),
    (
        "HERALDED_PAULI_CHANNEL_1",
        "X 0\nHERALDED_PAULI_CHANNEL_1(0.05, 0.05, 0.05, 0.05) 0\nM 0\n",
        "heralded_pauli_channel_1_unsupported",
    ),
    (
        "CORRELATED_ERROR",
        "CORRELATED_ERROR(0.1) X0 X1\nM 0 1\n",
        "correlated_error_unsupported",
    ),
    (
        "ELSE_CORRELATED_ERROR",
        "CORRELATED_ERROR(0.1) X0\nELSE_CORRELATED_ERROR(0.1) X1\nM 0 1\n",
        "else_correlated_error_unsupported",
    ),
]


def run() -> int:
    paths = core.CorpusPaths.default()
    failures: list[str] = []
    written = 0
    for instruction, source, name in UNSUPPORTED_FIXTURES:
        meta = core.FixtureMeta(
            name=name, category="unsupported", source=source, test_num_shots=0,
        )
        try:
            core.write_unsupported_fixture(
                meta, paths, awaiting_phase2_instruction=instruction
            )
            written += 1
        except Exception as e:
            failures.append(f"{name}: {e}")
    print(f"regen-stim unsupported: wrote {written} fixtures")
    if failures:
        print("regen-stim unsupported: failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0
```

- [ ] **Step 2: Run, verify, document**

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim unsupported
cd ../../../..
cargo test -p ppvm-stim --test stim_corpus
```

Expected: ~31 fixtures (overwriting the 2 existing `swap_unsupported`, `mx_unsupported` is fine — same content). All harness tests green; the harness asserts `awaiting_phase2_instruction` matches `NormalizeError::Unsupported.name`, which hits all 31.

If a fixture errors with "stim refused to parse" or similar, the source template is malformed for that instruction — open Stim's docs and adjust the template. (Some instructions like `MPP` use Pauli-target syntax that's easy to get wrong.)

Note: when the Stim source is pasted into ppvm's parser, ppvm should reject it with `NormalizeError::Unsupported { name: <instruction> }`. If the parser rejects it earlier (e.g. unknown gate name), that means the parser doesn't know about the instruction at all — we'd need to add it to ppvm's `GateName`/`MeasureName` enum first. The Phase-1 parser already accepts every Stim instruction listed above (per `crates/ppvm-stim/src/parser/ast.rs`). If a fixture fails, double-check by parsing it manually:

```bash
cargo run --example parse_one --release -- crates/ppvm-stim/tests/data/unsupported/<name>.stim
```

(If the example doesn't exist, just run `cargo test -p ppvm-stim --test stim_corpus` and read the failure message.)

Write `crates/ppvm-stim/tests/data/unsupported/README.md`:

```markdown
# unsupported/

One fixture per phase-1-unsupported Stim instruction. Each fixture flips from
`mode: "unsupported"` to `mode: "distribution"` (or `deterministic`) when
phase-2 implements that instruction.

Phase-2 flip workflow:

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim refresh ../data/unsupported/<name>.stim
```

`refresh` reads the existing JSON, sees that phase-2 now supports the
instruction (because ppvm's normalize::to_tableau succeeds), runs the
seed-search loop against the pre-recorded `stim_bit_means`, and writes a
new JSON with `mode: "distribution"`. The `.stim` source itself never
changes.

Generated by:

```bash
uv run regen-stim unsupported
```

The `regen-stim codes` subcommand also auto-routes fixtures here when
`stim gen` emits a phase-1-unsupported instruction.
```

- [ ] **Step 3: Commit**

```bash
git add crates/ppvm-stim/tests/regen-stim/src/regen_stim/unsupported.py
git add crates/ppvm-stim/tests/data/unsupported/
git commit -m "test(stim): generate unsupported/ corpus with one fixture per phase-1 gap"
```

---

### Task 12: Documentation — top-level + per-category

Most of the per-category READMEs already landed in their respective tasks. This task tightens up the top-level README and AGENTS.md mention.

**Files:**
- Modify: `crates/ppvm-stim/tests/data/README.md` (final cleanup, link to regen-stim README).
- Modify: `AGENTS.md` (workspace-level — add a one-liner about the corpus + regen workflow).

- [ ] **Step 1: Cross-reference the regen workflow in AGENTS.md**

Run:

```bash
grep -n "stim_corpus\|regen-stim\|tests/data" AGENTS.md
```

If there's no mention of the corpus harness, add a paragraph under whatever section talks about ppvm-stim. Insert at the appropriate location:

```markdown
### `crates/ppvm-stim` test corpus

Tests under `crates/ppvm-stim/tests/data/` are committed `.stim` + `.expected.json`
pairs consumed by `tests/stim_corpus.rs`. The harness asserts ppvm's output
matches the committed reference bit-for-bit. Cross-check against `quantumlib/Stim`
happens at regen time, not at test time:

```
cd crates/ppvm-stim/tests/regen-stim
uv sync
uv run regen-stim all
```

When phase-2 lifts a restriction, `uv run regen-stim refresh
../data/unsupported/<name>.stim` flips that fixture from "expected to fail
normalize" to "expected to match Stim's pre-recorded distribution".
```

If an `AGENTS.md` section already covers ppvm-stim testing, just add a one-line cross-reference instead of duplicating.

- [ ] **Step 2: Final top-level corpus README polish**

Re-read `crates/ppvm-stim/tests/data/README.md` (written in Task 2). Ensure:

- Each category subdir is mentioned with its `README.md`.
- The link to the spec is correct (`docs/superpowers/specs/2026-04-28-ppvm-stim-test-corpus-design.md`).
- The "Regenerating" section points at `crates/ppvm-stim/tests/regen-stim/README.md` for setup details (don't duplicate setup instructions).

If anything's missing, edit it now.

- [ ] **Step 3: Run the full workspace test suite**

```bash
cargo test --workspace
```

Expected: every test green. The corpus harness adds ~30–40 seconds per spec.

- [ ] **Step 4: Commit**

```bash
git add crates/ppvm-stim/tests/data/README.md AGENTS.md
git commit -m "docs(stim): document corpus layout and regen workflow"
```

---

## Final verification

- [ ] **Run the full Rust test suite**

```bash
cargo test --workspace
```

Expected: all green; no flakes.

- [ ] **Run the regen-stim unit tests**

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run pytest test/
```

Expected: all green.

- [ ] **Smoke-run `regen-stim verify` on a handful of distribution fixtures**

```bash
uv run regen-stim verify ../data/edge_cases/depolarize_smoke.stim
uv run regen-stim verify ../data/generated/noise_sweeps/depolarize1_n1_p001.stim
uv run regen-stim verify ../data/generated/codes/surface_code_rotated_memory_z_d3_r3_p001.stim
```

Expected: each prints `verify: <path> OK`. Cross-check fidelity holds.

- [ ] **Inspect total fixture count**

```bash
cd ../../../..
find crates/ppvm-stim/tests/data -name "*.stim" | wc -l
find crates/ppvm-stim/tests/data -name "*.expected.json" | wc -l
```

Expected: equal counts in the 200–250 range. (If the random/random_walk subsampling was kept tight, the lower bound is closer to 180; if all 8 random seeds are used, closer to 250.)

- [ ] **Workspace runtime check**

```bash
cargo test -p ppvm-stim --release --test stim_corpus -- --report-time
```

Expected: total runtime ~30–40 seconds (release build). If it's significantly higher (>2 min), reconsider per-category shot counts per the spec's "Per-Category Shot Defaults" guidance.

- [ ] **Push the branch**

```bash
git push origin david/44-ppvm-stim
```

(Pause here and confirm with the user before opening the PR.)

---

## Implementer notes

- **Editable `ppvm` install in regen-stim.** `[tool.uv.sources] ppvm = { path = "../../../../ppvm-python", editable = true }` makes `uv sync` use the in-tree ppvm-python source. The native extension (`ppvm-python-native` Rust crate) must be built first via `uv run --project ppvm-python --group dev maturin develop --uv` from the repo root. If `uv sync` errors with "ppvm-python-native: cannot find module", run that maturin command and re-sync. This indirection is intentional: regen-stim depends on `ppvm`, ppvm depends on `ppvm-python-native`, the native extension is a Rust crate built by maturin.
- **Bit-exact f64 compare relies on stable RNG semantics.** `GeneralizedTableau::new_with_seed` must always produce the same sample stream from the same `(seed, num_shots, program)` triple. If a future ppvm refactor changes RNG ordering inside the executor, every distribution-mode fixture will fail loudly — that's a feature, not a bug. The fix is to run `uv run regen-stim refresh ../data/...` per-fixture and commit the new means.
- **Stim version drift.** When Stim ships a new release that changes its RNG, `regen-stim verify` will start failing on existing fixtures. The fix is `uv run regen-stim all` and commit the regenerated JSONs (the `stim_version` field in each JSON makes the cause obvious).
- **Subsampling in `random_walk.py`.** The spec asks for ~30–40 fixtures; the full sweep is 288. Keep `QUBIT_COUNTS=[2,4,8]`, `LENGTHS=[10,50]`, `SEEDS=[0,1]` to land at 36. If you need more diversity later, expand the seed list — adding seeds is the cheapest knob.
- **`HERALDED_*` and Pauli-target syntax in `unsupported/`.** Some Stim instructions use exotic target syntax (`MPP X0*X1`, `CORRELATED_ERROR(p) X0 Y1`). The `unsupported/` fixtures need each instruction to (a) parse in ppvm's parser and (b) be rejected by `normalize::to_tableau`. If ppvm's parser doesn't accept a particular target form, that's a parser bug and lives outside this plan. The fallback is to comment out the offending instruction in `UNSUPPORTED_FIXTURES`, file an issue, and revisit.
- **Two-tier verification rationale.** ppvm and Stim use different internal RNGs, so bit-exact comparison between them is impossible. Statistical comparison at test time would introduce flake risk. The two-tier scheme — Stim cross-check at regen time, bit-exact ppvm-vs-self compare at test time — keeps both signals: "did our committed reference drift from Stim?" (regen) and "did ppvm's behavior change?" (test).
- **Why `serde_json` over `toml` for fixture metadata.** Considered TOML to avoid the dev-dep, but the regen script is Python and `json` is the path of least friction there. The dev-dep cost is one tiny crate; the upside is no impedance mismatch between the regen tool and the test harness.
- **`max_qubit_in_source` in `core.py` reads via `stim.Circuit`.** Using Stim's parser is intentional — Stim is the spec for `.stim` syntax, and the regen tool's job is to cross-check against it. The function's only consumer is `run_ppvm`, which uses it to size the tableau.
- **Loss-mode fixtures are out of scope.** Per spec Non-Goals: Stim doesn't simulate `I_ERROR[loss]` / `I_ERROR[correlated_loss]`, so loss has no oracle. Loss is exercised by hand-written tests in `crates/ppvm-stim/tests/executor.rs`. The harness rejects `None` measurement results in distribution and deterministic modes by design.
