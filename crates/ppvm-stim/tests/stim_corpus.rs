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
use ppvm_stim::{NormalizeError, ParseError, execute, normalize, parse, sample};
use ppvm_tableau::prelude::*;
use serde::Deserialize;
use walkdir::WalkDir;

type Tab = GeneralizedTableau<ByteFxHashF64<8>, usize>;

/// Storage `<8>` (64 bits per Pauli word) caps the tableau at 64 qubits.
/// `regen-stim codes` skips circuits that exceed this; bumping the budget
/// requires changing the storage parameter and regenerating every committed
/// `ppvm_bit_means` (because the bit-exact compare is sensitive to tableau
/// shape).
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
        #[allow(dead_code)]
        stim_seed: Option<u64>,
        #[serde(default)]
        #[allow(dead_code)]
        stim_num_shots: Option<usize>,
        #[serde(default)]
        #[allow(dead_code)]
        stim_bit_means: Option<Vec<f64>>,
        #[serde(default)]
        #[allow(dead_code)]
        tolerance_sigma_at_regen: Option<f64>,
        #[serde(default)]
        #[allow(dead_code)]
        stim_version: Option<String>,
    },
    Unsupported {
        awaiting_phase2_instruction: String,
        // Pre-recorded for phase-2 flip; harness does not use them.
        #[serde(default)]
        #[allow(dead_code)]
        stim_seed: Option<u64>,
        #[serde(default)]
        #[allow(dead_code)]
        stim_num_shots: Option<usize>,
        #[serde(default)]
        #[allow(dead_code)]
        stim_bit_means: Option<Vec<f64>>,
        #[serde(default)]
        #[allow(dead_code)]
        tolerance_sigma_at_regen: Option<f64>,
        #[serde(default)]
        #[allow(dead_code)]
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
    // Unsupported mode: rejection may come from either parse (unknown
    // instruction) or normalize. Both are "phase-1 doesn't support this"
    // signals; phase-2 will lift either.
    if let Expected::Unsupported {
        awaiting_phase2_instruction,
        ..
    } = expected
    {
        match parse(src) {
            Err(ParseError::UnknownInstruction { name, .. }) => {
                if name == *awaiting_phase2_instruction {
                    return Ok(());
                }
                return Err(format!(
                    "{label}: expected Unsupported({awaiting_phase2_instruction}), parser rejected '{name}'"
                ));
            }
            // Other ParseError variants (Syntax, ArgCount, TargetCount) are
            // accepted as phase-1 rejection — they fire for instructions
            // whose name the parser knows but whose operand syntax it
            // doesn't yet support (e.g. MPP's `X0*X1` Pauli-string targets).
            Err(_) => return Ok(()),
            Ok(prog) => match normalize::to_tableau(&prog) {
                Err(NormalizeError::Unsupported { name, .. }) => {
                    if name != *awaiting_phase2_instruction {
                        return Err(format!(
                            "{label}: expected Unsupported({awaiting_phase2_instruction}), got Unsupported({name})"
                        ));
                    }
                    return Ok(());
                }
                Err(other) => {
                    return Err(format!("{label}: expected Unsupported, got {other:?}"));
                }
                Ok(_) => {
                    return Err(format!(
                        "{label}: expected Unsupported({awaiting_phase2_instruction}), but normalize succeeded"
                    ));
                }
            },
        }
    }

    let prog = parse(src).map_err(|e| format!("{label}: parse failed: {e}"))?;

    match expected {
        Expected::Unsupported { .. } => unreachable!("handled above"),

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
            // sample() builds a fresh tableau per shot via this closure; the
            // closure increments `next_seed` so shot k uses `ppvm_seed + k`.
            // The regen tool MUST replicate this exact sequence (same start,
            // same +1 step, same per-bit `sum / num_shots` summation order) —
            // any deviation makes the bit-exact f64 compare below fail.
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
