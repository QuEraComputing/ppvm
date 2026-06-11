// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::path::PathBuf;

use ppvm_runtime::config::fxhash::ByteF64;
use ppvm_stim::{ExecError, execute, parse_extended};
use ppvm_tableau::prelude::*;

#[derive(Debug, Clone, Copy)]
enum Expect {
    /// File parses, prepares, executes.
    Ok,
    /// File parses, but prepare/execute must fail with `Unsupported(name)`.
    ExecUnsupported(&'static str),
    /// File should fail at parse time (e.g. uses `rec[-k]` targets).
    ParseFails,
}

const CASES: &[(&str, Expect)] = &[
    ("x_only.stim", Expect::Ok),
    ("bell_pair.stim", Expect::Ok),
    ("ghz.stim", Expect::Ok),
    ("repeat_block.stim", Expect::Ok),
    ("depolarize_smoke.stim", Expect::Ok),
    ("swap_unsupported.stim", Expect::ExecUnsupported("SWAP")),
    ("mx_unsupported.stim", Expect::ExecUnsupported("MX")),
    // The file uses `rec[-k]` targets on `DETECTOR` / `OBSERVABLE_INCLUDE`.
    // Phase-1 cannot represent measurement-record targets, but the parser
    // tolerates non-numeric tokens on annotations (which are no-ops in our
    // pipeline) so the file parses, prepares, and executes cleanly.
    ("repetition_code_d3_r3.stim", Expect::Ok),
    // `CX rec[-1] 1` is classically-controlled feedback. Gates do not
    // tolerate non-numeric targets, so the parser rejects this file.
    ("feedback_cx_unsupported.stim", Expect::ParseFails),
];

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
}

fn read(name: &str) -> String {
    std::fs::read_to_string(data_dir().join(name))
        .unwrap_or_else(|e| panic!("missing fixture {name}: {e}"))
}

#[test]
fn corpus_table_covers_every_file() {
    let dir = data_dir();
    let mut found: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| {
            let e = e.ok()?;
            let p = e.path();
            if p.extension().and_then(|s| s.to_str()) == Some("stim") {
                Some(p.file_name().unwrap().to_string_lossy().into_owned())
            } else {
                None
            }
        })
        .collect();
    found.sort();

    let mut declared: Vec<String> = CASES.iter().map(|(n, _)| (*n).to_string()).collect();
    declared.sort();

    assert_eq!(found, declared, "every .stim file must appear in CASES");
}

#[test]
fn corpus_obeys_expectations() {
    type Tab = GeneralizedTableau<ByteF64<8>, usize>;

    for (name, expect) in CASES {
        let src = read(name);
        let parsed = parse_extended(&src);
        match (expect, parsed) {
            (Expect::ParseFails, Ok(_)) => {
                panic!("{name}: expected parse failure, but parse succeeded");
            }
            (Expect::ParseFails, Err(_)) => continue,
            (Expect::Ok, Err(e)) | (Expect::ExecUnsupported(_), Err(e)) => {
                panic!("{name}: parse failed unexpectedly: {e}");
            }
            (Expect::Ok, Ok(prog)) => {
                let mut tab: Tab = GeneralizedTableau::new(64, 1e-10);
                execute(&prog, &mut tab).unwrap_or_else(|e| panic!("{name}: execute failed: {e}"));
            }
            (Expect::ExecUnsupported(expected_name), Ok(prog)) => {
                let mut tab: Tab = GeneralizedTableau::new(64, 1e-10);
                match execute(&prog, &mut tab) {
                    Err(ExecError::Unsupported { name: n, .. }) => {
                        assert_eq!(n, *expected_name, "{name}: wrong unsupported name");
                    }
                    Err(other) => {
                        panic!("{name}: expected Unsupported, got {other:?}")
                    }
                    Ok(_) => panic!("{name}: expected Unsupported, but execute succeeded"),
                }
            }
        }
    }
}
