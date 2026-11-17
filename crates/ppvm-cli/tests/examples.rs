// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Every shipped `.sst` example must compile and run to completion.

use std::path::Path;

#[test]
fn every_example_compiles_and_runs() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples");
    let mut ran = 0;
    for entry in std::fs::read_dir(&dir).expect("examples directory") {
        let path = entry.expect("directory entry").path();
        if path.extension().is_none_or(|ext| ext != "sst") {
            continue;
        }
        let name = path.display();
        let module = ppvm_vihaco::load_module_file(path.to_str().expect("utf-8 path"))
            .unwrap_or_else(|e| panic!("compile {name}: {e:?}"));
        let shots = ppvm_vihaco::shots::run_shots(&module, 4, Some(1))
            .unwrap_or_else(|e| panic!("run {name}: {e:?}"));
        assert_eq!(shots.len(), 4, "{name}");
        ran += 1;
    }
    assert!(ran > 0, "no examples found in {}", dir.display());
}
