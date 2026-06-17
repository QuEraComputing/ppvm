// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_runtime::config::indexmap::ByteFxHashF64;
use ppvm_stim::{execute, parse_extended};
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<1>, usize>;

#[test]
fn run_populates_measurement_record() {
    let prog = parse_extended("X 0\nM 0 1").unwrap();
    let mut t: Tab = GeneralizedTableau::new(2, 1e-12);
    let results = execute(&prog, &mut t).unwrap();
    assert_eq!(results.len(), 2);
    // The executor's returned results must equal the tableau's record exactly:
    // `measure` records once per measurement and the executor does not
    // separately push, so there is no double- or under-recording.
    assert_eq!(t.current_measurement_record(), results.as_slice());
    assert_eq!(results, vec![Some(true), Some(false)]);
}
