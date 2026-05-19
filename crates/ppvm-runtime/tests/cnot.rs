// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use approx::assert_relative_eq;
use ppvm_runtime::prelude::*;

#[test]
fn test_cnot() {
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<4>> =
        PauliSum::builder().n_qubits(2).build();
    let pat: PauliPattern = "Z?*".into();
    state += ("ZZ", 1.0);

    state.rz(0, 1.1);
    state.ry(0, 2.1);
    state.rz(0, 1.1);

    state.rz(1, 1.1);
    state.ry(1, 2.1);
    state.rz(1, 1.1);

    state.cnot(0, 1);

    state.rx(0, 2.1);
    state.rx(1, 2.1);

    assert_relative_eq!(state.trace(&pat), 0.18803675917759355)
}
