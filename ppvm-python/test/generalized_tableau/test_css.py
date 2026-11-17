# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

import numpy as np
import pytest

from ppvm import MeasurementResult, StimProgram, sample_stim

# Loss, leakage, a Hadamard-conjugated CZ, feedback and readout noise.
CIRCUIT = StimProgram.parse("""
RX 0
R 1 2 3
I_ERROR[loss](0.1) 0 1 2 3
CX 0 1
H 2
CZ 2 3
H 3
I_ERROR[leakage](0.05, 0.05) 1 3
I_ERROR[correlated_loss](0.02, 0.03, 0.2) 0 2
DEPOLARIZE2(0.05) 0 1
M(0.02) 0 1
CX rec[-1] 2
MX 2
MR 3
M 3
""")


def test_is_css():
    assert CIRCUIT.is_css
    assert not StimProgram.parse("S 0\nM 0").is_css
    assert not StimProgram.parse("CZ 0 1\nCZ 1 2\nCZ 0 2").is_css


def test_css_matches_tableau():
    n = 40_000
    a = sample_stim(CIRCUIT, num_shots=n, seed=1, as_numpy=True, method="css")
    b = sample_stim(CIRCUIT, num_shots=n, seed=2, as_numpy=True, method="tableau")
    # One-hot outcome / lost features; compare every frequency and pairwise joint frequency.
    fa = np.concatenate([a == 1, a == 2], axis=1).astype(float)
    fb = np.concatenate([b == 1, b == 2], axis=1).astype(float)
    ca, cb = fa.T @ fa / n, fb.T @ fb / n
    sigma = np.sqrt((ca * (1 - ca) + cb * (1 - cb)) / n) + 1 / n
    assert np.abs(ca - cb).max() > 0
    assert (np.abs(ca - cb) / sigma).max() < 5


def test_auto_uses_css_when_possible():
    css = sample_stim(CIRCUIT, num_shots=300, seed=3, as_numpy=True, method="css")
    auto = sample_stim(CIRCUIT, num_shots=300, seed=3, as_numpy=True)
    np.testing.assert_array_equal(auto, css)
    t = StimProgram.parse("H 0\nS 0\nH 0\nM 0")
    tableau = sample_stim(t, num_shots=300, seed=3, as_numpy=True, method="tableau")
    np.testing.assert_array_equal(sample_stim(t, num_shots=300, seed=3, as_numpy=True), tableau)


def test_css_output_format():
    bits = sample_stim(CIRCUIT, num_shots=1000, seed=4, as_numpy=True, method="css")
    assert bits.dtype == np.int8
    assert bits.shape == (1000, 5)
    assert bits.flags.writeable
    assert set(np.unique(bits)) == {0, 1, 2}
    shots = sample_stim(CIRCUIT, num_shots=1000, seed=4, method="css")
    assert shots == [[MeasurementResult(int(v)) for v in row] for row in bits]
    assert sample_stim(CIRCUIT, num_shots=0, as_numpy=True, method="css").shape == (0, 5)
    no_meas = StimProgram.parse("H 0")
    assert sample_stim(no_meas, num_shots=4, as_numpy=True, method="css").shape == (4, 0)


def test_css_seeds():
    a = sample_stim(CIRCUIT, num_shots=2000, seed=5, as_numpy=True, method="css")
    np.testing.assert_array_equal(
        a, sample_stim(CIRCUIT, num_shots=2000, seed=5, as_numpy=True, method="css")
    )
    assert (a != sample_stim(CIRCUIT, num_shots=2000, seed=6, as_numpy=True, method="css")).any()
    unseeded = sample_stim(CIRCUIT, num_shots=2000, as_numpy=True, method="css")
    assert (a != unseeded).any()


def test_css_beyond_the_tableau_qubit_limit():
    qubits = " ".join(map(str, range(3000)))
    prog = StimProgram.parse(f"X {qubits}\nI_ERROR[loss](1) 2999\nM {qubits}")
    bits = sample_stim(prog, num_shots=3, as_numpy=True)
    assert (bits[:, :-1] == 1).all()
    assert (bits[:, -1] == 2).all()


def test_css_method_rejects_other_circuits():
    with pytest.raises(ValueError, match="CSS"):
        sample_stim(StimProgram.parse("S 0\nM 0"), num_shots=1, method="css")
    with pytest.raises(ValueError, match="method"):
        sample_stim(CIRCUIT, num_shots=1, method="frame")  # ty: ignore[no-matching-overload]
