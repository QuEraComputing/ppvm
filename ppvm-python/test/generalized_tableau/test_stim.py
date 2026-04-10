import os
import tempfile

import pytest

from ppvm import GeneralizedTableau
from ppvm.generalized_tableau import MeasurementResult

# --- run_stim_string ---


def test_run_stim_string_single_measurement_zero():
    # Fresh |0> state → M returns ZERO
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("M 0")
    assert results == [MeasurementResult.ZERO]


def test_run_stim_string_x_then_measure():
    # X|0> = |1> → M returns ONE
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("X 0\nM 0")
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_measurement_order_follows_circuit():
    # M 0 1: qubit 0 was flipped, qubit 1 was not
    # Results must be [ONE, ZERO] in that order (circuit order, STIM convention)
    tab = GeneralizedTableau(2)
    results = tab.run_stim_string("X 0\nM 0 1")
    assert results == [MeasurementResult.ONE, MeasurementResult.ZERO]


def test_run_stim_string_multiple_m_instructions_appended_in_order():
    # Two separate M instructions: measurements from the first instruction
    # appear before those from the second, following STIM circuit order.
    tab = GeneralizedTableau(2)
    results = tab.run_stim_string("X 0\nM 1\nM 0")
    assert results == [MeasurementResult.ZERO, MeasurementResult.ONE]


def test_run_stim_string_double_measurement_same_qubit():
    # Measuring the same qubit twice appends two entries
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("X 0\nM 0\nM 0")
    assert results == [MeasurementResult.ONE, MeasurementResult.ONE]


def test_run_stim_string_no_measurements_returns_empty():
    tab = GeneralizedTableau(2)
    results = tab.run_stim_string("X 0\nX 1")
    assert results == []


def test_run_stim_string_reset():
    # R resets to |0>, so subsequent measurement is ZERO
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("X 0\nR 0\nM 0")
    assert results == [MeasurementResult.ZERO]


def test_run_stim_string_mr_resets_after_measure():
    # MR measures and then resets; second M should always be ZERO
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("X 0\nMR 0\nM 0")
    assert results == [MeasurementResult.ONE, MeasurementResult.ZERO]


def test_run_stim_string_cnot_entanglement():
    # H(0) + CX(0,1) creates a Bell state; both qubits must agree
    tab = GeneralizedTableau(2)
    results = tab.run_stim_string("H 0\nCX 0 1\nM 0 1")
    assert len(results) == 2
    assert results[0] == results[1]


def test_run_stim_string_cz_gate():
    tab = GeneralizedTableau(2)
    results = tab.run_stim_string("X 0\nX 1\nCZ 0 1\nM 0 1")
    assert results == [MeasurementResult.ONE, MeasurementResult.ONE]


def test_run_stim_string_h_xz_alias():
    # H_XZ is an alias for H
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("X 0\nH_XZ 0\nH_XZ 0\nM 0")
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_s_gate():
    # S|0> = |0>
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("S 0\nM 0")
    assert results == [MeasurementResult.ZERO]


def test_run_stim_string_t_gate_tag():
    # S[T] is the T gate; T|0> = |0>
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("S[T] 0\nM 0")
    assert results == [MeasurementResult.ZERO]


def test_run_stim_string_s_dag_t_tag():
    # S_DAG[T] is the T† gate; T T† |1> = |1>
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("X 0\nS[T] 0\nS_DAG[T] 0\nM 0")
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_sqrt_x():
    # SQRT_X SQRT_X = X
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("SQRT_X 0\nSQRT_X 0\nM 0")
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_sqrt_y():
    # SQRT_Y SQRT_Y = Y (up to global phase); Y|0> = i|1> → measures as 1
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("SQRT_Y 0\nSQRT_Y 0\nM 0")
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_i_rotation_x_pi():
    # I[R_X(theta=1.0*pi)]|0> = |1>
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("I[R_X(theta=1.0*pi)] 0\nM 0")
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_i_rotation_y_pi():
    # I[R_Y(theta=1.0*pi)]|0> = |1>
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("I[R_Y(theta=1.0*pi)] 0\nM 0")
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_i_rotation_z_pi():
    # I[R_Z(theta=1.0*pi)]|0> = |0> (Z rotation only adds phase)
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("I[R_Z(theta=1.0*pi)] 0\nM 0")
    assert results == [MeasurementResult.ZERO]


def test_run_stim_string_loss_channel():
    # I_ERROR[loss](1.0) loses the qubit with probability 1
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string("I_ERROR[loss](1.0) 0\nM 0")
    assert results == [MeasurementResult.LOST]


def test_run_stim_string_comments_and_blank_lines_ignored():
    # Comments (#) and blank lines must not affect execution
    circuit = """
# Prepare and measure
X 0

M 0
"""
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string(circuit)
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_tick_and_annotations_are_noops():
    # TICK, DETECTOR, OBSERVABLE_INCLUDE, QUBIT_COORDS, SHIFT_COORDS,
    # and MPAD are annotation-only and must not affect circuit outcomes.
    circuit = (
        "QUBIT_COORDS(0, 0) 0\nX 0\nTICK\nM 0\nDETECTOR rec[-1]\nOBSERVABLE_INCLUDE(0) rec[-1]\n"
    )
    tab = GeneralizedTableau(1)
    results = tab.run_stim_string(circuit)
    assert results == [MeasurementResult.ONE]


# --- run_stim_file ---


def test_run_stim_file_basic():
    circuit = "X 0\nM 0 1"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".stim", delete=False) as f:
        f.write(circuit)
        path = f.name
    try:
        tab = GeneralizedTableau(2)
        results = tab.run_stim_file(path)
        assert results == [MeasurementResult.ONE, MeasurementResult.ZERO]
    finally:
        os.unlink(path)


def test_run_stim_file_measurement_order_follows_circuit():
    # Verify that file execution also respects STIM circuit order
    circuit = "X 1\nM 0 1"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".stim", delete=False) as f:
        f.write(circuit)
        path = f.name
    try:
        tab = GeneralizedTableau(2)
        results = tab.run_stim_file(path)
        assert results == [MeasurementResult.ZERO, MeasurementResult.ONE]
    finally:
        os.unlink(path)


def test_run_stim_file_missing_file_raises():
    tab = GeneralizedTableau(1)
    with pytest.raises(BaseException):  # noqa: B017
        tab.run_stim_file("/nonexistent/path/circuit.stim")


def test_run_stim_file_matches_run_stim_string():
    # run_stim_file and run_stim_string must produce identical results
    # for the same deterministic circuit.
    circuit = "X 0\nH 1\nCX 1 2\nM 0 1 2"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".stim", delete=False) as f:
        f.write(circuit)
        path = f.name
    try:
        tab_str = GeneralizedTableau(3, seed=0)
        tab_file = GeneralizedTableau(3, seed=0)
        results_str = tab_str.run_stim_string(circuit)
        results_file = tab_file.run_stim_file(path)
        assert results_str == results_file
    finally:
        os.unlink(path)
