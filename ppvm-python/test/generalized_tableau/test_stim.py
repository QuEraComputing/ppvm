import os
import tempfile
import textwrap

import pytest

from ppvm import GeneralizedTableau, StimProgram, sample_stim
from ppvm.generalized_tableau import MeasurementResult

# --- run (via StimProgram.parse) ---


def test_run_stim_string_single_measurement_zero():
    # Fresh |0> state → M returns ZERO
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("M 0"))
    assert results == [MeasurementResult.ZERO]


def test_run_stim_string_x_then_measure():
    # X|0> = |1> → M returns ONE
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("X 0\nM 0"))
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_measurement_order_follows_circuit():
    # M 0 1: qubit 0 was flipped, qubit 1 was not
    # Results must be [ONE, ZERO] in that order (circuit order, STIM convention)
    tab = GeneralizedTableau(2)
    results = tab.run(StimProgram.parse("X 0\nM 0 1"))
    assert results == [MeasurementResult.ONE, MeasurementResult.ZERO]


def test_run_stim_string_multiple_m_instructions_appended_in_order():
    # Two separate M instructions: measurements from the first instruction
    # appear before those from the second, following STIM circuit order.
    tab = GeneralizedTableau(2)
    results = tab.run(StimProgram.parse("X 0\nM 1\nM 0"))
    assert results == [MeasurementResult.ZERO, MeasurementResult.ONE]


def test_run_stim_string_double_measurement_same_qubit():
    # Measuring the same qubit twice appends two entries
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("X 0\nM 0\nM 0"))
    assert results == [MeasurementResult.ONE, MeasurementResult.ONE]


def test_run_stim_string_no_measurements_returns_empty():
    tab = GeneralizedTableau(2)
    results = tab.run(StimProgram.parse("X 0\nX 1"))
    assert results == []


def test_run_stim_string_reset():
    # R resets to |0>, so subsequent measurement is ZERO
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("X 0\nR 0\nM 0"))
    assert results == [MeasurementResult.ZERO]


def test_run_stim_string_mr_resets_after_measure():
    # MR measures and then resets; second M should always be ZERO
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("X 0\nMR 0\nM 0"))
    assert results == [MeasurementResult.ONE, MeasurementResult.ZERO]


def test_run_stim_string_cnot_entanglement():
    # H(0) + CX(0,1) creates a Bell state; both qubits must agree
    tab = GeneralizedTableau(2)
    results = tab.run(StimProgram.parse("H 0\nCX 0 1\nM 0 1"))
    assert len(results) == 2
    assert results[0] == results[1]


def test_run_stim_string_cz_gate():
    tab = GeneralizedTableau(2)
    results = tab.run(StimProgram.parse("X 0\nX 1\nCZ 0 1\nM 0 1"))
    assert results == [MeasurementResult.ONE, MeasurementResult.ONE]


def test_run_stim_string_h_xz_alias():
    # H_XZ is an alias for H
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("X 0\nH_XZ 0\nH_XZ 0\nM 0"))
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_s_gate():
    # S|0> = |0>
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("S 0\nM 0"))
    assert results == [MeasurementResult.ZERO]


def test_run_stim_string_t_gate_tag():
    # S[T] is the T gate; T|0> = |0>
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("S[T] 0\nM 0"))
    assert results == [MeasurementResult.ZERO]


def test_run_stim_string_s_dag_t_tag():
    # S_DAG[T] is the T† gate; T T† |1> = |1>
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("X 0\nS[T] 0\nS_DAG[T] 0\nM 0"))
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_sqrt_x():
    # SQRT_X SQRT_X = X
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("SQRT_X 0\nSQRT_X 0\nM 0"))
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_sqrt_y():
    # SQRT_Y SQRT_Y = Y (up to global phase); Y|0> = i|1> → measures as 1
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("SQRT_Y 0\nSQRT_Y 0\nM 0"))
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_i_rotation_x_pi():
    # I[R_X(theta=1.0*pi)]|0> = |1>
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("I[R_X(theta=1.0*pi)] 0\nM 0"))
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_i_rotation_y_pi():
    # I[R_Y(theta=1.0*pi)]|0> = |1>
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("I[R_Y(theta=1.0*pi)] 0\nM 0"))
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_i_rotation_z_pi():
    # I[R_Z(theta=1.0*pi)]|0> = |0> (Z rotation only adds phase)
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("I[R_Z(theta=1.0*pi)] 0\nM 0"))
    assert results == [MeasurementResult.ZERO]


def test_run_stim_string_loss_channel():
    # I_ERROR[loss](1.0) loses the qubit with probability 1
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("I_ERROR[loss](1.0) 0\nM 0"))
    assert results == [MeasurementResult.LOST]


def test_run_stim_string_comments_and_blank_lines_ignored():
    # Comments (#) and blank lines must not affect execution
    circuit = """
# Prepare and measure
X 0

M 0
"""
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse(circuit))
    assert results == [MeasurementResult.ONE]


def test_run_stim_string_tick_and_annotations_are_noops():
    # TICK, DETECTOR, OBSERVABLE_INCLUDE, QUBIT_COORDS, SHIFT_COORDS,
    # are annotation-only and must not affect circuit outcomes.
    # Note: rec[-1] targets on annotations are silently dropped during
    # parsing; the annotation itself is preserved as a no-op.
    circuit = (
        "QUBIT_COORDS(0, 0) 0\nX 0\nTICK\nM 0\nDETECTOR rec[-1]\nOBSERVABLE_INCLUDE(0) rec[-1]\n"
    )
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse(circuit))
    assert results == [MeasurementResult.ONE]


# --- run (via StimProgram.from_file) ---


def test_run_stim_file_basic():
    circuit = "X 0\nM 0 1"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".stim", delete=False) as f:
        f.write(circuit)
        path = f.name
    try:
        tab = GeneralizedTableau(2)
        results = tab.run(StimProgram.from_file(path))
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
        results = tab.run(StimProgram.from_file(path))
        assert results == [MeasurementResult.ZERO, MeasurementResult.ONE]
    finally:
        os.unlink(path)


def test_run_stim_file_missing_file_raises():
    with pytest.raises(OSError):
        StimProgram.from_file("/nonexistent/path/circuit.stim")


def test_run_stim_file_matches_run_stim_string():
    # from_file and parse must produce identical results
    # for the same deterministic circuit.
    circuit = "X 0\nH 1\nCX 1 2\nM 0 1 2"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".stim", delete=False) as f:
        f.write(circuit)
        path = f.name
    try:
        tab_str = GeneralizedTableau(3, seed=0)
        tab_file = GeneralizedTableau(3, seed=0)
        results_str = tab_str.run(StimProgram.parse(circuit))
        results_file = tab_file.run(StimProgram.from_file(path))
        assert results_str == results_file
    finally:
        os.unlink(path)


# --- sample_stim and GeneralizedTableau.sample ---


def test_sample_stim_returns_list_of_lists():
    prog = StimProgram.parse("X 0\nM 0")
    shots = sample_stim(prog, n_qubits=1, num_shots=3, seed=0)
    assert shots == [
        [MeasurementResult.ONE],
        [MeasurementResult.ONE],
        [MeasurementResult.ONE],
    ]


def test_generalized_tableau_sample_classmethod_equivalent():
    prog = StimProgram.parse("X 0\nM 0")
    a = GeneralizedTableau.sample(prog, 1, num_shots=3, seed=0)
    b = sample_stim(prog, n_qubits=1, num_shots=3, seed=0)
    assert a == b


def test_sample_stim_zero_shots_returns_empty():
    prog = StimProgram.parse("X 0\nM 0")
    assert sample_stim(prog, n_qubits=1, num_shots=0) == []


def test_sample_stim_seeded_is_reproducible_for_large_batches():
    # A randomising circuit at a large shot count. Per-shot seeds are derived
    # from the shot index, so two runs must agree exactly regardless of whether
    # the batch runs serially or in parallel (the serial/parallel cutoff depends
    # on rayon's pool size, which we don't control here) and independent of how
    # rayon schedules the shots.
    prog = StimProgram.parse("H 0\nM 0")
    a = sample_stim(prog, n_qubits=1, num_shots=512, seed=7)
    b = sample_stim(prog, n_qubits=1, num_shots=512, seed=7)
    assert a == b
    assert len(a) == 512
    # Not degenerate: an H gate should produce a mix of outcomes.
    flat = {shot[0] for shot in a}
    assert flat == {MeasurementResult.ZERO, MeasurementResult.ONE}


def test_run_propagates_parse_error_as_value_error():
    with pytest.raises(ValueError):
        StimProgram.parse("FROBNICATE 0")


def test_sample_many_qubits():
    stim_str = textwrap.dedent("""
        X 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99
        M 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99
        """)
    prog = StimProgram.parse(stim_str)
    result = sample_stim(prog, n_qubits=100, num_shots=1, seed=0)
    assert result == [[MeasurementResult.ONE] * 100]


# --- X/Y-basis measurement & reset ---


def test_run_stim_rx_then_mx_reads_zero():
    # RX prepares |+>, the +1 eigenstate of X.
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("RX 0\nMX 0"))
    assert results == [MeasurementResult.ZERO]


def test_run_stim_z_flips_x_basis_outcome():
    # Z|+> = |->, the -1 eigenstate of X.
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("RX 0\nZ 0\nMX 0"))
    assert results == [MeasurementResult.ONE]


def test_run_stim_ry_then_my_reads_zero():
    # RY prepares |i>, the +1 eigenstate of Y.
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("RY 0\nMY 0"))
    assert results == [MeasurementResult.ZERO]


def test_run_stim_mrx_resets_to_plus():
    # MRX records the outcome and resets to |+>, so the next MX reads ZERO.
    tab = GeneralizedTableau(1)
    results = tab.run(StimProgram.parse("RX 0\nZ 0\nMRX 0\nMX 0"))
    assert results == [MeasurementResult.ONE, MeasurementResult.ZERO]


# --- measurement-record controlled feed-forward ---


def test_run_stim_cx_rec_applies_when_bit_set():
    # q0 measured as 1 -> CX rec[-1] 1 applies X to q1.
    tab = GeneralizedTableau(2)
    results = tab.run(StimProgram.parse("X 0\nM 0\nCX rec[-1] 1\nM 1"))
    assert results == [MeasurementResult.ONE, MeasurementResult.ONE]


def test_run_stim_cx_rec_noop_when_bit_clear():
    # q0 measured as 0 -> CX rec[-1] 1 does nothing.
    tab = GeneralizedTableau(2)
    results = tab.run(StimProgram.parse("M 0\nCX rec[-1] 1\nM 1"))
    assert results == [MeasurementResult.ZERO, MeasurementResult.ZERO]


def test_run_stim_record_target_rejected():
    # The Pauli target may never be a measurement record.
    with pytest.raises(ValueError):
        StimProgram.parse("M 0\nCX 1 rec[-1]")


# --- MPP multi-qubit Pauli-product measurement ---


def test_run_stim_mpp_zz_measures_parity():
    # |01> has odd Z0*Z1 parity -> -1 -> ONE.
    tab = GeneralizedTableau(2)
    results = tab.run(StimProgram.parse("X 0\nMPP Z0*Z1"))
    assert results == [MeasurementResult.ONE]


def test_run_stim_mpp_xx_on_bell_state():
    # Bell state is a +1 eigenstate of X0*X1 -> ZERO.
    tab = GeneralizedTableau(2)
    results = tab.run(StimProgram.parse("H 0\nCX 0 1\nMPP X0*X1"))
    assert results == [MeasurementResult.ZERO]


def test_run_stim_mpp_multiple_products():
    # Two space-separated products -> two results.
    tab = GeneralizedTableau(2)
    results = tab.run(StimProgram.parse("X 0\nMPP Z0 Z1"))
    assert results == [MeasurementResult.ONE, MeasurementResult.ZERO]
