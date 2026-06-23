# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

import pytest

from ppvm import GeneralizedTableau, MeasurementResult, PauliSum, StimProgram


def test_single_and_broadcast_gates():
    t = GeneralizedTableau(3)
    t.h(0)
    t.h(1, 2)
    rec = t.measure_many(0, 1, 2)
    assert len(rec) == 3


def test_two_qubit_pair_broadcast():
    t = GeneralizedTableau(4)
    t.h(0, 2)
    t.cnot(0, 1, 2, 3)
    assert t.measure(0) == t.measure(1)


def test_renamed_methods_exist():
    t = GeneralizedTableau(1)
    for name in [
        "s_dag",
        "sqrt_x_dag",
        "sqrt_y_dag",
        "t_dag",
        "depolarize1",
        "x_error",
        "y_error",
        "z_error",
        "reset_x",
        "reset_y",
        "reset_z",
        "cx",
        "zcx",
        "zcz",
        "zcy",
    ]:
        assert hasattr(t, name), name


def test_measurement_record():
    t = GeneralizedTableau(1)
    t.x(0)
    t.measure(0)
    assert t.current_measurement_record() == [MeasurementResult.ONE]


def test_noise_keyword_p():
    t = GeneralizedTableau(2)
    t.x_error(0, 1, p=0.0)
    t.depolarize2(0, 1, p=0.0)


def test_trailing_parameters_accept_positional_or_keyword():
    t = GeneralizedTableau(2)
    t.rx(0, 0.0)
    t.rxx(0, 1, 0.0)
    t.x_error(0, 0.0)
    t.pauli_error(0, [0.0, 0.0, 0.0])
    t.depolarize2(0, 1, 0.0)
    t.two_qubit_pauli_error(0, 1, [0.0] * 15)

    ps = PauliSum.new(2, "ZI")
    ps.rx(0, 0.0)
    ps.rxx(0, 1, 0.0, False)
    ps.x_error(0, 0.0)
    ps.pauli_error(0, [0.0, 0.0, 0.0])
    ps.depolarize2(0, 1, 0.0)
    ps.two_qubit_pauli_error(0, 1, [0.0] * 15)


def test_odd_two_qubit_targets_raise_value_error():
    t = GeneralizedTableau(3)
    with pytest.raises(ValueError, match="even number"):
        t.cnot(0, 1, 2)

    ps = PauliSum.new(3, "ZII")
    with pytest.raises(ValueError, match="even number"):
        ps.rzz(0, 1, 2, 0.0)


# --- targets as a single sequence (list / tuple / range / ndarray), mirroring
# --- stim's ``Circuit.append`` convention, on top of the variadic form.


def _x_flips(apply):
    """Apply X via ``apply`` to a fresh 3-qubit state, return the bit pattern."""
    t = GeneralizedTableau(3)
    apply(t)
    return [t.measure(q) for q in range(3)]


def test_single_qubit_gate_accepts_varargs_and_sequences():
    import numpy as np

    baseline = _x_flips(lambda t: t.x(0, 2))  # variadic ints
    assert baseline == [MeasurementResult.ONE, MeasurementResult.ZERO, MeasurementResult.ONE]
    assert _x_flips(lambda t: t.x([0, 2])) == baseline  # list
    assert _x_flips(lambda t: t.x((0, 2))) == baseline  # tuple
    assert _x_flips(lambda t: t.x(range(0, 3, 2))) == baseline  # range
    assert _x_flips(lambda t: t.x(np.array([0, 2]))) == baseline  # ndarray
    assert _x_flips(lambda t: t.x(0)) == [  # bare int still scalar
        MeasurementResult.ONE,
        MeasurementResult.ZERO,
        MeasurementResult.ZERO,
    ]


def test_numpy_integer_scalar_is_a_single_target():
    import numpy as np

    # A numpy int scalar is not iterable, so it must be treated as one target.
    assert _x_flips(lambda t: t.x(np.int64(2))) == [
        MeasurementResult.ZERO,
        MeasurementResult.ZERO,
        MeasurementResult.ONE,
    ]


def test_two_qubit_gate_accepts_flat_sequence():
    import numpy as np

    for targets in ([0, 1, 2, 3], (0, 1, 2, 3), np.array([0, 1, 2, 3])):
        t = GeneralizedTableau(4)
        t.h([0, 2])
        t.cnot(targets)  # pairs (0,1), (2,3)
        assert t.measure(0) == t.measure(1)
        assert t.measure(2) == t.measure(3)


def test_rotation_and_noise_accept_sequence_with_trailing_param():
    import numpy as np

    t = GeneralizedTableau(4)
    t.rx([0, 1, 2], theta=0.0)  # sequence + kw theta
    t.rx(np.array([0, 1]), 0.0)  # ndarray + positional theta
    t.rxx([0, 1, 2, 3], theta=0.0)  # sequence pairs (0,1),(2,3) + kw theta
    t.x_error([0, 1, 2], p=0.0)  # sequence + kw p
    t.pauli_error(np.array([0, 1]), p=[0.0, 0.0, 0.0])  # ndarray + kw p


def test_measure_many_accepts_sequence():
    import numpy as np

    for targets in ([0, 1, 2], (0, 1, 2), range(3), np.array([0, 1, 2])):
        t = GeneralizedTableau(3)
        t.x([0, 2])
        assert t.measure_many(targets) == [
            MeasurementResult.ONE,
            MeasurementResult.ZERO,
            MeasurementResult.ONE,
        ]


def test_pausisum_truncating_gate_accepts_sequence_and_truncate():
    ps = PauliSum.new(3, "ZII")
    ps.rx([0, 1, 2], 0.0, False)  # sequence + positional theta + truncate
    ps.rx([0, 1], theta=0.0, truncate=False)  # sequence + kw theta + kw truncate
    ps.rxx([0, 1], theta=0.0)  # sequence pairs + kw theta


def test_odd_two_qubit_sequence_raises_value_error():
    t = GeneralizedTableau(3)
    with pytest.raises(ValueError, match="even number"):
        t.cnot([0, 1, 2])


# --- StimProgram pretty-printing / round-trip ---------------------------------


@pytest.mark.parametrize(
    "src",
    [
        "H 0\nCX 0 1\nM 0 1\n",
        "REPEAT 3 {\n    X 0\n    M 0\n}\n",
        "S[T] 0\nI[R_X(theta=0.25)] 1\nI_ERROR[loss](0.01) 2\n",
        "MR 0\nCX rec[-1] 0\n",
        "MPP X0*Y1\nM 2\n",
    ],
)
def test_stim_program_print_is_a_fixpoint(src):
    # str(prog) is canonical: parse -> print -> parse -> print reaches a
    # byte-identical fixpoint, so a parsed program can be serialized via
    # str()/print() and re-parsed losslessly (modulo canonical normalization
    # of comments/whitespace).
    printed = str(StimProgram.parse(src))
    assert str(StimProgram.parse(printed)) == printed


def test_stim_program_print_normalizes_comments_and_whitespace():
    prog = StimProgram.parse("H 0  # flip\nCX  0   1\n")
    assert str(prog) == "H 0\nCX 0 1\n"


def test_stim_program_repr_html_is_highlighted():
    prog = StimProgram.parse("H 0\nM 0\n")
    html = prog._repr_html_()
    assert html.startswith("<pre")
    assert "<span" in html  # tokens are wrapped in coloured spans
    assert "</pre>" in html


def test_parse_error_points_at_the_source():
    with pytest.raises(ValueError) as exc:
        StimProgram.parse("H 0\nCX 0 1\nM 0 X\n")
    msg = str(exc.value)
    assert "M 0 X" in msg  # the offending source line is shown
    assert "invalid target" in msg  # ...with the diagnostic message
    assert "3" in msg  # ...located at line 3
