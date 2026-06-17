# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

import ppvm
from ppvm import GeneralizedTableau, MeasurementResult


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
    for name in ["s_dag", "sqrt_x_dag", "sqrt_y_dag", "t_dag",
                 "depolarize1", "x_error", "y_error", "z_error",
                 "reset_x", "reset_y", "reset_z", "cx", "zcx", "zcz", "zcy"]:
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
