# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

from copy import copy

import pytest

from ppvm import GeneralizedTableau, LossyPauliSum, MeasurementResult, PauliSum
from ppvm.generalized_tableau_sum import GeneralizedTableauSum


def ccz_with_t_gates(tab, a, b, c):
    """Reference: 4abc = a+b+c - (a^b) - (a^c) - (b^c) + (a^b^c)."""
    tab.t(a, b, c)
    tab.cnot(a, b)
    tab.t_dag(b)
    tab.cnot(a, b)
    tab.cnot(a, c)
    tab.t_dag(c)
    tab.cnot(b, c)
    tab.t(c)
    tab.cnot(a, c)
    tab.t_dag(c)
    tab.cnot(b, c)


@pytest.mark.parametrize("backend", [GeneralizedTableau, GeneralizedTableauSum])
def test_ccz_batch_matches_t_gates(backend):
    tab = backend(5)
    tab.h(range(5))
    tab.s(1, 3)
    tab.ry(2, theta=0.37)
    tab.cnot(0, 4)

    # Native CCZ followed by its T-gate inverse must preserve the input state.
    tab.ccz([4, 1, 3, 2, 4, 0])
    ccz_with_t_gates(tab, 2, 4, 0)
    ccz_with_t_gates(tab, 4, 1, 3)

    tab.cnot(0, 4)
    tab.ry(2, theta=-0.37)
    tab.s_dag(1, 3)
    tab.h(range(5))
    for q in range(5):
        if isinstance(tab, GeneralizedTableauSum):
            assert tab.measure(q)[MeasurementResult.ZERO] == pytest.approx(1.0)
        else:
            assert tab.expectation("I" * q + "Z" + "I" * (4 - q)) == pytest.approx(1.0)


@pytest.mark.parametrize("backend", [PauliSum, LossyPauliSum])
def test_ccz_can_defer_truncation(backend):
    observable = backend.new(3, "XII", min_abs_coeff=0.6)
    truncated = copy(observable)
    truncated.ccz(0, 1, 2)
    assert len(truncated) == 0

    observable.ccz([0, 1, 2], truncate=False)
    observable.ccz(0, 1, 2)
    assert observable.terms == [("XII", 1.0)]


@pytest.mark.parametrize("inactive", range(3))
@pytest.mark.parametrize("channel", ["loss", "leakage"])
def test_ccz_skips_inactive_operands(inactive, channel):
    tab = GeneralizedTableau(3)
    if channel == "loss":
        tab.loss_channel(inactive, 1.0)
    else:
        tab.leakage_channel(inactive, 0.0, 1.0)
    active = [q for q in range(3) if q != inactive]
    tab.h(active)
    tab.ccz(0, 1, 2)
    tab.h(active)
    assert tab.measure_many(active) == [MeasurementResult.ZERO] * 2


@pytest.mark.parametrize(
    "backend", [PauliSum, LossyPauliSum, GeneralizedTableau, GeneralizedTableauSum]
)
@pytest.mark.parametrize("targets", [[0, 1], [0, 1, 1], [0, 1, 3]])
def test_ccz_rejects_invalid_targets(backend, targets):
    tab = backend.new(3, "XII") if issubclass(backend, PauliSum) else backend(3)
    with pytest.raises(ValueError):
        tab.ccz(targets)
