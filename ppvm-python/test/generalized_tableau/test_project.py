import math

import pytest

from ppvm import GeneralizedTableau
from ppvm.generalized_tableau import MeasurementResult

ZERO, ONE, LOST = MeasurementResult.ZERO, MeasurementResult.ONE, MeasurementResult.LOST


def test_project_returns_probability_and_records():
    theta = 1.1
    tab = GeneralizedTableau(1)
    tab.ry(0, theta=theta)
    assert tab.project(0, ONE) == pytest.approx(math.sin(theta / 2) ** 2, abs=1e-10)
    assert tab.expectation("Z") == pytest.approx(-1.0, abs=1e-10)
    assert tab.current_measurement_record() == [ONE]


def test_impossible_outcome_raises():
    tab = GeneralizedTableau(1)
    with pytest.raises(ValueError, match="probability"):
        tab.project(0, ONE)


def test_tol_keyword():
    tab = GeneralizedTableau(1)
    tab.ry(0, theta=2 * math.asin(1e-4))  # P(1) = 1e-8
    with pytest.raises(ValueError, match="probability"):
        tab.project(0, ONE, tol=1e-6)
    assert tab.project(0, ONE) == pytest.approx(1e-8, rel=1e-6)


@pytest.mark.parametrize("tol", [0.0, -1.0, math.inf, math.nan])
def test_invalid_tol_raises(tol):
    with pytest.raises(ValueError, match="tol"):
        GeneralizedTableau(1).project(0, ZERO, tol=tol)


def test_lost_raises():
    tab = GeneralizedTableau(1)
    with pytest.raises(ValueError, match="ZERO or ONE"):
        tab.project(0, LOST)
    tab.loss_channel(0, 1.0)
    with pytest.raises(ValueError, match="lost"):
        tab.project(0, ZERO)
