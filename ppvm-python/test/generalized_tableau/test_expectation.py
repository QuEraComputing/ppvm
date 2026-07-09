import math

import pytest

from ppvm import GeneralizedTableau


def _bell() -> GeneralizedTableau:
    tab = GeneralizedTableau(2)
    tab.h(0)
    tab.cnot(0, 1)
    return tab


def test_expectation_z_on_zero_state():
    tab = GeneralizedTableau(1)
    assert tab.expectation("Z") == pytest.approx(1.0, abs=1e-10)
    assert tab.expectation("X") == pytest.approx(0.0, abs=1e-10)
    assert tab.expectation("I") == pytest.approx(1.0, abs=1e-10)


def test_expectation_x_on_plus_state():
    tab = GeneralizedTableau(1)
    tab.h(0)
    assert tab.expectation("X") == pytest.approx(1.0, abs=1e-10)
    assert tab.expectation("Z") == pytest.approx(0.0, abs=1e-10)


def test_bell_state_pauli_expectations():
    tab = _bell()
    assert tab.expectation("II") == pytest.approx(1.0, abs=1e-10)
    assert tab.expectation("ZZ") == pytest.approx(1.0, abs=1e-10)
    assert tab.expectation("XX") == pytest.approx(1.0, abs=1e-10)
    assert tab.expectation("YY") == pytest.approx(-1.0, abs=1e-10)
    assert tab.expectation("IZ") == pytest.approx(0.0, abs=1e-10)
    assert tab.expectation("ZI") == pytest.approx(0.0, abs=1e-10)


def test_ry_rotation_z_expectation_is_cos_theta():
    tab = GeneralizedTableau(1)
    theta = math.pi / 3
    tab.ry(0, theta)
    assert tab.expectation("Z") == pytest.approx(math.cos(theta), abs=1e-10)
    assert tab.expectation("X") == pytest.approx(math.sin(theta), abs=1e-10)


def test_trace_z_or_identity_pattern_on_bell():
    tab = _bell()
    # ⟨II⟩ + ⟨IZ⟩ + ⟨ZI⟩ + ⟨ZZ⟩ = 1 + 0 + 0 + 1 = 2
    assert tab.trace("Z?{2}") == pytest.approx(2.0, abs=1e-10)


def test_trace_y_or_identity_pattern_on_bell_is_zero():
    tab = _bell()
    # ⟨II⟩ + ⟨IY⟩ + ⟨YI⟩ + ⟨YY⟩ = 1 + 0 + 0 + (-1) = 0
    assert tab.trace("Y?{2}") == pytest.approx(0.0, abs=1e-10)


def test_trace_positional_pattern_matches_single_pauli():
    tab = _bell()
    assert tab.trace("Z0Z1") == pytest.approx(tab.expectation("ZZ"), abs=1e-10)
