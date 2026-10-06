import itertools
import math

import pytest

from ppvm import GeneralizedTableau
from ppvm.generalized_tableau import MeasurementResult

ZERO, ONE, LOST = MeasurementResult.ZERO, MeasurementResult.ONE, MeasurementResult.LOST


def _magic_circuit() -> GeneralizedTableau:
    tab = GeneralizedTableau(3, min_abs_coeff=1e-12)
    tab.h(0)
    tab.t(0)
    tab.cnot(0, 1)
    tab.h(2)
    tab.t(2)
    tab.cnot(1, 2)
    tab.ry(1, 0.7)
    tab.t(1)
    tab.h(1)
    tab.rx(2, 0.3)
    return tab


def _prob_by_projection(tab: GeneralizedTableau, bits: tuple[int, ...]) -> float:
    t = tab.fork(seed=0)
    prob = 1.0
    for q, b in enumerate(bits):
        try:
            prob *= t.project(q, MeasurementResult(b))
        except ValueError:
            return 0.0
    return prob


def _prob_by_expectation(tab: GeneralizedTableau, bits: tuple[int, ...]) -> float:
    # P(b) = 2^-n Σ_T (-1)^{b·T} ⟨Z_T⟩
    n = len(bits)
    total = 0.0
    for mask in itertools.product((0, 1), repeat=n):
        word = "".join("Z" if m else "I" for m in mask)
        sign = (-1) ** sum(m * b for m, b in zip(mask, bits))
        total += sign * tab.expectation(word)
    return total / 2**n


def test_project_zero_state():
    tab = GeneralizedTableau(1)
    assert tab.project(0, ZERO) == pytest.approx(1.0, abs=1e-12)
    assert tab.current_measurement_record() == [ZERO]


def _assert_zero_probability_preserves_state(tab: GeneralizedTableau, addr0: int, value):
    paulis = ["".join(p) for p in itertools.product("IXYZ", repeat=tab.n_qubits)]
    before_str = str(tab)
    before = [tab.expectation(p) for p in paulis]
    before_coeffs = tab.coefficients()
    before_record = tab.current_measurement_record()
    with pytest.raises(ValueError, match="probability"):
        tab.project(addr0, value)
    assert str(tab) == before_str
    assert tab.coefficients() == before_coeffs
    assert [tab.expectation(p) for p in paulis] == pytest.approx(before, abs=1e-12)
    assert tab.current_measurement_record() == before_record


def test_project_zero_probability_z_stabilizer():
    # |0⟩: Z is a stabilizer (case b).
    tab = GeneralizedTableau(1)
    _assert_zero_probability_preserves_state(tab, 0, ONE)

    # Projecting |+⟩·T onto 0 makes Z a stabilizer; projecting onto 1 is then impossible.
    tab = GeneralizedTableau(1)
    tab.h(0)
    tab.t(0)
    tab.project(0, ZERO)
    _assert_zero_probability_preserves_state(tab, 0, ONE)


def test_project_zero_probability_z_not_stabilizer():
    # H then RY(-π/2) is |0⟩, but the stabilizer frame still holds X, so the
    # projection takes the case-a path with P(1) = 0.
    tab = GeneralizedTableau(1)
    tab.h(0)
    tab.ry(0, -math.pi / 2)
    _assert_zero_probability_preserves_state(tab, 0, ONE)
    assert tab.project(0, ZERO) == pytest.approx(1.0, abs=1e-12)

    # 3 qubits: qubit 0 as above, next to an entangled non-Clifford pair.
    tab = GeneralizedTableau(3, min_abs_coeff=1e-12)
    tab.h(0, 1, 2)
    tab.t(1)
    tab.cz(1, 2)
    tab.ry(0, -math.pi / 2)
    assert tab.num_coefficients() > 1
    _assert_zero_probability_preserves_state(tab, 0, ONE)
    assert tab.project(0, ZERO) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("value", [ZERO, ONE])
def test_project_plus_state(value):
    tab = GeneralizedTableau(1)
    tab.h(0)
    assert tab.project(0, value) == pytest.approx(0.5, abs=1e-12)
    expected_z = 1.0 if value == ZERO else -1.0
    assert tab.expectation("Z") == pytest.approx(expected_z, abs=1e-12)
    assert tab.project(0, value) == pytest.approx(1.0, abs=1e-12)
    assert tab.current_measurement_record() == [value, value]


@pytest.mark.parametrize("theta", [0.1, 0.7, 1.3, 2.0, 2.9])
def test_project_ry_matches_cos_squared(theta):
    p0 = math.cos(theta / 2) ** 2
    for value, expected in [(ZERO, p0), (ONE, 1 - p0)]:
        tab = GeneralizedTableau(1)
        tab.ry(0, theta)
        assert tab.project(0, value) == pytest.approx(expected, abs=1e-10)


def test_project_matches_measure():
    base = _magic_circuit()
    paulis = ["".join(p) for p in itertools.product("IXYZ", repeat=3)]
    for q in range(3):
        for value in (ZERO, ONE):
            projected = base.fork(seed=0)
            try:
                projected.project(q, value)
            except ValueError:
                continue
            measured = next(t for s in range(1000) if (t := base.fork(seed=s)).measure(q) == value)
            for p in paulis:
                assert projected.expectation(p) == pytest.approx(measured.expectation(p), abs=1e-9)


def test_project_chain_rule_gives_bitstring_probabilities():
    tab = _magic_circuit()
    total = 0.0
    for bits in itertools.product((0, 1), repeat=3):
        p = _prob_by_projection(tab, bits)
        assert p == pytest.approx(_prob_by_expectation(tab, bits), abs=1e-9)
        total += p
    assert total == pytest.approx(1.0, abs=1e-9)


def test_project_ghz():
    tab = GeneralizedTableau(3)
    tab.h(0)
    tab.cnot(0, 1)
    tab.cnot(1, 2)
    for bits in itertools.product((0, 1), repeat=3):
        expected = 0.5 if len(set(bits)) == 1 else 0.0
        assert _prob_by_projection(tab, bits) == pytest.approx(expected, abs=1e-12)


def test_project_lost_value_not_implemented():
    tab = GeneralizedTableau(1)
    with pytest.raises(NotImplementedError):
        tab.project(0, LOST)


def test_project_lost_qubit_not_implemented():
    tab = GeneralizedTableau(2)
    tab.loss_channel(0, 1.0)
    with pytest.raises(NotImplementedError, match="lost"):
        tab.project(0, ZERO)
    assert tab.current_measurement_record() == []
    assert tab.project(1, ZERO) == pytest.approx(1.0, abs=1e-12)
