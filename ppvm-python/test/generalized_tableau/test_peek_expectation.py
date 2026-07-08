import math

import pytest

from ppvm import ExpectationResult, GeneralizedTableau


def bell() -> GeneralizedTableau:
    """The Stim docstring's Bell pair: H 0; CNOT 0 1 (forward/Schrödinger)."""
    tab = GeneralizedTableau(3, seed=0)
    tab.h(0)
    tab.cnot(0, 1)
    return tab


def peek(tab: GeneralizedTableau, obs: str) -> float:
    result = tab.peek_observable_expectation(obs)
    assert not result.is_lost
    return float(result)


def test_matches_stim_docstring():
    tab = bell()
    assert peek(tab, "X0*X1") == pytest.approx(1.0)
    assert peek(tab, "Y0*Y1") == pytest.approx(-1.0)
    assert peek(tab, "Z0*Z1") == pytest.approx(1.0)
    assert peek(tab, "-Z0*Z1") == pytest.approx(-1.0)
    assert peek(tab, "Z0") == pytest.approx(0.0)
    assert peek(tab, "Z2") == pytest.approx(1.0)


def test_dense_and_sparse_agree():
    tab = bell()
    assert peek(tab, "ZZI") == pytest.approx(peek(tab, "Z0*Z1"))
    assert peek(tab, "XXI") == pytest.approx(peek(tab, "X0*X1"))


def test_identity_is_plus_one():
    tab = bell()
    assert peek(tab, "III") == pytest.approx(1.0)
    assert peek(tab, "-III") == pytest.approx(-1.0)


def test_non_clifford_continuous_value():
    tab = GeneralizedTableau(1, seed=0)
    tab.ry(0, 0.7)
    assert peek(tab, "Z0") == pytest.approx(math.cos(0.7))


def test_lost_support_qubit_returns_lost():
    tab = GeneralizedTableau(2, seed=0)
    tab.h(0)
    tab.cnot(0, 1)
    tab.loss_channel(1, 1.0)
    assert tab.is_lost(1)

    result = tab.peek_observable_expectation("Z0*Z1")
    assert result.is_lost
    assert result is ExpectationResult.LOST or result.value is None
    with pytest.raises(ValueError):
        float(result)

    # An observable that avoids the lost qubit still returns a value.
    assert not tab.peek_observable_expectation("Z0").is_lost


def test_peek_does_not_disturb_state():
    tab = bell()
    before = tab.coefficients()
    tab.peek_observable_expectation("X0*X1")
    tab.peek_observable_expectation("Z0*Z1")
    assert tab.coefficients() == before


def test_malformed_observable_raises():
    tab = GeneralizedTableau(3, seed=0)
    with pytest.raises(ValueError):
        tab.peek_observable_expectation("Z5")  # out of range
    with pytest.raises(ValueError):
        tab.peek_observable_expectation("Z0*Z0")  # repeated qubit
    with pytest.raises(ValueError):
        tab.peek_observable_expectation("Q0")  # bad token


def test_expectation_result_value_object():
    r = ExpectationResult(0.5)
    assert r.value == 0.5
    assert not r.is_lost
    assert float(r) == 0.5
    assert ExpectationResult.LOST.is_lost
    assert ExpectationResult.LOST.value is None


def test_large_qubit_count_uses_wide_index():
    # n=70 selects the u128-indexed native interface (GeneralizedTableau2),
    # exercising the decomposition's bit-shifts beyond 64 bits.
    n = 70
    tab = GeneralizedTableau(n, seed=0)
    tab.h(0)
    tab.cnot(0, 69)  # Bell pair on the extreme qubits
    assert peek(tab, "Z0*Z69") == pytest.approx(1.0)
    assert peek(tab, "X0*X69") == pytest.approx(1.0)
    assert peek(tab, "Z0") == pytest.approx(0.0)
    assert peek(tab, f"Z{n - 1}") == pytest.approx(0.0)


def test_cross_validate_against_stim():
    """For Clifford circuits, ppvm's peek must agree with Stim's exactly.

    Stim's TableauSimulator is the Schrödinger-picture reference (gates forward,
    same as GeneralizedTableau) and returns the exact {-1, 0, +1} eigenvalue.
    """
    import random

    stim = pytest.importorskip("stim")

    rng = random.Random(20260626)
    n = 5
    single = ["h", "s", "x", "y", "z", "sqrt_x"]
    two = ["cnot", "cz"]

    for _ in range(25):
        tab = GeneralizedTableau(n, seed=0)
        sim = stim.TableauSimulator()
        for _ in range(40):
            if rng.random() < 0.3:
                g = rng.choice(two)
                a = rng.randrange(n)
                b = rng.randrange(n)
                while b == a:
                    b = rng.randrange(n)
                getattr(tab, g)(a, b)
                getattr(sim, g)(a, b)
            else:
                g = rng.choice(single)
                q = rng.randrange(n)
                getattr(tab, g)(q)
                getattr(sim, g)(q)

        for _ in range(12):
            dense = "".join(rng.choice("IXYZ") for _ in range(n))
            sign = rng.choice(["+", "-"])
            expected = sim.peek_observable_expectation(stim.PauliString(sign + dense))
            got = float(tab.peek_observable_expectation(sign + dense))
            assert got == pytest.approx(expected, abs=1e-9), (
                f"observable {sign + dense!r}: ppvm {got} vs stim {expected}"
            )
