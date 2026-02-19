import pytest

from ppvm import PauliSum


def test_basics():
    state = PauliSum(initial_terms=["ZZ"], coefficients=[1.0])  # ZZ

    state.cnot(0, 1)
    state.h(0)

    assert str(state) == "1.000 * IZ"
    assert state.overlap_with_zero() == 1.0


def test_noise():
    state = PauliSum(initial_terms=["IZ"], coefficients=[1.0])  # |00><00|

    error_probs = {"ZZ": 0.1, "XX": 0.2}
    error_probs_list = state.two_qubit_pauli_error_probabilities(error_probs)
    state.two_qubit_pauli_error(0, 1, error_probs_list)


def test_large_state():
    n = 200
    weight = 80

    terms = ["".join(["Z" if i == j else "I" for i in range(n)]) for j in range(n)]
    large_state = PauliSum(max_pauli_weight=weight, initial_terms=terms)

    for i in reversed(range(1, n)):
        large_state.cnot(i - 1, i)

    large_state.h(0)

    assert large_state.overlap_with_zero() == 0.0


def test_copy():

    state = PauliSum(initial_terms=["ZZ"], coefficients=[1.0])  # ZZ

    state.cnot(0, 1)
    state.h(0)

    assert str(state) == "1.000 * IZ"
    assert state.overlap_with_zero() == 1.0

    tmp = state.copy()
    assert tmp == state

    assert len(tmp) == len(state) == 1

    assert tmp.terms == [("IZ", 1.0)]


def test_weights():
    state = PauliSum(initial_terms=["IZ"])

    assert state.current_max_weight() == 1

    state.cnot(0, 1)
    state.h(0)

    assert state.current_max_weight() == 2

    state2 = PauliSum(initial_terms=["ZX", "IY"])
    weights = state2.weights()

    weights.sort(key=lambda w: w[1])
    assert weights == [('IY', 1), ('ZX', 2)]


def test_overlap():
    state = PauliSum(initial_terms=["IZ"])

    assert state.overlap(state) == 1.0

    state2 = PauliSum(initial_terms=["IX"])
    assert state.overlap(state2) == 0.0


def test_new():
    s = PauliSum.new(2, "IX")
    assert len(s) == 1
    assert s.terms == [("IX", 1.0)]

    s = PauliSum.new(2, ("IX", 0.25))
    assert len(s) == 1
    assert s.terms == [("IX", 0.25)]

    s = PauliSum.new(2, "X1")
    assert s == PauliSum.new(2, "IX")

    s = PauliSum.new(3, [("Y1", 0.1), "ZIZ"])
    assert len(s) == 2
    assert set(s.terms) == {("IYI", 0.1), ("ZIZ", 1.0)}

    n = 17
    s = PauliSum.new(n, [f"Z{i}" for i in range(n)])
    assert len(s) == n
    assert set(s.terms) == {("".join(["Z" if i==j else "I" for i in range(n)]), 1.0) for j in range(n)}

    with pytest.raises(ValueError, match="out of range"):
        PauliSum.new(2, "X2")

    n = 5
    terms = ["Z0", "Z1Z2"]
    s = PauliSum.new(n, terms)
