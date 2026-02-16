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
