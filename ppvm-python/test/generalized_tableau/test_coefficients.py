import math

from ppvm import GeneralizedTableau


def test_fresh_state_has_single_unit_coefficient():
    # |0...0> is represented by a single branch with amplitude 1.
    tab = GeneralizedTableau(2)
    assert tab.coefficients() == {0: complex(1.0, 0.0)}


def test_num_coefficients_fresh_state():
    tab = GeneralizedTableau(2)
    assert tab.num_coefficients() == 1


def test_coefficients_returns_plain_dict():
    tab = GeneralizedTableau(2)
    assert isinstance(tab.coefficients(), dict)


def test_clifford_gates_do_not_branch():
    # Clifford gates act on the tableau frame, leaving the coefficient
    # vector a single branch (the Bell state stays one branch).
    tab = GeneralizedTableau(2)
    tab.h(0)
    tab.cnot(0, 1)
    assert tab.num_coefficients() == 1


def test_non_clifford_gate_branches_coefficient_vector():
    # T on |+> is a magic state -> needs more than one stabilizer branch.
    tab = GeneralizedTableau(2)
    tab.h(0)
    tab.t(0)
    coeffs = tab.coefficients()
    assert tab.num_coefficients() > 1
    assert tab.num_coefficients() == len(coeffs)


def test_coefficient_vector_stays_normalized():
    # Unitary evolution preserves the L2 norm of the coefficient vector.
    tab = GeneralizedTableau(3)
    tab.h(0)
    tab.t(0)
    tab.rx(1, theta=math.pi / 3)
    tab.cnot(0, 2)
    norm_sq = sum(abs(c) ** 2 for c in tab.coefficients().values())
    assert math.isclose(norm_sq, 1.0, abs_tol=1e-9)


def test_coefficients_is_a_snapshot_not_a_live_view():
    # Mutating the returned dict must not change the tableau's state.
    tab = GeneralizedTableau(2)
    snapshot = tab.coefficients()
    snapshot[0] = complex(42.0, 42.0)
    snapshot[999] = complex(1.0, 0.0)
    assert tab.coefficients() == {0: complex(1.0, 0.0)}


def test_coefficients_lossless_for_wide_index_types():
    # Past 128 qubits the Rust index is a bnum type; the keys must still
    # round-trip as exact Python ints and the vector stay normalized.
    tab = GeneralizedTableau(200)
    tab.h(0)
    tab.t(0)
    coeffs = tab.coefficients()
    assert all(isinstance(k, int) for k in coeffs)
    assert all(isinstance(v, complex) for v in coeffs.values())
    assert tab.num_coefficients() == len(coeffs)
    norm_sq = sum(abs(c) ** 2 for c in coeffs.values())
    assert math.isclose(norm_sq, 1.0, abs_tol=1e-9)
