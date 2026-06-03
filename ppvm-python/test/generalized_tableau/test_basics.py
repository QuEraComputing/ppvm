import math

from ppvm import GeneralizedTableau


def test_construction():
    tab = GeneralizedTableau(4)
    assert tab.n_qubits == 4
    assert str(tab) is not None


def test_measure_zero_state():
    # Fresh state is |0...0> — all measurements deterministically False
    tab = GeneralizedTableau(4)
    for i in range(4):
        assert not tab.measure(i)


def test_x_gate():
    tab = GeneralizedTableau(2)
    tab.x(0)
    assert tab.measure(0)
    assert not tab.measure(1)


def test_x_twice_is_identity():
    tab = GeneralizedTableau(2)
    tab.x(0)
    tab.x(0)
    assert not tab.measure(0)


def test_z_on_zero_state():
    # Z|0> = |0>, Z-measurement unchanged
    tab = GeneralizedTableau(2)
    tab.z(0)
    assert not tab.measure(0)


def test_h_twice_is_identity():
    # H²|1> = |1>
    tab = GeneralizedTableau(2)
    tab.x(0)
    tab.h(0)
    tab.h(0)
    assert tab.measure(0)


def test_s_on_zero_state():
    # S|0> = |0>, Z-measurement unchanged
    tab = GeneralizedTableau(2)
    tab.s(0)
    assert not tab.measure(0)


def test_cnot_on_zero_state():
    # CNOT|00> = |00>
    tab = GeneralizedTableau(2)
    tab.cnot(0, 1)
    assert not tab.measure(0)
    assert not tab.measure(1)


def test_cnot_flips_target():
    # CNOT|10> = |11>
    tab = GeneralizedTableau(2)
    tab.x(0)
    tab.cnot(0, 1)
    assert tab.measure(0)
    assert tab.measure(1)


def test_cz_on_zero_state():
    tab = GeneralizedTableau(2)
    tab.cz(0, 1)
    assert not tab.measure(0)
    assert not tab.measure(1)


def test_bell_state_correlated_measurements():
    # H(0), CNOT(0,1) → (|00> + |11>) / sqrt(2)
    # Measuring both qubits must give the same outcome
    tab = GeneralizedTableau(2)
    tab.h(0)
    tab.cnot(0, 1)
    result0 = tab.measure(0)
    result1 = tab.measure(1)
    assert result0 == result1


def test_ghz_state():
    # H(0), CNOT(0,1), CNOT(0,2) → (|000> + |111>) / sqrt(2)
    tab = GeneralizedTableau(3)
    tab.h(0)
    tab.cnot(0, 1)
    tab.cnot(0, 2)
    result0 = tab.measure(0)
    result1 = tab.measure(1)
    result2 = tab.measure(2)
    assert result0 == result1 == result2


def test_reset():
    # After X, reset should bring qubit back to |0>
    tab = GeneralizedTableau(2)
    tab.x(0)
    tab.reset(0)
    assert not tab.measure(0)


def test_reset_on_zero_state():
    tab = GeneralizedTableau(2)
    tab.reset(0)
    assert not tab.measure(0)


def test_rx_pi_equals_x():
    # rx(π) ≈ -iX, so |0> → |1>
    tab = GeneralizedTableau(2)
    tab.rx(0, math.pi)
    assert tab.measure(0)


def test_rz_on_zero_state():
    # rz(θ)|0> = e^{-iθ/2}|0>, Z-measurement unchanged for any θ
    tab = GeneralizedTableau(2)
    tab.rz(0, math.pi / 4)
    assert not tab.measure(0)


def test_ry_pi_equals_y():
    # ry(π) ≈ -iY, so |0> → |1>
    tab = GeneralizedTableau(2)
    tab.ry(0, math.pi)
    assert tab.measure(0)


def test_r_axis_zero_equals_rx():
    # r(axis_angle=0, θ=π) = rx(π), so |0> → |1>
    tab = GeneralizedTableau(2)
    tab.r(0, 0.0, math.pi)
    assert tab.measure(0)


def test_r_axis_half_pi_equals_ry():
    # r(axis_angle=π/2, θ=π) = ry(π), so |0> → |1>
    tab = GeneralizedTableau(2)
    tab.r(0, math.pi / 2, math.pi)
    assert tab.measure(0)


def test_clifford_extensions_sqrt_x():
    # sqrt_x followed by sqrt_x_adj is identity
    tab = GeneralizedTableau(2)
    tab.x(0)
    tab.sqrt_x(0)
    tab.sqrt_x_adj(0)
    assert tab.measure(0)


def test_clifford_extensions_sqrt_y():
    # sqrt_y followed by sqrt_y_adj is identity
    tab = GeneralizedTableau(2)
    tab.x(0)
    tab.sqrt_y(0)
    tab.sqrt_y_adj(0)
    assert tab.measure(0)


def test_s_adj():
    # S S† = identity
    tab = GeneralizedTableau(2)
    tab.x(0)
    tab.s(0)
    tab.s_adj(0)
    assert tab.measure(0)


def test_t_gate_on_zero_state():
    # T|0> = |0>, Z-measurement unchanged
    tab = GeneralizedTableau(2)
    tab.t(0)
    assert not tab.measure(0)


def test_t_gate_four_times_equals_z():
    # T^4 = Z, so T^4|1> = Z|1> = -|1>, still measures as 1
    tab = GeneralizedTableau(2)
    tab.x(0)
    for _ in range(4):
        tab.t(0)
    assert tab.measure(0)


def test_t_adj_inverse_of_t():
    # T T† = identity
    tab = GeneralizedTableau(2)
    tab.x(0)
    tab.t(0)
    tab.t_adj(0)
    assert tab.measure(0)
