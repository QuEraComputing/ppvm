from ppvm import GeneralizedTableau
from ppvm.generalized_tableau import MeasurementResult


def test_is_lost_initially_false():
    tab = GeneralizedTableau(n_qubits=3)
    for i in range(3):
        assert not tab.is_lost(i)


def test_is_lost_after_loss_channel():
    # With p=1 the qubit is always lost
    tab = GeneralizedTableau(n_qubits=2, seed=0)
    tab.loss_channel(0, 1.0)
    assert tab.is_lost(0)
    assert not tab.is_lost(1)


def test_is_lost_reset_loss_channel():
    tab = GeneralizedTableau(n_qubits=1, seed=0)
    tab.loss_channel(0, 1.0)
    assert tab.is_lost(0)
    tab.reset_loss_channel(0)
    assert not tab.is_lost(0)


def test_loss_values_initially_all_false():
    n = 4
    tab = GeneralizedTableau(n_qubits=n)
    assert tab.loss_values() == [False] * n


def test_loss_values_length():
    n = 5
    tab = GeneralizedTableau(n_qubits=n)
    assert len(tab.loss_values()) == n


def test_loss_values_after_loss_channel():
    tab = GeneralizedTableau(n_qubits=3, seed=0)
    tab.loss_channel(1, 1.0)
    values = tab.loss_values()
    assert values[0] is False
    assert values[1] is True
    assert values[2] is False


def test_measure_zero():
    tab = GeneralizedTableau(n_qubits=1)
    assert tab.measure(0) == MeasurementResult.ZERO


def test_measure_one():
    tab = GeneralizedTableau(n_qubits=1)
    tab.x(0)
    assert tab.measure(0) == MeasurementResult.ONE


def test_measure_lost():
    tab = GeneralizedTableau(n_qubits=1, seed=0)
    tab.loss_channel(0, 1.0)
    assert tab.measure(0) == MeasurementResult.LOST


# === CorrelatedLossChannel ===


def test_correlated_loss_p0_no_loss():
    tab = GeneralizedTableau(n_qubits=2)
    tab.correlated_loss_channel(0, 1, [0.0, 0.0, 0.0])
    assert not tab.is_lost(0)
    assert not tab.is_lost(1)


def test_correlated_loss_p0_both_lost():
    tab = GeneralizedTableau(n_qubits=2)
    tab.correlated_loss_channel(0, 1, [1.0, 0.0, 0.0])
    assert tab.is_lost(0)
    assert tab.is_lost(1)


def test_correlated_loss_p1_exactly_one_lost():
    # p[1]=1 → exactly one qubit lost in every trial.
    for seed in range(200):
        tab = GeneralizedTableau(n_qubits=2, seed=seed)
        tab.correlated_loss_channel(0, 1, [0.0, 1.0, 0.0])
        assert tab.is_lost(0) ^ tab.is_lost(1), f"Expected exactly one lost qubit (seed {seed})"


def test_correlated_loss_p1_both_qubits_chosen_equally():
    # With p[1]=1 the 50/50 coin flip should lose addr0 and addr1 equally.
    trials = 1000
    addr0_lost = sum(
        1
        for seed in range(trials)
        if (
            tab := GeneralizedTableau(n_qubits=2, seed=seed),
            tab.correlated_loss_channel(0, 1, [0.0, 1.0, 0.0]),
            tab.is_lost(0),
        )[-1]
    )
    fraction = addr0_lost / trials
    assert abs(fraction - 0.5) < 0.08, f"Expected ~0.5, got {fraction:.3f}"


def test_correlated_loss_both_lost_resets_to_zero():
    # Lost qubits should be reset to |0⟩.
    tab = GeneralizedTableau(n_qubits=2, seed=0)
    tab.x(0)
    tab.x(1)
    tab.correlated_loss_channel(0, 1, [1.0, 0.0, 0.0])
    assert tab.is_lost(0)
    assert tab.is_lost(1)
    tab.reset_loss_channel(0)
    tab.reset_loss_channel(1)
    assert tab.measure(0) == MeasurementResult.ZERO
    assert tab.measure(1) == MeasurementResult.ZERO


def test_correlated_loss_single_lost_resets_to_zero():
    # The lost qubit should be reset to |0⟩.
    for seed in range(1000):
        tab = GeneralizedTableau(n_qubits=2, seed=seed)
        tab.x(0)
        tab.correlated_loss_channel(0, 1, [0.0, 1.0, 0.0])
        if tab.is_lost(0):
            tab.reset_loss_channel(0)
            assert tab.measure(0) == MeasurementResult.ZERO, "Lost qubit should be reset to |0⟩"
            return
    raise AssertionError("addr0 was never chosen as the lost qubit in 1000 trials")


def test_correlated_loss_addr0_already_lost_applies_p2_to_addr1():
    tab = GeneralizedTableau(n_qubits=2, seed=0)
    tab.loss_channel(0, 1.0)
    tab.correlated_loss_channel(0, 1, [0.0, 0.0, 1.0])
    assert tab.is_lost(0)
    assert tab.is_lost(1)


def test_correlated_loss_addr1_already_lost_applies_p2_to_addr0():
    tab = GeneralizedTableau(n_qubits=2, seed=0)
    tab.loss_channel(1, 1.0)
    tab.correlated_loss_channel(0, 1, [0.0, 0.0, 1.0])
    assert tab.is_lost(0)
    assert tab.is_lost(1)


def test_correlated_loss_addr0_already_lost_p2_zero_addr1_survives():
    tab = GeneralizedTableau(n_qubits=2, seed=0)
    tab.loss_channel(0, 1.0)
    tab.correlated_loss_channel(0, 1, [0.0, 0.0, 0.0])
    assert not tab.is_lost(1)


def test_correlated_loss_statistics_both():
    # P(both lost) should converge to p[0].
    p_both = 0.3
    trials = 1000
    both_lost = sum(
        1
        for seed in range(trials)
        if (
            tab := GeneralizedTableau(n_qubits=2, seed=seed),
            tab.correlated_loss_channel(0, 1, [p_both, 0.0, 0.0]),
            tab.is_lost(0) and tab.is_lost(1),
        )[-1]
    )
    fraction = both_lost / trials
    assert abs(fraction - p_both) < 0.07, f"Expected ~{p_both:.2f}, got {fraction:.3f}"


def test_correlated_loss_statistics_single():
    # P(exactly one lost) should converge to p[1].
    p_single = 0.4
    trials = 1000
    one_lost = sum(
        1
        for seed in range(trials)
        if (
            tab := GeneralizedTableau(n_qubits=2, seed=seed),
            tab.correlated_loss_channel(0, 1, [0.0, p_single, 0.0]),
            tab.is_lost(0) ^ tab.is_lost(1),
        )[-1]
    )
    fraction = one_lost / trials
    assert abs(fraction - p_single) < 0.08, f"Expected ~{p_single:.2f}, got {fraction:.3f}"
