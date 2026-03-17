from ppvm import GeneralizedTableau


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
