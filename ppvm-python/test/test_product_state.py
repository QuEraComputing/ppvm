import pytest

from ppvm import PauliSum, ProductState


def test_all_zero_n_qubits():
    assert ProductState.all_zero(4).n_qubits == 4


def test_bitstring_encoding():
    ps = ProductState.bitstring("01")
    assert ps._bloch == [0.0, 0.0, 1.0, 0.0, 0.0, -1.0]


def test_invalid_bitstring():
    with pytest.raises(ValueError):
        ProductState.bitstring("2")


def test_bloch_warning():
    with pytest.warns(UserWarning, match=r"\|b\|²"):
        ProductState.bloch_vectors([(2.0, 0.0, 0.0)])


def test_expectation_all_zero():
    # rho0 = |00⟩, observable = ZI + IZ (each coefficient 1).
    # Tr(|00⟩⟨00| ⊗ |0⟩⟨0| · ZI) = (+1)(1) = +1.
    # Tr(|00⟩⟨00| ⊗ |0⟩⟨0| · IZ) = (1)(+1) = +1.
    # Total = 2.
    rho0 = ProductState.all_zero(2)
    obs = PauliSum.new(2, ["ZI", "IZ"])
    result = rho0.expectation(obs)
    assert abs(result - 2.0) < 1e-12
