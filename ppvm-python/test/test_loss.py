from ppvm import LossyPauliSum


def test_ghz():
    ps = LossyPauliSum.new(n_qubits=2, terms=["ZZ"])

    # loss at the end of the circuit
    p = 0.1
    ps.reset_loss_channel(0)
    ps.loss_channel(0, p)

    ps.cnot(0, 1)
    ps.h(0)

    # losing one atom at the end of the circuit in GHZ
    # means that ZZ will be anti-correlated in p cases of the 50% |11> state
    # since there is a factor of two from adding up flipped signs, but a factor
    # of 1/2 since it's only the |11> part, we should arrive at 1 - p
    assert abs(ps.overlap_with_zero() - (1 - p)) < 1e-9


