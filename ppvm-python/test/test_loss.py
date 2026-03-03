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


def test_correlated_loss_ghz():
    # With only correlated loss (p[0]=0.1, p[1]=0, p[2]=0), both qubits are
    # always lost together, so measurement outcomes stay perfectly correlated.
    p = [0.1, 0.0, 0.0]

    ps = LossyPauliSum.new(n_qubits=2, terms=["ZZ"])
    ps.reset_loss_channel(0)
    ps.reset_loss_channel(1)
    ps.correlated_loss_channel(0, 1, p)
    ps.cnot(0, 1)
    ps.h(0)
    assert abs(ps.overlap_with_zero() - 1.0) < 1e-9

    # XX is reduced by the correlated loss probability (both qubits lost → no
    # XX contribution), leaving a factor of (1 - p[0]).
    ps = LossyPauliSum.new(n_qubits=2, terms=["XX"])
    ps.reset_loss_channel(0)
    ps.reset_loss_channel(1)
    ps.correlated_loss_channel(0, 1, p)
    ps.cnot(0, 1)
    ps.h(0)
    assert abs(ps.overlap_with_zero() - (1 - p[0])) < 1e-9
