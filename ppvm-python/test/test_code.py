# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

import numpy as np
import pytest

from ppvm import Code

FIVE = ("+XZZXI", "+IXZZX", "+XIXZZ", "+ZXIXZ")
STEANE_H = [[1, 1, 1, 1, 0, 0, 0], [0, 1, 1, 0, 1, 1, 0], [0, 0, 1, 1, 1, 0, 1]]


def test_non_css_tableau_and_target():
    xs = [[p in "XY" for p in row[1:]] for row in FIVE]
    zs = [[p in "ZY" for p in row[1:]] for row in FIVE]
    code = Code.from_tableau(xs, zs, logical_xs=["XXXXX"], logical_zs=["ZZZZZ"])
    assert (code.n, code.k, code.stabilizers) == (5, 1, FIVE)
    assert code.state_stabilizers() == (*FIVE, "+ZZZZZ")
    assert code.state_stabilizers([1]) == (*FIVE, "-ZZZZZ")


def test_steane_and_empty_checks():
    code = Code.from_css(STEANE_H, STEANE_H)
    assert (code.n, code.k) == (7, 1)
    assert len(code.state_stabilizers()) == 7
    empty = np.zeros((0, 3), dtype=bool)
    code = Code.from_tableau(empty, empty)
    assert (code.n, code.k) == (3, 3)


def test_y_sign_and_redundant_rows():
    code = Code.from_tableau([[1, 1], [0, 0], [1, 1]], [[0, 0], [1, 1], [1, 1]], [0, 0, 1])
    assert code.stabilizers == ("+XX", "+ZZ", "-YY")
    assert code.state_stabilizers() == ("+XX", "+ZZ")
    with pytest.raises(ValueError, match="-I"):
        Code.from_stabilizers(["XX", "ZZ", "YY"])


@pytest.mark.parametrize(
    ("args", "error"),
    [
        (([[2]], [[0]]), "binary"),
        (([[1]], [[1, 1]]), "shape"),
        (([[1]], [[1]], []), "sign"),
        (([[1], [0]], [[0], [1]]), "anticommute"),
    ],
)
def test_bad_tableau(args, error):
    with pytest.raises(ValueError, match=error):
        Code.from_tableau(*args)


def test_signed_logical_and_bad_inputs():
    code = Code.from_stabilizers(["-ZZ"], logical_xs=["XX"], logical_zs=["-ZI"])
    assert code.state_stabilizers([1]) == ("-ZZ", "+ZI")
    for xs, zs in [(["XI"], ["ZI"]), (["ZZ"], ["IZ"]), ([], [])]:
        with pytest.raises(ValueError):
            Code.from_stabilizers(["ZZ"], logical_xs=xs, logical_zs=zs)
    with pytest.raises(ValueError, match="invalid signed Pauli"):
        Code.from_stabilizers(["+iX"])


def test_non_css_logical_one_preparation():
    import stim

    code = Code.from_stabilizers(FIVE, logical_xs=["XXXXX"], logical_zs=["ZZZZZ"])
    circuit = code.synthesize_state_prep([1])
    sim = stim.TableauSimulator()
    sim.do_circuit(circuit)
    for row in code.state_stabilizers([1]):
        assert sim.peek_observable_expectation(stim.PauliString(row)) == 1


def test_terminal_readout_combinations_and_offsets():
    # Neither original check is directly measurable in Y, but their signed
    # product -YY is deterministic.
    bell = Code.from_stabilizers(["XX", "ZZ"])
    assert bell.terminal_readout("YY").detectors == (((0, 1), True),)
    odd = Code.from_stabilizers(["-ZZ"], logical_xs=["XX"], logical_zs=["ZI"])
    readout = odd.terminal_readout("XZ", order=[1, 0])
    # ZI is equivalent to -IZ modulo the signed stabilizer.
    assert readout.detectors == ()
    assert readout.observables == (((0,), True),)
    with pytest.raises(ValueError, match="inaccessible"):
        odd.terminal_readout("XX")
    with pytest.raises(ValueError, match="permutation"):
        odd.terminal_readout("ZZ", order=[0, 0])
    one = Code.from_stabilizers([], n_qubits=1, logical_xs=["-X"], logical_zs=["Z"])
    assert one.terminal_readout("Y", logical_bases="Y").observables == (((0,), True),)


def test_encoder_tableau_exact_logical_basis():
    import stim

    code = Code.from_stabilizers(["-ZZI"], logical_xs=["-XXI", "IIY"], logical_zs=["ZII", "-IIZ"])
    tab = code.encoder_tableau()
    assert [tab.z_output(i) for i in range(3)] == [
        stim.PauliString(s) for s in ("-ZZI", "ZII", "-IIZ")
    ]
    assert [tab.x_output(i) for i in (1, 2)] == [stim.PauliString(s) for s in code.logical_xs]


def test_encoder_and_logical_gadget():
    import stim

    code = Code.from_stabilizers(["-ZZ"], logical_xs=["XX"], logical_zs=["ZI"])
    encoder = stim.Tableau.from_circuit(code.synthesize_encoder())
    assert encoder == code.encoder_tableau()
    gadget = code.synthesize_logical_clifford(stim.Tableau.from_named_gate("S"))
    physical = stim.Tableau.from_circuit(gadget)
    assert physical(stim.PauliString("-ZZ")) == stim.PauliString("-ZZ")
    z = stim.PauliString("ZI")
    x = stim.PauliString("XX")
    assert physical(z) == z
    assert physical(x) == 1j * x * z


def test_gadget_preserves_group_and_logical_action():
    import stim

    code = Code.from_stabilizers([stim.PauliString("-ZZ")], logical_xs=["XX"], logical_zs=["ZI"])
    assert code.contains_stabilizer("-ZZ")
    assert not code.contains_stabilizer("ZZ")
    identity = stim.Tableau(1)
    assert code.implements_logical_clifford(stim.Tableau(2), identity)
    # SWAP preserves the odd-parity code but interchanges its two words.
    swap = stim.Tableau.from_named_gate("SWAP")
    assert not code.implements_logical_clifford(swap, identity)
    assert code.implements_logical_clifford(swap, stim.Tableau.from_named_gate("X"))


def test_two_steane_blocks_terminal_readout():
    h = np.array(STEANE_H)
    two = np.kron(np.eye(2, dtype=int), h)
    support = "XXIIIXI"
    xs = [support + "I" * 7, "I" * 7 + support]
    code = Code.from_css(two, two, logical_xs=xs, logical_zs=[x.replace("X", "Z") for x in xs])
    layout = list(range(7)) + list(range(13, 6, -1))
    readout = code.terminal_readout("Z" * 14, order=layout)
    assert [p.indices for p in readout.detectors] == [
        tuple(sorted(layout.index(q) for q in row.nonzero()[0])) for row in two
    ]
    assert readout.observables == (((0, 1, 5), False), ((8, 12, 13), False))


def test_idle_wires_survive_synthesis():
    import stim

    code = Code.from_stabilizers([], n_qubits=2, logical_xs=["XI", "IX"], logical_zs=["ZI", "IZ"])
    assert code.synthesize_state_prep().num_qubits == 2
    assert code.synthesize_encoder().to_tableau() == stim.Tableau(2)
