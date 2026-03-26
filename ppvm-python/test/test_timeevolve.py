"""Tests for ppvm.timeevolve Python bindings."""

import math

import pytest
from ppvm.timeevolve import LadderOp, LindbladOp, solve

from ppvm import PauliSum

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_2q_lower_lindblad(gamma: float = 1.0) -> LindbladOp:
    return LindbladOp(
        jump_ops=[
            LadderOp(qubit=0, direction="lower"),
            LadderOp(qubit=1, direction="lower"),
        ],
        rates=[gamma, gamma],
    )


def _excited_2q() -> PauliSum:
    # (I+Z)/2 ⊗ (I+Z)/2 — positive Z coefficients for both qubits.
    # Under lowering operators the Z terms decay toward zero.
    return PauliSum.new(2, [("II", 0.25), ("IZ", 0.25), ("ZI", 0.25), ("ZZ", 0.25)])


# ---------------------------------------------------------------------------
# Test 1: decay state snapshots
# ---------------------------------------------------------------------------


def test_decay_state_snapshots():
    state = _excited_2q()
    lindblad = _make_2q_lower_lindblad(gamma=1.0)
    save_at = [0.5, 1.0, 2.0, 4.0]

    times, states = solve(
        state=state, lindblad=lindblad, t_span=(0.0, 5.0), save_at=save_at
    )

    assert len(times) == 4
    assert len(states) == 4

    # qubit-0 Z coefficient should decay monotonically under lowering operators
    z0_vals = [s.trace("Z0") for s in states]
    for i in range(len(z0_vals) - 1):
        assert (
            z0_vals[i] > z0_vals[i + 1]
        ), f"Expected monotonic decay but got {z0_vals}"


# ---------------------------------------------------------------------------
# Test 2: scalar observable mode matches state-snapshot mode
# ---------------------------------------------------------------------------


def test_decay_scalar_observable():
    state = _excited_2q()
    lindblad = _make_2q_lower_lindblad(gamma=1.0)
    save_at = [0.5, 1.0, 2.0, 4.0]

    _, states = solve(
        state=state, lindblad=lindblad, t_span=(0.0, 5.0), save_at=save_at
    )
    _, values = solve(
        state=state,
        lindblad=lindblad,
        t_span=(0.0, 5.0),
        save_at=save_at,
        observable="trace:Z0",
    )

    assert isinstance(values, list)
    assert len(values) == len(save_at)
    for v, s in zip(values, states):
        assert (
            abs(v - s.trace("Z0")) < 1e-9
        ), f"scalar {v} != state snapshot trace {s.trace('Z0')}"


# ---------------------------------------------------------------------------
# Test 3: no Hamiltonian runs without error
# ---------------------------------------------------------------------------


def test_no_hamiltonian():
    state = _excited_2q()
    lindblad = _make_2q_lower_lindblad()
    times, states = solve(
        state=state,
        lindblad=lindblad,
        t_span=(0.0, 1.0),
        save_at=[0.5, 1.0],
        hamiltonian=None,
    )
    assert len(times) == 2
    assert len(states) == 2


# ---------------------------------------------------------------------------
# Test 4: Hamiltonian-driven Rabi oscillation on 1 qubit
# ---------------------------------------------------------------------------


def test_with_hamiltonian():
    # Initial state: Z0 coefficient = 0.5 (fully in +Z eigenstate)
    state = PauliSum.new(1, [("I", 0.5), ("Z", 0.5)])
    # H = pi/2 * X  → Rabi flip: Z0 goes from 0.5 → ~0 at t=0.5 → ~-0.5 at t=1
    ham = PauliSum.new(1, [("X", math.pi / 2)])
    lindblad = LindbladOp(jump_ops=[], rates=[])

    save_at = [0.25, 0.5, 0.75, 1.0]
    _, values = solve(
        state=state,
        lindblad=lindblad,
        t_span=(0.0, 1.0),
        save_at=save_at,
        hamiltonian=ham,
        observable="trace:Z0",
    )

    # At t=1, Z0 coefficient should be near -0.5 (flipped)
    assert values[-1] < -0.45, f"Expected Z0 ≈ -0.5 at t=1 but got {values[-1]}"
    # At t=0.5, Z0 should be near 0 (halfway through flip)
    assert abs(values[1]) < 0.01, f"Expected Z0 ≈ 0 at t=0.5 but got {values[1]}"


# ---------------------------------------------------------------------------
# Test 5-7: save_at validation
# ---------------------------------------------------------------------------


def test_save_at_validation_empty():
    state = _excited_2q()
    lindblad = _make_2q_lower_lindblad()
    with pytest.raises(ValueError, match="save_at must be non-empty"):
        solve(state=state, lindblad=lindblad, t_span=(0.0, 1.0), save_at=[])


def test_save_at_validation_unsorted():
    state = _excited_2q()
    lindblad = _make_2q_lower_lindblad()
    with pytest.raises(ValueError, match="sorted"):
        solve(state=state, lindblad=lindblad, t_span=(0.0, 3.0), save_at=[2.0, 1.0])


def test_save_at_validation_out_of_bounds():
    state = _excited_2q()
    lindblad = _make_2q_lower_lindblad()
    with pytest.raises(ValueError, match="t_span"):
        solve(state=state, lindblad=lindblad, t_span=(0.0, 1.0), save_at=[2.0])


# ---------------------------------------------------------------------------
# Test 8: type mismatch raises TypeError
# ---------------------------------------------------------------------------


def test_type_mismatch():
    state = PauliSum.new(2, [("II", 1.0)])  # 2-qubit  → N=1 interface
    ham = PauliSum.new(16, [("I" * 16, 1.0)])  # 16-qubit → N=4 interface
    lindblad = _make_2q_lower_lindblad()
    with pytest.raises(TypeError):
        solve(
            state=state,
            lindblad=lindblad,
            t_span=(0.0, 1.0),
            save_at=[0.5],
            hamiltonian=ham,
        )


# ---------------------------------------------------------------------------
# Test 9: returned states are PauliSum with working trace
# ---------------------------------------------------------------------------


def test_returned_states_are_paulisum():
    state = _excited_2q()
    lindblad = _make_2q_lower_lindblad()

    _, states = solve(state=state, lindblad=lindblad, t_span=(0.0, 1.0), save_at=[0.5])
    s = states[0]

    assert isinstance(s, PauliSum)
    # trace("Z?*") sums all Z/I coefficients, conserved at 1.0
    tr = s.trace("Z?*")
    assert isinstance(tr, float)
    assert abs(tr - 1.0) < 1e-6


# ---------------------------------------------------------------------------
# Test 10: dense rate matrix gives same result as diagonal
# ---------------------------------------------------------------------------


def test_dense_rate_matrix():
    state = _excited_2q()
    gamma = 0.5
    save_at = [0.5, 1.0]

    diag_lindblad = LindbladOp(
        jump_ops=[
            LadderOp(qubit=0, direction="lower"),
            LadderOp(qubit=1, direction="lower"),
        ],
        rates=[gamma, gamma],
    )
    dense_lindblad = LindbladOp(
        jump_ops=[
            LadderOp(qubit=0, direction="lower"),
            LadderOp(qubit=1, direction="lower"),
        ],
        rates=[[gamma, 0.0], [0.0, gamma]],
    )

    _, vals_diag = solve(
        state=state,
        lindblad=diag_lindblad,
        t_span=(0.0, 2.0),
        save_at=save_at,
        observable="trace:Z0",
    )
    _, vals_dense = solve(
        state=state,
        lindblad=dense_lindblad,
        t_span=(0.0, 2.0),
        save_at=save_at,
        observable="trace:Z0",
    )

    for vd, vn in zip(vals_diag, vals_dense):
        assert (
            abs(vd - vn) < 1e-9
        ), f"diagonal={vd} vs dense={vn} differ by {abs(vd - vn)}"
