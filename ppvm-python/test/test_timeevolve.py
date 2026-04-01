"""Tests for ppvm.timeevolve Python bindings (Heisenberg-picture API).

All tests use the Heisenberg picture: the observable is propagated under
dO/dt = i[H, O] + L†(O). Expectation values are computed via ProductState.
"""

import math

import pytest
from ppvm.timeevolve import LadderOp, LindbladOp, solve

from ppvm import PauliSum, ProductState

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


def _z0_observable() -> PauliSum:
    """Observable O = Z0 (single-qubit Pauli Z on qubit 0 of a 2-qubit system)."""
    return PauliSum.new(2, ["ZI"])


def _excited_2q_state() -> ProductState:
    """ρ₀ = |1⟩⊗|1⟩: bz = -1 for both qubits (excited state)."""
    return ProductState.bitstring("11")


# ---------------------------------------------------------------------------
# Test 1: decay state snapshots (Heisenberg picture returns PauliSum list)
# ---------------------------------------------------------------------------


def test_decay_state_snapshots():
    observable = _z0_observable()
    lindblad = _make_2q_lower_lindblad(gamma=1.0)
    save_at = [0.5, 1.0, 2.0, 4.0]

    times, states = solve(
        observable=observable, lindblad=lindblad, t_span=(0.0, 5.0), save_at=save_at
    )

    assert len(times) == 4
    assert len(states) == 4
    # States are PauliSum objects.
    assert all(isinstance(s, PauliSum) for s in states)


# ---------------------------------------------------------------------------
# Test 2: expectation value path matches raw-snapshot path
# ---------------------------------------------------------------------------


def test_decay_expectation_matches_snapshot():
    observable = _z0_observable()
    lindblad = _make_2q_lower_lindblad(gamma=1.0)
    save_at = [0.5, 1.0, 2.0, 4.0]
    rho0 = _excited_2q_state()

    _, states = solve(
        observable=observable, lindblad=lindblad, t_span=(0.0, 5.0), save_at=save_at
    )
    _, values = solve(
        observable=observable,
        lindblad=lindblad,
        t_span=(0.0, 5.0),
        save_at=save_at,
        initial_state=rho0,
    )

    assert isinstance(values, list)
    assert len(values) == len(save_at)
    for v, s in zip(values, states):
        expected = rho0.expectation(s)
        assert (
            abs(v - expected) < 1e-9
        ), f"expectation path {v} != snapshot path {expected}"


# ---------------------------------------------------------------------------
# Test 3: no Hamiltonian runs without error
# ---------------------------------------------------------------------------


def test_no_hamiltonian():
    observable = _z0_observable()
    lindblad = _make_2q_lower_lindblad()
    times, states = solve(
        observable=observable,
        lindblad=lindblad,
        t_span=(0.0, 1.0),
        save_at=[0.5, 1.0],
        hamiltonian=None,
    )
    assert len(times) == 2
    assert len(states) == 2


# ---------------------------------------------------------------------------
# Test 4: Hamiltonian-driven Larmor precession on 1 qubit
# ---------------------------------------------------------------------------


def test_with_hamiltonian():
    # Observable O = Z (Heisenberg picture); H = π/2 · X drives Z → -Z at t=1.
    observable = PauliSum.new(1, ["Z"])
    ham = PauliSum.new(1, [("X", math.pi / 2)])
    lindblad = LindbladOp(jump_ops=[], rates=[])
    # ρ₀ = |0⟩: bz = +1, so ⟨Z(0)⟩ = +1.
    rho0 = ProductState.all_zero(1)

    save_at = [0.25, 0.5, 0.75, 1.0]
    _, values = solve(
        observable=observable,
        lindblad=lindblad,
        t_span=(0.0, 1.0),
        save_at=save_at,
        hamiltonian=ham,
        initial_state=rho0,
    )

    # At t=1, ⟨Z⟩ ≈ -1 (full Rabi flip).
    assert values[-1] < -0.9, f"Expected ⟨Z⟩ ≈ -1 at t=1 but got {values[-1]}"
    # At t=0.5, ⟨Z⟩ ≈ 0 (halfway).
    assert abs(values[1]) < 0.01, f"Expected ⟨Z⟩ ≈ 0 at t=0.5 but got {values[1]}"


# ---------------------------------------------------------------------------
# Test 5-7: save_at validation
# ---------------------------------------------------------------------------


def test_save_at_validation_empty():
    observable = _z0_observable()
    lindblad = _make_2q_lower_lindblad()
    with pytest.raises(ValueError, match="save_at must be non-empty"):
        solve(observable=observable, lindblad=lindblad, t_span=(0.0, 1.0), save_at=[])


def test_save_at_validation_unsorted():
    observable = _z0_observable()
    lindblad = _make_2q_lower_lindblad()
    with pytest.raises(ValueError, match="sorted"):
        solve(
            observable=observable,
            lindblad=lindblad,
            t_span=(0.0, 3.0),
            save_at=[2.0, 1.0],
        )


def test_save_at_validation_out_of_bounds():
    observable = _z0_observable()
    lindblad = _make_2q_lower_lindblad()
    with pytest.raises(ValueError, match="t_span"):
        solve(
            observable=observable,
            lindblad=lindblad,
            t_span=(0.0, 1.0),
            save_at=[2.0],
        )


# ---------------------------------------------------------------------------
# Test 8: type mismatch raises TypeError
# ---------------------------------------------------------------------------


def test_type_mismatch():
    observable = PauliSum.new(2, [("II", 1.0)])  # 2-qubit → N=1 interface
    ham = PauliSum.new(16, [("I" * 16, 1.0)])  # 16-qubit → N=4 interface
    lindblad = _make_2q_lower_lindblad()
    with pytest.raises(TypeError):
        solve(
            observable=observable,
            lindblad=lindblad,
            t_span=(0.0, 1.0),
            save_at=[0.5],
            hamiltonian=ham,
        )


# ---------------------------------------------------------------------------
# Test 9: returned observables are PauliSum with working expectation
# ---------------------------------------------------------------------------


def test_returned_states_are_paulisum():
    observable = _z0_observable()
    lindblad = _make_2q_lower_lindblad()

    _, states = solve(
        observable=observable, lindblad=lindblad, t_span=(0.0, 1.0), save_at=[0.5]
    )
    s = states[0]

    assert isinstance(s, PauliSum)
    # The observable (ZI) evolves — verify we can trace it with a wildcard pattern.
    tr = s.trace("Z?*")
    assert isinstance(tr, float)


# ---------------------------------------------------------------------------
# Test 10: dense rate matrix gives same result as diagonal
# ---------------------------------------------------------------------------


def test_dense_rate_matrix():
    observable = _z0_observable()
    gamma = 0.5
    save_at = [0.5, 1.0]
    rho0 = _excited_2q_state()

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
        observable=observable,
        lindblad=diag_lindblad,
        t_span=(0.0, 2.0),
        save_at=save_at,
        initial_state=rho0,
    )
    _, vals_dense = solve(
        observable=observable,
        lindblad=dense_lindblad,
        t_span=(0.0, 2.0),
        save_at=save_at,
        initial_state=rho0,
    )

    for vd, vn in zip(vals_diag, vals_dense):
        assert (
            abs(vd - vn) < 1e-9
        ), f"diagonal={vd} vs dense={vn} differ by {abs(vd - vn)}"
