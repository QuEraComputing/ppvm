"""Superradiance integration test — Heisenberg-picture framing.

The observable O = Σᵢ Zᵢ is propagated under the adjoint master equation
dO/dt = L†(O), where L† uses collective raising operators.  The initial state
ρ₀ = |0⟩^N (all qubits in |0⟩, bz=+1) is provided to compute expectation values.

At t=0 the exact result is:
    ⟨O(0)⟩ = Tr(|0⟩⟨0|^N · Σᵢ Zᵢ) = N  (each Zᵢ has eigenvalue +1 in |0⟩)
"""

import pytest
from ppvm.timeevolve import LadderOp, LindbladOp, solve

from ppvm import PauliSum, ProductState

N = 4
GAMMA = 1.0
TMAX = 2.0


def test_superradiance_expectation():
    # Observable O(0) = Σᵢ Zᵢ — propagated under dO/dt = L†(O).
    observable = PauliSum.new(N, [f"Z{i}" for i in range(N)])

    jump_ops = [LadderOp(i, direction="lower") for i in range(N)]
    lindblad = LindbladOp(jump_ops=jump_ops, rates=[GAMMA] * N)

    # ρ₀ = |0⟩^N — static; used only at save checkpoints.
    rho0 = ProductState.all_zero(N)

    save_at = [0.0, 0.5, 1.0, TMAX]
    times, results = solve(
        observable=observable,
        lindblad=lindblad,
        t_span=(0.0, TMAX),
        save_at=save_at,
        initial_state=rho0,
    )

    # Results must be a list of floats, not PauliSum objects.
    assert isinstance(results, list)
    assert all(isinstance(v, float) for v in results)

    # At t=0: ⟨Σᵢ Zᵢ⟩ = N (each qubit in |0⟩ gives +1).
    assert results[0] == pytest.approx(N, abs=1e-9)
