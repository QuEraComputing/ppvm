from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

import ppvm_python_native

if TYPE_CHECKING:
    from .paulisum import PauliSum
    from .product_state import ProductState


@dataclass
class LadderOp:
    qubit: int
    direction: str  # "raise" or "lower"


@dataclass
class LindbladOp:
    jump_ops: list[LadderOp]
    rates: list[float] | list[list[float]]


@dataclass
class SolverConfig:
    rtol: float = 1e-6
    atol: float = 1e-9
    h0: float | None = None
    hmin: float = 1e-12
    hmax: float = float("inf")


def _wrap_native(native_obj) -> "PauliSum":
    """Wrap a raw Rust PauliSumIndexMapFxHashN object in a PauliSum shell."""
    from .paulisum import PauliSum

    new = object.__new__(PauliSum)
    object.__setattr__(new, "initial_terms", [])
    object.__setattr__(new, "n_qubits", None)
    object.__setattr__(new, "coefficients", ())
    object.__setattr__(new, "min_abs_coeff", 1e-10)
    object.__setattr__(new, "max_pauli_weight", None)
    object.__setattr__(new, "max_loss_weight", None)
    object.__setattr__(new, "_interface", native_obj)
    return new


def solve(
    observable: "PauliSum",
    lindblad: LindbladOp,
    t_span: tuple[float, float],
    save_at: Sequence[float],
    *,
    hamiltonian: "PauliSum | None" = None,
    initial_state: "ProductState | None" = None,
    config: SolverConfig | None = None,
) -> "tuple[list[float], list]":
    """Solve the Heisenberg-picture adjoint master equation.

    Propagates an observable O under dO/dt = i[H, O] + L†(O).
    To obtain expectation values, supply an `initial_state` ρ₀; the solver then
    returns ⟨O(t)⟩ = Tr(ρ₀ O(t)) at each save point without cloning the full state.

    Args:
        observable: Initial observable O (PauliSum) in the Heisenberg picture.
        lindblad: Dissipation operator (jump ops + rate matrix).
        t_span: (t_start, t_end) integration interval.
        save_at: Times at which to record results. Must be non-empty,
            sorted ascending, and within t_span.
        hamiltonian: Optional coherent Hamiltonian (same type as observable).
        initial_state: ρ₀ for computing ⟨O(t)⟩ = Tr(ρ₀ O(t)).
        config: ODE solver parameters (tolerances, step sizes).

    Returns:
        (times, results) where:
            - `initial_state` given  → results is list[float] (expectation values)
            - neither given          → results is list[PauliSum] (raw snapshots)
    """
    # --- validation ---
    t0, t1 = t_span
    if t0 >= t1:
        raise ValueError(f"t_span must satisfy t_span[0] < t_span[1], got {t_span}")
    save_list = list(save_at)
    if not save_list:
        raise ValueError("save_at must be non-empty")
    if save_list != sorted(save_list):
        raise ValueError("save_at must be sorted in ascending order")
    if save_list[0] < t0 or save_list[-1] > t1:
        raise ValueError(
            f"All save_at times must be within t_span={t_span}, "
            f"got [{save_list[0]}, {save_list[-1]}]"
        )
    for op in lindblad.jump_ops:
        if op.direction not in ("raise", "lower"):
            raise ValueError(
                f"invalid direction {op.direction!r}: expected 'raise' or 'lower'"
            )

    if config is None:
        config = SolverConfig()

    # --- build native args ---
    native_observable = observable._interface
    native_ham = hamiltonian._interface if hamiltonian is not None else None
    ops_list = [(op.qubit, op.direction) for op in lindblad.jump_ops]
    rates = lindblad.rates

    kwargs = dict(
        lindblad_ops=ops_list,
        rates=rates,
        t_span_start=t0,
        t_span_end=t1,
        save_at=save_list,
        hamiltonian=native_ham,
        rtol=config.rtol,
        atol=config.atol,
        h0=config.h0,
        hmin=config.hmin,
        hmax=config.hmax,
    )

    # --- dispatch ---
    if initial_state is not None:
        # Fast path: no state clone, single f64 per save point computed in Rust.
        times, results = ppvm_python_native.solve_timeevolve_expectation(
            observable=native_observable, bloch_vectors=initial_state._bloch, **kwargs
        )
        return times, results

    # Default: raw PauliSum snapshots.
    times, raw_states = ppvm_python_native.solve_timeevolve_states(
        observable=native_observable, **kwargs
    )
    return times, [_wrap_native(s) for s in raw_states]
