from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

import ppvm_python_native

if TYPE_CHECKING:
    from .paulisum import PauliSum


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
    state: "PauliSum",
    lindblad: LindbladOp,
    t_span: tuple[float, float],
    save_at: Sequence[float],
    *,
    hamiltonian: "PauliSum | None" = None,
    observable: "str | list[str] | None" = None,
    config: SolverConfig | None = None,
) -> "tuple[list[float], list]":
    """Solve the Lindblad master equation.

    Args:
        state: Initial density-matrix state as a PauliSum.
        lindblad: Dissipation operator (jump ops + rate matrix).
        t_span: (t_start, t_end) integration interval.
        save_at: Times at which to record results. Must be non-empty,
            sorted ascending, and within t_span.
        hamiltonian: Optional coherent Hamiltonian (same type as state).
        observable: Controls what is returned at each save point.
            - None → list[PauliSum] (full state snapshots)
            - "trace:<pattern>" → list[float] (single scalar)
            - ["trace:<p1>", "trace:<p2>", ...] → list[list[float]]
        config: ODE solver parameters (tolerances, step sizes).

    Returns:
        (times, results) where times are the actual save times and
        results depend on the observable mode.
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
    native_state = state._interface
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
    if observable is None:
        times, raw_states = ppvm_python_native.solve_timeevolve_states(
            state=native_state, **kwargs
        )
        return times, [_wrap_native(s) for s in raw_states]

    # scalar observable mode
    single = isinstance(observable, str)
    obs_list = [observable] if single else list(observable)
    patterns = []
    for obs in obs_list:
        if not obs.startswith("trace:"):
            raise ValueError(
                f"unsupported observable {obs!r}: only 'trace:<pattern>' is supported"
            )
        patterns.append(obs[len("trace:") :])

    times, results = ppvm_python_native.solve_timeevolve_observables(
        state=native_state, patterns=patterns, **kwargs
    )
    if single:
        # unwrap inner list: list[list[float]] → list[float]
        return times, [row[0] for row in results]
    return times, results
