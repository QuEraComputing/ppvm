# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

import enum
from dataclasses import InitVar, dataclass, field

import ppvm_python_native
from ppvm_python_native import StimProgram

from .mixins import (
    CliffordExtensionMixin,
    CliffordMixin,
    LossMixin,
    NoiseMixin,
    RotationsMixin,
)
from .types import GeneralizedTableauInterface

MAX_N_QUBITS = 2048
"""Maximum number of qubits supported by the Python bindings.

The native module pre-compiles a fixed set of tableau interfaces; beyond this
limit, use the Rust crate directly.
"""


def _native_tableau_cls(n_qubits: int):
    if n_qubits < 1:
        raise ValueError(
            f"n_qubits must be between 1 and {MAX_N_QUBITS} (got {n_qubits})."
        )
    if n_qubits > MAX_N_QUBITS:
        raise ValueError(
            f"n_qubits must be between 1 and {MAX_N_QUBITS} (got {n_qubits}); "
            "to simulate more qubits, use the ppvm-tableau Rust crate directly."
        )
    N_interface = (n_qubits + 63) // 64
    return getattr(ppvm_python_native, f"GeneralizedTableau{N_interface}")


class MeasurementResult(enum.IntEnum):
    """A measurement outcome, which accounts for a qubit being potentially lost."""

    ZERO = 0
    ONE = 1
    LOST = 2


@dataclass(frozen=True)
class GeneralizedTableau(
    CliffordMixin,
    CliffordExtensionMixin,
    RotationsMixin,
    NoiseMixin,
    LossMixin,
):
    """Generalized stabilizer tableau for quantum circuit simulation.

    Represents an arbitrary quantum state in the basis spanned by the
    stabilizer tableau. It supports Clifford gates, arbitrary single- and two-qubit rotations,
    noise channels, and mid-circuit measurement.
    The coefficient vector grows exponentially with the
    number of non-Clifford operations applied.

    Attributes:
        n_qubits: The number of qubits.
        min_abs_coeff: Coefficient threshold - coefficients smaller than this
            are pruned from the sparse coefficient vector.
        seed: Optional RNG seed for reproducible simulations. If ``None``
            (the default), the RNG is seeded from OS entropy.
    """

    n_qubits: int
    min_abs_coeff: float = 1e-10
    seed: InitVar[int | None] = None

    _interface: GeneralizedTableauInterface = field(init=False, repr=False)

    def __post_init__(self, seed: int | None):
        native_cls = _native_tableau_cls(self.n_qubits)
        object.__setattr__(
            self,
            "_interface",
            native_cls(self.n_qubits, self.min_abs_coeff, seed),
        )

    def fork(self, seed: int | None = None) -> "GeneralizedTableau":
        """Fork this tableau into an independent simulation branch.

        Clones all quantum state but reinitializes the RNG, so the returned
        tableau evolves independently from this one. If ``seed`` is provided
        the new RNG is seeded deterministically; otherwise it is seeded from
        OS entropy.

        Use this when branching a simulation into independent trajectories.
        To preserve the RNG state exactly (e.g. for checkpointing), use
        ``copy.copy()`` or ``copy.deepcopy()`` instead.

        Args:
            seed: Optional integer seed for the forked RNG.

        Returns:
            A new ``GeneralizedTableau`` with the same quantum state but an
            independent RNG.
        """
        forked = GeneralizedTableau(self.n_qubits, self.min_abs_coeff)
        object.__setattr__(forked, "_interface", self._interface.fork(seed))
        return forked

    def __copy__(self) -> "GeneralizedTableau":
        """Return a copy of this tableau, including its RNG state.

        Both the original and the copy will produce identical random sequences
        from this point forward. To get an independent copy with a fresh RNG,
        use :meth:`fork` instead.
        """
        copied = GeneralizedTableau(self.n_qubits, self.min_abs_coeff)
        object.__setattr__(copied, "_interface", self._interface.__copy__())
        return copied

    def __deepcopy__(self, memo: dict) -> "GeneralizedTableau":
        """Return a deep copy of this tableau, including its RNG state.

        Both the original and the copy will produce identical random sequences
        from this point forward. To get an independent copy with a fresh RNG,
        use :meth:`fork` instead.
        """
        copied = GeneralizedTableau(self.n_qubits, self.min_abs_coeff)
        object.__setattr__(copied, "_interface", self._interface.__deepcopy__(memo))
        return copied

    def __str__(self) -> str:
        """Return a human-readable representation of the tableau state."""
        return self._interface.__str__()

    def t(self, *targets: int) -> None:
        """Apply a T gate (π/8 rotation) to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.t(list(targets))

    def t_dag(self, *targets: int) -> None:
        """Apply a T adjoint gate (negative π/8 rotation) to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.t_dag(list(targets))

    def measure(self, addr0: int) -> MeasurementResult:
        """Measure the specified qubit in the Z basis.

        Args:
            addr0: The index of the target qubit.

        Returns:
            The measurement outcome as a ``MeasurementResult``, which is
            ``LOST`` if the qubit has been lost, ``ZERO`` or ``ONE`` otherwise.
        """
        return MeasurementResult(self._interface.measure(addr0))

    def measure_many(self, *targets: int) -> list[MeasurementResult]:
        """Measure several qubits in the Z basis.

        Args:
            *targets: The indices of the target qubits.

        Returns:
            A list of ``MeasurementResult`` outcomes, one per target.
        """
        return [MeasurementResult(v) for v in self._interface.measure_many(list(targets))]

    def current_measurement_record(self) -> list[MeasurementResult]:
        """Return all measurement outcomes recorded so far.

        Returns:
            A list of ``MeasurementResult`` outcomes in measurement order.
        """
        return [MeasurementResult(v) for v in self._interface.current_measurement_record()]

    def u3(self, addr0: int, theta: float, phi: float, lam: float):
        """Apply the U3 gate to the specified qubit.

        U3(θ, φ, λ) = RZ(φ) · RY(θ) · RZ(λ)

        The corresponding unitary matrix is:

            [ cos(θ/2)            -e^(iλ)·sin(θ/2)       ]
            [ e^(iφ)·sin(θ/2)     e^(i(φ+λ))·cos(θ/2)   ]

        Args:
            addr0: The address of the qubit to apply the gate to.
            theta: Rotation angle θ (in radians) for the RY component.
            phi: Rotation angle φ (in radians) for the first RZ component.
            lam: Rotation angle λ (in radians) for the second RZ component.
        """
        self._interface.u3(addr0, theta, phi, lam)

    def reset(self, *targets: int) -> None:
        """Reset each target qubit to the |0> state.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.reset(list(targets))

    def reset_x(self, *targets: int) -> None:
        """Reset each target qubit to the |+> state.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.reset_x(list(targets))

    def reset_y(self, *targets: int) -> None:
        """Reset each target qubit to the |+i> state.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.reset_y(list(targets))

    def reset_z(self, *targets: int) -> None:
        """Reset each target qubit to the |0> state.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.reset_z(list(targets))

    def reset_loss_channel(self, addr0: int) -> None:
        """Reset a lost qubit to being active again.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.reset_loss_channel(addr0)

    def is_lost(self, addr0: int) -> bool:
        """Check whether a qubit has been lost.

        Args:
            addr0: The index of the qubit.

        Returns:
            True if the qubit is lost, False otherwise.
        """
        return self._interface.is_lost(addr0)

    def loss_values(self) -> list[bool]:
        """Return the loss state of all qubits.

        Returns:
            A list of booleans of length ``n_qubits``, where each entry is
            True if the corresponding qubit is lost and False otherwise.
        """
        return self._interface.loss_values()

    def run(self, prog: StimProgram) -> list[MeasurementResult]:
        """Execute a parsed Stim program against this tableau (single shot).

        .. note::
            This **mutates** the tableau in place. For independent shots use
            :meth:`fork` or the :func:`ppvm.sample_stim` / :meth:`sample`
            helpers (which build a fresh tableau per shot).
        """
        raw = self._interface.run(prog)
        return [MeasurementResult(x) for x in raw]

    # stim familiarity alias
    do = run

    @classmethod
    def sample(
        cls,
        prog: StimProgram,
        n_qubits: int,
        min_abs_coeff: float = 1e-10,
        num_shots: int = 1,
        seed: int | None = None,
    ) -> list[list[MeasurementResult]]:
        """Run ``num_shots`` shots of ``prog`` and return all measurement results.

        Each shot starts from a fresh tableau, so this is the right entry
        point for multi-shot sampling.
        """
        native_cls = _native_tableau_cls(n_qubits)
        raw = native_cls.sample(prog, n_qubits, min_abs_coeff, num_shots, seed)
        return [[MeasurementResult(x) for x in shot] for shot in raw]


def sample_stim(
    prog: StimProgram,
    n_qubits: int,
    min_abs_coeff: float = 1e-10,
    num_shots: int = 1,
    seed: int | None = None,
) -> list[list[MeasurementResult]]:
    """Multi-shot sampling — module-level alias for ``GeneralizedTableau.sample``."""
    return GeneralizedTableau.sample(
        prog, n_qubits, min_abs_coeff=min_abs_coeff, num_shots=num_shots, seed=seed
    )
