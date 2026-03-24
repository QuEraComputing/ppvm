import enum
import math
from dataclasses import InitVar, dataclass, field
from typing import Optional, Sequence

import ppvm_python_native

from .clifford import (
    CliffordExtensionMixin,
    CliffordMixin,
    NoiseMixin,
    NonCliffordMixin,
)
from .types import GeneralizedTableauInterface


class MeasurementResult(enum.IntEnum):
    """A measurement outcome, which accounts for a qubit being potentially lost."""

    ZERO = 0
    ONE = 1
    LOST = 2

    @staticmethod
    def _from_raw(result: bool | None) -> "MeasurementResult":
        if result is None:
            return MeasurementResult.LOST
        if result:
            return MeasurementResult.ONE
        return MeasurementResult.ZERO


@dataclass(frozen=True)
class GeneralizedTableau(
    CliffordMixin, CliffordExtensionMixin, NonCliffordMixin, NoiseMixin
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
    seed: InitVar[Optional[int]] = None

    _interface: GeneralizedTableauInterface = field(init=False, repr=False)

    def __post_init__(self, seed: Optional[int]):
        N_interface = math.ceil(self.n_qubits / 64.0)
        object.__setattr__(
            self,
            "_interface",
            getattr(ppvm_python_native, f"GeneralizedTableau{N_interface}")(
                self.n_qubits, self.min_abs_coeff, seed
            ),
        )

    def fork(self, seed: Optional[int] = None) -> "GeneralizedTableau":
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

    def t(self, addr0: int) -> None:
        """Apply a T gate (π/8 rotation) to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.t(addr0)

    def t_adj(self, addr0: int) -> None:
        """Apply a T adjoint gate (negative π/8 rotation) to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.t_adj(addr0)

    # additional noise methods
    def depolarize(self, addr0: int, p: float) -> None:
        """Apply a depolarizing channel to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            p: The depolarizing probability.
        """
        self._interface.depolarize(addr0, p)

    def depolarize2(self, addr0: int, addr1: int, p: float) -> None:
        """Apply a two-qubit depolarizing channel to the specified qubits.

        Args:
            addr0: The index of the first target qubit.
            addr1: The index of the second target qubit.
            p: The depolarizing probability.
        """
        self._interface.depolarize2(addr0, addr1, p)

    def loss_channel(self, addr0: int, p: float) -> None:
        """Apply a loss channel to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            p: The loss probability.
        """
        self._interface.loss_channel(addr0, p)

    def correlated_loss_channel(
        self, addr0: int, addr1: int, p: Sequence[float]
    ) -> None:
        """Apply a correlated loss channel to two qubits.

        Args:
            addr0: The index of the first target qubit.
            addr1: The index of the second target qubit.
            p: A list of three probabilities:

                - ``p[0]``: probability of losing both qubits simultaneously
                  when both are in the qubit subspace.
                - ``p[1]``: probability of losing exactly one qubit when both
                  are in the qubit subspace (which qubit is lost is 50/50 random).
                - ``p[2]``: probability of losing the remaining active qubit
                  when the other has already been lost prior to this channel.
        """
        self._interface.correlated_loss_channel(addr0, addr1, p)

    def measure(self, addr0: int) -> MeasurementResult:
        """Measure the specified qubit in the Z basis.

        Args:
            addr0: The index of the target qubit.

        Returns:
            The measurement outcome as a ``MeasurementResult``, which is
            ``LOST`` if the qubit has been lost, ``ZERO`` or ``ONE`` otherwise.
        """
        m = self._interface.measure(addr0)
        return MeasurementResult._from_raw(m)

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

    def reset(self, addr0: int) -> None:
        """Reset the specified qubit to the |0> state.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.reset(addr0)

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

    def run_stim_string(self, circuit: str) -> list[MeasurementResult]:
        """Execute a STIM circuit given as a string and return all measurement results.

        .. note::
            This method **mutates** the tableau in place. To run multiple
            independent shots, call :meth:`fork` to obtain a fresh copy before
            each run.

        Parses and runs a STIM-format circuit, applying each instruction to
        this tableau in sequence. Only the squin subset of the STIM instruction
        set is supported; control flow is not supported. Measurements are
        collected in the order they appear in the circuit, following the STIM
        convention: each measurement instruction appends its results
        left-to-right as the qubits are listed, and later instructions append
        after earlier ones.

        Args:
            circuit: A multi-line string containing a STIM circuit.

        Returns:
            A list of ``MeasurementResult`` values, one per measured qubit,
            in circuit order. Each value is ``ZERO``, ``ONE``, or ``LOST``
            (if the qubit had been lost prior to measurement).
        """
        results = self._interface.run_stim_string(circuit)
        return list(map(MeasurementResult._from_raw, results))

    def run_stim_file(self, file_path: str) -> list[MeasurementResult]:
        """Execute a STIM circuit from a file and return all measurement results.

        .. note::
            This method **mutates** the tableau in place. To run multiple
            independent shots, call :meth:`fork` to obtain a fresh copy before
            each run.

        Reads the circuit from ``file_path`` and runs it identically to
        :meth:`run_stim_string`. Only the squin subset of the STIM instruction
        set is supported; control flow is not supported. Measurements are
        collected in the order they appear in the circuit, following the STIM
        convention: each measurement instruction appends its results
        left-to-right as the qubits are listed, and later instructions append
        after earlier ones.

        Args:
            file_path: Path to a ``.stim`` file containing a STIM circuit.

        Returns:
            A list of ``MeasurementResult`` values, one per measured qubit,
            in circuit order. Each value is ``ZERO``, ``ONE``, or ``LOST``
            (if the qubit had been lost prior to measurement).

        Raises:
            pyo3_runtime.PanicException: If the file cannot be read.
        """
        results = self._interface.run_stim_file(file_path)
        return list(map(MeasurementResult._from_raw, results))
