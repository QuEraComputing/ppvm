import math
from dataclasses import InitVar, dataclass, field
from typing import Optional

import ppvm_python_native

from .clifford import (
    CliffordExtensionMixin,
    CliffordMixin,
    NoiseMixin,
    NonCliffordMixin,
)
from .types import GeneralizedTableauInterface


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
        N_interface = math.ceil(self.n_qubits / 8.0)
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

    def measure(self, addr0: int) -> bool:
        """Measure the specified qubit in the Z basis.

        Args:
            addr0: The index of the target qubit.

        Returns:
            The measurement outcome (False = 0, True = 1).
        """
        return self._interface.measure(addr0)

    def reset(self, addr0: int) -> None:
        """Reset the specified qubit to the |0> state.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.reset(addr0)

    def reset_loss_channel(self, addr0: int) -> None:
        """Reset the loss channel for the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.reset_loss_channel(addr0)
