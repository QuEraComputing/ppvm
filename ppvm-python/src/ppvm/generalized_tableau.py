import math
from dataclasses import dataclass, field

import ppvm_python_native

from .clifford import CliffordExtensionMixin, CliffordMixin, NoiseMixin
from .types import GeneralizedTableauInterface


@dataclass(frozen=True)
class GeneralizedTableau(CliffordMixin, CliffordExtensionMixin, NoiseMixin):
    n_qubits: int
    min_abs_coeff: float = 1e-10

    _interface: GeneralizedTableauInterface = field(init=False, repr=False)

    def __post_init__(self):
        N_interface = math.ceil(self.n_qubits / 8.0)
        object.__setattr__(
            self,
            "_interface",
            getattr(ppvm_python_native, f"GeneralizedTableau{N_interface}")(
                self.n_qubits, self.min_abs_coeff
            ),
        )

    def __str__(self) -> str:
        return self._interface.__str__()

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
