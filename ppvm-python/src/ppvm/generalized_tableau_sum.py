from dataclasses import InitVar, dataclass, field
from typing import cast

import ppvm_python_native

from .generalized_tableau import MeasurementResult
from .mixins import (
    CliffordExtensionMixin,
    CliffordMixin,
    LossMixin,
    NoiseMixin,
    RotationsMixin,
)
from .types import GeneralizedTableauSumInterface, TableauSumSamplerInterface

# Indexed by integer outcome value (0/1/2) to reuse the singleton enum members.
# This is much faster than calling ``MeasurementResult(i)`` per element: the
# IntEnum constructor dominates large shot batches, while a tuple index just
# bumps a refcount.
_BY_VALUE = (MeasurementResult.ZERO, MeasurementResult.ONE, MeasurementResult.LOST)


@dataclass(frozen=True)
class GeneralizedTableauSum(
    CliffordMixin,
    CliffordExtensionMixin,
    RotationsMixin,
    NoiseMixin,
    LossMixin,
):
    """A sum over generalized tableaus, representing a density matrix that is a
    weighted sum over pure state projectors (classical mixture). Each projector
    is represented by a generalized tableau. The respective weights are the
    probabilities with which the system is in the corresponding state.
    """

    n_qubits: int
    min_abs_coeff: float = 1e-10
    sum_cutoff: float = 1e-8
    seed: InitVar[int | None] = None

    _interface: GeneralizedTableauSumInterface = field(init=False, repr=False)

    def __post_init__(self, seed: int | None):
        N_interface = (self.n_qubits + 63) // 64
        object.__setattr__(
            self,
            "_interface",
            getattr(ppvm_python_native, f"GeneralizedTableauSum{N_interface}")(
                self.n_qubits, self.min_abs_coeff, self.sum_cutoff, seed
            ),
        )

    def __copy__(self) -> "GeneralizedTableauSum":
        """Return a copy of this tableau sum, including its RNG state.

        Both the original and the copy will produce identical random sequences
        from this point forward.
        """
        copied = GeneralizedTableauSum(self.n_qubits, self.min_abs_coeff, self.sum_cutoff)
        object.__setattr__(copied, "_interface", self._interface.__copy__())
        return copied

    def __deepcopy__(self, memo: dict) -> "GeneralizedTableauSum":
        """Return a deep copy of this tableau sum, including its RNG state.

        Both the original and the copy will produce identical random sequences
        from this point forward.
        """
        copied = GeneralizedTableauSum(self.n_qubits, self.min_abs_coeff, self.sum_cutoff)
        object.__setattr__(copied, "_interface", self._interface.__deepcopy__(memo))
        return copied

    def __str__(self) -> str:
        """Return a human-readable representation of the tableau state."""
        return self._interface.__str__()

    def __len__(self) -> int:
        """Return the number of branches currently in the sum."""
        return len(self._interface)

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

    def measure(self, addr0: int) -> dict[MeasurementResult, float]:
        """Branch on a mid-circuit measurement and return probabilities for outcomes

        Args:
            addr0: The index of the target qubit.

        Returns:
            
        """
        p0, p1, pl = self._interface.measure(addr0)
        return {
            MeasurementResult.ZERO: p0,
            MeasurementResult.ONE: p1,
            MeasurementResult.LOST: pl,
        }

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
    
    def sampler(self) -> "TableauSumSampler":
        base_sampler = self._interface.sampler()
        return TableauSumSampler(
            cast(TableauSumSamplerInterface, base_sampler)
        )



@dataclass(frozen=True)
class TableauSumSampler:
    """Sample a `GeneralizedTableauSum`.

    Construct via `GeneralizedTableauSum.sampler()`.
    """

    _interface: TableauSumSamplerInterface = field(repr=False)

    def sample(self) -> list[MeasurementResult]:
        return [_BY_VALUE[i] for i in self._interface.sample()]

    def raw_sample(self) -> list[int]:
        return self._interface.sample()

    def sample_shots(self, num_shots: int) -> list[list[MeasurementResult]]:
        raw_samples = self._interface.sample_shots(num_shots=num_shots)
        return [[_BY_VALUE[i] for i in ints] for ints in raw_samples]

    def raw_shots(self, num_shots: int) -> list[list[int]]:
        return self._interface.sample_shots(num_shots=num_shots)
