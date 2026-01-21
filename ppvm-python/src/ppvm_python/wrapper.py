import math
from dataclasses import dataclass, field
from typing import Sequence, Union

import ppvm_python_native

T = Union[
    ppvm_python_native.PauliSumIndexMapFxHash0,
    ppvm_python_native.PauliSumIndexMapFxHash1,
    ppvm_python_native.PauliSumIndexMapFxHash2,
    ppvm_python_native.PauliSumIndexMapFxHash3,
    ppvm_python_native.PauliSumIndexMapFxHash4,
    ppvm_python_native.PauliSumIndexMapFxHash5,
    ppvm_python_native.PauliSumIndexMapFxHash6,
    ppvm_python_native.PauliSumIndexMapFxHash7,
    ppvm_python_native.PauliSumIndexMapFxHash8,
    ppvm_python_native.PauliSumIndexMapFxHash9,
    ppvm_python_native.PauliSumIndexMapFxHash10,
    ppvm_python_native.PauliSumIndexMapFxHash11,
    ppvm_python_native.PauliSumIndexMapFxHash12,
    ppvm_python_native.PauliSumIndexMapFxHash13,
    ppvm_python_native.PauliSumIndexMapFxHash14,
    ppvm_python_native.PauliSumIndexMapFxHash15,
]


@dataclass(frozen=True)
class PauliSum:
    terms: Sequence[str]
    n_qubits: int | None = None
    coefficients: Sequence[float] = ()
    min_abs_coeff: float = 1e-10
    max_pauli_weight: int | None = None

    _interface: T = field(init=False, repr=False)

    def __post_init__(self):
        object.__setattr__(
            self,
            "_interface",
            self._init_ppvm_interface(),
        )

    def _init_ppvm_interface(
        self,
    ):

        n_qubits = self.n_qubits
        terms = self.terms
        coefficients = self.coefficients

        if not terms:
            raise ValueError(
                "At least one term must be provided to initialize PauliSum."
            )

        if n_qubits is None:
            n_qubits = len(terms[0])

        for term in terms:
            if len(term) != n_qubits:
                raise ValueError(
                    "All terms must have the same length! Expected length "
                    f"{n_qubits}, but got term of length {len(term)}: {term!r}"
                )

        # number of bytes we need
        N = math.ceil(n_qubits / 8.0)

        # number of bytes we have
        possible_interfaces = range(15)
        N_interface = next(n for n in possible_interfaces if 2**n > N)

        interface = getattr(ppvm_python_native, f"PauliSumIndexMapFxHash{N_interface}")

        if terms and not coefficients:
            coefficients = (1.0,) * len(terms)

        if self.max_pauli_weight is None:
            # NOTE: let rust handle the default setting for max_pauli_weight
            return interface(
                n_qubits,
                min_abs_coeff=self.min_abs_coeff,
                terms=terms,
                coefficients=coefficients,
            )
        else:
            return interface(
                n_qubits,
                min_abs_coeff=self.min_abs_coeff,
                max_pauli_weight=self.max_pauli_weight,
                terms=terms,
                coefficients=coefficients,
            )

    @staticmethod
    def from_str(s: str) -> "PauliSum":
        s = s.strip()
        # Validate the string: must only contain I, X, Y, Z
        allowed = set("IXYZ")
        if not set(s).issubset(allowed):
            raise ValueError(
                f"Invalid Pauli string: {s!r}. Only 'I', 'X', 'Y', 'Z' are allowed."
            )
        n_qubits = len(s)
        terms = [s]
        coefficients = [1.0]
        return PauliSum(n_qubits=n_qubits, terms=terms, coefficients=coefficients)

    def __str__(self) -> str:
        return self._interface.__str__()

    # Getting results
    def overlap_with_zero(self) -> float:
        return self._interface.trace("Z?*")

    def trace(self, pattern: str) -> float:
        return self._interface.trace(pattern)

    # Clifford operations
    def x(self, addr0: int) -> None:
        self._interface.x(addr0)

    def y(self, addr0: int) -> None:
        self._interface.y(addr0)

    def z(self, addr0: int) -> None:
        self._interface.z(addr0)

    def h(self, addr0: int) -> None:
        self._interface.h(addr0)

    def s(self, addr0: int) -> None:
        self._interface.s(addr0)

    def cnot(self, addr0: int, addr1: int) -> None:
        self._interface.cnot(addr0, addr1)

    def cz(self, addr0: int, addr1: int) -> None:
        self._interface.cz(addr0, addr1)

    # Rotations
    def rx(self, addr0: int, theta: float) -> None:
        self._interface.rx(addr0, theta)

    def ry(self, addr0: int, theta: float) -> None:
        self._interface.ry(addr0, theta)

    def rz(self, addr0: int, theta: float) -> None:
        self._interface.rz(addr0, theta)

    # Two qubit rotations
    def rxx(self, addr0: int, addr1: int, theta: float) -> None:
        self._interface.rxx(addr0, addr1, theta)

    def ryy(self, addr0: int, addr1: int, theta: float) -> None:
        self._interface.ryy(addr0, addr1, theta)

    def rzz(self, addr0: int, addr1: int, theta: float) -> None:
        self._interface.rzz(addr0, addr1, theta)

    # Noise operations
    def pauli_error(self, addr0: int, p: Sequence[float]) -> None:
        self._interface.pauli_error(addr0, p)

    @staticmethod
    def two_qubit_pauli_error_probabilities(
        error_probabilities: dict[str, float],
    ) -> list[float]:
        """Convert a dictionary of two-qubit Pauli error probabilities to a list.

        Convenience method to convert a dictionary mapping two-qubit Pauli strings ensuring correct order.
        """
        keys = (
            "IX",
            "IT",
            "IZ",
            "XI",
            "XX",
            "XY",
            "XZ",
            "YI",
            "YX",
            "YY",
            "YZ",
            "ZI",
            "ZX",
            "ZY",
            "ZZ",
        )
        return [error_probabilities.get(key, 0.0) for key in keys]

    def two_qubit_pauli_error(self, addr0: int, addr1: int, p: Sequence[float]) -> None:
        self._interface.two_qubit_pauli_error(addr0, addr1, p)
