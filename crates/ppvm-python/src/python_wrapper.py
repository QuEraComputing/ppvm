import math
from typing import Sequence, Any

import ppvm_python


class PauliSum:
    n_qubits: int
    min_abs_coeff: float
    max_pauli_weight: int | None
    terms: Sequence[str]
    coefficients: Sequence[float]

    _interface: Any

    def __init__(
        self,
        n_qubits: int,
        terms: Sequence[str],
        min_abs_coeff: float = 1e-10,
        max_pauli_weight: int | None = None,
        coefficients: Sequence[float] = (),
    ):
        self.n_qubits = n_qubits
        self.min_abs_coeff = min_abs_coeff
        self.max_pauli_weight = max_pauli_weight
        self.terms = terms
        self.coefficients = coefficients
        self._interface = self._init_ppvm_interface(terms, coefficients)

    def _init_ppvm_interface(self, terms: Sequence[str], coefficients: Sequence[float]):
        n_qubits = self.n_qubits

        # number of bytes we need
        N = math.ceil(n_qubits / 8.0)

        # number of bytes we have
        possible_interfaces = range(15)
        N_interface = next(n for n in possible_interfaces if 2**n > N)

        interface = getattr(ppvm_python, f"PauliSumIndexMapFxHash{N_interface}")

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

    def two_qubit_pauli_error(self, addr0: int, addr1: int, p: Sequence[float]) -> None:
        self._interface.two_qubit_pauli_error(addr0, addr1, p)
