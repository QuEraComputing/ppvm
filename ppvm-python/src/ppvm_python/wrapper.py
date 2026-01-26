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
    """A weighted sum of Pauli strings for quantum simulation.

    PauliSum represents a linear combination of Pauli operators, commonly used
    to represent quantum observables or Hamiltonians. It provides methods for
    applying quantum gates (Clifford operations and rotations) and computing
    expectation values via the trace operation.

    Attributes:
        terms: Pauli strings, each containing only 'I', 'X', 'Y', 'Z' characters.
            All terms must have the same length (number of qubits).
        n_qubits: Number of qubits. If None, inferred from the length of the
            first term.
        coefficients: Coefficients for each Pauli term. If empty, all terms
            are assigned coefficient 1.0.
        min_abs_coeff: Minimum absolute coefficient value. Terms with smaller
            coefficients are dropped for efficiency.
        max_pauli_weight: Maximum number of non-identity Paulis allowed per term.
            If None, uses the backend default.

    Note:
        Gates must be applied in reverse circuit order. This is because PauliSum
        evolves observables in the Heisenberg picture rather than states in the
        Schrödinger picture.

    Example:
        Basic usage with a simple Pauli sum:

        ```python
        # Create a simple Pauli sum: 0.5 * ZZ + 0.3 * XI
        ps = PauliSum(terms=["ZZ", "XI"], coefficients=[0.5, 0.3])
        # For a circuit: RZ(0.5) on qubit 1, then H on qubit 0
        # Apply in reverse order:
        ps.rz(1, 0.5)
        ps.h(0)
        # Compute overlap with |0...0> state
        result = ps.overlap_with_zero()
        ```

        Simulating a 3-qubit GHZ state preparation circuit:

        ```python
        # Start with ZZZ observable (measures all qubits in Z basis)
        ps = PauliSum.from_str("ZZZ")
        # GHZ circuit: H(0), CNOT(0,1), CNOT(1,2)
        # Apply in reverse order:
        ps.cnot(1, 2)
        ps.cnot(0, 1)
        ps.h(0)
        # Expectation value of ZZZ for GHZ state is 0
        result = ps.overlap_with_zero()
        ```
    """

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
        """Create a PauliSum from a single Pauli string with coefficient 1.0.

        Args:
            s: A Pauli string containing only 'I', 'X', 'Y', 'Z' characters.

        Returns:
            A PauliSum with a single term and coefficient 1.0.

        Raises:
            ValueError: If the string contains invalid characters.
        """
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
        """Compute the overlap with the all-zeros computational basis state.

        Returns:
            The expectation value of the Pauli sum with respect to |0...0>.
        """
        return self._interface.trace("Z?*")

    def trace(self, pattern: str) -> float:
        """Compute the trace using a pattern string.

        Args:
            pattern: A pattern specifying which terms to include in the trace.
                Use 'Z' to project onto |0>, '?' for any single character,
                and '*' to match zero or more characters.

        Returns:
            The trace result.
        """
        return self._interface.trace(pattern)

    # Clifford operations
    def x(self, addr0: int) -> None:
        """Apply a Pauli X gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.x(addr0)

    def y(self, addr0: int) -> None:
        """Apply a Pauli Y gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.y(addr0)

    def z(self, addr0: int) -> None:
        """Apply a Pauli Z gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.z(addr0)

    def h(self, addr0: int) -> None:
        """Apply a Hadamard gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.h(addr0)

    def s(self, addr0: int) -> None:
        """Apply an S gate (sqrt(Z)) to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.s(addr0)

    def cnot(self, addr0: int, addr1: int) -> None:
        """Apply a CNOT (controlled-X) gate.

        Args:
            addr0: The index of the control qubit.
            addr1: The index of the target qubit.
        """
        self._interface.cnot(addr0, addr1)

    def cz(self, addr0: int, addr1: int) -> None:
        """Apply a CZ (controlled-Z) gate.

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
        """
        self._interface.cz(addr0, addr1)

    # Rotations
    def rx(self, addr0: int, theta: float) -> None:
        """Apply an RX rotation gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            theta: The rotation angle in radians.
        """
        self._interface.rx(addr0, theta)

    def ry(self, addr0: int, theta: float) -> None:
        """Apply an RY rotation gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            theta: The rotation angle in radians.
        """
        self._interface.ry(addr0, theta)

    def rz(self, addr0: int, theta: float) -> None:
        """Apply an RZ rotation gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            theta: The rotation angle in radians.
        """
        self._interface.rz(addr0, theta)

    # Two qubit rotations
    def rxx(self, addr0: int, addr1: int, theta: float) -> None:
        """Apply an RXX (Ising XX) rotation gate to two qubits.

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            theta: The rotation angle in radians.
        """
        self._interface.rxx(addr0, addr1, theta)

    def ryy(self, addr0: int, addr1: int, theta: float) -> None:
        """Apply an RYY (Ising YY) rotation gate to two qubits.

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            theta: The rotation angle in radians.
        """
        self._interface.ryy(addr0, addr1, theta)

    def rzz(self, addr0: int, addr1: int, theta: float) -> None:
        """Apply an RZZ (Ising ZZ) rotation gate to two qubits.

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            theta: The rotation angle in radians.
        """
        self._interface.rzz(addr0, addr1, theta)

    # Noise operations
    def pauli_error(self, addr0: int, p: Sequence[float]) -> None:
        """Apply a single-qubit Pauli error channel.

        Args:
            addr0: The index of the target qubit.
            p: Error probabilities [p_x, p_y, p_z] for X, Y, Z errors.
                The identity probability is implicitly 1 - sum(p).
        """
        self._interface.pauli_error(addr0, p)

    @staticmethod
    def two_qubit_pauli_error_probabilities(
        error_probabilities: dict[str, float],
    ) -> list[float]:
        """Convert a dictionary of two-qubit Pauli error probabilities to a list.

        Convenience method to convert a dictionary mapping two-qubit Pauli
        strings (e.g., "IX", "ZZ") to their probabilities into the ordered
        list format required by two_qubit_pauli_error.

        Args:
            error_probabilities: Dictionary mapping two-qubit Pauli strings
                to their error probabilities. Missing keys default to 0.0.

        Returns:
            A list of 15 probabilities in the canonical order (excludes "II").
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
        """Apply a two-qubit Pauli error channel.

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            p: Error probabilities for the 15 non-identity two-qubit Pauli
                operators. Use two_qubit_pauli_error_probabilities to convert
                from a dictionary format.
        """
        self._interface.two_qubit_pauli_error(addr0, addr1, p)
