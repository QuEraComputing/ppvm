from typing import Any, Sequence


# TODO: also use this in PauliSum
class CliffordMixin:
    _interface: Any

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


class RotationsMixin:
    _interface: Any

    # Rotations
    def rx(self, addr0: int, theta: float) -> None:
        """Apply an RX rotation gate to the specified qubit.

        ```math
        R_X(\\theta) = e^{-i \\frac{\\theta}{2} X} = \\cos\\frac{\\theta}{2} I - i \\sin\\frac{\\theta}{2} X
        ```

        Args:
            addr0: The index of the target qubit.
            theta: The rotation angle in radians.
        """
        self._interface.rx(addr0, theta)

    def ry(self, addr0: int, theta: float) -> None:
        """Apply an RY rotation gate to the specified qubit.

        ```math
        R_Y(\\theta) = e^{-i \\frac{\\theta}{2} Y} = \\cos\\frac{\\theta}{2} I - i \\sin\\frac{\\theta}{2} Y
        ```

        Args:
            addr0: The index of the target qubit.
            theta: The rotation angle in radians.
        """
        self._interface.ry(addr0, theta)

    def rz(self, addr0: int, theta: float) -> None:
        """Apply an RZ rotation gate to the specified qubit.

        ```math
        R_Z(\\theta) = e^{-i \\frac{\\theta}{2} Z} = \\cos\\frac{\\theta}{2} I - i \\sin\\frac{\\theta}{2} Z
        ```

        Args:
            addr0: The index of the target qubit.
            theta: The rotation angle in radians.
        """
        self._interface.rz(addr0, theta)

    # Two qubit rotations
    def rxx(self, addr0: int, addr1: int, theta: float) -> None:
        """Apply an RXX (Ising XX) rotation gate to two qubits.

        ```math
        R_{XX}(\\theta) = e^{-i \\frac{\\theta}{2} X \\otimes X}
        ```

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            theta: The rotation angle in radians.
        """
        self._interface.rxx(addr0, addr1, theta)

    def ryy(self, addr0: int, addr1: int, theta: float) -> None:
        """Apply an RYY (Ising YY) rotation gate to two qubits.

        ```math
        R_{YY}(\\theta) = e^{-i \\frac{\\theta}{2} Y \\otimes Y}
        ```

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            theta: The rotation angle in radians.
        """
        self._interface.ryy(addr0, addr1, theta)

    def rzz(self, addr0: int, addr1: int, theta: float) -> None:
        """Apply an RZZ (Ising ZZ) rotation gate to two qubits.

        ```math
        R_{ZZ}(\\theta) = e^{-i \\frac{\\theta}{2} Z \\otimes Z}
        ```

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            theta: The rotation angle in radians.
        """
        self._interface.rzz(addr0, addr1, theta)


class CliffordExtensionMixin:
    _interface: Any

    def s_adj(self, addr0: int) -> None:
        """Apply an S adjoint gate (sqrt(Z) dagger) to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.s_adj(addr0)

    def sqrt_x(self, addr0: int) -> None:
        """Apply a sqrt(X) gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.sqrt_x(addr0)

    def sqrt_y(self, addr0: int) -> None:
        """Apply a sqrt(Y) gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.sqrt_y(addr0)

    def sqrt_x_adj(self, addr0: int) -> None:
        """Apply a sqrt(X) adjoint gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.sqrt_x_adj(addr0)

    def sqrt_y_adj(self, addr0: int) -> None:
        """Apply a sqrt(Y) adjoint gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
        """
        self._interface.sqrt_y_adj(addr0)

    def cy(self, addr0: int, addr1: int) -> None:
        """Apply a controlled-Y gate.

        Args:
            addr0: The index of the control qubit.
            addr1: The index of the target qubit.
        """
        self._interface.cy(addr0, addr1)


class NoiseMixin:
    _interface: Any

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


class LossMixin:
    _interface: Any

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
