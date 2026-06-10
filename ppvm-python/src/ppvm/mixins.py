# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

from collections.abc import Sequence
from typing import Any

# Gate vocabulary shared across backends comes in two flavours:
#
# * The plain mixins below (``CliffordMixin``, ``RotationsMixin``, ...) forward
#   straight to the native interface with no ``truncate`` kwarg. Truncation is
#   not a per-gate concept for the generalized stabilizer tableau — that
#   representation is exact — so ``GeneralizedTableau`` uses these directly.
# * The ``Truncating*`` variants further down add an optional
#   ``truncate: bool = True`` kwarg and are used by ``PauliSum``. When ``True``
#   (the default), the configured truncation strategy fires immediately after
#   the gate — historical behaviour. Pass ``truncate=False`` to defer the cut;
#   the user is then responsible for calling :meth:`PauliSum.truncate`
#   explicitly when the next truncation point is reached. This is the supported
#   way to chain commuting gates (e.g. ``rxx + ryy`` on the same edge for a
#   U(1)-conserving exchange step) without losing a conserved-charge component
#   to intermediate truncation.


class CliffordMixin:
    """Clifford gates without per-gate truncation control."""

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
    """Rotation gates without per-gate truncation control."""

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

    def r(self, addr0: int, axis_angle: float, theta: float) -> None:
        """Apply a rotation about an axis in the X-Y plane to the specified qubit.

        ```math
        R(\\phi, \\theta) = e^{-i \\frac{\\theta}{2} (\\cos\\phi\\, X + \\sin\\phi\\, Y)}
            = R_Z(\\phi) R_X(\\theta) R_Z(-\\phi)
        ```

        Args:
            addr0: The index of the target qubit.
            axis_angle: The angle ``φ`` (in radians) of the rotation axis
                within the X-Y plane, measured from the X-axis.
            theta: The rotation angle in radians.
        """
        self._interface.r(addr0, axis_angle, theta)

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
    """Additional Clifford gates without per-gate truncation control."""

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
    """Noise channels without per-gate truncation control."""

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

    def two_qubit_pauli_error(
        self,
        addr0: int,
        addr1: int,
        p: Sequence[float],
    ) -> None:
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
    """Loss channels without per-gate truncation control."""

    _interface: Any

    def loss_channel(self, addr0: int, p: float) -> None:
        """Apply a loss channel to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            p: The loss probability.
        """
        self._interface.loss_channel(addr0, p)

    def correlated_loss_channel(
        self,
        addr0: int,
        addr1: int,
        p: Sequence[float],
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


# Truncating variants used by ``PauliSum``: same gates, plus a per-gate
# ``truncate`` kwarg. Each subclasses its plain counterpart, so members that
# have no truncation concept (``cy``, ``two_qubit_pauli_error_probabilities``)
# are inherited unchanged and only the branching gates are overridden.


class TruncatingCliffordMixin(CliffordMixin):
    """Clifford gates with a per-gate ``truncate`` kwarg (used by PauliSum)."""

    def x(self, addr0: int, *, truncate: bool = True) -> None:
        """Apply a Pauli X gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            truncate: If ``True`` (default), run the configured truncation
                strategy after the gate; if ``False``, leave the map
                untruncated so the next gate sees the full unpruned state.
        """
        self._interface.x(addr0, truncate=truncate)

    def y(self, addr0: int, *, truncate: bool = True) -> None:
        """Apply a Pauli Y gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            truncate: See :meth:`x`.
        """
        self._interface.y(addr0, truncate=truncate)

    def z(self, addr0: int, *, truncate: bool = True) -> None:
        """Apply a Pauli Z gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            truncate: See :meth:`x`.
        """
        self._interface.z(addr0, truncate=truncate)

    def h(self, addr0: int, *, truncate: bool = True) -> None:
        """Apply a Hadamard gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            truncate: See :meth:`x`.
        """
        self._interface.h(addr0, truncate=truncate)

    def s(self, addr0: int, *, truncate: bool = True) -> None:
        """Apply an S gate (sqrt(Z)) to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            truncate: See :meth:`x`.
        """
        self._interface.s(addr0, truncate=truncate)

    def cnot(self, addr0: int, addr1: int, *, truncate: bool = True) -> None:
        """Apply a CNOT (controlled-X) gate.

        Args:
            addr0: The index of the control qubit.
            addr1: The index of the target qubit.
            truncate: See :meth:`x`.
        """
        self._interface.cnot(addr0, addr1, truncate=truncate)

    def cz(self, addr0: int, addr1: int, *, truncate: bool = True) -> None:
        """Apply a CZ (controlled-Z) gate.

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            truncate: See :meth:`x`.
        """
        self._interface.cz(addr0, addr1, truncate=truncate)


class TruncatingRotationsMixin(RotationsMixin):
    """Rotation gates with a per-gate ``truncate`` kwarg (used by PauliSum)."""

    def rx(self, addr0: int, theta: float, *, truncate: bool = True) -> None:
        """Apply an RX rotation gate to the specified qubit.

        See :meth:`RotationsMixin.rx` for the gate definition.

        Args:
            addr0: The index of the target qubit.
            theta: The rotation angle in radians.
            truncate: If ``True`` (default), run the configured truncation
                strategy after the gate; if ``False``, defer it.
        """
        self._interface.rx(addr0, theta, truncate=truncate)

    def ry(self, addr0: int, theta: float, *, truncate: bool = True) -> None:
        """Apply an RY rotation gate to the specified qubit.

        See :meth:`RotationsMixin.ry` for the gate definition.

        Args:
            addr0: The index of the target qubit.
            theta: The rotation angle in radians.
            truncate: See :meth:`rx`.
        """
        self._interface.ry(addr0, theta, truncate=truncate)

    def rz(self, addr0: int, theta: float, *, truncate: bool = True) -> None:
        """Apply an RZ rotation gate to the specified qubit.

        See :meth:`RotationsMixin.rz` for the gate definition.

        Args:
            addr0: The index of the target qubit.
            theta: The rotation angle in radians.
            truncate: See :meth:`rx`.
        """
        self._interface.rz(addr0, theta, truncate=truncate)

    def r(
        self, addr0: int, axis_angle: float, theta: float, *, truncate: bool = True
    ) -> None:
        """Apply a rotation about an axis in the X-Y plane to the specified qubit.

        See :meth:`RotationsMixin.r` for the gate definition.

        Args:
            addr0: The index of the target qubit.
            axis_angle: The angle ``φ`` (in radians) of the rotation axis
                within the X-Y plane, measured from the X-axis.
            theta: The rotation angle in radians.
            truncate: See :meth:`rx`.
        """
        self._interface.r(addr0, axis_angle, theta, truncate=truncate)

    # Two qubit rotations
    def rxx(self, addr0: int, addr1: int, theta: float, *, truncate: bool = True) -> None:
        """Apply an RXX (Ising XX) rotation gate to two qubits.

        See :meth:`RotationsMixin.rxx` for the gate definition.

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            theta: The rotation angle in radians.
            truncate: If ``True`` (default), run the configured truncation
                strategy after the gate. Set to ``False`` to compose a
                U(1)-conserving step like ``rxx + ryy`` on the same edge
                without dropping the conserved-charge component between
                them — then call :meth:`PauliSum.truncate` once at the
                end of the composition.
        """
        self._interface.rxx(addr0, addr1, theta, truncate=truncate)

    def ryy(self, addr0: int, addr1: int, theta: float, *, truncate: bool = True) -> None:
        """Apply an RYY (Ising YY) rotation gate to two qubits.

        See :meth:`RotationsMixin.ryy` for the gate definition.

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            theta: The rotation angle in radians.
            truncate: See :meth:`rxx`.
        """
        self._interface.ryy(addr0, addr1, theta, truncate=truncate)

    def rzz(self, addr0: int, addr1: int, theta: float, *, truncate: bool = True) -> None:
        """Apply an RZZ (Ising ZZ) rotation gate to two qubits.

        See :meth:`RotationsMixin.rzz` for the gate definition.

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            theta: The rotation angle in radians.
            truncate: See :meth:`rxx`.
        """
        self._interface.rzz(addr0, addr1, theta, truncate=truncate)


class TruncatingCliffordExtensionMixin(CliffordExtensionMixin):
    """Additional Clifford gates with a per-gate ``truncate`` kwarg.

    ``cy`` is inherited unchanged from :class:`CliffordExtensionMixin` (it does
    not branch the Pauli sum, so it takes no ``truncate`` kwarg).
    """

    def s_adj(self, addr0: int, *, truncate: bool = True) -> None:
        """Apply an S adjoint gate (sqrt(Z) dagger) to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            truncate: See :meth:`TruncatingCliffordMixin.x`.
        """
        self._interface.s_adj(addr0, truncate=truncate)

    def sqrt_x(self, addr0: int, *, truncate: bool = True) -> None:
        """Apply a sqrt(X) gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            truncate: See :meth:`TruncatingCliffordMixin.x`.
        """
        self._interface.sqrt_x(addr0, truncate=truncate)

    def sqrt_y(self, addr0: int, *, truncate: bool = True) -> None:
        """Apply a sqrt(Y) gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            truncate: See :meth:`TruncatingCliffordMixin.x`.
        """
        self._interface.sqrt_y(addr0, truncate=truncate)

    def sqrt_x_adj(self, addr0: int, *, truncate: bool = True) -> None:
        """Apply a sqrt(X) adjoint gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            truncate: See :meth:`TruncatingCliffordMixin.x`.
        """
        self._interface.sqrt_x_adj(addr0, truncate=truncate)

    def sqrt_y_adj(self, addr0: int, *, truncate: bool = True) -> None:
        """Apply a sqrt(Y) adjoint gate to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            truncate: See :meth:`TruncatingCliffordMixin.x`.
        """
        self._interface.sqrt_y_adj(addr0, truncate=truncate)


class TruncatingNoiseMixin(NoiseMixin):
    """Noise channels with a per-gate ``truncate`` kwarg (used by PauliSum).

    ``two_qubit_pauli_error_probabilities`` is inherited unchanged from
    :class:`NoiseMixin`.
    """

    def pauli_error(
        self, addr0: int, p: Sequence[float], *, truncate: bool = True
    ) -> None:
        """Apply a single-qubit Pauli error channel.

        Args:
            addr0: The index of the target qubit.
            p: Error probabilities [p_x, p_y, p_z] for X, Y, Z errors.
                The identity probability is implicitly 1 - sum(p).
            truncate: If ``True`` (default), run the configured truncation
                strategy after the channel; if ``False``, defer it.
        """
        self._interface.pauli_error(addr0, p, truncate=truncate)

    def two_qubit_pauli_error(
        self,
        addr0: int,
        addr1: int,
        p: Sequence[float],
        *,
        truncate: bool = True,
    ) -> None:
        """Apply a two-qubit Pauli error channel.

        Args:
            addr0: The index of the first qubit.
            addr1: The index of the second qubit.
            p: Error probabilities for the 15 non-identity two-qubit Pauli
                operators. Use two_qubit_pauli_error_probabilities to convert
                from a dictionary format.
            truncate: See :meth:`pauli_error`.
        """
        self._interface.two_qubit_pauli_error(addr0, addr1, p, truncate=truncate)

    def depolarize(self, addr0: int, p: float, *, truncate: bool = True) -> None:
        """Apply a depolarizing channel to the specified qubit.

        Args:
            addr0: The index of the target qubit.
            p: The depolarizing probability.
            truncate: See :meth:`pauli_error`.
        """
        self._interface.depolarize(addr0, p, truncate=truncate)

    def depolarize2(
        self, addr0: int, addr1: int, p: float, *, truncate: bool = True
    ) -> None:
        """Apply a two-qubit depolarizing channel to the specified qubits.

        Args:
            addr0: The index of the first target qubit.
            addr1: The index of the second target qubit.
            p: The depolarizing probability.
            truncate: See :meth:`pauli_error`.
        """
        self._interface.depolarize2(addr0, addr1, p, truncate=truncate)
