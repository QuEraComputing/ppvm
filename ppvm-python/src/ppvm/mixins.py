# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

from collections.abc import Iterable, Sequence
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
#   the user is then responsible for calling `PauliSum.truncate`
#   explicitly when the next truncation point is reached. This is the supported
#   way to chain commuting gates (e.g. ``rxx + ryy`` on the same edge for a
#   U(1)-conserving exchange step) without losing a conserved-charge component
#   to intermediate truncation.
#
# Gate methods broadcast over targets and accept them in two interchangeable
# forms, mirroring stim's two surfaces:
#
# * Variadic indices, like stim's ``TableauSimulator`` — ``x(0, 1, 2)``,
#   ``cnot(0, 1, 2, 3)``.
# * A single sequence (list / tuple / range / numpy array), like stim's
#   ``Circuit.append`` — ``x([0, 1, 2])``, ``cnot(np.array([0, 1, 2, 3]))``.
#
# Single-qubit gates apply to each index; two-qubit gates consume consecutive
# pairs (the target count must be even). Trailing scalar parameters (``theta``
# for rotations, ``p`` for noise) may be passed either as keywords or as the
# final positional argument. (``measure`` stays scalar and never broadcasts —
# use ``measure_many`` for multiple qubits — matching stim, which keeps
# ``measure`` single to avoid the ``if sim.measure(q):`` returns-a-list trap.)


def _is_sequence(obj: Any) -> bool:
    """True for a target *collection* (list / tuple / range / ndarray, ...).

    Any iterable counts except ``str``/``bytes`` (iterable but never targets).
    A bare ``int`` — including a numpy integer scalar, which is not iterable —
    is not a sequence, so it falls through to the variadic path.
    """
    # Concrete-type fast paths first: ``list``/``tuple`` (the overwhelmingly
    # common splatted form) and bare ``int`` short-circuit before the slow ABC
    # ``isinstance(obj, Iterable)`` dispatch, which is run only for the rare
    # range / ndarray / generator cases. ``str``/``bytes`` are iterable but are
    # never targets, so they report False.
    if isinstance(obj, (list, tuple)):
        return True
    if isinstance(obj, (int, str, bytes)):
        return False
    return isinstance(obj, Iterable)


def _normalize_targets(args: tuple[Any, ...]) -> Sequence[int]:
    """Resolve gate targets passed either as variadic indices (``x(0, 1, 2)``)
    or as a single sequence (``x([0, 1, 2])``, ``x(np.array([0, 1, 2]))``).

    Returns the targets as-is — the single sequence, or the variadic ``args``
    tuple. The native layer extracts a ``Vec<usize>`` directly (PyO3 handles
    Python ints, numpy integer scalars, ranges and ndarrays), so there is no
    need to rebuild the list with a per-element ``int()`` on the hot path."""
    if len(args) == 1 and _is_sequence(args[0]):
        return args[0]
    return args


def _split_targets_parameter(
    args: tuple[Any, ...],
    value: Any | None,
    name: str,
) -> tuple[Sequence[int], Any]:
    """Split ``(*targets, value)`` accepting ``value=...`` and a single leading
    sequence of targets (``([0, 1, 2], theta)`` as well as ``(0, 1, 2, theta)``).

    Targets are returned as-is (sequence or tuple slice); the native layer does
    the ``Vec<usize>`` extraction, so no per-element ``int()`` rebuild is needed."""
    if args and _is_sequence(args[0]):
        targets, rest = args[0], args[1:]
    elif value is None:
        if not args:
            raise TypeError(f"missing required argument: {name!r}")
        targets, rest = args[:-1], args[-1:]
    else:
        targets, rest = args, ()
    if value is None:
        if not rest:
            raise TypeError(f"missing required argument: {name!r}")
        value = rest[-1]
    return targets, value


def _split_targets_parameter_truncate(
    args: tuple[Any, ...],
    value: Any | None,
    name: str,
    truncate: bool,
) -> tuple[Sequence[int], Any, bool]:
    """Split ``(*targets, value[, truncate])`` for PauliSum methods, also
    accepting a single leading sequence of targets.

    Targets are returned as-is (sequence or tuple/list slice); the native layer
    does the ``Vec<usize>`` extraction, so no per-element ``int()`` is needed."""
    if args and _is_sequence(args[0]):
        targets = args[0]
        rest = list(args[1:])
        if rest and isinstance(rest[-1], bool):
            truncate = rest.pop()
        if value is None:
            if not rest:
                raise TypeError(f"missing required argument: {name!r}")
            value = rest[-1]
        return targets, value, truncate
    if value is None:
        if not args:
            raise TypeError(f"missing required argument: {name!r}")
        args_list = list(args)
        if len(args_list) >= 2 and isinstance(args_list[-1], bool):
            truncate = args_list.pop()
        return args_list[:-1], args_list[-1], truncate
    return args, value, truncate


class CliffordMixin:
    """Clifford gates without per-gate truncation control."""

    _interface: Any

    # Clifford operations
    def x(self, *targets: int | Iterable[int]) -> None:
        """Apply a Pauli X gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.x(_normalize_targets(targets))

    def y(self, *targets: int | Iterable[int]) -> None:
        """Apply a Pauli Y gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.y(_normalize_targets(targets))

    def z(self, *targets: int | Iterable[int]) -> None:
        """Apply a Pauli Z gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.z(_normalize_targets(targets))

    def h(self, *targets: int | Iterable[int]) -> None:
        """Apply a Hadamard gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.h(_normalize_targets(targets))

    def s(self, *targets: int | Iterable[int]) -> None:
        """Apply an S gate (sqrt(Z)) to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.s(_normalize_targets(targets))

    def cnot(self, *targets: int | Iterable[int]) -> None:
        """Apply CNOT (controlled-X) gates over consecutive control/target pairs.

        Args:
            *targets: A flat list of qubit indices, broadcast as
                ``(control, target)`` pairs.
        """
        self._interface.cnot(_normalize_targets(targets))

    def cz(self, *targets: int | Iterable[int]) -> None:
        """Apply CZ (controlled-Z) gates over consecutive qubit pairs.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
        """
        self._interface.cz(_normalize_targets(targets))

    # stim aliases
    cx = cnot
    zcx = cnot
    zcz = cz


class RotationsMixin:
    """Rotation gates without per-gate truncation control."""

    _interface: Any

    # Rotations
    def rx(self, *args: Any, theta: float | None = None) -> None:
        """Apply an RX rotation gate to each target qubit.

        ```math
        R_X(\\theta) = e^{-i \\frac{\\theta}{2} X} = \\cos\\frac{\\theta}{2} I - i \\sin\\frac{\\theta}{2} X
        ```

        Args:
            *targets: The indices of the target qubits.
            theta: The rotation angle in radians.
        """
        targets, theta = _split_targets_parameter(args, theta, "theta")
        self._interface.rx(targets, theta)

    def ry(self, *args: Any, theta: float | None = None) -> None:
        """Apply an RY rotation gate to each target qubit.

        ```math
        R_Y(\\theta) = e^{-i \\frac{\\theta}{2} Y} = \\cos\\frac{\\theta}{2} I - i \\sin\\frac{\\theta}{2} Y
        ```

        Args:
            *targets: The indices of the target qubits.
            theta: The rotation angle in radians.
        """
        targets, theta = _split_targets_parameter(args, theta, "theta")
        self._interface.ry(targets, theta)

    def rz(self, *args: Any, theta: float | None = None) -> None:
        """Apply an RZ rotation gate to each target qubit.

        ```math
        R_Z(\\theta) = e^{-i \\frac{\\theta}{2} Z} = \\cos\\frac{\\theta}{2} I - i \\sin\\frac{\\theta}{2} Z
        ```

        Args:
            *targets: The indices of the target qubits.
            theta: The rotation angle in radians.
        """
        targets, theta = _split_targets_parameter(args, theta, "theta")
        self._interface.rz(targets, theta)

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
    def rxx(self, *args: Any, theta: float | None = None) -> None:
        """Apply RXX (Ising XX) rotation gates over consecutive qubit pairs.

        ```math
        R_{XX}(\\theta) = e^{-i \\frac{\\theta}{2} X \\otimes X}
        ```

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            theta: The rotation angle in radians.
        """
        targets, theta = _split_targets_parameter(args, theta, "theta")
        self._interface.rxx(targets, theta)

    def ryy(self, *args: Any, theta: float | None = None) -> None:
        """Apply RYY (Ising YY) rotation gates over consecutive qubit pairs.

        ```math
        R_{YY}(\\theta) = e^{-i \\frac{\\theta}{2} Y \\otimes Y}
        ```

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            theta: The rotation angle in radians.
        """
        targets, theta = _split_targets_parameter(args, theta, "theta")
        self._interface.ryy(targets, theta)

    def rzz(self, *args: Any, theta: float | None = None) -> None:
        """Apply RZZ (Ising ZZ) rotation gates over consecutive qubit pairs.

        ```math
        R_{ZZ}(\\theta) = e^{-i \\frac{\\theta}{2} Z \\otimes Z}
        ```

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            theta: The rotation angle in radians.
        """
        targets, theta = _split_targets_parameter(args, theta, "theta")
        self._interface.rzz(targets, theta)


class CliffordExtensionMixin:
    """Additional Clifford gates without per-gate truncation control."""

    _interface: Any

    def s_dag(self, *targets: int | Iterable[int]) -> None:
        """Apply an S adjoint gate (sqrt(Z) dagger) to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.s_dag(_normalize_targets(targets))

    def sqrt_x(self, *targets: int | Iterable[int]) -> None:
        """Apply a sqrt(X) gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.sqrt_x(_normalize_targets(targets))

    def sqrt_y(self, *targets: int | Iterable[int]) -> None:
        """Apply a sqrt(Y) gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.sqrt_y(_normalize_targets(targets))

    def sqrt_x_dag(self, *targets: int | Iterable[int]) -> None:
        """Apply a sqrt(X) adjoint gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.sqrt_x_dag(_normalize_targets(targets))

    def sqrt_y_dag(self, *targets: int | Iterable[int]) -> None:
        """Apply a sqrt(Y) adjoint gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
        """
        self._interface.sqrt_y_dag(_normalize_targets(targets))

    def cy(self, *targets: int | Iterable[int]) -> None:
        """Apply controlled-Y gates over consecutive control/target pairs.

        Args:
            *targets: A flat list of qubit indices, broadcast as
                ``(control, target)`` pairs.
        """
        self._interface.cy(_normalize_targets(targets))

    # stim alias
    zcy = cy


class NoiseMixin:
    """Noise channels without per-gate truncation control."""

    _interface: Any

    # Noise operations
    def pauli_error(self, *args: Any, p: Sequence[float] | None = None) -> None:
        """Apply a single-qubit Pauli error channel to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            p: Error probabilities [p_x, p_y, p_z] for X, Y, Z errors.
                The identity probability is implicitly 1 - sum(p).
        """
        targets, p = _split_targets_parameter(args, p, "p")
        self._interface.pauli_error(targets, p)

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

    def two_qubit_pauli_error(self, *args: Any, p: Sequence[float] | None = None) -> None:
        """Apply a two-qubit Pauli error channel over consecutive qubit pairs.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            p: Error probabilities for the 15 non-identity two-qubit Pauli
                operators. Use two_qubit_pauli_error_probabilities to convert
                from a dictionary format.
        """
        targets, p = _split_targets_parameter(args, p, "p")
        self._interface.two_qubit_pauli_error(targets, p)

    # additional noise methods
    def x_error(self, *args: Any, p: float | None = None) -> None:
        """Apply an X error channel to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            p: The probability of applying an X error.
        """
        targets, p = _split_targets_parameter(args, p, "p")
        self._interface.x_error(targets, p)

    def y_error(self, *args: Any, p: float | None = None) -> None:
        """Apply a Y error channel to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            p: The probability of applying a Y error.
        """
        targets, p = _split_targets_parameter(args, p, "p")
        self._interface.y_error(targets, p)

    def z_error(self, *args: Any, p: float | None = None) -> None:
        """Apply a Z error channel to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            p: The probability of applying a Z error.
        """
        targets, p = _split_targets_parameter(args, p, "p")
        self._interface.z_error(targets, p)

    def depolarize1(self, *args: Any, p: float | None = None) -> None:
        """Apply a single-qubit depolarizing channel to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            p: The depolarizing probability.
        """
        targets, p = _split_targets_parameter(args, p, "p")
        self._interface.depolarize1(targets, p)

    def depolarize2(self, *args: Any, p: float | None = None) -> None:
        """Apply a two-qubit depolarizing channel over consecutive qubit pairs.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            p: The depolarizing probability.
        """
        targets, p = _split_targets_parameter(args, p, "p")
        self._interface.depolarize2(targets, p)


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
# have no truncation concept (``two_qubit_pauli_error_probabilities``) are
# inherited unchanged and only the branching gates are overridden. Aliases are
# re-declared at class level because they must point at the overridden methods.


class TruncatingCliffordMixin(CliffordMixin):
    """Clifford gates with a per-gate ``truncate`` kwarg (used by PauliSum)."""

    def x(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply a Pauli X gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            truncate: If ``True`` (default), run the configured truncation
                strategy after the gate; if ``False``, leave the map
                untruncated so the next gate sees the full unpruned state.
        """
        self._interface.x(_normalize_targets(targets), truncate)

    def y(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply a Pauli Y gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            truncate: See `x`.
        """
        self._interface.y(_normalize_targets(targets), truncate)

    def z(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply a Pauli Z gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            truncate: See `x`.
        """
        self._interface.z(_normalize_targets(targets), truncate)

    def h(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply a Hadamard gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            truncate: See `x`.
        """
        self._interface.h(_normalize_targets(targets), truncate)

    def s(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply an S gate (sqrt(Z)) to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            truncate: See `x`.
        """
        self._interface.s(_normalize_targets(targets), truncate)

    def cnot(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply CNOT (controlled-X) gates over consecutive control/target pairs.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            truncate: See `x`.
        """
        self._interface.cnot(_normalize_targets(targets), truncate)

    def cz(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply CZ (controlled-Z) gates over consecutive qubit pairs.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            truncate: See `x`.
        """
        self._interface.cz(_normalize_targets(targets), truncate)

    # stim aliases
    cx = cnot
    zcx = cnot
    zcz = cz


class TruncatingRotationsMixin(RotationsMixin):
    """Rotation gates with a per-gate ``truncate`` kwarg (used by PauliSum)."""

    def rx(self, *args: Any, theta: float | None = None, truncate: bool = True) -> None:
        """Apply an RX rotation gate to each target qubit.

        See `RotationsMixin.rx` for the gate definition.

        Args:
            *targets: The indices of the target qubits.
            theta: The rotation angle in radians.
            truncate: If ``True`` (default), run the configured truncation
                strategy after the gate; if ``False``, defer it.
        """
        targets, theta, truncate = _split_targets_parameter_truncate(args, theta, "theta", truncate)
        self._interface.rx(targets, theta, truncate)

    def ry(self, *args: Any, theta: float | None = None, truncate: bool = True) -> None:
        """Apply an RY rotation gate to each target qubit.

        See `RotationsMixin.ry` for the gate definition.

        Args:
            *targets: The indices of the target qubits.
            theta: The rotation angle in radians.
            truncate: See `rx`.
        """
        targets, theta, truncate = _split_targets_parameter_truncate(args, theta, "theta", truncate)
        self._interface.ry(targets, theta, truncate)

    def rz(self, *args: Any, theta: float | None = None, truncate: bool = True) -> None:
        """Apply an RZ rotation gate to each target qubit.

        See `RotationsMixin.rz` for the gate definition.

        Args:
            *targets: The indices of the target qubits.
            theta: The rotation angle in radians.
            truncate: See `rx`.
        """
        targets, theta, truncate = _split_targets_parameter_truncate(args, theta, "theta", truncate)
        self._interface.rz(targets, theta, truncate)

    def r(self, addr0: int, axis_angle: float, theta: float, *, truncate: bool = True) -> None:
        """Apply a rotation about an axis in the X-Y plane to the specified qubit.

        See `RotationsMixin.r` for the gate definition.

        Args:
            addr0: The index of the target qubit.
            axis_angle: The angle ``φ`` (in radians) of the rotation axis
                within the X-Y plane, measured from the X-axis.
            theta: The rotation angle in radians.
            truncate: See `rx`.
        """
        self._interface.r(addr0, axis_angle, theta, truncate=truncate)

    # Two qubit rotations
    def rxx(self, *args: Any, theta: float | None = None, truncate: bool = True) -> None:
        """Apply RXX (Ising XX) rotation gates over consecutive qubit pairs.

        See `RotationsMixin.rxx` for the gate definition.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            theta: The rotation angle in radians.
            truncate: If ``True`` (default), run the configured truncation
                strategy after the gate. Set to ``False`` to compose a
                U(1)-conserving step like ``rxx + ryy`` on the same edge
                without dropping the conserved-charge component between
                them — then call `PauliSum.truncate` once at the
                end of the composition.
        """
        targets, theta, truncate = _split_targets_parameter_truncate(args, theta, "theta", truncate)
        self._interface.rxx(targets, theta, truncate)

    def ryy(self, *args: Any, theta: float | None = None, truncate: bool = True) -> None:
        """Apply RYY (Ising YY) rotation gates over consecutive qubit pairs.

        See `RotationsMixin.ryy` for the gate definition.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            theta: The rotation angle in radians.
            truncate: See `rxx`.
        """
        targets, theta, truncate = _split_targets_parameter_truncate(args, theta, "theta", truncate)
        self._interface.ryy(targets, theta, truncate)

    def rzz(self, *args: Any, theta: float | None = None, truncate: bool = True) -> None:
        """Apply RZZ (Ising ZZ) rotation gates over consecutive qubit pairs.

        See `RotationsMixin.rzz` for the gate definition.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            theta: The rotation angle in radians.
            truncate: See `rxx`.
        """
        targets, theta, truncate = _split_targets_parameter_truncate(args, theta, "theta", truncate)
        self._interface.rzz(targets, theta, truncate)


class TruncatingCliffordExtensionMixin(CliffordExtensionMixin):
    """Additional Clifford gates with a per-gate ``truncate`` kwarg."""

    def s_dag(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply an S adjoint gate (sqrt(Z) dagger) to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            truncate: See `TruncatingCliffordMixin.x`.
        """
        self._interface.s_dag(_normalize_targets(targets), truncate)

    def sqrt_x(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply a sqrt(X) gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            truncate: See `TruncatingCliffordMixin.x`.
        """
        self._interface.sqrt_x(_normalize_targets(targets), truncate)

    def sqrt_y(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply a sqrt(Y) gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            truncate: See `TruncatingCliffordMixin.x`.
        """
        self._interface.sqrt_y(_normalize_targets(targets), truncate)

    def sqrt_x_dag(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply a sqrt(X) adjoint gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            truncate: See `TruncatingCliffordMixin.x`.
        """
        self._interface.sqrt_x_dag(_normalize_targets(targets), truncate)

    def sqrt_y_dag(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply a sqrt(Y) adjoint gate to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            truncate: See `TruncatingCliffordMixin.x`.
        """
        self._interface.sqrt_y_dag(_normalize_targets(targets), truncate)

    def cy(self, *targets: int | Iterable[int], truncate: bool = True) -> None:
        """Apply controlled-Y gates over consecutive control/target pairs.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            truncate: See `TruncatingCliffordMixin.x`.
        """
        self._interface.cy(_normalize_targets(targets), truncate)

    # stim alias
    zcy = cy


class TruncatingNoiseMixin(NoiseMixin):
    """Noise channels with a per-gate ``truncate`` kwarg (used by PauliSum).

    ``two_qubit_pauli_error_probabilities`` is inherited unchanged from
    `NoiseMixin`.
    """

    def pauli_error(
        self,
        *args: Any,
        p: Sequence[float] | None = None,
        truncate: bool = True,
    ) -> None:
        """Apply a single-qubit Pauli error channel to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            p: Error probabilities [p_x, p_y, p_z] for X, Y, Z errors.
                The identity probability is implicitly 1 - sum(p).
            truncate: If ``True`` (default), run the configured truncation
                strategy after the channel; if ``False``, defer it.
        """
        targets, p, truncate = _split_targets_parameter_truncate(args, p, "p", truncate)
        self._interface.pauli_error(targets, p, truncate)

    def two_qubit_pauli_error(
        self,
        *args: Any,
        p: Sequence[float] | None = None,
        truncate: bool = True,
    ) -> None:
        """Apply a two-qubit Pauli error channel over consecutive qubit pairs.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            p: Error probabilities for the 15 non-identity two-qubit Pauli
                operators. Use two_qubit_pauli_error_probabilities to convert
                from a dictionary format.
            truncate: See `pauli_error`.
        """
        targets, p, truncate = _split_targets_parameter_truncate(args, p, "p", truncate)
        self._interface.two_qubit_pauli_error(targets, p, truncate)

    def x_error(self, *args: Any, p: float | None = None, truncate: bool = True) -> None:
        """Apply an X error channel to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            p: The probability of applying an X error.
            truncate: See `pauli_error`.
        """
        targets, p, truncate = _split_targets_parameter_truncate(args, p, "p", truncate)
        self._interface.x_error(targets, p, truncate)

    def y_error(self, *args: Any, p: float | None = None, truncate: bool = True) -> None:
        """Apply a Y error channel to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            p: The probability of applying a Y error.
            truncate: See `pauli_error`.
        """
        targets, p, truncate = _split_targets_parameter_truncate(args, p, "p", truncate)
        self._interface.y_error(targets, p, truncate)

    def z_error(self, *args: Any, p: float | None = None, truncate: bool = True) -> None:
        """Apply a Z error channel to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            p: The probability of applying a Z error.
            truncate: See `pauli_error`.
        """
        targets, p, truncate = _split_targets_parameter_truncate(args, p, "p", truncate)
        self._interface.z_error(targets, p, truncate)

    def depolarize1(self, *args: Any, p: float | None = None, truncate: bool = True) -> None:
        """Apply a single-qubit depolarizing channel to each target qubit.

        Args:
            *targets: The indices of the target qubits.
            p: The depolarizing probability.
            truncate: See `pauli_error`.
        """
        targets, p, truncate = _split_targets_parameter_truncate(args, p, "p", truncate)
        self._interface.depolarize1(targets, p, truncate)

    def depolarize2(self, *args: Any, p: float | None = None, truncate: bool = True) -> None:
        """Apply a two-qubit depolarizing channel over consecutive qubit pairs.

        Args:
            *targets: A flat list of qubit indices, broadcast as pairs.
            p: The depolarizing probability.
            truncate: See `pauli_error`.
        """
        targets, p, truncate = _split_targets_parameter_truncate(args, p, "p", truncate)
        self._interface.depolarize2(targets, p, truncate)
