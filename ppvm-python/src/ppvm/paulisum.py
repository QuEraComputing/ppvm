# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

import math
import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Self, Union

import ppvm_python_native

from .mixins import CliffordExtensionMixin, CliffordMixin, NoiseMixin, RotationsMixin

_COMPACT_RE = re.compile(r"^([IXYZ]\d+)+$")
_COMPACT_TOKEN_RE = re.compile(r"([IXYZ])(\d+)")


def _parse_term(term: "str | tuple[str, float]", n_qubits: int) -> "tuple[str, float]":
    if isinstance(term, tuple):
        s, coeff = term
    else:
        s, coeff = term, 1.0

    if _COMPACT_RE.match(s):
        chars = ["I"] * n_qubits
        for pauli, idx_str in _COMPACT_TOKEN_RE.findall(s):
            idx = int(idx_str)
            if idx >= n_qubits:
                raise ValueError(f"Qubit index {idx} out of range for {n_qubits}-qubit system.")
            chars[idx] = pauli
        s = "".join(chars)

    return s, coeff


PauliSumInterface = Union[
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

LossyPauliSumInterface = Union[
    ppvm_python_native.PauliSumLossIndexMapFxHash0,
    ppvm_python_native.PauliSumLossIndexMapFxHash1,
    ppvm_python_native.PauliSumLossIndexMapFxHash2,
    ppvm_python_native.PauliSumLossIndexMapFxHash3,
    ppvm_python_native.PauliSumLossIndexMapFxHash4,
    ppvm_python_native.PauliSumLossIndexMapFxHash5,
    ppvm_python_native.PauliSumLossIndexMapFxHash6,
    ppvm_python_native.PauliSumLossIndexMapFxHash7,
    ppvm_python_native.PauliSumLossIndexMapFxHash8,
    ppvm_python_native.PauliSumLossIndexMapFxHash9,
    ppvm_python_native.PauliSumLossIndexMapFxHash10,
    ppvm_python_native.PauliSumLossIndexMapFxHash11,
    ppvm_python_native.PauliSumLossIndexMapFxHash12,
    ppvm_python_native.PauliSumLossIndexMapFxHash13,
    ppvm_python_native.PauliSumLossIndexMapFxHash14,
    ppvm_python_native.PauliSumLossIndexMapFxHash15,
]


@dataclass(frozen=True)
class PauliSum(
    CliffordExtensionMixin,
    CliffordMixin,
    NoiseMixin,
    RotationsMixin,
):
    """A weighted sum of Pauli strings for quantum simulation.

    PauliSum represents a linear combination of Pauli operators, commonly used
    to represent quantum observables or Hamiltonians. It provides methods for
    applying quantum gates (Clifford operations and rotations) and computing
    expectation values via the trace operation.

    Attributes:
        initial_terms: Pauli strings, each containing only 'I', 'X', 'Y', 'Z' characters.
            All terms must have the same length (number of qubits).
        n_qubits: Number of qubits.
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
        ps = PauliSum(n_qubits = 2, initial_terms=["ZZ", "XI"], coefficients=[0.5, 0.3])
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
        ps = PauliSum.new(3, "ZZZ")
        # GHZ circuit: H(0), CNOT(0,1), CNOT(1,2)
        # Apply in reverse order:
        ps.cnot(1, 2)
        ps.cnot(0, 1)
        ps.h(0)
        # Expectation value of ZZZ for GHZ state is 0
        result = ps.overlap_with_zero()
        ```
    """

    initial_terms: Sequence[str]
    n_qubits: int
    coefficients: Sequence[float] = ()
    min_abs_coeff: float = 1e-10
    max_pauli_weight: int | None = None
    max_loss_weight: int | None = None
    preserve_strings: Sequence[str] | None = None

    _interface: PauliSumInterface = field(init=False, repr=False)

    def __post_init__(self):
        object.__setattr__(
            self,
            "_interface",
            self._init_ppvm_interface(),
        )

    def _get_interface(self, n_interface: int):
        return getattr(ppvm_python_native, f"PauliSumIndexMapFxHash{n_interface}")

    def _init_ppvm_interface(
        self,
    ):

        n_qubits = self.n_qubits
        terms = self.initial_terms
        coefficients = self.coefficients

        if not terms:
            raise ValueError("At least one term must be provided to initialize PauliSum.")

        for term in terms:
            if len(term) != n_qubits:
                raise ValueError(
                    "All terms must have the same length! Expected length "
                    f"{n_qubits}, but got term of length {len(term)}: {term!r}"
                )

        # number of bytes we need
        N = math.ceil(n_qubits / 8.0)

        # number of bytes we have
        possible_interfaces = range(16)
        N_interface = next(n for n in possible_interfaces if 2**n > N)
        interface = self._get_interface(N_interface)

        if terms and not coefficients:
            coefficients = (1.0,) * len(terms)

        # set the kwargs for the interface
        options = {
            "min_abs_coeff": self.min_abs_coeff,
            "terms": terms,
            "coefficients": coefficients,
        }

        # these are just set to the defaults on the rust side if None
        if self.max_pauli_weight is not None:
            options["max_pauli_weight"] = self.max_pauli_weight

        if self.max_loss_weight is not None:
            options["max_loss_weight"] = self.max_loss_weight

        if self.preserve_strings:
            preserve_list = list(self.preserve_strings)
            for s in preserve_list:
                if len(s) != n_qubits:
                    raise ValueError(
                        "All preserve strings must have length n_qubits "
                        f"({n_qubits}); got {len(s)}: {s!r}"
                    )
            options["preserve_strings"] = preserve_list

        return interface(
            n_qubits,
            **options,
        )

    def __len__(self) -> int:
        """Get the number of terms in the PauliSum.

        Returns:
            The number of Pauli terms.
        """
        return len(self._interface)

    @classmethod
    def new(
        cls,
        n_qubits: int,
        terms: "str | tuple | list",
        min_abs_coeff: float = 1e-10,
        max_pauli_weight: int | None = None,
        max_loss_weight: int | None = None,
        preserve_strings: Sequence[str] | None = None,
    ) -> Self:
        """Create a PauliSum from one or more terms with flexible input formats.

        Args:
            n_qubits: Number of qubits.
            terms: A single term or list of terms. Each term is either:
                - A full Pauli string (e.g. ``"IX"``), with coefficient 1.0.
                - A compact string ``"P{i}"`` (e.g. ``"X1"``), placing Pauli P
                  at 0-based qubit index i with coefficient 1.0.
                - A tuple ``(str, float)`` pairing either of the above with an
                  explicit coefficient.
            min_abs_coeff: Terms with absolute coefficient below this threshold
                are dropped. Defaults to 1e-10.
            max_pauli_weight: Maximum number of non-identity Paulis per term.
                If None, truncation is disabled.
            max_loss_weight: Maximum loss weight per term (only used by
                LossyPauliSum). If None, truncation is disabled.
                Note, that this should usually be chosen to be quite low, since
                e.g. 10 would correspond to keeping terms that contribute if
                up to 10 qubits are lost simultaneously.
            preserve_strings: Pauli strings (length ``n_qubits`` each) that
                truncation must never drop. Empty by default.

        Returns:
            A new instance of the class this method is called on.

        Raises:
            ValueError: If a compact qubit index is out of range for n_qubits.

        Example:
            Full Pauli string with implicit coefficient 1.0::

                ps = PauliSum.new(2, "IX")

            Full string with explicit coefficient::

                ps = PauliSum.new(2, ("IX", 0.5))

            Compact notation — X on qubit 1 in a 3-qubit system::

                ps = PauliSum.new(3, "X1")  # equivalent to PauliSum.new(3, "IXI")

            Multiple terms mixing both notations::

                ps = PauliSum.new(3, [("Y1", 0.1), "ZIZ"])

            Building a Z-basis observable for each qubit::

                n = 5
                ps = PauliSum.new(n, [f"Z{i}" for i in range(n)])

            Creating a lossy PauliSum with a loss weight cutoff::

                ps = LossyPauliSum.new(3, "ZZZ", max_loss_weight=1)
                ps.loss_channel(0, 0.01)
        """
        if isinstance(terms, (str, tuple)):
            terms = [terms]
        parsed = [_parse_term(t, n_qubits) for t in terms]
        return cls(
            n_qubits=n_qubits,
            initial_terms=[s for s, _ in parsed],
            coefficients=[c for _, c in parsed],
            min_abs_coeff=min_abs_coeff,
            max_pauli_weight=max_pauli_weight,
            max_loss_weight=max_loss_weight,
            preserve_strings=preserve_strings,
        )

    def __str__(self) -> str:
        return self._interface.__str__()

    def __copy__(self) -> Self:
        new = object.__new__(type(self))
        object.__setattr__(new, "initial_terms", self.initial_terms)
        object.__setattr__(new, "n_qubits", self.n_qubits)
        object.__setattr__(new, "coefficients", self.coefficients)
        object.__setattr__(new, "min_abs_coeff", self.min_abs_coeff)
        object.__setattr__(new, "max_pauli_weight", self.max_pauli_weight)
        object.__setattr__(new, "max_loss_weight", self.max_loss_weight)
        object.__setattr__(new, "preserve_strings", self.preserve_strings)
        object.__setattr__(new, "_interface", self._interface.__copy__())
        return new

    def copy(self) -> Self:
        """Create a copy of the PauliSum instance.

        Returns:
            A new PauliSum instance that is a copy of the current one.
        """
        return self.__copy__()

    @property
    def terms(self) -> list[tuple[str, float]]:
        """Get the list of Pauli terms and their coefficients.

        Returns:
            A list of tuples, each containing a Pauli string and its coefficient.
        """
        return self._interface.terms()

    def weights(self) -> list[tuple[str, int]]:
        """Get the weight of each Pauli term.

        Returns:
            A list of tuples, each containing a Pauli string and its weight.
        """
        return self._interface.weights()

    def current_max_weight(self) -> int:
        """Get the current maximum weight of the Pauli sum.

        Returns:
            The weight as integer.
        """
        return self._interface.current_max_weight()

    # Getting results
    def overlap_with_zero(self) -> float:
        """Compute the overlap with the all-zeros computational basis state.

        Returns:
            The expectation value of the Pauli sum with respect to |0...0>.
        """
        return self._interface.trace("Z?*")

    def overlap(self, other: "PauliSum") -> float:
        """Compute the overlap of the current PauliSum with another PauliSum.
        The overlap is defined as

        ```math
        \\text{tr}(A^\\dagger \\cdot B),
        ```

        where
            self -> A
            other -> B

        Returns:
            The trace overlap of the PauliSums.
        """
        return self._interface.overlap(other._interface)

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

    def symmetry_merge(self, group) -> None:
        """Merge entries into orbit-representative form under a translation group.

        Each Pauli word in the sum is replaced by its canonical (lex-min)
        representative under the action of ``group``; coefficients of words
        that collapse to the same representative are summed. Entry count
        reduces by up to ``|group|×`` for translation-invariant operators.

        For a translation-invariant dynamics that you apply between
        merging steps, this preserves all ``group``-invariant expectation
        values (Theorem 1 of Teng et al., arXiv:2512.12094). Plain
        real-coefficient merge — handles the trivial (``k=0``) momentum
        sector only.

        Args:
            group: A :class:`ppvm_python_native.TranslationGroup`
                (use ``TranslationGroup.chain_1d(n)``, ``.torus_2d``,
                ``.torus_3d``, ``.ladder``, or ``.from_generators``).
        """
        self._interface.symmetry_merge(group)

    def amplitude_damping(self, addr0: int, gamma: float, *, truncate: bool = True):
        """Apply an amplitude-damping channel.

        Args:
            addr0: The index of the target qubit.
            gamma: The damping rate. `X` and `Y` are damped with `sqrt(1 - gamma)`,
                whereas `Z` branches into `gamma * I + (1 - gamma) * Z`.
            truncate: If ``True`` (default), run the configured truncation
                strategy after the channel; if ``False``, defer it (use
                :meth:`truncate` to fire the cut later).
        """

        # TODO: move to mixins once also implemented for GeneralizedTableau

        self._interface.amplitude_damping(addr0, gamma, truncate=truncate)

    def truncate(self) -> None:
        """Run the configured truncation strategy (``min_abs_coeff`` and/or
        ``max_pauli_weight``) on the current state.

        Useful when gates were called with ``truncate=False`` to chain a
        composition of commuting operations (e.g. ``rxx + ryy`` on the
        same edge, an exchange-like step that conserves total ``Z``), so
        that intermediate truncation does not drop conserved-charge
        components. Call :meth:`truncate` once at the end of the
        composition to apply the cut to the combined result.
        """
        self._interface.truncate()


@dataclass(frozen=True)
class LossyPauliSum(PauliSum):
    """A PauliSum that supports modelling qubit loss.

    This is achieved by extending the set of Pauli basis operators to include
    an addition operator ``{I, X, Y, Z, L}``, where ``L`` is the projector on
    a third leakage state. This basis effectively allows simulating qutrits,
    where we neglect any coherences between the qubit subspace and the leakage
    state.

    In addition to the new channels, there is also another truncation strategy:
    Since Pauli Strings that have an `L` at multiple positions contribute only
    minimally, these can get truncated by setting an appropriate `max_loss_weight`.
    The truncation is similar to how `max_pauli_weight` truncates strings, but only
    counting `L`s.
    """

    _interface: LossyPauliSumInterface = field(init=False, repr=False)

    def _get_interface(self, n_interface: int):
        return getattr(ppvm_python_native, f"PauliSumLossIndexMapFxHash{n_interface}")

    # NOTE: purposely not using mixin for better docstrings here

    def loss_channel(self, addr0: int, p: float, *, truncate: bool = True) -> None:
        """Apply a single-qubit loss channel.

        Reduces the trace of qubit-subspace operators by `(1 - p)`.
        Adds back population into `I` or `Z` if the Pauli string has an `L`
        at `addr0`. This can only occur if a reset channel has been applied before
        and accounts for the fact of falsely counting a lost qubit as 0 in a
        measurement.

        Args:
            addr0: The index of the target qubit.
            p: Loss probability in [0, 1].
            truncate: If ``True`` (default), run the configured truncation
                strategy after the channel; if ``False``, defer it.
        """
        self._interface.loss_channel(addr0, p, truncate=truncate)

    def correlated_loss_channel(
        self,
        addr0: int,
        addr1: int,
        p: Sequence[float],
        *,
        truncate: bool = True,
    ) -> None:
        """Apply a correlated loss channel.

        This applies a correlated loss channel to the qubits at `addr0` and `addr1`.
        The channel accepts 3 probabilities as argument:
            * `p[0]`: The probability of losing both qubits, when they are originally
                in the qubit subspace.
            * `p[1]`: The probability of losing a single qubit, when both qubits
                are originally in the qubit subspace.
            * `p[2]`: The probability of losing one qubit when the other one
                has already been lost prior to applying the channel. This is to
                account for the fact that when one qubit is missing during e.g.
                a controlled gate, the remaining qubit undergoes a different dynamic.
                We account for this difference with this distinct probability.
            truncate: If ``True`` (default), run the configured truncation
                strategy after the channel; if ``False``, defer it.
        """
        self._interface.correlated_loss_channel(addr0, addr1, p, truncate=truncate)

    def reset_loss_channel(self, addr0: int, *, truncate: bool = True) -> None:
        """Reset a lost qubit to the 0 state. Usually, you want to apply
        this channel at the end of the circuit, i.e. at the beginning when
        propagating backwards.

        **NOTE**: This channel causes exponential branching in `I` and `Z`.
        Make sure to set an appropriate `max_loss_weight` to truncate.

        Args:
            addr0: The index of the qubit to reset.
            truncate: If ``True`` (default), run the configured truncation
                strategy after the channel; if ``False``, defer it.
        """
        self._interface.reset_loss_channel(addr0, truncate=truncate)
