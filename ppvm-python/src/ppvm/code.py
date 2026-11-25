# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Qubit codes specified by signed stabilizer generators."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from operator import index
from typing import NamedTuple

import numpy as np
import stim
from numpy.typing import ArrayLike

from ._core import CodeSpec


def _bits(data: ArrayLike, name: str, ndim: int) -> np.ndarray:
    a = np.asarray(data)
    if a.ndim != ndim:
        raise ValueError(f"{name} must be {ndim}-dimensional")
    if a.size and (a.dtype.kind not in "biu" or not np.all((a == 0) | (a == 1))):
        raise ValueError(f"{name} must contain binary integers or booleans")
    return a.astype(np.bool_, copy=False)


def _strings(rows: Iterable[str | stim.PauliString]) -> list[str]:
    if isinstance(rows, str):
        raise TypeError("supply a sequence of signed Pauli strings, e.g. ['+ZZ']")
    strings = []
    for row in rows:
        if not isinstance(row, str):
            if not isinstance(row, stim.PauliString):
                raise TypeError("expected a signed Pauli string or stim.PauliString")
            row = str(row)
        strings.append(row)
    return strings


def _synthesize_clifford(target: stim.Tableau) -> stim.Circuit:
    circuit = target.to_circuit()
    if circuit.to_tableau() != target:
        raise RuntimeError("synthesized circuit does not implement the requested Clifford")
    return circuit


class Parity(NamedTuple):
    """XOR the indicated record bits and the offset to evaluate the parity."""

    indices: tuple[int, ...]
    offset: bool


class TerminalReadout(NamedTuple):
    detectors: tuple[Parity, ...]
    observables: tuple[Parity, ...]


@dataclass(frozen=True, slots=True, eq=False)
class Code:
    """Joint +1 eigenspace of signed, commuting qubit Pauli generators.

    Column zero is qubit zero. Redundant, consistent checks retain their input
    order. The default logical basis is deterministic; supply both ordered
    logical X and Z lists to use a particular convention.
    """

    _spec: CodeSpec

    @classmethod
    def from_stabilizers(
        cls,
        stabilizers: Iterable[str | stim.PauliString],
        *,
        n_qubits: int | None = None,
        logical_xs: Iterable[str | stim.PauliString] | None = None,
        logical_zs: Iterable[str | stim.PauliString] | None = None,
    ) -> Code:
        """Construct from rows such as '+XZZXI' or '-YY'; I and _ mean identity."""
        rows = _strings(stabilizers)
        if n_qubits is None:
            if not rows:
                raise ValueError("specify n_qubits for an empty generator set")
            n_qubits = len(rows[0].removeprefix("+").removeprefix("-"))
        n = index(n_qubits)
        if n < 0:
            raise ValueError("n_qubits must be nonnegative")
        return cls(
            CodeSpec(
                n,
                rows,
                None if logical_xs is None else _strings(logical_xs),
                None if logical_zs is None else _strings(logical_zs),
            )
        )

    @classmethod
    def from_tableau(
        cls,
        xs: ArrayLike,
        zs: ArrayLike,
        signs: ArrayLike | None = None,
        *,
        logical_xs: Iterable[str | stim.PauliString] | None = None,
        logical_zs: Iterable[str | stim.PauliString] | None = None,
    ) -> Code:
        """Each row is one check. (x,z)=(1,1) is Y; sign 1 denotes minus."""
        x = _bits(xs, "xs", 2)
        z = _bits(zs, "zs", 2)
        if x.shape != z.shape:
            raise ValueError("xs and zs must have identical shapes")
        m = x.shape[0]
        negative = np.zeros(m, dtype=bool) if signs is None else _bits(signs, "signs", 1)
        if len(negative) != m:
            raise ValueError("expected one sign bit per generator")
        return cls(
            CodeSpec.from_tableau(
                x,
                z,
                negative,
                None if logical_xs is None else _strings(logical_xs),
                None if logical_zs is None else _strings(logical_zs),
            )
        )

    @classmethod
    def from_css(
        cls,
        hx: ArrayLike,
        hz: ArrayLike,
        *,
        logical_xs: Iterable[str | stim.PauliString] | None = None,
        logical_zs: Iterable[str | stim.PauliString] | None = None,
    ) -> Code:
        """Construct from positive X and Z checks with the same physical-qubit order."""
        x = _bits(hx, "hx", 2)
        z = _bits(hz, "hz", 2)
        if x.shape[1] != z.shape[1]:
            raise ValueError("hx and hz must have the same number of columns")
        return cls.from_tableau(
            np.vstack((x, np.zeros_like(z))),
            np.vstack((np.zeros_like(x), z)),
            logical_xs=logical_xs,
            logical_zs=logical_zs,
        )

    @property
    def n(self) -> int:
        """Number of physical qubits."""
        return self._spec.n

    @property
    def k(self) -> int:
        """Number of encoded qubits (n minus stabilizer rank)."""
        return self._spec.k

    @property
    def stabilizers(self) -> tuple[str, ...]:
        """Original signed checks, including redundant rows."""
        return tuple(self._spec.stabilizers)

    @property
    def logical_xs(self) -> tuple[str, ...]:
        """Ordered logical X representatives in physical-qubit order."""
        return tuple(self._spec.logical_xs)

    @property
    def logical_zs(self) -> tuple[str, ...]:
        """Ordered logical Z representatives paired with logical_xs."""
        return tuple(self._spec.logical_zs)

    def state_stabilizers(self, bits: Sequence[int] | None = None) -> tuple[str, ...]:
        """Independent signed generators of a logical computational-basis state.

        Defaults to all logical zeros. bits[i] sets the eigenvalue of logical_zs[i].
        """
        signs = np.zeros(self.k, dtype=bool) if bits is None else _bits(bits, "bits", 1)
        if len(signs) != self.k:
            raise ValueError(f"expected {self.k} logical state bits")
        checks = self._spec.independent_stabilizers
        for z, toggle in zip(self.logical_zs, signs, strict=True):
            checks.append(("-" if z.startswith("+") else "+") + z[1:] if toggle else z)
        return tuple(checks)

    def synthesize_state_prep(self, bits: Sequence[int] | None = None) -> stim.Circuit:
        """Prepare and verify the signed logical basis state. Returns a Stim circuit."""
        rows = [stim.PauliString(p) for p in self.state_stabilizers(bits)]
        return _synthesize_clifford(stim.Tableau.from_stabilizers(rows))

    def terminal_readout(
        self,
        bases: str,
        *,
        order: Sequence[int] | None = None,
        logical_bases: str | None = None,
    ) -> TerminalReadout:
        """Derive parities for destructive single-qubit Pauli measurements.

        bases[q] is the physical X/Y/Z axis on qubit q. order[j] is the qubit
        measured into record j (default: identity). logical_bases chooses X, Y or
        Z for each logical qubit, or I to skip it; defaults to all Z. Detectors
        form an independent basis of measurable stabilizers. Observables are
        returned in logical order, using equivalent representatives when needed.
        """
        measured = range(self.n) if order is None else [index(q) for q in order]
        if sorted(measured) != list(range(self.n)):
            raise ValueError("order must be a permutation of the physical qubits")
        records = {q: j for j, q in enumerate(measured)}
        dets, obs = self._spec.terminal_readout(
            bases, "Z" * self.k if logical_bases is None else logical_bases
        )

        def remap(rows):
            return tuple(
                Parity(tuple(sorted(records[q] for q in support)), sign) for support, sign in rows
            )

        return TerminalReadout(remap(dets), remap(obs))

    def contains_stabilizer(self, pauli: str | stim.PauliString) -> bool:
        """Test membership in the signed stabilizer group, including the phase."""
        return self._spec.contains_stabilizer(_strings([pauli])[0])

    def encoder_tableau(self) -> stim.Tableau:
        """Forward Stim Clifford tableau. The last k inputs carry logical data.

        The first n-k inputs are zero ancillas. Their Z images are the independent
        checks, and both Pauli images of every data input use our logical basis.
        """
        checks = [stim.PauliString(s) for s in self._spec.independent_stabilizers]
        logical_zs = [stim.PauliString(s) for s in self.logical_zs]
        logical_xs = [stim.PauliString(s) for s in self.logical_xs]
        zs = checks + logical_zs
        completion = stim.Tableau.from_stabilizers(zs)
        xs = []
        for i in range(len(checks)):
            dual = completion.x_output(i)
            # Remove overlap with our chosen logical X representatives while
            # preserving the dual's commutation with every ancilla check.
            for x, z in zip(logical_xs, logical_zs, strict=True):
                if not dual.commutes(x):
                    dual *= z
            xs.append(dual)
        return stim.Tableau.from_conjugated_generators(xs=xs + logical_xs, zs=zs)

    def synthesize_encoder(self) -> stim.Circuit:
        """Return a verified ideal encoding circuit using Stim."""
        return _synthesize_clifford(self.encoder_tableau())

    def synthesize_logical_clifford(self, logical: stim.Tableau) -> stim.Circuit:
        """Return a verified ideal circuit acting as logical on the encoded data.

        Its target is logical decoding, the requested Clifford, then re-encoding.
        Fault tolerance requires a separate circuit-level analysis.
        """
        return _synthesize_clifford(self._logical_target(logical))

    def _logical_target(self, logical: stim.Tableau) -> stim.Tableau:
        if not isinstance(logical, stim.Tableau) or len(logical) != self.k:
            raise ValueError(f"expected a {self.k}-qubit logical stim.Tableau")
        encoder = self.encoder_tableau()
        lifted = stim.Tableau(self.n)
        lifted.append(logical, list(range(self.n - self.k, self.n)))
        return encoder.inverse().then(lifted).then(encoder)

    def implements_logical_clifford(self, physical: stim.Tableau, logical: stim.Tableau) -> bool:
        """Check a physical Clifford's action on the codespace up to global phase."""
        if not isinstance(physical, stim.Tableau) or len(physical) != self.n:
            raise ValueError(f"expected a {self.n}-qubit physical stim.Tableau")
        target = self._logical_target(logical)
        for check in self._spec.independent_stabilizers:
            if not self.contains_stabilizer(physical(stim.PauliString(check))):
                return False
        for rep in (*self.logical_xs, *self.logical_zs):
            p = stim.PauliString(rep)
            residual = physical(p) * target(p)
            if residual.sign.imag or not self.contains_stabilizer(residual):
                return False
        return True
