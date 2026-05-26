# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Direct Pauli-Lindbladian time evolution on an adaptive Pauli-string basis.

Given a Hermitian Pauli Hamiltonian H = Σ c_i P_i and Hermitian Pauli jump
operators L_k with rates γ_k ≥ 0, this module exposes three primitives
needed for adaptive Heisenberg-picture evolution:

- ``action(p)`` / ``action_arr(p)``: L*(p) for one Pauli string p
- ``leakage(basis, coeffs)`` / ``leakage_arr(...)``: off-basis component of
  L*(Σ c_j p_j), driving basis expansion
- ``generator(basis)``: scipy CSC matrix M such that L* restricted to
  ``basis`` is ``M @ coeffs``

The ``*_arr`` variants pass Pauli strings as ``(N, n_qubits)`` ``uint8``
arrays of Pauli codes (``0=I, 1=X, 2=Z, 3=Y``) and skip string
construction entirely — at ~10^5 basis rows per evolution step, per-row
``str.join`` dominates wall time.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import scipy.sparse as sp
from ppvm_python_native import LindbladSpec as _LindbladSpec

_PAULI_CODE = {"I": 0, "X": 1, "Z": 2, "Y": 3}
# Lookup table mapping code -> ASCII byte for vectorised string output.
_CODE_TO_ASCII = np.array([ord("I"), ord("X"), ord("Z"), ord("Y")], dtype=np.uint8)


def string_to_codes(s: str, n_qubits: int) -> np.ndarray:
    """Encode a Pauli string ``"IXYZ..."`` as a length-``n_qubits`` uint8 array."""
    if len(s) != n_qubits:
        raise ValueError(f"Pauli string {s!r} has length {len(s)} != n_qubits {n_qubits}")
    return np.array([_PAULI_CODE[c] for c in s], dtype=np.uint8)


def codes_to_string(codes: np.ndarray) -> str:
    """Decode one length-``n_qubits`` row of Pauli codes back to a string."""
    return _CODE_TO_ASCII[codes].tobytes().decode("ascii")


def basis_to_codes(basis: Sequence[str], n_qubits: int) -> np.ndarray:
    """Stack a sequence of Pauli strings into an ``(N, n_qubits)`` uint8 array."""
    arr = np.zeros((len(basis), n_qubits), dtype=np.uint8)
    for i, s in enumerate(basis):
        arr[i] = string_to_codes(s, n_qubits)
    return arr


def codes_to_basis(arr: np.ndarray) -> list[str]:
    """Inverse of :func:`basis_to_codes`. One call into C per row."""
    bytes_per_row = _CODE_TO_ASCII[arr].tobytes()
    n = arr.shape[1]
    return [bytes_per_row[i * n : (i + 1) * n].decode("ascii") for i in range(arr.shape[0])]


class Lindbladian:
    """Pre-compiled adjoint Pauli-Lindbladian acting on Pauli strings.

    Parameters
    ----------
    n_qubits:
        Number of qubits.
    h_terms:
        Iterable of ``(pauli_string, coefficient)`` pairs for the
        Hermitian Hamiltonian ``H = Σ c_i P_i``. Each ``pauli_string`` is
        a length-``n_qubits`` ``str`` over ``"IXYZ"``.
    jump_terms:
        Iterable of ``(pauli_string, rate)`` pairs for the Hermitian Pauli
        jump operators ``L_k`` with non-negative rates ``γ_k``.
    """

    def __init__(
        self,
        n_qubits: int,
        h_terms: Iterable[tuple[str, float]],
        jump_terms: Iterable[tuple[str, float]] = (),
    ):
        self.n_qubits = int(n_qubits)
        h_strs: list[str] = []
        h_coeffs: list[float] = []
        for s, c in h_terms:
            h_strs.append(s)
            h_coeffs.append(float(c))
        j_strs: list[str] = []
        j_rates: list[float] = []
        for s, g in jump_terms:
            j_strs.append(s)
            j_rates.append(float(g))
        self._spec = _LindbladSpec(self.n_qubits, h_strs, h_coeffs, j_strs, j_rates)

    @property
    def num_h_terms(self) -> int:
        return self._spec.num_h_terms

    @property
    def num_jump_terms(self) -> int:
        return self._spec.num_jump_terms

    @property
    def cache_size(self) -> int:
        return self._spec.cache_size

    def clear_cache(self) -> None:
        self._spec.clear_cache()

    # ── Pure-ndarray hot path ──

    def action_arr(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Apply ``L*`` to a single Pauli string given as uint8 codes.

        Returns ``(out_basis, out_coeffs)``: a ``(M, n_qubits)`` uint8
        array and a length-``M`` float64 array.
        """
        return self._spec.action(np.ascontiguousarray(p, dtype=np.uint8))

    def leakage_arr(
        self,
        basis_arr: np.ndarray,
        coeffs: np.ndarray,
        protected_arr: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Off-basis component of ``L*( Σ_j coeffs[j] basis[j] )``.

        ``basis_arr``: ``(N, n_qubits)`` uint8. ``coeffs``: length-N float64.
        ``protected_arr``: optional ``(K, n_qubits)`` uint8 of Pauli strings
        that must NEVER appear in the leakage output.

        Returns ``(out_basis, out_coeffs)`` packed the same way as
        :meth:`action_arr`.
        """
        n = self.n_qubits
        if protected_arr is None:
            protected_arr = np.zeros((0, n), dtype=np.uint8)
        return self._spec.leakage(
            np.ascontiguousarray(basis_arr, dtype=np.uint8),
            np.ascontiguousarray(coeffs, dtype=np.float64),
            np.ascontiguousarray(protected_arr, dtype=np.uint8),
        )

    def generator_arr(self, basis_arr: np.ndarray) -> sp.csc_matrix:
        """Sparse generator matrix in CSC form, basis given as uint8 codes."""
        n_basis = basis_arr.shape[0]
        rows, cols, vals = self._spec.generator(np.ascontiguousarray(basis_arr, dtype=np.uint8))
        return sp.coo_matrix((vals, (rows, cols)), shape=(n_basis, n_basis)).tocsc()

    # ── String-keyed convenience API (slower; for tests / display) ──

    def action(self, p: str) -> dict[str, float]:
        """Apply ``L*`` to a single Pauli string ``p`` (string-keyed dict)."""
        codes = string_to_codes(p, self.n_qubits)
        out_basis, out_coeffs = self._spec.action(codes)
        keys = codes_to_basis(out_basis)
        return {k: float(v) for k, v in zip(keys, out_coeffs) if v != 0.0}

    def leakage(
        self,
        basis: Sequence[str],
        coeffs: np.ndarray,
        protected: Sequence[str] | None = None,
    ) -> dict[str, float]:
        """Off-basis leakage as a ``dict[str, float]`` (slower API)."""
        n = self.n_qubits
        basis_arr = basis_to_codes(basis, n)
        protected_arr = (
            basis_to_codes(list(protected), n) if protected else np.zeros((0, n), dtype=np.uint8)
        )
        out_basis, out_coeffs = self._spec.leakage(
            basis_arr,
            np.ascontiguousarray(coeffs, dtype=np.float64),
            protected_arr,
        )
        keys = codes_to_basis(out_basis)
        return {k: float(v) for k, v in zip(keys, out_coeffs) if v != 0.0}

    def generator(self, basis: Sequence[str]) -> sp.csc_matrix:
        """Sparse generator matrix in CSC form, basis given as strings."""
        n = self.n_qubits
        basis_arr = basis_to_codes(basis, n)
        return self.generator_arr(basis_arr)
