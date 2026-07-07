# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Direct Pauli-Lindbladian time evolution on an adaptive Pauli-string basis.

Given a Hermitian Pauli Hamiltonian H = Σ c_i P_i and jump operators
L_k = Σ_a λ_{k,a} P_{k,a} (each a complex linear combination of Pauli
strings) with rates γ_k ≥ 0, this module exposes three primitives needed
for adaptive Heisenberg-picture evolution:

- ``action(p)`` / ``action_arr(p)``: L*(p) for one Pauli string p
- ``leakage(basis, coeffs)`` / ``leakage_arr(...)``: off-basis component of
  L*(Σ c_j p_j), driving basis expansion
- ``generator(basis)``: COO triples ``(rows, cols, vals)`` for the generator
  matrix M such that L* restricted to ``basis`` is ``M @ coeffs``. Users
  wanting a sparse matrix can wrap them — e.g.
  ``scipy.sparse.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsc()``

The ``*_arr`` variants pass Pauli strings as ``(N, n_qubits)`` ``uint8``
arrays of Pauli codes (``0=I, 1=X, 2=Z, 3=Y``) and skip string
construction entirely — at ~10^5 basis rows per evolution step, per-row
``str.join`` dominates wall time.

Each jump term can be either:

- a single Hermitian Pauli (`("ZZII", γ)`), routed to a fast diagonal path,
  or
- a complex Pauli sum (`([("XIII", 0.5+0j), ("YIII", 0+0.5j)], γ)`) to
  describe e.g. amplitude-damping (`σ⁻`) and excitation (`σ⁺`) operators.

For the general case the shim evaluates
``γ ( L† p L − ½ {L†L, p} )`` directly; the L†L Pauli expansion is
precomputed once at construction.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Union

import numpy as np
from ._core import LindbladSpec as _LindbladSpec

_PAULI_CODE = {"I": 0, "X": 1, "Z": 2, "Y": 3}
# Lookup table mapping code -> ASCII byte for vectorised string output.
_CODE_TO_ASCII = np.array([ord("I"), ord("X"), ord("Z"), ord("Y")], dtype=np.uint8)

# A jump operator is either a Hermitian Pauli (single string) or a complex
# linear combination of Pauli strings.
PauliLincomb = Iterable[tuple[str, complex]]
JumpSpec = Union[tuple[str, float], tuple[PauliLincomb, float]]


def _string_to_codes(s: str, n_qubits: int) -> np.ndarray:
    """Encode a Pauli string ``"IXYZ..."`` as a length-``n_qubits`` uint8 array.

    Underscores in the input are ignored, matching the Rust parser
    (``parse_pauli_string`` in `ppvm-lindblad`) so users can write
    ``"X_Y_Z"`` for readability.
    """
    s_clean = s.replace("_", "")
    if len(s_clean) != n_qubits:
        raise ValueError(
            f"Pauli string {s!r} has length {len(s_clean)} (after stripping '_') "
            f"!= n_qubits {n_qubits}"
        )
    try:
        return np.array([_PAULI_CODE[c] for c in s_clean], dtype=np.uint8)
    except KeyError as exc:
        bad = exc.args[0]
        raise ValueError(
            f"Pauli string {s!r} contains invalid character {bad!r}; "
            f"expected one of 'I', 'X', 'Y', 'Z' (and '_' is allowed as a separator)"
        ) from None


def _codes_to_string(codes: np.ndarray) -> str:
    """Decode one length-``n_qubits`` row of Pauli codes back to a string."""
    return _CODE_TO_ASCII[codes].tobytes().decode("ascii")


def _basis_to_codes(basis: Sequence[str], n_qubits: int) -> np.ndarray:
    """Stack a sequence of Pauli strings into an ``(N, n_qubits)`` uint8 array."""
    arr = np.zeros((len(basis), n_qubits), dtype=np.uint8)
    for i, s in enumerate(basis):
        arr[i] = _string_to_codes(s, n_qubits)
    return arr


def _codes_to_basis(arr: np.ndarray) -> list[str]:
    """Inverse of :func:`_basis_to_codes`. One call into C per row."""
    bytes_per_row = _CODE_TO_ASCII[arr].tobytes()
    n = arr.shape[1]
    return [bytes_per_row[i * n : (i + 1) * n].decode("ascii") for i in range(arr.shape[0])]


def sigma_plus(site: int, n_qubits: int) -> list[tuple[str, complex]]:
    """``σ⁺_q = (X_q + i Y_q) / 2``. Use as a Lindblad jump for excitation."""
    if not 0 <= site < n_qubits:
        raise ValueError(f"site {site} out of range for n_qubits={n_qubits}")
    x_str = "I" * site + "X" + "I" * (n_qubits - site - 1)
    y_str = "I" * site + "Y" + "I" * (n_qubits - site - 1)
    return [(x_str, 0.5 + 0.0j), (y_str, 0.0 + 0.5j)]


def sigma_minus(site: int, n_qubits: int) -> list[tuple[str, complex]]:
    """``σ⁻_q = (X_q − i Y_q) / 2``. Use as a Lindblad jump for amplitude damping."""
    if not 0 <= site < n_qubits:
        raise ValueError(f"site {site} out of range for n_qubits={n_qubits}")
    x_str = "I" * site + "X" + "I" * (n_qubits - site - 1)
    y_str = "I" * site + "Y" + "I" * (n_qubits - site - 1)
    return [(x_str, 0.5 + 0.0j), (y_str, 0.0 - 0.5j)]


def _normalize_jump(jump_op: object) -> list[tuple[str, float, float]]:
    """Convert a user-supplied jump operator to ``[(pauli_str, re, im), ...]``.

    Accepts either a single Pauli string (treated as a Hermitian-Pauli jump
    with coefficient 1) or an iterable of ``(pauli_str, complex_coeff)``
    pairs.
    """
    if isinstance(jump_op, str):
        return [(jump_op, 1.0, 0.0)]
    out: list[tuple[str, float, float]] = []
    for term in jump_op:
        s, c = term
        cc = complex(c)
        out.append((str(s), float(cc.real), float(cc.imag)))
    if not out:
        raise ValueError("jump operator lincomb must contain at least one Pauli term")
    return out


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
        Iterable of ``(jump_op, rate)`` pairs. ``jump_op`` is either a
        Pauli string ``"XYZI..."`` (treated as a Hermitian-Pauli jump
        with coefficient 1, hitting the fast path) or an iterable of
        ``(pauli_string, complex_coeff)`` pairs for a general complex
        Pauli linear combination such as :func:`sigma_plus` or
        :func:`sigma_minus`. ``rate`` is the non-negative GKSL rate
        ``γ_k``.

    Examples
    --------
    Dephasing (Hermitian Pauli):

    >>> Lindbladian(2, [("XX", 1.0)], [("ZI", 0.3), ("IZ", 0.3)])

    Amplitude damping on site 0 (non-Hermitian):

    >>> jumps = [(sigma_minus(0, 2), 0.5)]
    >>> Lindbladian(2, [("XX", 1.0)], jumps)
    """

    def __init__(
        self,
        n_qubits: int,
        h_terms: Iterable[tuple[str, float]],
        jump_terms: Iterable[tuple[object, float]] = (),
    ):
        self.n_qubits = int(n_qubits)
        h_strs: list[str] = []
        h_coeffs: list[float] = []
        for s, c in h_terms:
            h_strs.append(s)
            h_coeffs.append(float(c))
        j_lincombs: list[list[tuple[str, float, float]]] = []
        j_rates: list[float] = []
        for jump_op, rate in jump_terms:
            j_lincombs.append(_normalize_jump(jump_op))
            j_rates.append(float(rate))
        self._spec = _LindbladSpec(self.n_qubits, h_strs, h_coeffs, j_lincombs, j_rates)

    @property
    def num_h_terms(self) -> int:
        return self._spec.num_h_terms

    @property
    def num_jump_terms(self) -> int:
        return self._spec.num_jump_terms

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

    def pc_step_arr(
        self,
        basis_arr: np.ndarray,
        coeffs: np.ndarray,
        dt: float,
        max_basis: int,
        drop_tol: float = 1e-12,
        protected_arr: np.ndarray | None = None,
        num_threads: int | None = None,
        admit_basis: int | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """One predictor-corrector adaptive step.

        All work — leakage expansion, matrix-exponential step, second-hop
        re-expansion, corrector — runs in Rust; SciPy is not required.
        The matrix-exponential action is computed matrix-free via the external
        ``quspin-expm`` crate (Al-Mohy & Higham scaling-and-squaring).

        Truncation. ``max_basis`` is a hard rank cap on the live basis:
        enrichment adds at most ``max_basis - len(basis)`` of the largest
        leakage strings, and the post-step basis is trimmed to the
        top-``max_basis`` entries by ``|coeff|`` (``protected`` words always
        kept). Pass a large value (e.g. ``10_000_000``) for the near-exact,
        uncapped case. ``drop_tol`` additionally prunes basis entries whose
        absolute coefficient is below the threshold after the corrector
        (unless the word is ``protected``).

        ``num_threads``, when set, pins this call to a freshly-built rayon
        pool of that size — useful for benchmarking parallel scaling.

        ``admit_basis``, when set (must be >= ``max_basis``), bounds the
        enriched working set during the step instead of ``max_basis``: the
        step may hold up to ``admit_basis`` strings transiently, and the
        final truncation keeps the top-``max_basis`` by evolved ``|coeff|``
        over the whole union (retained + newly admitted) — genuine rank
        displacement, so no ``drop_tol`` is needed to sustain membership
        turnover. Default ``None`` reproduces the historical behaviour
        (admission bounded by ``max_basis``; turnover requires
        ``drop_tol > 0``).

        Returns ``(new_basis_arr, new_coeffs)``; the basis may have grown
        (or shrunk, if ``max_basis`` / ``drop_tol`` pruned entries).
        """
        n = self.n_qubits
        if protected_arr is None:
            protected_arr = np.zeros((0, n), dtype=np.uint8)
        return self._spec.pc_step(
            np.ascontiguousarray(basis_arr, dtype=np.uint8),
            np.ascontiguousarray(coeffs, dtype=np.float64),
            float(dt),
            int(max_basis),
            float(drop_tol),
            np.ascontiguousarray(protected_arr, dtype=np.uint8),
            None if num_threads is None else int(num_threads),
            None if admit_basis is None else int(admit_basis),
        )

    def pc_step_orbit_rep(
        self,
        basis_arr: np.ndarray,
        coeffs: np.ndarray,
        dt: float,
        max_basis: int,
        group,
        momentum: np.ndarray,
        drop_tol: float = 1e-12,
        protected_arr: np.ndarray | None = None,
        canonicalize_first: bool = False,
        admit_basis: int | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Per-step orbit-representative pc evolution.

        State lives entirely in orbit-rep form throughout: ``basis_arr``
        contains only canonical translation-orbit representatives,
        ``coeffs`` are complex. Phase-aware action + complex CSR. Basis
        is ~``|group|×`` smaller than the equivalent full-basis complex
        evolution, and the reduction persists across every step.

        Truncation. ``max_basis`` is a hard rank cap on the live orbit-rep
        basis: enrichment adds at most ``max_basis - len(basis)`` of the
        largest leakage reps, and the post-step basis is trimmed to the
        top-``max_basis`` reps by ``|c|`` (``protected`` reps always kept).
        Pass a large value (e.g. ``10_000_000``) for the near-exact,
        uncapped case. ``drop_tol`` additionally prunes reps whose absolute
        coefficient is below the threshold after the corrector.

        ``admit_basis``, when set (>= ``max_basis``), bounds the enriched
        working set instead of ``max_basis``: the step may hold up to
        ``admit_basis`` reps transiently and the final truncation keeps the
        top-``max_basis`` by evolved ``|c|`` over the whole union — the
        displacement scheme, matching the real-space ``pc_step_arr``.

        ``basis_arr`` is assumed to contain canonical reps only. Pass
        ``canonicalize_first=True`` to rewrite each row to its canonical
        rep on entry (coefficients unchanged).
        """
        n = self.n_qubits
        if protected_arr is None:
            protected_arr = np.zeros((0, n), dtype=np.uint8)
        return self._spec.pc_step_orbit_rep(
            np.ascontiguousarray(basis_arr, dtype=np.uint8),
            np.ascontiguousarray(coeffs, dtype=np.complex128),
            float(dt),
            int(max_basis),
            group,
            np.ascontiguousarray(momentum, dtype=np.int32),
            float(drop_tol),
            np.ascontiguousarray(protected_arr, dtype=np.uint8),
            bool(canonicalize_first),
            None if admit_basis is None else int(admit_basis),
        )

    def pc_step(
        self,
        basis: Sequence[str],
        coeffs: np.ndarray,
        dt: float,
        max_basis: int,
        drop_tol: float = 1e-12,
        protected: Sequence[str] | None = None,
        num_threads: int | None = None,
    ) -> tuple[list[str], np.ndarray]:
        """String-keyed variant of :meth:`pc_step_arr`."""
        n = self.n_qubits
        basis_arr = _basis_to_codes(basis, n)
        protected_arr = (
            _basis_to_codes(list(protected), n) if protected else np.zeros((0, n), dtype=np.uint8)
        )
        new_basis_arr, new_coeffs = self.pc_step_arr(
            basis_arr,
            coeffs,
            dt,
            max_basis,
            drop_tol,
            protected_arr,
            num_threads,
        )
        return _codes_to_basis(new_basis_arr), new_coeffs

    def generator_arr(
        self, basis_arr: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generator matrix as COO triples ``(rows, cols, vals)``.

        Basis given as uint8 codes. To get a SciPy sparse matrix:

        >>> import scipy.sparse as sp
        >>> rows, cols, vals = L_op.generator_arr(basis_arr)
        >>> M = sp.coo_matrix(
        ...     (vals, (rows, cols)), shape=(len(basis_arr), len(basis_arr))
        ... ).tocsc()
        """
        return self._spec.generator(np.ascontiguousarray(basis_arr, dtype=np.uint8))

    # ── String-keyed convenience API (slower; for tests / display) ──

    def action(self, p: str) -> dict[str, float]:
        """Apply ``L*`` to a single Pauli string ``p`` (string-keyed dict)."""
        codes = _string_to_codes(p, self.n_qubits)
        out_basis, out_coeffs = self._spec.action(codes)
        keys = _codes_to_basis(out_basis)
        return {k: float(v) for k, v in zip(keys, out_coeffs) if v != 0.0}

    def leakage(
        self,
        basis: Sequence[str],
        coeffs: np.ndarray,
        protected: Sequence[str] | None = None,
    ) -> dict[str, float]:
        """Off-basis leakage as a ``dict[str, float]`` (slower API)."""
        n = self.n_qubits
        basis_arr = _basis_to_codes(basis, n)
        protected_arr = (
            _basis_to_codes(list(protected), n) if protected else np.zeros((0, n), dtype=np.uint8)
        )
        out_basis, out_coeffs = self._spec.leakage(
            basis_arr,
            np.ascontiguousarray(coeffs, dtype=np.float64),
            protected_arr,
        )
        keys = _codes_to_basis(out_basis)
        return {k: float(v) for k, v in zip(keys, out_coeffs) if v != 0.0}

    def generator(
        self, basis: Sequence[str]
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generator matrix as COO triples ``(rows, cols, vals)``,
        basis given as strings. See :meth:`generator_arr` for the conversion
        to a SciPy sparse matrix."""
        n = self.n_qubits
        basis_arr = _basis_to_codes(basis, n)
        return self.generator_arr(basis_arr)
