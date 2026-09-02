# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Translation-symmetry merging of Pauli sums.

A `TranslationGroup` is a finite abelian group acting on qubit positions by
permutation. Every Pauli word then belongs to a translation orbit, and
dynamics that commutes with the group can be tracked using **one canonical
representative per orbit** instead of all ``|G|`` members — cutting per-step
memory and compute by up to ``|G|×`` (Teng et al., arXiv:2512.12094).

Two representations are supported:

- `ppvm.PauliSum.symmetry_merge` / `ppvm.PauliSum.momentum_merge` for the
  dictionary representation used by gate-based Trotter evolution.
- the ``*_arr`` functions here for the ``(basis_arr, coeffs)`` array
  representation used by `ppvm.Lindbladian.pc_step_arr` and
  `ppvm.Lindbladian.pc_step_orbit_rep`, where ``basis_arr`` is an
  ``(N, n_qubits)`` uint8 array with the encoding ``0=I, 1=X, 2=Z, 3=Y``.

These are thin wrappers over `ppvm._core` that coerce their arguments to the
dtypes the compiled entry points require (uint8 basis, float64 / complex128
coefficients, int32 momentum), so plain Python lists and default-dtype numpy
arrays work.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from . import _core
from ._core import TranslationGroup as TranslationGroup

__all__ = [
    "TranslationGroup",
    "canonicalize_basis_arr",
    "canonicalize_basis_arr_complex",
    "check_momentum_sector_arr",
]


def _momentum(momentum: npt.ArrayLike) -> np.ndarray:
    return np.ascontiguousarray(momentum, dtype=np.int32)


def _basis(basis_arr: npt.ArrayLike) -> np.ndarray:
    return np.ascontiguousarray(basis_arr, dtype=np.uint8)


def canonicalize_basis_arr(
    basis_arr: npt.ArrayLike,
    coeffs: npt.ArrayLike,
    group: TranslationGroup,
) -> tuple[np.ndarray, np.ndarray]:
    """Merge a real-coefficient ``(basis_arr, coeffs)`` Pauli sum into
    orbit-representative form.

    Each row of ``basis_arr`` is replaced by its canonical representative
    under ``group``; coefficients of rows collapsing to the same
    representative are **summed**. The output is no longer than the input.

    This is the trivial (``k=0``) symmetry sector. For dynamics that commutes
    with ``group`` and a ``group``-invariant initial state, it preserves every
    ``group``-invariant expectation value. Use `canonicalize_basis_arr_complex`
    for non-trivial momentum sectors.

    Args:
        basis_arr: ``(N, n_qubits)`` array of Pauli codes.
        coeffs: length-``N`` real coefficients.
        group: the symmetry group to merge under.

    Returns:
        ``(merged_basis_arr, merged_coeffs)``.
    """
    return _core.canonicalize_basis_arr(
        _basis(basis_arr),
        np.ascontiguousarray(coeffs, dtype=np.float64),
        group,
    )


def canonicalize_basis_arr_complex(
    basis_arr: npt.ArrayLike,
    coeffs: npt.ArrayLike,
    group: TranslationGroup,
    momentum: npt.ArrayLike,
) -> tuple[np.ndarray, np.ndarray]:
    """Phase-aware merge of a complex-coefficient ``(basis_arr, coeffs)``
    Pauli sum into orbit-representative form, projected onto momentum sector
    ``momentum``.

    Coefficients on each orbit's distinct members are **averaged** with the
    character weight ``χ_k(g)`` — a ``1/|orbit|`` normalization that
    `canonicalize_basis_arr` (which sums) does not apply. Orbits whose
    stabilizer cannot carry ``momentum`` project to zero and are dropped.

    If the input does not actually lie in sector ``momentum``, the projection
    silently discards the other components; call `check_momentum_sector_arr`
    first to validate.

    Args:
        basis_arr: ``(N, n_qubits)`` array of Pauli codes.
        coeffs: length-``N`` complex coefficients.
        group: the symmetry group to merge under.
        momentum: one integer mode index per group generator. The wavenumber
            along generator ``g`` is ``2π · momentum[g] / order_g``;
            ``[0, ...]`` is the trivial sector.

    Returns:
        ``(merged_basis_arr, merged_coeffs)`` with complex coefficients.
    """
    return _core.canonicalize_basis_arr_complex(
        _basis(basis_arr),
        np.ascontiguousarray(coeffs, dtype=np.complex128),
        group,
        _momentum(momentum),
    )


def check_momentum_sector_arr(
    basis_arr: npt.ArrayLike,
    coeffs: npt.ArrayLike,
    group: TranslationGroup,
    momentum: npt.ArrayLike,
    tol: float = 1e-8,
) -> None:
    """Verify that a ``(basis_arr, complex_coeffs)`` Pauli sum lies entirely
    in momentum sector ``momentum``.

    For every orbit represented in the basis, all members must satisfy
    ``c_{g·r} = χ_k(g)⁻¹ · c_r``. Orbit members absent from ``basis_arr``
    count as zero rather than being ignored, so a partially-populated orbit
    fails.

    Run this on a user-supplied initial state before feeding it to
    `canonicalize_basis_arr_complex` or
    `ppvm.Lindbladian.pc_step_orbit_rep` — silently projecting a
    wrongly-typed input throws away meaningful physics.

    Args:
        basis_arr: ``(N, n_qubits)`` array of Pauli codes.
        coeffs: length-``N`` complex coefficients.
        group: the symmetry group.
        momentum: one integer mode index per group generator.
        tol: relative tolerance on the coefficient comparison.

    Raises:
        ValueError: if the input is not in the sector, naming the offending
            orbit representative with its expected and actual coefficient.
    """
    return _core.check_momentum_sector_arr(
        _basis(basis_arr),
        np.ascontiguousarray(coeffs, dtype=np.complex128),
        group,
        _momentum(momentum),
        tol,
    )
