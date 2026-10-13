# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Coercion of array-like Pauli codes and momentum modes to the dtypes the core
expects."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt


def as_pauli_codes(codes: npt.ArrayLike) -> np.ndarray:
    """Return `codes` as a C-contiguous uint8 array of Pauli codes
    (``0=I, 1=X, 2=Z, 3=Y``).

    uint8 input passes through unchanged (the core range-checks it). Any
    other input must be integer-typed with every value in ``0..=3``, so an
    out-of-range code raises ``ValueError`` instead of wrapping on the cast
    (e.g. ``256`` silently becoming ``0 = I``).
    """
    arr = np.asarray(codes)
    if arr.dtype != np.uint8 and arr.size:
        if not np.issubdtype(arr.dtype, np.integer):
            raise ValueError(f"Pauli codes must be integers, got dtype {arr.dtype}")
        lo, hi = arr.min(), arr.max()
        if lo < 0 or hi > 3:
            bad = lo if lo < 0 else hi
            raise ValueError(f"Pauli code must be 0 (I), 1 (X), 2 (Z), or 3 (Y); got {bad}")
    return np.ascontiguousarray(arr, dtype=np.uint8)


def as_momentum(momentum: npt.ArrayLike) -> np.ndarray:
    """Return `momentum` as a C-contiguous int32 array of mode indices.

    Input must be integer-typed with every value in the int32 range, so a
    fractional mode (``0.5``) or an out-of-range one raises ``ValueError``
    instead of being truncated or wrapped into a different sector.
    """
    arr = np.asarray(momentum)
    if arr.dtype != np.int32 and arr.size:
        if not np.issubdtype(arr.dtype, np.integer):
            raise ValueError(f"momentum modes must be integers, got dtype {arr.dtype}")
        info = np.iinfo(np.int32)
        lo, hi = arr.min(), arr.max()
        if lo < info.min or hi > info.max:
            bad = lo if lo < info.min else hi
            raise ValueError(f"momentum mode {bad} is out of int32 range")
    return np.ascontiguousarray(arr, dtype=np.int32)
