# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Tests for momentum-sector (k != 0) symmetry merging of real PauliSum pairs.

A complex operator O = O_re + i·O_im is carried as a pair of real PauliSums.
``PauliSum.momentum_merge`` folds the pair onto translation-orbit
representatives in momentum sector k, generalizing ``symmetry_merge`` (k=0).

These checks compare against *exact* references — the projector definition,
idempotency, and exact diagonalization of the dynamics — NOT against any
other propagation scheme.
"""

import cmath
import math

import numpy as np
import pytest

from ppvm import PauliSum
from ppvm._core import TranslationGroup

# ── dense Pauli helpers (exact references) ───────────────────────────────────
_I = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_P = {"I": _I, "X": _X, "Y": _Y, "Z": _Z}


def dense(pauli_str):
    m = np.array([[1]], dtype=complex)
    for ch in pauli_str:
        m = np.kron(m, _P[ch])
    return m


def zstr(n, q):
    return "".join("Z" if i == q else "I" for i in range(n))


def chain_bonds(n):
    return [(i, (i + 1) % n, 1.0) for i in range(n)]


# ── helpers shared with the k-resolved Trotter driver ────────────────────────
def _seed_pair(n, k):
    a = np.arange(n)
    re = np.cos(2 * np.pi * k * a / n)
    im = -np.sin(2 * np.pi * k * a / n)  # e^{-2πi k a/n} = cos - i sin
    Z = [zstr(n, q) for q in range(n)]
    PA = PauliSum.new(
        n, [(Z[q], float(re[q])) for q in range(n)], min_abs_coeff=0.0, max_pauli_weight=n
    )
    PB = PauliSum.new(
        n, [(Z[q], float(im[q])) for q in range(n)], min_abs_coeff=0.0, max_pauli_weight=n
    )
    return PA, PB


def _to_complex_dict(PA, PB):
    d = {}
    for s, c in PA.terms:
        d[s] = d.get(s, 0j) + c
    for s, c in PB.terms:
        d[s] = d.get(s, 0j) + 1j * c
    return {s: v for s, v in d.items() if v != 0j}


def _ovl(sA, sB, oA, oB):
    re = sA.overlap(oA) + sB.overlap(oB)
    im = sA.overlap(oB) - sB.overlap(oA)
    return complex(re, im)


# =============================================================================
# 1. The merge is an exact sector projector: idempotent, and it leaves a
#    genuine momentum-k eigenoperator unchanged.
# =============================================================================
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_momentum_merge_idempotent(k):
    n = 4
    g = TranslationGroup.chain_1d(n)
    PA, PB = _seed_pair(n, k)  # S^z_k is exactly in sector k
    PA.momentum_merge(PB, g, [k])
    once = _to_complex_dict(PA, PB)
    PA.momentum_merge(PB, g, [k])  # merging again must be a no-op
    twice = _to_complex_dict(PA, PB)
    keys = set(once) | set(twice)
    assert max(abs(once.get(x, 0j) - twice.get(x, 0j)) for x in keys) < 1e-12


def test_momentum_merge_projects_out_other_sectors():
    """Merging a pure sector-k operator in sector k' != k gives ~zero."""
    n = 4
    g = TranslationGroup.chain_1d(n)
    PA, PB = _seed_pair(n, 1)  # operator lives in k=1
    PA.momentum_merge(PB, g, [2])  # project onto k=2
    d = _to_complex_dict(PA, PB)
    assert all(abs(v) < 1e-12 for v in d.values()), d


# =============================================================================
# 2. End-to-end: k-resolved, symmetry-compressed Trotter reproduces the
#    EXACT (dense-diagonalization) operator autocorrelator as dt -> 0.
# =============================================================================
def _ed_autocorr(n, bonds, k, ts):
    """C_k(t) = Tr[O0^dagger O(t)] / Tr[O0^dagger O0], O0 = S^z_k, exact."""
    H = np.zeros((2**n, 2**n), dtype=complex)
    for i, j, J in bonds:
        for q in "XY":
            s = ["I"] * n
            s[i] = q
            s[j] = q
            H += J * dense("".join(s))
    O0 = np.zeros((2**n, 2**n), dtype=complex)
    for a in range(n):
        O0 += cmath.exp(-2j * math.pi * k * a / n) * dense(zstr(n, a))
    E, V = np.linalg.eigh(H)
    out = []
    with np.errstate(all="ignore"):  # silence spurious macOS-Accelerate matmul warnings
        norm = np.trace(O0.conj().T @ O0).real
        for t in ts:
            U = (V * np.exp(-1j * E * t)) @ V.conj().T
            Ot = U.conj().T @ O0 @ U
            out.append(np.trace(O0.conj().T @ Ot) / norm)
    return np.array(out)


def _ctrotter_autocorr(n, bonds, k, dt, steps):
    g = TranslationGroup.chain_1d(n)
    PA, PB = _seed_pair(n, k)
    PA.momentum_merge(PB, g, [k])
    refA, refB = PA.copy(), PB.copy()
    C0 = _ovl(refA, refB, PA, PB)
    out = [1.0 + 0j]
    for _ in range(steps):
        for i, j, J in bonds:  # Strang: forward then reversed
            PA.rxx(i, j, J * dt, truncate=False)
            PA.ryy(i, j, J * dt, truncate=False)
            PB.rxx(i, j, J * dt, truncate=False)
            PB.ryy(i, j, J * dt, truncate=False)
        for i, j, J in reversed(bonds):
            PA.rxx(i, j, J * dt, truncate=False)
            PA.ryy(i, j, J * dt, truncate=False)
            PB.rxx(i, j, J * dt, truncate=False)
            PB.ryy(i, j, J * dt, truncate=False)
        PA.momentum_merge(PB, g, [k])
        out.append(_ovl(refA, refB, PA, PB) / C0)
    return np.array(out)


@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_k_resolved_trotter_converges_to_exact(k):
    n, T = 4, 0.3
    bonds = chain_bonds(n)
    # exact reference at the matching times for two step sizes
    err = {}
    for dt in (0.04, 0.02):
        steps = round(T / dt)
        ts = np.arange(steps + 1) * dt
        c = _ctrotter_autocorr(n, bonds, k, dt, steps)
        ed = _ed_autocorr(n, bonds, k, ts)
        err[dt] = np.max(np.abs(c - ed))

    assert abs(_ctrotter_autocorr(n, bonds, k, 0.02, 1)[0] - 1.0) < 1e-12  # C_k(0)=1
    if k == 0:
        # total Z is conserved -> exact in every sector-0 step
        assert err[0.02] < 1e-10
    else:
        assert err[0.02] < 5e-3  # close to exact at dt=0.02
        assert err[0.02] < err[0.04]  # converges toward exact as dt->0


def test_compressed_matches_uncompressed_evolution():
    """Merging must not change observables beyond the O(dt^2) equivariance
    error: compressed (merge each step) vs the same gates with no merge."""
    n, k, dt, steps = 4, 2, 0.02, 10
    bonds = chain_bonds(n)
    g = TranslationGroup.chain_1d(n)

    # uncompressed: evolve the full real pair, project only at readout
    PA, PB = _seed_pair(n, k)
    rA, rB = _seed_pair(n, k)
    rA.momentum_merge(rB, g, [k])
    C0 = _ovl(rA, rB, *_merged_copy(PA, PB, g, k))
    comp = _ctrotter_autocorr(n, bonds, k, dt, steps)
    unc = []
    for _ in range(steps):
        for i, j, J in bonds:
            PA.rxx(i, j, J * dt, truncate=False)
            PA.ryy(i, j, J * dt, truncate=False)
            PB.rxx(i, j, J * dt, truncate=False)
            PB.ryy(i, j, J * dt, truncate=False)
        for i, j, J in reversed(bonds):
            PA.rxx(i, j, J * dt, truncate=False)
            PA.ryy(i, j, J * dt, truncate=False)
            PB.rxx(i, j, J * dt, truncate=False)
            PB.ryy(i, j, J * dt, truncate=False)
        mA, mB = _merged_copy(PA, PB, g, k)
        unc.append(_ovl(rA, rB, mA, mB) / C0)
    unc = np.array([1.0 + 0j, *unc])
    assert np.max(np.abs(comp - unc)) < 5e-3  # only O(dt^2) equivariance


def _merged_copy(PA, PB, g, k):
    a, b = PA.copy(), PB.copy()
    a.momentum_merge(b, g, [k])
    return a, b
