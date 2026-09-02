# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Orbit-representative predictor-corrector evolution through the Python
binding (:meth:`Lindbladian.pc_step_orbit_rep`).

The state lives entirely in orbit-rep form: the basis holds only canonical
translation-orbit representatives and the coefficients are complex. For a
translation-invariant Lindbladian, evolving in orbit-rep form and projecting
onto the momentum sector at the *end* of a full-space evolution must agree —
that is the content of the projection theorem, and it is what the first test
checks against a dense numpy matrix exponential, independent of the Rust
Al-Mohy & Higham implementation.
"""

from __future__ import annotations

import cmath

import numpy as np
import pytest

from ppvm import Lindbladian, TranslationGroup, canonicalize_basis_arr_complex

from ._helpers import all_strings

_CODE = {"I": 0, "X": 1, "Z": 2, "Y": 3}
_CHAR = {v: k for k, v in _CODE.items()}


def basis_arr(strings, n):
    arr = np.zeros((len(strings), n), dtype=np.uint8)
    for i, s in enumerate(strings):
        arr[i] = [_CODE[c] for c in s]
    return arr


def string(row):
    return "".join(_CHAR[int(c)] for c in row)


def to_dict(basis, coeffs):
    return {string(w): c for w, c in zip(basis, coeffs, strict=True)}


def momentum(*modes):
    """A plain tuple: both the wrappers and `Lindbladian.pc_step_orbit_rep`
    coerce momentum to the int32 the compiled code needs."""
    return modes


def xy_chain_pbc(n, gamma):
    """Translation-invariant XY chain with PBC plus uniform Z dephasing."""
    h_terms = []
    for j in range(n):
        nxt = (j + 1) % n
        for op in "XY":
            s = ["I"] * n
            s[j] = op
            s[nxt] = op
            h_terms.append(("".join(s), 1.0))
    jumps = [("I" * j + "Z" + "I" * (n - j - 1), gamma) for j in range(n)]
    return Lindbladian(n, h_terms, jumps)


def z_momentum_seed(n, k):
    """``O_k = Σ_a e^{-2πi k a / n} Z_a`` as ``(basis_arr, complex coeffs)``."""
    words = ["I" * a + "Z" + "I" * (n - a - 1) for a in range(n)]
    coeffs = np.array([cmath.exp(-2j * cmath.pi * k * a / n) for a in range(n)])
    return basis_arr(words, n), coeffs


def _dense_expm(A, terms=40):
    """``exp(A)`` for a real matrix, by Taylor series with scaling and squaring.

    The full-space Lindbladian is not diagonalizable in general (numpy's
    ``eig`` returns a singular eigenvector matrix here), and the test suite
    deliberately has no scipy dependency — so the reference is this direct
    series, independent of the Rust Al-Mohy & Higham implementation.
    """
    norm = np.abs(A).sum(axis=0).max()
    s = max(0, int(np.ceil(np.log2(norm))) + 1) if norm > 0 else 0
    B = A / 2**s
    total = np.eye(A.shape[0])
    term = np.eye(A.shape[0])
    for k in range(1, terms + 1):
        term = term @ B / k
        total = total + term
    for _ in range(s):
        total = total @ total
    return total


@pytest.mark.parametrize("k", [0, 1, 2])
def test_orbit_rep_matches_dense_full_space_then_project(k):
    """Orbit-rep evolution == full-space evolution projected at the end.

    Full space is all 4^n Pauli strings (n=3 -> 64), exponentiated densely
    with numpy. The orbit-rep side runs with a huge ``max_basis`` so its rank
    cap never binds and the only remaining difference would be a bug in the
    phase-aware action.
    """
    n = 3
    dt = 0.02
    n_steps = 3
    op = xy_chain_pbc(n, gamma=0.3)
    group = TranslationGroup.chain_1d(n)
    k_arr = momentum(k)

    # --- dense full-space reference ---
    full = all_strings(n)
    generator = np.zeros((len(full), len(full)), dtype=float)
    rows, cols, vals = op.generator(full)
    generator[rows, cols] = vals
    seed_basis, seed_coeffs = z_momentum_seed(n, k)
    index = {s: i for i, s in enumerate(full)}
    v = np.zeros(len(full), dtype=complex)
    for w, c in zip(seed_basis, seed_coeffs, strict=True):
        v[index[string(w)]] = c
    # `generator` is real, so exp(dt·G) is too: evolve the real and imaginary
    # parts of the coefficient vector separately.
    step = _dense_expm(dt * generator)
    v_re, v_im = v.real.copy(), v.imag.copy()
    for _ in range(n_steps):
        v_re = step @ v_re
        v_im = step @ v_im
    v = v_re + 1j * v_im
    expected = to_dict(*canonicalize_basis_arr_complex(basis_arr(full, n), v, group, k_arr))

    # --- orbit-rep evolution ---
    rep_basis, rep_coeffs = canonicalize_basis_arr_complex(seed_basis, seed_coeffs, group, k_arr)
    for _ in range(n_steps):
        rep_basis, rep_coeffs = op.pc_step_orbit_rep(
            rep_basis, rep_coeffs, dt, 10_000_000, group, k_arr, drop_tol=0.0
        )
    got = to_dict(rep_basis, rep_coeffs)

    # Compare on the union; the dense side keeps exact zeros the orbit-rep
    # side never admits, so only nonzero entries must match.
    for word in set(expected) | set(got):
        e = expected.get(word, 0.0)
        g = got.get(word, 0.0)
        assert abs(e - g) < 1e-9, f"k={k}: rep {word} dense {e} vs orbit-rep {g}"
    assert any(abs(c) > 1e-6 for c in got.values()), "orbit-rep state decayed away"


def test_canonicalize_first_accepts_non_canonical_input():
    """The same physical state seeded on a non-canonical orbit member gives
    the same evolution once ``canonicalize_first=True`` normalizes it."""
    n = 3
    dt = 0.02
    op = xy_chain_pbc(n, gamma=0.3)
    group = TranslationGroup.chain_1d(n)
    k_arr = momentum(0)

    seed_basis, seed_coeffs = z_momentum_seed(n, 0)
    canonical, coeffs = canonicalize_basis_arr_complex(seed_basis, seed_coeffs, group, k_arr)

    ref_basis, ref_coeffs = op.pc_step_orbit_rep(canonical, coeffs, dt, 10_000_000, group, k_arr)
    # Feed a shifted (non-canonical) representative of the same orbit.
    shifted = np.array([[_CODE[c] for c in "IZI"]], dtype=np.uint8)
    got_basis, got_coeffs = op.pc_step_orbit_rep(
        shifted, coeffs, dt, 10_000_000, group, k_arr, canonicalize_first=True
    )
    ref = to_dict(ref_basis, ref_coeffs)
    got = to_dict(got_basis, got_coeffs)
    assert ref.keys() == got.keys()
    for w in ref:
        assert abs(ref[w] - got[w]) < 1e-12


def test_max_basis_caps_the_live_basis():
    n = 4
    op = xy_chain_pbc(n, gamma=0.1)
    group = TranslationGroup.chain_1d(n)
    k_arr = momentum(0)
    seed_basis, seed_coeffs = z_momentum_seed(n, 0)
    basis, coeffs = canonicalize_basis_arr_complex(seed_basis, seed_coeffs, group, k_arr)
    for _ in range(4):
        basis, coeffs = op.pc_step_orbit_rep(basis, coeffs, 0.05, 6, group, k_arr)
        assert basis.shape[0] <= 6
        assert basis.shape == (len(coeffs), n)


def test_protected_reps_are_never_dropped():
    n = 4
    op = xy_chain_pbc(n, gamma=0.1)
    group = TranslationGroup.chain_1d(n)
    k_arr = momentum(0)
    seed_basis, seed_coeffs = z_momentum_seed(n, 0)
    basis, coeffs = canonicalize_basis_arr_complex(seed_basis, seed_coeffs, group, k_arr)
    protected = basis.copy()
    keep = {string(w) for w in protected}
    for _ in range(3):
        # max_basis=1 with a huge drop_tol would wipe everything unprotected.
        basis, coeffs = op.pc_step_orbit_rep(
            basis, coeffs, 0.05, 1, group, k_arr, drop_tol=1e3, protected_arr=protected
        )
        assert keep <= {string(w) for w in basis}


def test_pc_step_orbit_rep_validates_inputs():
    n = 3
    op = xy_chain_pbc(n, gamma=0.0)
    group = TranslationGroup.chain_1d(n)
    basis, coeffs = z_momentum_seed(n, 0)
    with pytest.raises(ValueError, match="momentum has 2 entries but group has 1 generators"):
        op.pc_step_orbit_rep(basis, coeffs, 0.01, 100, group, momentum(0, 0))
    with pytest.raises(ValueError, match="coeffs has length 2 but basis has 3 rows"):
        op.pc_step_orbit_rep(basis, coeffs[:2], 0.01, 100, group, momentum(0))


def test_returns_complex_arrays_of_matching_shape():
    n = 3
    op = xy_chain_pbc(n, gamma=0.2)
    group = TranslationGroup.chain_1d(n)
    k_arr = momentum(1)
    seed_basis, seed_coeffs = z_momentum_seed(n, 1)
    basis, coeffs = canonicalize_basis_arr_complex(seed_basis, seed_coeffs, group, k_arr)
    out_basis, out_coeffs = op.pc_step_orbit_rep(basis, coeffs, 0.01, 500, group, k_arr)
    assert out_basis.dtype == np.uint8
    assert out_coeffs.dtype == np.complex128
    assert out_basis.shape == (len(out_coeffs), n)
