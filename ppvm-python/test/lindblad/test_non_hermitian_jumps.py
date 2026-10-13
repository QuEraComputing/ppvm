# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Non-Hermitian dissipators: cross-check :meth:`Lindbladian.action` /
:meth:`Lindbladian.generator` / :meth:`Lindbladian.leakage` against the dense
2^L × 2^L Liouvillian reference. Only viable for L ≤ 3.
"""

from __future__ import annotations

import numpy as np

from ppvm import Lindbladian, sigma_minus, sigma_plus

from ._helpers import (
    SIGMA_MINUS_MAT,
    SIGMA_PLUS_MAT,
    all_strings,
    coo_to_dense,
    dense_action,
    embed_op,
    pauli_mat,
    random_pauli_str,
)


def test_amplitude_damping_action():
    L = 3
    gamma = 0.5
    h_terms = [("XXI", 1.0), ("IXX", 0.7), ("ZIZ", 0.3)]
    jump_terms = [(sigma_minus(i, L), gamma) for i in range(L)]
    L_op = Lindbladian(L, h_terms, jump_terms)

    H = sum(c * pauli_mat(s) for s, c in h_terms)
    jumps_dense = [(embed_op(SIGMA_MINUS_MAT, i, L), gamma) for i in range(L)]

    rng = np.random.default_rng(11)
    strings = ["III", "ZII", "IZI", "IIZ", "XYZ", "YYZ"]
    strings += [random_pauli_str(rng, L) for _ in range(10)]

    for p in strings:
        got = L_op.action(p)
        want = dense_action(H, jumps_dense, p, L)
        for k in set(got) | set(want):
            diff = abs(got.get(k, 0.0) - want.get(k, 0.0))
            assert diff < 1e-10, (
                f"sigma_minus action mismatch at p={p!r} k={k!r}: "
                f"shim={got.get(k, 0.0)} ref={want.get(k, 0.0)} diff={diff}"
            )


def test_thermal_excitation_damping_action():
    """σ⁺ + σ⁻ jumps together (thermal bath at finite temperature)."""
    L = 2
    h_terms = [("XX", 1.0), ("ZI", 0.2), ("IZ", 0.1)]
    jump_terms = [(sigma_minus(i, L), 0.4) for i in range(L)] + [
        (sigma_plus(i, L), 0.1) for i in range(L)
    ]
    L_op = Lindbladian(L, h_terms, jump_terms)

    H = sum(c * pauli_mat(s) for s, c in h_terms)
    jumps_dense = [(embed_op(SIGMA_MINUS_MAT, i, L), 0.4) for i in range(L)] + [
        (embed_op(SIGMA_PLUS_MAT, i, L), 0.1) for i in range(L)
    ]

    for p in all_strings(L):
        got = L_op.action(p)
        want = dense_action(H, jumps_dense, p, L)
        for k in set(got) | set(want):
            diff = abs(got.get(k, 0.0) - want.get(k, 0.0))
            assert diff < 1e-10, (
                f"sigma_plus/sigma_minus action mismatch at p={p!r} k={k!r}: "
                f"shim={got.get(k, 0.0)} ref={want.get(k, 0.0)}"
            )


def test_amplitude_damping_generator_and_leakage():
    L = 3
    gamma = 0.3
    h_terms = [("XXI", 0.5), ("IXX", 0.5), ("ZII", 0.2), ("IZI", 0.2), ("IIZ", 0.2)]
    jump_terms = [(sigma_minus(0, L), gamma), (sigma_minus(2, L), gamma)]
    L_op = Lindbladian(L, h_terms, jump_terms)

    basis = ["III", "ZII", "IZI", "IIZ", "ZZI", "IZZ"]
    coeffs = np.array([0.1, 0.5, -0.3, 0.4, 0.2, -0.1])

    H = sum(c * pauli_mat(s) for s, c in h_terms)
    jumps_dense = [(embed_op(SIGMA_MINUS_MAT, i, L), gamma) for i in (0, 2)]

    M_shim = coo_to_dense(L_op.generator(basis), len(basis))
    M_ref = np.zeros((len(basis), len(basis)))
    idx = {p: i for i, p in enumerate(basis)}
    leak_ref: dict = {}
    for col, p in enumerate(basis):
        action_p = dense_action(H, jumps_dense, p, L)
        for q, v in action_p.items():
            if q in idx:
                M_ref[idx[q], col] += v
            else:
                leak_ref[q] = leak_ref.get(q, 0.0) + v * coeffs[col]
    assert np.max(np.abs(M_shim - M_ref)) < 1e-10

    leak_shim = L_op.leakage(basis, coeffs)
    leak_ref = {k: v for k, v in leak_ref.items() if abs(v) > 1e-14}
    for k in set(leak_shim) | set(leak_ref):
        diff = abs(leak_shim.get(k, 0.0) - leak_ref.get(k, 0.0))
        assert diff < 1e-10, f"leakage mismatch at k={k!r}: diff={diff}"
