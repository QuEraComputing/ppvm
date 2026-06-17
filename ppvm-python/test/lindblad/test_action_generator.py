# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Cross-checks for :meth:`Lindbladian.action` / :meth:`Lindbladian.generator`
/ :meth:`Lindbladian.leakage` against the Hermitian-Pauli reference built from
the single-qubit Pauli multiplication table.
"""

from __future__ import annotations

import numpy as np
import pytest

from ppvm import Lindbladian

from ._helpers import (
    _reference_action,
    _tuple_to_str,
    assert_action_matches,
    coo_to_dense,
    expm_mv_dense,
    random_pauli_str,
    tfim_xdeph,
    xy_dephasing,
)


def test_action_xy_dephasing():
    L = 8
    h_terms, jump_terms = xy_dephasing(L, alpha=1.0, gamma=0.3)
    L_op = Lindbladian(L, h_terms, jump_terms)
    rng = np.random.default_rng(42)
    strings = ["I" * L]
    strings += ["I" * i + "Z" + "I" * (L - i - 1) for i in range(L)]
    strings += [random_pauli_str(rng, L) for _ in range(20)]
    assert_action_matches(L_op, h_terms, jump_terms, strings)


def test_action_tfim_xdephasing():
    L = 6
    h_terms, jump_terms = tfim_xdeph(L, J=0.7, h=0.4, gamma=0.2)
    L_op = Lindbladian(L, h_terms, jump_terms)
    rng = np.random.default_rng(7)
    strings = [random_pauli_str(rng, L) for _ in range(30)]
    assert_action_matches(L_op, h_terms, jump_terms, strings)


def test_generator_leakage_and_expm():
    L = 5
    dt = 0.1
    h_terms, jump_terms = xy_dephasing(L, alpha=1.0, gamma=0.4)
    L_op = Lindbladian(L, h_terms, jump_terms)
    basis = ["I" * i + "Z" + "I" * (L - i - 1) for i in range(L)]
    basis += ["YIYII", "IYIYI", "ZZIII"]
    coeffs = np.array([0.5, -0.3, 0.2, 0.4, -0.1, 0.6, 0.2, 0.1])

    M_shim = coo_to_dense(L_op.generator(basis), len(basis))
    M_ref = np.zeros((len(basis), len(basis)))
    index = {p: i for i, p in enumerate(basis)}
    for col, p in enumerate(basis):
        for r_tuple, v in _reference_action(p, h_terms, jump_terms).items():
            r = _tuple_to_str(r_tuple)
            if r in index:
                M_ref[index[r], col] += v
    assert np.max(np.abs(M_shim - M_ref)) < 1e-12

    shim_leak = L_op.leakage(basis, coeffs)
    ref_leak: dict = {}
    for p, cf in zip(basis, coeffs):
        for r_tuple, v in _reference_action(p, h_terms, jump_terms).items():
            r = _tuple_to_str(r_tuple)
            if r not in index:
                ref_leak[r] = ref_leak.get(r, 0.0) + v * cf
    ref_leak = {kk: v for kk, v in ref_leak.items() if v}
    for kk in set(shim_leak) | set(ref_leak):
        assert abs(shim_leak.get(kk, 0.0) - ref_leak.get(kk, 0.0)) < 1e-12

    c_shim = expm_mv_dense(dt * M_shim, coeffs)
    c_ref = expm_mv_dense(dt * M_ref, coeffs)
    assert np.allclose(c_shim, c_ref, atol=1e-13)


def test_generator_rejects_duplicate_basis():
    """Duplicate basis rows would silently overwrite each other in the
    row-index map and produce an incorrect sparse generator. The user-facing
    entry point must reject them with a clear ValueError instead.
    """
    L = 4
    h_terms, jump_terms = xy_dephasing(L, alpha=1.0, gamma=0.3)
    L_op = Lindbladian(L, h_terms, jump_terms)
    basis = ["ZIII", "IZII", "ZIII"]  # duplicate at rows 0 and 2
    with pytest.raises(ValueError, match=r"duplicate Pauli word at row 0 and row 2"):
        L_op.generator(basis)
    # pc_step also builds the row index and must reject duplicates.
    with pytest.raises(ValueError, match=r"duplicate Pauli word"):
        L_op.pc_step(basis, np.ones(len(basis)), 0.01, 1e-10)


def test_protected_strings_suppressed():
    L = 4
    h_terms, jump_terms = xy_dephasing(L, alpha=1.0, gamma=0.0)
    L_op = Lindbladian(L, h_terms, jump_terms)
    basis = ["ZIII"]
    coeffs = np.array([1.0])
    leak = L_op.leakage(basis, coeffs)
    assert leak, "expected some leakage"
    protected_key = next(iter(leak))
    leak2 = L_op.leakage(basis, coeffs, protected=[protected_key])
    assert protected_key not in leak2
