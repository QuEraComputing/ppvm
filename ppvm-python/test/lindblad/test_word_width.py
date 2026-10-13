# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Registers wider than 128 qubits (variable Pauli-word width)."""

import numpy as np
import pytest

from ppvm import Lindbladian
from ppvm._core import TranslationGroup, canonicalize_basis_arr_complex
from ppvm.lindblad import _basis_to_codes


def _word(n, ops):
    s = ["I"] * n
    for q, p in ops:
        s[q] = p
    return "".join(s)


def _ring(n):
    terms = []
    for i in range(n):
        j = (i + 1) % n
        terms += [
            (_word(n, [(i, "X"), (j, "X")]), 1.0),
            (_word(n, [(i, "Y"), (j, "Y")]), 1.0),
            (_word(n, [(i, "Z"), (j, "Z")]), 0.5),
            (_word(n, [(i, "Z")]), 0.3),
        ]
    return terms


def _orbit_run(n, steps=2):
    lind = Lindbladian(n, _ring(n), [])
    group = TranslationGroup.chain_1d(n)
    mom = np.array([0], dtype=np.int32)
    seed = _basis_to_codes([_word(n, [(q, "X")]) for q in range(n)], n)
    basis, co = canonicalize_basis_arr_complex(seed, np.ones(n, dtype=np.complex128), group, mom)
    for _ in range(steps):
        basis, co = lind.pc_step_orbit_rep(
            basis, co, 0.1, 10**6, group=group, momentum=mom, drop_tol=0.0
        )
    return basis, co


@pytest.mark.parametrize("n", [130, 256, 300, 512])
def test_wide_lindbladian_constructs_and_steps(n):
    # 130 qubits used to fail with "LindbladSpec supports n_qubits ≤ 128".
    lind = Lindbladian(
        n,
        [(_word(n, [(0, "Z"), (n - 1, "Z")]), 0.7), (_word(n, [(n - 1, "X")]), 1.3)],
        [(_word(n, [(n - 1, "Z")]), 0.05)],
    )
    assert lind.n_qubits == n
    basis, _ = lind.pc_step([_word(n, [(n - 1, "Z")])], np.array([1.0]), 0.05, 1000)
    assert len(basis) > 1
    assert all(len(b) == n for b in basis)
    # The top qubit is live: something acts on it.
    assert any(b[n - 1] == "Y" for b in basis)


def test_more_than_512_qubits_is_rejected():
    with pytest.raises(ValueError, match="512"):
        Lindbladian(513, [(_word(513, [(0, "Z")]), 1.0)], [])


def test_orbit_rep_step_beyond_128_qubits():
    basis, co = _orbit_run(130)
    assert basis.shape[1] == 130
    assert len(basis) > 5
    assert np.all(np.isfinite(co))


def test_orbit_rep_step_matches_across_the_width_boundary():
    # Rings of 120 and 136 sites use 128- and 256-qubit words. Compare the
    # k=0 autocorrelation of M_x, which is ring-size independent here.
    def autocorr(n):
        basis, co = _orbit_run(n, steps=3)
        # The single-X rep, whatever position the canonical form puts it at.
        m = np.where(((basis == 1).sum(axis=1) == 1) & ((basis != 0).sum(axis=1) == 1))[0]
        assert m.size == 1
        return co[m[0]].real

    # With nearest-neighbour H and 3 short steps the operator front is far
    # smaller than either ring, so the value is ring-size independent.
    assert autocorr(120) == pytest.approx(autocorr(136), rel=1e-12, abs=1e-14)
