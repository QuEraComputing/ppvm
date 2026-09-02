# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Tests for the ``TranslationGroup`` binding and ``PauliSum.symmetry_merge``.

``symmetry_merge`` is the plain real-coefficient (``k=0``) merge: every Pauli
word is replaced by its canonical translation-orbit representative and
coefficients of colliding words are summed. See ``test_momentum_merge.py``
for the phase-aware (``k != 0``) counterpart.
"""

import numpy as np
import pytest

from ppvm import PauliSum
from ppvm._core import TranslationGroup

_CODE = {"I": 0, "X": 1, "Z": 2, "Y": 3}
_CHAR = {v: k for k, v in _CODE.items()}


def codes(s):
    return np.array([_CODE[c] for c in s], dtype=np.uint8)


def string(arr):
    return "".join(_CHAR[int(c)] for c in arr)


def psum(n, terms):
    return PauliSum.new(n, terms, min_abs_coeff=0.0, max_pauli_weight=n)


# ── TranslationGroup constructors and properties ─────────────────────────────
@pytest.mark.parametrize(
    "group, n_qubits, n_generators, order",
    [
        (TranslationGroup.chain_1d(6), 6, 1, 6),
        (TranslationGroup.torus_2d(3, 2), 6, 2, 6),
        (TranslationGroup.torus_3d(2, 2, 2), 8, 3, 8),
        (TranslationGroup.ladder(3, 2), 6, 1, 3),
    ],
)
def test_group_shapes(group, n_qubits, n_generators, order):
    assert group.n_qubits == n_qubits
    assert group.n_generators == n_generators
    assert group.order == order


def test_from_generators_matches_chain_1d():
    n = 4
    shift = [(i + 1) % n for i in range(n)]
    g = TranslationGroup.from_generators(n, [shift], [n])
    ref = TranslationGroup.chain_1d(n)
    assert (g.n_qubits, g.n_generators, g.order) == (ref.n_qubits, ref.n_generators, ref.order)
    for s in ["ZIII", "IZII", "XYII", "IXYI"]:
        assert string(g.canonicalize(codes(s))) == string(ref.canonicalize(codes(s)))


@pytest.mark.parametrize(
    "perms, orders, message",
    [
        ([[1, 0, 2, 3]], [4, 4], "same length"),
        ([[1, 0, 2]], [2], "permutation length"),
        ([[1, 0, 2, 9]], [2], "out of range"),
        ([[1, 1, 2, 3]], [2], "duplicate target"),
    ],
)
def test_from_generators_validates(perms, orders, message):
    with pytest.raises(ValueError, match=message):
        TranslationGroup.from_generators(4, perms, orders)


def test_canonicalize_is_orbit_invariant():
    g = TranslationGroup.chain_1d(4)
    shifts = ["IIXY", "IXYI", "XYII", "YIIX"]
    reps = {string(g.canonicalize(codes(s))) for s in shifts}
    assert len(reps) == 1, "all cyclic shifts must share one representative"
    # The rep is itself a member of the orbit (lex-min is over the internal
    # (xbits, zbits) ordering, which isn't observable from Python).
    assert reps.pop() in shifts


def test_canonicalize_rejects_wrong_length():
    g = TranslationGroup.chain_1d(4)
    with pytest.raises(ValueError, match="length 3 but group expects 4"):
        g.canonicalize(codes("IXY"))


# ── PauliSum.symmetry_merge ──────────────────────────────────────────────────
def test_symmetry_merge_sums_one_orbit():
    """Σ_j Z_j on a 4-chain is a single free orbit: 4 entries -> 1 with c=4."""
    n = 4
    g = TranslationGroup.chain_1d(n)
    p = psum(n, [("I" * j + "Z" + "I" * (n - j - 1), 1.0) for j in range(n)])
    assert len(p.terms) == n
    p.symmetry_merge(g)
    assert len(p.terms) == 1
    (word, coeff) = p.terms[0]
    assert coeff == pytest.approx(4.0)
    assert string(g.canonicalize(codes(word))) == word


def test_symmetry_merge_keeps_distinct_orbits_and_weights():
    n = 4
    g = TranslationGroup.chain_1d(n)
    terms = [("I" * j + "Z" + "I" * (n - j - 1), 1.0) for j in range(n)]
    terms += [("I" * j + "X" + "I" * (n - j - 1), 0.25) for j in range(n)]
    p = psum(n, terms)
    p.symmetry_merge(g)
    coeffs = sorted(c for _, c in p.terms)
    assert coeffs == pytest.approx([1.0, 4.0])


def test_symmetry_merge_is_idempotent():
    """A merged sum is already in orbit-rep form, so re-merging is a no-op."""
    n = 4
    g = TranslationGroup.chain_1d(n)
    p = psum(n, [("I" * j + "Z" + "I" * (n - j - 1), 1.0) for j in range(n)])
    p.symmetry_merge(g)
    once = sorted(p.terms)
    p.symmetry_merge(g)
    assert sorted(p.terms) == once


def test_symmetry_merge_preserves_translation_invariant_trace():
    """Merging conserves Σ_p c_p, hence any orbit-summed observable."""
    n = 4
    g = TranslationGroup.chain_1d(n)
    rng = np.random.default_rng(7)
    words = ["ZIII", "IZII", "IIZI", "IIIZ", "XXII", "IXXI", "IIXX", "XIIX"]
    coeffs = rng.normal(size=len(words))
    p = psum(n, list(zip(words, coeffs, strict=True)))
    total = sum(c for _, c in p.terms)
    p.symmetry_merge(g)
    assert sum(c for _, c in p.terms) == pytest.approx(total)


def test_symmetry_merge_rejects_qubit_count_mismatch():
    p = psum(4, [("ZIII", 1.0)])
    with pytest.raises(ValueError, match="4 qubits but the TranslationGroup acts on 3"):
        p.symmetry_merge(TranslationGroup.chain_1d(3))
