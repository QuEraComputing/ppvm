# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""In-place Hamiltonian-coefficient setter: ``Lindbladian.set_h_coeffs``.

Contract: after ``L.set_h_coeffs(c)`` every output of ``L`` is bit-identical
(``np.array_equal``) to a freshly constructed ``Lindbladian`` with
coefficients ``c``; jump / Kossakowski data are untouched.
"""

import itertools
import time

import numpy as np
import pytest

from ppvm import Lindbladian, sigma_minus, sigma_plus
from ppvm.lindblad import _string_to_codes


def _local_terms(n: int, max_weight: int) -> list[str]:
    """All Pauli strings on n qubits of weight 1..max_weight (nearest-neighbour for weight 2)."""
    out = []
    for q in range(n):
        for p in "XYZ":
            s = ["I"] * n
            s[q] = p
            out.append("".join(s))
    if max_weight >= 2:
        for q in range(n - 1):
            for a, b in itertools.product("XYZ", repeat=2):
                s = ["I"] * n
                s[q], s[q + 1] = a, b
                out.append("".join(s))
    return out


def _build(n, words, coeffs, rng_seed):
    rng = np.random.default_rng(rng_seed)
    jumps = []
    # Hermitian single-Pauli jumps (fast path).
    for q in range(n):
        s = ["I"] * n
        s[q] = "Z"
        jumps.append(("".join(s), float(rng.uniform(0.01, 0.2))))
    # Non-Hermitian jumps.
    jumps.append((sigma_minus(0, n), float(rng.uniform(0.05, 0.3))))
    jumps.append((sigma_plus(n - 1, n), float(rng.uniform(0.05, 0.3))))
    # Multi-term complex lincomb jump.
    xs = ["I"] * n
    xs[1] = "X"
    ys = ["I"] * n
    ys[2] = "Y"
    jumps.append(([("".join(xs), 0.7 + 0.1j), ("".join(ys), -0.3j)], 0.11))
    # Kossakowski block on sigma_minus ops with a random PSD K.
    ops = [sigma_minus(j, n) for j in range(n)]
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    K = 0.05 * (a @ a.conj().T) / n
    h_terms = list(zip(words, [float(c) for c in coeffs]))
    return Lindbladian(n, h_terms, jumps, kossakowski=(ops, K))


def _assert_same_tuple(a, b):
    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert x.dtype == y.dtype
        assert np.array_equal(x, y)


def _start_basis(n):
    words = []
    for q in range(n):
        s = ["I"] * n
        s[q] = "Z"
        words.append("".join(s))
    s = ["I"] * n
    s[0] = "X"
    words.append("".join(s))
    basis = np.stack([_string_to_codes(w, n) for w in words])
    coeffs = np.linspace(1.0, 0.3, len(words))
    return basis, coeffs


def _compare_all(L_set, L_fresh, n, rng):
    words_all = _local_terms(n, 2)
    # action on several strings
    for w in rng.choice(words_all, size=8, replace=False):
        p = _string_to_codes(str(w), n)
        _assert_same_tuple(L_set.action_arr(p), L_fresh.action_arr(p))
    # leakage + generator on a random basis
    basis_words = list(rng.choice(words_all, size=min(12, len(words_all)), replace=False))
    basis = np.stack([_string_to_codes(str(w), n) for w in basis_words])
    c = rng.normal(size=len(basis_words))
    _assert_same_tuple(L_set.leakage_arr(basis, c), L_fresh.leakage_arr(basis, c))
    prot = basis[:2]
    _assert_same_tuple(L_set.leakage_arr(basis, c, prot), L_fresh.leakage_arr(basis, c, prot))
    _assert_same_tuple(L_set.generator_arr(basis), L_fresh.generator_arr(basis))
    # several pc_steps with truncation, two truncation modes
    for kwargs in (
        {"max_basis": 24, "drop_tol": 1e-10},
        {"max_basis": 16, "drop_tol": 0.0, "admit_basis": 40},
    ):
        b1, c1 = _start_basis(n)
        b2, c2 = b1.copy(), c1.copy()
        for _ in range(6):
            b1, c1 = L_set.pc_step_arr(b1, c1, 0.05, **kwargs)
            b2, c2 = L_fresh.pc_step_arr(b2, c2, 0.05, **kwargs)
            assert np.array_equal(b1, b2)
            assert np.array_equal(c1, c2)
        assert len(c1) <= kwargs["max_basis"]


@pytest.mark.parametrize("n", [3, 5])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_set_h_coeffs_bit_identical_to_fresh(n, seed):
    rng = np.random.default_rng(1000 * n + seed)
    words = _local_terms(n, 2)
    c0 = rng.normal(size=len(words))
    c1 = rng.normal(size=len(words))
    c1[rng.integers(len(words))] = 0.0  # a zero coefficient is still a term
    L = _build(n, words, c0, rng_seed=seed)
    L.set_h_coeffs(c1)
    L_fresh = _build(n, words, c1, rng_seed=seed)
    _compare_all(L, L_fresh, n, rng)
    # Set again (repeated updates, then back to the original).
    L.set_h_coeffs(c0)
    _compare_all(L, _build(n, words, c0, rng_seed=seed), n, rng)


def test_getter_round_trip_and_copy():
    n = 3
    words = _local_terms(n, 2)
    rng = np.random.default_rng(7)
    c0 = rng.normal(size=len(words))
    L = _build(n, words, c0, rng_seed=0)
    got = L.h_coeffs
    assert got.dtype == np.float64 and got.shape == (len(words),)
    assert np.array_equal(got, c0)
    got[:] = 99.0  # a copy: mutating it must not change L
    assert np.array_equal(L.h_coeffs, c0)
    c1 = rng.normal(size=len(words))
    L.set_h_coeffs(c1)
    assert np.array_equal(L.h_coeffs, c1)
    # list input and non-contiguous views are accepted
    L.set_h_coeffs(list(c0))
    assert np.array_equal(L.h_coeffs, c0)
    big = np.zeros(2 * len(words))
    big[::2] = c1
    L.set_h_coeffs(big[::2])
    assert np.array_equal(L.h_coeffs, c1)
    assert L.num_h_terms == len(words)
    assert L.num_jump_terms > 0


def test_wrong_shapes_raise():
    n = 3
    words = _local_terms(n, 2)
    c0 = np.arange(len(words), dtype=float)
    L = _build(n, words, c0, rng_seed=0)
    for bad in (
        np.zeros(len(words) - 1),
        np.zeros(len(words) + 1),
        np.zeros((len(words), 1)),
        np.zeros((1, len(words))),
        np.zeros(0),
        np.float64(1.0),
    ):
        with pytest.raises(ValueError):
            L.set_h_coeffs(bad)
        with pytest.raises(ValueError):
            L._spec.set_h_coeffs(bad)
    # The raw spec also rejects a non-float64 dtype rather than casting.
    with pytest.raises(ValueError):
        L._spec.set_h_coeffs(np.zeros(len(words), dtype=np.int64))
    assert np.array_equal(L.h_coeffs, c0)  # unchanged after failures


def test_setter_timing_54_terms():
    n = 6
    words = _local_terms(n, 2)[:54]  # 18 one-body + the first 36 nearest-neighbour two-body
    assert len(words) == 54
    rng = np.random.default_rng(3)
    L = _build(n, words, rng.normal(size=54), rng_seed=0)
    cs = rng.normal(size=(10_000, 54))
    t0 = time.perf_counter()
    for c in cs:
        L.set_h_coeffs(c)
    t_wrap = time.perf_counter() - t0
    t0 = time.perf_counter()
    spec = L._spec
    for c in cs:
        spec.set_h_coeffs(c)
    t_raw = time.perf_counter() - t0
    print(
        f"\n10000 set_h_coeffs calls on 54 terms: Lindbladian {t_wrap * 1e6 / 1e4:.2f} us/call, "
        f"_LindbladSpec {t_raw * 1e6 / 1e4:.2f} us/call"
    )
    assert np.array_equal(L.h_coeffs, cs[-1])
    assert t_wrap / 1e4 < 1e-3  # well under a millisecond per call
