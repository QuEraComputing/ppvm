# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Cross-checks for the adjoint Pauli-Lindbladian shim (``ppvm.Lindbladian``).

Each test compares the Rust shim's ``action`` / ``leakage`` / ``generator``
against an independent, dense reference built from the single-qubit Pauli
multiplication table, on a few representative Lindbladians (XY + Z dephasing,
TFIM + X dephasing).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse.linalg as spla

from ppvm import Lindbladian

# --- independent reference via the dense Pauli multiplication table ---------
I, X, Y, Z = range(4)  # noqa: E741  (I is standard Pauli notation here)
MUL = {
    (I, I): (1, I),
    (I, X): (1, X),
    (I, Y): (1, Y),
    (I, Z): (1, Z),
    (X, I): (1, X),
    (X, X): (1, I),
    (X, Y): (1j, Z),
    (X, Z): (-1j, Y),
    (Y, I): (1, Y),
    (Y, X): (-1j, Z),
    (Y, Y): (1, I),
    (Y, Z): (1j, X),
    (Z, I): (1, Z),
    (Z, X): (1j, Y),
    (Z, Y): (-1j, X),
    (Z, Z): (1, I),
}
CODE = {I: "I", X: "X", Y: "Y", Z: "Z"}
ICODE = {"I": I, "X": X, "Y": Y, "Z": Z}


def _str_to_tuple(s):
    return tuple(ICODE[ch] for ch in s)


def _tuple_to_str(t):
    return "".join(CODE[ch] for ch in t)


def _mul_pauli(p, q):
    """Return (phase, p·q) for two Pauli strings given as code tuples."""
    phase = 1
    r = list(p)
    for i, (pi, qi) in enumerate(zip(p, q)):
        ph, rr = MUL[(pi, qi)]
        phase *= ph
        r[i] = rr
    return phase, tuple(r)


def _reference_action(p_str, h_terms, jump_terms):
    """L*(p) = i[H, p] + sum_k gamma_k (L_k p L_k - p), term by term."""
    p = _str_to_tuple(p_str)
    out: dict = {}
    for h_str, coeff_h in h_terms:
        h = _str_to_tuple(h_str)
        ph_pp, pp = _mul_pauli(h, p)
        ph_pq, _ = _mul_pauli(p, h)
        # i[H, p] = i (Hp - pH) = i (ph_pp - ph_pq) r ; real for Hermitian H, p.
        coeff = (1j * coeff_h * (ph_pp - ph_pq)).real
        if coeff:
            out[pp] = out.get(pp, 0.0) + coeff
    for j_str, gamma in jump_terms:
        j = _str_to_tuple(j_str)
        ph_pp, _ = _mul_pauli(j, p)
        # Hermitian Pauli L, p: Lp has imaginary phase iff they anti-commute,
        # and then L p L = -p, contributing -2 gamma p.
        if abs(ph_pp.imag) > 0.5:
            out[p] = out.get(p, 0.0) + (-2.0 * gamma)
    return {kk: v for kk, v in out.items() if v}


def _to_str_dict(d):
    return {_tuple_to_str(kk): v for kk, v in d.items()}


def _xy_dephasing(L, alpha, gamma):
    pairs = [
        (a, b, 1.0 / min(b - a, L - b + a) ** alpha) for a in range(L) for b in range(a + 1, L)
    ]
    kac = sum(j for _, _, j in pairs) / L
    pairs = [(a, b, j / kac) for a, b, j in pairs]
    h_terms = []
    for a, b, j in pairs:
        for q in "XY":
            term = ["I"] * L
            term[a] = term[b] = q
            h_terms.append(("".join(term), j))
    jump_terms = []
    for i in range(L):
        term = ["I"] * L
        term[i] = "Z"
        jump_terms.append(("".join(term), gamma))
    return h_terms, jump_terms


def _tfim_xdeph(L, J, h, gamma):
    h_terms = []
    for i in range(L - 1):
        term = ["I"] * L
        term[i] = term[i + 1] = "Z"
        h_terms.append(("".join(term), J))
    for i in range(L):
        term = ["I"] * L
        term[i] = "X"
        h_terms.append(("".join(term), h))
    jump_terms = []
    for i in range(L):
        term = ["I"] * L
        term[i] = "X"
        jump_terms.append(("".join(term), gamma))
    return h_terms, jump_terms


def _random_pauli_str(rng, L):
    chars = ["I"] * L
    positions = rng.choice(L, size=int(rng.integers(1, L + 1)), replace=False)
    for q in positions:
        chars[q] = rng.choice(["X", "Y", "Z"])
    return "".join(chars)


def _assert_action_matches(L_op, h_terms, jump_terms, strings):
    for p in strings:
        got = L_op.action(p)
        want = _to_str_dict(_reference_action(p, h_terms, jump_terms))
        for kk in set(got) | set(want):
            assert abs(got.get(kk, 0.0) - want.get(kk, 0.0)) < 1e-12, (
                f"action mismatch at p={p!r} k={kk!r}: "
                f"shim={got.get(kk, 0.0)} ref={want.get(kk, 0.0)}"
            )


def test_action_xy_dephasing():
    L = 8
    h_terms, jump_terms = _xy_dephasing(L, alpha=1.0, gamma=0.3)
    L_op = Lindbladian(L, h_terms, jump_terms)
    rng = np.random.default_rng(42)
    strings = ["I" * L]
    strings += ["I" * i + "Z" + "I" * (L - i - 1) for i in range(L)]
    strings += [_random_pauli_str(rng, L) for _ in range(20)]
    _assert_action_matches(L_op, h_terms, jump_terms, strings)


def test_action_tfim_xdephasing():
    L = 6
    h_terms, jump_terms = _tfim_xdeph(L, J=0.7, h=0.4, gamma=0.2)
    L_op = Lindbladian(L, h_terms, jump_terms)
    rng = np.random.default_rng(7)
    strings = [_random_pauli_str(rng, L) for _ in range(30)]
    _assert_action_matches(L_op, h_terms, jump_terms, strings)


def test_generator_leakage_and_expm():
    L = 5
    dt = 0.1
    h_terms, jump_terms = _xy_dephasing(L, alpha=1.0, gamma=0.4)
    L_op = Lindbladian(L, h_terms, jump_terms)
    basis = ["I" * i + "Z" + "I" * (L - i - 1) for i in range(L)]
    basis += ["YIYII", "IYIYI", "ZZIII"]
    coeffs = np.array([0.5, -0.3, 0.2, 0.4, -0.1, 0.6, 0.2, 0.1])

    M_shim = L_op.generator(basis).toarray()
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

    c_shim = spla.expm_multiply(dt * L_op.generator(basis), coeffs)
    c_ref = spla.expm_multiply(dt * M_ref, coeffs)
    assert np.allclose(c_shim, c_ref, atol=1e-13)


def test_protected_strings_suppressed():
    L = 4
    h_terms, jump_terms = _xy_dephasing(L, alpha=1.0, gamma=0.0)
    L_op = Lindbladian(L, h_terms, jump_terms)
    basis = ["ZIII"]
    coeffs = np.array([1.0])
    leak = L_op.leakage(basis, coeffs)
    assert leak, "expected some leakage"
    protected_key = next(iter(leak))
    leak2 = L_op.leakage(basis, coeffs, protected=[protected_key])
    assert protected_key not in leak2
