# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Shared helpers for the Lindbladian tests.

Two reference kernels live here:

- :func:`_reference_action` builds `L*(p)` for **Hermitian-Pauli** jumps via
  the single-qubit Pauli multiplication table (cheap; only depends on `p`'s
  weight).
- :func:`_dense_action` builds `L*(p)` for **arbitrary jump operators** by
  constructing the full 2^L × 2^L dense Liouvillian. Only viable for L ≤ 3,
  but it makes no assumption about the shape of the jumps.

The bilinear NN-XY + Z-dephasing reference :func:`_bilinear_nn_xy_z_dephasing_obc`
gives an exact closed-form answer that the predictor-corrector tests converge
toward as dt → 0.
"""

from __future__ import annotations

import numpy as np

# --- Pauli multiplication table for the Hermitian-Pauli reference -----------
I, X, Y, Z = range(4)  # noqa: E741  (I is standard Pauli notation here)
MUL = {
    (I, I): (1, I), (I, X): (1, X), (I, Y): (1, Y), (I, Z): (1, Z),
    (X, I): (1, X), (X, X): (1, I), (X, Y): (1j, Z), (X, Z): (-1j, Y),
    (Y, I): (1, Y), (Y, X): (-1j, Z), (Y, Y): (1, I), (Y, Z): (1j, X),
    (Z, I): (1, Z), (Z, X): (1j, Y), (Z, Y): (-1j, X), (Z, Z): (1, I),
}
CODE = {I: "I", X: "X", Y: "Y", Z: "Z"}
ICODE = {"I": I, "X": X, "Y": Y, "Z": Z}


def _str_to_tuple(s):
    return tuple(ICODE[ch] for ch in s)


def _tuple_to_str(t):
    return "".join(CODE[ch] for ch in t)


def _mul_pauli(p, q):
    """Return ``(phase, p·q)`` for two Pauli strings given as code tuples."""
    phase = 1
    r = list(p)
    for i, (pi, qi) in enumerate(zip(p, q)):
        ph, rr = MUL[(pi, qi)]
        phase *= ph
        r[i] = rr
    return phase, tuple(r)


def _reference_action(p_str, h_terms, jump_terms):
    """``L*(p) = i[H, p] + sum_k gamma_k (L_k p L_k - p)``, term by term."""
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


# --- model builders ---------------------------------------------------------


def xy_dephasing(L, alpha, gamma):
    """Long-range XY model with PBC + Z dephasing."""
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
    jump_terms = [("I" * i + "Z" + "I" * (L - i - 1), gamma) for i in range(L)]
    return h_terms, jump_terms


def tfim_xdeph(L, J, h, gamma):
    """TFIM (ZZ + transverse X) with X dephasing."""
    h_terms = []
    for i in range(L - 1):
        term = ["I"] * L
        term[i] = term[i + 1] = "Z"
        h_terms.append(("".join(term), J))
    for i in range(L):
        term = ["I"] * L
        term[i] = "X"
        h_terms.append(("".join(term), h))
    jump_terms = [("I" * i + "X" + "I" * (L - i - 1), gamma) for i in range(L)]
    return h_terms, jump_terms


def nn_xy_z_dephasing_obc(L, J, gamma):
    """Nearest-neighbour XY (OBC) + per-site Z dephasing."""
    h_terms = []
    for i in range(L - 1):
        a, b = i, i + 1
        xs = ["I"] * L
        xs[a] = xs[b] = "X"
        ys = ["I"] * L
        ys[a] = ys[b] = "Y"
        h_terms += [("".join(xs), J), ("".join(ys), J)]
    jump_terms = [("I" * j + "Z" + "I" * (L - j - 1), gamma) for j in range(L)]
    return h_terms, jump_terms


def random_pauli_str(rng, L):
    chars = ["I"] * L
    positions = rng.choice(L, size=int(rng.integers(1, L + 1)), replace=False)
    for q in positions:
        chars[q] = rng.choice(["X", "Y", "Z"])
    return "".join(chars)


def assert_action_matches(L_op, h_terms, jump_terms, strings):
    """Compare `L_op.action(p)` against :func:`_reference_action` for each `p`."""
    for p in strings:
        got = L_op.action(p)
        want = _to_str_dict(_reference_action(p, h_terms, jump_terms))
        for kk in set(got) | set(want):
            assert abs(got.get(kk, 0.0) - want.get(kk, 0.0)) < 1e-12, (
                f"action mismatch at p={p!r} k={kk!r}: "
                f"shim={got.get(kk, 0.0)} ref={want.get(kk, 0.0)}"
            )


# --- dense Liouvillian reference (for non-Hermitian jumps) ------------------

_DENSE_PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}
SIGMA_MINUS_MAT = np.array([[0, 0], [1, 0]], dtype=complex)
SIGMA_PLUS_MAT = np.array([[0, 1], [0, 0]], dtype=complex)


def pauli_mat(s):
    """Dense matrix for Pauli string ``s`` (leftmost char = leftmost factor)."""
    M = np.array([[1.0]], dtype=complex)
    for c in s:
        M = np.kron(M, _DENSE_PAULI[c])
    return M


def all_strings(L):
    if L == 0:
        return [""]
    sub = all_strings(L - 1)
    return [c + s for c in "IXYZ" for s in sub]


def dense_action(H, jumps, p_str, L):
    """`L*(p)` computed densely, returned as a real Pauli-coefficient dict."""
    p_mat = pauli_mat(p_str)
    # macOS Accelerate emits spurious divide warnings on exact zeros.
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        out_mat = 1j * (H @ p_mat - p_mat @ H)
        for Lop, gamma in jumps:
            Ld = Lop.conj().T
            out_mat += gamma * (
                Ld @ p_mat @ Lop - 0.5 * (Ld @ Lop @ p_mat + p_mat @ Ld @ Lop)
            )
    d = 2**L
    out = {}
    for q_str in all_strings(L):
        coef = np.trace(pauli_mat(q_str) @ out_mat) / d
        assert abs(coef.imag) < 1e-9, f"non-real coefficient for {q_str}: {coef}"
        if abs(coef.real) > 1e-12:
            out[q_str] = coef.real
    return out


def embed_op(op_1q, site, L):
    """Embed a single-qubit dense operator at ``site`` of an L-qubit register."""
    eye = _DENSE_PAULI["I"]
    M = np.array([[1.0]], dtype=complex)
    for j in range(L):
        M = np.kron(M, op_1q if j == site else eye)
    return M


# --- closed bilinear evolution for NN-XY + Z-dephasing (OBC) ----------------


def bilinear_nn_xy_z_dephasing_obc(L, J, gamma, times, site0):
    """Closed bilinear evolution of `C_j(t)` for the NN XY + Z-dephasing model
    with open boundary conditions.

    `F_{mn}(t) = 2^{-L} Tr[B_{mn}(t) Z_i]` evolves as
    `∂_t F_{mn} = i·2J·(F_{m+1,n}+F_{m-1,n}-F_{m,n+1}-F_{m,n-1}) - 4γ(1-δ_{mn})F_{mn}`
    on the L×L lattice; OBC means edge terms (m=0 or m=L-1, etc.) drop.
    Z_j = I - 2 n_j gives C_j = -2 F_{jj}.
    """
    dim = L * L

    def idx(m, n):
        return m * L + n

    gen = np.zeros((dim, dim), dtype=complex)
    for m in range(L):
        for n in range(L):
            row = idx(m, n)
            if m + 1 < L:
                gen[row, idx(m + 1, n)] += 1j * 2 * J
            if m - 1 >= 0:
                gen[row, idx(m - 1, n)] += 1j * 2 * J
            if n + 1 < L:
                gen[row, idx(m, n + 1)] += -1j * 2 * J
            if n - 1 >= 0:
                gen[row, idx(m, n - 1)] += -1j * 2 * J
            if m != n:
                gen[row, row] += -4 * gamma

    evals, evecs = np.linalg.eig(gen)
    evecs_inv = np.linalg.inv(evecs)
    f0 = np.zeros(dim, dtype=complex)
    f0[idx(site0, site0)] = -0.5
    coeffs = evecs_inv @ f0

    corr = np.empty((len(times), L))
    diag = [idx(j, j) for j in range(L)]
    for nt, t in enumerate(times):
        ft = evecs @ (np.exp(evals * t) * coeffs)
        corr[nt] = np.real(-2 * ft[diag])
    return corr


# --- numpy-only matrix exponential reference --------------------------------


def coo_to_dense(triples, n_basis):
    """Build a dense ``(n_basis, n_basis)`` array from COO triples."""
    rows, cols, vals = triples
    M = np.zeros((n_basis, n_basis), dtype=float)
    M[rows, cols] = vals
    return M


def expm_mv_dense(M, v):
    """``exp(M) @ v`` via numpy eigendecomposition. Independent of the Rust
    Al-Mohy & Higham implementation; small bases only.

    The Lindbladian is generally diagonalizable, so
    ``M = V diag(λ) V^{-1}`` and ``exp(M) v = V (exp(λ) ⊙ (V^{-1} v))``.
    """
    evals, evecs = np.linalg.eig(M)
    rhs = np.linalg.solve(evecs, v.astype(complex))
    return np.real(evecs @ (np.exp(evals) * rhs))


# --- adaptive PC evolution reference (used by test_adaptive_pc.py) ----------


def _generator_dense(L_op, basis):
    """`L_op.generator(basis)` returns COO triples; convert to dense."""
    return coo_to_dense(L_op.generator(basis), len(basis))


def adaptive_z_correlator(L_op, L, site0, dt, n_steps, tau_add):
    """Adaptive Heisenberg-picture evolution of Z_{site0} on a growing basis.

    First-hop only: each step adds the strings from `L_op.leakage(...)`,
    then matrix-exponentiates the (small, dense) restricted generator via
    numpy eigendecomposition. Local truncation is O(dt²) per step.
    """
    z_strings = ["I" * j + "Z" + "I" * (L - j - 1) for j in range(L)]
    basis = [z_strings[site0]]
    coeffs = np.array([1.0])
    protected = [z_strings[site0]]

    corr = np.zeros((n_steps + 1, L))
    corr[0, site0] = 1.0

    for step in range(n_steps):
        leak = L_op.leakage(basis, coeffs, protected=protected)
        new = [k for k, v in leak.items() if abs(v) > tau_add]
        if new:
            basis = basis + new
            coeffs = np.concatenate([coeffs, np.zeros(len(new))])
        M = _generator_dense(L_op, basis)
        coeffs = expm_mv_dense(dt * M, coeffs)
        index = {s: i for i, s in enumerate(basis)}
        for j in range(L):
            if z_strings[j] in index:
                corr[step + 1, j] = coeffs[index[z_strings[j]]]
    return corr


def adaptive_z_correlator_pc(L_op, L, site0, dt, n_steps, tau_add):
    """Same as :func:`adaptive_z_correlator` but with predictor-corrector
    basis expansion: predict, then add the second-hop leakage strings before
    redoing the step. Lifts the per-step truncation from O(dt²) to O(dt³).
    """
    z_strings = ["I" * j + "Z" + "I" * (L - j - 1) for j in range(L)]
    basis = [z_strings[site0]]
    coeffs = np.array([1.0])
    protected = [z_strings[site0]]

    corr = np.zeros((n_steps + 1, L))
    corr[0, site0] = 1.0

    for step in range(n_steps):
        leak = L_op.leakage(basis, coeffs, protected=protected)
        new = [k for k, v in leak.items() if abs(v) > tau_add]
        if new:
            basis = basis + new
            coeffs = np.concatenate([coeffs, np.zeros(len(new))])
        coeffs_pre = coeffs.copy()
        M = _generator_dense(L_op, basis)
        coeffs_predict = expm_mv_dense(dt * M, coeffs)
        leak2 = L_op.leakage(basis, coeffs_predict, protected=protected)
        new2 = [k for k, v in leak2.items() if abs(v) > tau_add]
        if new2:
            basis = basis + new2
            coeffs_pre = np.concatenate([coeffs_pre, np.zeros(len(new2))])
        M = _generator_dense(L_op, basis)
        coeffs = expm_mv_dense(dt * M, coeffs_pre)
        index = {s: i for i, s in enumerate(basis)}
        for j in range(L):
            if z_strings[j] in index:
                corr[step + 1, j] = coeffs[index[z_strings[j]]]
    return corr
