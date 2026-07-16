# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Shared references for the collective-decay (Kossakowski) tests.

Vendored from the superradiant-burst study (free-space photon-mediated
Lindbladian, arXiv:2309.11376 Eqs. 6-8): emitter couplings from the dyadic
Green's function, the ppvm model builder, and the exact Lindblad reference
via the excitation-number cascade. `exact_rate_sector` is geometry-agnostic
(it takes the J and Gamma matrices), so the same reference covers chains
and rings.

Conventions: excited = |0> (Z = +1), sigma^- = (X - i Y)/2;
H = sum_{n<m} J_nm (s+_n s-_m + h.c.);
D[rho] = sum_{nm} Gamma_nm (s-_m rho s+_n - 1/2 {s+_n s-_m, rho});
observable R(t) = sum_nm Gamma_nm <s+_n s-_m>(t), the photon emission rate.
"""

import itertools

import numpy as np

LAM = 1.0
K0 = 2 * np.pi / LAM
G0 = 1.0
POL = np.array([1.0, 1.0j, 0.0]) / np.sqrt(2.0)


def greens(r_vec):
    """Free-space dyadic Green's function G(r, omega0)."""
    r = np.linalg.norm(r_vec)
    rh = np.outer(r_vec, r_vec) / r**2
    kr = K0 * r
    pref = np.exp(1j * kr) / (4 * np.pi * K0**2 * r**3)
    return pref * ((kr**2 + 1j * kr - 1) * np.eye(3) - (kr**2 + 3j * kr - 3) * rh)


def couplings(pos):
    """J_nm, Gamma_nm from emitter positions; J_nn = 0, Gamma_nn = Gamma_0."""
    n = len(pos)
    J = np.zeros((n, n))
    Gam = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            if a == b:
                Gam[a, b] = G0
                continue
            g = POL.conj() @ greens(pos[a] - pos[b]) @ POL
            J[a, b] = -3 * np.pi * G0 / K0 * g.real
            Gam[a, b] = 6 * np.pi * G0 / K0 * g.imag
    return J, Gam


def chain_positions(n, d):
    return [np.array([j * d, 0.0, 0.0]) for j in range(n)]


def ring_positions(n, d):
    """n emitters on a circle with nearest-neighbour arc spacing ~d.

    Chord-based radius so that adjacent emitters sit exactly d apart.
    The resulting J/Gamma matrices are circulant (exact C_n symmetry).
    """
    radius = d / (2 * np.sin(np.pi / n))
    return [
        np.array([radius * np.cos(2 * np.pi * j / n), radius * np.sin(2 * np.pi * j / n), 0.0])
        for j in range(n)
    ]


def pstr(n, **sites):
    s = ["I"] * n
    for k, v in sites.items():
        s[int(k[1:])] = v
    return "".join(s)


def hamiltonian_terms(n, J):
    """H = sum_{a<b} (J_ab/2)(X_a X_b + Y_a Y_b) as ppvm term pairs."""
    h_terms = []
    for a in range(n):
        for b in range(a + 1, n):
            if abs(J[a, b]) > 1e-14:
                h_terms.append((pstr(n, **{f"q{a}": "X", f"q{b}": "X"}), J[a, b] / 2))
                h_terms.append((pstr(n, **{f"q{a}": "Y", f"q{b}": "Y"}), J[a, b] / 2))
    return h_terms


def eigenmode_jumps(ops, K):
    """Jump list equivalent to the Kossakowski pair (ops, K).

    K = V diag(g) V^dagger; L_nu = sqrt(g_nu) sum_j conj(V_j_nu) A_j with
    rate g_nu (the sqrt is absorbed into the rate as g_nu).
    """
    g_nu, V = np.linalg.eigh(np.asarray(K, dtype=complex))
    jumps = []
    for nu in range(len(ops)):
        if g_nu[nu] < 1e-12:
            continue
        lin = []
        for j, op in enumerate(ops):
            v = np.conj(V[j, nu])
            if abs(v) > 1e-14:
                for p, c in op:
                    lin.append((p, complex(c) * v))
        jumps.append((lin, float(g_nu[nu])))
    return jumps


def rate_observable(n, Gam):
    """O = sum_nm Gamma_nm s+_n s-_m as {pauli_string: real_coeff}."""
    obs = {pstr(n): n * G0 / 2}
    for a in range(n):
        obs[pstr(n, **{f"q{a}": "Z"})] = G0 / 2
        for b in range(a + 1, n):
            obs[pstr(n, **{f"q{a}": "X", f"q{b}": "X"})] = Gam[a, b] / 2
            obs[pstr(n, **{f"q{a}": "Y", f"q{b}": "Y"})] = Gam[a, b] / 2
    return obs


def exact_rate_sector(n, J, Gam, T_run, dt_out, dt_inner=2e-3):
    """Exact Lindblad R(t) via the excitation-number cascade (dense blocks).

    From the fully inverted state, H conserves the excitation number M and
    every jump lowers M symmetrically on both sides of rho, so rho(t) is a
    direct sum of C(n, M)-sized blocks, evolved here with RK4.
    """
    sectors = []
    for M in range(n, -1, -1):
        confs = [frozenset(c) for c in itertools.combinations(range(n), M)]
        sectors.append({c: i for i, c in enumerate(confs)})

    def block(i, mat):
        idx = sectors[i]
        out = np.zeros((len(idx), len(idx)), dtype=complex)
        for conf, b in idx.items():
            for m in conf:
                for nn in range(n):
                    if nn == m:
                        out[b, b] += mat[m, m]
                    elif nn not in conf:
                        out[idx[(conf - {m}) | {nn}], b] += mat[nn, m]
        return out

    H = [block(i, J) for i in range(n + 1)]
    A = [block(i, Gam) for i in range(n + 1)]
    add = []
    for i in range(1, n + 1):
        idx, idx_up = sectors[i], sectors[i - 1]
        amap = np.full((n, len(idx)), -1, dtype=np.int64)
        for conf, b in idx.items():
            for m in range(n):
                if m not in conf:
                    amap[m, b] = idx_up[conf | {m}]
        add.append(amap)

    def feed(i, rho_up):
        amap = add[i - 1]
        out = np.zeros((len(sectors[i]), len(sectors[i])), dtype=complex)
        padded = np.pad(rho_up, ((0, 1), (0, 1)))
        for m in range(n):
            for nn in range(n):
                if abs(Gam[nn, m]) < 1e-14:
                    continue
                sub = padded[amap[m][:, None], amap[nn][None, :]]
                mask = (amap[m][:, None] >= 0) & (amap[nn][None, :] >= 0)
                out += Gam[nn, m] * np.where(mask, sub, 0.0)
        return out

    def rhs(blocks):
        out = []
        for i, r in enumerate(blocks):
            d = -1j * (H[i] @ r - r @ H[i]) - 0.5 * (A[i] @ r + r @ A[i])
            if i > 0:
                d += feed(i, blocks[i - 1])
            out.append(d)
        return out

    blocks = [np.zeros((len(s), len(s)), dtype=complex) for s in sectors]
    blocks[0][0, 0] = 1.0
    n_out = round(T_run / dt_out)
    sub = max(1, int(np.ceil(dt_out / dt_inner)))
    h = dt_out / sub
    R = np.zeros(n_out + 1)
    R[0] = sum(np.trace(r @ a).real for r, a in zip(blocks, A))
    for k in range(n_out):
        for _ in range(sub):
            k1 = rhs(blocks)
            k2 = rhs([r + h / 2 * d for r, d in zip(blocks, k1)])
            k3 = rhs([r + h / 2 * d for r, d in zip(blocks, k2)])
            k4 = rhs([r + h * d for r, d in zip(blocks, k3)])
            blocks = [
                r + h / 6 * (a + 2 * b + 2 * c + e) for r, a, b, c, e in zip(blocks, k1, k2, k3, k4)
            ]
        R[k + 1] = sum(np.trace(r @ a).real for r, a in zip(blocks, A))
    return np.arange(n_out + 1) * dt_out, R
