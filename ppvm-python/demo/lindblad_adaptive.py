# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Adaptive Pauli-Lindbladian time evolution on a growing Pauli-string basis.

Demonstrates ``ppvm.Lindbladian``: direct Heisenberg-picture evolution of a
transport observable under a generic Pauli-Lindbladian -- here an all-to-all
XY model (Kac-normalised 1/r^alpha) with single-site Z dephasing -- *without*
Trotterisation. The observable is kept as a sum over a finite set of Pauli
strings (the "basis"), and each step:

  1. measures the leakage  (1 - P_B) L*(O)  out of the current basis B;
  2. adds the leakage strings above ``add_tol`` to B;
  3. advances the coefficient vector by  exp(dt * M),  where
     M = P_B L* P_B  is the Lindbladian restricted to B
     (scipy.sparse.linalg.expm_multiply on the sparse generator);
  4. prunes coefficients below ``drop_tol`` (a protected set -- the target
     observable's own support -- is never dropped).

The matrix exponential is exact in dt within the basis, so there is no
Trotter splitting error; the only approximation is the finite basis, whose
adequacy is monitored by the cumulative discarded weight.

``--predictor_corrector 1`` additionally enriches B with the strings the
predictor step flows into and redoes the step, lifting the adaptive-
integration error from O(dt^2) to O(dt^3) at the cost of one extra matrix
exponential per step.

The observable is the single-Z Fourier mode O_k = sum_j cos(k x_j) Z_j; its
decay C_k(t) = <O_k, O_k(t)> / <O_k, O_k> is the (infinite-temperature)
spin transport coefficient at momentum k.

Run, e.g.::

    python lindblad_adaptive.py --L 12 --kmax 2 --steps 20 --dt 0.1
    python lindblad_adaptive.py --L 12 --kmax 2 --predictor_corrector 1
"""

from argparse import ArgumentParser

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ppvm import Lindbladian

CONSTANTS = dict(
    L=8, alpha=3.0, gamma=0.1, dt=0.05, steps=20, kmax=3,
    add_tol=1e-8, drop_tol=1e-10, max_pauli_weight=0,
    max_basis=0, predictor_corrector=0, out="",
)
p = ArgumentParser(description=__doc__.splitlines()[0])
for key, val in CONSTANTS.items():
    p.add_argument(f"--{key}", type=type(val), default=val)
c = p.parse_args()
if c.max_pauli_weight == 0:
    c.max_pauli_weight = c.L

pairs = [(a, b, min(b - a, c.L - b + a)) for a in range(c.L) for b in range(a + 1, c.L)]
kac = sum(1 / r**c.alpha for _, _, r in pairs) / c.L
pairs = [(a, b, 1 / kac / r**c.alpha) for a, b, r in pairs]
times = np.arange(c.steps + 1) * c.dt
n = np.arange(1, c.kmax + 1)
k_modes = 2 * np.pi * n / c.L
x = (np.arange(c.L) - c.L // 2 + c.L // 2) % c.L - c.L // 2


def zterm_codes(j):
    """Z on site j as a length-L array of Pauli codes (0=I, 1=X, 2=Z, 3=Y)."""
    row = np.zeros(c.L, dtype=np.uint8)
    row[j] = 2  # Z
    return row


def build_hamiltonian_and_jumps():
    """All-to-all XY Hamiltonian H = sum_{a<b} j_ab (X_a X_b + Y_a Y_b),
    with single-site Z dephasing jump operators sqrt(gamma) Z_i."""
    h_terms: list[tuple[str, float]] = []
    for a, b, j in pairs:
        for q in ("X", "Y"):
            s = ["I"] * c.L
            s[a] = q
            s[b] = q
            h_terms.append(("".join(s), j))
    jump_terms: list[tuple[str, float]] = []
    for i in range(c.L):
        s = ["I"] * c.L
        s[i] = "Z"
        jump_terms.append(("".join(s), c.gamma))
    return h_terms, jump_terms


h_terms, jump_terms = build_hamiltonian_and_jumps()
L_op = Lindbladian(c.L, h_terms, jump_terms)


def weights_of(basis_arr):
    """Per-row weight = number of non-identity Paulis."""
    return (basis_arr != 0).sum(axis=1).astype(np.int64)


def dedup_against(candidate_arr, existing_keys):
    """Boolean mask of rows of ``candidate_arr`` not present in ``existing_keys``."""
    return np.array([row.tobytes() not in existing_keys for row in candidate_arr])


def compress(basis_arr, coeff, protected_keys, keep_mask):
    basis_arr = basis_arr[keep_mask]
    coeff = coeff[keep_mask]
    keys = [row.tobytes() for row in basis_arr]
    index = {kk: i for i, kk in enumerate(keys)}
    protected_keys = {pk for pk in protected_keys if pk in index}
    return basis_arr, coeff, index, protected_keys


def cap_basis(basis_arr, coeff, index, protected_keys):
    if not c.max_basis or len(basis_arr) <= c.max_basis:
        return basis_arr, coeff, index, protected_keys, 0
    slots = max(c.max_basis - len(protected_keys), 0)
    keys = [row.tobytes() for row in basis_arr]
    is_protected = np.array([kk in protected_keys for kk in keys])
    cand_idx = np.where(~is_protected)[0]
    order = sorted(cand_idx.tolist(), key=lambda i: abs(coeff[i]), reverse=True)
    keep_extra = order[:slots]
    keep = is_protected.copy()
    keep[keep_extra] = True
    dropped = len(basis_arr) - int(keep.sum())
    return (*compress(basis_arr, coeff, protected_keys, keep), dropped)


Ck = np.empty((c.steps + 1, c.kmax))
Ck0 = np.empty(c.kmax)
n_basis = np.empty_like(Ck, dtype=np.int64)
max_weight = np.empty_like(Ck, dtype=np.int64)
# Cumulative discarded l2 weight / ||O_target|| is an a posteriori error
# estimate. The "total" (trapezoidal) variant adds the post-step leakage so it
# also bounds the O(dt^2) finite-dt adaptive-integration error, not just the
# truncation error.
discarded_cum = np.zeros_like(Ck)
discarded_cum_total = np.zeros_like(Ck)

for m, kk in enumerate(k_modes):
    # Initial basis: Z_j with cosine projection coefficients (the target).
    target_rows, target_coeffs = [], []
    for j, xj in enumerate(x):
        coef = float(np.cos(kk * xj))
        if abs(coef) > c.drop_tol:
            target_rows.append(zterm_codes(j))
            target_coeffs.append(coef)
    basis_arr = np.array(target_rows, dtype=np.uint8)
    coeff = np.array(target_coeffs, dtype=float)
    target_arr = basis_arr.copy()
    target_coeff_arr = coeff.copy()
    protected_arr = basis_arr.copy()  # never dropped, never emitted as leakage
    protected_keys = {row.tobytes() for row in protected_arr}
    index = {row.tobytes(): i for i, row in enumerate(basis_arr)}
    Ck0[m] = float(np.dot(coeff, coeff))
    norm_target = np.sqrt(Ck0[m])

    L_op.clear_cache()  # action cache is per-k since each k starts a fresh basis

    for nt in range(c.steps + 1):
        ck_t = 0.0
        for trow, tc in zip(target_arr, target_coeff_arr):
            idx = index.get(trow.tobytes())
            if idx is not None:
                ck_t += coeff[idx] * tc
        Ck[nt, m] = ck_t / Ck0[m]
        n_basis[nt, m] = len(basis_arr)
        max_weight[nt, m] = int(weights_of(basis_arr).max())
        if nt == c.steps:
            break

        # Leakage out of the basis; the part below add_tol is irrecoverably
        # discarded this step (rate_before), the rest is added to the basis.
        leak_basis, leak_coeffs = L_op.leakage_arr(basis_arr, coeff, protected_arr)
        if len(leak_coeffs):
            lw = weights_of(leak_basis)
            add_mask = (np.abs(leak_coeffs) > c.add_tol) & (lw <= c.max_pauli_weight)
            rate_before = float(np.linalg.norm(leak_coeffs[~add_mask]))
            if add_mask.any():
                cand = leak_basis[add_mask]
                cand_abs = np.abs(leak_coeffs[add_mask])
                cand = cand[np.argsort(cand_abs)[::-1]]
                cand = cand[dedup_against(cand, index)]
                if c.max_basis:
                    cand = cand[: max(c.max_basis - len(basis_arr), 0)]
                if len(cand):
                    n0 = len(basis_arr)
                    basis_arr = np.vstack([basis_arr, cand])
                    coeff = np.r_[coeff, np.zeros(len(cand))]
                    for i, row in enumerate(cand):
                        index[row.tobytes()] = n0 + i
        else:
            rate_before = 0.0

        # Predictor: advance the restricted generator (exact in dt).
        coeff_t = coeff.copy()
        coeff = expm_multiply(c.dt * L_op.generator_arr(basis_arr), coeff)

        if c.predictor_corrector:
            # Enrich with the predictor's second-hop strings and redo the step
            # from the time-t state on the enriched basis (dt^2 -> dt^3).
            lb2, lc2 = L_op.leakage_arr(basis_arr, coeff, protected_arr)
            if len(lc2):
                w2 = weights_of(lb2)
                am2 = (np.abs(lc2) > c.add_tol) & (w2 <= c.max_pauli_weight)
                if am2.any():
                    cand2 = lb2[am2]
                    cand2 = cand2[dedup_against(cand2, index)]
                    if c.max_basis:
                        cand2 = cand2[: max(c.max_basis - len(basis_arr), 0)]
                    if len(cand2):
                        n0 = len(basis_arr)
                        basis_arr = np.vstack([basis_arr, cand2])
                        coeff_t = np.r_[coeff_t, np.zeros(len(cand2))]
                        for i, row in enumerate(cand2):
                            index[row.tobytes()] = n0 + i
                        coeff = expm_multiply(
                            c.dt * L_op.generator_arr(basis_arr), coeff_t
                        )

        w_pre = float(np.dot(coeff, coeff))
        _, leak_after = L_op.leakage_arr(basis_arr, coeff, protected_arr)
        rate_after = float(np.linalg.norm(leak_after)) if len(leak_after) else 0.0

        keep = np.array(
            [(row.tobytes() in protected_keys) or (abs(v) >= c.drop_tol)
             for row, v in zip(basis_arr, coeff)]
        )
        basis_arr, coeff, index, protected_keys = compress(
            basis_arr, coeff, protected_keys, keep
        )
        basis_arr, coeff, index, protected_keys, _ = cap_basis(
            basis_arr, coeff, index, protected_keys
        )
        # Pruning/capping only remove entries, so the discarded l2 weight is
        # sqrt(||c||^2_pre - ||c||^2_post).
        dropped_w = np.sqrt(max(w_pre - float(np.dot(coeff, coeff)), 0.0))
        d_lead = c.dt * rate_before + dropped_w
        d_total = 0.5 * c.dt * (rate_before + rate_after) + dropped_w
        discarded_cum[nt + 1, m] = discarded_cum[nt, m] + d_lead / norm_target
        discarded_cum_total[nt + 1, m] = (
            discarded_cum_total[nt, m] + d_total / norm_target
        )

print(f"L={c.L}  gamma={c.gamma}  alpha={c.alpha}  dt={c.dt}  steps={c.steps}  "
      f"predictor_corrector={c.predictor_corrector}")
for m, kk in enumerate(k_modes):
    print(f"\nk = 2*pi*{n[m]}/{c.L} = {kk:.4f}")
    print(f"  {'t':>6} {'C_k(t)':>10} {'n_basis':>8} {'max_wt':>7} {'err_est':>10}")
    for nt in range(c.steps + 1):
        print(f"  {times[nt]:>6.2f} {Ck[nt, m]:>10.6f} {n_basis[nt, m]:>8d} "
              f"{max_weight[nt, m]:>7d} {discarded_cum_total[nt, m]:>10.2e}")

if c.out:
    np.savez(
        c.out, times=times, k=k_modes, Ck=Ck, Ck0=Ck0, n_basis=n_basis,
        max_weight=max_weight, discarded_cum=discarded_cum,
        discarded_cum_total=discarded_cum_total,
    )
    print(f"\nsaved -> {c.out}")
