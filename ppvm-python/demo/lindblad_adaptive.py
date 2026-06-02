# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: ppvm (3.12.12)
#     language: python
#     name: python3
# ---

# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

# %% [markdown]
# # Adaptive Pauli-Lindbladian time evolution
#
# Direct Heisenberg-picture evolution of a transport observable on a growing
# Pauli-string basis, without Trotterisation. Composes the three
# `ppvm.Lindbladian` primitives — `leakage`, `generator`, and the matrix
# exponential — into a predictor-corrector integrator on an all-to-all XY
# model (Kac-normalised $1/r^\alpha$ couplings) with single-site Z dephasing.
#
# The observable is kept as a finite sum over Pauli strings (the *basis*),
# and each step:
#
# 1. measures the leakage $(\mathbf 1 - P_B)\,\mathcal L^\dagger(\mathcal O)$
#    out of the current basis $B$;
# 2. adds the leakage strings above `add_tol` to $B$;
# 3. advances the coefficient vector by $\exp(dt\,M)$, where
#    $M = P_B \mathcal L^\dagger P_B$ is the Lindbladian restricted to $B$;
# 4. prunes coefficients below `drop_tol`; the protected set (the target
#    observable's own support) is never dropped.
#
# The matrix exponential is exact in `dt` within the basis, so there is no
# Trotter splitting error; the only approximation is the finite basis,
# monitored by the cumulative discarded weight.
#
# With `predictor_corrector=True`, the predicted state is fed back as a
# second leakage probe — enriching $B$ with the strings the predictor flows
# into before re-running the step from the pre-step state. This lifts the
# per-step adaptive-integration error from $\mathcal O(dt^2)$ to $\mathcal O(dt^3)$.
#
# The observable is the single-Z Fourier mode
# $\mathcal O_k = \sum_j \cos(k x_j) Z_j$; its decay
# $C_k(t)/C_k(0)$ is the (infinite-temperature) spin transport coefficient.

# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

from ppvm import Lindbladian

# %% [markdown]
# ## Parameters

# %%
L = 8
alpha = 3.0
gamma = 0.1
dt = 0.05
steps = 20
kmax = 3
add_tol = 1e-8
drop_tol = 1e-10
max_pauli_weight = L
max_basis = 0  # 0 = no cap
predictor_corrector = True

times = np.arange(steps + 1) * dt
k_indices = np.arange(1, kmax + 1)
k_modes = 2 * np.pi * k_indices / L
x = (np.arange(L) - L // 2 + L // 2) % L - L // 2


# %% [markdown]
# ## Model
#
# All-to-all XY Hamiltonian $H = \sum_{a<b} J_{ab} (X_aX_b + Y_aY_b)$ with
# $J_{ab} = (1/r_{ab}^\alpha)/k_\mathrm{ac}$ (Kac-normalised) and per-site Z
# dephasing jumps $\sqrt{\gamma} Z_i$.

# %%
def build_xy_dephasing(L, alpha, gamma):
    pairs = [(a, b, min(b - a, L - b + a)) for a in range(L) for b in range(a + 1, L)]
    kac = sum(1 / r**alpha for _, _, r in pairs) / L
    pairs = [(a, b, 1 / kac / r**alpha) for a, b, r in pairs]
    h_terms = []
    for a, b, j in pairs:
        for q in ("X", "Y"):
            s = ["I"] * L
            s[a] = s[b] = q
            h_terms.append(("".join(s), j))
    jump_terms = []
    for i in range(L):
        s = ["I"] * L
        s[i] = "Z"
        jump_terms.append(("".join(s), gamma))
    return h_terms, jump_terms


h_terms, jump_terms = build_xy_dephasing(L, alpha, gamma)
L_op = Lindbladian(L, h_terms, jump_terms)


# %% [markdown]
# ## Adaptive-evolution helpers
#
# `init_mode(kk)` seeds the basis with the target observable
# $\sum_j \cos(k x_j) Z_j$, dropping coefficients below `drop_tol`.
# `add_leakage_to_basis` extends the basis with above-threshold leakage from
# a probe coefficient vector. `prune_and_cap` drops below-threshold
# coefficients and (optionally) caps the basis size.

# %%
def zterm_codes(j):
    """Z on site j as a length-L array of Pauli codes (0=I, 1=X, 2=Z, 3=Y)."""
    row = np.zeros(L, dtype=np.uint8)
    row[j] = 2
    return row


def weights_of(basis_arr):
    return (basis_arr != 0).sum(axis=1).astype(np.int64)


def generator_sparse(L_op, basis_arr):
    """COO triples → scipy CSC. `generator_arr` itself has no scipy dep; this
    convenience wrapper is local to the demo."""
    rows, cols, vals = L_op.generator_arr(basis_arr)
    n = basis_arr.shape[0]
    return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()


def init_mode(kk):
    rows, coeffs = [], []
    for j, xj in enumerate(x):
        c = float(np.cos(kk * xj))
        if abs(c) > drop_tol:
            rows.append(zterm_codes(j))
            coeffs.append(c)
    basis_arr = np.array(rows, dtype=np.uint8)
    coeff = np.array(coeffs, dtype=float)
    index = {row.tobytes(): i for i, row in enumerate(basis_arr)}
    protected_keys = set(index)
    return basis_arr, coeff, index, protected_keys


def add_leakage_to_basis(basis_arr, probe_coeff, extend_coeffs, protected_arr, index):
    """Compute leakage from `probe_coeff`, add above-threshold strings to
    `basis_arr`, and pad each vector in `extend_coeffs` with zeros for the
    new rows. Returns ``(basis_arr, [extended...], index, rate_below)`` where
    `rate_below` is the l2 norm of the leakage too small to add."""
    leak_basis, leak_coeffs = L_op.leakage_arr(basis_arr, probe_coeff, protected_arr)
    if not len(leak_coeffs):
        return basis_arr, list(extend_coeffs), index, 0.0
    weights = weights_of(leak_basis)
    add_mask = (np.abs(leak_coeffs) > add_tol) & (weights <= max_pauli_weight)
    rate_below = float(np.linalg.norm(leak_coeffs[~add_mask]))
    if not add_mask.any():
        return basis_arr, list(extend_coeffs), index, rate_below
    cand = leak_basis[add_mask]
    cand = cand[np.argsort(np.abs(leak_coeffs[add_mask]))[::-1]]
    cand = cand[np.array([row.tobytes() not in index for row in cand])]
    if max_basis:
        cand = cand[: max(max_basis - len(basis_arr), 0)]
    if not len(cand):
        return basis_arr, list(extend_coeffs), index, rate_below
    n0 = len(basis_arr)
    basis_arr = np.vstack([basis_arr, cand])
    extended = [np.r_[c, np.zeros(len(cand))] for c in extend_coeffs]
    for i, row in enumerate(cand):
        index[row.tobytes()] = n0 + i
    return basis_arr, extended, index, rate_below


def prune_and_cap(basis_arr, coeff, protected_keys):
    """Drop below-`drop_tol` coefficients; if `max_basis` is set, cap by
    keeping protected rows + largest-magnitude others."""
    keep = np.array(
        [(row.tobytes() in protected_keys) or (abs(v) >= drop_tol)
         for row, v in zip(basis_arr, coeff)]
    )
    basis_arr = basis_arr[keep]
    coeff = coeff[keep]
    keys = [row.tobytes() for row in basis_arr]
    index = {kk: i for i, kk in enumerate(keys)}
    protected_keys = {pk for pk in protected_keys if pk in index}
    if max_basis and len(basis_arr) > max_basis:
        slots = max(max_basis - len(protected_keys), 0)
        is_protected = np.array([kk in protected_keys for kk in keys])
        order = sorted(np.where(~is_protected)[0].tolist(),
                       key=lambda i: abs(coeff[i]), reverse=True)
        keep2 = is_protected.copy()
        keep2[order[:slots]] = True
        basis_arr = basis_arr[keep2]
        coeff = coeff[keep2]
        keys = [row.tobytes() for row in basis_arr]
        index = {kk: i for i, kk in enumerate(keys)}
        protected_keys = {pk for pk in protected_keys if pk in index}
    return basis_arr, coeff, index, protected_keys


# %% [markdown]
# ## Run one $k$-mode
#
# Returns per-step $C_k(t)$, basis size, max weight, and cumulative discarded
# weight (trapezoidal a-posteriori error estimate normalised by
# $\|\mathcal O_k\|$).

# %%
def run_mode(kk):
    basis_arr, coeff, index, protected_keys = init_mode(kk)
    # `basis_arr` rebinds to a new ndarray on every vstack/slice; nothing
    # mutates the initial array in place, so the target/protected views can
    # simply alias it without a copy.
    target_arr = basis_arr
    target_coeff = coeff
    protected_arr = basis_arr
    norm0 = float(np.dot(coeff, coeff))
    norm_target = np.sqrt(norm0)

    ck = np.empty(steps + 1)
    n_basis_t = np.empty(steps + 1, dtype=np.int64)
    max_w_t = np.empty(steps + 1, dtype=np.int64)
    discarded_cum = np.zeros(steps + 1)

    L_op.clear_cache()
    for nt in range(steps + 1):
        # Overlap with the (fixed) initial target.
        c_t = sum(coeff[index[trow.tobytes()]] * tc
                  for trow, tc in zip(target_arr, target_coeff)
                  if trow.tobytes() in index)
        ck[nt] = c_t / norm0
        n_basis_t[nt] = len(basis_arr)
        max_w_t[nt] = int(weights_of(basis_arr).max())
        if nt == steps:
            break

        # Predictor: enrich basis with leakage from current state, then
        # advance by exp(dt · M).
        basis_arr, [coeff], index, rate_before = add_leakage_to_basis(
            basis_arr, coeff, [coeff], protected_arr, index
        )
        coeff_pre = coeff.copy()
        coeff = expm_multiply(dt * generator_sparse(L_op, basis_arr), coeff)

        if predictor_corrector:
            # Probe leakage with the predicted state; extend the pre-step
            # vector with zeros for the new rows and re-run on the enlarged
            # basis. Lifts O(dt²) -> O(dt³).
            basis_arr, [coeff_pre], index, _ = add_leakage_to_basis(
                basis_arr, coeff, [coeff_pre], protected_arr, index
            )
            coeff = expm_multiply(dt * generator_sparse(L_op, basis_arr), coeff_pre)

        # Post-step leakage rate, for the trapezoidal error estimate.
        norm2_pre_prune = float(np.dot(coeff, coeff))
        _, leak_after = L_op.leakage_arr(basis_arr, coeff, protected_arr)
        rate_after = float(np.linalg.norm(leak_after)) if len(leak_after) else 0.0

        basis_arr, coeff, index, protected_keys = prune_and_cap(
            basis_arr, coeff, protected_keys
        )
        # Pruning only removes entries, so the discarded l2 weight is
        # sqrt(‖c‖²_pre − ‖c‖²_post).
        dropped_w = np.sqrt(max(norm2_pre_prune - float(np.dot(coeff, coeff)), 0.0))
        d_total = 0.5 * dt * (rate_before + rate_after) + dropped_w
        discarded_cum[nt + 1] = discarded_cum[nt] + d_total / norm_target

    return ck, n_basis_t, max_w_t, discarded_cum


# %% [markdown]
# ## Run all $k$-modes and plot

# %%
Ck = np.empty((steps + 1, kmax))
n_basis = np.empty((steps + 1, kmax), dtype=np.int64)
max_weight = np.empty((steps + 1, kmax), dtype=np.int64)
discarded_cum = np.empty((steps + 1, kmax))
for m, kk in enumerate(k_modes):
    Ck[:, m], n_basis[:, m], max_weight[:, m], discarded_cum[:, m] = run_mode(kk)

# %%
fig, ax = plt.subplots()
for m in range(kmax):
    ax.plot(times, Ck[:, m], "o-", ms=3, label=rf"$k = 2\pi\cdot{k_indices[m]}/L$")
ax.set_xlabel("$t$")
ax.set_ylabel(r"$C_k(t)/C_k(0)$")
ax.set_title(f"Adaptive Lindbladian evolution  L={L}  γ={gamma}  α={alpha}")
ax.legend()
plt.tight_layout()
plt.show()

# %%
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for m in range(kmax):
    ax[0].plot(times, n_basis[:, m], "o-", ms=3, label=rf"$k_{k_indices[m]}$")
    ax[1].semilogy(times, discarded_cum[:, m] + 1e-16, "o-", ms=3,
                   label=rf"$k_{k_indices[m]}$")
ax[0].set(xlabel="$t$", ylabel="|basis|", title="Basis-size growth")
ax[1].set(xlabel="$t$", ylabel=r"cum. discarded / $||O_k||$",
          title="A-posteriori error estimate")
for a in ax:
    a.legend()
plt.tight_layout()
plt.show()
