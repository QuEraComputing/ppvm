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
# # `pc_step` parallel scaling
#
# End-to-end wall-time scaling of the pure-Rust predictor-corrector step
# (`ppvm.Lindbladian.pc_step`) with rayon thread count. The entire `pc_step`
# body — both leakage calls, the generator build, and both matrix
# exponentials — runs inside a rayon pool of the requested size:
#
# * leakage and generator parallelise over basis elements;
# * the matrix exponential parallelises over SpMV rows.
#
# So the speedup numbers reflect overall PC throughput, not just SpMV.

# %%
from statistics import median
import time

import matplotlib.pyplot as plt
import numpy as np

from ppvm import Lindbladian


# %% [markdown]
# ## Parameters

# %%
L = 51
J = 1.0
gamma = 1.0
alpha = 1.0
dt = 0.05
n_steps = 20
max_basis = 10_000_000  # large: rank cap never binds (full enrichment)
max_cores = 4
warmup_steps = 4
model = "long-range"  # "nn" or "long-range"


# %% [markdown]
# ## Model
#
# All-to-all XY with $1/r^\alpha$ couplings (Kac-normalised) and per-site Z
# dephasing. Long-range activates every bond every step, giving a basis
# size that meaningfully exercises parallel scaling.

# %%
def build_nn_xy_dephasing(L, J, gamma):
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


def build_long_range_xy_dephasing(L, J, alpha, gamma):
    pairs = [
        (a, b, 1.0 / min(b - a, L - b + a) ** alpha)
        for a in range(L)
        for b in range(a + 1, L)
    ]
    kac = sum(j for _, _, j in pairs) / L
    pairs = [(a, b, j / kac) for a, b, j in pairs]
    h_terms = []
    for a, b, j in pairs:
        for q in "XY":
            term = ["I"] * L
            term[a] = term[b] = q
            h_terms.append(("".join(term), J * j))
    jump_terms = [("I" * j + "Z" + "I" * (L - j - 1), gamma) for j in range(L)]
    return h_terms, jump_terms


if model == "nn":
    h_terms, jump_terms = build_nn_xy_dephasing(L, J, gamma)
else:
    h_terms, jump_terms = build_long_range_xy_dephasing(L, J, alpha, gamma)
L_op = Lindbladian(L, h_terms, jump_terms)


# %% [markdown]
# ## Timing harness
#
# Each call to `run_pc_steps` runs `n_steps` consecutive PC steps from
# $Z_{L/2}$, returning the per-step wall times and the final basis size.
# The `num_threads` kwarg pins this call to a freshly-built rayon pool of
# that size, isolating thread-count effects from JIT cache state.

# %%
def run_pc_steps(L_op, L, site0, dt, n_steps, max_basis, num_threads):
    z_strings = ["I" * j + "Z" + "I" * (L - j - 1) for j in range(L)]
    basis = [z_strings[site0]]
    coeffs = np.array([1.0])
    protected = [z_strings[site0]]
    times = []
    for _ in range(n_steps):
        t0 = time.perf_counter()
        basis, coeffs = L_op.pc_step(
            basis,
            coeffs,
            dt,
            max_basis,
            protected=protected,
            num_threads=num_threads,
        )
        times.append(time.perf_counter() - t0)
    return times, len(basis)


# %% [markdown]
# ## Warmup
#
# Each thread count pre-builds its rayon pool and amortises one-time setup
# before the real timing pass.

# %%
site0 = L // 2
for n in range(1, max_cores + 1):
    run_pc_steps(L_op, L, site0, dt, warmup_steps, max_basis, n)


# %% [markdown]
# ## Scaling sweep

# %%
results = []
for n in range(1, max_cores + 1):
    times, basis_size = run_pc_steps(L_op, L, site0, dt, n_steps, max_basis, n)
    first = times[0] * 1000.0
    steady = median(times[1:]) * 1000.0
    results.append({"threads": n, "first_ms": first, "steady_ms": steady,
                    "basis": basis_size})

baseline = results[0]["steady_ms"]
print(f"{'threads':>8s}  {'first-step (ms)':>16s}  {'steady (ms)':>12s}  "
      f"{'speedup':>9s}  {'|basis|':>8s}")
for r in results:
    speedup = baseline / r["steady_ms"]
    print(f"{r['threads']:>8d}  {r['first_ms']:>16.1f}  {r['steady_ms']:>12.2f}  "
          f"{speedup:>8.2f}x  {r['basis']:>8d}")


# %% [markdown]
# ## Plot
#
# Steady-state speedup vs thread count, with the linear-scaling reference.

# %%
threads = np.array([r["threads"] for r in results])
steady_ms = np.array([r["steady_ms"] for r in results])
speedup = steady_ms[0] / steady_ms

fig, ax = plt.subplots()
ax.plot(threads, speedup, "o-", label="measured")
ax.plot(threads, threads, "k--", alpha=0.4, label="linear")
ax.set_xlabel("threads")
ax.set_ylabel("speedup (vs 1 thread)")
ax.set_title(f"pc_step parallel scaling  L={L}  model={model}  "
             f"|basis|={results[-1]['basis']}")
ax.legend()
plt.tight_layout()
plt.show()
