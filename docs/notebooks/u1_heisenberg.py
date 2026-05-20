# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
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

# %% [markdown]
# # U(1)-Conserving Trotter Dynamics: XY / Heisenberg Chain
#
# This example shows how to use ppvm's U(1)-symmetric gate helpers to simulate
# z-magnetization-conserving spin dynamics. The Hamiltonian we evolve is the
# XXZ Heisenberg chain with an on-site Z field,
#
# $$ H = \sum_{(i,j)} J_{ij} \big( X_i X_j + Y_i Y_j \big) + \sum_{(i,j)} \Delta_{ij} Z_i Z_j + \sum_i h_i Z_i . $$
#
# Each term commutes with the total magnetization $\sum_k Z_k$, so the
# dynamics preserve the U(1) symmetry sector of any starting observable that
# already lives in that sector. ppvm exposes three pieces of API for this:
#
# 1. `PauliSum.exchange(a, b, theta)` — fused $\mathrm{exp}\!\big({-i\,\theta/2\,(X_a X_b + Y_a Y_b)}\big)$.
# 2. `PauliSum.xyzz(a, b, theta_xy, theta_zz)` — combined XY + ZZ interaction.
# 3. `PauliSum.apply_u1_trotter_step(edges, theta_xy, theta_zz, fields_z)` —
#    a single Trotter slice expressed in terms of the above.
#
# All three act in the **Heisenberg picture**, like every other PauliSum gate.

# %%
from ppvm import PauliSum

# %% [markdown]
# ## Parameters

# %%
n = 4
edges = [(i, i + 1) for i in range(n - 1)]
J_xy = 0.5
Delta = 0.2
field = 0.1
dt = 0.1
n_steps = 5

# %% [markdown]
# ## Observable: total magnetization
#
# We initialise the observable $O = \sum_i Z_i$. Because the XXZ Hamiltonian
# commutes with $\sum_i Z_i$, the propagated observable must stay numerically
# identical to its starting form (up to truncation noise) — a sharp built-in
# check of U(1) symmetry.

# %%
state = PauliSum.new(
    n_qubits=n,
    terms=[f"Z{i}" for i in range(n)],
    min_abs_coeff=1e-12,
)
initial_terms = dict(state.terms)

# %% [markdown]
# ## Trotterized evolution
#
# Each `apply_u1_trotter_step` call applies, for every edge $(i,j)$:
# `xyzz(i, j, J_xy * dt, Delta * dt)`, then a uniform $R_z(field \cdot dt)$ on
# every site. The full circuit is one Heisenberg-picture sweep per call.

# %%
for _ in range(n_steps):
    state.apply_u1_trotter_step(
        edges=edges,
        theta_xy=J_xy * dt,
        theta_zz=Delta * dt,
        fields_z=[field * dt] * n,
    )

# %% [markdown]
# ## Symmetry check
#
# Total magnetization is conserved, so the final and initial PauliSums are
# the same map of Pauli strings to coefficients.

# %%
final_terms = dict(state.terms)
for term, coeff in initial_terms.items():
    assert abs(final_terms.get(term, 0.0) - coeff) < 1e-10, term
print("Total Z preserved across", n_steps, "Trotter steps.")

# %% [markdown]
# ## A starting observable outside the conserved sector
#
# If we instead start from $Z_0 - Z_1$ — the magnetization *difference* — the
# dynamics rotate it inside the two-dimensional U(1) sector spanned by
# $Z_0 - Z_1$ and the current-like operator $Y_0 X_1 - X_0 Y_1$. This is the
# operator-level manifestation of "spin currents flowing across an XY bond".

# %%
diff = PauliSum.new(2, [("Z0", 0.5), ("Z1", -0.5)], min_abs_coeff=1e-12)
diff.exchange(0, 1, theta=0.4)
for name, coeff in diff.terms:
    print(f"  {name}: {coeff:+.4f}")
