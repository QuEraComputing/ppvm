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
# # Conserving Magnetization in the XY Ladder with Protected Strings
#
# The two-leg XY ladder
#
# $$H = J \sum_{\langle i,j\rangle} \left(X_i X_j + Y_i Y_j\right)$$
#
# (nearest-neighbour bonds along both legs and across the rungs) conserves the
# total magnetization $M = \sum_i Z_i$. In Heisenberg-picture Pauli propagation
# this is a sum rule on the propagated observable $O(t)$: the single-site
# magnetizations $a_i(t) = \langle Z_i, O(t)\rangle$ satisfy
# $\sum_i a_i(t) = \langle M, O(t)\rangle = \text{const}$.
#
# Truncation can break this. By orthonormality of the Pauli basis, the only words
# overlapping $M$ are the single-site $Z_i$ themselves, so the conserved charge
# leaks precisely when a $Z_i$ at the small-coefficient front is dropped below the
# threshold. ppvm's `preserve_strings` exempts chosen Pauli words from truncation;
# protecting the $\{Z_i\}$ pins the magnetization exactly, at the cost of $n$ words.

# %%
import matplotlib.pyplot as plt

from ppvm import PauliSum

# two-leg ladder: L rungs, N = 2L qubits; qubit (rung j, leg a) -> j + a*L
L, dt, steps, delta = 10, 0.1, 20, 1e-3
N = 2 * L
site = lambda j, a: j + a * L
bonds = [(site(j, a), site((j + 1) % L, a)) for a in (0, 1) for j in range(L)]  # legs (periodic)
bonds += [(site(j, 0), site(j, 1)) for j in range(L)]                          # rungs

M = PauliSum.new(n_qubits=N, terms=[f"Z{q}" for q in range(N)])  # total magnetization
z_strings = ["I" * q + "Z" + "I" * (N - 1 - q) for q in range(N)]  # the single-site Z_i

# %% [markdown]
# One first-order Trotter step. Each bond gate $e^{-i\theta(X_iX_j+Y_iY_j)}$ commutes
# with $M$, but its `rxx` factor alone does not — so we truncate **only after** the
# full gate is applied: `rxx` runs with `truncate=False`, and the following `ryy`
# performs the truncation. We seed the observable with a unit of magnetization on
# the central rung and track $\sum_i a_i$ with and without protecting the $\{Z_i\}$.

# %%
def run(preserve):
    o = PauliSum.new(
        n_qubits=N,
        terms=[(f"Z{site(L // 2, 0)}", 0.5), (f"Z{site(L // 2, 1)}", 0.5)],  # seed: centre rung
        min_abs_coeff=delta,
        preserve_strings=z_strings if preserve else None,
    )
    total_z = [o.overlap(M)]
    for _ in range(steps):
        for a, b in bonds:
            o.rxx(a, b, dt, truncate=False)  # no truncation between the two factors
            o.ryy(a, b, dt)                  # ... truncate only once the bond gate is complete
        total_z.append(o.overlap(M))
    return total_z


times = [i * dt for i in range(steps + 1)]
mz_unprotected = run(preserve=False)
mz_protected = run(preserve=True)
print(f"total magnetization at t = {steps * dt:.1f}:")
print(f"  preserve = 0  ->  {mz_unprotected[-1]:.6f}")
print(f"  preserve = 1  ->  {mz_protected[-1]:.6f}")

# %% [markdown]
# Without protection the magnetization leaks as the spreading front is truncated;
# protecting the $\{Z_i\}$ holds it at its initial value to machine precision.

# %%
fig, ax = plt.subplots()
ax.plot(times, mz_unprotected, "--", label="preserve = 0")
ax.plot(times, mz_protected, "-", label="preserve = 1")
ax.axhline(1.0, ls=":", c="gray")
ax.set_xlabel("Time $t$")
ax.set_ylabel(r"$\sum_i \langle Z_i, O(t)\rangle$")
ax.set_title("XY ladder — magnetization conservation under truncation")
ax.legend()
plt.tight_layout()
plt.show()
