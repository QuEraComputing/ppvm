# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Trotterized Time Evolution of the XZZ Ising Chain
#
# In this example we simulate the time evolution of the **XZZ Ising chain Hamiltonian**
#
# $$H = h \sum_i X_i + J \sum_i Z_i Z_{i+1}$$
#
# using first-order Trotterized time evolution and **Heisenberg-picture Pauli propagation**.
#
# In the Heisenberg picture we propagate *observables* backwards through the circuit rather
# than evolving a state vector forwards. ppvm represents the observable as a weighted sum of
# Pauli strings and applies each gate layer analytically.
#
# One Trotter step of duration $\delta t$ decomposes the unitary as
#
# $$U(\delta t) \approx \prod_i e^{-i h \delta t X_i} \prod_i e^{-i J \delta t Z_i Z_{i+1}}$$
#
# We also include a local depolarising noise channel after the single-qubit layer and a
# two-qubit depolarising channel after each two-qubit gate to model realistic hardware.

# %%
import matplotlib.pyplot as plt

from ppvm import PauliSum

# %% [markdown]
# ## Parameters
#
# | Symbol | Variable | Description |
# |--------|----------|-------------|
# | $n$ | `n` | Number of qubits |
# | $h$ | `h` | Transverse-field strength |
# | $J$ | `j` | Ising interaction strength ($J = h/8$) |
# | $\delta t$ | `dt` | Trotter step size |
# | $T$ | `time` | Total simulation time |
# | $\varepsilon$ | `min_abs_coeff` | Truncation threshold: Pauli strings with coefficient below this value are discarded after every gate |

# %%
n = 20
h = 1.0
j = h / 8.0
dt = 0.1 / h
time = 1.0 / h

# %% [markdown]
# ## Initial observable
#
# We want to compute $\langle \sum_i Z_i \rangle$ on the all-zeros state.
# In the Heisenberg picture, we initialise the observable as
#
# $$O = \sum_{i=0}^{n-1} Z_i$$
#
# using ppvm's compact notation where `"Z3"` means $Z$ on qubit 3 and $I$ everywhere else.

# %%
state = PauliSum.new(
    n_qubits=n,
    terms=[f"Z{i}" for i in range(n)],
    min_abs_coeff=1e-6,
)
print(state)

# %% [markdown]
# ## Trotter step
#
# Each call to `trotter_step` applies one Trotter step in reverse (Heisenberg picture):
#
# 1. **RX** on every qubit — implements $e^{-i h \delta t X_i}$
# 2. **Single-qubit depolarising noise** on every qubit
# 3. **RZZ** on every neighbouring pair — implements $e^{-i J \delta t Z_i Z_{i+1}}$
# 4. **Two-qubit depolarising noise** on every pair
#
# The `min_abs_coeff` threshold is applied automatically after every gate, so no explicit
# truncation call is needed.

# %%
# Single-qubit depolarising noise: equal probability for X, Y, Z errors
noise_1q = [1e-4, 1e-4, 1e-4]

# Two-qubit depolarising noise: symmetric over all 15 non-identity two-qubit Pauli operators
p_2q = 1e-4
noise_2q = [p_2q / 15.0] * 15


def trotter_step(state, n, theta_x, theta_zz):
    for i in range(n):
        state.rx(i, theta_x)
        state.pauli_error(i, noise_1q)

    for i in range(n - 1):
        state.rzz(i, i + 1, theta_zz)
        state.two_qubit_pauli_error(i, i + 1, noise_2q)


# %% [markdown]
# ## Time evolution
#
# We run the Trotter loop and record $\langle \sum_i Z_i \rangle$ (via `overlap_with_zero`)
# at each time step.

# %%
steps = int(time / dt)
theta_x = dt * h
theta_zz = dt * j

times = [i * dt for i in range(steps + 1)]
expectation_values = []

for step in range(steps):
    expectation_values.append(state.overlap_with_zero())
    trotter_step(state, n, theta_x, theta_zz)

expectation_values.append(state.overlap_with_zero())

print(f"Maximum Pauli weight at final time: {state.current_max_weight()}")

# %% [markdown]
# ## Results
#
# The plot below shows how $\langle \sum_i Z_i \rangle / n$ decays from 1 (all qubits in the
# $|0\rangle$ state) as the system evolves under the combined effect of the transverse field,
# the Ising interaction, and depolarising noise.

# %%
fig, ax = plt.subplots()
ax.plot(times, [ev / n for ev in expectation_values])
ax.set_xlabel("Time $t$")
ax.set_ylabel(r"$\langle \sum_i Z_i \rangle / n$")
ax.set_title("XZZ Ising chain — Trotterized evolution")
plt.tight_layout()
plt.show()
