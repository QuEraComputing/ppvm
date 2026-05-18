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
# We run two simulations: a noiseless baseline and a noisy version with qubit loss,
# using ppvm's `LossyPauliSum`. Note that the noiseless simulation is still subject
# to truncation and therefore an approximate solution.

# %%
import matplotlib.pyplot as plt

from ppvm import LossyPauliSum, PauliSum

# %% [markdown]
# ## Parameters

# %%
n = 6
h = 1.0
j = 1.5 * h
dt = 0.1 / h
time = 3.0 / h

# %% [markdown]
# ## Noiseless simulation
#
# We want to compute $\langle \sum_i Z_i \rangle$ on the all-zeros state.
# In the Heisenberg picture, we initialise the observable as
#
# $$O = \sum_{i=0}^{n-1} Z_i$$
#
# using ppvm's compact notation where `"Z3"` means $Z$ on qubit 3 and $I$ everywhere else.
#
# Each call to `trotter_step` applies one Trotter step in reverse (Heisenberg picture):
#
# 1. **RX** on every qubit — implements $e^{-i h \delta t X_i}$
# 2. **RZZ** on every neighbouring pair — implements $e^{-i J \delta t Z_i Z_{i+1}}$
#
# Truncation applies automatically after every gate and channel.
# There are currently two possible ways to truncate in a loss-less simulation:
# * Coefficient truncation: every Pauli string with a leading coefficient, whose absolute value is smaller than `min_abs_coeff` is truncated.
# * Max Pauli Weight truncation: every Pauli string with more than `max_pauli_weight` non-identity Paulis gets truncated.

# %%
state = PauliSum.new(
    n_qubits=n,
    terms=[f"Z{i}" for i in range(n)],
    min_abs_coeff=1e-6,
    max_pauli_weight=8,
)

theta_x = dt * h
theta_zz = dt * j


def trotter_step(state, n, theta_x, theta_zz):
    for i in range(n):
        state.rx(i, theta_x)
    for i in range(n - 1):
        state.rzz(i, i + 1, theta_zz)


# %% [markdown]
# ## Noiseless time evolution

# %%
steps = int(time / dt)
times = [i * dt for i in range(steps + 1)]
ev_noiseless = []

for _ in range(steps):
    ev_noiseless.append(state.overlap_with_zero())
    trotter_step(state, n, theta_x, theta_zz)

ev_noiseless.append(state.overlap_with_zero())
print(f"Max Pauli weight (noiseless): {state.current_max_weight()}")

# %% [markdown]
# ## Noisy simulation with qubit loss
#
# We repeat the simulation using `LossyPauliSum`, which extends the Pauli basis with a
# loss operator $L$ to track qubits that have left the computational subspace.
#
# After each gate layer we apply:
# - a **single-qubit depolarising channel** (`pauli_error`) or
#   **two-qubit depolarising channel** (`two_qubit_pauli_error`) to model gate errors
# - a **loss channel** (`loss_channel`) to model qubit loss at the same locations
#
# In addition to the two truncation strategies, we can now also truncate using
# `max_loss_weight`, which removes any Pauli String which has more `L`s than
# that thresholds. This is justified since these Pauli strings only contribute
# little to the final average value.

# %%
noise_1q = [
    1e-3,
    1e-3,
    1e-3,
]  # symmetric single-qubit depolarising: equal p for X, Y, Z
noise_2q = [
    1e-3 / 15.0
] * 15  # symmetric two-qubit depolarising over all 15 non-identity Paulis
p_loss = 1e-3  # loss probability per gate location

noisy_state = LossyPauliSum.new(
    n_qubits=n,
    terms=[f"Z{i}" for i in range(n)],
    min_abs_coeff=1e-6,
    max_pauli_weight=8,
    max_loss_weight=2,
)

# Reset the loss register on every qubit before propagation: this ensures that qubits
# which become lost during the circuit are counted as |0⟩ when computing overlap_with_zero.
# **NOTE**: without truncation, this scales exponentially
for i in range(n):
    noisy_state.reset_loss_channel(i)


def noisy_trotter_step(state, n, theta_x, theta_zz):
    for i in range(n):
        state.rx(i, theta_x)
        state.pauli_error(i, noise_1q)
        state.loss_channel(i, p_loss)

    for i in range(n - 1):
        state.rzz(i, i + 1, theta_zz)
        state.two_qubit_pauli_error(i, i + 1, noise_2q)
        state.loss_channel(i, p_loss)
        state.loss_channel(i + 1, p_loss)


# %% [markdown]
# ## Noisy time evolution

# %%
ev_noisy = []

for _ in range(steps):
    ev_noisy.append(noisy_state.overlap_with_zero())
    noisy_trotter_step(noisy_state, n, theta_x, theta_zz)

ev_noisy.append(noisy_state.overlap_with_zero())
print(f"Max Pauli weight (noisy): {noisy_state.current_max_weight()}")

# %% [markdown]
# ## Results
#
# Both curves start at 1 (all qubits in $|0\rangle$). In the interaction-dominated regime
# ($J > h$) the magnetisation oscillates before decaying. The noisy simulation decays
# faster due to depolarisation and qubit loss.

# %%
fig, ax = plt.subplots()
ax.plot(times, [ev / n for ev in ev_noiseless], label="noiseless")
ax.plot(times, [ev / n for ev in ev_noisy], label="noisy + loss")
ax.set_xlabel("Time $t$")
ax.set_ylabel(r"$\langle \sum_i Z_i \rangle / n$")
ax.set_title("XZZ Ising chain — Trotterized evolution")
ax.legend()
plt.tight_layout()
plt.show()
