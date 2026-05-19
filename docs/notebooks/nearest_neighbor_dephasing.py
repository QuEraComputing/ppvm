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
# # Dephased Nearest-Neighbour XY Chain
#
# In this example we simulate the infinite-temperature spin correlator of the
# **nearest-neighbour XY chain with local dephasing**
#
# $$H = J \sum_i (X_i X_{i+1} + Y_i Y_{i+1})$$
#
# with periodic boundary conditions.  We work in the Heisenberg picture, where
# observables evolve under the adjoint Lindblad equation
#
# $$\dot O=i[H,O]+\gamma\sum_j(Z_jOZ_j-O).$$
#
# The observable is
#
# $$C_j(t)=2^{-L}\operatorname{Tr}[Z_j(t)Z_i],$$
#
# where $i$ is the middle site of the chain.  The goal is to extract transport
# coefficients from the long-wavelength relaxation of this correlator.  In one
# dimension, diffusion implies that the local autocorrelator decays as
#
# $$C_i(t)\sim t^{-1/2}.$$
#
# It is often cleaner to analyze the same dynamics in momentum space, with
#
# $$C_k(t)=\sum_j e^{-ikx_j}C_j(t),$$
#
# where $x_j$ is the periodic distance from the initially perturbed site.  For a
# diffusive mode one expects
#
# $$C_k(t)\sim e^{-Dk^2t}$$
#
# at small momentum, so the decay rate $\lambda(k)$ gives the diffusion
# constant through $\lambda(k)/k^2\to D$.
#
# We compute the correlator in two ways:
#
# 1. Trotterized Heisenberg-picture Pauli propagation with ppvm.
# 2. An exact finite-size benchmark obtained by solving the closed bilinear
#    operator equation for $B_{mn}=c_m^\dagger c_n$.
#
# The ppvm result is approximate because of Trotter error and Pauli-string
# truncation.  The bilinear benchmark is exact for the nearest-neighbour
# dephased free-fermion problem.

# %%
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from ppvm import PauliSum

# %% [markdown]
# ## Parameters

# %%
L = 51
J = 1.0
gamma = 1.0
dt = 0.01
time = 4.0
min_abs_coeff = 1e-8
max_pauli_weight = L+1 # no weight truncation

site0 = L // 2
steps = int(time / dt)
times = np.arange(steps + 1) * dt

# %% [markdown]
# ## Trotterized ppvm simulation
#
# We initialise the Heisenberg observable as $Z_i$.  During each Trotter step
# we apply local $Z$-dephasing and then reverse the nearest-neighbour
# $R_{XX}$, $R_{YY}$ gate layers.  The bonds are split into even and odd
# layers.  Bonds within one layer are disjoint for even $L$, so this ordering is
# closer to the usual parallel nearest-neighbour Trotter circuit than sweeping
# through all bonds one by one.
#
# The dephasing part of the adjoint Lindblad equation
#
# $$\dot O=\gamma\sum_j(Z_jOZ_j-O)$$
#
# damps a single $X$ or $Y$ operator by $e^{-2\gamma\delta t}$.  Therefore
# the corresponding Pauli-error probability per Trotter step is
#
# $$p_Z=\frac{1-e^{-2\gamma\delta t}}{2}.$$

# %%
observable = PauliSum.new(
    L,
    f"Z{site0}",
    min_abs_coeff=min_abs_coeff,
    max_pauli_weight=max_pauli_weight,
)
z_ops = [PauliSum.new(L, f"Z{j}") for j in range(L)]

theta = 2 * J * dt
p_z = (1 - np.exp(-2 * gamma * dt)) / 2
bond_layers = [
    [(j, (j + 1) % L) for j in range(0, L, 2)],
    [(j, (j + 1) % L) for j in range(1, L, 2)],
]


def trotter_step(obs):
    for q in range(L):
        obs.pauli_error(q, [0.0, 0.0, p_z])
    for layer in reversed(bond_layers):
        for a, b in reversed(layer):
            obs.ryy(a, b, theta)
            obs.rxx(a, b, theta)


# %%
ppvm_corr = np.empty((steps + 1, L))

for n, _ in tqdm(enumerate(times), total=len(times)):
    ppvm_corr[n] = [observable.overlap(z) for z in z_ops]
    if n == steps:
        break
    trotter_step(observable)

print(f"Current max Pauli weight: {observable.current_max_weight()}")

# %% [markdown]
# ## Exact bilinear benchmark
#
# For nearest-neighbour interactions the Jordan-Wigner fermions remain closed
# under the adjoint Lindblad equation.  The bilinears
#
# $$B_{mn}=c_m^\dagger c_n$$
#
# obey
#
# $$
# \partial_t B_{mn}
# =
# i\,2J(B_{m+1,n}+B_{m-1,n}-B_{m,n+1}-B_{m,n-1})
# -4\gamma(1-\delta_{mn})B_{mn}.
# $$
#
# By linearity, the quantities
#
# $$F_{mn}(t)=2^{-L}\operatorname{Tr}[B_{mn}(t)Z_i]$$
#
# obey the same equation with initial condition
#
# $$F_{mn}(0)=-\frac{1}{2}\delta_{mi}\delta_{ni}.$$
#
# The spin correlator is recovered from the diagonal entries:
#
# $$C_j(t)=-2F_{jj}(t).$$

# %%
def exact_bilinear_correlator(L, J, gamma, times, site0):
    dim = L * L

    def idx(m, n):
        return (m % L) * L + (n % L)

    generator = np.zeros((dim, dim), dtype=complex)
    for m in range(L):
        for n in range(L):
            row = idx(m, n)
            generator[row, idx(m + 1, n)] += 1j * 2 * J
            generator[row, idx(m - 1, n)] += 1j * 2 * J
            generator[row, idx(m, n + 1)] += -1j * 2 * J
            generator[row, idx(m, n - 1)] += -1j * 2 * J
            if m != n:
                generator[row, row] += -4 * gamma

    evals, evecs = np.linalg.eig(generator)
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


exact_corr = exact_bilinear_correlator(L, J, gamma, times, site0)

# %% [markdown]
# ## Results
#
# The autocorrelator checks the time dependence at the initially perturbed site.
# We plot it on log-log axes and then show the relative error to the exact
# bilinear solution.  The final-time profile checks spatial propagation around
# the periodic chain.

# %%
fig, ax = plt.subplots()
ax.loglog(times[1:], ppvm_corr[1:, site0], "o", ms=3, label="ppvm")
ax.loglog(times[1:], exact_corr[1:, site0], "-", label="exact bilinear")
ax.set_xlabel("Time $t$")
ax.set_ylabel(r"$C_i(t)$")
ax.set_title("Nearest-neighbour XY chain with dephasing")
ax.legend()
plt.tight_layout()
plt.show()

# %%
rel_err = np.abs(ppvm_corr[:, site0] - exact_corr[:, site0]) / np.maximum(
    np.abs(exact_corr[:, site0]), 1e-15
)

fig, ax = plt.subplots()
ax.semilogy(times, rel_err, "o-", ms=3)
ax.set_xlabel("Time $t$")
ax.set_ylabel("Relative error")
ax.set_title("Autocorrelator error relative to exact bilinear solution")
plt.tight_layout()
plt.show()

# %%
x = (np.arange(L) - site0 + L // 2) % L - L // 2
order = np.argsort(x)

fig, ax = plt.subplots()
ax.plot(x[order], ppvm_corr[-1, order], "o", ms=4, label="ppvm")
ax.plot(x[order], exact_corr[-1, order], "-", label="exact bilinear")
ax.set_xlabel("Distance from initial site")
ax.set_ylabel(r"$C_j(t)$")
ax.set_title(f"Profile at $t={times[-1]:.2f}$")
ax.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Fourier-space decay
#
# We now inspect low-momentum relaxation by Fourier transforming the spatial correlator,
#
# $$C_k(t)=\sum_j e^{-ikx_j}C_j(t).$$
#
# For diffusive transport the low-$k$ decay rate scales as
# $\lambda(k)\propto k^2$.

# %%
k = 2 * np.pi * np.arange(1, min(5, L // 2) + 1) / L
phase = np.exp(-1j * np.outer(x, k))
ppvm_ck = ppvm_corr @ phase
exact_ck = exact_corr @ phase

fig, ax = plt.subplots()
for n, kk in enumerate(k):
    ax.plot(times, np.real(ppvm_ck[:, n] / ppvm_ck[0, n]), "o", ms=3, label=rf"ppvm $k={kk:.2f}$")
    ax.plot(times, np.real(exact_ck[:, n] / exact_ck[0, n]), "-", label=rf"exact $k={kk:.2f}$")
ax.set_xlabel("Time $t$")
ax.set_ylabel(r"$C_k(t)/C_k(0)$")
ax.set_yscale("log")
ax.set_title("Low-momentum mode decay")
ax.legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Diffusive rate fit
#
# We fit each Fourier mode to
#
# $$\log |C_k(t)| = a_k - \lambda(k)t$$
#
# in the window $2\leq t\leq4$.  For diffusion, $\lambda(k)/k^2$ should approach
# a constant.  In the infinite-system hydrodynamic limit of the exact
# nearest-neighbour solution, that constant is
#
# $$D_\infty=\frac{2J^2}{\gamma}.$$

# %%
fit_tmin = 2.0
fit_tmax = 4.0
fit_mask = (times >= fit_tmin) & (times <= fit_tmax)

lambda_k = []
for n in range(len(k)):
    y = np.abs(exact_ck[:, n])
    mask = fit_mask & (y > 0)
    slope, intercept = np.polyfit(times[mask], np.log(y[mask]), 1)
    lambda_k.append(-slope)
lambda_k = np.array(lambda_k)

D_exact = 2 * J**2 / gamma
print(f"D_infinity = {D_exact:.6g}")

fig, ax = plt.subplots()
ax.plot(k, lambda_k / k**2, "o-", label=r"$\lambda(k)/k^2$")
ax.axhline(D_exact, ls="--", color="black", label=rf"$D_\infty={D_exact:.3g}$")
ax.set_xlabel(r"$k$")
ax.set_ylabel(r"$\lambda(k)/k^2$")
ax.set_title(r"Diffusion estimate from $2\leq t\leq4$")
ax.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Appendix: Errors in ppvm simulation
#
# Here we discuss the interplay between Trotter and truncation errors. We vary only $\gamma$, keeping all other parameters fixed, and compare
# the ppvm autocorrelator to the exact bilinear result.  The plotted error is
# the maximum relative error over all times.  We repeat the comparison with the
# original Trotter step and with a ten times smaller step size.

# %%
gamma_values = np.array([1e-4,0.001, 0.01,0.1, 1.0])
dt_values = [dt, dt / 10]
max_rel_errors = {dt_value: [] for dt_value in dt_values}


def ppvm_autocorrelator(gamma_value, dt_value):
    steps_value = int(time / dt_value)
    times_value = np.arange(steps_value + 1) * dt_value
    theta_value = 2 * J * dt_value
    p_z_value = (1 - np.exp(-2 * gamma_value * dt_value)) / 2
    obs = PauliSum.new(
        L,
        f"Z{site0}",
        min_abs_coeff=min_abs_coeff,
        max_pauli_weight=max_pauli_weight,
    )
    z0 = PauliSum.new(L, f"Z{site0}")
    ppvm_auto = np.empty(len(times_value))

    for n, _ in enumerate(times_value):
        ppvm_auto[n] = obs.overlap(z0)
        if n == steps_value:
            break
        for q in range(L):
            obs.pauli_error(q, [0.0, 0.0, p_z_value])
        for layer in reversed(bond_layers):
            for a, b in reversed(layer):
                obs.ryy(a, b, theta_value)
                obs.rxx(a, b, theta_value)
    return times_value, ppvm_auto


for dt_value in dt_values:
    for gamma_value in tqdm(gamma_values, desc=f"dt={dt_value:g}"):
        times_value, ppvm_auto = ppvm_autocorrelator(gamma_value, dt_value)
        exact_auto = exact_bilinear_correlator(L, J, gamma_value, times_value, site0)[:, site0]
        rel = np.abs(ppvm_auto - exact_auto) / np.maximum(np.abs(exact_auto), 1e-15)
        max_rel_errors[dt_value].append(np.max(rel))

fig, ax = plt.subplots()
for dt_value in dt_values:
    ax.loglog(gamma_values, max_rel_errors[dt_value], "o-", label=rf"$\delta t={dt_value:g}$")
ax.set_xlabel(r"$\gamma$")
ax.set_ylabel("Maximum relative autocorrelator error")
ax.set_title("Maximum ppvm error versus dephasing strength")
ax.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
#
# While smaller dt values improve accuracy for larger $\gamma$, this is not true for smaller ones - this is because we kept the truncation error fixed and at small time step, more truncations happen which dominate at small $\gamma$. We therefore recommend the following convergence procedure: fix $\delta t$, converge with respect to truncation. Then decrease $\delta t$ and repeat until time-step error negligible.
