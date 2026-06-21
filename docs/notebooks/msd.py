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
# # Magic State Distillation with the Generalized Stabilizer Tableau
#
# Here we simulate an 85-qubit MSD circuit (5 code blocks of 17 qubits each) using ppvm's
# `GeneralizedTableau`. It is loosely based on [the TSIM example](https://bloqade.quera.com/latest/digital/examples/tsim/magic_state_distillation/),
# but with fewer details.

# %%
import time

from ppvm import GeneralizedTableau, MeasurementResult

QUBITS_PER_CODE_BLOCK = 17


def encode(tab: GeneralizedTableau, qubits: list[int]) -> None:
    """Apply the 17-qubit surface code encoding circuit."""
    for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16]:
        tab.sqrt_y(qubits[i])

    for i, j in [[1, 3], [7, 10], [12, 14], [13, 16]]:
        tab.cz(qubits[i], qubits[j])
    for i in [7, 16]:
        tab.sqrt_y_adj(qubits[i])
    for i, j in [[4, 7], [8, 10], [11, 14], [15, 16]]:
        tab.cz(qubits[i], qubits[j])
    for i in [4, 10, 14, 16]:
        tab.sqrt_y_adj(qubits[i])
    for i, j in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]]:
        tab.cz(qubits[i], qubits[j])
    for i in [3, 6, 9, 10, 12, 13]:
        tab.sqrt_y(qubits[i])
    for i, j in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]]:
        tab.cz(qubits[i], qubits[j])
    for i in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14]:
        tab.sqrt_y(qubits[i])
    for i, j in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]]:
        tab.cz(qubits[i], qubits[j])
    for i in [0, 2, 5, 6, 8, 10, 12]:
        tab.sqrt_y_adj(qubits[i])


def build_msd(tab) -> None:
    """Build the full 85-qubit MSD circuit (no measurement).

    Works on any tableau exposing the standard gate methods — both
    ``GeneralizedTableau`` and ``GeneralizedTableauSum``.
    """
    n_qubits = QUBITS_PER_CODE_BLOCK * 5
    qubit_addrs = list(range(n_qubits))

    # Split into 5 code blocks
    ql = [
        qubit_addrs[i * QUBITS_PER_CODE_BLOCK : (i + 1) * QUBITS_PER_CODE_BLOCK]
        for i in range(5)
    ]

    # Prepare magic state in each block: H + T on encoding qubit, then encode
    for q in ql:
        encoding_qubit = q[7]
        tab.h(encoding_qubit)
        tab.t(encoding_qubit)
        encode(tab, q)

    # Cross-block entangling operations
    for i in [0, 1, 4]:
        for q in ql[i]:
            tab.sqrt_x(q)

    for control, target in zip(ql[0], ql[1]):
        tab.cz(control, target)
    for control, target in zip(ql[2], ql[3]):
        tab.cz(control, target)

    for q in ql[0]:
        tab.sqrt_y(q)
    for q in ql[3]:
        tab.sqrt_y(q)

    for control, target in zip(ql[0], ql[2]):
        tab.cz(control, target)
    for control, target in zip(ql[3], ql[4]):
        tab.cz(control, target)

    for q in ql[0]:
        tab.sqrt_x_adj(q)

    for control, target in zip(ql[0], ql[4]):
        tab.cz(control, target)
    for control, target in zip(ql[1], ql[3]):
        tab.cz(control, target)

    for i in range(5):
        for q in ql[i]:
            tab.sqrt_x_adj(q)


def msd_circuit(tab: GeneralizedTableau) -> list[MeasurementResult]:
    """Build and measure the full 85-qubit MSD circuit."""
    build_msd(tab)
    n_qubits = QUBITS_PER_CODE_BLOCK * 5
    return [tab.measure(i) for i in range(n_qubits)]


# %% [markdown]
# ## Running the circuit
#
# Each shot requires its own copy of the initial tableau since measurement mutates the state.
# We use `fork()` to create independent copies with separate RNG streams.

# %%
n_qubits = QUBITS_PER_CODE_BLOCK * 5
n_shots = 1000

tab = GeneralizedTableau(n_qubits)

start = time.perf_counter()

results = []
for shot in range(n_shots):
    tab_shot = tab.fork(seed=shot)
    results.append(msd_circuit(tab_shot))

elapsed = time.perf_counter() - start

print(f"Simulated {n_shots} shots of the {n_qubits}-qubit MSD circuit")
print(f"Total time: {elapsed:.2f} s ({elapsed / n_shots * 1e3:.2f} ms per shot)")

# %% [markdown]
# Let's look at the measurement outcomes. Each shot produces an
# 85-bit string (one bit per measured qubit). We summarise rather
# than dump every shot — printing all `n_shots` strings would balloon
# the rendered notebook page without telling the reader anything
# new.

# %%
from collections import Counter

bitstrings = [
    "".join("1" if r == MeasurementResult.ONE else "0" for r in shot)
    for shot in results
]
counts = Counter(bitstrings)

print(f"Distinct outcomes: {len(counts)} / {n_shots}")
print("Most common 5 patterns:")
for pattern, count in counts.most_common(5):
    print(f"  {count:>4}×  {pattern}")

# %% [markdown]
# ## Sampling a noisy circuit with `GeneralizedTableauSum`
#
# The approach above re-runs the entire circuit for every shot. When the
# circuit is *noisy*, `GeneralizedTableauSum` offers an alternative: it tracks
# the full **mixed state** as a probability-weighted sum of stabilizer branches
# (one per error outcome), so the circuit is built **once**. A `Sampler` then
# draws as many shots as you like from the resulting distribution, in parallel
# across cores.
#
# We reuse the exact same `build_msd` circuit — the summed tableau exposes the
# same gate methods — and add a depolarizing channel on every qubit so the sum
# branches into a genuine mixture. `sum_cutoff` drops branches whose probability
# weight falls below the threshold, keeping the mixture tractable.

# %%
from ppvm.generalized_tableau_sum import GeneralizedTableauSum

p_depolarize = 1e-4

tab_sum = GeneralizedTableauSum(n_qubits, sum_cutoff=1e-7, seed=0)
build_msd(tab_sum)  # same gates, on the summed tableau

for q in range(n_qubits):
    tab_sum.depolarize(q, p_depolarize)

print(f"Mixed state holds {len(tab_sum)} branches")

# %% [markdown]
# Compiling a `Sampler` snapshots the current state; we can then draw all the
# shots in a single call.

# %%
start = time.perf_counter()

sampler = tab_sum.sampler()
shots = sampler.sample_shots(n_shots)  # list[list[MeasurementResult]]

elapsed = time.perf_counter() - start
print(f"Drew {n_shots} shots from the mixture in {elapsed * 1e3:.1f} ms")

# %% [markdown]
# The shots come from the same distribution as the trajectory approach, but the
# expensive circuit construction happened only once. Summarising the same way:

# %%
sum_bitstrings = [
    "".join("1" if r == MeasurementResult.ONE else "0" for r in shot)
    for shot in shots
]
sum_counts = Counter(sum_bitstrings)

print(f"Distinct outcomes: {len(sum_counts)} / {n_shots}")
print("Most common 5 patterns:")
for pattern, count in sum_counts.most_common(5):
    print(f"  {count:>4}×  {pattern}")
