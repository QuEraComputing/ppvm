"""
Stim-style sampling.

A common workflow for researchers: write a Stim program describing a small
QEC circuit, sample many shots, and inspect the per-shot measurement
outcomes. Loss-aware measurement results come out of ppvm as
``MeasurementResult`` values (``ZERO``, ``ONE``, or ``LOST``).
"""

from collections import Counter

from ppvm import StimProgram, sample_stim

# A two-qubit GHZ-preparation followed by measurement on both qubits.
prog = StimProgram.parse(
    """
    H 0
    CX 0 1
    M 0 1
    """
)

shots = sample_stim(prog, n_qubits=2, num_shots=200, seed=42)

# Each shot is a list of MeasurementResult; cast to int for tallying.
patterns = Counter(tuple(int(r) for r in shot) for shot in shots)
for pattern, count in sorted(patterns.items()):
    print(f"{pattern}: {count}")

# A GHZ state should only produce (0, 0) or (1, 1) — never (0, 1) or (1, 0).
correlated_fraction = sum(c for p, c in patterns.items() if p[0] == p[1]) / len(shots)
print(f"correlated fraction: {correlated_fraction}")
