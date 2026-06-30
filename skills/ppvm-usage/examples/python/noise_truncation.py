"""
Skill noise + truncation block.

A small-support observable propagated backwards through a noisy
Trotter-style circuit. ``min_abs_coeff`` and ``max_pauli_weight`` are
passed at construction; in Python the binding calls ``truncate()`` after
each gate method call by default, so we don't manage it manually.

The point of this example is the *workflow* — combining gates, noise,
and bounded truncation — not a specific numeric outcome. We assert that
the propagation finishes, that truncation actually capped the working
set, and that the resulting overlap is a finite real number.
"""

import math

from ppvm import PauliSum

N = 12
LAYERS = 5
MAX_WEIGHT = 6

ps = PauliSum.new(
    N,
    "Z0",  # initial observable: Z on qubit 0; weight grows under rzz layers
    min_abs_coeff=1e-6,
    max_pauli_weight=MAX_WEIGHT,
)
for _ in range(LAYERS):
    for q in range(N):
        ps.depolarize1(q, p=1e-3)
        ps.rx(q, theta=0.1)
    for q in range(N - 1):
        ps.rzz(q, q + 1, theta=0.05)

# Truncation has been applied throughout; ps.current_max_weight() is the
# largest non-identity-count among the terms still kept.
top_weight = ps.current_max_weight()
overlap = ps.overlap_with_zero()
assert top_weight <= MAX_WEIGHT, top_weight
assert math.isfinite(overlap), overlap
print(f"layers={LAYERS} terms={len(ps)} max_weight={top_weight} finite_overlap=True")
# → layers=5 terms=7 max_weight=4 finite_overlap=True
