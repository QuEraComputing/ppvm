"""Pedagogical example: Trotter evolution of the XY two-leg ladder with ppvm.

We track an observable in the Heisenberg picture as a `PauliSum` (a weighted
sum of Pauli strings). Seeding with a localized Z at the centre rung and
applying second-order Trotter steps of the XX+YY ladder Hamiltonian, we read off
the spin profile a_q(t) = <Z_q, O(t)> each step and reduce it to a
mean-squared displacement (MSD) along the chain — the signature of spin
transport. Truncating tiny coefficients (`min_abs_coeff`) is what keeps the
growing operator tractable.

    python trotter_ladder.py
"""
import numpy as np
from ppvm import PauliSum

# ── parameters ───────────────────────────────────────────────────────────
L = 10              # rungs; N = 2L qubits (two legs)
dt = 0.05           # Trotter step
steps = 20
min_abs_coeff = 1e-4   # drop Pauli terms smaller than this (accuracy knob)

import sys  # optional override for benchmarking:  python trotter_ladder.py L dt steps min_abs_coeff
if len(sys.argv) == 5:
    L, dt, steps, min_abs_coeff = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])

N = 2 * L
j0 = L // 2                       # seed at the centre rung
site = lambda j, a: j + a * L     # qubit index of rung j on leg a (0 or 1)

# Brick-wall bond layering: group bonds into vertex-disjoint layers so adjacent
# gates within a sweep commute (no within-sweep propagation). The ladder is
# degree-3, so it takes up to 4 layers -- even leg-bonds, odd leg-bonds, the
# wrap-around seam (an odd ring is not 2-edge-colourable), and the rungs.
leg = lambda j, a: (site(j, a), site((j + 1) % L, a))
bonds = ([leg(j, a) for a in (0, 1) for j in range(0, L - 1, 2)]   # even leg bonds
       + [leg(j, a) for a in (0, 1) for j in range(1, L - 1, 2)]   # odd leg bonds
       + [leg(L - 1, a) for a in (0, 1)]                           # wrap seam
       + [(site(j, 0), site(j, 1)) for j in range(L)])             # rungs

# Observable O(0) = (Z_{j0,0} + Z_{j0,1}) / 2.  Σ_q a_q = <total Z> is conserved.
o = PauliSum.new(N, [(f"Z{site(j0, 0)}", 0.5), (f"Z{site(j0, 1)}", 0.5)],
                 min_abs_coeff=min_abs_coeff, max_pauli_weight=N)
z_targets = [PauliSum.new(N, f"Z{q}") for q in range(N)]   # for the profile readout


def trotter_step():
    """One second-order Strang step (O(dt^3)): a forward sweep over the bonds
    then a reversed sweep, each bond evolved by rxx(dt) then ryy(dt). rxx and ryy
    commute on a bond, so the reversed bond order makes the step a palindrome ->
    O(dt^3). rxx is applied WITHOUT truncation and ryy WITH it, on purpose: rxx
    turns Z_a into cos(.)Z_a + sin(.)Y_aX_b -- the hopping intermediate -- which
    the very next ryy on the SAME bond rotates into Z_b. Dropping terms only
    after that pair (never the bare intermediate) keeps the Z-hopping flux
    intact, so total Z stays conserved even under aggressive truncation."""
    for a, b in bonds:
        o.rxx(a, b, theta=dt, truncate=False); o.ryy(a, b, theta=dt)
    for a, b in reversed(bonds):
        o.rxx(a, b, theta=dt, truncate=False); o.ryy(a, b, theta=dt)


def msd():
    """MSD along the chain: Σ_q (Δj_q)^2 a_q / Σ_q a_q, with a_q = <Z_q, O(t)>,
    Δj the min-image chain displacement from j0 (both legs summed via q % L)."""
    a = np.array([o.overlap(z) for z in z_targets])
    dj = ((np.arange(N) % L) - j0 + L // 2) % L - L // 2
    return float((dj**2 * a).sum() / a.sum())


print(" t       MSD     n_terms")
for n in range(steps + 1):
    print(f"{n*dt:5.2f}  {msd():8.4f}  {len(o):7d}")
    if n < steps:
        trotter_step()
