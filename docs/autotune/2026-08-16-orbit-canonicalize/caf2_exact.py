import sys
sys.path.insert(0, "/Users/alexschuckert/dev/26_ppvm/CTPP Figures/fig16_caf2_bulk_fid")
import numpy as np, model as M
# Uncapped (cap far above the reached basis) => truncation-free: orbit and real
# must agree to round-off, which tests the rotation end to end.
L, NS = 3, 2
t, g_orb, nb, N, w1 = M.fid_orbit(L, 20_000_000, dt_us=2.0, nsteps=NS)
t, g_real, nbr, N, w2 = M.fid_real(L, 20_000_000, dt_us=2.0, nsteps=NS)
print(f"orbit peak {nb.max():,} reps ({w1:.1f}s) | real peak {nbr.max():,} strings ({w2:.1f}s)"
      f" | compression {nbr.max()/nb.max():.1f}x")
print(f"max |G_orbit - G_real| = {np.abs(g_orb-g_real).max():.3e}")
