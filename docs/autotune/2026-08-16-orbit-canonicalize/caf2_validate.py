import sys
sys.path.insert(0, "/Users/alexschuckert/dev/26_ppvm/CTPP Figures/fig16_caf2_bulk_fid")
import numpy as np, model as M
L, B, NS = 3, 512, 10
t, g_orb, nb, N, w1 = M.fid_orbit(L, B, dt_us=2.0, nsteps=NS)
t, g_real, nbr, N, w2 = M.fid_real(L, B * L**3, dt_us=2.0, nsteps=NS)
print(f"L={L} N={N}  orbit: {nb.max():,} reps in {w1:.1f}s   real: {nbr.max():,} strings in {w2:.1f}s")
print(f"max |G_orbit - G_real| = {np.abs(g_orb-g_real).max():.3e}")
print("t(us) ", np.round(t[:6], 1))
print("orbit ", np.round(g_orb[:6], 12))
print("real  ", np.round(g_real[:6], 12))
