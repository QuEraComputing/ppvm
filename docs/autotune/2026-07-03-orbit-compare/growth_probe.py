"""Basis-growth feasibility probe for the orbit-rep expm path at scale.

Prints one line per step: step, t, basis size, step wall, cumulative wall, RSS.
Usage: growth_probe.py L k T dt drop max_basis
"""
import sys, time, resource
import numpy as np
from ppvm import Lindbladian
from ppvm.lindblad import _basis_to_codes
from ppvm._core import TranslationGroup, canonicalize_basis_arr_complex


def ladder(L):
    N = 2 * L
    site = lambda j, leg: leg * L + j
    bonds = []
    for leg in (0, 1):
        for j in range(L):
            bonds.append((site(j, leg), site((j + 1) % L, leg)))
    for j in range(L):
        bonds.append((site(j, 0), site(j, 1)))
    return N, bonds


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


L, k, T, dt, drop, mb = (int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]),
                         float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]))
N, bonds = ladder(L)
phi = 2 * np.pi * k / L
steps = round(T / dt)

h = []
for (p, q) in bonds:
    for P in "XY":
        s = ["I"] * N
        s[p] = P
        s[q] = P
        h.append(("".join(s), 1.0))
Lop = Lindbladian(N, h, [])
g = TranslationGroup.ladder(L, 2)
mom = np.array([k], dtype=np.int32)
zs = ["I" * q + "Z" + "I" * (N - q - 1) for q in range(N)]
zb = _basis_to_codes(zs, N)
seed = np.array([np.exp(-1j * phi * (q % L)) for q in range(N)], dtype=np.complex128)
bo, coo = canonicalize_basis_arr_complex(zb, seed, g, mom)
pr = bo.copy()

print(f"# probe: ladder L={L} (N={N}) k={k} T={T} dt={dt} drop={drop:.0e} mb={mb}", flush=True)
print(f"{'step':>4} {'t':>6} {'basis':>9} {'dwall_s':>8} {'wall_s':>8} {'RSS_mb':>8}", flush=True)
t0 = time.time()
for st in range(steps):
    ts = time.time()
    bo, coo = Lop.pc_step_orbit_rep(bo, coo, dt, mb, g, mom, drop, pr)
    now = time.time()
    print(f"{st+1:>4} {(st+1)*dt:>6.2f} {len(bo):>9} {now-ts:>8.1f} {now-t0:>8.1f} {rss_mb():>8.0f}", flush=True)
