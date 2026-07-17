"""Exact reference for the XY-ladder REAL-SPACE profile + MSD, exploiting
total-magnetization conservation (the quspin-style sector reduction).

a_x(t) = (1/D) Tr[Z_x Z_seed(t)], seed = (Z_{j0,0}+Z_{j0,1})/2, j0 = L//2.
H = XX+YY on all ladder bonds conserves total Z, so H, U(t), and the diagonal
Z observables are block-diagonal in the magnetization sectors (dim <= C(N,N/2),
e.g. 3432 at N=14 instead of 16384). Per sector: dense eigh (trivial at these
dims), U_m = V e^{iEt} V^T, and
   a_x(t) += (1/D) d_x[m]^T |U_m(t)|^2 d_seed[m].
Validated against the full-space |U|^2 bilinear at L=5 (which itself was
validated against kron ED).

Usage: exact_msd_L7.py L T M out.npz     (times t_m = m*T/M)
"""
import sys, time
import numpy as np

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

L = int(sys.argv[1])
T, M = float(sys.argv[2]), int(sys.argv[3])
out = sys.argv[4]
N, bonds = ladder(L)
dim = 2 ** N
j0 = L // 2

t0 = time.time()
idx = np.arange(dim)
pop = np.array([bin(s).count("1") for s in range(dim)])
d = np.array([1 - 2 * ((idx >> (N - 1 - q)) & 1) for q in range(N)], dtype=float)  # (N, dim)
dseed = 0.5 * d[j0] + 0.5 * d[L + j0]

ts = np.array([m * T / M for m in range(M + 1)])
a = np.zeros((M + 1, N))
for nup in range(N + 1):
    sec = np.flatnonzero(pop == nup)
    dm = len(sec)
    pos = {int(s): i for i, s in enumerate(sec)}
    Hm = np.zeros((dm, dm))
    for (p, q) in bonds:
        bp = (sec >> (N - 1 - p)) & 1
        bq = (sec >> (N - 1 - q)) & 1
        differ = np.flatnonzero(bp != bq)
        for i in differ:
            dst = int(sec[i]) ^ (1 << (N - 1 - p)) ^ (1 << (N - 1 - q))
            Hm[pos[dst], i] += 2.0
    E, V = np.linalg.eigh(Hm)
    Vt = np.ascontiguousarray(V.T)
    ds_m = dseed[sec]
    dx_m = d[:, sec]                              # (N, dm)
    for m, t in enumerate(ts):
        ph = np.exp(1j * E * t)
        P = (V * ph.real) @ Vt
        P = P ** 2
        P += ((V * ph.imag) @ Vt) ** 2
        a[m] += dx_m @ (P @ ds_m) / dim
    print(f"sector Nup={nup} dim={dm} done ({time.time()-t0:.0f}s)", flush=True)

jcoord = np.arange(N) % L
dj = (jcoord - j0) % L
dj = np.where(dj > L // 2, dj - L, dj).astype(float)
msd = (a * dj[None, :] ** 2).sum(axis=1) / a.sum(axis=1)
np.savez(out, ts=ts, a=a, msd=msd, j0=j0, L=L)
print("sum_x a_x(T) =", a[-1].sum(), "(should be 1)", flush=True)
print("MSD(t):", np.round(msd, 5), flush=True)
print(f"total {time.time()-t0:.0f}s", flush=True)
