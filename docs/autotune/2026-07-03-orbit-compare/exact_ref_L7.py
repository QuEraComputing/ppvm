"""Exact ED reference for the L=7 (N=14) XY ladder momentum-k Z autocorrelation.

C_k(t_m) = <O_k, e^{tL} O_k> / <O_k, O_k>,  O_k = sum_j e^{-i phi (j%L)} Z_j,
phi = 2 pi k / L, t_m = m*T/M.  Same |U|^2 bilinear form as the L=5 driver in
orbit_bench.py (verified there vs both methods), with two speed tricks:
- H for the XY ladder is REAL symmetric -> real dsyevd (fast, multithreaded),
- U_block = (V_b * e^{iEt}) @ V^T computed as two real dgemms (re/im parts),
  blocked over rows so peak RAM stays ~ 2|V| + one block.

Usage: exact_ref_L7.py k T M out.npz
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


L = 7
k, T, M = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
out = sys.argv[4]
N, bonds = ladder(L)
phi = 2 * np.pi * k / L
dim = 2 ** N

t0 = time.time()
# XX+YY on a bond: flips the two bits where they differ, matrix element +2 (real).
H = np.zeros((dim, dim))
idx = np.arange(dim)
for (p, q) in bonds:
    bp = (idx >> (N - 1 - p)) & 1
    bq = (idx >> (N - 1 - q)) & 1
    src = idx[bp != bq]
    dst = src ^ (1 << (N - 1 - p)) ^ (1 << (N - 1 - q))
    H[dst, src] += 2.0
print(f"H built ({time.time()-t0:.0f}s), symmetric: {np.allclose(H, H.T)}", flush=True)

t1 = time.time()
E, V = np.linalg.eigh(H)   # real symmetric -> dsyevd
del H
print(f"eigh done ({time.time()-t1:.0f}s)", flush=True)

d = np.zeros(dim, complex)
for q in range(N):
    d += np.exp(-1j * phi * (q % L)) * (1 - 2 * ((idx >> (N - 1 - q)) & 1))
den = (np.conj(d) @ d).real

ts = np.array([m * T / M for m in range(M + 1)])
res = np.zeros(len(ts), complex)
Vt = np.ascontiguousarray(V.T)
BLK = 2048
t2 = time.time()
for lo in range(0, dim, BLK):
    hi = min(lo + BLK, dim)
    Vb = V[lo:hi, :]
    for m, t in enumerate(ts):
        ph = np.exp(1j * E * t)
        Wre = Vb * ph.real
        Wim = Vb * ph.imag
        P = (Wre @ Vt) ** 2
        P += (Wim @ Vt) ** 2          # P = |U_block|^2, real
        res[m] += np.conj(d[lo:hi]) @ (P @ d)
    print(f"rows {hi}/{dim} ({time.time()-t2:.0f}s)", flush=True)
res /= den
np.savez(out, ts=ts, ref=res, k=k, T=T, L=L)
print("REF real:", np.round(res.real, 6), flush=True)
print("REF |imag| max:", np.abs(res.imag).max(), flush=True)
print(f"total {time.time()-t0:.0f}s", flush=True)
