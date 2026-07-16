"""Exact test of the revival-carrier hypothesis at N=10 (L=5 ladder).

Protocol: evolve O(0) = (Z_{j0,0}+Z_{j0,1})/2 EXACTLY to t0, decompose
O(t0) over the full 4^N Pauli basis (fast per-qubit transform), split into
the top-B strings by |coefficient| (kept) and the remainder (dust), then
evolve both parts EXACTLY (no further truncation) and read off the MSD
profile a_q(t) = Tr(Z_q O(t))/2^N, which is linear in O. Prediction: the
kept part reproduces the smooth plateau; the revival wiggles live in the
dust (individually tiny, collectively coherent coefficients).

Usage: revival_test_N10.py t0 T out_prefix
"""
import sys, time
import numpy as np

L = 5
N = 2 * L
dim = 2 ** N
t0, T = float(sys.argv[1]), float(sys.argv[2])
pref = sys.argv[3]

def ladder(L):
    site = lambda j, leg: leg * L + j
    bonds = []
    for leg in (0, 1):
        for j in range(L):
            bonds.append((site(j, leg), site((j + 1) % L, leg)))
    for j in range(L):
        bonds.append((site(j, 0), site(j, 1)))
    return bonds

idx = np.arange(dim)
H = np.zeros((dim, dim))
for (p, q) in ladder(L):
    bp = (idx >> (N - 1 - p)) & 1
    bq = (idx >> (N - 1 - q)) & 1
    src = idx[bp != bq]
    dst = src ^ (1 << (N - 1 - p)) ^ (1 << (N - 1 - q))
    H[dst, src] += 2.0
E, V = np.linalg.eigh(H)

dz = np.array([1 - 2 * ((idx >> (N - 1 - q)) & 1) for q in range(N)], dtype=float)
j0 = L // 2
O0 = np.diag(0.5 * dz[j0] + 0.5 * dz[L + j0]).astype(complex)

P4 = np.stack([np.eye(2), np.array([[0, 1], [1, 0]]),
               np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]).astype(complex)

def pauli_coeffs(O):
    """(2^N x 2^N) matrix -> (4^N,) Pauli coefficients, qubit-1 slowest."""
    Tm = O.reshape((2,) * (2 * N))
    order = []
    for k in range(N):
        order += [k, N + k]                      # pairs (i_k, j_k)
    Tm = np.transpose(Tm, order).reshape((-1,))
    Tm = Tm.reshape((1, -1))
    for k in range(N):
        # leading done-paulis x (2,2) x rest  ->  contract this qubit's (i,j)
        Tm = Tm.reshape((4 ** k, 2, 2, -1))
        c0 = (Tm[:, 0, 0, :] + Tm[:, 1, 1, :]) / 2
        c1 = (Tm[:, 0, 1, :] + Tm[:, 1, 0, :]) / 2
        c2 = 1j * (Tm[:, 0, 1, :] - Tm[:, 1, 0, :]) / 2
        c3 = (Tm[:, 0, 0, :] - Tm[:, 1, 1, :]) / 2
        Tm = np.stack([c0, c1, c2, c3], axis=1).reshape((4 ** (k + 1), -1))
    return Tm.reshape(-1)

def pauli_matrix(c):
    """inverse of pauli_coeffs."""
    Tm = c.reshape((4,) * N)
    for _ in range(N):
        Tm = np.tensordot(Tm, P4, axes=([0], [0]))   # consume leading pauli axis, append (i,j)
    order = [2 * k for k in range(N)] + [2 * k + 1 for k in range(N)]
    return np.transpose(Tm, order).reshape((dim, dim))

def evolve(O, t):
    ph = np.exp(1j * E * t)
    W = (V * ph) @ V.T
    return W.conj().T @ O @ W

t_start = time.time()
# round-trip validation on O0 evolved a bit
Otest = evolve(O0, 0.3)
cc = pauli_coeffs(Otest)
err = np.abs(pauli_matrix(cc) - Otest).max()
assert err < 1e-10, f"transform round-trip failed: {err}"
print(f"transform round-trip OK ({err:.1e})", flush=True)

Ot0 = evolve(O0, t0)
c = pauli_coeffs(Ot0)
fro = np.linalg.norm(Ot0) ** 2
csum = (np.abs(c) ** 2).sum() * dim
print(f"Parseval check: Tr(O^2)={fro:.6f}  2^N*sum|c|^2={csum:.6f}", flush=True)

mag = np.abs(c)
order = np.argsort(mag)[::-1]
pw = np.zeros(4 ** N, dtype=np.int8)
tmp = np.arange(4 ** N)
for _ in range(N):
    pw += (tmp % 4 != 0).astype(np.int8)
    tmp //= 4
tot2 = (mag ** 2).sum()
n_nonzero = int((mag > 1e-14).sum())
print(f"t0={t0}: {n_nonzero:,} strings with |c|>1e-14 (of {4**N:,})", flush=True)
for B in (100, 1000, 10000, 100000):
    kept = order[:B]
    frac = (mag[kept] ** 2).sum() / tot2
    print(f"  B={B:>7,}: |c|^2 fraction {frac:.5f}  min|c| kept {mag[kept].min():.2e}  "
          f"mean weight kept {pw[kept].mean():.2f} vs dust {pw[order[B:min(B+200000,len(order))]].mean():.2f}", flush=True)

jcoord = np.arange(N) % L
dj2 = ((jcoord - j0) % L)
dj2 = np.where(dj2 > L // 2, dj2 - L, dj2).astype(float) ** 2

ts = np.arange(t0, T + 1e-9, 0.1)
Bs = (1000, 10000, 100000)
mats = {"full": Ot0}
for B in Bs:
    ck = c.copy()
    ck[order[B:]] = 0.0
    mats[f"kept{B}"] = pauli_matrix(ck)
res = {"ts": ts}
prof = {k: np.zeros((len(ts), N)) for k in mats}
for m, t in enumerate(ts):
    s = t - t0
    ph = np.exp(1j * E * s)
    W = (V * ph) @ V.T
    Wd = W.conj().T
    for k, M in mats.items():
        d = np.real(np.diagonal(Wd @ M @ W))
        prof[k][m] = (dz @ d) / dim
for k in prof:
    a = prof[k]
    res[f"msd_{k}"] = (a * dj2[None, :]).sum(axis=1)
    res[f"suma_{k}"] = a.sum(axis=1)
np.savez(pref + ".npz", **res)
print(f"saved {pref}.npz ({time.time()-t_start:.0f}s)", flush=True)
fu = res["msd_full"]
late = ts >= t0 + 1.0
wig_f = fu[late] - fu[late].mean()
for B in Bs:
    ke = res[f"msd_kept{B}"]
    wig_k = ke[late] - ke[late].mean()
    corr = np.corrcoef(wig_k, wig_f)[0, 1]
    print(f"B={B:>7,}: late <MSD_kept>={ke[late].mean():.4f} vs exact {fu[late].mean():.4f};  "
          f"wiggle corr {corr:+.3f}, amp ratio {wig_k.std()/wig_f.std():.3f}", flush=True)
