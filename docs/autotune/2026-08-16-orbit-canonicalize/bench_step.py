"""Per-rep cost of one pc_step: momentum (orbit-rep, k=1) vs real space,
at matched live basis size B. Prints us per basis entry."""
import sys, time
import numpy as np
from ppvm import Lindbladian
from ppvm.lindblad import _basis_to_codes
from ppvm._core import TranslationGroup, canonicalize_basis_arr_complex

tag = sys.argv[1]
B = 30000
DT = 0.1

def _time(fn):
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


def build(L):
    N = 2 * L
    site = lambda j, a: j + a * L
    def pstr(pairs):
        s = ["I"] * N
        for q, o in pairs:
            s[q] = o
        return "".join(s)
    h = []
    for a in (0, 1):
        for j in range(L):
            for o in "XY":
                h.append((pstr([(site(j, a), o), (site((j + 1) % L, a), o)]), 1.0))
    for j in range(L):
        for o in "XY":
            h.append((pstr([(site(j, 0), o), (site(j, 1), o)]), 1.0))
    Z = [pstr([(q, "Z")]) for q in range(N)]
    return N, h, Z

for L in (11, 16):
    N, h, Z = build(L)
    op = Lindbladian(N, h, [])
    group = TranslationGroup.ladder(L, 2)

    # ---- momentum, k=1 -------------------------------------------------
    mom = np.array([1], dtype=np.int32)
    ph = np.array([np.exp(-2j * np.pi * (q % L) / L) for q in range(N)])
    basis, coeffs = canonicalize_basis_arr_complex(_basis_to_codes(Z, N), ph, group, mom)
    for _ in range(60):
        if len(basis) >= B:
            break
        basis, coeffs = op.pc_step_orbit_rep(basis, coeffs, DT, B, group=group,
                                             momentum=mom, drop_tol=0.0)
    nb_k = len(basis)
    tk = min(_time(lambda: op.pc_step_orbit_rep(basis, coeffs, DT, B, group=group,
                                                momentum=mom, drop_tol=0.0)) for _ in range(3))

    # ---- real space ----------------------------------------------------
    # single-Z seed: the uniform sum is conserved (L*(ΣZ)=0) and never grows
    rb = _basis_to_codes([Z[0]], N)
    rc = np.ones(1, dtype=np.float64)
    for _ in range(60):
        if len(rb) >= B:
            break
        rb, rc = op.pc_step_arr(rb, rc, DT, B, drop_tol=0.0)
    nb_r = len(rb)
    tr = min(_time(lambda: op.pc_step_arr(rb, rc, DT, B, drop_tol=0.0)) for _ in range(3))

    print(f"{tag} L={L:<3} momentum k=1: basis={nb_k:6d}  step={tk:6.3f}s -> {tk/nb_k*1e6:7.2f} us/rep", flush=True)
    print(f"{tag} L={L:<3} real space  : basis={nb_r:6d}  step={tr:6.3f}s -> {tr/nb_r*1e6:7.2f} us/term", flush=True)
    print(f"{tag} L={L:<3} momentum/real per-entry ratio: {(tk/nb_k)/(tr/nb_r):.2f}x", flush=True)
