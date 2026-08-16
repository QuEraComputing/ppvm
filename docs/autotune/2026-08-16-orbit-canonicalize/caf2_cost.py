"""Per-entry cost with the REAL CaF2 dipolar H (all pairs, min image)."""
import sys, time
sys.path.insert(0, "/Users/alexschuckert/dev/26_ppvm/CTPP Figures/fig16_caf2_bulk_fid")
import numpy as np
import model as M
from ppvm import Lindbladian
from ppvm.lindblad import _basis_to_codes
from ppvm._core import TranslationGroup, canonicalize_basis_arr_complex

B, DT = 2048, 2.0 / M.TAU_US
print(f"{'L':>3} {'N':>5} {'H terms':>9} {'path':>12} {'entries':>8} {'ms/step':>9} {'us/entry':>9}", flush=True)
for L in (3, 4):
    T, N = M.terms(L, "100", None)
    lind = Lindbladian(N, T, [])
    g = TranslationGroup.torus_3d(L, L, L)
    mom = np.zeros(3, dtype=np.int32)
    X = ["I" * q + "X" + "I" * (N - q - 1) for q in range(N)]

    b, c = canonicalize_basis_arr_complex(_basis_to_codes(X, N),
                                          np.ones(N, dtype=np.complex128), g, mom)
    for _ in range(25):
        if len(b) >= B:
            break
        b, c = lind.pc_step_orbit_rep(b, c, DT, B, group=g, momentum=mom,
                                      drop_tol=0.0, admit_basis=3 * B)
    t0 = time.perf_counter()
    lind.pc_step_orbit_rep(b, c, DT, B, group=g, momentum=mom, drop_tol=0.0, admit_basis=3 * B)
    tk = time.perf_counter() - t0
    print(f"{L:>3} {N:>5} {len(T):>9,} {'orbit k=0':>12} {len(b):>8} {tk*1e3:>9.1f} {tk/len(b)*1e6:>9.1f}", flush=True)

    rb = _basis_to_codes(X, N); rc = np.ones(N)
    for _ in range(25):
        if len(rb) >= B:
            break
        rb, rc = lind.pc_step_arr(rb, rc, DT, B, drop_tol=0.0, admit_basis=3 * B)
    t0 = time.perf_counter()
    lind.pc_step_arr(rb, rc, DT, B, drop_tol=0.0, admit_basis=3 * B)
    tr = time.perf_counter() - t0
    print(f"{L:>3} {N:>5} {len(T):>9,} {'real space':>12} {len(rb):>8} {tr*1e3:>9.1f} "
          f"{tr/len(rb)*1e6:>9.1f}   ratio {(tk/len(b))/(tr/len(rb)):5.1f}x  |G|={L**3}", flush=True)
