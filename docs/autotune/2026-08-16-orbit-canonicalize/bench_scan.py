import sys, time, numpy as np
from ppvm import Lindbladian
from ppvm.lindblad import _basis_to_codes
from ppvm._core import TranslationGroup, canonicalize_basis_arr_complex
tag, B, DT = sys.argv[1], 30000, 0.1
def build(L):
    N = 2*L; site = lambda j,a: j + a*L
    def pstr(pairs):
        s = ["I"]*N
        for q,o in pairs: s[q]=o
        return "".join(s)
    h=[]
    for a in (0,1):
        for j in range(L):
            for o in "XY": h.append((pstr([(site(j,a),o),(site((j+1)%L,a),o)]),1.0))
    for j in range(L):
        for o in "XY": h.append((pstr([(site(j,0),o),(site(j,1),o)]),1.0))
    return N,h,[pstr([(q,"Z")]) for q in range(N)]
def t3(fn):
    return min((lambda: (lambda t0: (fn(), time.perf_counter()-t0)[1])(time.perf_counter()))() for _ in range(3))
print(f"{'L':>4} {'|G|':>4} {'mom us/rep':>11} {'real us/term':>13} {'ratio':>7}", flush=True)
for L in (8, 11, 16, 24, 32):
    N,h,Z = build(L); op = Lindbladian(N,h,[]); g = TranslationGroup.ladder(L,2)
    mom = np.array([1], dtype=np.int32)
    ph = np.array([np.exp(-2j*np.pi*(q%L)/L) for q in range(N)])
    b,c = canonicalize_basis_arr_complex(_basis_to_codes(Z,N), ph, g, mom)
    for _ in range(60):
        if len(b) >= B: break
        b,c = op.pc_step_orbit_rep(b,c,DT,B,group=g,momentum=mom,drop_tol=0.0)
    tk = t3(lambda: op.pc_step_orbit_rep(b,c,DT,B,group=g,momentum=mom,drop_tol=0.0)) / len(b)
    rb, rc = _basis_to_codes([Z[0]],N), np.ones(1)
    for _ in range(60):
        if len(rb) >= B: break
        rb,rc = op.pc_step_arr(rb,rc,DT,B,drop_tol=0.0)
    tr = t3(lambda: op.pc_step_arr(rb,rc,DT,B,drop_tol=0.0)) / len(rb)
    print(f"{L:>4} {L:>4} {tk*1e6:>11.2f} {tr*1e6:>13.2f} {tk/tr:>6.2f}x", flush=True)
