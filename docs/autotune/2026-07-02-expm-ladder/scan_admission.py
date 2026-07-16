"""Fair admission-control scan: admit-all vs K-filter (PPVM_K_LEAKAGE) vs
max_basis cap, on TWO models (2-leg XY ladder AND long-range XY chain, alpha=1),
exact-ED referenced. Prints (rel_err, wall, RSS, peak) per method so the Pareto
fronts (wall-vs-rel, rss-vs-rel) can be compared.

worker:  scan_admission.py <model> <N> <T> <dt> <max_basis> <drop>   (K from env PPVM_K_LEAKAGE)
driver:  scan_admission.py                                            (runs the sweep)
"""
import sys, os, json, time, resource
from functools import reduce
import numpy as np, scipy.sparse as sp
from ppvm import Lindbladian
from ppvm.lindblad import _basis_to_codes

SX=np.array([[0,1],[1,0]],complex); SY=np.array([[0,-1j],[1j,0]],complex)

def ladder_bonds(N):
    L=N//2; site=lambda j,a:j+a*L; b=[]
    for a in(0,1):
        for j in range(L): b.append((site(j,a),site((j+1)%L,a),1.0))
    for j in range(L): b.append((site(j,0),site(j,1),1.0))
    return b, site(L//2,0), site(L//2,1)          # bonds(with coupling), two seed sites

def lrxy_bonds(N, alpha=1.0):
    b=[]
    for i in range(N):
        for j in range(i+1,N):
            r=min(abs(i-j), N-abs(i-j))            # periodic chain distance
            b.append((i,j,1.0/r**alpha))
    c=N//2
    return b, c, c                                 # seed = middle site (both = c)

def model(name, N):
    if name=="ladder": return ladder_bonds(N)
    if name=="lrxy":   return lrxy_bonds(N, 1.0)
    raise ValueError(name)

def sop(s,q,N): return reduce(lambda x,y:sp.kron(x,y,format='csr'),[sp.csr_matrix(s) if k==q else sp.identity(2,format='csr',dtype=complex) for k in range(N)])
def exact(name,N,T):
    bonds,s0,s1=model(name,N); dim=2**N; H=sp.csr_matrix((dim,dim),dtype=complex)
    for(p,q,J)in bonds: H=H+J*(sop(SX,p,N)@sop(SX,q,N)+sop(SY,p,N)@sop(SY,q,N))
    E,V=np.linalg.eigh(H.toarray()); U=(V*np.exp(1j*E*T))@V.conj().T
    idx=np.arange(dim)
    o0=0.5*((1-2*((idx>>(N-1-s0))&1))+(1-2*((idx>>(N-1-s1))&1)))
    dO=(np.abs(U)**2)@o0
    return np.array([((1-2*((idx>>(N-1-q))&1))*dO).sum()/dim for q in range(N)])

def run_expm(name,N,T,dt,max_basis,drop):
    bonds,s0,s1=model(name,N); steps=round(T/dt); h=[]
    for(p,q,J)in bonds:
        for P in"XY": s=["I"]*N; s[p]=P; s[q]=P; h.append(("".join(s),J))
    Lop=Lindbladian(N,h,[]); zs=["I"*q+"Z"+"I"*(N-q-1) for q in range(N)]
    zb=_basis_to_codes(zs,N); co=np.zeros(N); co[s0]+=0.5; co[s1]+=0.5
    ba=zb.copy(); pr=zb.copy(); peak=len(ba)
    for _ in range(steps):
        ba,co=Lop.pc_step_arr(ba,co,dt,max_basis=max_basis,drop_tol=drop,protected_arr=pr); peak=max(peak,len(ba))
    keys=[zb[q].tobytes() for q in range(N)]; idx={ba[i].tobytes():i for i in range(len(ba))}
    return np.array([co[idx[k]] if k in idx else 0.0 for k in keys]),peak

if len(sys.argv)>1:  # worker
    name=sys.argv[1]; N=int(sys.argv[2]); T=float(sys.argv[3]); dt=float(sys.argv[4]); mb=int(sys.argv[5]); dr=float(sys.argv[6])
    t0=time.time(); a,pk=run_expm(name,N,T,dt,mb,dr)
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024*1024)
    print(json.dumps(dict(profile=a.tolist(),wall=time.time()-t0,peak=pk,rss=rss))); sys.exit(0)

import subprocess
W=os.path.abspath(__file__); T=2.0; DT=0.05; N=10
INF=10_000_000
def call(name,mb,dr,K):
    env={**os.environ, "PPVM_K_LEAKAGE":str(K)}
    best=None
    for _ in range(2):
        r=subprocess.run([sys.executable,W,name,str(N),str(T),str(DT),str(mb),str(dr)],capture_output=True,text=True,env=env)
        try: d=json.loads(r.stdout.strip().splitlines()[-1])
        except Exception: continue
        if best is None or d["wall"]<best["wall"]: best=d
    return best

# (label, max_basis, K, [drop values])  — each row traces part of a method's front
PLAN=[
 ("admit-all",      INF,   0.0, [3e-3,1e-3,3e-4,1e-4]),
 ("K=0.3",          INF,   0.3, [1e-3,3e-4,1e-4,3e-5]),
 ("K=1",            INF,   1.0, [1e-3,3e-4,1e-4,3e-5]),
 ("K=3",            INF,   3.0, [3e-4,1e-4,3e-5,1e-5]),
 ("K=5",            INF,   5.0, [3e-4,1e-4,3e-5,1e-5]),
 ("mb=30k",         30000, 0.0, [3e-4,1e-4]),
 ("mb=60k",         60000, 0.0, [3e-4,1e-4]),
]
for name in ("ladder","lrxy"):
    ex=exact(name,N,T); en=np.linalg.norm(ex)
    print(f"\n################ model={name}  N={N} T={T} dt={DT} (min of 2) ################",flush=True)
    print(f"{'method':10} {'max_basis':>9} {'K':>4} {'drop':>7} {'rel_err':>9} {'wall_s':>7} {'RSS_mb':>7} {'peak':>8}",flush=True)
    for label,mb,K,drops in PLAN:
        for dr in drops:
            d=call(name,mb,dr,K)
            if d is None or d["profile"] is None:
                print(f"{label:10} {mb:>9} {K:>4} {dr:>7.0e}   FAILED",flush=True); continue
            rel=float(np.linalg.norm(np.array(d["profile"])-ex)/en)
            print(f"{label:10} {mb:>9} {K:>4} {dr:>7.0e} {rel:>9.2e} {d['wall']:>7.2f} {d['rss']:>7.0f} {d['peak']:>8}",flush=True)
