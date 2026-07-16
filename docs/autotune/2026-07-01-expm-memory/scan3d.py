"""Wide (dt, drop_tol, max_basis) scan for expm -> optimal params per accuracy.

Uses /tmp/charac.py as the subprocess worker (fresh RSS per run). Exact-ED ref.
"""
import sys, json, subprocess
from functools import reduce
import numpy as np, scipy.sparse as sp

SX=np.array([[0,1],[1,0]],complex); SY=np.array([[0,-1j],[1j,0]],complex)
def ladder(L):
    N=2*L; site=lambda j,a:j+a*L; b=[]
    for a in(0,1):
        for j in range(L): b.append((site(j,a),site((j+1)%L,a)))
    for j in range(L): b.append((site(j,0),site(j,1)))
    return N,site,b
def sop(s,q,N): return reduce(lambda x,y:sp.kron(x,y,format='csr'),[sp.csr_matrix(s) if k==q else sp.identity(2,format='csr',dtype=complex) for k in range(N)])
def exact(L,T):
    N,site,bd=ladder(L); dim=2**N; H=sp.csr_matrix((dim,dim),dtype=complex)
    for(p,q)in bd: H=H+sop(SX,p,N)@sop(SX,q,N)+sop(SY,p,N)@sop(SY,q,N)
    E,V=np.linalg.eigh(H.toarray()); U=(V*np.exp(1j*E*T))@V.conj().T
    idx=np.arange(dim); j0=L//2; a0,b0=site(j0,0),site(j0,1)
    o0=0.5*((1-2*((idx>>(N-1-a0))&1))+(1-2*((idx>>(N-1-b0))&1))); dO=(np.abs(U)**2)@o0
    return np.array([((1-2*((idx>>(N-1-q))&1))*dO).sum()/dim for q in range(N)])

W="/tmp/charac.py"; L=5; T=2.0; REPS=2
ex=exact(L,T); en=np.linalg.norm(ex); N=2*L
def run(dt,mb,dr):
    best=None
    for _ in range(REPS):
        r=subprocess.run([sys.executable,W,"expm",str(L),str(T),str(dt),str(mb),str(dr)],capture_output=True,text=True)
        try: d=json.loads(r.stdout.strip().splitlines()[-1])
        except Exception: continue
        d["rel"]=float(np.linalg.norm(np.array(d["profile"])-ex)/en)
        if best is None or d["wall"]<best["wall"]: best=d
    return best or dict(rel=float('nan'),wall=float('nan'),rss=float('nan'),peak=-1)

DTS=[0.05,0.1,0.15,0.2,0.25,0.3]
DROPS=[3e-3,1e-3,3e-4,1e-4,3e-5]
print(f"\n################ N={N} expm (dt x drop_tol), max_basis=inf, min of {REPS} ################",flush=True)
print(f"{'dt':>5} {'drop':>7} {'rel_err':>9} {'wall_s':>7} {'RSS_mb':>7} {'peak':>8}",flush=True)
rows=[]
for dt in DTS:
    for dr in DROPS:
        d=run(dt,-1,dr); rows.append((dt,dr,d["rel"],d["wall"],d["rss"],d["peak"]))
        print(f"{dt:>5} {dr:>7.0e} {d['rel']:>9.2e} {d['wall']:>7.1f} {d['rss']:>7.0f} {d['peak']:>8}",flush=True)

print(f"\n== optimal (dt,drop) per accuracy target (max_basis=inf) ==",flush=True)
best_cfg={}
for tgt in (3e-2,1e-2,3e-3,1e-3):
    ok=[r for r in rows if r[2]<=tgt]
    if not ok: print(f"  rel<= {tgt:.0e}: none reached",flush=True); continue
    bw=min(ok,key=lambda r:r[3]); br=min(ok,key=lambda r:r[4])
    best_cfg[tgt]=bw
    print(f"  rel<= {tgt:.0e}: MIN-WALL dt={bw[0]} drop={bw[1]:.0e} -> {bw[3]:.1f}s {bw[4]:.0f}MB (rel {bw[2]:.1e}) | "
          f"MIN-RSS dt={br[0]} drop={br[1]:.0e} -> {br[4]:.0f}MB {br[3]:.1f}s (rel {br[2]:.1e})",flush=True)

# max_basis interaction on the min-wall config for a moderate & a tight target
print(f"\n== max_basis interaction (does capping help RAM without losing accuracy?) ==",flush=True)
for tgt in (1e-2,3e-3):
    if tgt not in best_cfg: continue
    dt,dr=best_cfg[tgt][0],best_cfg[tgt][1]
    print(f"  base config dt={dt} drop={dr:.0e} (target rel {tgt:.0e}):",flush=True)
    print(f"    {'max_basis':>10} {'rel_err':>9} {'wall_s':>7} {'RSS_mb':>7} {'peak':>8}",flush=True)
    for mb in (20000,50000,100000,200000,-1):
        d=run(dt,mb,dr); print(f"    {mb:>10} {d['rel']:>9.2e} {d['wall']:>7.1f} {d['rss']:>7.0f} {d['peak']:>8}",flush=True)
