"""Trotter tight-accuracy scan (rel 1e-4, 1e-5) at N=10, to compare vs expm.
Trotter (2nd-order Strang) error is O(dt^2) + truncation(min_abs_coeff)."""
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
def run(dt,mac):
    best=None
    for _ in range(REPS):
        r=subprocess.run([sys.executable,W,"trotter",str(L),str(T),str(dt),str(mac)],capture_output=True,text=True)
        try: d=json.loads(r.stdout.strip().splitlines()[-1])
        except Exception: continue
        d["rel"]=float(np.linalg.norm(np.array(d["profile"])-ex)/en)
        if best is None or d["wall"]<best["wall"]: best=d
    return best or dict(rel=float('nan'),wall=float('nan'),rss=float('nan'),peak=-1)
DTS=[0.00625,0.0125,0.025,0.05]; MACS=[1e-5,1e-6,1e-7]
print(f"\n### N={N} TROTTER tight (dt x min_abs_coeff), min of {REPS} ###",flush=True)
print(f"{'dt':>8} {'mac':>7} {'rel_err':>9} {'wall_s':>7} {'RSS_mb':>7} {'peak':>8}",flush=True)
rows=[]
for dt in DTS:
    for mc in MACS:
        d=run(dt,mc); rows.append((dt,mc,d["rel"],d["wall"],d["rss"],d["peak"]))
        print(f"{dt:>8} {mc:>7.0e} {d['rel']:>9.2e} {d['wall']:>7.1f} {d['rss']:>7.0f} {d['peak']:>8}",flush=True)
for tgt in (1e-4,1e-5):
    ok=[r for r in rows if r[2]<=tgt]
    if not ok: print(f"  rel<= {tgt:.0e}: NOT REACHED (finest {min(r[2] for r in rows):.1e})",flush=True); continue
    bw=min(ok,key=lambda r:r[3]); print(f"  TROTTER rel<= {tgt:.0e}: dt={bw[0]} mac={bw[1]:.0e} -> {bw[3]:.1f}s {bw[4]:.0f}MB {bw[5]} terms (rel {bw[2]:.1e})",flush=True)
