"""Orbit-preserving (momentum-k) Trotter vs expm(cache/stream), matched rel~1e-3.
worker: orbit_bench.py <expm|trotter> L k T dt knob max_basis  (PPVM_EXPM_STREAM=1 for stream)
driver: orbit_bench.py driver
Exact-ED referenced at L=5 (N=10)."""
import sys, os, json, time, resource
from functools import reduce
import numpy as np
from ppvm import Lindbladian, PauliSum
from ppvm.lindblad import _basis_to_codes
from ppvm._core import TranslationGroup, canonicalize_basis_arr_complex

def ladder(L):
    N=2*L; site=lambda j,leg: leg*L+j; bonds=[]
    for leg in(0,1):
        for j in range(L): bonds.append((site(j,leg),site((j+1)%L,leg)))
    for j in range(L): bonds.append((site(j,0),site(j,1)))
    return N,bonds
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024*1024)

def run_expm(L,k,T,dt,drop,mb):
    N,bonds=ladder(L); phi=2*np.pi*k/L; steps=round(T/dt); snap=max(1,steps//10)
    h=[]
    for(p,q)in bonds:
        for P in"XY": s=["I"]*N; s[p]=P; s[q]=P; h.append(("".join(s),1.0))
    Lop=Lindbladian(N,h,[]); g=TranslationGroup.ladder(L,2); mom=np.array([k],dtype=np.int32)
    zs=["I"*q+"Z"+"I"*(N-q-1) for q in range(N)]; zb=_basis_to_codes(zs,N)
    seed=np.array([np.exp(-1j*phi*(q%L)) for q in range(N)],dtype=np.complex128)
    bo,coo=canonicalize_basis_arr_complex(zb,seed,g,mom); pr=bo.copy()
    sb,sc=bo.copy(),coo.copy(); den=np.sum(np.conj(sc)*sc)
    def Ck(bo,coo):
        idxm={bo[i].tobytes():i for i in range(len(bo))}; num=0j
        for i in range(len(sb)):
            key=sb[i].tobytes()
            if key in idxm: num+=np.conj(sc[i])*coo[idxm[key]]
        return num/den
    curve=[Ck(bo,coo)]; peak=len(bo); t0=time.time()
    for st in range(steps):
        bo,coo=Lop.pc_step_orbit_rep(bo,coo,dt,mb,g,mom,drop,pr); peak=max(peak,len(bo))
        if (st+1)%snap==0: curve.append(Ck(bo,coo))
    return curve,time.time()-t0,peak

def run_trotter(L,k,T,dt,mac,mb):
    N,bonds=ladder(L); phi=2*np.pi*k/L; steps=round(T/dt); snap=max(1,steps//10)
    g=TranslationGroup.ladder(L,2); mom=[k]
    def mk(fn): return PauliSum.new(N,[(f"Z{q}",float(fn(phi*(q%L)))) for q in range(N)],min_abs_coeff=mac,max_pauli_weight=N)
    R=mk(np.cos); Im=mk(lambda x:-np.sin(x)); Rs=mk(np.cos); Ims=mk(lambda x:-np.sin(x))
    Rs.momentum_merge(Ims,g,mom); norm=Rs.overlap(Rs)+Ims.overlap(Ims)
    def Ck(R,Im):
        re=Rs.overlap(R)+Ims.overlap(Im); im=Rs.overlap(Im)-Ims.overlap(R); return (re+1j*im)/norm
    def snapCk():
        Rc=R.copy(); Ic=Im.copy(); Rc.momentum_merge(Ic,g,mom); return Ck(Rc,Ic)
    curve=[snapCk()]; peak=len(R); t0=time.time()
    for st in range(steps):
        for a,b in bonds: R.rxx(a,b,dt,truncate=False); R.ryy(a,b,dt,truncate=False)
        for a,b in reversed(bonds): R.rxx(a,b,dt,truncate=False); R.ryy(a,b,dt,truncate=False)
        for a,b in bonds: Im.rxx(a,b,dt,truncate=False); Im.ryy(a,b,dt,truncate=False)
        for a,b in reversed(bonds): Im.rxx(a,b,dt,truncate=False); Im.ryy(a,b,dt,truncate=False)
        R.truncate(); Im.truncate(); peak=max(peak,len(R))
        if (st+1)%snap==0: curve.append(snapCk())
    return curve,time.time()-t0,peak

if len(sys.argv)>1 and sys.argv[1]!="driver":
    m=sys.argv[1];L=int(sys.argv[2]);k=int(sys.argv[3]);T=float(sys.argv[4]);dt=float(sys.argv[5]);knob=float(sys.argv[6]);mb=int(sys.argv[7])
    curve,wall,peak=(run_expm if m=="expm" else run_trotter)(L,k,T,dt,knob,mb)
    print(json.dumps(dict(curve=[[c.real,c.imag] for c in curve],wall=wall,peak=peak,rss=rss()))); sys.exit(0)

# ---- driver (L=5, exact-ED reference) ----
import subprocess
W=os.path.abspath(__file__); PY=sys.executable
L,k,T=5,2,2.0
# exact C_k(t_m), t_m = m*T/10
SX=np.array([[0,1],[1,0]],complex);SY=np.array([[0,-1j],[1j,0]],complex);I2=np.eye(2)
N,bonds=ladder(L); phi=2*np.pi*k/L
def op(s,q):return reduce(np.kron,[s if kk==q else I2 for kk in range(N)])
H=np.zeros((2**N,2**N),complex)
for p,q in bonds: H+=op(SX,p)@op(SX,q)+op(SY,p)@op(SY,q)
E,V=np.linalg.eigh(H); idx=np.arange(2**N)
d=np.zeros(2**N,complex)
for q in range(N): d+=np.exp(-1j*phi*(q%L))*(1-2*((idx>>(N-1-q))&1))
REF=np.array([ (np.conj(d)@((np.abs((V*np.exp(1j*E*(m*T/10)))@V.conj().T)**2)@d))/(np.conj(d)@d) for m in range(11)])
def call(method,dt,knob,mb,stream=False):
    env={**os.environ}
    if stream: env["PPVM_EXPM_STREAM"]="1"
    else: env.pop("PPVM_EXPM_STREAM",None)
    r=subprocess.run([PY,W,method,str(L),str(k),str(T),str(dt),str(knob),str(mb)],capture_output=True,text=True,env=env)
    try: d=json.loads(r.stdout.strip().splitlines()[-1])
    except Exception: return None
    d["curve"]=np.array([c[0]+1j*c[1] for c in d["curve"]]); return d
def rel(c):
    n=min(len(c),len(REF)); return float(np.linalg.norm(c[:n]-REF[:n])/np.linalg.norm(REF[:n]))
print(f"ORBIT k-RESOLVED BENCH: ladder L={L}(N={2*L}) k={k} T={T}  (exact-ED reference)",flush=True)
print(f"{'method':14} {'dt':>5} {'knob':>7} {'max_basis':>9} {'rel_err':>9} {'wall_s':>7} {'RSS_mb':>7} {'peak':>7}",flush=True)
PLAN=[
 ("expm-cache",0.1,1e-4,10**7,0),("expm-cache",0.05,1e-4,10**7,0),("expm-cache",0.025,1e-4,10**7,0),
 ("trotter",0.1,1e-4,10**7,0),("trotter",0.05,1e-4,10**7,0),("trotter",0.025,1e-4,10**7,0),
 ("expm-cache",0.05,3e-4,10**7,0),("trotter",0.05,3e-4,10**7,0),
]
for name,dt,knob,mb,stream in PLAN:
    method="trotter" if name=="trotter" else "expm"
    d=call(method,dt,knob,mb,bool(stream))
    if d is None: print(f"{name:14} {dt:>5} {knob:>7.0e} {mb:>9} FAILED",flush=True); continue
    print(f"{name:14} {dt:>5} {knob:>7.0e} {mb:>9} {rel(d['curve']):>9.2e} {d['wall']:>7.1f} {d['rss']:>7.0f} {d['peak']:>7}",flush=True)
