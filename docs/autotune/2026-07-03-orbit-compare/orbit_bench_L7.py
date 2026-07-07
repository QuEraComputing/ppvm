"""Orbit-preserving (momentum-k) Trotter vs expm at scale: ladder L=7 (N=14), k=3, T=2.
Reference: exact ED curve from exact_ref_L7.py (npz next to this file).

worker: orbit_bench_L7.py <expm|trotter> L k T dt knob max_basis   (PPVM_EXPM_STREAM=1 for stream)
driver: orbit_bench_L7.py driver [row ...]   rows like expm:0.1:1e-3:10000000[:stream]
        (no rows -> default plan)

NOTE (2026-07-07): PPVM_K_LEAKAGE / PPVM_EXPM_STREAM were removed from
ppvm; the K / stream row tokens in this legacy driver are now inert. Use
scan_xy_mid.py / scan_realspace_msd.py (explicit --tau_add) instead.
"""
import sys, os, json, time, resource
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

def run_trotter(L,k,T,dt,mac,mb,perbond=False):
    # perbond=True is the scheme of the verified real-space baseline
    # (demo/trotter_ladder.py): rxx untruncated, each ryy truncates. The legacy
    # perbond=False variant defers all truncation to the end of the step, which
    # blows up the intra-step transient at large N.
    N,bonds=ladder(L); phi=2*np.pi*k/L; steps=round(T/dt); snap=max(1,steps//10)
    g=TranslationGroup.ladder(L,2); mom=[k]
    def mk(fn): return PauliSum.new(N,[(f"Z{q}",float(fn(phi*(q%L)))) for q in range(N)],min_abs_coeff=mac,max_pauli_weight=N)
    R=mk(np.cos); Im=mk(lambda x:-np.sin(x)); Rs=mk(np.cos); Ims=mk(lambda x:-np.sin(x))
    Rs.momentum_merge(Ims,g,mom); norm=Rs.overlap(Rs)+Ims.overlap(Ims)
    def Ck(R,Im):
        re=Rs.overlap(R)+Ims.overlap(Im); im=Rs.overlap(Im)-Ims.overlap(R); return (re+1j*im)/norm
    def snapCk():
        Rc=R.copy(); Ic=Im.copy(); Rc.momentum_merge(Ic,g,mom); return Ck(Rc,Ic)
    def sweep(o):
        if perbond:
            for a,b in bonds: o.rxx(a,b,dt,truncate=False); o.ryy(a,b,dt)
            for a,b in reversed(bonds): o.rxx(a,b,dt,truncate=False); o.ryy(a,b,dt)
        else:
            for a,b in bonds: o.rxx(a,b,dt,truncate=False); o.ryy(a,b,dt,truncate=False)
            for a,b in reversed(bonds): o.rxx(a,b,dt,truncate=False); o.ryy(a,b,dt,truncate=False)
            o.truncate()
    curve=[snapCk()]; peak=len(R); t0=time.time()
    for st in range(steps):
        sweep(R); sweep(Im)
        peak=max(peak,max(len(R),len(Im)))
        if (st+1)%snap==0: curve.append(snapCk())
    return curve,time.time()-t0,peak

if len(sys.argv)>1 and sys.argv[1]!="driver":
    m=sys.argv[1];L=int(sys.argv[2]);k=int(sys.argv[3]);T=float(sys.argv[4]);dt=float(sys.argv[5]);knob=float(sys.argv[6]);mb=int(sys.argv[7])
    if m=="expm": curve,wall,peak=run_expm(L,k,T,dt,knob,mb)
    else: curve,wall,peak=run_trotter(L,k,T,dt,knob,mb,perbond=(m=="trotter-pb"))
    print(json.dumps(dict(curve=[[c.real,c.imag] for c in curve],wall=wall,peak=peak,rss=rss()))); sys.exit(0)

# ---- driver: exact-ED-referenced from npz ----
import subprocess
W=os.path.abspath(__file__); PY=sys.executable; D=os.path.dirname(W)
L,k,T=7,int(os.environ.get("PPVM_BENCH_K","3")),2.0
npz=np.load(os.path.join(D,f"exact_ref_L7_k{k}_T2.npz"))
REF=npz["ref"]
def call(method,dt,knob,mb,stream=False,kleak=None):
    env={**os.environ}
    if stream: env["PPVM_EXPM_STREAM"]="1"
    else: env.pop("PPVM_EXPM_STREAM",None)
    if kleak is not None: env["PPVM_K_LEAKAGE"]=str(kleak)
    else: env.pop("PPVM_K_LEAKAGE",None)
    r=subprocess.run([PY,W,method,str(L),str(k),str(T),str(dt),str(knob),str(mb)],capture_output=True,text=True,env=env)
    try: d=json.loads(r.stdout.strip().splitlines()[-1])
    except Exception: return None
    d["curve"]=np.array([c[0]+1j*c[1] for c in d["curve"]]); return d
def rel(c):
    n=min(len(c),len(REF)); return float(np.linalg.norm(c[:n]-REF[:n])/np.linalg.norm(REF[:n]))
print(f"ORBIT k-RESOLVED BENCH AT SCALE: ladder L={L}(N={2*L}) k={k} T={T}  (exact-ED reference)",flush=True)
print(f"{'method':14} {'dt':>6} {'knob':>7} {'K':>4} {'max_basis':>9} {'rel_err':>9} {'wall_s':>7} {'RSS_mb':>7} {'peak':>8}",flush=True)
if len(sys.argv)>2:
    PLAN=[]
    for row in sys.argv[2:]:
        f=row.split(":")
        stream=1 if "stream" in f[4:] else 0
        kleak=next((tok[1:] for tok in f[4:] if tok.startswith("K")),None)
        PLAN.append((f[0],float(f[1]),float(f[2]),int(f[3]),stream,kleak))
else:
    PLAN=[
     ("expm-cache",0.1,1e-3,10**7,0,None),("trotter",0.1,1e-3,10**7,0,None),
     ("expm-cache",0.05,1e-3,10**7,0,None),("trotter",0.05,1e-3,10**7,0,None),
    ]
for name,dt,knob,mb,stream,kleak in PLAN:
    method=name if name.startswith("trotter") else "expm"
    d=call(method,dt,knob,mb,bool(stream),kleak)
    ktag=kleak if kleak is not None else "-"
    if d is None: print(f"{name:14} {dt:>6} {knob:>7.0e} {ktag:>4} {mb:>9} FAILED",flush=True); continue
    print(f"{name:14} {dt:>6} {knob:>7.0e} {ktag:>4} {mb:>9} {rel(d['curve']):>9.2e} {d['wall']:>7.1f} {d['rss']:>7.0f} {d['peak']:>8}",flush=True)
