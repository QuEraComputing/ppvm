import numpy as np
from functools import reduce
from ppvm import Lindbladian
from ppvm.lindblad import _basis_to_codes
from ppvm._core import TranslationGroup, canonicalize_basis_arr_complex

N=8; k=1; dt=0.02; nsnap=5; steps_per=3; drop=1e-10; MB=10_000_000
phi=2*np.pi*k/N
# --- 1D XY chain PBC ---
h=[]
for j in range(N):
    nx=(j+1)%N
    for P in "XY":
        s=["I"]*N; s[j]=P; s[nx]=P; h.append(("".join(s),1.0))
Lop=Lindbladian(N,h,[])
g=TranslationGroup.chain_1d(N); mom=np.array([k],dtype=np.int32)

# --- exact ED ref ---
SX=np.array([[0,1],[1,0]],complex);SY=np.array([[0,-1j],[1j,0]],complex);SZ=np.array([[1,0],[0,-1]],complex);I2=np.eye(2)
def op(s,q):return reduce(np.kron,[s if k2==q else I2 for k2 in range(N)])
H=np.zeros((2**N,2**N),complex)
for j in range(N):
    nx=(j+1)%N; H+=op(SX,j)@op(SX,nx)+op(SY,j)@op(SY,nx)
idx=np.arange(2**N)
d=np.zeros(2**N,complex)
for j in range(N): d+=np.exp(-1j*phi*j)*(1-2*((idx>>(N-1-j))&1))
E,V=np.linalg.eigh(H)
def exact_Ck(t):
    U=(V*np.exp(1j*E*t))@V.conj().T; U2=np.abs(U)**2
    return (np.conj(d)@(U2@d))/(np.conj(d)@d)

# --- seed (real-space k-mode) ---
zs=["I"*j+"Z"+"I"*(N-j-1) for j in range(N)]; zb=_basis_to_codes(zs,N)
seed=np.array([np.exp(-1j*phi*j) for j in range(N)],dtype=np.complex128)

# --- real-space complex expm ---
ba=zb.copy(); co=seed.copy(); pr=zb.copy()
def Ck_realspace(ba,co):
    idxm={ba[i].tobytes():i for i in range(len(ba))}
    s=0j
    for j in range(N):
        key=zb[j].tobytes()
        if key in idxm: s+=np.exp(1j*phi*j)*co[idxm[key]]
    return s/N

# --- orbit-rep expm ---
bo,coo=canonicalize_basis_arr_complex(zb.copy(),seed.copy(),g,mom)
# the Z-orbit rep row(s): track total k-mode Z amplitude via overlap with seed's rep
seed_bo, seed_coo = bo.copy(), coo.copy()
def Ck_orbit(bo,coo):
    idxm={bo[i].tobytes():i for i in range(len(bo))}
    num=0j
    for i in range(len(seed_bo)):
        key=seed_bo[i].tobytes()
        if key in idxm: num+=np.conj(seed_coo[i])*coo[idxm[key]]
    den=np.sum(np.conj(seed_coo)*seed_coo)
    return num/den
pr_o=seed_bo.copy()

print(f"N={N} k={k} dt={dt}  (near-exact: drop={drop})")
print(f"{'t':>6} {'exact':>22} {'realspace':>22} {'orbit-rep':>22}")
bo2,coo2=bo.copy(),coo.copy()
for si in range(nsnap+1):
    t=si*steps_per*dt
    ex=exact_Ck(t); rs=Ck_realspace(ba,co); orb=Ck_orbit(bo2,coo2)
    print(f"{t:6.3f} {ex.real:+.4f}{ex.imag:+.4f}j  {rs.real:+.4f}{rs.imag:+.4f}j  {orb.real:+.4f}{orb.imag:+.4f}j")
    if si<nsnap:
        for _ in range(steps_per):
            ba,co=Lop.pc_step_complex(ba,co,dt,0.0,drop,pr,None,g,mom)
            bo2,coo2=Lop.pc_step_orbit_rep(bo2,coo2,dt,MB,g,mom,drop,pr_o)
