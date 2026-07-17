import numpy as np
from functools import reduce
from ppvm import PauliSum
from ppvm._core import TranslationGroup

N=8; k=1; dt=0.02; nsnap=5; steps_per=3; mac=1e-10
phi=2*np.pi*k/N
g=TranslationGroup.chain_1d(N); mom=[k]
bonds=[(j,(j+1)%N) for j in range(N)]

# exact ED (same as before)
SX=np.array([[0,1],[1,0]],complex);SY=np.array([[0,-1j],[1j,0]],complex);I2=np.eye(2)
def op(s,q):return reduce(np.kron,[s if k2==q else I2 for k2 in range(N)])
H=np.zeros((2**N,2**N),complex)
for j,nx in bonds: H+=op(SX,j)@op(SX,nx)+op(SY,j)@op(SY,nx)
idx=np.arange(2**N); d=np.zeros(2**N,complex)
for j in range(N): d+=np.exp(-1j*phi*j)*(1-2*((idx>>(N-1-j))&1))
E,V=np.linalg.eigh(H)
def exact_Ck(t):
    U=(V*np.exp(1j*E*t))@V.conj().T; return (np.conj(d)@((np.abs(U)**2)@d))/(np.conj(d)@d)

# real-pair Trotter, merged each step (orbit-preserving)
def mk():
    R=PauliSum.new(N,[(f"Z{j}", float(np.cos(phi*j))) for j in range(N)],min_abs_coeff=mac,max_pauli_weight=N)
    Im=PauliSum.new(N,[(f"Z{j}", float(-np.sin(phi*j))) for j in range(N)],min_abs_coeff=mac,max_pauli_weight=N)
    return R,Im
R,Im=mk(); Rs,Ims=mk()   # Rs,Ims = seed copies for overlap
Rs.momentum_merge(Ims,g,mom)
norm=Rs.overlap(Rs)+Ims.overlap(Ims)
def Ck(R,Im):
    Rm=R.copy() if hasattr(R,'copy') else R
    # overlap <O_k_seed, O(t)> = (Rs.R + Ims.Im) + i(Rs.Im - Ims.R)
    re=Rs.overlap(R)+Ims.overlap(Im); im=Rs.overlap(Im)-Ims.overlap(R)
    return (re+1j*im)/norm

print(f"N={N} k={k} dt={dt} (merged Trotter, near-exact mac={mac})")
print(f"{'t':>6} {'exact':>20} {'trotter-merged':>22}")
for si in range(nsnap+1):
    t=si*steps_per*dt
    Rc=R.copy(); Imc=Im.copy()
    Rc.momentum_merge(Imc,g,mom)
    ex=exact_Ck(t); tr=Ck(Rc,Imc)
    print(f"{t:6.3f} {ex.real:+.4f}{ex.imag:+.4f}j  {tr.real:+.4f}{tr.imag:+.4f}j")
    if si<nsnap:
        for _ in range(steps_per):
            for a,b in bonds: R.rxx(a,b,2*dt,truncate=False); R.ryy(a,b,2*dt,truncate=False)
            for a,b in bonds: Im.rxx(a,b,2*dt,truncate=False); Im.ryy(a,b,2*dt,truncate=False)
            R.truncate(); Im.truncate()
