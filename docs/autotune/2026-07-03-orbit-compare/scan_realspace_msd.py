"""Real-space MSD comparison at L=7, gamma=0, T=2: Trotter vs adaptive (pec),
driving xy-experiments/main_realspace_ladder.py.  Metric: median over
t = 0.2..2.0 (t=0 excluded, MSD=0) of |MSD(t)-MSD_exact(t)|/|MSD_exact(t)|.
Reference: exact_msd_L7_T2.npz (blocked |U|^2 ED, validated at L=5).
Peak RSS via /usr/bin/time -l; wall from the run's own wall_s dataset.
K (PPVM_K_LEAKAGE) applies to the adaptive mode only.

Usage: scan_realspace_msd.py [row ...]
       rows like  adaptive:0.05:1e-3[:K1]   trotter:0.05:3e-4
       (knob = drop_tol for adaptive, min_abs_coeff for trotter)
"""
import os, re, subprocess, sys
import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
XY = "/Users/alexschuckert/dev/26_ppvm/xy-experiments"
sys.path.insert(0, XY)
from msd import msd_from_profile

L = int(os.environ.get("SCAN_L", "7"))
T = float(os.environ.get("SCAN_T", "2.0"))
try:
    npz = np.load(os.path.join(HERE, f"exact_msd_L{L}_T{T:g}.npz"))
    tR, MSDR = npz["ts"], npz["msd"]      # t = 0 .. T
    REF_DT = float(tR[1] - tR[0])
except FileNotFoundError:
    tR = MSDR = None                       # no exact reference (large L)
    REF_DT = None

def run_cell(mode, dt, knob, kleak=None, stream=False, mb=None, ab=None):
    steps = round(T / dt)
    out = (f"data/msd{'' if L==7 else f'_L{L}'}{'' if T==2.0 else f'_T{T:g}'}_{mode}_dt{dt}_knob{knob:g}" + (f"_K{kleak:g}" if kleak is not None else "")
           + (f"_M{mb:g}" if mb is not None else "") + (f"_A{ab:g}" if ab is not None else "")
           + ("_stream" if stream else "") + ".h5")
    cmd = ["/usr/bin/time", "-l", "./run", "run", "python", "main_realspace_ladder.py",
           "--mode", mode, "--L", str(L), "--gamma", "0.0", "--dt", str(dt),
           "--steps", str(steps), "--pbc", "1", "--preserve", "1", "--out", out]
    cmd += ["--min_abs_coeff", str(knob)] if mode == "trotter" else ["--drop_tol", str(knob)]
    if mb is not None:
        cmd += ["--max_basis", str(int(mb))]
    if ab is not None:
        cmd += ["--admit_basis", str(int(ab))]
    if kleak is not None and mode != "trotter":
        cmd += ["--tau_add", str(kleak * knob / dt)]   # K token -> explicit tau_add
    r = subprocess.run(cmd, cwd=XY, capture_output=True, text=True)
    m = re.search(r"(\d+)\s+maximum resident set size", r.stderr)
    rss_mb = int(m.group(1)) / (1024 * 1024) if m else float("nan")
    try:
        t, msd = msd_from_profile(os.path.join(XY, out))
        with h5py.File(os.path.join(XY, out)) as h:
            wall = float(h["wall_s"][-1])
            peak = int(h["n_basis"][:].max())
        if MSDR is None:
            return dict(med=float("nan"), mx=float("nan"), wall=wall, rss=rss_mb, peak=peak)
        if dt <= REF_DT:                      # run finer than ref: subsample run
            stride = max(1, round(REF_DT / dt))
            msd_s = msd[::stride][:len(MSDR)]
            ref = MSDR[:len(msd_s)]
        else:                                  # run coarser than ref: subsample ref
            rstride = max(1, round(dt / REF_DT))
            ref = MSDR[::rstride][:len(msd)]
            msd_s = msd[:len(ref)]
        rel = np.abs(msd_s[1:] - ref[1:]) / np.abs(ref[1:])
        return dict(med=float(np.median(rel)), mx=float(rel.max()), wall=wall, rss=rss_mb, peak=peak)
    except Exception as e:
        return dict(fail=f"{e} / {r.stdout[-150:]}{r.stderr[-150:]}")

print(f"REAL-SPACE MSD SCAN: ladder L={L}(N={2*L}) gamma=0 T={T} (exact-ED ref; median rel over t=0.2..2)", flush=True)
print(f"{'mode':9} {'dt':>6} {'knob':>7} {'K':>4} {'median_rel':>10} {'max_rel':>10} {'wall_s':>7} {'peakRSS_mb':>10} {'peak':>9}", flush=True)
for row in sys.argv[1:]:
    f = row.split(":")
    mode, dt, knob = f[0], float(f[1]), float(f[2])
    kleak = next((float(t[1:]) for t in f[3:] if t.startswith("K")), None)
    mb = next((float(t[1:]) for t in f[3:] if t.startswith("M")), None)
    ab = next((float(t[1:]) for t in f[3:] if t.startswith("A")), None)
    stream = "stream" in f[3:]
    d = run_cell(mode, dt, knob, kleak, stream, mb, ab)
    ktag = ("-" if kleak is None else f"{kleak:g}") + ("s" if stream else "") + ("" if mb is None else f"/M{mb:g}") + ("" if ab is None else f"/A{ab:g}")
    if "fail" in d:
        print(f"{mode:9} {dt:>6} {knob:>7.0e} {ktag:>4} FAILED {d['fail'][:120]}", flush=True)
        continue
    print(f"{mode:9} {dt:>6} {knob:>7.0e} {ktag:>4} {d['med']:>10.2e} {d.get('mx',float('nan')):>10.2e} {d['wall']:>7.1f} {d['rss']:>10.0f} {d['peak']:>9}", flush=True)
