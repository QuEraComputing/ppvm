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

L, T = 7, 2.0
npz = np.load(os.path.join(HERE, "exact_msd_L7_T2.npz"))
tR, MSDR = npz["ts"], npz["msd"]          # t = 0, .2, ..., 2.0

def run_cell(mode, dt, knob, kleak=None):
    steps = round(T / dt)
    out = f"data/msd_{mode}_dt{dt}_knob{knob:g}" + (f"_K{kleak:g}" if kleak is not None else "") + ".h5"
    cmd = ["/usr/bin/time", "-l", "./run", "run", "python", "main_realspace_ladder.py",
           "--mode", mode, "--L", str(L), "--gamma", "0.0", "--dt", str(dt),
           "--steps", str(steps), "--pbc", "1", "--preserve", "1", "--out", out]
    cmd += ["--min_abs_coeff", str(knob)] if mode == "trotter" else ["--drop_tol", str(knob)]
    env = {**os.environ}
    env.pop("PPVM_K_LEAKAGE", None)
    if kleak is not None:
        env["PPVM_K_LEAKAGE"] = str(kleak)
    r = subprocess.run(cmd, cwd=XY, capture_output=True, text=True, env=env)
    m = re.search(r"(\d+)\s+maximum resident set size", r.stderr)
    rss_mb = int(m.group(1)) / (1024 * 1024) if m else float("nan")
    try:
        t, msd = msd_from_profile(os.path.join(XY, out))
        with h5py.File(os.path.join(XY, out)) as h:
            wall = float(h["wall_s"][-1])
            peak = int(h["n_basis"][:].max())
        stride = round(0.2 / dt)
        msd_s = msd[::stride]                 # -> t = 0, .2, ..., 2.0
        rel = np.abs(msd_s[1:] - MSDR[1:]) / np.abs(MSDR[1:])   # exclude t=0
        return dict(med=float(np.median(rel)), wall=wall, rss=rss_mb, peak=peak)
    except Exception as e:
        return dict(fail=f"{e} / {r.stdout[-150:]}{r.stderr[-150:]}")

print(f"REAL-SPACE MSD SCAN: ladder L={L}(N={2*L}) gamma=0 T={T} (exact-ED ref; median rel over t=0.2..2)", flush=True)
print(f"{'mode':9} {'dt':>6} {'knob':>7} {'K':>4} {'median_rel':>10} {'wall_s':>7} {'peakRSS_mb':>10} {'peak':>9}", flush=True)
for row in sys.argv[1:]:
    f = row.split(":")
    mode, dt, knob = f[0], float(f[1]), float(f[2])
    kleak = next((float(t[1:]) for t in f[3:] if t.startswith("K")), None)
    d = run_cell(mode, dt, knob, kleak)
    ktag = "-" if kleak is None else f"{kleak:g}"
    if "fail" in d:
        print(f"{mode:9} {dt:>6} {knob:>7.0e} {ktag:>4} FAILED {d['fail'][:120]}", flush=True)
        continue
    print(f"{mode:9} {dt:>6} {knob:>7.0e} {ktag:>4} {d['med']:>10.2e} {d['wall']:>7.1f} {d['rss']:>10.0f} {d['peak']:>9}", flush=True)
