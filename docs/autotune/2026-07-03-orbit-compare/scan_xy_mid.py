"""Mid-precision dt x drop_tol (x K for pec) scan at L=7, k=1, T=2, driving the
clean xy-experiments harnesses (main_k_pec_ladder.py / main_k_xy_ladder.py)
instead of orbit_bench_L7.py.  Differences vs orbit_bench_L7 that matter:
- Trotter evolves the momentum-MERGED pair (momentum_merge each step, ~L x
  fewer terms live) -- the proper symmetry-compressed baseline; orbit_bench
  evolved the full real-space pair and merged only at readout.
- pec: same pc_step_orbit_rep kernel (no protected_arr; validated 1.86e-3
  at the L=5 exact-ED gate, matching the ledger).
Peak RSS via /usr/bin/time -l (their psutil attr is end-RSS, not peak).
rel vs exact ED npz (t = 0, .2, ..., 2.0).

Usage: scan_xy_mid.py [row ...]   rows like  pec:0.1:3e-3[:K1]  trot:0.05:3e-4
"""
import os, re, subprocess, sys, time
import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
XY = "/Users/alexschuckert/dev/26_ppvm/xy-experiments"
L, K_MODE, T = 7, 1, 2.0
REF = np.load(os.path.join(HERE, "exact_ref_L7_k1_T2.npz"))["ref"]

def run_cell(kind, dt, drop, kleak=None):
    steps = round(T / dt)
    stride = round(0.2 / dt)
    out = f"data/scan_{kind}_dt{dt}_drop{drop:g}" + (f"_K{kleak:g}" if kleak else "") + ".h5"
    script = "main_k_pec_ladder.py" if kind == "pec" else "main_k_xy_ladder.py"
    cmd = ["/usr/bin/time", "-l", "./run", "run", "python", script,
           "--L", str(L), "--dt", str(dt), "--steps", str(steps),
           "--drop_tol", str(drop), "--ks", str(K_MODE), "--out", out]
    # PPVM_K_LEAKAGE was removed from ppvm (2026-07-07): K tokens now convert
    # to the explicit --tau_add flag (pec cells only).
    if kleak is not None and kind == "pec":
        cmd += ["--tau_add", str(kleak * drop / dt)]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=XY, capture_output=True, text=True)
    wall = time.time() - t0
    m = re.search(r"(\d+)\s+maximum resident set size", r.stderr)
    rss_mb = int(m.group(1)) / (1024 * 1024) if m else float("nan")
    try:
        with h5py.File(os.path.join(XY, out)) as h:
            C = (h["C_re"][0] + 1j * h["C_im"][0])[::stride]
            peak = int(h["n_basis"][0].max())
    except Exception as e:
        return dict(fail=str(e) + " / " + r.stdout[-200:] + r.stderr[-200:])
    rel = float(np.linalg.norm(C - REF) / np.linalg.norm(REF))
    return dict(rel=rel, wall=wall, rss=rss_mb, peak=peak)

rows = sys.argv[1:]
print(f"XY-EXPERIMENTS MID SCAN: ladder L={L}(N={2*L}) k={K_MODE} T={T} (exact-ED ref)", flush=True)
print(f"{'method':8} {'dt':>7} {'drop':>7} {'K':>4} {'rel_err':>9} {'wall_s':>7} {'peakRSS_mb':>10} {'peak':>8}", flush=True)
for row in rows:
    f = row.split(":")
    kind, dt, drop = f[0], float(f[1]), float(f[2])
    kleak = next((float(t[1:]) for t in f[3:] if t.startswith("K")), None)
    d = run_cell(kind, dt, drop, kleak)
    ktag = "-" if kleak is None else f"{kleak:g}"
    if "fail" in d:
        print(f"{kind:8} {dt:>7} {drop:>7.0e} {ktag:>4} FAILED {d['fail'][:120]}", flush=True)
        continue
    print(f"{kind:8} {dt:>7} {drop:>7.0e} {ktag:>4} {d['rel']:>9.2e} {d['wall']:>7.1f} {d['rss']:>10.0f} {d['peak']:>8}", flush=True)
