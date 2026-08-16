"""Figures for the orbit-canonicalizer work: canonicalization cost vs group
order (a), per-entry `pc_step` cost (b), and the momentum-vs-real-space
per-entry ratio across sizes (c).

Three implementations:
  original  — rebuild each group element from the identity, O(|G|^2 N)
  odometer  — incremental group walk,                       O(|G| N)
  booth     — staged least-rotation (Booth/Duval) scan,      O(N)

Data are the measurements recorded in log.md — same machine, same commit,
back-to-back builds (M4 MacBook Air, release, min-of-3).
"""

import matplotlib.pyplot as plt
import numpy as np

L = np.array([4, 8, 16, 32, 64, 96])
ORIGINAL = np.array([0.52, 0.99, 4.05, 27.03, 198.90, 657.40])   # us/call
ODOMETER = np.array([0.23, 0.34, 0.84, 2.74, 9.10, 19.25])
BOOTH = np.array([0.21, 0.24, 0.32, 0.47, 0.84, 1.20])

# one pc_step, B = 30_000 live entries, dt = 0.1, us per basis entry
STEP_L = [11, 16]
MOM = {"original": [23.11, 61.37], "odometer": [8.08, 12.99], "booth": [3.77, 4.27]}
REAL = [1.42, 1.42]

# momentum/real per-entry ratio vs size (booth), plus the two original points
RATIO_L = np.array([8, 11, 16, 24, 32])
RATIO_BOOTH = np.array([2.58, 2.62, 3.01, 3.45, 3.96])
RATIO_ORIG_L = np.array([11, 16])
RATIO_ORIG = np.array([15.64, 35.13])

C = {"original": "#c1272d", "odometer": "#e08214", "booth": "#0b6e4f", "real": "#3b6ea5"}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

a = ax[0]
a.loglog(L, ORIGINAL, "o-", color=C["original"], label=r"original  $O(|G|^2 N)$")
a.loglog(L, ODOMETER, "o-", color=C["odometer"], label=r"odometer  $O(|G| N)$")
a.loglog(L, BOOTH, "o-", color=C["booth"], label=r"Booth/Duval  $O(N)$")
for x, b, f in zip(L, ORIGINAL, BOOTH):
    a.annotate(f"{b / f:.0f}×", (x, np.sqrt(b * f)), ha="center", va="center",
               fontsize=8, color="0.35")
a.set_xlabel("ladder length $L$   ($|G| = L$, $N = 2L$)")
a.set_ylabel(r"$\mu$s per canonicalization")
a.set_title("(a) orbit-rep canonicalization")
a.legend(fontsize=8, frameon=False)
a.grid(alpha=0.25, which="both", lw=0.4)

a = ax[1]
w, xs = 0.2, np.arange(len(STEP_L))
for i, key in enumerate(["original", "odometer", "booth"]):
    a.bar(xs + (i - 1.5) * w, MOM[key], w, color=C[key], label=f"momentum $k{{=}}1$, {key}")
    for x, v in zip(xs + (i - 1.5) * w, MOM[key]):
        a.text(x, v + 1.2, f"{v:.1f}", ha="center", fontsize=7)
a.bar(xs + 1.5 * w, REAL, w, color=C["real"], label="real space")
for x, v in zip(xs + 1.5 * w, REAL):
    a.text(x, v + 1.2, f"{v:.1f}", ha="center", fontsize=7)
a.set_xticks(xs)
a.set_xticklabels([f"$L={x}$" for x in STEP_L])
a.set_ylabel(r"$\mu$s per basis entry, one pc_step")
a.set_title("(b) pc_step at $B = 30{,}000$")
a.legend(fontsize=8, frameon=False)
a.grid(alpha=0.25, axis="y", lw=0.4)

a = ax[2]
a.plot(RATIO_L, RATIO_BOOTH, "o-", color=C["booth"], label="after (Booth/Duval)")
a.plot(RATIO_ORIG_L, RATIO_ORIG, "o-", color=C["original"], label="original")
a.axhline(1.0, color="0.5", lw=0.8, ls=":")
a.text(30, 1.15, "per-entry parity", fontsize=8, color="0.4", ha="right")
# a single k-mode holds |G|x fewer entries, so the break-even vs a full
# real-space run sits at ratio = |G|
a.plot(RATIO_L, RATIO_L, "--", color="0.6", lw=1,
       label=r"break-even vs full real-space run ($=|G|$)")
a.set_yscale("log")
a.set_xlabel("ladder length $L$")
a.set_ylabel("momentum / real space, per entry")
a.set_title("(c) remaining per-entry gap")
a.legend(fontsize=8, frameon=False, loc="upper left")
a.grid(alpha=0.25, which="both", lw=0.4)

fig.tight_layout()
fig.savefig(__file__.replace("plot_canon.py", "canonicalize_fix.png"), dpi=160)
print("wrote canonicalize_fix.png")
