"""Figure for the orbit-canonicalizer fix: measured `TranslationGroup`
canonicalization cost vs group order (left) and the resulting per-entry
`pc_step` cost, momentum vs real space (right).

Data are the measurements recorded in log.md — same machine, same commit,
back-to-back before/after builds (M4 MacBook Air, release, min-of-3).
"""

import matplotlib.pyplot as plt
import numpy as np

L = np.array([4, 8, 16, 32, 64, 96])
BEFORE = np.array([0.52, 0.99, 4.05, 27.03, 198.90, 657.40])   # us/call
AFTER = np.array([0.23, 0.34, 0.84, 2.74, 9.10, 19.25])

# pc_step, B = 30_000 live entries, dt = 0.1, us per basis entry
STEP_L = [11, 16]
STEP = {                     # (momentum before, momentum after, real space)
    11: (23.11, 8.08, 1.38),
    16: (61.37, 12.99, 1.49),
}

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

a = ax[0]
a.loglog(L, BEFORE, "o-", color="#c1272d", label="before (rebuild per element)")
a.loglog(L, AFTER, "o-", color="#0b6e4f", label="after (odometer walk)")
ref = L.astype(float) ** 3
a.loglog(L, ref / ref[2] * AFTER[2] * 8, "--", color="#c1272d", alpha=0.4, lw=1,
         label=r"$\propto |G|^2 N$")
ref2 = L.astype(float) ** 2
a.loglog(L, ref2 / ref2[2] * AFTER[2], "--", color="#0b6e4f", alpha=0.4, lw=1,
         label=r"$\propto |G|\, N$")
for x, b, f in zip(L, BEFORE, AFTER):
    a.annotate(f"{b / f:.0f}×", (x, np.sqrt(b * f)), ha="center", va="center",
               fontsize=8, color="0.35")
a.set_xlabel("ladder length $L$   ($|G| = L$, $N = 2L$)")
a.set_ylabel(r"$\mu$s per canonicalization")
a.set_title("(a) orbit-rep canonicalization")
a.legend(fontsize=8, frameon=False)
a.grid(alpha=0.25, which="both", lw=0.4)

a = ax[1]
w, xs = 0.26, np.arange(len(STEP_L))
mb = [STEP[x][0] for x in STEP_L]
ma = [STEP[x][1] for x in STEP_L]
rs = [STEP[x][2] for x in STEP_L]
a.bar(xs - w, mb, w, color="#c1272d", label="momentum $k{=}1$, before")
a.bar(xs, ma, w, color="#0b6e4f", label="momentum $k{=}1$, after")
a.bar(xs + w, rs, w, color="#3b6ea5", label="real space")
for x, v in zip(xs - w, mb):
    a.text(x, v + 1.5, f"{v:.0f}", ha="center", fontsize=8)
for x, v in zip(xs, ma):
    a.text(x, v + 1.5, f"{v:.1f}", ha="center", fontsize=8)
for x, v in zip(xs + w, rs):
    a.text(x, v + 1.5, f"{v:.1f}", ha="center", fontsize=8)
a.set_xticks(xs)
a.set_xticklabels([f"$L={x}$" for x in STEP_L])
a.set_ylabel(r"$\mu$s per basis entry, one pc_step")
a.set_title("(b) pc_step at $B = 30{,}000$")
a.legend(fontsize=8, frameon=False)
a.grid(alpha=0.25, axis="y", lw=0.4)

fig.tight_layout()
fig.savefig(__file__.replace("plot_canon.py", "canonicalize_fix.png"), dpi=160)
print("wrote canonicalize_fix.png")
