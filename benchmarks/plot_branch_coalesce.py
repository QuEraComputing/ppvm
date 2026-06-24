# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Plot the branch-coalesce scaling study: sort-merge (PR #154) vs the
pre-#154 FxHashMap coalesce, as a function of the branch count ``m``.

Reads the JSON criterion writes for the ``branch-coalesce-scaling`` bench
(``crates/ppvm-tableau/benches/branch-coalesce-scaling.rs``) — no CSV step.
Run the bench first, then this script:

    cargo bench -p ppvm-tableau --bench branch-coalesce-scaling
    uv run --with matplotlib python benchmarks/plot_branch_coalesce.py \
        --out /tmp/branch_coalesce_scaling.png

The left panel is the raw time-vs-``m`` scaling (log-log); the right panel is
the sort-merge speedup ``t_hashmap / t_sortmerge`` (>1 → sort-merge wins,
<1 → hash wins), with the crossover line and the "hash wins" band shaded. The
two regimes (doubling = fresh branching, merge = collision-heavy) tell opposite
stories, which is the whole point.
"""

import argparse
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CRITERION_DIR = "target/criterion"

# regime -> (label, color, marker)
REGIMES = {
    "doubling": ("doubling (fresh bit · output 2m)", "#5b3fb8", "o"),
    "merge": ("merge (closed set · output m)", "#d08700", "s"),
}
# algo -> (label, color, linestyle)
ALGOS = {
    "hashmap": ("FxHashMap coalesce (pre-#154)", "#c0392b", "--"),
    "sortmerge": ("sort-merge (this PR)", "#5b3fb8", "-"),
}
# packed-path cutoff in the sort-merge: m > 65535 drops to the generic fallback.
PACKED_CUTOFF = 65535


def load_series(regime, algo):
    """Return ([m...], [seconds...]) sorted by m, from criterion JSON."""
    base = os.path.join(CRITERION_DIR, f"branch-coalesce-{regime}", algo)
    pts = []
    for est in glob.glob(os.path.join(base, "*", "new", "estimates.json")):
        m = int(os.path.basename(os.path.dirname(os.path.dirname(est))))
        with open(est) as f:
            ns = json.load(f)["median"]["point_estimate"]
        pts.append((m, ns * 1e-9))
    if not pts:
        raise SystemExit(
            f"no criterion data under {base} — run the bench first:\n"
            "  cargo bench -p ppvm-tableau --bench branch-coalesce-scaling"
        )
    pts.sort()
    return [m for m, _ in pts], [s for _, s in pts]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fig, (ax_t, ax_s) = plt.subplots(1, 2, figsize=(13, 5.2))

    # ---- left: absolute time vs m (log-log) -----------------------------
    for regime, (rlabel, rcolor, marker) in REGIMES.items():
        for algo, (alabel, acolor, ls) in ALGOS.items():
            m, t = load_series(regime, algo)
            ax_t.plot(
                m,
                t,
                ls,
                color=acolor,
                marker=marker,
                ms=5,
                lw=1.7,
                alpha=0.95 if regime == "doubling" else 0.6,
                label=f"{algo} · {regime}",
            )
    ax_t.set_xscale("log", base=2)
    ax_t.set_yscale("log")
    ax_t.set_xlabel("branch count  m  (= 2^k for k T gates)")
    ax_t.set_ylabel("coalesce time per T gate (s)")
    ax_t.set_title("Coalesce cost vs branch count", fontsize=11)
    ax_t.axvline(PACKED_CUTOFF, color="0.5", ls=":", lw=1)
    ax_t.text(
        PACKED_CUTOFF, ax_t.get_ylim()[0], "  packed→generic",
        rotation=90, va="bottom", ha="left", fontsize=7.5, color="0.45",
    )
    ax_t.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)
    ax_t.legend(frameon=False, fontsize=8.5, loc="upper left")

    # ---- right: sort-merge speedup vs m ---------------------------------
    ax_s.axhline(1.0, color="0.3", lw=1)
    ymax_band = 4.5
    ax_s.axhspan(0, 1.0, color="#c0392b", alpha=0.06, lw=0)
    ax_s.text(
        0.98, 0.04, "hash wins", transform=ax_s.transAxes,
        ha="right", va="bottom", fontsize=9, color="#c0392b",
    )
    ax_s.text(
        0.02, 0.96, "sort-merge wins", transform=ax_s.transAxes,
        ha="left", va="top", fontsize=9, color="#5b3fb8",
    )
    for regime, (rlabel, rcolor, marker) in REGIMES.items():
        m, t_h = load_series(regime, "hashmap")
        _, t_s = load_series(regime, "sortmerge")
        speedup = [h / s for h, s in zip(t_h, t_s)]
        ax_s.plot(m, speedup, "-", color=rcolor, marker=marker, ms=5, lw=1.8, label=rlabel)
    ax_s.set_xscale("log", base=2)
    ax_s.set_ylim(0, ymax_band)
    ax_s.axvline(PACKED_CUTOFF, color="0.5", ls=":", lw=1)
    ax_s.set_xlabel("branch count  m")
    ax_s.set_ylabel("sort-merge speedup  (t_hash / t_sortmerge)")
    ax_s.set_title("Where each coalesce wins", fontsize=11)
    ax_s.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)
    ax_s.legend(frameon=False, fontsize=9, loc="upper center")

    fig.suptitle(
        "ppvm-tableau branch coalesce: sort-merge vs FxHashMap  (80 qubits, u128 index)",
        fontsize=12.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
