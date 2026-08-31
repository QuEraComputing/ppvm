# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Render the cross-library benchmark from `run_xbench.py`'s `results.csv`.

One row per model, three panels: runtime vs qubit count (log-y), runtime
relative to `ppvm`, and the size of the propagated operator — which is the
workload, and the thing that makes the runtimes comparable in the first place.

    uv run --no-project --with matplotlib python3 plot_xbench.py \
        --csv target/xbench/results.csv --out target/xbench/xbench.png
"""

from __future__ import annotations

import argparse
import csv
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Slots from the validated default palette, in this fixed order — it clears
# every hard gate on the adjacent pairlist that line charts use (worst CVD
# ΔE 9.2, worst normal-vision ΔE 27.5). Colour follows the entity, so a run with
# a subset of `--libs` keeps every survivor's hue. Magenta and aqua sit below
# 3:1 on the light surface, which obligates relief — `results.csv`, the driver's
# summary table, and the per-series markers are it.
BASELINE = "ppvm"
SERIES = {
    "ppvm": ("ppvm (this repo)", "#2a78d6", "o"),
    "pauli-propagation-jl": ("PauliPropagation.jl", "#eb6834", "s"),
    "pauli-strings-jl": ("PauliStrings.jl", "#1baf7a", "^"),
    "pauli-prop": ("pauli-prop (Qiskit)", "#4a3aa7", "D"),
    "monoprop": ("monoprop", "#e87ba4", "v"),
}
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#d8d7d3"
SURFACE = "#fcfcfb"
MODEL_TITLES = {
    "tfim": "TFIM Trotter — ⟨0|Σ Z_i(t)|0⟩",
    "heisenberg": "Heisenberg correlations — tr[Z₀ Z₀(t)]/2ⁿ",
}


def load(path: Path):
    """`{model: {library: {n: row}}}`."""
    out: dict[str, dict[str, dict[int, dict[str, str]]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    with path.open() as fh:
        for row in csv.DictReader(fh):
            out[row["model"]][row["library"]][int(row["qubits"])] = row
    return out


def style(ax) -> None:
    """Recessive grid and axes; the data carries the ink."""
    ax.set_facecolor(SURFACE)
    ax.grid(True, which="major", color=GRID, linewidth=0.8, alpha=0.9)
    ax.grid(True, which="minor", color=GRID, linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=9)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    data = load(args.csv)
    models = [m for m in ("tfim", "heisenberg") if m in data]
    fig, axes = plt.subplots(
        len(models),
        3,
        figsize=(15.5, 4.4 * len(models)),
        squeeze=False,
        facecolor=SURFACE,
    )

    for r, model in enumerate(models):
        per_lib = data[model]
        ax_t, ax_s, ax_n = axes[r]
        widths = sorted({n for lib in per_lib.values() for n in lib})
        term_series: dict[str, tuple[int, ...]] = {}
        for key, (label, color, marker) in SERIES.items():
            if key not in per_lib:
                continue
            ns = sorted(per_lib[key])
            times = [float(per_lib[key][n]["time_s"]) for n in ns]
            ax_t.plot(
                ns,
                times,
                color=color,
                marker=marker,
                markersize=6,
                linewidth=2,
                label=label,
            )
            terms = [int(per_lib[key][n]["terms"]) for n in ns]
            ax_n.plot(
                ns,
                terms,
                color=color,
                marker=marker,
                markersize=6,
                linewidth=2,
                label=label,
            )
            term_series[key] = tuple(terms)
            if key != BASELINE and BASELINE in per_lib:
                shared = [n for n in ns if n in per_lib[BASELINE]]
                ratios = [
                    float(per_lib[key][n]["time_s"])
                    / float(per_lib[BASELINE][n]["time_s"])
                    for n in shared
                ]
                ax_s.plot(
                    shared,
                    ratios,
                    color=color,
                    marker=marker,
                    markersize=6,
                    linewidth=2,
                    label=label,
                )

        # Engines whose support is identical draw the same line, so the ones
        # underneath are invisible. Say so rather than let the panel imply that
        # only the top series was measured.
        if BASELINE in term_series:
            same = [k for k, v in term_series.items() if v == term_series[BASELINE]]
            if len(same) > 1:
                ax_n.annotate(
                    "identical (exactly):\n"
                    + "\n".join(SERIES[k][0].split(" (")[0] for k in same),
                    # The support rises left-to-right, so the low-right corner
                    # is the one that stays clear of the marks.
                    xy=(0.97, 0.04),
                    xycoords="axes fraction",
                    fontsize=8,
                    color=TEXT_SECONDARY,
                    ha="right",
                    va="bottom",
                )

        ax_s.axhline(1.0, color=SERIES[BASELINE][1], linewidth=2, linestyle=(0, (4, 3)))
        ax_s.annotate(
            f"{SERIES[BASELINE][0]} = 1",
            xy=(0.02, 1.0),
            xycoords=("axes fraction", "data"),
            va="bottom",
            fontsize=8,
            color=TEXT_SECONDARY,
        )

        for ax, ylabel, title in (
            (ax_t, "runtime (s, min of repeats)", "runtime"),
            (ax_s, f"× {BASELINE}", "relative runtime"),
            (ax_n, "terms in the propagated operator", "workload size"),
        ):
            style(ax)
            ax.set_yscale("log")
            # Qubit counts are integers; let matplotlib interpolate ticks and it
            # invents 6.25 qubits.
            ax.set_xticks(widths)
            ax.set_xticklabels([str(n) for n in widths])
            ax.set_xlabel("qubits", color=TEXT_SECONDARY, fontsize=9)
            ax.set_ylabel(ylabel, color=TEXT_SECONDARY, fontsize=9)
            ax.set_title(title, color=TEXT_PRIMARY, fontsize=10, loc="left")

        ax_t.set_ylabel("runtime (s, min of repeats)", color=TEXT_SECONDARY, fontsize=9)
        ax_t.annotate(
            MODEL_TITLES.get(model, model),
            xy=(0, 1.16),
            xycoords="axes fraction",
            fontsize=12,
            fontweight="bold",
            color=TEXT_PRIMARY,
        )

    # One figure-level legend above the panels: repeating it per row wastes
    # space, and inside the axes it lands on the fastest series.
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=False,
        fontsize=9.5,
        labelcolor=TEXT_SECONDARY,
        ncol=len(labels),
        loc="upper left",
        bbox_to_anchor=(0.006, 0.995),
    )

    sample = next(iter(next(iter(data.values())).values()))
    row = next(iter(sample.values()))
    n_engines = len({lib for per_lib in data.values() for lib in per_lib})
    caption = (
        f"first-order Trotter, steps={row['steps']}, dt={row['dt']}, "
        f"truncation |c| < {row['atol']}; single-threaded; min of repeats. "
        f"All {n_engines} engines are validated to propagate the identical operator "
        "term-for-term before timing, but pauli-prop and monoprop prune at "
        "branch-creation time rather than after accumulation, so at this "
        "truncation they carry a different support — compare the workload panel "
        "before reading their ratios. monoprop is capped to one thread "
        "(monoprop_NUM_THREADS=1); uncapped it takes one partition per core. "
        "Per-point numbers in results.csv."
    )
    if args.title:
        fig.suptitle(
            args.title, fontsize=14, color=TEXT_PRIMARY, x=0.006, y=0.995, ha="left"
        )
        # Legend sits just under the title when there is one.
        fig.legends[0].set_bbox_to_anchor((0.006, 0.962))
    wrapped = textwrap.fill(caption, width=178)
    fig.text(
        0.008,
        0.004,
        wrapped,
        fontsize=8.5,
        color=TEXT_SECONDARY,
        ha="left",
        va="bottom",
    )
    top = 0.90 if args.title else 0.945
    fig.tight_layout(rect=(0, 0.018 * (wrapped.count("\n") + 1) + 0.012, 1, top))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=160, facecolor=SURFACE)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
