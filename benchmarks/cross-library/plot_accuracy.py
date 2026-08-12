# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Render the truncation-error scaling from `xbench_accuracy.py`'s CSVs.

Four panels, because four separate questions were asked of this data and each one
has a different answer:

1. **Error vs `atol`** — every engine falls one decade per decade of `atol`, so
   the two truncation rules share a convergence order and differ only in the
   constant. This panel is why the rules are not qualitatively apart.
2. **Error ratio vs depth** — the loss compounds per gate, so the gap grows with
   circuit depth.
3. **Error ratio vs rotation angle** — the ratio *crosses 1.0*. Large angles
   favour monoprop, because that is where accumulate-then-truncate has its own
   failure mode: `truncate()` after each gate permanently deletes a term whose
   merged coefficient transiently cancelled. Any summary quoting one ratio is
   hiding this panel.
4. **The scrambling instance** — where the rules genuinely separate, beside a
   `θ = π/4` control on the same circuit family that reverses the sign.

Ratios are monoprop / ppvm on the `L2` distance to a converged reference, so
above 1.0 means monoprop is less accurate. The scalar `observable` is
deliberately not plotted: its truncation error changes sign, so ratios built from
it are dominated by whichever engine happened to land near a zero crossing (see
`README.md`).

    uv run --no-project --with matplotlib python3 plot_accuracy.py \
        --csv accuracy.csv --divergence accuracy_divergence.csv \
        --out target/xbench/accuracy.png
"""

from __future__ import annotations

import argparse
import csv
import statistics
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_xbench import SERIES, SURFACE, TEXT_PRIMARY, TEXT_SECONDARY, style

BASELINE = "ppvm"
PAIR = ("ppvm", "monoprop")
# Ratio panels plot a *derived* quantity, not an engine, so they stay off the
# engine palette — those hues are reserved for identity throughout the harness.
RATIO_INK = ["#0b0b0b", "#8a8986"]
RATIO_MARKS = ["v", "o"]


def load(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh))


def relative(rows, library: str) -> dict[float, list[float]]:
    """`{atol: [relative L2 error, ...]}` for one engine, pooled over seeds."""
    out: dict[float, list[float]] = defaultdict(list)
    for r in rows:
        if r["library"] == library:
            out[float(r["atol"])].append(float(r["l2_err"]) / float(r["ref_norm"]))
    return out


def ratio_vs(rows, key: str) -> dict[float, list[float]]:
    """`{x: [monoprop/ppvm error ratio, ...]}` as `key` is swept.

    A ratio is only formed inside one fully-specified cell — same model, width,
    depth, angle, seed and `atol` — so sweeping one axis never divides across
    another. Remaining axes (seeds) pool into the list.
    """
    cells: dict[tuple, dict[str, float]] = defaultdict(dict)
    for r in rows:
        if r["library"] not in PAIR:
            continue
        cell = (r["model"], r["qubits"], r["steps"], r["dt"], r["seed"], r["atol"])
        cells[cell][r["library"]] = float(r["l2_err"])
        cells[cell]["x"] = float(r[key])
    out: dict[float, list[float]] = defaultdict(list)
    for libs in cells.values():
        if len(libs) < 3 or not libs.get(BASELINE):
            continue
        out[libs["x"]].append(libs["monoprop"] / libs[BASELINE])
    return dict(sorted(out.items()))


def band(ax, data, color, marker, label) -> None:
    """Median line, with a min/max ribbon where an axis pooled several seeds."""
    xs = list(data)
    ax.plot(
        xs,
        [statistics.median(data[x]) for x in xs],
        color=color,
        marker=marker,
        markersize=5,
        linewidth=1.6,
        label=label,
    )
    if any(len(v) > 1 for v in data.values()):
        ax.fill_between(
            xs,
            [min(data[x]) for x in xs],
            [max(data[x]) for x in xs],
            color=color,
            alpha=0.16,
            linewidth=0,
        )


def panel_atol(ax, rows) -> None:
    """Absolute error against atol, every engine, both Trotter models."""
    for model, dash in (("tfim", "-"), ("heisenberg", "--")):
        for lib, (label, color, marker) in SERIES.items():
            pts = sorted(
                (float(r["atol"]), float(r["l2_err"]) / float(r["ref_norm"]))
                for r in rows
                if r["model"] == model and r["library"] == lib
            )
            if pts:
                ax.plot(
                    *zip(*pts),
                    dash,
                    color=color,
                    marker=marker,
                    markersize=4,
                    linewidth=1.4,
                    alpha=0.9,
                    label=label if model == "tfim" else None,
                )
    # Slope-1 guide, anchored on the baseline's own worst point so it sits beside
    # the data rather than floating below it -- the eye is comparing gradients.
    anchor = max(
        (float(r["atol"]), float(r["l2_err"]) / float(r["ref_norm"]))
        for r in rows
        if r["model"] == "tfim" and r["library"] == BASELINE
    )
    lo = min(float(r["atol"]) for r in rows)
    ax.plot(
        [lo, anchor[0]],
        [anchor[1] * lo / anchor[0], anchor[1]],
        ":",
        color=TEXT_SECONDARY,
        linewidth=1.3,
        label="slope 1  (error ∝ atol)",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("truncation threshold   atol")
    ax.set_ylabel("relative L2 error   ‖c − c*‖ / ‖c*‖")
    ax.set_title(
        "1 · Every rule converges at the same order\n"
        "solid TFIM, dashed Heisenberg — n=8, 10 steps, dt=0.1",
        fontsize=10,
        color=TEXT_PRIMARY,
        loc="left",
    )
    ax.legend(fontsize=7, frameon=False, labelcolor=TEXT_SECONDARY, loc="lower left")


def panel_ratio(ax, rows, key, xlabel, title) -> None:
    """monoprop / ppvm error ratio against one swept axis, one line per atol.

    The swept values are logarithmically spaced and uneven, so the axis is log
    and ticks are pinned to the values actually run.
    """
    ax.set_xscale("log")
    ax.minorticks_off()
    # Loosest threshold first: that is the series with the largest effect, and
    # `sorted` on the raw strings would put 1e-5 ahead of 1e-3.
    xs: list[float] = []
    for i, atol in enumerate(
        sorted({r["atol"] for r in rows}, key=float, reverse=True)
    ):
        data = ratio_vs([r for r in rows if r["atol"] == atol], key)
        if data:
            band(
                ax,
                data,
                RATIO_INK[i % len(RATIO_INK)],
                RATIO_MARKS[i % len(RATIO_MARKS)],
                f"atol = {atol}",
            )
            xs = list(data)
    # Set after the scale, which would otherwise install its own locator.
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{x:g}" for x in xs], fontsize=8)
    ax.axhline(1.0, color=SERIES["ppvm"][1], linewidth=1.2, linestyle=":")
    ax.annotate(
        "equal accuracy",
        (0.03, 1.0),
        xycoords=("axes fraction", "data"),
        fontsize=7,
        color=SERIES["ppvm"][1],
        va="bottom",
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel("monoprop error / ppvm error")
    ax.set_title(title, fontsize=10, color=TEXT_PRIMARY, loc="left")
    ax.legend(fontsize=7, frameon=False, labelcolor=TEXT_SECONDARY)


def panel_scramble(ax, rows) -> None:
    """The scrambling instance, small angle against the theta = pi/4 control."""
    for dt, dash, tag in (
        ("0.05", "-", "dt=0.05, 3200 gates"),
        ("0.7853981633974483", "--", "dt=π/4, 80 gates"),
    ):
        sub = [r for r in rows if r["dt"] == dt]
        for lib in PAIR:
            data = relative(sub, lib)
            if data:
                label, color, marker = SERIES[lib]
                ax.plot(
                    sorted(data),
                    [statistics.median(data[x]) for x in sorted(data)],
                    dash,
                    color=color,
                    marker=marker,
                    markersize=5,
                    linewidth=1.6,
                    label=f"{label} — {tag}",
                )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("truncation threshold   atol")
    ax.set_ylabel("relative L2 error")
    ax.set_title(
        "4 · Scrambling instance: where they separate\n"
        "n=8, all 65535 Pauli words populated, median of 5 seeds",
        fontsize=10,
        color=TEXT_PRIMARY,
        loc="left",
    )
    ax.legend(fontsize=7, frameon=False, labelcolor=TEXT_SECONDARY, loc="lower left")


CAPTION = (
    "Error is the L2 distance of the whole propagated coefficient vector from a "
    "converged reference (ppvm at atol=1e-16), relative to ‖c*‖. ppvm, "
    "PauliPropagation.jl and PauliStrings.jl accumulate every contribution to a term "
    "and then threshold the merged coefficient; pauli-prop and monoprop threshold each "
    "branch as it is emitted. A Pauli rotation sends Q to exactly two words, so a "
    "merged value can exceed the larger contribution by at most 2x — which bounds the "
    "divergence and ties both errors to the same atol. Panels 2-4: monoprop needs small "
    "sin θ (a wide rejection band), depth (to populate and compound it) and a scrambled "
    "operator (so the loss does not average out) simultaneously before it separates; "
    "drop any one and the gap collapses or reverses. Ribbons are min/max over seeds."
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, required=True, help="the 5-engine atol sweep")
    ap.add_argument(
        "--divergence",
        type=Path,
        required=True,
        help="the depth / angle / scramble sweeps",
    )
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    div = load(args.divergence)
    # The two Heisenberg sweeps are told apart by width, which is how they were
    # run: depth at n=6 (cheap enough for 320 steps), angle at n=8.
    depth = [r for r in div if r["model"] == "heisenberg" and r["qubits"] == "6"]
    angle = [r for r in div if r["model"] == "heisenberg" and r["qubits"] == "8"]
    scramble = [r for r in div if r["model"] == "scramble"]
    for name, rows in (("depth", depth), ("angle", angle), ("scramble", scramble)):
        if not rows:
            raise SystemExit(f"no {name} rows in {args.divergence}")

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.6), facecolor=SURFACE)
    for ax in axes.flat:
        style(ax)

    panel_atol(axes[0][0], load(args.csv))
    panel_ratio(
        axes[0][1],
        depth,
        "steps",
        "Trotter steps",
        "2 · The gap compounds with depth\nHeisenberg n=6, dt=0.05 (sin θ ≈ 0.1)",
    )
    panel_ratio(
        axes[1][0],
        angle,
        "dt",
        "dt      (θ = 2·dt, so dt = π/8 gives sin θ = 0.707)",
        "3 · …and reverses at large angles\nHeisenberg n=8, 10 steps",
    )
    panel_scramble(axes[1][1], scramble)

    fig.text(
        0.008,
        0.012,
        textwrap.fill(CAPTION, 168),
        fontsize=7.5,
        color=TEXT_SECONDARY,
        va="bottom",
    )
    fig.tight_layout(rect=(0, 0.072, 1, 1))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200, facecolor=SURFACE)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
