# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Plot TFIM Trotter runtime vs qubit count, log-y, for the ppvm hashers and
PauliPropagation.jl.

Inputs are the CSVs produced by
``crates/ppvm-runtime/examples/trotter_qubit_sweep.rs`` and
``julia-benchmarks/benches/trotter_sweep.jl`` (columns:
``qubits,hasher,bytes,time_s,terms``).

    uv run --with matplotlib python benchmarks/plot_tfim_sweep.py \
        --ppvm /tmp/tfim_sweep/ppvm.csv \
        --pp   /tmp/tfim_sweep/pp.csv \
        --out  /tmp/tfim_sweep/tfim_trotter_scaling.png
"""

import argparse
import csv
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# series key -> (label, color, marker, linestyle)
SERIES = {
    "fxhash_nofold": ("ppvm — fxhash, no fold (pre-PR)", "#c0392b", "o", "-"),
    "fxhash": ("ppvm — fxhash + fold (this PR)", "#5b3fb8", "s", "-"),
    "gxhash": ("ppvm — gxhash", "#1f9e6e", "^", "-"),
    "pauli_propagation_jl": ("PauliPropagation.jl", "#7f8c8d", "D", "--"),
}


def load(path):
    """Return {hasher: ([qubits], [time_s])} sorted by qubits."""
    rows = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f):
            rows[r["hasher"]].append((int(r["qubits"]), float(r["time_s"])))
    return {k: tuple(zip(*sorted(v))) for k, v in rows.items()}


def storage_bytes(n):
    need = -(-n // 8)
    k = 0
    while (1 << k) <= need:
        k += 1
    return 1 << k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ppvm", required=True)
    ap.add_argument("--pp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.ppvm)
    data.update(load(args.pp))

    fig, ax = plt.subplots(figsize=(9, 5.5))

    # Shade storage-tier bands (the bump is tied to these) using ppvm's qubits.
    all_q = sorted({q for s in data.values() for q in s[0]})
    qmin, qmax = min(all_q), max(all_q)
    boundaries = [n for n in range(qmin, qmax + 2) if storage_bytes(n) != storage_bytes(n - 1)]
    band_edges = [qmin] + boundaries + [qmax + 1]
    for lo, hi in zip(band_edges, band_edges[1:]):
        b = storage_bytes(lo)
        ax.axvspan(lo - 0.5, hi - 0.5, color="0.5", alpha=0.05, lw=0)
        ax.text(
            (lo + min(hi, qmax) - 1) / 2,
            0.93,
            f"[u8;{b}]",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=8,
            color="0.45",
        )

    for key, (label, color, marker, ls) in SERIES.items():
        if key not in data:
            continue
        q, t = data[key]
        ax.plot(q, t, ls, color=color, marker=marker, ms=5, lw=1.7, label=label)

    ax.set_yscale("log")
    ax.set_xlabel("number of qubits")
    ax.set_ylabel("runtime per Trotter run (s)")
    ax.set_title(
        "TFIM Trotter scaling: fxhash bucket cliff vs gxhash\n"
        "J=1.0, h=1, dt=0.1, 20 steps, truncation 1e-6, depolarizing 1e-4",
        fontsize=11,
    )
    ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
