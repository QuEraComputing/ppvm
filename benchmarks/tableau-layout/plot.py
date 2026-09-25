# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Plot recorded summaries; no benchmarks run here."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Directory containing summary.csv")
    parser.add_argument("--word-size-n", type=int, default=1024)
    args = parser.parse_args()
    with (args.directory / "summary.csv").open() as file:
        rows = list(csv.DictReader(file))
    metadata = args.directory / "metadata.json"
    parameters = (
        json.loads(metadata.read_text())["arguments"] if metadata.exists() else {}
    )
    caption = (
        f"Median of {parameters['launches']} launches; {parameters['samples']} samples per launch"
        if parameters
        else "Median across process launches"
    )
    series = [
        ("ppvm-native-kernel", "column", "ppvm column kernel", "#171717", "-"),
        ("ppvm-public", "column", "ppvm public API", "#171717", "--"),
        ("ppvm-native-kernel", "row", "ppvm row kernel", "#777777", "-"),
        (
            "ppvm-native-kernel",
            "column-transpose-batch",
            "ppvm batch + transpose",
            "#777777",
            ":",
        ),
        ("quantumclifford", "fastrow", "QC fastrow", "#2468a0", "-"),
        ("quantumclifford", "fastcolumn", "QC fastcolumn", "#2468a0", "--"),
        ("candidate", "qubit_bits", "Rust qubit bits", "#b64d20", "-"),
        ("candidate", "generator_bits", "Rust generator bits", "#b64d20", "--"),
    ]
    operations = [
        ("h", "Hadamard"),
        ("cnot", "CNOT"),
        ("comm", "Row commutation"),
        ("mul", "Signed row multiplication"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for ax, (operation, title) in zip(axes.flat, operations):
        for implementation, layout, label, color, style in series:
            selected = sorted(
                (
                    r
                    for r in rows
                    if r["implementation"] == implementation
                    and r["layout"] == layout
                    and r["operation"] == operation
                    and int(r["word_bits"]) == 64
                ),
                key=lambda r: int(r["n"]),
            )
            if selected:
                ax.plot(
                    [int(r["n"]) for r in selected],
                    [float(r["ns_per_op"]) for r in selected],
                    color=color,
                    linestyle=style,
                    marker=".",
                    label=label,
                    linewidth=1.3,
                )
        ax.set(
            xscale="log",
            yscale="log",
            title=title,
            xlabel="Qubits (2n rows)",
            ylabel="ns per operation",
        )
        ax.grid(alpha=0.15)
    handles, labels = {}, {}
    for ax in axes.flat:
        for handle, label in zip(*ax.get_legend_handles_labels()):
            handles[label], labels[label] = handle, label
    fig.legend(
        handles.values(),
        labels.values(),
        loc="outside lower center",
        ncol=4,
        fontsize=8,
    )
    fig.suptitle(
        f"Tableau operations with 64-bit words\n{caption}",
        fontsize=13,
    )
    fig.savefig(args.directory / "sizes.png", dpi=160)
    fig.savefig(args.directory / "sizes.svg")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for ax, (operation, title) in zip(axes.flat, operations):
        for implementation, layout, label, color, style in series[4:]:
            selected = sorted(
                (
                    r
                    for r in rows
                    if r["implementation"] == implementation
                    and r["layout"] == layout
                    and r["operation"] == operation
                    and int(r["n"]) == args.word_size_n
                ),
                key=lambda r: int(r["word_bits"]),
            )
            ax.plot(
                [int(r["word_bits"]) for r in selected],
                [float(r["ns_per_op"]) for r in selected],
                color=color,
                linestyle=style,
                marker="o",
                label=label,
                linewidth=1.3,
            )
        ax.set_xscale("log", base=2)
        ax.set_xticks([8, 16, 32, 64, 128], ["8", "16", "32", "64", "128"])
        ax.set(yscale="log", title=title, xlabel="Word bits", ylabel="ns per operation")
        ax.grid(alpha=0.15)
    fig.legend(
        *axes.flat[0].get_legend_handles_labels(),
        loc="outside lower center",
        ncol=4,
        fontsize=9,
    )
    fig.suptitle(
        f"Word sizes at n={args.word_size_n}\nQC UInt128 fastrow multiplication is unsupported",
        fontsize=13,
    )
    fig.savefig(args.directory / "words.png", dpi=160)
    fig.savefig(args.directory / "words.svg")


if __name__ == "__main__":
    main()
