"""Summarize paired old/new Criterion medians from terminal output."""

from __future__ import annotations

import argparse
import itertools
import re
import statistics
from pathlib import Path

TIME = re.compile(
    r"^\s*time:\s+\[[^ ]+\s+\S+\s+(?P<median>[0-9.]+)\s+"
    r"(?P<unit>ps|ns|us|µs|ms|s)\s+"
)
TARGET = re.compile(r"Running benches/(?P<name>[^ ]+)\.rs ")
SIDE = re.compile(r"(?:(?<=/)|(?<=_))(old|new)(?=/|_|$)")
SCALE = {
    "ps": 0.001,
    "ns": 1.0,
    "us": 1_000.0,
    "µs": 1_000.0,
    "ms": 1_000_000.0,
    "s": 1_000_000_000.0,
}


def parse(path: Path) -> tuple[dict[str, dict[str, list[float]]], list[str]]:
    lines = path.read_text().splitlines()
    measured: dict[tuple[str, int], dict[str, float]] = {}
    launches: dict[str, int] = {}
    current = ("unknown", 0)
    for label, timing in itertools.pairwise(lines):
        target_match = TARGET.search(label)
        if target_match is not None:
            name = target_match.group("name")
            launches[name] = launches.get(name, 0) + 1
            current = (name, launches[name])
            continue
        if not label or label[0].isspace() or label.startswith(("Benchmarking ", "Found ")):
            continue
        time_match = TIME.match(timing)
        if time_match is None:
            continue
        median_ns = float(time_match.group("median")) * SCALE[time_match.group("unit")]
        measured.setdefault(current, {})[label] = median_ns

    pairs: dict[str, dict[str, list[float]]] = {}
    unpaired: list[str] = []
    for run, run_measurements in measured.items():
        paired_labels: set[str] = set()
        for label in run_measurements:
            if label in paired_labels:
                continue
            for match in SIDE.finditer(label):
                key = f"{label[: match.start(1)]}{{side}}{label[match.end(1) :]}"
                old_label = key.format(side="old")
                new_label = key.format(side="new")
                if old_label in run_measurements and new_label in run_measurements:
                    target = pairs.setdefault(key, {"old": [], "new": []})
                    target["old"].append(run_measurements[old_label])
                    target["new"].append(run_measurements[new_label])
                    paired_labels.update((old_label, new_label))
                    break
        unpaired.extend(
            f"{run[0]}#{run[1]}:{label}"
            for label in run_measurements
            if label not in paired_labels
        )
    return pairs, unpaired


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--only-regressions", action="store_true")
    args = parser.parse_args()

    pairs, unpaired = parse(args.output)
    rows = []
    paired_count = 0
    for key, values in pairs.items():
        if values.keys() >= {"old", "new"}:
            paired_count += 1
            process_ratios = [
                new / old for old, new in zip(values["old"], values["new"], strict=False)
            ]
            ratio = statistics.median(process_ratios)
            if not args.only_regressions or ratio > 1.0:
                rows.append(
                    (
                        ratio,
                        key,
                        statistics.median(values["old"]),
                        statistics.median(values["new"]),
                        min(process_ratios),
                        max(process_ratios),
                        len(process_ratios),
                    )
                )
    rows.sort(reverse=True)

    print("ratio\tmin\tmax\truns\told_ns\tnew_ns\tbenchmark")
    for ratio, key, old, new, minimum, maximum, runs in rows:
        print(
            f"{ratio:.4f}\t{minimum:.4f}\t{maximum:.4f}\t{runs}\t"
            f"{old:.3f}\t{new:.3f}\t{key}"
        )
    print(
        f"\npaired={paired_count} shown={len(rows)} "
        f"unpaired_measurements={len(unpaired)}"
    )


if __name__ == "__main__":
    main()
