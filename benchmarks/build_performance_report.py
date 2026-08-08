"""Build the checked-in core old/new benchmark report."""

from __future__ import annotations

import argparse
import statistics
from pathlib import Path

from summarize_criterion_output import parse


def human(ns: float) -> str:
    if ns >= 1_000_000:
        return f"{ns / 1_000_000:.3f} ms"
    if ns >= 1_000:
        return f"{ns / 1_000:.3f} µs"
    return f"{ns:.3f} ns"


def measurements(
    values: dict[str, list[float]],
) -> tuple[float, float, float, float, float, int]:
    ratios = [new / old for old, new in zip(values["old"], values["new"], strict=False)]
    return (
        statistics.median(values["old"]),
        statistics.median(values["new"]),
        statistics.median(ratios),
        min(ratios),
        max(ratios),
        len(ratios),
    )


def status(ratio: float, runs: int) -> str:
    if ratio > 1.03:
        return "confirmed regression" if runs >= 3 else "provisional regression"
    if ratio < 0.97:
        return "improvement"
    return "parity"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--screening", type=Path, action="append", default=[], required=True
    )
    parser.add_argument(
        "--confirmation", type=Path, action="append", default=[], required=True
    )
    args = parser.parse_args()

    screening: dict[str, dict[str, list[float]]] = {}
    unpaired: list[str] = []
    for path in args.screening:
        measured, missing = parse(path)
        unpaired.extend(missing)
        for key, sides in measured.items():
            target = screening.setdefault(key, {})
            for side, values in sides.items():
                target.setdefault(side, []).extend(values)
    confirmation: dict[str, dict[str, list[float]]] = {}
    for path in args.confirmation:
        measured, _ = parse(path)
        for key, sides in measured.items():
            target = confirmation.setdefault(key, {})
            for side, values in sides.items():
                target.setdefault(side, []).extend(values)
    rows = []
    for key, initial in screening.items():
        if initial.keys() < {"old", "new"}:
            continue
        source = confirmation.get(key, initial)
        if source.keys() < {"old", "new"}:
            source = initial
        old, new, ratio, minimum, maximum, runs = measurements(source)
        rows.append((key, old, new, ratio, minimum, maximum, runs, status(ratio, runs)))

    regressions = sorted(
        (row for row in rows if "regression" in row[-1]),
        key=lambda row: row[3],
        reverse=True,
    )
    improvements = sum(row[-1] == "improvement" for row in rows)
    parity = sum(row[-1] == "parity" for row in rows)
    confirmed = sum(row[-1] == "confirmed regression" for row in rows)
    provisional = sum(row[-1] == "provisional regression" for row in rows)

    out = [
        "# Core old/new benchmark report — 2026-08-07",
        "",
        "This report compares the latest `-2` core crates with the old reference crates.",
        "Every comparable public operation is represented; operations without a valid old",
        "semantic twin are listed separately rather than assigned an artificial ratio.",
        "",
        "## Method",
        "",
        "- Platform: Darwin, release profile, Criterion 0.7.",
        "- Full screening: every comparative target, 10 samples, 0.2 s warm-up, 0.5 s measurement.",
        "- Confirmation: targeted screening regressions and complete surface suites were",
        "  rerun in fresh Cargo-launched processes, 20 samples, 1 s warm-up,",
        "  2 s measurement. The third complete-surface pass was stopped during the",
        "  pathological mixture branch-scaling case; no fourth full pass began.",
        "  Rows with fewer than three measurements remain provisional.",
        "- Ratio is `new / old`; below 0.97 is improvement, 0.97–1.03 is parity.",
        "  Above 1.03 is confirmed with at least three processes and provisional otherwise.",
        "- Setup is excluded unless construction or cloning is itself the target.",
        "",
        "## Summary",
        "",
        f"- Comparable benchmark pairs: **{len(rows)}**",
        f"- Improvements: **{improvements}**",
        f"- Parity: **{parity}**",
        f"- Confirmed regressions: **{confirmed}**",
        f"- Provisional regressions: **{provisional}**",
        "",
        "## Regressions",
        "",
        "| status | ratio | process range | runs | old | new | benchmark |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for key, old, new, ratio, minimum, maximum, runs, result in regressions:
        out.append(
            f"| {result} | **{ratio:.3f}×** | {minimum:.3f}–{maximum:.3f}× | {runs} | "
            f"{human(old)} | {human(new)} | `{key}` |"
        )

    out.extend(
        [
            "",
            "## Attribution summary",
            "",
            "- **Lossy branch keys (1.6–1.8×):** new clone-then-toggle copies three",
            "  atomic caches and performs guarded invalidation; old copies one plain",
            "  digest and uses unchecked setters. The benchmark excludes the later hash.",
            "- **Disabled truncation (1.5–1.7×):** both sides return immediately;",
            "  the roughly 1 ns absolute delta is confirmed but remains unattributed.",
            "- **Pauli-sum construction/clone:** new preallocates persistent scratch and",
            "  clones atomic key caches; this explains much of the 1.27–1.34× delta.",
            "- **Pauli-sum Clifford families and qubit scaling:** drift is localized to",
            "  repeated bijective re-keying. Existing controls exclude reserve and show",
            "  only a small drain/clone contribution; the remaining mechanism is unknown.",
            "- **Lossy-sum integration:** loss, reset, Clifford and rotation stages remain",
            "  slower; the component-cache representation is a source difference, but no",
            "  controlled ablation proves it is the cause.",
            "- **Symbolic propagation/trace:** matched fixtures confirm the regression,",
            "  but the remaining engine/trace mechanism is unattributed. Symbolic eval",
            "  is slower on one-use variables because the new angle cache initializes 32",
            "  entries and computes both sine and cosine on each miss.",
            "- **Tiny tableau observation helpers (2.3–4.0×):** implementations are",
            "  effectively identical and absolute deltas are 3–6 ns; this is attributed",
            "  to inlining/code placement or struct-layout effects, not algorithmic work.",
            "- **Tableau many-gate batches (about 1.06–1.12×):** new adds lazy-hash",
            "  invalidation and uses a different inner loop shape; full scaling workloads",
            "  are otherwise parity or faster.",
            "- **Mixture clone/measurement:** new persistently clones fingerprint buckets",
            "  and measurement eagerly dirties/rebuilds fingerprints. Two-qubit mixture",
            "  noise repeatedly scans rows for each of 15 Pauli branches.",
            "",
            "## Complete comparison table",
            "",
            "| status | ratio | runs | old | new | benchmark |",
            "|---|---:|---:|---:|---:|---|",
        ]
    )
    for key, old, new, ratio, _, _, runs, result in sorted(rows):
        out.append(
            f"| {result} | {ratio:.3f}× | {runs} | {human(old)} | "
            f"{human(new)} | `{key}` |"
        )

    out.extend(
        [
            "",
            "## No-old-twin measurements and exclusions",
            "",
            "- New-only Pauli-sum `reduce`, Hermitian overlap, and sum×sum multiplication:",
            "  old sum-RHS multiplication is uninstantiable even for a singleton RHS;",
            "  `multiply_into` has no old callable semantic twin.",
            "- Mixed/non-unit projection is new-only because Lean adjudicated the old behavior",
            "  as incorrect; only the common I/Z unit-coefficient subset is compared.",
            "- Lossy Pauli words intentionally have no native product.",
            "- Complex symbolic evaluation and exact Gaussian-ring operations are new-only.",
            "- Symbolic projection is blocked because the new symbolic term deliberately",
            "  lacks `Halvable`; symbolic amplitude damping lacks `Float`; neither engine",
            "  implements `PauliErrorAll` for symbolic sums.",
            "- Direct tableau sampling has no common old/new API; mixture sampling is compared.",
            "",
            f"Unpaired measured benchmark IDs in the screening output: **{len(set(unpaired))}**.",
        ]
    )
    args.output.write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
