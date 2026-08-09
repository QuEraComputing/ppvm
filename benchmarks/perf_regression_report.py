# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Run the `ppvm-conformance-2` head-to-head benchmarks and report old-vs-`-2`
regressions.

Every benchmark target in `crates/ppvm-conformance-2` holds *both* engines in
one binary, so each `<name>/old` + `<name>/new` pair is a same-build ratio and
the code-layout bias cancels. This script drives those targets, pairs the
medians (reusing `summarize_criterion_output.parse`), classifies each pair
against the audit gate, and fails if a regression survives.

Ratios are `new / old`:

    < 0.97      improvement
    0.97-1.03   parity
    > 1.03      regression; *robust* when the slowest process is also above the
                gate, which is the only class that fails the run

Run it with `mise run perf-report` — see that task for the common invocations.
"""

from __future__ import annotations

import argparse
import fnmatch
import importlib.util
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CONFORMANCE = ROOT / "crates" / "ppvm-conformance-2"
PACKAGE = "ppvm-conformance-2"
DEFAULT_ALLOWLIST = ROOT / "benchmarks" / "perf-allowlist.txt"

# The audit protocol: short screening runs, repeated in fresh processes. A
# single long process cannot separate a real regression from an executable
# layout artifact, which is why the gate needs `--launches > 1` to bite.
SCREEN = {"sample_size": 20, "warm_up_time": 1.0, "measurement_time": 2.0}

IMPROVED = 0.97
GATE = 1.03


def load_parser():
    """Import `parse()` from the sibling summarizer so both agree on units."""
    path = Path(__file__).with_name("summarize_criterion_output.py")
    spec = importlib.util.spec_from_file_location("summarize_criterion_output", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.parse


def discover_targets() -> list[str]:
    """Every `[[bench]]` in the conformance manifest, in declaration order."""
    manifest = (CONFORMANCE / "Cargo.toml").read_text()
    names = re.findall(r"\[\[bench\]\]\s*\nname\s*=\s*\"([^\"]+)\"", manifest)
    # `harness_smoke` only checks that the generators still link; it has no
    # old/new pair to compare and only costs wall time here.
    return [name for name in names if name != "harness_smoke"]


def build(env: dict[str, str]) -> None:
    subprocess.run(
        ["cargo", "build", "--release", "-p", PACKAGE, "--benches"],
        cwd=ROOT,
        env=env,
        check=True,
    )


def bench_command(target: str, args) -> list[str]:
    command = ["cargo", "bench", "-q", "-p", PACKAGE, "--bench", target, "--"]
    command += ["--noplot", "--discard-baseline"]
    if not args.full:
        command += [
            "--sample-size", str(args.sample_size),
            "--warm-up-time", str(args.warm_up_time),
            "--measurement-time", str(args.measurement_time),
        ]
    return command


def count_cases(target: str, args, env: dict[str, str]) -> int:
    """How many Criterion cases `target` will run under the current filter.

    `--list` enumerates them without measuring, which is what lets the progress
    bar know its denominator — and the ETA its horizon — before the first
    sample. It still executes each target's fixture setup, so it is not free,
    just cheap relative to the sweep it is sizing.
    """
    command = bench_command(target, args)
    command.remove("--noplot")
    command.remove("--discard-baseline")
    command.append("--list")
    if args.filter:
        command.append(args.filter)
    result = subprocess.run(
        command, cwd=ROOT, env=env, capture_output=True, text=True, check=False
    )
    return sum(1 for line in result.stdout.splitlines() if line.endswith(": benchmark"))


# One Criterion case prints exactly one of these when it stops measuring, which
# makes it the progress tick. The id is everything between the prefix and the
# suffix.
ANALYZING = re.compile(r"^Benchmarking (?P<id>.+): Analyzing$")


class Progress:
    """A single-line bar over `total` Criterion cases, with a measured ETA.

    The remaining time is extrapolated from the cases already finished rather
    than from the nominal per-case cost, so per-target fixture setup and the
    targets that override the screening durations are absorbed automatically —
    the estimate is rough at the start and tightens as the sweep proceeds.
    """

    def __init__(self, total: int, nominal: float) -> None:
        self.total = total
        self.nominal = nominal
        self.done = 0
        self.start = time.monotonic()
        self.tty = sys.stdout.isatty()
        self.last_line = 0.0

    @staticmethod
    def clock(seconds: float) -> str:
        seconds = max(0, int(seconds))
        hours, seconds = divmod(seconds, 3600)
        minutes, seconds = divmod(seconds, 60)
        return f"{hours}h{minutes:02d}m" if hours else f"{minutes}m{seconds:02d}s"

    def eta(self) -> float:
        elapsed = time.monotonic() - self.start
        remaining = self.total - self.done
        if self.done < 3:
            return remaining * self.nominal
        return remaining * (elapsed / self.done)

    def render(self, note: str) -> None:
        elapsed = time.monotonic() - self.start
        fraction = self.done / self.total if self.total else 1.0
        width = shutil.get_terminal_size((100, 24)).columns
        bar_width = 24
        filled = int(bar_width * fraction)
        bar = "█" * filled + "░" * (bar_width - filled)
        head = (
            f"[{bar}] {fraction * 100:3.0f}% │ {self.done}/{self.total} │ "
            f"{self.clock(elapsed)} elapsed │ ~{self.clock(self.eta())} left │ "
        )
        if self.tty:
            line = (head + note)[: max(width - 1, 20)]
            print(f"\r\033[K{line}", end="", flush=True)
            return
        # Non-interactive (CI, `| tee`): no carriage returns, and throttled so
        # the log stays readable.
        now = time.monotonic()
        if now - self.last_line >= 30 or self.done in (0, self.total):
            self.last_line = now
            print(head + note, flush=True)

    def tick(self, note: str) -> None:
        self.done += 1
        self.render(note)

    def finish(self) -> None:
        if self.tty:
            print("\r\033[K", end="", flush=True)


def run_target(target: str, launch: int, args, log, env: dict[str, str], bar: Progress) -> None:
    """One fresh process for one target, its output appended to `log`.

    The `Running benches/<target>.rs` header is written by us rather than by
    cargo: `-q` suppresses cargo's own line, and the parser keys per-process
    measurements off that header, so synthesising it is what keeps repeated
    launches from overwriting each other.
    """
    command = bench_command(target, args)
    if args.filter:
        command.append(args.filter)

    log.write(f"Running benches/{target}.rs (launch {launch})\n")
    log.flush()
    label = f"{target} ({launch}/{args.launches})"
    bar.render(label)
    process = subprocess.Popen(
        command,
        cwd=ROOT,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    for line in process.stdout:
        log.write(line)
        match = ANALYZING.match(line.rstrip("\n"))
        if match:
            bar.tick(f"{label} {match.group('id')}")
    process.wait()
    log.flush()
    if process.returncode != 0:
        bar.finish()
        print(f"warning: {target} exited {process.returncode}", file=sys.stderr)


def load_allowlist(path: Path) -> list[str]:
    if not path.exists():
        return []
    patterns = []
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            patterns.append(line)
    return patterns


def allowlisted(key: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatch(key, pattern) for pattern in patterns)


def classify(ratio: float, minimum: float, runs: int, gate: float) -> str:
    """`robust` needs the slowest of *several* processes to stay above the gate.

    One process cannot supply that evidence — its minimum is its median — so a
    single-launch sweep tops out at `regression`, which reports but never
    blocks. That is the whole reason the gate wants `--launches >= 2`.
    """
    if ratio > gate:
        return "robust" if runs >= 2 and minimum > gate else "regression"
    if ratio < IMPROVED:
        return "improved"
    return "parity"


def summarize(pairs: dict, gate: float, allowlist: list[str]) -> list[dict]:
    rows = []
    for key, values in pairs.items():
        ratios = [new / old for old, new in zip(values["old"], values["new"])]
        if not ratios:
            continue
        median = statistics.median(ratios)
        rows.append(
            {
                "key": key,
                "ratio": median,
                "min": min(ratios),
                "max": max(ratios),
                "runs": len(ratios),
                "old_ns": statistics.median(values["old"]),
                "new_ns": statistics.median(values["new"]),
                "class": classify(median, min(ratios), len(ratios), gate),
                "allowlisted": allowlisted(key, allowlist),
            }
        )
    rows.sort(key=lambda row: row["ratio"], reverse=True)
    return rows


def table(rows: list[dict]) -> list[str]:
    lines = [
        "| ratio | min | max | runs | old | new | benchmark |",
        "|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        note = " *(allowlisted)*" if row["allowlisted"] else ""
        lines.append(
            f"| {row['ratio']:.3f}× | {row['min']:.3f}× | {row['max']:.3f}× "
            f"| {row['runs']} | {row['old_ns']:.1f} ns | {row['new_ns']:.1f} ns "
            f"| `{row['key']}`{note} |"
        )
    return lines


def write_report(path: Path, rows: list[dict], unpaired: list[str], args, gate: float) -> None:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=False
    ).stdout.strip()
    counts = {name: 0 for name in ("robust", "regression", "parity", "improved")}
    for row in rows:
        counts[row["class"]] += 1
    blocking = [row for row in rows if row["class"] == "robust" and not row["allowlisted"]]
    protocol = (
        "each target's own Criterion config"
        if args.full
        else f"{args.sample_size} samples, {args.warm_up_time} s warm-up, "
        f"{args.measurement_time} s measurement"
    )

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        f"# ppvm `-2` vs legacy performance report — {stamp}",
        "",
        f"- Commit: `{commit}`",
        f"- Host: {platform.platform()}",
        f"- Protocol: {protocol}; {args.launches} launch(es) per target, fresh process each",
        f"- Targets: {', '.join(args.bench)}",
        f"- Filter: `{args.filter}`" if args.filter else "- Filter: (none)",
        f"- Gate: `new/old` median and process minimum both above {gate:.2f}×",
        "",
        "## Summary",
        "",
        "| class | pairs |",
        "|---|---:|",
        f"| robust regression | {counts['robust']} |",
        f"| regression (process-unstable) | {counts['regression']} |",
        f"| parity | {counts['parity']} |",
        f"| improved | {counts['improved']} |",
        f"| **total paired** | **{len(rows)}** |",
        "",
        f"Unpaired measurements (new-only, old-only, or non-`{{side}}` labels): {len(unpaired)}",
        "",
    ]

    if args.launches < 2:
        lines += [
            (
                "> Screening run (one launch per target). Nothing here can be "
                "*robust*: with a single process the minimum is the median, so no "
                "row blocks. Rerun the rows below with `--launches 4` to decide "
                "whether they are real regressions or executable-layout artifacts."
            ),
            "",
        ]

    lines += ["## Blocking regressions", ""]
    lines += table(blocking) if blocking else ["None.", ""]

    unstable = [row for row in rows if row["class"] == "regression"]
    heading = (
        "## Above gate — needs confirmation"
        if args.launches < 2
        else "## Above gate but process-unstable"
    )
    lines += ["", heading, ""]
    lines += table(unstable) if unstable else ["None.", ""]

    lines += ["", "## All pairs, slowest first", ""] + table(rows) + [""]

    if unpaired:
        lines += ["## Unpaired measurements", "", "```"] + unpaired + ["```", ""]

    path.write_text("\n".join(lines) + "\n")


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w") as handle:
        handle.write("ratio\tmin\tmax\truns\told_ns\tnew_ns\tclass\tbenchmark\n")
        for row in rows:
            handle.write(
                f"{row['ratio']:.4f}\t{row['min']:.4f}\t{row['max']:.4f}\t{row['runs']}\t"
                f"{row['old_ns']:.3f}\t{row['new_ns']:.3f}\t{row['class']}\t{row['key']}\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--bench", action="append", metavar="TARGET",
        help="Benchmark target to run; repeatable. Default: every conformance target.",
    )
    parser.add_argument(
        "--filter", metavar="SUBSTRING",
        help="Criterion filter, e.g. 'noise' — narrows every target to matching ids.",
    )
    parser.add_argument(
        "--launches", type=int, default=1, metavar="N",
        help="Fresh processes per target. Use >=4 to confirm a regression (default: 1).",
    )
    parser.add_argument("--sample-size", type=int, default=SCREEN["sample_size"])
    parser.add_argument("--warm-up-time", type=float, default=SCREEN["warm_up_time"])
    parser.add_argument("--measurement-time", type=float, default=SCREEN["measurement_time"])
    parser.add_argument(
        "--full", action="store_true",
        help="Use each target's own Criterion config instead of the screening protocol.",
    )
    parser.add_argument(
        "--out", type=Path, default=ROOT / "target" / "perf-report", metavar="DIR",
        help="Output directory (default: target/perf-report).",
    )
    parser.add_argument("--gate", type=float, default=GATE, metavar="RATIO")
    parser.add_argument("--allowlist", type=Path, default=DEFAULT_ALLOWLIST, metavar="PATH")
    parser.add_argument(
        "--reuse", type=Path, metavar="RAW",
        help="Skip benchmarking and re-summarize an existing raw log.",
    )
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument(
        "--no-gate", action="store_true",
        help="Always exit 0; report only.",
    )
    args = parser.parse_args()

    if shutil.which("cargo") is None:
        print("error: cargo not found on PATH", file=sys.stderr)
        return 2

    args.bench = args.bench or discover_targets()
    args.out.mkdir(parents=True, exist_ok=True)
    raw = args.reuse or (args.out / "raw.txt")

    if args.reuse is None:
        # Criterion's own runtime is the measurement, so keep cargo from
        # rebuilding mid-sweep and keep rayon out of the ratio unless a target
        # asked for it.
        env = dict(os.environ)
        env.setdefault("CARGO_TERM_COLOR", "never")
        if not args.no_build:
            print("Building benchmarks (release)…", flush=True)
            build(env)

        print(f"Counting cases in {len(args.bench)} target(s)…", flush=True)
        counts_per_target = {target: count_cases(target, args, env) for target in args.bench}
        empty = [target for target, count in counts_per_target.items() if count == 0]
        for target in empty:
            print(f"  note: {target} has no case matching the filter; skipping", flush=True)
            args.bench.remove(target)
        total = sum(counts_per_target.values()) * args.launches
        if total == 0:
            print("error: no benchmark case matched", file=sys.stderr)
            return 2

        # Warm-up + measurement is the floor; the rest is Criterion's analysis
        # and the harness's own per-case setup. Only the first few cases use it,
        # after which the bar extrapolates from real timings.
        nominal = (
            8.6 if args.full else args.warm_up_time + args.measurement_time + 0.6
        )
        bar = Progress(total, nominal)
        print(
            f"{total} case(s) across {len(args.bench)} target(s) × {args.launches} launch(es); "
            f"rough estimate {Progress.clock(total * nominal)} → {raw}",
            flush=True,
        )
        with raw.open("w") as log:
            for launch in range(1, args.launches + 1):
                for target in args.bench:
                    run_target(target, launch, args, log, env, bar)
        bar.finish()
        print(f"Benchmarking finished in {Progress.clock(time.monotonic() - bar.start)}.")

    pairs, unpaired = load_parser()(raw)
    allowlist = load_allowlist(args.allowlist)
    rows = summarize(pairs, args.gate, allowlist)

    report = args.out / "report.md"
    tsv = args.out / "pairs.tsv"
    write_report(report, rows, unpaired, args, args.gate)
    write_tsv(tsv, rows)

    blocking = [row for row in rows if row["class"] == "robust" and not row["allowlisted"]]
    counts = {name: sum(row["class"] == name for row in rows) for name in
              ("robust", "regression", "parity", "improved")}
    above = "unconfirmed" if args.launches < 2 else "process-unstable"
    print(
        f"\n{len(rows)} paired: {counts['improved']} improved, {counts['parity']} parity, "
        f"{counts['regression']} above gate ({above}), "
        f"{counts['robust']} robust regression(s)."
    )
    for row in rows:
        if row["ratio"] > args.gate:
            print(f"  {row['ratio']:.3f}× ({row['min']:.3f}–{row['max']:.3f}) {row['key']}")
    print(f"\nReport: {report}\nTSV:    {tsv}\nRaw:    {raw}")
    if args.launches < 2 and counts["regression"]:
        print(
            "\nSingle-launch screening: rerun the rows above with "
            "`--launches 4` (add `--filter`) before treating them as regressions."
        )
    sys.stdout.flush()

    if blocking and not args.no_gate:
        print(
            f"\nFAIL: {len(blocking)} robust regression(s) above {args.gate:.2f}×.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
