# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Drive the cross-library Pauli-propagation benchmark and merge the results.

Runs the same two circuits — TFIM Trotter and Heisenberg-model correlations —
through `ppvm-*-2`, PauliPropagation.jl, PauliStrings.jl and Qiskit's
`pauli-prop`, all from one parameter contract, and writes one tidy CSV.

Before it will report a single timing it **validates**: every engine dumps its
propagated support at a small width and the driver diffs them term-for-term
against `ppvm-2`. A cross-library benchmark whose engines are quietly computing
different things is worse than no benchmark, and the two bugs this check caught
while the harness was being written (a reversed circuit in the `pauli-prop`
runner, duplicate Paulis silently collapsing in its readout) both produced
plausible-looking numbers that agreed to 9 digits on the observable.

    uv run --no-project --with pauli-prop python3 run_xbench.py --help
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
JULIA_PROJECT = REPO / "julia-benchmarks"
CSV_COLUMNS = [
    "model",
    "library",
    "qubits",
    "steps",
    "dt",
    "atol",
    "time_s",
    "terms",
    "observable",
]


@dataclass(frozen=True)
class Runner:
    """One engine: how to invoke it, and whether it is available here."""

    name: str
    argv: list[str]
    needs: str  # the executable that must exist on PATH

    def available(self) -> bool:
        return shutil.which(self.needs) is not None


RUNNERS = {
    "ppvm-2": Runner(
        "ppvm-2",
        [
            "cargo",
            "run",
            "--release",
            "-q",
            "-p",
            "ppvm-pauli-sum-2",
            "--example",
            "xbench",
        ],
        "cargo",
    ),
    "ppvm-1": Runner(
        "ppvm-1",
        [
            "cargo",
            "run",
            "--release",
            "-q",
            "-p",
            "ppvm-pauli-sum",
            "--example",
            "xbench_v1",
        ],
        "cargo",
    ),
    "pauli-propagation-jl": Runner(
        "pauli-propagation-jl",
        [
            "julia",
            f"--project={JULIA_PROJECT}",
            "-t1",
            str(JULIA_PROJECT / "benches" / "xbench_pp.jl"),
        ],
        "julia",
    ),
    "pauli-strings-jl": Runner(
        "pauli-strings-jl",
        [
            "julia",
            f"--project={JULIA_PROJECT}",
            "-t1",
            str(JULIA_PROJECT / "benches" / "xbench_ps.jl"),
        ],
        "julia",
    ),
    "pauli-prop": Runner(
        "pauli-prop",
        [
            "uv",
            "run",
            "--no-project",
            "--with",
            "pauli-prop",
            "python3",
            str(Path(__file__).with_name("xbench_qiskit.py")),
        ],
        "uv",
    ),
}


def invoke(runner: Runner, env_extra: dict[str, str], quiet: bool) -> str:
    """Run one engine and return its stdout, streaming its stderr as progress."""
    env = os.environ.copy()
    env.update(env_extra)
    proc = subprocess.run(
        runner.argv,
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"{runner.name} failed with exit code {proc.returncode}")
    if not quiet and proc.stderr.strip():
        sys.stderr.write(proc.stderr)
    return proc.stdout


def parse_dump(text: str) -> dict[str, float]:
    """Parse a `DUMP=1` support listing into `{word: coefficient}`."""
    terms: dict[str, float] = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        word, coeff = line.split()
        terms[word] = float(coeff)
    return terms


def validate(
    libs: list[str], params: dict[str, str], atol_cmp: float, quiet: bool
) -> None:
    """Assert every engine propagates the *same* operator, term for term."""
    print("validating: all engines must agree term-for-term", file=sys.stderr)
    for model in ("tfim", "heisenberg"):
        env = dict(params)
        env.update(
            {"MODEL": model, "QUBITS": "4", "STEPS": "3", "ATOL": "1e-14", "DUMP": "1"}
        )
        reference: dict[str, float] | None = None
        for lib in libs:
            terms = parse_dump(invoke(RUNNERS[lib], env, quiet=True))
            if reference is None:
                reference, ref_lib = terms, lib
                print(
                    f"  {model}: {lib} — {len(terms)} terms (reference)",
                    file=sys.stderr,
                )
                continue
            missing = set(reference) - set(terms)
            extra = set(terms) - set(reference)
            worst = max(
                (abs(reference[w] - terms[w]) for w in set(reference) & set(terms)),
                default=0.0,
            )
            if missing or extra or worst > atol_cmp:
                raise SystemExit(
                    f"{model}: {lib} disagrees with {ref_lib} — "
                    f"{len(terms)} vs {len(reference)} terms, "
                    f"{len(missing)} missing, {len(extra)} extra, max|Δc|={worst:.3e} "
                    f"(tolerance {atol_cmp:.0e}). Refusing to report timings."
                )
            print(
                f"  {model}: {lib} — {len(terms)} terms, max|Δc|={worst:.1e} OK",
                file=sys.stderr,
            )
    print("validation passed\n", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--model", default="both", choices=["tfim", "heisenberg", "both"])
    ap.add_argument(
        "--qubits", default="8,12,16,20,24,28,32", help="applies to every model"
    )
    ap.add_argument(
        "--qubits-tfim",
        default=None,
        help="override --qubits for the TFIM sweep (its support grows linearly in n, "
        "so it reaches far wider systems than the Heisenberg one at equal cost)",
    )
    ap.add_argument(
        "--qubits-heisenberg", default=None, help="override --qubits for Heisenberg"
    )
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--j", type=float, default=1.0, help="bond coupling J")
    ap.add_argument("--h", type=float, default=1.0, help="field strength h")
    ap.add_argument(
        "--atol", type=float, default=1e-6, help="coefficient truncation threshold"
    )
    ap.add_argument(
        "--iters", type=int, default=3, help="timed repeats; the minimum is reported"
    )
    ap.add_argument("--libs", default=",".join(RUNNERS))
    ap.add_argument("--out", default="target/xbench", type=Path)
    ap.add_argument("--skip-validate", action="store_true")
    ap.add_argument(
        "--validate-tol",
        type=float,
        default=1e-10,
        help="absolute coefficient tolerance for the cross-engine agreement check",
    )
    ap.add_argument(
        "--reuse",
        type=Path,
        default=None,
        help="re-print the summary from an existing results.csv",
    )
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()

    if args.reuse is not None:
        with args.reuse.open() as fh:
            summarize(list(csv.DictReader(fh)))
        return

    libs = [lib.strip() for lib in args.libs.split(",") if lib.strip()]
    for lib in libs:
        if lib not in RUNNERS:
            raise SystemExit(f"unknown library {lib!r}; known: {', '.join(RUNNERS)}")
    missing = [lib for lib in libs if not RUNNERS[lib].available()]
    if missing:
        print(
            f"skipping (missing {', '.join(RUNNERS[m].needs for m in missing)}): "
            f"{', '.join(missing)}",
            file=sys.stderr,
        )
        libs = [lib for lib in libs if lib not in missing]
    if not libs:
        raise SystemExit("no runnable engines")

    params = {
        "DT": repr(args.dt),
        "JCOUP": repr(args.j),
        "HFIELD": repr(args.h),
        "ITERS": str(args.iters),
    }
    if not args.skip_validate:
        validate(libs, params, args.validate_tol, args.quiet)

    models = ["tfim", "heisenberg"] if args.model == "both" else [args.model]
    per_model_qubits = {
        "tfim": args.qubits_tfim or args.qubits,
        "heisenberg": args.qubits_heisenberg or args.qubits,
    }
    rows: list[dict[str, str]] = []
    for model in models:
        for lib in libs:
            env = dict(params)
            env.update(
                {
                    "MODEL": model,
                    "QUBITS": per_model_qubits[model],
                    "STEPS": str(args.steps),
                    "ATOL": repr(args.atol),
                }
            )
            print(f"running {model} / {lib}", file=sys.stderr)
            out = invoke(RUNNERS[lib], env, args.quiet)
            reader = csv.DictReader(line for line in out.splitlines() if line.strip())
            rows.extend(dict(row) for row in reader)

    out_dir = REPO / args.out if not args.out.is_absolute() else args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "results.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    summarize(rows)
    print(f"\nwrote {csv_path}", file=sys.stderr)


def summarize(rows: list[dict[str, str]]) -> None:
    """Print a time table plus the speedup of every engine relative to `ppvm-2`."""
    by_model: dict[str, dict[int, dict[str, dict[str, str]]]] = {}
    for row in rows:
        by_model.setdefault(row["model"], {}).setdefault(int(row["qubits"]), {})[
            row["library"]
        ] = row

    for model, per_n in by_model.items():
        libs = sorted({lib for cells in per_n.values() for lib in cells})
        print(f"\n=== {model} ===")
        header = f"{'n':>4}  {'terms':>9}  " + "  ".join(f"{lib:>22}" for lib in libs)
        print(header)
        print("-" * len(header))
        for n in sorted(per_n):
            cells = per_n[n]
            base = cells.get("ppvm-2")
            terms = base["terms"] if base else next(iter(cells.values()))["terms"]
            line = f"{n:>4}  {terms:>9}  "
            parts = []
            for lib in libs:
                row = cells.get(lib)
                if row is None:
                    parts.append(f"{'—':>22}")
                    continue
                t = float(row["time_s"])
                if base is not None and lib != "ppvm-2":
                    ratio = t / float(base["time_s"])
                    parts.append(f"{t:>11.4f}s ({ratio:5.2f}x)")
                else:
                    parts.append(f"{t:>11.4f}s {'(1.00x)':>9}")
            print(line + "  ".join(parts))

        # A runtime ratio only means something at equal work. Any engine whose
        # truncation rule leaves it carrying a materially different support is
        # called out here rather than left to the reader to spot in the plot.
        drift: dict[str, list[str]] = {}
        for n in sorted(per_n):
            base = per_n[n].get("ppvm-2")
            if base is None:
                continue
            for lib, row in per_n[n].items():
                if lib == "ppvm-2":
                    continue
                rel = int(row["terms"]) / int(base["terms"]) - 1.0
                if abs(rel) > 0.02:
                    drift.setdefault(lib, []).append(f"n={n}: {rel:+.0%}")
        for lib, notes in drift.items():
            print(f"  note: {lib} support vs ppvm-2 — {', '.join(notes)}")


if __name__ == "__main__":
    main()
