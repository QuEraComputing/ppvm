# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Measure what each engine's truncation rule costs in accuracy, not in time.

`run_xbench.py` answers "how fast", holding `atol` fixed. It cannot answer "how
wrong", because the engines do not all mean the same thing by `atol`: `ppvm`,
PauliPropagation.jl and PauliStrings.jl accumulate every contribution to a term
and then drop the merged coefficient, while `pauli-prop` and monoprop test a
branch's prospective coefficient before they emit it. Same threshold, different
rule, so the fair question is how far each lands from the untruncated answer.

This sweeps `atol` at a small width and diffs each engine's **whole coefficient
vector** against a converged reference. The vector norm is the point. Judging
this on the scalar `observable` is actively misleading, because its truncation
error changes sign — `ppvm` on TFIM walks 7.7e-4, 1.8e-3, 1.0e-4, 3.0e-7, 1.2e-6
over `atol` 1e-3 … 1e-7, and that accidental near-zero at 1e-6 manufactures a 27x
gap against monoprop where the norm shows 6 %. Both are recorded so the
discrepancy stays visible rather than being something you have to rediscover.

The reference is the baseline engine at `--ref-atol`. Widths whose reachable
sector saturates give an exact reference (Heisenberg `n=8` tops out at 16 384
terms); elsewhere check it has converged before trusting a small ratio.

    uv run --no-project python3 xbench_accuracy.py --help
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import sys
from pathlib import Path

from run_xbench import BASELINE, RUNNERS, invoke, parse_dump

CSV_COLUMNS = [
    "model",
    "library",
    "qubits",
    "steps",
    "dt",
    "seed",
    "atol",
    "ref_atol",
    "ref_terms",
    "ref_norm",
    "terms",
    "l2_err",
    "l1_err",
    "max_coeff_err",
    "lost_mass",
    "kept_subthreshold",
    "dropped_above_atol",
    "observable",
    "ref_observable",
    "obs_abs_err",
]


def norm(terms: dict[str, float]) -> float:
    """`‖c‖₂`. Unitary conjugation preserves it, so it is the scale the absolute
    errors should be read against."""
    return math.sqrt(sum(v * v for v in terms.values()))


def measure(approx: dict[str, float], ref: dict[str, float], atol: float) -> dict:
    """Compare one propagated support against the converged one."""
    words = set(approx) | set(ref)
    diffs = [abs(approx.get(w, 0.0) - ref.get(w, 0.0)) for w in words]
    return {
        "terms": len(approx),
        "l2_err": math.sqrt(sum(d * d for d in diffs)),
        "l1_err": sum(diffs),
        "max_coeff_err": max(diffs, default=0.0),
        # Weight the engine discarded that the reference says was there.
        "lost_mass": sum(abs(c) for w, c in ref.items() if w not in approx),
        # Rows held whose converged value is under the threshold they were given:
        # the signature of a branch that cleared `atol` when it was emitted and
        # was never re-tested after later contributions cancelled it.
        "kept_subthreshold": sum(1 for w in approx if abs(ref.get(w, 0.0)) < atol),
        # The opposite failure: mass thrown away that belonged above `atol`.
        "dropped_above_atol": sum(
            1 for w, c in ref.items() if w not in approx and abs(c) >= atol
        ),
    }


def observable(text: str) -> float:
    """Pull the `observable` column out of a runner's one-row CSV."""
    return float(
        next(csv.DictReader(l for l in text.splitlines() if l.strip()))["observable"]
    )


def run(lib: str, env: dict[str, str], quiet: bool) -> tuple[dict[str, float], float]:
    """Get one engine's support and its readout at these parameters."""
    dump = parse_dump(invoke(RUNNERS[lib], {**env, "DUMP": "1"}, quiet))
    return dump, observable(invoke(RUNNERS[lib], env, quiet))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--models", default="tfim,heisenberg")
    ap.add_argument("--qubits", type=int, default=8)
    ap.add_argument(
        "--steps",
        default="10",
        help="comma-separated depths; the reference is recomputed per depth",
    )
    ap.add_argument("--dt", default="0.1", help="comma-separated angles scales")
    ap.add_argument("--atols", default="1e-3,1e-4,1e-5,1e-6,1e-7")
    ap.add_argument(
        "--seeds",
        default="12345",
        help="comma-separated circuit seeds; only `scramble` reads them, and it "
        "needs several because a random instance varies",
    )
    ap.add_argument(
        "--ref-atol",
        default="1e-16",
        help="threshold for the reference run; must be converged at this width",
    )
    ap.add_argument("--libs", default=",".join(RUNNERS))
    ap.add_argument("--out", type=Path, default=Path("target/xbench/accuracy.csv"))
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    libs = [lib for lib in args.libs.split(",") if lib.strip()]
    if unknown := set(libs) - set(RUNNERS):
        raise SystemExit(f"unknown libraries: {sorted(unknown)}")
    if missing := [lib for lib in libs if not RUNNERS[lib].available()]:
        print(f"skipping (toolchain absent): {', '.join(missing)}", file=sys.stderr)
        libs = [lib for lib in libs if lib not in missing]
    if BASELINE not in libs:
        raise SystemExit(f"{BASELINE} is the reference and cannot be skipped")

    base = {"QUBITS": str(args.qubits), "ITERS": "1", "JCOUP": "1.0", "HFIELD": "1.0"}
    grid = itertools.product(
        args.models.split(","),
        args.steps.split(","),
        args.dt.split(","),
        args.seeds.split(","),
    )
    rows = []
    for model, steps, dt, seed in grid:
        env = {**base, "MODEL": model, "STEPS": steps, "DT": dt, "SEED": seed}
        # The reference depends on every axis, so it is recomputed per cell.
        ref, ref_obs = run(BASELINE, {**env, "ATOL": args.ref_atol}, args.quiet)
        print(
            f"{model} n={args.qubits} steps={steps} dt={dt} seed={seed}: reference is"
            f" {BASELINE} at atol={args.ref_atol} — {len(ref)} terms,"
            f" observable {ref_obs:.12g}",
            file=sys.stderr,
        )
        for atol in args.atols.split(","):
            for lib in libs:
                approx, obs = run(lib, {**env, "ATOL": atol}, args.quiet)
                stats = measure(approx, ref, float(atol))
                rows.append(
                    {
                        "model": model,
                        "library": lib,
                        "qubits": args.qubits,
                        "steps": steps,
                        "dt": dt,
                        "seed": seed,
                        "atol": atol,
                        "ref_atol": args.ref_atol,
                        "ref_terms": len(ref),
                        "ref_norm": f"{norm(ref):.6e}",
                        "observable": f"{obs:.15g}",
                        "ref_observable": f"{ref_obs:.15g}",
                        "obs_abs_err": f"{abs(obs - ref_obs):.6e}",
                        **{
                            k: f"{v:.6e}" if isinstance(v, float) else v
                            for k, v in stats.items()
                        },
                    }
                )
                print(
                    f"  atol={atol} {lib:>20}: {stats['terms']:>7} terms,"
                    f" L2 {stats['l2_err']:.3e}, obs err {rows[-1]['obs_abs_err']}",
                    file=sys.stderr,
                )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
