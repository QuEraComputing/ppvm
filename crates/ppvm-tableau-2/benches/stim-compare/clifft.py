"""clifft-bench workloads, noise included, on ppvm-tableau-2's `run_clifft` and
on reference binaries with the same CLI. Follows clifft-bench's probe contract
loosely: per tool and workload, a calibration run (also the warmup) picks the
shots for about `--seconds`, then `--probes` timed runs alternate over the
tools; the median is reported. A probe is at least one shot, so on slow
workloads it runs longer than `--seconds`. Writes `clifft.csv` to
`target/stim-compare/`.

    python clifft.py --reference NAME=PATH... [--workloads DIR] [--case NAME]...
        [--skip NAME]... [--seconds 1] [--probes 3]

Workloads are clifft-bench's `.stim` files (`scripts/fetch_clifft_bench.sh`
fetches them into `target/clifft-bench/`). Every tool reads the original file:
the ppvm runners rewrite Clifft's rotation and `U3` lines themselves.

Stim merges consecutive instructions with the same gate and arguments, and
shares setup across an instruction's targets; ppvm's `stim-parser` does not, so
ppvm runs one instruction per line. See the design doc's Benchmarks section.
"""

import argparse
import csv
import pathlib
import statistics

import common


def calibrated_shots(cmd: list[str], circuit: pathlib.Path, seconds: float) -> int:
    """Shots for about `seconds`: one shot, then a ~10% run if that is fast."""
    estimate = common.time_per_shot(cmd, circuit, 1)
    if estimate * 1e-6 < seconds / 20:
        shots = max(1, int(seconds / 10 * 1e6 / estimate))
        estimate = common.time_per_shot(cmd, circuit, shots)
    return max(1, int(seconds * 1e6 / estimate))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", action="append", default=[])
    parser.add_argument(
        "--workloads", type=pathlib.Path, default=common.REPO / "target" / "clifft-bench"
    )
    parser.add_argument("--case", action="append", help="workload name, e.g. coherent_d3_r1")
    parser.add_argument("--skip", action="append", default=[], help="workload name to leave out")
    parser.add_argument("--seconds", type=float, default=1.0, help="target time per probe")
    parser.add_argument("--probes", type=int, default=3)
    args = parser.parse_args()

    common.build_ppvm("run_clifft")
    tools = common.tools(args.reference, "run_clifft")
    references = list(tools)[1:]
    paths = sorted(args.workloads.glob("*.stim"))
    if args.case:
        paths = [p for p in paths if p.stem in args.case]
    paths = [p for p in paths if p.stem not in args.skip]
    if not paths:
        raise SystemExit(f"no workloads in {args.workloads}")
    rows = []
    for path in paths:
        shots = {name: calibrated_shots(cmd, path, args.seconds) for name, cmd in tools.items()}
        times = {name: [] for name in tools}
        for _ in range(args.probes):
            for name, cmd in tools.items():
                times[name].append(common.time_per_shot(cmd, path, shots[name]))
        row = {"workload": path.stem}
        row |= {f"{name}_us": round(statistics.median(ts), 3) for name, ts in times.items()}
        row |= {f"{name}_shots": shots[name] for name in tools}
        for ref in references:
            row[f"ppvm_over_{ref}"] = round(row["ppvm_us"] / row[f"{ref}_us"], 3)
        rows.append(row)
        cells = "  ".join(f"{name} {row[f'{name}_us']:12.1f} µs" for name in tools)
        ratios = "  ".join(f"ppvm/{ref} {row[f'ppvm_over_{ref}']:.2f}" for ref in references)
        print(f"{path.stem:32} {cells}  {ratios}", flush=True)
    common.OUT.mkdir(parents=True, exist_ok=True)
    with open(common.OUT / "clifft.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {common.OUT / 'clifft.csv'}")


if __name__ == "__main__":
    main()
