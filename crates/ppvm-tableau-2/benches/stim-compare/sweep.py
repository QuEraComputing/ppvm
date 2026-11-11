"""Size sweep: Stim-generated Clifford memory circuits (see `circuits.py`) on
ppvm-tableau-2's `run_stim` and on reference binaries, alternating tools over a
few rounds; median per tool. Writes the circuits and `sweep.csv` to
`target/stim-compare/`.

    uv run --no-project --with stim python sweep.py --reference stim=PATH
        [--case surface:7]... [--family surface] [--max-qubits N] [--rounds 3]

Python Stim only generates the circuits; timing is against the references
(e.g. `stim_driver`, see README.md). Times are per shot, noise-free,
gates/measurements/resets only; see the design doc's Benchmarks section for how
instruction batching differs between the tools.
"""

import argparse
import csv
import statistics

import circuits
import common


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--max-qubits", type=int, default=10**9)
    parser.add_argument("--family", action="append")
    parser.add_argument(
        "--case", action="append", help="family:distance, e.g. surface:7"
    )
    parser.add_argument("--reference", action="append", default=[])
    parser.add_argument(
        "--seconds", type=float, default=0.3, help="target runtime per measurement"
    )
    args = parser.parse_args()

    common.build_ppvm()
    tools = common.tools(args.reference)
    reference = list(tools)[1]
    cases = {(f, int(d)) for f, d in (c.split(":") for c in args.case or [])}
    (common.OUT / "circuits").mkdir(parents=True, exist_ok=True)
    rows = []
    for family, task, distances in circuits.FAMILIES:
        if args.family and family not in args.family:
            continue
        for d in distances:
            if cases and (family, d) not in cases:
                continue
            text = circuits.generate(task, d)
            n = circuits.num_qubits(text)
            if n > args.max_qubits:
                continue
            path = common.OUT / "circuits" / f"{family}_d{d}.stim"
            path.write_text(text)
            shots = {
                name: common.calibrated_shots(cmd, path, args.seconds)
                for name, cmd in tools.items()
            }
            times = {name: [] for name in tools}
            for _ in range(args.rounds):
                for name, cmd in tools.items():
                    times[name].append(common.time_per_shot(cmd, path, shots[name]))
            row = {"family": family, "distance": d, "qubits": n}
            row |= {
                f"{name}_us": round(statistics.median(ts), 3)
                for name, ts in times.items()
            }
            row[f"ppvm_over_{reference}"] = round(
                row["ppvm_us"] / row[f"{reference}_us"], 3
            )
            rows.append(row)
            cells = "  ".join(f"{name} {row[f'{name}_us']:10.1f} µs" for name in tools)
            ratio = row[f"ppvm_over_{reference}"]
            print(
                f"{family:10} d={d:<4} {n:5} q  {cells}  ppvm/{reference} {ratio:.2f}",
                flush=True,
            )
    with open(common.OUT / "sweep.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {common.OUT / 'sweep.csv'}")


if __name__ == "__main__":
    main()
