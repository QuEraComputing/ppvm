"""Size sweep: Stim-generated Clifford memory circuits (see `circuits.py`) on
ppvm-tableau-2's `run_stim` and on official Stim, alternating tools over a few
rounds; median per tool. Writes the circuits and `sweep.csv` to
`target/stim-compare/`.

    uv run --no-project --with stim python sweep.py [--rounds 3] [--max-qubits N]
        [--family surface] [--extra NAME=PATH]...

`--extra` adds a reference binary with `run_stim`'s CLI. Times are per shot,
noise-free, gates/measurements/resets only; see the design doc's Benchmarks
section for how instruction batching differs between the tools.
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
    parser.add_argument("--extra", action="append", default=[])
    parser.add_argument(
        "--seconds", type=float, default=0.3, help="target runtime per measurement"
    )
    args = parser.parse_args()

    common.build_ppvm()
    tools = common.tools(args.extra)
    (common.OUT / "circuits").mkdir(parents=True, exist_ok=True)
    rows = []
    for family, task, distances in circuits.FAMILIES:
        if args.family and family not in args.family:
            continue
        for d in distances:
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
            row["ppvm_over_stim"] = round(row["ppvm_us"] / row["stim_us"], 3)
            rows.append(row)
            others = "  ".join(
                f"{name} {row[f'{name}_us']:.1f}"
                for name in tools
                if name not in ("ppvm", "stim")
            )
            print(
                f"{family:10} d={d:<4} {n:5} q  ppvm {row['ppvm_us']:10.1f} µs  stim {row['stim_us']:10.1f} µs"
                f"  ppvm/stim {row['ppvm_over_stim']:.2f}  {others}",
                flush=True,
            )
    with open(common.OUT / "sweep.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {common.OUT / 'sweep.csv'}")


if __name__ == "__main__":
    main()
