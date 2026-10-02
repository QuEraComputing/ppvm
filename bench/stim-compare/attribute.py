"""Per-commit attribution: time harness binaries built at each commit (see
README) against official Stim on a fixed workload set, alternating over rounds.

Usage: `attribute.py <bin-dir> <commit>... [--rounds N]`. Writes `attribution.csv`.
"""

import argparse
import csv
import pathlib
import re
import statistics
import subprocess

HERE = pathlib.Path(__file__).resolve().parent
PYTHON = "/Users/david/git/ppvm/ppvm-python/.venv/bin/python"
# (label, harness args after the binary, stim_bench.py args)
WORKLOADS = [
    ("surface_d30 full", ["20", "ppvm", "full"], ["20", "full"]),
    ("surface_d30 1 round", ["20", "ppvm", "full-r1"], ["20", "full-r1"]),
    ("surface_d30 no measure", ["20", "ppvm", "no-measure"], ["20", "no-measure"]),
    *(
        (name, ["file", str(HERE / f"circuits/{name}.stim"), shots], ["file", str(HERE / f"circuits/{name}.stim"), shots])
        for name, shots in [
            ("surface_d7", "500"),
            ("surface_d11", "500"),
            ("surface_d19", "200"),
            ("repetition_d75", "500"),
            ("color_d31", "20"),
        ]
    ),
]


def execute_us(out: str) -> float:
    m = re.search(r"execute median\s+([\d.]+)(ns|µs|ms|s)\b", out)
    return float(m.group(1)) * {"ns": 1e-3, "µs": 1.0, "ms": 1e3, "s": 1e6}[m.group(2)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("bin_dir")
    ap.add_argument("commits", nargs="+")
    ap.add_argument("--rounds", type=int, default=3)
    args = ap.parse_args()
    tools = {c: [str(pathlib.Path(args.bin_dir) / c)] for c in args.commits}
    tools["stim"] = [PYTHON, str(HERE / "stim_bench.py")]

    rows = []
    for label, harness, stim_args in WORKLOADS:
        times = {t: [] for t in tools}
        for _ in range(args.rounds):
            for tool, cmd in tools.items():
                extra = stim_args if tool == "stim" else harness
                out = subprocess.run([*cmd, *extra], capture_output=True, text=True, check=True).stdout
                times[tool].append(execute_us(out))
        med = {t: statistics.median(v) for t, v in times.items()}
        rows.append({"workload": label, **{f"{t}_us": round(v, 1) for t, v in med.items()}})
        print(label, "  ".join(f"{t}={v:.1f}" for t, v in med.items()), flush=True)
    with open(HERE / "attribution.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
