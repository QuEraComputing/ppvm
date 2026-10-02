"""Size sweep: Stim-generated Clifford memory circuits, timed on ppvm (this
branch and PR 204) and on official Stim, alternating tools over a few rounds.

Usage: `sweep.py [rounds]`. Writes `circuits/*.stim` and `sweep.csv`.
"""

import csv
import pathlib
import re
import statistics
import subprocess
import sys

import stim

HERE = pathlib.Path(__file__).resolve().parent
WORKTREES = HERE.parents[2]
PYTHON = "/Users/david/git/ppvm/ppvm-python/.venv/bin/python"
TOOLS = {
    "branch": [str(HERE / "target/release/stim-compare"), "file"],
    "pr204": [str(WORKTREES / "pr204-baseline/bench/stim-compare/target/release/stim-compare"), "file"],
    "stim": [PYTHON, str(HERE / "stim_bench.py"), "file"],
}
FAMILIES = [
    ("surface", "surface_code:rotated_memory_z", [3, 5, 7, 9, 11, 15, 19, 23, 27, 31]),
    ("repetition", "repetition_code:memory", [3, 9, 25, 75, 225, 675]),
    ("color", "color_code:memory_xyz", [3, 5, 7, 9, 13, 17, 21, 25, 31, 37, 43]),
]
CLIFFT_D7 = pathlib.Path.home() / "git/clifft-bench/workloads/circuits/pure_surface_d7_r7_p1e-3.stim"


def rewrite_c_xyz(text: str) -> str:
    """ppvm-stim has no `C_XYZ`; it equals `H` then `SQRT_X_DAG`."""

    def expand(m: re.Match) -> str:
        indent, targets = m.group(1), m.group(2)
        return f"{indent}H{targets}\n{indent}SQRT_X_DAG{targets}"

    return re.sub(r"^(\s*)C_XYZ(\s.*)$", expand, text, flags=re.M)


def circuits() -> list[tuple[str, int, pathlib.Path]]:
    out_dir = HERE / "circuits"
    out_dir.mkdir(exist_ok=True)
    out = [("clifft_surface", 7, CLIFFT_D7)]
    for family, task, distances in FAMILIES:
        for d in distances:
            c = stim.Circuit.generated(
                task,
                distance=d,
                rounds=d,
                after_clifford_depolarization=1e-3,
                before_round_data_depolarization=1e-3,
                before_measure_flip_probability=1e-3,
                after_reset_flip_probability=1e-3,
            )
            path = out_dir / f"{family}_d{d}.stim"
            path.write_text(rewrite_c_xyz(str(c)) + "\n")
            out.append((family, d, path))
    return out


def run(tool: str, path: pathlib.Path, shots: int) -> tuple[float, float, int, int]:
    """(execute median in µs, mean ones, qubits, measurements)."""
    out = subprocess.run(
        [*TOOLS[tool], str(path), str(shots)], capture_output=True, text=True, check=True
    ).stdout
    head = re.search(r"\[file\] (\d+) qubits, (\d+) measurements", out)
    m = re.search(r"execute median\s+([\d.]+)(ns|µs|ms|s)\b", out)
    scale = {"ns": 1e-3, "µs": 1.0, "ms": 1e3, "s": 1e6}[m.group(2)]
    ones = float(re.search(r"mean ones/shot ([\d.]+)", out).group(1))
    return float(m.group(1)) * scale, ones, int(head.group(1)), int(head.group(2))


def main() -> None:
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    rows = []
    for family, d, path in circuits():
        pilot, *_ = run("branch", path, 3)
        shots = max(10, min(500, int(1.0e6 / max(pilot, 1.0))))
        times = {tool: [] for tool in TOOLS}
        ones = {}
        for _ in range(rounds):
            for tool in TOOLS:
                t, o, n_qubits, n_meas = run(tool, path, shots)
                times[tool].append(t)
                ones[tool] = o
        med = {tool: statistics.median(ts) for tool, ts in times.items()}
        row = {
            "family": family,
            "d": d,
            "qubits": n_qubits,
            "measurements": n_meas,
            "shots": shots,
            **{f"{tool}_us": round(med[tool], 1) for tool in TOOLS},
            "branch_over_stim": round(med["branch"] / med["stim"], 3),
            "pr204_over_stim": round(med["pr204"] / med["stim"], 3),
            **{f"{tool}_ones": ones[tool] for tool in TOOLS},
        }
        rows.append(row)
        print(
            f"{family:14s} d={d:3d} n={n_qubits:5d} meas={n_meas:7d} shots={shots:3d}  "
            f"branch {med['branch']:10.1f}µs  pr204 {med['pr204']:10.1f}µs  stim {med['stim']:10.1f}µs  "
            f"branch/stim {row['branch_over_stim']:.2f}x  pr204/stim {row['pr204_over_stim']:.2f}x",
            flush=True,
        )
    with open(HERE / "sweep.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
