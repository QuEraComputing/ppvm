"""Regression check on non-Clifford programs (ppvm only; Stim can't run them):
time harness binaries built at two commits, alternating over rounds.

Usage: `nonclifford.py <bin-dir> <commit>... [--rounds N]`. Writes `nonclifford.csv`.
"""

import argparse
import csv
import pathlib
import re
import statistics
import subprocess

HERE = pathlib.Path(__file__).resolve().parent
CLIFFT = pathlib.Path.home() / "git/clifft-bench/workloads/circuits"
PROGRAMS = [
    HERE / "../../crates/ppvm-stim/tests/data/cultivation_d5.stim",
    *(CLIFFT / f"{name}.stim" for name in [
        "msc_d3_inject_cultivate_p1e-3",
        "msc_d5_inject_cultivate_p1e-3",
        "distillation",
        "coherent_d3_r1",
        "coherent_d3_r3",
        "quantum_volume_q10_seed42",
    ]),
]
ANNOTATIONS = {"DETECTOR", "OBSERVABLE_INCLUDE", "QUBIT_COORDS", "SHIFT_COORDS", "TICK", "REPEAT"}


def to_ppvm_dialect(text: str) -> str:
    """Clifft's `R_X(a)` / `U3(t, p, l)` (half-turns) as ppvm's `I[...]` tags (radians)."""

    def rot(m: re.Match) -> str:
        return f"I[R_{m.group(1)}(theta={m.group(2).strip()}*pi)]"

    def u3(m: re.Match) -> str:
        t, p, l = (a.strip() for a in m.group(1).split(","))
        return f"I[U3(theta={t}*pi, phi={p}*pi, lambda={l}*pi)]"

    text = re.sub(r"\bR_([XYZ])\(([^)]*)\)", rot, text)
    return re.sub(r"\bU3\(([^)]*)\)", u3, text)


def n_qubits(text: str) -> int:
    """Largest qubit target plus one, without a Stim parse (Stim rejects `T`)."""
    top = 0
    for line in text.splitlines():
        line = line.split("#")[0].strip()
        name = re.split(r"[\s(]", line, maxsplit=1)[0]
        if not line or name in ANNOTATIONS or line == "}":
            continue
        targets = re.sub(r"\([^)]*\)", " ", line[len(name):])
        for token in re.split(r"[\s*]+", targets):
            token = token.lstrip("!").lstrip("XYZ")
            if token.isdigit():
                top = max(top, int(token) + 1)
    return max(top, 1)


def execute_us(out: str) -> float:
    m = re.search(r"execute median\s+([\d.]+)(ns|µs|ms|s)\b", out)
    return float(m.group(1)) * {"ns": 1e-3, "µs": 1.0, "ms": 1e3, "s": 1e6}[m.group(2)]


def run(binary: pathlib.Path, program: pathlib.Path, shots: int, n: int) -> float:
    out = subprocess.run(
        [str(binary), "file", str(program), str(shots), "ppvm", str(n)],
        capture_output=True, text=True, check=True,
    ).stdout
    return execute_us(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("bin_dir")
    ap.add_argument("commits", nargs="+")
    ap.add_argument("--rounds", type=int, default=3)
    args = ap.parse_args()
    bins = {c: pathlib.Path(args.bin_dir) / c for c in args.commits}

    out_dir = HERE / "circuits/nonclifford"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for source in PROGRAMS:
        program = out_dir / source.name
        program.write_text(to_ppvm_dialect(source.read_text()))
        n = n_qubits(program.read_text())
        pilot = run(bins[args.commits[-1]], program, 2, n)
        shots = max(3, min(500, int(1.0e6 / max(pilot, 1.0))))
        times = {c: [] for c in bins}
        for _ in range(args.rounds):
            for c, b in bins.items():
                times[c].append(run(b, program, shots, n))
        med = {c: statistics.median(v) for c, v in times.items()}
        first, last = args.commits[0], args.commits[-1]
        row = {"program": program.stem, "qubits": n, "shots": shots,
               **{f"{c}_us": round(v, 1) for c, v in med.items()},
               "ratio": round(med[last] / med[first], 3)}
        rows.append(row)
        print(f"{program.stem:34s} n={n:3d} shots={shots:3d}  "
              + "  ".join(f"{c} {v:10.1f}µs" for c, v in med.items())
              + f"  {last}/{first} {row['ratio']:.2f}x", flush=True)
    with open(HERE / "nonclifford.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
