"""Running the timing tools: ppvm's `run_stim` and reference binaries with
the same CLI (`<tool> <circuit> <shots>`, printing `<x> µs/shot` or
`<x> us/shot`), such as `stim_driver` built against a Stim checkout."""

import pathlib
import re
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[3]
OUT = REPO / "target" / "stim-compare"
RUN_STIM = REPO / "target" / "release" / "examples" / "run_stim"


def build_ppvm() -> None:
    subprocess.run(
        [
            "cargo",
            "build",
            "--release",
            "-q",
            "-p",
            "ppvm-tableau-2",
            "--example",
            "run_stim",
        ],
        cwd=REPO,
        check=True,
    )


def tools(references: list[str]) -> dict[str, list[str]]:
    """`ppvm`, then the `NAME=PATH` references in order."""
    if not references:
        sys.exit("give at least one --reference NAME=PATH (see README.md)")
    found = {"ppvm": [str(RUN_STIM)]}
    for reference in references:
        name, path = reference.split("=", 1)
        found[name] = [path]
    return found


def time_per_shot(cmd: list[str], circuit: pathlib.Path, shots: int) -> float:
    """µs per shot."""
    out = subprocess.run(
        [*cmd, str(circuit), str(shots)], capture_output=True, text=True, check=True
    ).stdout
    match = re.search(r"([\d.]+) (?:µs|us)/shot", out)
    if not match:
        raise RuntimeError(f"no timing in output of {cmd}: {out!r}")
    return float(match.group(1))


def calibrated_shots(cmd: list[str], circuit: pathlib.Path, seconds: float) -> int:
    """Shots for about `seconds` of runtime."""
    estimate = time_per_shot(cmd, circuit, 3)
    return max(3, min(50_000, int(seconds * 1e6 / max(estimate, 1e-3))))
