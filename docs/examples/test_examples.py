"""
Documentation doctests.

Each example shipped in ``docs/examples/*.py`` is executed end-to-end as a
subprocess. The expected stdout, recorded below, must match exactly. If you
edit an example file (or the API it exercises) you must also update the
expected output here — that is the point: docs and code cannot drift.

Run from the repository root:

    uv run --project ppvm-python --group dev pytest docs/examples
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).resolve().parent

# (filename, expected_stdout) pairs. Expected output is matched exactly
# after stripping trailing whitespace from each line.
EXAMPLES: list[tuple[str, str]] = [
    (
        "paulisum_ghz.py",
        """\
1.000 * IZ
1.0
""",
    ),
    (
        "tableau_ghz.py",
        """\
qubit 0: 1, qubit 1: 1
correlated: True
""",
    ),
    (
        "stim_program.py",
        """\
single-run ok
sampled 16 shots
first shot: [<MeasurementResult.ONE: 1>, <MeasurementResult.ONE: 1>]
all shots correlated: True
""",
    ),
    (
        "stim_sampling.py",
        """\
(0, 0): 102
(1, 1): 98
correlated fraction: 1.0
""",
    ),
]


def _run(path: Path) -> str:
    result = subprocess.run(
        [sys.executable, str(path)],
        capture_output=True,
        text=True,
        check=False,
        cwd=path.parent,
    )
    assert result.returncode == 0, (
        f"example {path.name} failed (exit {result.returncode}):\n"
        f"--- stdout ---\n{result.stdout}\n"
        f"--- stderr ---\n{result.stderr}"
    )
    return result.stdout


def _normalize(text: str) -> str:
    return "\n".join(line.rstrip() for line in text.splitlines()).strip() + "\n"


@pytest.mark.parametrize(("filename", "expected"), EXAMPLES, ids=[e[0] for e in EXAMPLES])
def test_example_runs_and_matches_expected_output(filename: str, expected: str) -> None:
    path = EXAMPLES_DIR / filename
    assert path.exists(), f"missing example file: {path}"

    stdout = _run(path)

    assert _normalize(stdout) == _normalize(expected), (
        f"\nexample {filename} produced unexpected output.\n"
        f"--- expected ---\n{expected}"
        f"--- got ---\n{stdout}"
    )


def test_all_examples_are_covered() -> None:
    """Guard against forgetting to register a new example in EXAMPLES."""
    on_disk = {p.name for p in EXAMPLES_DIR.glob("*.py") if p.name not in {"__init__.py", "test_examples.py"}}
    registered = {name for name, _ in EXAMPLES}
    missing = on_disk - registered
    assert not missing, (
        f"example file(s) present on disk but not registered in EXAMPLES: {sorted(missing)}"
    )
