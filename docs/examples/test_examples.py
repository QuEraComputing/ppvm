"""
Documentation doctests.

Each example shipped in ``docs/examples/*.py`` (and the runnable Python
snippets shipped with the ``ppvm-usage`` skill) is executed end-to-end as
a subprocess. Expected output lives **inside the example file**, as
``# → <line>`` comments next to the ``print`` statements (or in a small
trailing block when one ``print`` is called from a loop). The test
extracts those markers in source order, joins them, and compares against
the captured stdout.

Format:

    print(f"qubit 0: {r0}")  # → qubit 0: 1
    for x in xs:
        print(x)
    # → first
    # → second

Example files are auto-discovered from the examples directories via
``_discover()``; there is no separate registration step. A per-file
precondition check still fails if a file has zero ``# → `` markers (we
don't want examples that silently pass without verifying any output).

Run from the repository root:

    uv run --project ppvm-python --group dev pytest docs/examples
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = Path(__file__).resolve().parent
SKILL_DIR = (REPO_ROOT / "skills/ppvm-usage/examples/python").resolve()

# Match any line containing "# → <text>" — the marker is the arrow
# (U+2192), preceded by "# " and an optional single space afterwards.
_MARKER = re.compile(r"#\s*→\s?(.*?)\s*$")


def _expected_from_source(path: Path) -> str:
    lines: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        m = _MARKER.search(raw)
        if m is not None:
            lines.append(m.group(1))
    return ("\n".join(lines) + "\n") if lines else ""


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


def _discover(directory: Path) -> list[Path]:
    skip = {"__init__.py", "test_examples.py"}
    return sorted(
        p
        for p in directory.glob("*.py")
        if p.name not in skip and "__pycache__" not in p.parts
    )


DOCS_EXAMPLES = _discover(EXAMPLES_DIR)
SKILL_EXAMPLES = _discover(SKILL_DIR)


@pytest.mark.parametrize("path", DOCS_EXAMPLES, ids=[p.name for p in DOCS_EXAMPLES])
def test_docs_example(path: Path) -> None:
    expected = _expected_from_source(path)
    assert expected, (
        f"{path.relative_to(REPO_ROOT)} has no `# → ` expected-output markers; "
        "add at least one next to a print() so CI verifies behavior."
    )
    stdout = _run(path)
    assert _normalize(stdout) == _normalize(expected), (
        f"\nexample {path.name} produced unexpected output.\n"
        f"--- expected (from `# → ` markers) ---\n{expected}"
        f"--- got ---\n{stdout}"
    )


@pytest.mark.parametrize(
    "path", SKILL_EXAMPLES, ids=[f"skill/{p.name}" for p in SKILL_EXAMPLES]
)
def test_skill_example(path: Path) -> None:
    expected = _expected_from_source(path)
    assert expected, (
        f"{path.relative_to(REPO_ROOT)} has no `# → ` expected-output markers; "
        "add at least one next to a print() so CI verifies behavior."
    )
    stdout = _run(path)
    assert _normalize(stdout) == _normalize(expected), (
        f"\nskill example {path.name} produced unexpected output.\n"
        f"--- expected (from `# → ` markers) ---\n{expected}"
        f"--- got ---\n{stdout}"
    )


def test_examples_directories_are_nonempty() -> None:
    """Guard against an accidentally empty examples directory hiding the matrix."""
    assert DOCS_EXAMPLES, f"no example files discovered under {EXAMPLES_DIR}"
    assert SKILL_EXAMPLES, f"no skill example files discovered under {SKILL_DIR}"
