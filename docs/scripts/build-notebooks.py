#!/usr/bin/env python3
"""Convert Jupytext notebooks under ``docs/notebooks/`` into executed HTML
fragments that the Astro site can embed.

Pipeline per notebook (``foo.py``):

1. Read the Jupytext ``percent`` script and turn it into an in-memory
   ``ipynb`` via the Jupytext Python API.
2. Run every code cell via ``nbclient.NotebookClient`` so that text
   outputs, tracebacks, and matplotlib figures are captured inline. We
   force inline matplotlib so plots get embedded as base64 PNGs.
3. Render to an HTML fragment with ``nbconvert.HTMLExporter`` using the
   ``basic`` template — no page chrome, no toolbar, just the cell
   stream. The site's stylesheet (`global.css`) carries our own theme
   for the ``.jp-*`` classes.
4. Write ``docs/src/generated/notebooks/<name>.html`` and a sidecar
   ``<name>.json`` with metadata (title, ordered headings) that the
   Astro index page can render without parsing HTML.

Designed to be invoked from CI as the step before ``npx astro build``.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import jupytext
import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter

HERE = Path(__file__).resolve().parent
DOCS = HERE.parent
SOURCE_DIR = DOCS / "notebooks"
OUTPUT_DIR = DOCS / "src" / "generated" / "notebooks"

# Force a non-interactive matplotlib backend before any cell runs so
# inline ``plt.show()`` calls produce embedded PNGs rather than trying
# to open a window. We prepend this as the first code cell of every
# notebook with the ``ppvm-hidden-setup`` tag, execute the notebook,
# and then drop the tagged cell entirely via ``drop_hidden_setup_cells``
# before rendering — so neither the input nor any output of this cell
# survives into the HTML fragment that the site embeds.
MATPLOTLIB_SETUP = (
    "import matplotlib\n"
    "matplotlib.use('Agg')\n"
    "from matplotlib import pyplot as _ppvm_plt\n"
    "_ppvm_plt.rcParams['figure.dpi'] = 110\n"
)


def slug_for(path: Path) -> str:
    return path.stem.lower().replace("_", "-")


def extract_title_and_headings(nb: nbformat.NotebookNode) -> tuple[str, list[str]]:
    """Pull the document title (first H1) and the rest of the headings
    out of the markdown cells. We only walk markdown cells; output
    cells from executed code aren't part of the document structure.
    """
    title: str | None = None
    headings: list[str] = []
    pattern = re.compile(r"^(#{1,3})\s+(.+?)\s*$", re.MULTILINE)
    for cell in nb.cells:
        if cell.cell_type != "markdown":
            continue
        for match in pattern.finditer(cell.source):
            level, text = match.group(1), match.group(2).strip()
            if level == "#" and title is None:
                title = text
            else:
                headings.append(text)
    return title or "Untitled", headings


def prepend_setup_cell(nb: nbformat.NotebookNode) -> None:
    cell = nbformat.v4.new_code_cell(MATPLOTLIB_SETUP)
    cell["metadata"] = {"tags": ["ppvm-hidden-setup"]}
    nb.cells.insert(0, cell)


def drop_hidden_setup_cells(nb: nbformat.NotebookNode) -> None:
    nb.cells = [
        c
        for c in nb.cells
        if "ppvm-hidden-setup" not in (c.get("metadata", {}) or {}).get("tags", [])
    ]


def execute(nb: nbformat.NotebookNode, source_label: str) -> None:
    client = NotebookClient(
        nb,
        timeout=600,
        kernel_name="python3",
        # Plain text/HTML outputs only; we don't want stale `In [N]`
        # prompts cluttering the rendered fragment.
        allow_errors=False,
        record_timing=False,
        resources={"metadata": {"path": str(SOURCE_DIR)}},
    )
    try:
        client.execute()
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"notebook execution failed for {source_label}: {exc}\n")
        raise


_BODY_PATTERN = re.compile(r"<body[^>]*>(.*)</body>", re.DOTALL)
_MAIN_PATTERN = re.compile(r"</?main[^>]*>", re.IGNORECASE)
_SCRIPT_PATTERN = re.compile(r"<script\b[^>]*>.*?</script>", re.DOTALL | re.IGNORECASE)
_EVENT_HANDLER_PATTERN = re.compile(r"\s+on[a-z]+\s*=\s*(?:\"[^\"]*\"|'[^']*'|[^\s>]+)", re.IGNORECASE)


def sanitise(fragment: str) -> str:
    """Defense-in-depth scrub of the rendered notebook HTML before we
    embed it via ``set:html``. Notebook authors are trusted (the source
    .py files live in this repo and run during CI), but nbconvert's
    ``basic`` template will faithfully include any HTML a cell emits as
    rich output — including ``<script>`` and inline event handlers.
    Drop both, plus the ``<main>`` wrapper nbconvert adds (the page
    already has its own ``<main>`` from ``Base.astro``; nesting two is
    invalid HTML).
    """
    fragment = _SCRIPT_PATTERN.sub("", fragment)
    fragment = _EVENT_HANDLER_PATTERN.sub("", fragment)
    fragment = _MAIN_PATTERN.sub("", fragment)
    return fragment.strip()


def render_html(nb: nbformat.NotebookNode) -> str:
    """Render the notebook to an HTML *fragment* — just the cell stream,
    none of the surrounding JupyterLab stylesheet. The site's own CSS
    in ``global.css`` handles the visual treatment.
    """
    exporter = HTMLExporter()
    exporter.template_name = "basic"
    exporter.exclude_input_prompt = True
    exporter.exclude_output_prompt = True
    body, _ = exporter.from_notebook_node(nb)
    match = _BODY_PATTERN.search(body)
    inner = match.group(1) if match else body
    return sanitise(inner)


def detect_language(nb: nbformat.NotebookNode) -> str:
    kernelspec = (nb.metadata or {}).get("kernelspec", {}) or {}
    lang = kernelspec.get("language") or (
        (nb.metadata or {}).get("language_info", {}) or {}
    ).get("name")
    return (lang or "python").lower()


def build_one(source: Path) -> dict:
    sys.stderr.write(f"[notebooks] building {source.name}\n")
    nb = jupytext.read(source, fmt="py:percent")
    title, headings = extract_title_and_headings(nb)
    slug = slug_for(source)
    language = detect_language(nb)

    prepend_setup_cell(nb)
    execute(nb, source.name)
    drop_hidden_setup_cells(nb)

    html = render_html(nb)
    (OUTPUT_DIR / f"{slug}.html").write_text(html, encoding="utf-8")
    meta = {
        "slug": slug,
        "title": title,
        "headings": headings,
        "language": language,
        "source": f"docs/notebooks/{source.name}",
    }
    (OUTPUT_DIR / f"{slug}.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8"
    )
    return meta


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    sources = sorted(SOURCE_DIR.glob("*.py"))
    if not sources:
        sys.stderr.write(f"[notebooks] no .py files found under {SOURCE_DIR}\n")
        return 0
    index = [build_one(s) for s in sources]
    (OUTPUT_DIR / "index.json").write_text(
        json.dumps(index, indent=2), encoding="utf-8"
    )
    sys.stderr.write(f"[notebooks] built {len(index)} notebook(s)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
