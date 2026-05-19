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

Executed outputs are also written to a content-addressed cache under
``docs/.notebook-cache/`` (or ``$PPVM_NOTEBOOK_CACHE_DIR``). The cache
key is ``sha256(CACHE_SCHEMA_VERSION + this script's source + notebook
source + Cargo.lock + Cargo.toml files + uv.lock)`` — so docs-only or
CSS-only PRs reuse already-executed notebooks, and any change to the
extractor itself (rendering, sanitiser, matplotlib setup) invalidates
every entry automatically. Bump ``CACHE_SCHEMA_VERSION`` if you need
to force a global invalidation for some other reason. GH Actions
persists this directory across runs via ``actions/cache``.

The fingerprint deliberately does not hash Rust ``.rs`` or Python
package sources, so a numerical change inside a crate that doesn't
touch a lockfile or ``Cargo.toml`` will not invalidate cached
outputs. Rely on ``cargo test --workspace`` / ``pytest`` to catch
those; if we add a scheduled full-rebuild workflow in the future, it
will serve as a second safety net.

Designed to be invoked from CI as the step before ``npx astro build``.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path

import bleach
import jupytext
import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter

HERE = Path(__file__).resolve().parent
DOCS = HERE.parent
REPO_ROOT = DOCS.parent
SOURCE_DIR = DOCS / "notebooks"
OUTPUT_DIR = DOCS / "src" / "generated" / "notebooks"
# Content-addressed cache for executed notebooks. Keyed by
# ``sha256(notebook source + shared runtime fingerprint)``; populated
# on every successful execution, read on subsequent runs. CI restores
# this directory via ``actions/cache`` so docs-only PRs reuse already-
# executed notebooks instead of re-running them. Override with
# ``PPVM_NOTEBOOK_CACHE_DIR=…`` (used by the GH Actions workflow to
# point at a stable cross-job location).
CACHE_DIR = Path(
    os.environ.get("PPVM_NOTEBOOK_CACHE_DIR", DOCS / ".notebook-cache")
).resolve()
# Set ``PPVM_NOTEBOOK_CACHE=0`` to force re-execution regardless of
# what's on disk (useful when investigating numerical drift).
CACHE_ENABLED = os.environ.get("PPVM_NOTEBOOK_CACHE", "1") != "0"
# Bump this to force-invalidate every cached notebook output. Useful
# when the cached artefact layout changes incompatibly (e.g. a new
# field in the sidecar JSON that downstream Astro pages depend on).
# Routine pipeline edits don't need a bump because the script's own
# source is part of the fingerprint via ``_shared_fingerprint_files``.
CACHE_SCHEMA_VERSION = "1"

# Switch matplotlib to the IPython "inline" backend before any cell
# runs so ``plt.show()`` triggers Jupyter's display hook and embeds
# the figure as a base64 PNG inside the cell's output. The plain
# ``Agg`` backend renders to a buffer but never reaches the
# notebook output stream, which is why an early version of this
# pipeline silently dropped every plot.
#
# We prepend this as the first code cell of every notebook with the
# ``ppvm-hidden-setup`` tag, execute the notebook, and then drop the
# tagged cell entirely via ``drop_hidden_setup_cells`` before
# rendering — so neither the input nor any output of this cell
# survives into the HTML fragment the site embeds.
MATPLOTLIB_SETUP = (
    "%matplotlib inline\n"
    "import matplotlib\n"
    "matplotlib.rcParams['figure.dpi'] = 110\n"
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

# Tag and attribute allow-lists for the rendered notebook HTML
# fragment. Notebook authors live in this repo and we trust their
# .py source, but nbconvert's ``basic`` template will faithfully
# emit any HTML a cell produces as rich output. We hand the
# fragment to ``bleach`` (html5lib-based, parses to a real DOM
# rather than string-matching) so a future notebook that prints
# an ``<iframe>`` or a ``javascript:`` URL can't smuggle active
# content into the page.
ALLOWED_TAGS = frozenset({
    # Block structure used by nbconvert's `basic` template
    "div", "section", "article", "header", "footer", "main",
    "p", "br", "hr",
    "h1", "h2", "h3", "h4", "h5", "h6",
    "ul", "ol", "li", "dl", "dt", "dd",
    "blockquote",
    "table", "thead", "tbody", "tfoot", "tr", "th", "td", "caption",
    # Inline prose
    "a", "em", "strong", "i", "b", "u", "s", "sub", "sup", "small", "mark", "code", "kbd", "samp",
    "span", "abbr", "cite", "q", "del", "ins",
    # Code rendering (Pygments wraps tokens in nested spans)
    "pre",
    # Images: matplotlib figures embed as `data:image/png;base64,...`
    "img",
    # Math from notebook LaTeX (rare, but `nbconvert` may pass through)
    "math", "mrow", "mi", "mn", "mo", "ms", "mtext", "mfrac", "msup", "msub",
    "msubsup", "msqrt", "mroot", "mfenced", "mtable", "mtr", "mtd",
    "annotation", "annotation-xml", "semantics",
})

ALLOWED_ATTRS = {
    "*": ["class", "id", "title", "aria-hidden", "aria-label", "role",
          "tabindex", "data-mime-type", "lang"],
    "a": ["href", "rel", "target"],
    "img": ["src", "alt", "width", "height"],
    "td": ["colspan", "rowspan", "headers", "scope"],
    "th": ["colspan", "rowspan", "headers", "scope"],
    "ol": ["start", "type"],
    "li": ["value"],
}

# Permitted URL schemes for `<a href>` and `<img src>`. `data:` is
# allowed so matplotlib figures (which arrive as
# `data:image/png;base64,...`) can survive — `bleach` enforces this
# at the attribute level.
ALLOWED_PROTOCOLS = frozenset({"http", "https", "mailto", "data"})


def sanitise(fragment: str) -> str:
    """Parse the rendered notebook HTML with bleach and re-serialise
    with only tags/attrs/URL schemes on the allow-list above.
    Anything else — `<script>`, `<iframe>`, inline event handlers,
    `javascript:` hrefs — is stripped at the DOM level, so we
    don't accidentally delete legitimate text that happens to
    mention `onload="…"` in prose or a code block.

    Drop nbconvert's `<main>` wrapper first because the surrounding
    page layout already provides its own `<main>` (one per page).
    """
    fragment = _MAIN_PATTERN.sub("", fragment)
    return bleach.clean(
        fragment,
        tags=ALLOWED_TAGS,
        attributes=ALLOWED_ATTRS,
        protocols=ALLOWED_PROTOCOLS,
        strip=True,
        strip_comments=False,
    ).strip()


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


# Files whose contents influence notebook outputs. Hashed once into
# ``_shared_fingerprint`` and combined with the notebook source to
# form each notebook's cache key.
#
# Includes this script itself so changes to the rendering / sanitiser
# / matplotlib-setup path invalidate cached outputs automatically —
# without that, a fix to (say) the bleach allow-list would silently
# keep serving the previous (potentially insecure or broken) HTML
# for every unchanged notebook source.
#
# We intentionally do NOT hash Rust ``.rs`` or Python package
# sources — that would blow up the fingerprint on cosmetic changes.
# Cargo.toml + Cargo.lock + uv.lock cover dependency bumps and
# version changes; ``cargo test`` / ``pytest`` catch the rest.
def _shared_fingerprint_files() -> list[Path]:
    files: list[Path] = [Path(__file__).resolve()]
    for name in ("Cargo.lock", "Cargo.toml"):
        p = REPO_ROOT / name
        if p.exists():
            files.append(p)
    uv_lock = REPO_ROOT / "ppvm-python" / "uv.lock"
    if uv_lock.exists():
        files.append(uv_lock)
    files.extend(sorted((REPO_ROOT / "crates").glob("*/Cargo.toml")))
    return files


def _compute_shared_fingerprint() -> bytes:
    h = hashlib.sha256()
    h.update(b"schema=")
    h.update(CACHE_SCHEMA_VERSION.encode("utf-8"))
    h.update(b"\0")
    for f in _shared_fingerprint_files():
        try:
            label = f.relative_to(REPO_ROOT).as_posix()
        except ValueError:
            label = f.name
        h.update(label.encode("utf-8"))
        h.update(b"\0")
        h.update(f.read_bytes())
        h.update(b"\0")
    return h.digest()


_shared_fingerprint: bytes | None = None


def shared_fingerprint() -> bytes:
    global _shared_fingerprint
    if _shared_fingerprint is None:
        _shared_fingerprint = _compute_shared_fingerprint()
    return _shared_fingerprint


def notebook_cache_key(source: Path) -> str:
    h = hashlib.sha256()
    h.update(shared_fingerprint())
    h.update(b"\0")
    h.update(source.read_bytes())
    return h.hexdigest()


def _try_restore_from_cache(source: Path, slug: str) -> dict | None:
    if not CACHE_ENABLED:
        return None
    key = notebook_cache_key(source)
    html_cached = CACHE_DIR / f"{key}.html"
    meta_cached = CACHE_DIR / f"{key}.json"
    if not (html_cached.exists() and meta_cached.exists()):
        return None
    try:
        meta = json.loads(meta_cached.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    # Slug may have changed since the cache entry was written (e.g. file
    # rename). Trust the current slug; rewrite meta and copy under the
    # current output name.
    meta["slug"] = slug
    shutil.copyfile(html_cached, OUTPUT_DIR / f"{slug}.html")
    (OUTPUT_DIR / f"{slug}.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8"
    )
    sys.stderr.write(f"[notebooks] cache hit  {source.name} ({key[:12]})\n")
    return meta


def _write_cache(source: Path, html: str, meta: dict) -> None:
    if not CACHE_ENABLED:
        return
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    key = notebook_cache_key(source)
    (CACHE_DIR / f"{key}.html").write_text(html, encoding="utf-8")
    (CACHE_DIR / f"{key}.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8"
    )


def build_one(source: Path) -> dict:
    slug = slug_for(source)
    cached = _try_restore_from_cache(source, slug)
    if cached is not None:
        return cached

    sys.stderr.write(f"[notebooks] executing {source.name}\n")
    nb = jupytext.read(source, fmt="py:percent")
    title, headings = extract_title_and_headings(nb)
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
    _write_cache(source, html, meta)
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
