# AGENTS.md

> **Read the Developer Guide first.** The canonical contributor reference for
> this repository — for both human and AI contributors — is the Astro page
> at [`docs/src/pages/develop.astro`](docs/src/pages/develop.astro) (rendered
> at `/develop/` on the live site). It covers project layout, build/test
> commands for the Rust workspace, the Python package, **and the docs site
> itself** (`§ 2 Build & test`), architecture, conventions, the Python
> binding pipeline, extension recipes, and the file-by-file "where to look
> for X" table.
>
> This file is a short pointer so agents can find that guide quickly. Do
> not add content here that belongs in the Developer Guide — keep the
> guide as the single source of truth.

## Install the ppvm-usage skill

If you have [ion](https://ion.rogerluo.dev) installed, install the
`ppvm-usage` skill before writing any ppvm code:

```bash
ion add QuEraComputing/ppvm/skills/ppvm-usage
```

The skill (`skills/ppvm-usage/SKILL.md` in this repo) covers the Heisenberg /
Schrödinger gate-order trap, `Config`-generic `PauliSum` usage, truncation
strategies, and Python / Rust call sites for both backends. Read it before
the Developer Guide if your task is *using* ppvm rather than modifying its
internals.

## TL;DR for agents

If you are an AI agent picking up a task in this repository:

1. Open [`docs/src/pages/develop.astro`](docs/src/pages/develop.astro) and
   read the sections relevant to your task. The "For AI agents" callout at
   the top tells you which sections are load-bearing.
2. Use `uv` for anything Python; never `pip`.
3. Use Conventional Commits: `<type>(<scope>): <description>`.
4. Build & test the relevant target:
   - **Rust workspace**: `cargo test --workspace`
   - **Python package**: `uv run --project ppvm-python --group dev pytest …`
   - **Docs site**: `cd docs && npm run build` (chains `extract:rust`,
     `extract:python`, `extract:notebooks`, then `astro build`).
     Per-step commands and the `astro:dev` / `astro:build` escape
     hatches are documented under `§ 2 → "This docs site"` in the
     Developer Guide and in [`docs/README.md`](docs/README.md).
   - **Docs notebooks (`docs/notebooks/*.py`)**: executed at build time
     by [`docs/scripts/build-notebooks.py`](docs/scripts/build-notebooks.py)
     and embedded into `/examples/<slug>`. Outputs are content-addressed
     cached under `docs/.notebook-cache/`, fingerprinted by
     `CACHE_SCHEMA_VERSION` + the extractor script itself + the
     notebook source + every `Cargo.toml` + `Cargo.lock` +
     `ppvm-python/uv.lock`. (The extractor is in the fingerprint so
     rendering/sanitiser edits invalidate cached outputs automatically.)
     CI persists the directory via `actions/cache` (see the "Restore
     executed-notebook cache" step in
     [`.github/workflows/docs.yml`](.github/workflows/docs.yml)).
     **If you change the fingerprint inputs**, update both
     `_shared_fingerprint_files()` in the script *and* the
     `hashFiles(...)` call in the workflow — they must stay in sync
     (both already list `docs/scripts/build-notebooks.py`). Force a
     clean re-execution with `PPVM_NOTEBOOK_CACHE=0`, or bump
     `CACHE_SCHEMA_VERSION` for a global invalidation that survives
     in the cache. Full rationale lives under `§ 2.1 Notebook
     execution & caching` in the Developer Guide.
5. Respect the `Config`-trait generics in `ppvm-traits`; do not introduce
   runtime dispatch where a compile-time bound suffices.
6. Pauli propagation runs **backwards** (Heisenberg picture). Reverse the
   gate order accordingly when writing tests.
7. **Python docstring cross-references** use Markdown backtick spans, not
   RST syntax. The docs pipeline (griffe → `marked`) renders docstrings as
   Markdown, so `:meth:` / `:func:` / `:class:` are never parsed as links —
   they appear as literal text. Use plain backticks instead:
   - ✅ `` `fork` `` or `` `GeneralizedTableau.sample` ``
   - ❌ `` :meth:`fork` `` or `` :func:`ppvm.sample_stim` ``

## Workspace at a glance

```
crates/ppvm-traits          # Trait system, Config bundle, Pauli alphabet, map impls
crates/ppvm-pauli-word      # Packed Pauli strings: PauliWord, phased, lossy, pattern
crates/ppvm-pauli-sum       # PauliSum engine, truncation strategy, concrete configs
crates/ppvm-tableau         # Stabilizer + generalized-tableau simulator
crates/ppvm-sym             # Symbolic (parametric) Pauli propagation
crates/ppvm-stim            # Stim program execution against the tableau
crates/stim-parser          # Standalone Stim parser
crates/ppvm-python-native   # PyO3 cdylib (maturin)
ppvm-python/                # Pure-Python wrapper (uv_build)
docs/                       # Astro docs site — includes the Developer Guide
```

Everything else — design patterns, extension recipes, the file-by-file
"where to look for X" table — is in the Developer Guide.

## Docs-site visual language

If you're modifying anything under `docs/`, respect these conventions. The
site is meant to read as research-grade documentation, not a marketing
landing — that intent guides every other call.

### Design vocabulary

**Swiss / flat-academic.** Modern neo-grotesque sans for UI and body
(`Inter`), an editorial serif kept in reserve for italic captions and
sidenotes (`Source Serif 4`), monospace for code (`JetBrains Mono`).
Hierarchy comes from size, weight (300 / 400 / 500), and whitespace —
not from decorative dividers. No boxed widgets, no drop shadows, no
gradient text fills. Hairline rules (`var(--rule)`) separate sections;
they never enclose them.

**One brand accent.** The page is monochrome except for `var(--brand)`,
which carries interactive affordances only — links, active nav, the
primary CTA, the active install tab's underline. Brand purple is
`#5b3fb8` on the off-white canvas and lifts to `#c2a8ff` on the dark
canvas for AA contrast.

**Two canvases, same gradient.** Light: `#fafaf7` (off-white, picks up
the Bloqade-docs precedent). Dark: `#0c0820` (the deep indigo
quera.com itself uses). Both modes share the same QuEra "dawn"
gradient stops (`#7f3eff` → `#fa82ec` → `#ff7c24`) — those colors are
saturated and read against either canvas.

**Text-link buttons.** No filled chrome. A `.btn` is a sans label with
a hairline underline (`var(--rule)`) and an animated `→` arrow that
slides on hover. Primary CTAs use `var(--accent)` for the underline;
secondary CTAs darken to `var(--ink)`.

### The atom-cloud signature

The landing hero (`docs/src/pages/index.astro`) has an inline SVG of
~11 overlapping spheres, all filled with the QuEra dawn linear
gradient, blurred (stdDeviation=14) and masked to the upper-right
corner so the reading column on the left stays clean canvas.

**Don't import the painterly hero `.webp`** that quera.com uses. The
SVG cloud is fully styleable, scales without artefacts, and stays
crisp in dark mode — the dawn stops are the same in both themes.

**Reserve the cloud for the landing page only.** Every other page
(`/quickstart`, `/develop`, `/api`) is strict monochrome. The cloud is
a signal, not wallpaper; repeating it cheapens it.

### Code blocks and syntax highlighting

`highlight.js` runs client-side via a CDN script in `Base.astro`. The
syntax theme is defined in `global.css` using the same design tokens
as the rest of the site:

- Keywords / decorators / TOML headers → `var(--brand)` (purple).
- Strings → muted dawn-orange (`#8c4a18` light, `#ffa56b` dark).
- Numbers → muted dawn-violet (`#6a3aa5` light, `#c8a8ff` dark).
- Comments → `var(--ink-faint)`, italic.
- Function / type names → plain ink, no extra color.

Five distinct tones — not a full IDE rainbow. **Highlighting is
opt-in**: only `<pre><code class="language-…">` blocks get themed.
ASCII project trees, agent prompts, and commit-message examples are
left without a language hint, so they read as plain monospace.

When adding new code blocks: tag with `language-python`, `language-rust`,
`language-toml`, or `language-bash`. The `<ExampleBlock>` component
infers the language from the file extension automatically.

### Things to keep doing

- Cross-reference Rust symbols in prose to their `/api/` anchor (see
  the `link.X` helper in `docs/src/pages/develop.astro`).
- Use the `<Base toc={…}>` prop on long pages — the layout renders a
  sticky left sidebar with scroll-spy. The landing page intentionally
  has no TOC (it isn't a long-read).
- Keep the install-tabs widget on the landing page synchronized:
  Python / Rust / Agents, in that order; Python active by default.

### Things to avoid

- Adding a new accent color. If you reach for orange or pink, that's a
  signal to use the dawn gradient (only in the atom cloud) or do
  without color entirely.
- Boxing widgets in a border-radius card. Hairline rules are the
  separator. The hero code card is the only "card" on the page.
- Repeating the atom cloud, the dawn gradient stamp, or any
  brand-coloured rule beyond the hero.
- Importing additional fonts. The three families above carry the
  whole site.
- `pip` in any code or doc example. Project policy is `uv` everywhere.

When in doubt, open one of the existing pages in your browser at
`http://127.0.0.1:4322/` (run `npx astro dev` in `docs/`) and match
the tone of that page rather than guessing.
