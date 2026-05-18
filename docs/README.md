# ppvm docs site

Astro site for the ppvm documentation. Deployed to gh-pages root by
`.github/workflows/docs.yml`; PRs that touch this directory get a
per-PR preview deploy.

## Quick reference

From this directory:

```bash
npm install                 # one-time
npm run dev                 # extract + astro dev → http://127.0.0.1:4321
npm run build               # extract + astro build → ./dist/
```

`npm run dev` and `npm run build` run every prerequisite extraction
step first, so a fresh clone has one command to remember. Each
extraction step is also available on its own when you're iterating on
a specific layer:

| Command                       | What it does |
|-------------------------------|-----|
| `npm run extract:rust`        | Re-run `cargo +nightly rustdoc -- --output-format json` for the three Rust crates and reshape into `src/data/rust-api.json`. Source: `scripts/extract-rust.mjs`. |
| `npm run extract:python`      | Re-run `griffe dump ppvm` and reshape into `src/data/python-api.json`. Source: `scripts/extract-python.mjs`. |
| `npm run extract:notebooks`   | Execute every Jupytext `.py` under `notebooks/` and emit HTML fragments + metadata to `src/generated/notebooks/`. Source: `scripts/build-notebooks.py`. |
| `npm run extract`             | All three of the above, in order. |
| `npm run astro:dev`           | `astro dev` without extracting anything — handy when you're editing only the Astro layer and trust the existing inputs. |
| `npm run astro:build`         | `astro build` without extracting anything. |

## Workflow patterns

- **First clone / fresh checkout:** `npm install && npm run dev`.
- **Iterating on the homepage / layout:** keep `npm run astro:dev` running; re-run `npm run extract:*` only when you change the underlying source.
- **Changed a Rust public API:** `npm run extract:rust` then refresh the browser. The site picks up the regenerated `src/data/rust-api.json` automatically.
- **Edited a notebook under `notebooks/`:** `npm run extract:notebooks`.
- **Added a *new* notebook:** drop it into `notebooks/<name>.py` as a Jupytext-percent file, then `npm run extract:notebooks`. The Examples landing page picks it up from the regenerated `src/generated/notebooks/index.json`.

## Layout

```
docs/
├── notebooks/                       # Jupytext .py — executed at build time
├── scripts/
│   ├── extract-rust.mjs
│   ├── extract-python.mjs
│   └── build-notebooks.py
├── src/
│   ├── data/                        # *.json from extract:rust / extract:python  (gitignored)
│   ├── generated/                   # notebooks/*.html + meta + index           (gitignored)
│   ├── components/, layouts/, pages/, styles/
│   └── ...
├── package.json
└── astro.config.mjs
```

`src/data/` and `src/generated/` are `.gitignore`'d; they're rebuilt on
every CI run via `npm run extract`.

## CI

`.github/workflows/docs.yml` builds the site via `npm run build` on
every PR that touches `docs/`, `crates/`, `ppvm-python/`, or itself.
PR previews land at `gh-pages/pr-preview/pr-<N>/`; the main deploy
publishes to `gh-pages` root on every push to `main`.

## Prerequisites

- Node ≥ 20 (Astro 5 requirement).
- `uv` (https://docs.astral.sh/uv/) — used to drive the Python
  notebook executor and the griffe-based Python extractor.
- Rust nightly — `extract-rust.mjs` invokes `cargo +nightly rustdoc
  --output-format json`. Install with `rustup toolchain install
  nightly`.
