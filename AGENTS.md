# AGENTS.md

> **Read the Developer Guide first.** The canonical contributor reference for
> this repository — for both human and AI contributors — lives in the docs
> site at `docs/src/pages/develop.astro` (rendered at `/develop/`). It covers
> project layout, build/test commands, architecture, conventions, the Python
> binding pipeline, and where to look for each subsystem.
>
> This file is a short pointer so agents can find that guide quickly. Do not
> add content here that belongs in the Developer Guide — keep the guide as
> the single source of truth.

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

1. Open `docs/src/pages/develop.astro` and read the sections relevant to your
   task. The "For AI agents" callout at the top tells you which sections are
   load-bearing.
2. Use `uv` for anything Python; never `pip`.
3. Use Conventional Commits: `<type>(<scope>): <description>`.
4. Build & test with the commands documented in the guide
   (`cargo test --workspace`, `uv run --project ppvm-python --group dev pytest …`).
5. Respect the `Config`-trait generics in `ppvm-runtime`; do not introduce
   runtime dispatch where a compile-time bound suffices.
6. Pauli propagation runs **backwards** (Heisenberg picture). Reverse the
   gate order accordingly when writing tests.

## Workspace at a glance

```
crates/ppvm-runtime         # Core: Pauli arithmetic, PauliSum, traits, Config
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
