# `ppvm-usage` examples

Every large code block in [`../SKILL.md`](../SKILL.md) has a runnable twin
here. The skill is the source of truth for what an agent should *write*;
these files are the source of truth for what actually *compiles and runs*
against the current ppvm API.

```
examples/
├── python/        # ran by docs/examples/test_examples.py
│   ├── verify.py
│   └── noise_truncation.py
└── rust/          # workspace member; built + tested by cargo
    ├── Cargo.toml
    ├── src/bin/
    │   ├── paulisum.rs
    │   ├── tableau.rs
    │   ├── stim_sample.rs
    │   └── sym.rs
    └── tests/
        └── skill_bins_run.rs
```

The whole directory is colocated with `SKILL.md` so that `ion add
QuEraComputing/ppvm/skills/ppvm-usage` brings the examples along with
the instructions — agents can use them as a copy-paste starting point.

## How CI catches drift

- **Python:** `uv run --project ppvm-python --group dev pytest
  docs/examples/` runs every script and matches its stdout against the
  expected output recorded in `docs/examples/test_examples.py`.
- **Rust:** the `rust/` directory is a workspace member named
  `ppvm-skill-examples`. `cargo build --workspace --all-targets` compiles
  every binary (catches signature changes), and `cargo test -p
  ppvm-skill-examples` actually invokes each binary (catches runtime
  regressions).

If you change `SKILL.md`, update the matching file here in the same PR.
If you change the public ppvm API in a way that breaks one of these
files, the build will fail before the skill ships to users.

## A note on `examples/` vs `scripts/`

`ion validate` emits two advisory warnings on this directory
("Script file outside `scripts/`") because by ion's convention any
executable content shipped with a skill lives in `scripts/`. We
deliberately use `examples/` instead because these files are *reference
examples and CI fixtures*, not helper scripts the skill instructs an
agent to invoke — the skill itself never says "run
`examples/python/verify.py`". The warnings don't block installation.
