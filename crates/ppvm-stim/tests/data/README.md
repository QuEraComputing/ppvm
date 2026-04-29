# Stim corpus

Fixtures consumed by `crates/ppvm-stim/tests/stim_corpus.rs`. Each fixture is
two committed files:

- `<name>.stim` — the source circuit.
- `<name>.expected.json` — declares the test's expected behavior in one of three
  modes: `deterministic`, `distribution`, or `unsupported`. See the spec
  (`docs/superpowers/specs/2026-04-28-ppvm-stim-test-corpus-design.md`) for
  the schema.

The harness walks the directory tree recursively and asserts that every
`.stim` has a sibling `.expected.json` (and vice versa). It then dispatches
on the JSON's `mode` field.

## Categories

| Subdir | Source | Purpose |
|---|---|---|
| `edge_cases/` | hand-written | Empty programs, REPEAT, every tag/Pi-expression form, dense/sparse measurement, comments/whitespace stress. |
| `noise_channels/` | hand-written | Boundary probabilities (p=0.0, p=1.0) and ordering corner cases. |
| `unsupported/` | hand-written + `regen-stim unsupported` | One fixture per phase-1-unsupported instruction. Flips to `distribution` in phase-2. |
| `generated/codes/` | `regen-stim codes` | `stim gen` sweeps over surface/repetition/color codes. |
| `generated/noise_sweeps/` | `regen-stim noise-sweeps` | Per-channel parameter sweeps. |
| `generated/dialect/` | `regen-stim dialect` | ppvm-specific `I[R_X(...)]`, `S[T]`, etc. |
| `generated/random/` | `regen-stim random` | Random-walk programs. |

## Provenance

- Hand-written fixtures (`edge_cases/`, `noise_channels/`, hand-written
  `unsupported/`) are authored by the ppvm team.
- Generated fixtures are produced by the `regen-stim` Python CLI in
  `crates/ppvm-stim/tests/regen-stim/`. See its README for invocation.

## Regenerating

`regen-stim` is a uv-managed Python tool. Cross-check against Stim happens at
regen time, not at test time:

```bash
cd crates/ppvm-stim/tests/regen-stim
uv sync
uv run regen-stim all          # regenerate everything
uv run regen-stim codes        # subcommand-specific
uv run regen-stim refresh ../data/<category>/<name>.stim
```

The committed `ppvm_bit_means` are bit-exact-compared against ppvm's output
at `cargo test` time. Bit drift here means ppvm's behavior changed — that's
the signal we want.
