# regen-stim

Generates ppvm-stim test corpus fixtures by cross-checking against
[`quantumlib/Stim`](https://github.com/quantumlib/Stim). Run by humans, not in CI.

## Setup

```bash
cd crates/ppvm-stim/tests/regen-stim
uv sync
```

`uv sync` resolves the editable `ppvm` Python package. The first sync requires
the `ppvm-python-native` extension built — run from the repo root once:

```bash
uv run --project ppvm-python --group dev maturin develop --uv
```

## Commands

```bash
uv run regen-stim codes        # generated/codes/  (stim gen sweeps)
uv run regen-stim noise-sweeps # generated/noise_sweeps/
uv run regen-stim dialect      # generated/dialect/
uv run regen-stim random       # generated/random/
uv run regen-stim unsupported  # unsupported/
uv run regen-stim refresh ../data/<category>/<name>.stim
uv run regen-stim verify  ../data/<category>/<name>.stim
uv run regen-stim all          # everything
```

## Dev

```bash
uv run pytest test/   # unit tests for the seed-search loop and tolerance math
```

## How it works

For distribution-mode fixtures, the seed-search loop in
`src/regen_stim/core.py` runs Stim at a high shot count to compute reference
per-bit means, then tries ppvm seeds in [0, 32) until ppvm's empirical means
(at the test-time `num_shots`) are within `tolerance_sigma * sigma` of Stim's.

Test-time runs use the committed seed at the committed shot count and
bit-exact-compare against the committed `ppvm_bit_means`. The cross-check
against Stim only happens here, at regen time.
