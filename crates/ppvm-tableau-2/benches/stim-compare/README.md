# ppvm-tableau-2 vs Stim

Development benchmarks comparing the stabilizer core of `ppvm-tableau-2`
(`InverseTableau`, through the `run_stim` example) with official Stim's
`TableauSimulator`. Ported from PR 236's `bench/stim-compare`. Possibly not
shipped; nothing in the crate depends on this directory.

```bash
cd crates/ppvm-tableau-2/benches/stim-compare
uv run --no-project --with stim python sweep.py            # sizes, three code families
uv run --no-project --with stim python gates.py            # CX and H cost over n
uv run --no-project --with stim python stim_bench.py <circuit.stim> <shots>
```

Outputs (circuits, `sweep.csv`) go to `target/stim-compare/`. `--extra
NAME=PATH` adds a reference binary with `run_stim`'s CLI
(`<tool> <circuit> <shots>`, printing `<x> µs/shot`).

## What is measured

- Time per shot of gates, measurements and resets on a fresh simulator per
  shot. Noise is removed (`circuits.py`), so both sides run the same
  noise-free circuit; noise sampling is the circuit runner's job.
- **sweep.py**: Stim-generated memory circuits (surface d = 3…31, repetition
  d = 3…675, color d = 3…43, rounds = d), alternating tools over rounds,
  median per tool. The color code's `C_XYZ`, `MX`, `MY` are rewritten into
  gates ppvm's harness supports, for both tools.
- **gates.py**: a brickwork `CX` layer and a full `H` layer, timed at two
  repeat counts; the difference gives ns per gate without per-shot overhead.

## Caveats

- **PyPI Stim** is a generic build and is called from Python once per shot,
  which adds a few µs per shot. Per-shot comparisons on small circuits favor
  ppvm; `gates.py` cancels the overhead. A natively built Stim
  (`-march=native`) is the stricter baseline; pass it with `--extra`.
- **Instruction batching differs**: Stim merges consecutive same-gate
  instructions and shares measurement setup per instruction; ppvm's
  `stim-parser` does not merge (issue 291). See the design doc's Benchmarks
  section.
- This compares the stabilizer core. Once `GeneralizedTableau` exists, its
  `k = 0` overhead (active-mask checks) should be compared separately.
