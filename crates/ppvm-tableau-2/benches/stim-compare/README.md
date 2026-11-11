# ppvm-tableau-2 vs Stim

Development benchmarks comparing the stabilizer core of `ppvm-tableau-2`
(`InverseTableau`, through the `run_stim` example) with official Stim's
`TableauSimulator`. Ported from PR 236's `bench/stim-compare`. Possibly not
shipped; nothing in the crate depends on this directory.

```bash
# A reference: stim_driver.cc built against any Stim checkout (C++, -march=native).
cmake -S <stim-checkout> -B target/stim-build -DCMAKE_BUILD_TYPE=Release
cmake --build target/stim-build --target libstim -j
clang++ -O3 -std=c++20 -march=native -fno-strict-aliasing -I <stim-checkout>/src \
    crates/ppvm-tableau-2/benches/stim-compare/stim_driver.cc \
    target/stim-build/out/libstim.a -o target/stim-compare/stim_driver

cd crates/ppvm-tableau-2/benches/stim-compare
uv run --no-project --with stim python sweep.py --reference stim=$PWD/../../../../target/stim-compare/stim_driver \
    --case surface:7 --case color:25            # or no --case for the full sweep
uv run --no-project --with stim python gates.py --reference stim=... --size 64 --size 2048
```

Python Stim (`uv run --with stim`) only generates circuits. Outputs (circuits,
`sweep.csv`) go to `target/stim-compare/`. Each `--reference NAME=PATH` is a
binary with `run_stim`'s CLI (`<tool> <circuit> <shots>`, printing
`<x> µs/shot`); ratios are against the first.

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

- Compare against a native C++ build, as above. PyPI Stim is a generic build
  and its per-shot Python call adds a few µs, which flatters ppvm on small
  circuits.
- **Instruction batching differs**: Stim merges consecutive same-gate
  instructions and shares measurement setup per instruction; ppvm's
  `stim-parser` does not merge (issue 291). See the design doc's Benchmarks
  section.
- This compares the stabilizer core. Once `GeneralizedTableau` exists, its
  `k = 0` overhead (active-mask checks) should be compared separately.
