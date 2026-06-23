<!--
SPDX-FileCopyrightText: 2026 The PPVM Authors
SPDX-License-Identifier: Apache-2.0
-->

# `measure_many` profiling: where the time goes (2026-06-22)

Question: beyond reusing a common `MeasureScratch`, is there an algorithmic win
for `measure_many`? This note records the profiling that answers it. **No
optimization has been applied yet** — these are findings to act on later.

## Harness

- `examples/profile_measure_many.rs` — builds a full-width (85-qubit) entangled
  Clifford state so every measurement takes the case-a (random) path, with `n_t`
  T-gates controlling the coefficient count. Three workloads bracket the regimes:
  `few` (~1 coefficient), `mid` (~128), `large` (~1024). `quick` mode prints
  timings + the achieved coefficient count; `flame` mode runs a sustained loop.
- `scripts/profile_measure_many.sh` — records all three with **samply** (no
  sudo on macOS, unlike dtrace/flamegraph) into `target/profiles/`, then prints
  a top-functions summary.
- `scripts/samply_top.py` — turns a samply `profile.json.gz` (+ presymbolicated
  `.syms.json` sidecar) into a top-functions-by-self/inclusive-time table,
  headless (no browser UI needed).

Reproduce: `./scripts/profile_measure_many.sh`

## Results (self time; per-shot from `quick` mode)

Per-shot `measure_many` cost: few 17.7 µs, mid 21 µs, large 48 µs (85 qubits).

| function (self)                          | few (1) | mid (128) | large (1024) |
|------------------------------------------|---------|-----------|--------------|
| `compute_decomposition` (tableau scan)   | 39%     | 33%       | 14%          |
| `measure_one_with_scratch` *(see note)*  | 56%     | 54%       | 49%          |
| `insert` (build case-a HashMap)          | ~1%     | 7%        | 23%          |
| `retain` (drain Vec↔HashMap)             | ~1%     | 3%        | 9%           |
| memset / memmove                         | 1.4%    | 1.3%      | 4%           |

*Note: `measure_with_scratch`, `project_case_a`, and
`update_tableau_according_to_outcome` are inlined into `measure_one_with_scratch`.*

## Conclusions

The cost is two distinct regimes:

1. **Tableau-bound floor (dominates few/mid).** `compute_decomposition` costs a
   *constant* ~6.8 µs regardless of coefficient count (39% at 1 coeff → 14% at
   1024, same absolute µs). With the inlined `update_tableau_according_to_outcome`,
   the tableau work is a fixed **~16 µs O(n²)-per-batch floor** (each of the n
   measured qubits scans all 2n stabilizer/destabilizer rows). This is ~95% of
   the cost on Clifford-heavy / low-coefficient states.

2. **Coefficient-bound `Vec`↔`HashMap` round-trip (dominates large).** Each
   case-a measurement drains `self.coefficients` (`Vec`) into the scratch
   HashMap, works, then drains back; the next qubit rebuilds it. At 1024
   coefficients `insert` + `retain` + their memset/memmove churn is **~35%
   (~17 µs)** of the cost, and ~0% at low coefficient counts. The shared scratch
   keeps the *allocation*, but not the *representation*, across qubits.

## Recommended next steps (not yet implemented)

1. **Keep case-a coefficients HashMap-resident across the batch** — convert
   `Vec`→`HashMap` once and back once *per batch* instead of per qubit.
   Estimated ~30–35% on coefficient-heavy states; well-scoped refactor of
   `measure_with_scratch` / `project_case_a` that leaves the RNG draw order
   unchanged (so the existing characterization tests still apply).
2. **Batch the O(n²) tableau floor** — a blocked/incremental measurement that
   avoids rescanning all 2n rows per qubit. Larger, riskier algorithmic change;
   payoff concentrated on cheap (low-coefficient) states. Treat separately.
