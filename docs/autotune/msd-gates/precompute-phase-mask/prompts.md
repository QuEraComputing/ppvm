# Approach: precompute-phase-mask

## Hypothesis
`compute_phase` is called once per coefficient entry during T gates and measurements.
It contains an O(n) loop over all destabilizers checking `destab.phase % 2 != 0`.
With 32 coefficients and 85 measurements, this O(n) loop runs thousands of times.

By precomputing a bitmask of destabilizers with odd phase (once per operation),
the O(n) loop becomes a single `(active & odd_phase_mask).count_ones()` — O(1).

## Target metrics
- msd/msd-0: baseline 181 µs
- measurement/generalized: baseline 1.04 µs

## Files to modify
- `crates/ppvm-tableau/src/data.rs` — modify `compute_phase` and callers
