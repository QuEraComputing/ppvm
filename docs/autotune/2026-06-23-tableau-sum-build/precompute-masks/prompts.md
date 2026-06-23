# Approach: precompute per-row masks (kill redundant splitmix)

## Hypothesis
`sign_mask(row)`/`imag_mask(row)`/`loss_mask(q)` are splitmix64 hashes of a pure
index. They are recomputed per-row, per-entry, on every op:
- `pauli_error`'s dx/dy/dz loop computes `sign_mask(row)` for all 2n rows of
  every entry on every depolarize.
- `phase_loss_hash` (called per entry in `rebuild_fingerprints_if_dirty`)
  recomputes sign/imag masks per set phase and loss masks per lost qubit.

Precompute the per-index mask tables ONCE per op and index them. Values are
identical, so all fingerprints (and the accuracy fingerprint) are unchanged.

## Target
`./target/release/examples/msd-noisy-bench`; baseline now build_median ~573ms.
Keep branches=2025, sum_p2=0.725135705447, top5[0]=0.8515413524292632, per_shot ~22us.
