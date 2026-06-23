# Approach: lazy branch materialization (loss_channel + pauli_error)

## Hypothesis
Build time is dominated by `fork` (deep-clone of a ~7KB `GeneralizedTableau`,
47% inclusive, 32% `_platform_memmove` self). Each `depolarize1` forks 3 full
tableaux per entry and `loss_channel` forks 1 per entry, but ~85% of those
branches are immediately merged into an existing entry or dropped below
`sum_cutoff`. Those clones are pure waste.

Key fact (verified in `ppvm-tableau/src/gates/clifford.rs`): applying X/Y/Z to a
`GeneralizedTableau` only flips per-row **sign bits** and leaves the Pauli
words, the `coefficients` vector, and `is_lost` identical to the parent. Loss
only sets one `is_lost` bit. So a branch's fingerprint and structural identity
are derivable from the parent **without cloning**. Materialize (clone+mutate)
only when a branch survives as a *new* entry.

Per-row sign-flip rule at column `addr0` (matches the gate code exactly):
- X flips sign of row iff `z[addr0] == 1`
- Y flips sign of row iff `x[addr0] ^ z[addr0] == 1`
- Z flips sign of row iff `x[addr0] == 1`
(Only phase bit 1 = sign; the imaginary bit 0 is untouched. So the phase/loss
hash delta is `XOR sign_mask(row)` over flipped rows — same as the existing
`pauli_branch_phase_loss`.)

## Target metric
`cargo build --release -p ppvm-tableau-sum --example msd-noisy-bench && ./target/release/examples/msd-noisy-bench`
Baseline: build_median ~2620ms, per_shot ~22.5us, branches=2025,
sum_p2=0.725135705447, top5[0]=0.8515413524292632. The math must be unchanged:
branches stays 2025 and the accuracy fingerprint must match to ~1e-9.

## Expected win
Replace ~3N depolarize clones + ~N loss clones with clones only for survivors
(~1.2N total entries). Target 25-40% build-time reduction.
