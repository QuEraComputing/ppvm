# noise_channels/

Hand-written noise corner cases. The per-channel breadth lives in
`generated/noise_sweeps/`; this category covers boundary probabilities
(p=0.0, p=1.0), ordering invariants, and noise-channel composition that
sweeps don't naturally produce.

Distribution-mode fixtures use `num_shots=4096` because the circuits are
tiny and statistical drift is the exact signal we lock down.

Loss tests are deliberately absent (per spec Non-Goals: Stim has no oracle
for `I_ERROR[loss]`). Loss is exercised by Rust unit tests in
`crates/ppvm-stim/tests/executor.rs`.
