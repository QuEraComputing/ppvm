# Approach: direct single-pass word_fingerprint (no gather, no thread_local)

## Hypothesis
`word_fingerprint` currently gathers all row words into a thread_local buffer and
calls `gxhash64` (the `with` frame is 22% self). The gather (memcpy ~2.7KB/entry),
the thread_local access, and the separate hash pass have overhead. Replace with a
direct single-pass scalar hash that reads the row words straight from the tableau
(no buffer, no thread_local, one pass). fxhash-style mixing is proven adequate
here (it's the existing wasm fallback) and `structurally_equal` resolves any
extra collisions (it's <1% of runtime, lots of headroom). Bonus: the hash becomes
portable (same on native + wasm), simplifying the cfg.

## Target
`./target/release/examples/msd-noisy-bench`; baseline now build_median ~552ms.
Keep branches=2025, sum_p2=0.725135705447, top5[0]=0.8515413524292632, per_shot ~22us.
