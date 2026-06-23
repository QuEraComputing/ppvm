# Approach: bulk word_fingerprint hashing

## Hypothesis
After lazy materialization, `rebuild_fingerprints_if_dirty` dominates (47% self,
61% inclusive); it re-hashes every entry's words after each clifford gate marks
them dirty. `word_fingerprint` currently does 2 small `Hash::hash` calls per row
(`xbits.data` then `zbits.data`) = ~340 hasher writes for 170 rows, with
per-call overhead. Gather the row words into one contiguous buffer and hash once
with `gxhash::gxhash64` (native) — far less per-call overhead, single SIMD pass.

## Target
`./target/release/examples/msd-noisy-bench`; baseline now build_median ~958ms.
Must keep branches=2025, sum_p2=0.725135705447, top5[0]=0.8515413524292632.
Fingerprints are transient dedup keys (resolved by structurally_equal), so the
hash VALUE may change freely as long as it's consistent within a build.
