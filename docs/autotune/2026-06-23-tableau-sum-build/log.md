# Log for 2026-06-23-tableau-sum-build

## 2026-06-23
- Architecture Notes / baseline profile (samply, 2733 samples):
- Target: examples/msd-noisy build time. Baseline build_median ~2620ms, per_shot ~22.5us, final branches=2025. Config Byte8F64<2> (storage [u64;2]=128bit), index u128, 85 qubits => tableau has 170 rows, each row ~32B word data.
- INCLUSIVE: for_each_mut_with_keys 85%; depolarize1 53%; fork(clone) 47%; loss_channel 38%; rebuild_fingerprints_if_dirty 23%; mimalloc alloc ~15-20%.
- SELF: _platform_memmove 32% (the tableau deep-clone in fork); rebuild_fingerprints_if_dirty 18% (re-hashes all words after every clifford gate marks dirty); for_each_mut_with_keys 11%; phase_loss_hash 5%; gates y/cz/sqrt_* ~10% total.
- Root causes: (1) noise branching deep-clones a full ~7KB tableau per branch; depolarize forks 3x/entry, ~85% of branches are then merged or truncated -> wasted clones. (2) every clifford gate marks all entries dirty -> next noise op re-hashes all words of all entries.
- Accuracy guard: branch count must stay 2025 (optimizations must not change the math). Cutoff fixed at 1e-7.
- Bench cmd: cargo build --release -p ppvm-tableau-sum --example msd-noisy-bench && ./target/release/examples/msd-noisy-bench
