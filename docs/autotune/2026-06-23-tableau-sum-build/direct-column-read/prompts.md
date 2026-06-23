# Approach: direct-word column reads in hot loops

## Hypothesis
`pauli_error`'s per-entry dx/dy/dz loop and `structurally_equal_mutated`'s Pauli
branch read the tableau column with `bitvec`'s generic `Index` (`pw.word.xbits[addr0]`),
which recomputes word/bit and bounds-checks per access — done for all 2n rows of
every entry on every depolarize (part of the 23% `for_each_mut_with_keys` self).
Replace with direct storage-word access: compute `word_idx = addr0 / bits_per_word`
and `bit = addr0 % bits_per_word` ONCE, then test `(data.as_raw_slice()[word_idx] >> bit) & 1`.
Same bit values (Lsb0, matches `Tableau::build_masks`), so branches stay 2025.

## Target
`./target/release/examples/msd-noisy-bench`; baseline now build_median ~542ms.
Keep branches=2025, sum_p2=0.725135705447, top5[0]=0.8515413524292632, per_shot ~22us.
