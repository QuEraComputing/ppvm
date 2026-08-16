# Log for 2026-08-16-orbit-canonicalize

## 2026-08-16

Task: remove the per-term tax that makes the momentum (orbit-rep) path
slower than real space, blocking the k-resolved D(k) benchmark.

### Diagnosis

`TranslationGroup::canonicalize` and `canonicalize_with_shift` enumerated
the group by *rebuilding each element from the identity*: for group-element
index `idx` they decoded a mixed-radix counter and applied generator `g`
exactly `c[g]` times starting from `w`. For a single-generator group
(chain / ladder) that is `Σ_{c<|G|} c = |G|²/2` calls to `apply_generator`,
each `O(N)` with a fresh `PauliWord` + `rehash`, plus a `Vec<u32>` counter
allocation per element.

So the real cost was **O(|G|² · N)**, not the `O(|G| × n_qubits)` both
docstrings claimed. `orbit()` had the same shape. Since the momentum step
canonicalizes every action output of every rep, four times per `pc_step`
(leakage, predictor `build_orbit_rep_cols`, second-hop leakage, corrector),
this was essentially the entire momentum overhead — the earlier
"momentum is 13.5× slower per term than real space" note
(2026-07-07, ledger 2026-07-03-orbit-compare) is explained by it, and by
the group-order growth, not by phase bookkeeping.

### Fix (this commit)

1. **Odometer walk.** Applying generator `g` once always advances digit `g`
   *cyclically*, so a roll-over is one more application, not a rebuild.
   `advance(cur, idx)` applies generator 0 always, generator `g` only when
   `idx` is a multiple of `orders[0..g]`. Total `O(|G|)` applications per
   canonicalization, allocation-free. Used by `canonicalize`,
   `canonicalize_with_index` and `orbit`.
2. **Character table.** New `canonicalize_with_index` returns the
   mixed-radix *index* of the group element instead of a counter `Vec`;
   `character_table(k_modes)` precomputes all `|G|` characters once per
   call. The hot loops (`build_orbit_rep_cols`, `leakage_orbit_rep`,
   `canonicalize_pauli_sum_complex`, `check_momentum_sector`) index the
   table instead of calling `Complex::from_polar` (a `sin`/`cos` pair) and
   allocating a counter per action term.
   `canonicalize_with_shift` is kept as a thin wrapper (public API, and the
   error path of the sector check still wants the counter).

### Measurements (M4 MacBook Air, release, min-of-3, back-to-back
before/after builds of the SAME commit — not vs. the stale July wheel)

`TranslationGroup::canonicalize`, ladder(L,2), weight-4 words, µs/call:

| L (=\|G\|) | 4 | 8 | 16 | 32 | 64 | 96 |
|---|---|---|---|---|---|---|
| before | 0.52 | 0.99 | 4.05 | 27.03 | 198.90 | 657.40 |
| after  | 0.23 | 0.34 | 0.84 |  2.74 |   9.10 |  19.25 |
| speedup | 2.3× | 2.9× | 4.8× | 9.9× | 21.9× | 34.2× |

Before tracks `|G|²N`, after tracks `|G|N`; the speedup is `≈ |G|/2.8`, so
it keeps growing with system size.

End-to-end `pc_step`, B = 30 000 live entries, dt = 0.1, µs per basis entry:

| | momentum k=1 before | momentum k=1 after | real space |
|---|---|---|---|
| L=11 | 23.11 | 8.08 | 1.38 |
| L=16 | 61.37 | 12.99 | 1.49 |

Momentum/real per-entry ratio: 15.6× → 5.9× (L=11), 35.1× → 8.7× (L=16).

**Consequence for D(k).** A single k-mode holds ~|G|× fewer reps than the
equivalent real-space run, so a one-mode run now costs roughly
`(B/|G|)·8.7 / (B·1.49) ≈ 0.36` of a full real-space run at L=16 — a ~2.8×
*win* where before the fix it was a ~2.2× loss. Full-MSD-via-all-|G|-modes
is still the wrong benchmark (the mode count cancels the compression);
single-k D(k) is the right one.

Validation: `cargo test -p ppvm-pauli-sum -p ppvm-lindblad` (incl.
`pc_step_orbit_rep_matches_full_basis_projection`,
`complex_full_matches_real_at_kzero`,
`pc_step_matches_symmetry_merged_on_small_chain`) and the 34 Python
lindblad/momentum tests all pass. Two new unit tests pin the fix:
`canonicalize_with_index_matches_shift_and_table` (index form agrees with
the counter form and the table with `character()`, on a 2-generator
mixed-radix group) and `odometer_walk_covers_the_whole_group` (the
incremental walk still enumerates each element exactly once, and every
member canonicalizes to the lex-min).

### Not done (next lever)

For single-generator groups the canonical rep is a least-rotation problem,
solvable in `O(N)` by Booth/Duval instead of `O(|G|·N)` — that would close
the remaining ~6-9× per-entry gap to real space, at the cost of changing
the rep convention (any consistent choice is valid). Worth it only if D(k)
turns out to be canonicalization-bound after this fix.

Note: `cargo clippy --all-targets -p ppvm-lindblad` fails on the two
benches (`drug_dipolar`, `kossakowski`) with "missing generics for type
alias `Word`" — pre-existing on `wide-pauli-words` (verified by stashing
this change), unrelated.
