# Phase 5: Fused Cross-Site Kernel with Weight-Filtered Generation

## Motivation

### The scaling problem

With a dense rate matrix (N operators, N×N rates), `LindbladOp::new` creates N² individual
`LindbladTerm::Ladder` entries. Each entry independently iterates over the full observable
in `apply()`. The per-step cost is O(N² × |P|), but the real problem is worse: **cross-site
terms cause the observable to grow rapidly in Pauli weight**.

A weight-1 term like Z_i, processed by cross-site pair (i,j), produces weight-2 terms at
qubits (i,j). Next step, those weight-2 terms spawn weight-3 terms, and so on. The cascade:

| Step k | New weight | Term count scale |
|--------|-----------|-----------------|
| 0      | 1         | O(N)            |
| 1      | 2         | O(N²)           |
| 2      | 3         | O(N³)           |
| k      | k+1       | O(N^{k+1})     |

With diagonal rates (Vector), on-site terms can only maintain or reduce weight — no growth
at all. The dense rate matrix makes both the per-step cost AND the observable size scale
poorly with N. For N=20 or beyond, the solve becomes impractical.

### The fix: weight-filtered generation in a fused cross-site kernel

With a Pauli weight cutoff `w_max`, we can **skip generating terms that would exceed the
weight limit** before doing any expensive work (set_new_2, hash map insertion). For an
observable term at the weight ceiling (weight = w_max), only cross-site pairs where **both
qubits are already non-identity** can contribute without exceeding the limit. That's at
most w_max × (w_max - 1) pairs instead of N × (N - 1).

The improvement factor is **(N / w_max)²** — e.g. **100×** for N=100, w_max=10.

This requires two structural changes:
1. **Fused kernel**: replace N² individual `LindbladTerm::Ladder` entries with a single
   data structure that stores the rate matrix directly and iterates over qubit pairs
   inline, enabling per-term weight checks.
2. **Observable-outer loop**: iterate over observable terms in the outer loop, caching
   per-term metadata (weight, active qubit set) for use across all qubit pairs.

### Additional optimisations (stacking)

Two further constant-factor optimisations apply on top of the fused kernel:

1. **Precomputed action tables**: for each (left_dir, right_dir) direction pair and each
   (p_qi, p_qj) input Pauli pair, precompute the output (out_qi, out_qj, factor) tuples.
   Eliminates `pauli_mul` (15% of step time), `re_phase` (4%), and
   `pauli_anticommutes` from the inner loop. 16 input pairs × 4 direction pairs = 64
   table entries, each with 0–4 outputs.

2. **Rate matrix symmetry**: rate matrices are always symmetric (γ_ij = γ_ji). Terms
   (i,j) and (j,i) produce identical output words with identical factors — this is
   guaranteed by the physics of the Lindblad superoperator. Processing only the upper
   triangle (i < j) with doubled weight halves the number of hash map insertions.

Combined estimated improvement for the cross-site kernel:

| Optimisation               | Factor          | Applies to           |
|----------------------------|-----------------|----------------------|
| Weight filter              | (N / w_max)²    | pair count           |
| Precomputed action tables  | ~1.3×           | per-pair arithmetic  |
| Symmetry (upper triangle)  | 2×              | pair count           |
| **Combined**               | **~2.6 × (N/w_max)²** | |

---

## Design

### New data structure: `CrossSiteLadder`

```rust
/// Fused cross-site ladder kernel data.
///
/// Replaces N(N-1) individual `LindbladTerm::Ladder` entries (qi != qj) with a single
/// structure that stores the rate matrix and operator metadata directly, enabling
/// weight-filtered generation and precomputed action tables.
pub(crate) struct CrossSiteLadder {
    /// (qubit_index, direction) for each operator, in input order.
    ops: Vec<(usize, LadderDirection)>,

    /// N×N rate matrix. Accessed as rates[i][j].
    /// Stored in full (not upper triangle) for flexible access patterns.
    rates: Vec<Vec<f64>>,

    /// Precomputed action tables indexed by [left_dir][right_dir][p_qi][p_qj].
    /// Each entry: list of (out_qi: Pauli, out_qj: Pauli, factor: f64).
    /// Factor includes the multiplicity and re_phase contributions.
    /// At most 4 entries per table cell; many cells have 0–2.
    tables: ActionTables,
}

/// Compact action table storage. Indexed by direction pair, then input Pauli pair.
struct ActionTables {
    /// Flat storage: 2 left_dirs × 2 right_dirs × 4 p_qi × 4 p_qj = 64 cells.
    /// Each cell stores up to MAX_ACTIONS entries inline.
    data: [ActionCell; 64],
}

const MAX_ACTIONS: usize = 4;

struct ActionCell {
    entries: [(Pauli, Pauli, f64); MAX_ACTIONS],
    count: u8,
}
```

### Changes to `LindbladOp`

```rust
pub struct LindbladOp<T: Config> {
    /// Generic terms + on-site Ladder terms (qi == qj). Unchanged from Phase 3.
    pub(crate) terms: Vec<LindbladTerm<T>>,

    /// Fused cross-site ladder kernel. Present when any Ladder-Ladder pairs have qi != qj.
    pub(crate) cross_site: Option<CrossSiteLadder>,
}
```

In `LindbladOp::new`, Ladder-Ladder pairs with `qi != qj` are **no longer** pushed to
`terms`. Instead, their operator metadata and rate matrix entries are collected into
`CrossSiteLadder`. On-site Ladder terms (qi == qj) and all Generic terms remain in `terms`.

### Where does `w_max` come from?

The weight cutoff is needed at `apply` time for the generation filter. The generation
filter and the truncation strategy must agree on `w_max` — if the filter skips generating
terms above `w_max` but truncation doesn't enforce that limit, results are silently wrong.
Conversely, if the user has to set `w_max` in two places, they will inevitably diverge.

**Design**: `w_max` comes from the `Strategy`, not from `SolverConfig`. The user sets it
in exactly one place: the `MaxPauliWeight` strategy (or a `CombinedStrategy` that includes
it). The kernel reads `w_max` from the observable's strategy at `apply` time.

This requires adding a `max_weight()` method to the `Strategy` trait in ppvm-runtime:

```rust
pub trait Strategy: Default + Clone + Copy {
    fn capacity(&self, n_qubits: usize) -> usize;
    fn truncate<S, V, H, M, W>(&self, map: &mut M) where ...;

    /// Maximum Pauli weight retained by this strategy.
    /// Returns `usize::MAX` if no weight limit is enforced (default).
    /// Used by the cross-site kernel for generation-time filtering.
    fn max_weight(&self) -> usize { usize::MAX }
}
```

Implementations:
- `NoStrategy`: returns `usize::MAX` (default).
- `CoefficientThreshold`: returns `usize::MAX` (default).
- `MaxPauliWeight(w)`: returns `w`.
- `Budget`: returns `usize::MAX` (default).
- `CombinedStrategy<S1, S2>`: returns `min(s1.max_weight(), s2.max_weight())`.

The value is read inside `rhs_into` from the result PauliSum's strategy (which has the
same strategy type as the observable being evolved). No `SolverConfig` change needed.

**Note**: This requires a change to ppvm-runtime (adding `max_weight()` to the trait).
This is a non-breaking addition (default implementation returns `usize::MAX`), but it
crosses the crate boundary. Task 32 handles the ppvm-runtime change.

### Apply kernel structure (final form, after symmetry — Task 33)

```rust
fn apply_cross_site(
    cs: &CrossSiteLadder,
    p: &PauliSum<T>,
    result: &mut PauliSum<T>,
    w_max: usize,
) {
    let n_ops = cs.ops.len();

    for (w_a, coeff_a) in p.data().iter() {
        let weight = w_a.weight();
        let headroom = w_max.saturating_sub(weight);

        // Collect active (non-identity) operator indices for weight filtering.
        let active: SmallVec<[usize; 32]> = (0..n_ops)
            .filter(|&idx| w_a.get(cs.ops[idx].0) != Pauli::I)
            .collect();

        match headroom {
            0 => {
                // At ceiling: only pairs where BOTH qubits are already active.
                for ai in 0..active.len() {
                    for aj in (ai+1)..active.len() {
                        process_pair(cs, w_a, coeff_a, active[ai], active[aj], result);
                    }
                }
            }
            1 => {
                // One unit of headroom: skip pairs where BOTH qubits are identity.
                for i in 0..n_ops {
                    let i_active = w_a.get(cs.ops[i].0) != Pauli::I;
                    for j in (i+1)..n_ops {
                        if !i_active && w_a.get(cs.ops[j].0) == Pauli::I { continue; }
                        process_pair(cs, w_a, coeff_a, i, j, result);
                    }
                }
            }
            _ => {
                // Ample headroom: all pairs.
                for i in 0..n_ops {
                    for j in (i+1)..n_ops {
                        process_pair(cs, w_a, coeff_a, i, j, result);
                    }
                }
            }
        }
    }
}

/// Processes the merged (i,j)+(j,i) contribution for a symmetric rate matrix.
/// Uses doubled weight and a single table lookup per output (symmetry verified
/// for all direction pairs — see Resolved Question 1).
#[inline]
fn process_pair(cs, w_a, coeff_a, i, j, result) {
    let gamma = 2.0 * cs.rates[i][j];  // doubled: merges (i,j) and (j,i)
    if gamma == 0.0 { return; }
    let qi = cs.ops[i].0;
    let qj = cs.ops[j].0;
    let left_dir = cs.ops[i].1.flip();
    let right_dir = cs.ops[j].1;
    let cell = cs.tables.get(left_dir, right_dir, w_a.get(qi), w_a.get(qj));
    for k in 0..cell.count {
        let (out_qi, out_qj, factor) = cell.entries[k];
        *result += (w_a.set_new_2(qi, out_qi, qj, out_qj),
                    (gamma * factor).into() * *coeff_a);
    }
}
```

Task 32 implements this without symmetry first (full `i != j` iteration, single table
lookup per pair). Task 33 switches to `i < j` with doubled weight.

### Parallel dispatch

The current `PAR_THRESHOLD = 200` check on `self.terms.len()` will undercount work after
moving cross-site terms out of `terms`. The fused kernel has its own parallelism
opportunity: **parallelize over observable terms** (the outer loop). This is a natural
split because each observable term's contributions are independent until accumulated into
`result`.

The parallel strategy for the fused kernel should be a separate fold/reduce pattern
(similar to `apply_par`), iterating over `p.data().par_iter()` and folding into
thread-local `PauliSum`s.

### Interaction with existing code

- **On-site Ladder terms**: unchanged, remain in `terms`, processed by existing `apply()`.
- **Generic terms**: unchanged, remain in `terms`.
- **Hamiltonian commutator**: unchanged, called before `apply()` in `rhs_into`.
- **Truncation**: called after `apply()` in `rhs_into`, as before. The weight filter in
  the generation loop is a *pre-filter* that avoids generating doomed terms; truncation
  still runs and handles coefficient-based pruning.
- **DOPRI5 stepper**: unchanged. `rhs_into` reads `w_max` from the result's strategy
  internally; its external signature does not change.
- **SolverCache**: unchanged.

---

## Task breakdown

### Task 30 — Precomputed action tables

**Goal**: Build and test the action table data structure. Pure precomputation, no changes
to the hot path.

**Background**: For a cross-site Ladder term with direction pair (left_dir, right_dir)
acting on observable word w_a, the current kernel runs a 2×2 inner loop over
(left_sub, right_sub) ∈ {X, Y} × {X, Y}, performing `pauli_anticommutes`,
`pauli_mul`, phase arithmetic, and `re_phase` per combination. The action table
precomputes this for all 16 input (p_qi, p_qj) pairs, storing only the non-zero outputs.

**Steps**:
1. Define `ActionCell` and `ActionTables` types in `src/lindblad.rs` (see Design section).
2. Implement `ActionTables::build() -> Self`:
   - For each of the 4 direction pairs (left_dir, right_dir):
     - Compute `l_y_phase` and `r_y_phase` from the directions.
     - For each of the 16 (p_qi, p_qj) input pairs:
       - Run the existing 2×2 loop logic (pauli_anticommutes, pauli_mul, re_phase).
       - Collect non-zero (out_qi, out_qj, factor) results into the cell.
   - Factor includes `multiplicity * 2.0 * re_phase(total_phase)`.
3. Implement `ActionTables::get(left_dir, right_dir, p_qi, p_qj) -> &ActionCell`.

**Unit tests**:
- `action_table_zz_lower_raise`: verify table entry for (Lower, Raise, Z, Z) matches
  hand-computed values: [(Y, Y, 4.0), (X, X, 4.0)].
- `action_table_zi_lower_raise`: verify (Lower, Raise, Z, I) → [(Y, Y, -2.0), (X, X, -2.0)].
- `action_table_ii_all_dirs`: verify all 4 direction pairs give empty cells for (I, I).
- `action_table_matches_inline_kernel`: for all 64 (dir_pair, p_qi, p_qj) combinations,
  run the existing inline 2×2 kernel and verify the table produces identical outputs.

**Review checklist**:
- [ ] `ActionTables::build()` uses the exact same algebra as the existing cross-site kernel.
- [ ] No changes to `apply()`, `LindbladOp::new`, or any hot path.
- [ ] `action_table_matches_inline_kernel` covers all 64 combinations exhaustively.
- [ ] Benchmark: `bench_rhs` unchanged (table construction is not in the hot path).

---

### Task 31 — `CrossSiteLadder` structure and `LindbladOp::new` routing

**Goal**: Introduce the `CrossSiteLadder` struct. Modify `LindbladOp::new` to separate
cross-site Ladder-Ladder pairs from `terms` and store them in a new `cross_site` field.
No changes to `apply()` yet — the cross-site field is populated but not consumed.

**Steps**:
1. Define `CrossSiteLadder` struct in `src/lindblad.rs` (see Design section).
2. Add `pub(crate) cross_site: Option<CrossSiteLadder>` field to `LindbladOp<T>`.
3. In `LindbladOp::new`:
   - Collect `(qubit, direction)` metadata from all `JumpOp::Ladder` ops.
   - For Ladder-Ladder pairs with `qi != qj`: continue pushing `LindbladTerm::Ladder`
     to `terms` (keep old path active), AND additionally accumulate the rate matrix
     entry for `CrossSiteLadder`.
   - After the loop: if any cross-site Ladder pairs were found, construct
     `CrossSiteLadder { ops, rates, tables: ActionTables::build() }` and store in
     `self.cross_site`.
   - On-site Ladder pairs (qi == qj) and all Generic terms go to `terms` as before.
   - The `cross_site` field is populated but not yet consumed. Task 32 activates the
     fused kernel and stops pushing cross-site terms to `terms`.

**Unit tests**:
- `cross_site_populated_for_dense_rate_matrix`: n=3 all-Ladder with 3×3 dense rates.
  Assert `lop.cross_site.is_some()`, `cs.ops.len() == 3`, `cs.rates` is 3×3.
- `cross_site_none_for_diagonal_rate_matrix`: n=3 all-Ladder with Vector rates.
  Assert `lop.cross_site.is_none()` (no cross-site pairs for diagonal rates).
- `cross_site_none_for_generic_ops`: n=2 with Generic CollapseOps + dense rates.
  Assert `lop.cross_site.is_none()` (Generic ops don't participate).
- `cross_site_mixed_ops`: n=3 with 2 Ladder + 1 Generic and dense rates.
  Assert `lop.cross_site.is_some()` with `cs.ops.len() == 2` (only the 2 Ladder ops).
  Assert the Ladder-Generic and Generic-Ladder pairs are in `terms` as `Generic`.
- All existing tests pass (old path still active).

**Review checklist**:
- [ ] `CrossSiteLadder` stores ops, rates, and precomputed tables.
- [ ] `LindbladOp::new` correctly separates cross-site Ladder pairs.
- [ ] Mixed JumpOp scenarios handled: only Ladder-Ladder cross-site pairs go to `cross_site`.
- [ ] All existing tests still pass (old code path not yet removed).
- [ ] Benchmark: `bench_rhs` may show slight regression from double work; record it.

---

### Task 32 — Fused cross-site apply kernel with weight filtering

**Goal**: Implement the fused cross-site kernel in `apply()` using the action tables and
weight-based generation filtering. Remove the old per-term cross-site path.

**Steps**:
1. Add `max_weight(&self) -> usize` to the `Strategy` trait in ppvm-runtime (default
   impl returns `usize::MAX`). Override in `MaxPauliWeight` (returns `self.0`) and
   `CombinedStrategy` (returns `min(s1.max_weight(), s2.max_weight())`). This is a
   non-breaking change to ppvm-runtime.

2. In `rhs_into`, read `w_max` from `result.strategy().max_weight()` and pass it to
   `LindbladOp::apply()` (and `apply_par()`). The external signatures of `rhs_into`,
   `rhs_into_par`, `step`, and `solve` do not change. Only the internal `apply()` and
   `apply_cross_site()` functions receive `w_max` as a parameter.

3. Implement `apply_cross_site()` (see Design section):
   - Outer loop over observable terms.
   - For each term: compute `weight` via `w_a.weight()`, build active qubit set.
   - Branch on `headroom = w_max - weight`:
     - `headroom == 0`: iterate only over active × active pairs.
     - `headroom == 1`: iterate all pairs, skip if both qubits inactive.
     - `headroom >= 2`: iterate all pairs.
   - For each pair: look up action table, emit contributions via `set_new_2` + `+=`.
   - Iterate all pairs `i != j` (not yet upper-triangle; Task 33 adds symmetry).

4. Call `apply_cross_site()` from `apply()` after processing `self.terms`:
   ```rust
   if let Some(cs) = &self.cross_site {
       apply_cross_site(cs, p, result, w_max);
   }
   ```

5. Remove cross-site `LindbladTerm::Ladder` entries from `terms` in `LindbladOp::new`
   (stop the double-push from Task 31).

6. Update `rhs_into` parallel dispatch: the `PAR_THRESHOLD` check on `terms.len()` no
   longer accounts for cross-site work. For now, keep the existing threshold for `terms`
   and run the fused kernel sequentially. Parallel dispatch for the fused kernel is
   Task 34.

**Unit tests**:
- `fused_kernel_matches_old_path_n3`: n=3 all-Ladder with dense rates. Compare output
  of the fused kernel against the old per-term kernel (computed via Generic expand path)
  for a representative observable (weight-1 and weight-2 terms). Assert
  coefficient-level equality (within f64 epsilon).
- `fused_kernel_matches_old_path_n6`: n=6 all-Ladder (superradiance fixture). Full
  solve comparison: run the same solve with old code and new fused kernel, assert results
  match within solver tolerance.
- `weight_filter_skips_correctly`: n=4, w_max=2. Start with a weight-1 observable.
  After one RHS evaluation, verify no weight-3 or higher terms appear in the result.
  Compare against unfiltered result (w_max=MAX) truncated to weight ≤ 2 — should match.
- `weight_filter_headroom_1`: n=4, w_max=3. Observable with weight-2 terms. Verify that
  cross-site pairs where both qubits are identity are skipped (no weight-4 terms), but
  pairs where one qubit is active do contribute weight-3 terms.
- `w_max_usize_max_matches_unfiltered`: verify that w_max=usize::MAX produces the same
  result as the old unfiltered kernel.

**Review checklist**:
- [ ] `max_weight()` added to `Strategy` trait in ppvm-runtime with default `usize::MAX`.
- [ ] `w_max` read from strategy in `rhs_into`, passed internally to `apply`. No public API changes.
- [ ] `apply_cross_site` uses `ActionTables::get` — no `pauli_mul` or `re_phase` calls.
- [ ] Weight filter correctly handles headroom 0, 1, and ≥2.
- [ ] `w_a.weight()` called once per observable term (not per pair).
- [ ] Old cross-site `LindbladTerm::Ladder` entries no longer generated in `LindbladOp::new`.
- [ ] `fused_kernel_matches_old_path_n6` passes (full solve-level agreement).
- [ ] **Performance**: Run `bench_rhs` and new fused-kernel bench for n=6. Expect improvement
     from action tables even without weight filtering (pauli_mul elimination). Record numbers.

---

### Task 33 — Symmetry exploitation (upper triangle)

**Goal**: For symmetric rate matrices, process only the upper triangle (i < j) and
combine contributions from (i,j) and (j,i) into a single hash map insertion per output.

**Background**: For a symmetric rate matrix (γ_ij = γ_ji), terms (i,j) and (j,i) produce
identical output words with identical coefficient factors for any input observable term.
This has been **verified algebraically for all 4 direction-pair combinations and all 16
input Pauli pairs** (64 total cases, zero failures). The symmetry holds because:
- Same-direction pairs (Lower/Lower, Raise/Raise): both terms get identical direction
  pairs, so the table is the same and input transposition maps to output transposition.
- Mixed-direction pairs (Lower/Raise, Raise/Lower): the direction pairs swap, giving
  different l_y_phase/r_y_phase assignments (e.g. (3,3) vs (1,1)), but the total phases
  after combining with pauli_mul phases always produce the same re_phase values.

Since the merge is universally valid, no runtime fallback is needed. The kernel can
unconditionally use upper-triangle iteration with doubled weight.

**Steps**:
1. In `apply_cross_site`, change the pair iteration from `for i in 0..n, for j in 0..n,
   i != j` to `for i in 0..n, for j in (i+1)..n`:
   - For each (i, j) pair, use `gamma = 2.0 * rates[i][j]` (doubled weight).
   - Single table lookup and single insertion per output.
   - No conditional merge logic or fallback needed.

2. Add a `symmetry_verified` unit test that exhaustively checks all 64 (direction pair ×
   input Pauli pair) combinations, confirming that the (i,j) and (j,i) table entries
   produce identical outputs after qubit swap. This serves as a regression guard.

**Unit tests**:
- `symmetry_verified_all_64_cases`: exhaustively check all 4 direction pairs × 16 input
  Pauli pairs, confirming (i,j) and (j,i) table entries produce identical outputs after
  qubit swap. Regression guard for the algebraic verification.
- `symmetry_all_same_direction`: all-Lower or all-Raise operators. Verify upper-triangle
  output matches the full-iteration (Task 32) path.
- `symmetry_mixed_directions`: mix of Raise and Lower operators. Verify upper-triangle
  output matches the full-iteration path.

**Review checklist**:
- [ ] Upper-triangle iteration: `j in (i+1)..n_ops`, `gamma = 2.0 * rates[i][j]`.
- [ ] No conditional merge logic — unconditional doubled weight.
- [ ] Output matches full-iteration path for both same-direction and mixed-direction cases.
- [ ] **Performance**: `bench_rhs` for n=6 superradiance. Expect measurable improvement
     from halved hash map insertions. Record before/after numbers.

---

### Task 34 — Parallel dispatch for the fused kernel

**Goal**: Add rayon parallelism to `apply_cross_site` for large observables, following the
existing cold/never-inline isolation pattern.

**Steps**:
1. Implement `apply_cross_site_par()`:
   - `#[cold] #[inline(never)]` to isolate rayon atomics from the sequential path.
   - `p.data().par_iter().fold(|| local_result, |local, (w_a, coeff_a)| { ... }).reduce(...)`.
   - Each thread processes a chunk of observable terms independently.
   - Final reduce merges thread-local PauliSums into `result`.
2. Add a threshold check in `apply_cross_site` (or in `rhs_into`): dispatch to the
   parallel path when `|P| × n_pairs` exceeds a tunable threshold.
   - The threshold should account for both observable size and pair count (the fused
     kernel's work is proportional to their product).
3. Update `rhs_into` to use the new parallel dispatch instead of (or in addition to) the
   existing `PAR_THRESHOLD` on `terms.len()`.

**Unit tests**:
- `parallel_fused_matches_sequential`: same observable + LindbladOp, verify parallel and
  sequential paths produce identical results.
- All existing parallel tests still pass.

**Review checklist**:
- [ ] `apply_cross_site_par` is `#[cold] #[inline(never)]`.
- [ ] Threshold is based on estimated work, not just term count.
- [ ] Sequential path has zero rayon code in scope.
- [ ] **Performance**: benchmark parallel vs sequential for n=6 (small, expect sequential
     wins) and n=10+ (larger, expect parallel wins). Record crossover point.

---

### Task 35 — Benchmarks, integration tests, and superradiance example update

**Goal**: Comprehensive validation and performance measurement.

**Steps**:
1. Add new Criterion benchmarks:
   - `bench_rhs_fused_n6`: n=6 dense rates, all-Ladder, w_max=MAX (no filter). Compare
     against pre-Phase-5 baseline.
   - `bench_rhs_fused_n6_wmax3`: same but w_max=3. Measure filter benefit.
   - `bench_rhs_fused_n10_wmax4`: n=10, w_max=4. Larger system, filter benefit scales.
   - `bench_solve_superradiance_n6`: full solve, compare old vs new.
2. Integration test: `solve_superradiance_fused_matches_baseline`:
   - Run the superradiance n=6 problem with the old (pre-Phase-5) code path and the new
     fused kernel. Assert results agree within solver tolerance.
3. Update `examples/superradiance_flame.rs`:
   - Use the fused kernel (it should activate automatically for Ladder+Dense).
   - Use a `CombinedStrategy<MaxPauliWeight, CoefficientThreshold>` to demonstrate
     weight-filtered generation. The `MaxPauliWeight` value is read automatically by
     the kernel via `Strategy::max_weight()`.
4. Record all benchmark numbers in the commit message. Compare against Phase 3/4 baselines.

**Review checklist**:
- [ ] Criterion benchmarks cover: no filter, with filter, scaling with N.
- [ ] Integration test confirms numerical equivalence with old code path.
- [ ] superradiance_flame example updated and runs.
- [ ] Benchmark numbers recorded with explicit before/after comparison.

---

### Task 36 — Optimise `PauliWordTrait::weight()` using `count_ones` (ppvm-runtime)

**Goal**: Replace the per-bit loop in `PauliWord::weight()` with a byte-level `count_ones`
fast path. This is a ppvm-runtime change, exempt from the Phase 5 guideline of only
touching ppvm-timeevolve.

**Background**: `weight()` currently iterates over all N qubits checking
`xbits[i] || zbits[i]`. For N=100 (13 bytes of storage), that's 100 branch-heavy
iterations. A qubit is non-identity iff its xbit or zbit is set, so
`weight = popcount(xbits | zbits)` over the raw storage bytes. For `[u8; K]` storage,
this is K OR operations + K `count_ones()` calls — typically 1–2 instructions per byte
on ARM/x86 with hardware popcount.

**Steps**:
1. In `ppvm-runtime/src/word/data.rs`, change the `weight()` implementation on
   `PauliWord<A, S>`:
   ```rust
   fn weight(&self) -> usize {
       // Each qubit occupies 2 bits (xbit, zbit). Non-identity iff either is set.
       // OR the raw storage bytes and popcount.
       self.xbits.as_raw_slice().iter()
           .zip(self.zbits.as_raw_slice().iter())
           .map(|(x, z)| (x | z).count_ones() as usize)
           .sum()
   }
   ```
   Adjust based on the actual storage layout (BitArray vs raw bytes). The key requirement
   is to use `count_ones()` on integer types, not per-bit iteration.

2. Verify that the new implementation matches the old one for all edge cases (all-I, all-X,
   mixed, single qubit).

**Unit tests**:
- `weight_matches_reference`: for a representative set of Pauli words (n=1 through n=20),
  compute weight via the old per-bit loop and via the new implementation. Assert equality.
- `weight_all_identity`: n=10 all-I word → weight 0.
- `weight_all_nonidentity`: n=10 all-X word → weight 10.

**Review checklist**:
- [ ] Uses `count_ones()` or equivalent, not per-bit iteration.
- [ ] All existing ppvm-runtime tests pass.
- [ ] Benchmark: measure `weight()` cost for n=6, n=20, n=100 before/after.

---

## Resolved questions

### 1. Symmetry merge — always valid? ✓

**Verified.** Exhaustive algebraic check of all 4 direction-pair combinations × 16 input
Pauli pairs (64 cases total): symmetry holds universally. No fallback needed. See Task 33
background for details.

### 2. `w_max` correctness contract ✓

**Resolved: single source of truth.** `w_max` is read from the `Strategy` trait's
`max_weight()` method, not from `SolverConfig`. The user sets `MaxPauliWeight(w)` (or a
`CombinedStrategy` that includes it) on their `PauliSum`, and the generation filter reads
it from there. No separate configuration, no divergence possible. See "Where does `w_max`
come from?" in the Design section.

### 3. `weight()` cost ✓

**Resolved: accept for now + dedicated optimisation task.** The per-bit iteration in
`weight()` is O(N) per call. For small storage (1–3 bytes) this is fast. For larger N,
Task 36 adds a `count_ones`-based fast path in ppvm-runtime. This is the only Phase 5
change that touches ppvm-runtime source (beyond the `max_weight()` trait addition).

### 4. Active qubit set representation ✓

**Resolved: decide based on profiling.** Start with `SmallVec<[usize; 32]>`. For
w_max ≤ 20 (the common regime), the nested loop has at most ~400 iterations and the
SmallVec fits in a cache line. If profiling reveals bottlenecks for larger w_max, switch
to a bitset. No upfront over-engineering.

---

## Open questions

### 5. Parallel dispatch threshold

The current `PAR_THRESHOLD = 200` was tuned for per-term parallelism. The fused kernel's
work per observable term is proportional to the number of pairs processed (which depends
on weight and w_max). A simple threshold on `|P|` might not capture this well. Consider:
`parallel_if |P| * estimated_pairs_per_term > THRESHOLD`. The
`estimated_pairs_per_term` could be a conservative upper bound (e.g. `n_ops²`) or a
sampled average.

### 6. Non-ladder operators in cross-site position

If only some operators are Ladders, the `CrossSiteLadder` includes only the Ladder subset.
The remaining Ladder-Generic and Generic-Generic cross-site pairs stay as
`LindbladTerm::Generic` in `terms`. This means the weight filter does NOT apply to Generic
cross-site terms. Is this acceptable? For all-Ladder systems (the primary use case for
dense rate matrices), this is not an issue. For mixed systems, the Generic terms dominate
cost anyway and the filter provides no benefit.

### 7. Rate matrix storage

The current design stores the full N×N rate matrix in `CrossSiteLadder::rates`. For large
N (>100), this is 80KB+ of f64s. An alternative is to store only the upper triangle
(N(N-1)/2 entries) and compute the flat index. This saves memory but adds index arithmetic
to the inner loop. For N < 200 the full matrix fits comfortably in L1/L2 cache; optimise
only if profiling shows cache pressure.
