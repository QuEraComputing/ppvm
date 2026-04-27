# Phase 5 Tasks

## Performance policy

Run `cargo bench bench_rhs` before and after every task and record both numbers in the
commit message. If a task introduces a regression, document it explicitly. No regression
may be left unresolved at the end of the phase.

Benchmark fixture: n=6 dense all-Raise superradiance (36 Lindblad terms).
Use `--release` for all benchmarks.

---

## Task 30 — Precomputed action tables

**Goal**: Build and test the `ActionCell` / `ActionTables` data structures. Pure
precomputation, no changes to the hot path.

**Steps**:

1. In `src/lindblad.rs`, below the existing `pauli_mul` function (~line 272), add:

   ```rust
   const MAX_ACTIONS: usize = 4;

   #[derive(Debug, Clone, Copy)]
   pub(crate) struct ActionCell {
       pub entries: [(Pauli, Pauli, f64); MAX_ACTIONS],
       pub count: u8,
   }

   pub(crate) struct ActionTables {
       /// Flat storage: [left_dir][right_dir][p_qi][p_qj].
       /// Index = left_dir_idx * 32 + right_dir_idx * 16 + p_qi_idx * 4 + p_qj_idx
       /// where Raise=0, Lower=1; I=0, X=1, Y=2, Z=3.
       data: [ActionCell; 64],
   }
   ```

2. Implement index helpers:

   ```rust
   impl LadderDirection {
       fn idx(&self) -> usize {
           match self { LadderDirection::Raise => 0, LadderDirection::Lower => 1 }
       }
   }

   fn pauli_idx(p: Pauli) -> usize {
       match p { Pauli::I => 0, Pauli::X => 1, Pauli::Y => 2, Pauli::Z => 3 }
   }

   impl ActionTables {
       fn index(left_dir: LadderDirection, right_dir: LadderDirection,
                p_qi: Pauli, p_qj: Pauli) -> usize {
           left_dir.idx() * 32 + right_dir.idx() * 16
               + pauli_idx(p_qi) * 4 + pauli_idx(p_qj)
       }

       pub fn get(&self, left_dir: LadderDirection, right_dir: LadderDirection,
                   p_qi: Pauli, p_qj: Pauli) -> &ActionCell {
           &self.data[Self::index(left_dir, right_dir, p_qi, p_qj)]
       }
   }
   ```

3. Implement `ActionTables::build() -> Self`:

   ```rust
   impl ActionTables {
       pub fn build() -> Self {
           let empty_cell = ActionCell {
               entries: [(Pauli::I, Pauli::I, 0.0); MAX_ACTIONS],
               count: 0,
           };
           let mut tables = ActionTables { data: [empty_cell; 64] };

           let dirs = [LadderDirection::Raise, LadderDirection::Lower];
           let paulis = [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z];

           for &left_dir in &dirs {
               let l_y_phase: u8 = match left_dir {
                   LadderDirection::Raise => 3,
                   LadderDirection::Lower => 1,
               };
               for &right_dir in &dirs {
                   let r_y_phase: u8 = match right_dir {
                       LadderDirection::Lower => 1,
                       LadderDirection::Raise => 3,
                   };
                   let left_subs = [(Pauli::X, 0u8), (Pauli::Y, l_y_phase)];
                   let right_subs = [(Pauli::X, 0u8), (Pauli::Y, r_y_phase)];

                   for &p_qi in &paulis {
                       for &p_qj in &paulis {
                           let idx = Self::index(left_dir, right_dir, p_qi, p_qj);
                           let cell = &mut tables.data[idx];
                           cell.count = 0;

                           for &(l_pauli, l_phase) in &left_subs {
                               for &(r_pauli, r_phase) in &right_subs {
                                   let p1 = pauli_anticommutes(l_pauli, p_qi) as u8;
                                   let p2 = pauli_anticommutes(p_qj, r_pauli) as u8;
                                   let mult = p1 + p2;
                                   if mult == 0 { continue; }

                                   let (out_qi, qi_phase) = pauli_mul(l_pauli, p_qi);
                                   let (out_qj, qj_phase) = pauli_mul(p_qj, r_pauli);
                                   let total_phase = (l_phase as u16 + qi_phase as u16
                                       + qj_phase as u16 + r_phase as u16) as u8 % 4;
                                   let s = re_phase(total_phase);
                                   if s == 0.0 { continue; }

                                   let factor = mult as f64 * 2.0 * s;
                                   let c = cell.count as usize;
                                   cell.entries[c] = (out_qi, out_qj, factor);
                                   cell.count += 1;
                               }
                           }
                       }
                   }
               }
           }
           tables
       }
   }
   ```

   This reuses the existing `pauli_anticommutes`, `pauli_mul`, and `re_phase` functions
   directly. The build logic mirrors the existing cross-site Ladder kernel inner loop
   (lines 512–529 of lindblad.rs).

4. Write unit tests in `lindblad.rs::tests`:

   - `action_table_zz_lower_raise`: build tables, look up `(Lower, Raise, Z, Z)`.
     Assert count == 2, entries contain `(Y, Y, 4.0)` and `(X, X, 4.0)` (order
     doesn't matter — sort by first Pauli before comparing).

   - `action_table_zi_lower_raise`: look up `(Lower, Raise, Z, I)`.
     Assert count == 2, entries contain `(Y, Y, -2.0)` and `(X, X, -2.0)`.

   - `action_table_ii_all_dirs`: for all 4 direction pairs, look up `(_, _, I, I)`.
     Assert count == 0 for all.

   - `action_table_matches_inline_kernel`: for all 64 combinations, run the existing
     inline 2×2 loop (extracted into a helper function) and compare against the table.
     This is the exhaustive regression guard.

     Helper function for the reference implementation:
     ```rust
     fn reference_cross_site_action(
         left_dir: LadderDirection, right_dir: LadderDirection,
         p_qi: Pauli, p_qj: Pauli,
     ) -> Vec<(Pauli, Pauli, f64)> {
         // Copy the 2×2 loop from the existing cross-site kernel
         // Return non-zero (out_qi, out_qj, factor) entries
     }
     ```

**Review checklist**:
- [ ] `ActionTables::build()` reuses `pauli_anticommutes`, `pauli_mul`, `re_phase` — no
      duplicated algebra.
- [ ] `action_table_matches_inline_kernel` covers all 64 cells exhaustively.
- [ ] No changes to `apply()`, `LindbladOp::new`, or any hot path.
- [ ] `cargo test -p ppvm-timeevolve` clean.
- [ ] Benchmark: `bench_rhs` unchanged (record numbers in commit message).

---

## Task 31 — `CrossSiteLadder` structure and `LindbladOp::new` routing

**Goal**: Define `CrossSiteLadder`, populate it in `LindbladOp::new`. Old code path
stays active; the new struct is built but not yet consumed.

**Steps**:

1. In `src/lindblad.rs`, define the struct (no generics — it stores indices, not Pauli
   words):

   ```rust
   /// Fused cross-site ladder kernel data.
   pub(crate) struct CrossSiteLadder {
       /// (qubit_index, direction) per operator. Only Ladder ops are included.
       pub ops: Vec<(usize, LadderDirection)>,
       /// Rates sub-matrix for the Ladder ops. rates[i][j] is the rate for
       /// (ops[i], ops[j]). Size: ops.len() × ops.len().
       pub rates: Vec<Vec<f64>>,
       /// Precomputed action tables.
       pub tables: ActionTables,
   }
   ```

2. Add field to `LindbladOp`:

   ```rust
   pub struct LindbladOp<T: Config> {
       pub(crate) terms: Vec<LindbladTerm<T>>,
       pub(crate) cross_site: Option<CrossSiteLadder>,
   }
   ```

3. Modify `LindbladOp::new` (currently at line 141):

   **Before the main loop**: scan `ops` to collect Ladder metadata.
   ```rust
   // Collect Ladder op metadata: (original_index, qubit, direction).
   // Build a mapping from original op index → ladder sub-index (None for Generic).
   let mut ladder_meta: Vec<(usize, LadderDirection)> = Vec::new();
   let mut op_to_ladder: Vec<Option<usize>> = vec![None; ops.len()];
   for (idx, op) in ops.iter().enumerate() {
       if let JumpOp::Ladder(l) = op {
           op_to_ladder[idx] = Some(ladder_meta.len());
           ladder_meta.push((l.qubit, l.direction));
       }
   }
   let n_ladder = ladder_meta.len();
   let mut ladder_rates = vec![vec![0.0f64; n_ladder]; n_ladder];
   let mut has_cross_site = false;
   ```

   **Inside the `(i, j)` loop**: when both ops are Ladder and `qi != qj`, record the
   rate in `ladder_rates`:
   ```rust
   (JumpOp::Ladder(li), JumpOp::Ladder(lj)) => {
       // Existing: push LindbladTerm::Ladder to terms (keep old path active)
       terms.push(LindbladTerm::Ladder { ... });

       // New: also record rate for CrossSiteLadder
       if li.qubit != lj.qubit {
           let li_idx = op_to_ladder[i].unwrap();
           let lj_idx = op_to_ladder[j].unwrap();
           ladder_rates[li_idx][lj_idx] = gamma_ij;
           has_cross_site = true;
       }
   }
   ```

   **After the loop**: construct `CrossSiteLadder` if applicable.
   ```rust
   let cross_site = if has_cross_site {
       Some(CrossSiteLadder {
           ops: ladder_meta,
           rates: ladder_rates,
           tables: ActionTables::build(),
       })
   } else {
       None
   };
   LindbladOp { terms, cross_site }
   ```

4. Write unit tests:

   - `cross_site_populated_for_dense_rate_matrix`:
     ```rust
     let ops: Vec<JumpOp<SB>> = (0..3)
         .map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower }))
         .collect();
     let rates = RateMatrix::Dense(vec![vec![1.0; 3]; 3]);
     let lop = LindbladOp::new(ops, rates);
     let cs = lop.cross_site.as_ref().unwrap();
     assert_eq!(cs.ops.len(), 3);
     assert_eq!(cs.rates.len(), 3);
     assert_eq!(cs.rates[0].len(), 3);
     // Diagonal entries should be 0 (on-site is handled by terms)
     // Wait — actually diagonal entries ARE recorded because the loop processes all (i,j).
     // But we only set has_cross_site when qi != qj, and we only record rates for qi != qj.
     // The diagonal entries in ladder_rates stay 0.0. This is fine.
     assert_eq!(cs.rates[0][1], 1.0);
     assert_eq!(cs.rates[0][0], 0.0);
     ```

   - `cross_site_none_for_diagonal_rate_matrix`:
     ```rust
     let ops = ...; // 3 Ladder ops
     let rates = RateMatrix::Vector(vec![1.0; 3]);
     let lop = LindbladOp::new(ops, rates);
     assert!(lop.cross_site.is_none());
     ```

   - `cross_site_none_for_generic_ops`:
     ```rust
     // 2 Generic CollapseOps with dense rates → no Ladder ops → cross_site is None
     let lop = LindbladOp::new(generic_ops, dense_rates);
     assert!(lop.cross_site.is_none());
     ```

   - `cross_site_mixed_ops`:
     ```rust
     // 2 Ladder + 1 Generic, dense 3×3 rates
     // cross_site should have ops.len() == 2, rates is 2×2
     let cs = lop.cross_site.as_ref().unwrap();
     assert_eq!(cs.ops.len(), 2);
     assert_eq!(cs.rates.len(), 2);
     // Verify Ladder-Generic pairs are in terms as Generic
     let generic_count = lop.terms.iter()
         .filter(|t| matches!(t, LindbladTerm::Generic { .. }))
         .count();
     assert!(generic_count > 0);
     ```

   - All existing tests pass unchanged.

**Review checklist**:
- [ ] `CrossSiteLadder` stores a sub-matrix (Ladder ops only), not the full rate matrix.
- [ ] `op_to_ladder` mapping correctly handles mixed Generic/Ladder ops.
- [ ] Old code path (cross-site `LindbladTerm::Ladder` in `terms`) is still active.
- [ ] `cargo test -p ppvm-timeevolve` clean. All existing tests pass.
- [ ] Benchmark: `bench_rhs` — record numbers. May show slight regression from
      `CrossSiteLadder` construction overhead (but this is outside the hot loop).

---

## Task 32 — Fused cross-site apply kernel with weight filtering

**Goal**: Implement the fused cross-site kernel, wire `w_max` from the Strategy trait,
remove the old per-term cross-site path.

This is the largest task. Three sub-parts: (A) ppvm-runtime trait change, (B) `w_max`
threading, (C) fused kernel implementation.

### Step 1 — ppvm-runtime: `Strategy::max_weight()` + `PauliSum::strategy()`

In `crates/ppvm-runtime/src/traits/strategy.rs`:

```rust
pub trait Strategy: Default + Clone + Copy {
    fn capacity(&self, n_qubits: usize) -> usize;
    fn truncate<S, V, H, M, W>(&self, map: &mut M) where ...;

    /// Maximum Pauli weight retained by this strategy.
    /// Default: `usize::MAX` (no weight limit).
    fn max_weight(&self) -> usize { usize::MAX }
}
```

In `crates/ppvm-runtime/src/strategy.rs`, override for `MaxPauliWeight`:

```rust
impl Strategy for MaxPauliWeight {
    // ... existing capacity + truncate ...

    fn max_weight(&self) -> usize {
        self.0
    }
}
```

Override for `CombinedStrategy`:

```rust
impl<S1: Strategy, S2: Strategy> Strategy for CombinedStrategy<S1, S2> {
    // ... existing capacity + truncate ...

    fn max_weight(&self) -> usize {
        self.0.max_weight().min(self.1.max_weight())
    }
}
```

All other Strategy impls (`NoStrategy`, `CoefficientThreshold`, `Budget`,
`MaxLossWeight`) use the default `usize::MAX` — no changes needed.

In `crates/ppvm-runtime/src/sum/data.rs`, add a public getter:

```rust
impl<T: Config> PauliSum<T> {
    pub fn strategy(&self) -> T::Strategy {
        self.strategy
    }
}
```

### Step 2 — Thread `w_max` into `apply`

In `src/lindblad.rs`:

- Change `LindbladOp::apply` signature:
  ```rust
  pub(crate) fn apply(&self, p: &PauliSum<T>, result: &mut PauliSum<T>, w_max: usize)
  ```

- Change `apply_par` signature similarly.

- In `rhs_into` (line 634):
  ```rust
  let w_max = result.strategy().max_weight();
  ```
  Then pass `w_max` to `lindblad.apply(p, result, w_max)`.

- In `rhs_into_par`, same: read `w_max` from `result.strategy().max_weight()`, pass it
  to `apply_par`.

- Add `T::Strategy: Strategy` bound where needed (it should already be implied by the
  `Config` trait, but verify).

### Step 3 — Implement `apply_cross_site`

Add a new free function in `src/lindblad.rs` (below `apply_par`, around line 402):

```rust
/// Fused cross-site ladder kernel with weight-filtered generation.
///
/// Iterates over observable terms (outer loop) and qubit pairs (inner loop),
/// using precomputed action tables instead of runtime pauli_mul/re_phase.
/// The `w_max` parameter controls generation-time filtering: pairs that would
/// increase the observable term's weight beyond `w_max` are skipped.
#[inline]
fn apply_cross_site<T: Config>(
    cs: &CrossSiteLadder,
    p: &PauliSum<T>,
    result: &mut PauliSum<T>,
    w_max: usize,
) where
    // Same bounds as LindbladOp::apply
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>,
    T::Coeff: std::ops::AddAssign + Copy + std::ops::Mul<Output = T::Coeff>,
    T::PauliWordType: PauliWordTrait + Clone,
    f64: Into<T::Coeff>,
{
    let n_ops = cs.ops.len();

    for (w_a, coeff_a) in p.data().iter() {
        let weight = w_a.weight();
        let headroom = w_max.saturating_sub(weight);

        // Build active set: operator indices where w_a has non-identity.
        // Only needed for headroom == 0; for headroom >= 2 we process all pairs.
        match headroom {
            0 => {
                // At weight ceiling. Only pairs where both qubits are already
                // non-identity can contribute without exceeding w_max.
                let active: SmallVec<[usize; 32]> = (0..n_ops)
                    .filter(|&idx| w_a.get(cs.ops[idx].0) != Pauli::I)
                    .collect();
                for ai in 0..active.len() {
                    let i = active[ai];
                    for aj in 0..active.len() {
                        let j = active[aj];
                        if i == j { continue; }
                        process_pair(cs, w_a, coeff_a, i, j, result);
                    }
                }
            }
            1 => {
                // One unit of headroom. Skip pairs where BOTH qubits are identity.
                for i in 0..n_ops {
                    let i_active = w_a.get(cs.ops[i].0) != Pauli::I;
                    for j in 0..n_ops {
                        if i == j { continue; }
                        let j_active = w_a.get(cs.ops[j].0) != Pauli::I;
                        if !i_active && !j_active { continue; }
                        process_pair(cs, w_a, coeff_a, i, j, result);
                    }
                }
            }
            _ => {
                // Ample headroom: all pairs.
                for i in 0..n_ops {
                    for j in 0..n_ops {
                        if i == j { continue; }
                        process_pair(cs, w_a, coeff_a, i, j, result);
                    }
                }
            }
        }
    }
}

#[inline(always)]
fn process_pair<T: Config>(
    cs: &CrossSiteLadder,
    w_a: &T::PauliWordType,
    coeff_a: &T::Coeff,
    i: usize,
    j: usize,
    result: &mut PauliSum<T>,
) where
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>,
    T::Coeff: std::ops::AddAssign + Copy + std::ops::Mul<Output = T::Coeff>,
    T::PauliWordType: PauliWordTrait + Clone,
    f64: Into<T::Coeff>,
{
    let gamma = cs.rates[i][j];
    if gamma == 0.0 { return; }
    let qi = cs.ops[i].0;
    let qj = cs.ops[j].0;
    let left_dir = cs.ops[i].1.flip();
    let right_dir = cs.ops[j].1;
    let cell = cs.tables.get(left_dir, right_dir, w_a.get(qi), w_a.get(qj));
    for k in 0..cell.count as usize {
        let (out_qi, out_qj, factor) = cell.entries[k];
        *result += (w_a.set_new_2(qi, out_qi, qj, out_qj),
                    (gamma * factor).into() * *coeff_a);
    }
}
```

Note: this is the pre-symmetry version — full `i != j` iteration. Task 33 changes to
`i < j` with doubled weight.

Add `use smallvec::SmallVec;` to the imports (add `smallvec` to `Cargo.toml` if not
already a dependency — check first).

### Step 4 — Wire into `apply()` and remove old path

In `LindbladOp::apply` (line 428), after the existing `for term in &self.terms { ... }`
loop, add:

```rust
if let Some(cs) = &self.cross_site {
    apply_cross_site(cs, p, result, w_max);
}
```

Do the same in `apply_par` (line 282).

### Step 5 — Remove old cross-site terms from `LindbladOp::new`

In `LindbladOp::new`, in the `(JumpOp::Ladder(li), JumpOp::Ladder(lj))` arm:
- **Keep** the `LindbladTerm::Ladder` push for `li.qubit == lj.qubit` (on-site).
- **Remove** the `LindbladTerm::Ladder` push for `li.qubit != lj.qubit` (cross-site).
  These are now handled by `cross_site`.

```rust
(JumpOp::Ladder(li), JumpOp::Ladder(lj)) => {
    if li.qubit != lj.qubit {
        // Cross-site: handled by CrossSiteLadder. Record rate only.
        let li_idx = op_to_ladder[i].unwrap();
        let lj_idx = op_to_ladder[j].unwrap();
        ladder_rates[li_idx][lj_idx] = gamma_ij;
        has_cross_site = true;
    } else {
        // On-site: keep in terms as before.
        terms.push(LindbladTerm::Ladder {
            qi: li.qubit, qj: lj.qubit,
            left_dir: li.direction.flip(),
            right_dir: lj.direction,
            weight: gamma_ij,
        });
    }
}
```

### Step 6 — Update parallel dispatch

In `rhs_into` (line 640), the check `lindblad.terms.len() >= 200` will now undercount
work (cross-site terms removed from `terms`). For now, keep the threshold on
`terms.len()` for the existing `terms`-based parallelism. The fused kernel runs
sequentially. Task 34 addresses parallelism for the fused kernel.

### Unit tests

- `fused_kernel_matches_old_path_n3`: build a reference `LindbladOp` using the old path
  (wrap Ladder ops as `JumpOp::Generic(ladder.expand(...))` to force all-Generic terms).
  Compare `rhs()` output for a representative observable against the new fused-kernel path.
  Assert coefficient-level equality (within 1e-14).

- `fused_kernel_matches_old_path_n6`: superradiance n=6, full `solve()` comparison.
  Run with both paths, compare results within solver tolerance.

- `weight_filter_skips_correctly`: n=4, `MaxPauliWeight(2)` strategy. Initial observable:
  weight-1 terms only. After one `rhs()` call, verify no terms with weight > 2 exist.
  Compare against `rhs()` with `usize::MAX` (no filter) followed by manual weight
  truncation — should match.

- `weight_filter_headroom_1`: n=4, `MaxPauliWeight(3)`. Observable with weight-2 terms.
  Verify weight-3 terms appear (headroom allows +1) but no weight-4 terms.

- `w_max_usize_max_matches_unfiltered`: `CoefficientThreshold` strategy (max_weight
  returns `usize::MAX`). Verify fused kernel output matches old per-term kernel exactly.

**Review checklist**:
- [ ] `Strategy::max_weight()` added with default impl in ppvm-runtime trait.
- [ ] `MaxPauliWeight::max_weight()` returns `self.0`.
- [ ] `CombinedStrategy::max_weight()` returns `min(...)`.
- [ ] `PauliSum::strategy()` getter added in ppvm-runtime.
- [ ] `rhs_into` reads `w_max` from `result.strategy().max_weight()`.
- [ ] `apply()` and `apply_par()` receive and pass `w_max`.
- [ ] `apply_cross_site` uses `ActionTables::get` — no `pauli_mul` or `re_phase`.
- [ ] Weight filter: headroom 0 → active×active, headroom 1 → skip both-identity,
      headroom ≥2 → all pairs.
- [ ] `w_a.weight()` called once per observable term, not per pair.
- [ ] Cross-site `LindbladTerm::Ladder` entries no longer in `terms`.
- [ ] On-site `LindbladTerm::Ladder` entries still in `terms`.
- [ ] `fused_kernel_matches_old_path_n6` passes.
- [ ] `cargo test -p ppvm-timeevolve` and `cargo test -p ppvm-runtime` both clean.
- [ ] `cargo clippy -p ppvm-timeevolve -- -D warnings` clean.
- [ ] **Performance**: bench_rhs before/after. Record numbers.

---

## Task 33 — Symmetry exploitation (upper triangle)

**Goal**: Change `apply_cross_site` from full `i != j` iteration to upper-triangle
`i < j` with doubled weight.

**Steps**:

1. In `process_pair`, change:
   ```rust
   let gamma = cs.rates[i][j];
   ```
   to:
   ```rust
   let gamma = 2.0 * cs.rates[i][j];  // Merged (i,j) + (j,i); verified for all dir pairs
   ```

2. In `apply_cross_site`, change all three headroom branches from `i != j` to `i < j`:

   **headroom == 0**:
   ```rust
   for ai in 0..active.len() {
       for aj in (ai+1)..active.len() {
           process_pair(cs, w_a, coeff_a, active[ai], active[aj], result);
       }
   }
   ```

   **headroom == 1**:
   ```rust
   for i in 0..n_ops {
       let i_active = w_a.get(cs.ops[i].0) != Pauli::I;
       for j in (i+1)..n_ops {
           let j_active = w_a.get(cs.ops[j].0) != Pauli::I;
           if !i_active && !j_active { continue; }
           process_pair(cs, w_a, coeff_a, i, j, result);
       }
   }
   ```

   **headroom >= 2**:
   ```rust
   for i in 0..n_ops {
       for j in (i+1)..n_ops {
           process_pair(cs, w_a, coeff_a, i, j, result);
       }
   }
   ```

3. Add `symmetry_verified_all_64_cases` test:
   ```rust
   #[test]
   fn symmetry_verified_all_64_cases() {
       let tables = ActionTables::build();
       let dirs = [LadderDirection::Raise, LadderDirection::Lower];
       let paulis = [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z];

       for &di in &dirs {
           for &dj in &dirs {
               let left_dir_ij = di.flip();
               let right_dir_ij = dj;
               let left_dir_ji = dj.flip();
               let right_dir_ji = di;

               for &p_qi in &paulis {
                   for &p_qj in &paulis {
                       let cell_ij = tables.get(left_dir_ij, right_dir_ij, p_qi, p_qj);
                       let cell_ji = tables.get(left_dir_ji, right_dir_ji, p_qj, p_qi);

                       assert_eq!(cell_ij.count, cell_ji.count,
                           "count mismatch for dirs ({:?},{:?}), paulis ({:?},{:?})",
                           di, dj, p_qi, p_qj);

                       for k in 0..cell_ij.count as usize {
                           let (oqi, oqj, f) = cell_ij.entries[k];
                           let (oqj_r, oqi_r, f_r) = cell_ji.entries[k];
                           assert_eq!(oqi, oqi_r);
                           assert_eq!(oqj, oqj_r);
                           assert!((f - f_r).abs() < 1e-15);
                       }
                   }
               }
           }
       }
   }
   ```

   Note: this test assumes entries appear in the same order for (i,j) and (j,i).
   This is guaranteed because `ActionTables::build()` iterates the 2×2 loop in the
   same order regardless of input. If this assumption breaks, sort entries before
   comparing.

4. Add integration tests:

   - `symmetry_all_same_direction`: n=4 all-Lower, dense rates. Compare `rhs()` output
     of upper-triangle kernel against Task 32's full-iteration kernel (revert temporarily
     or use a reference function). Assert exact match.

   - `symmetry_mixed_directions`: n=4, 2 Raise + 2 Lower, dense rates. Same comparison.

**Review checklist**:
- [ ] All three headroom branches use `i < j` (or `ai < aj` for active set).
- [ ] `gamma = 2.0 * cs.rates[i][j]` in `process_pair`.
- [ ] `symmetry_verified_all_64_cases` passes.
- [ ] Integration tests match full-iteration path.
- [ ] `cargo test -p ppvm-timeevolve` clean.
- [ ] **Performance**: bench_rhs before/after. Expect ~2× improvement on cross-site
      kernel portion. Record numbers.

---

## Task 34 — Parallel dispatch for the fused kernel

**Goal**: Add rayon parallelism to `apply_cross_site` for large observables.

**Steps**:

1. Add `apply_cross_site_par` function, following the isolation pattern:

   ```rust
   #[cold]
   #[inline(never)]
   fn apply_cross_site_par<T: Config>(
       cs: &CrossSiteLadder,
       p: &PauliSum<T>,
       result: &mut PauliSum<T>,
       w_max: usize,
   ) where
       // Same bounds as apply_cross_site + Send + Sync
   {
       let n = p.n_qubits();
       let combined = p.data()
           .par_iter()
           .fold(
               || PauliSum::<T>::builder().n_qubits(n).build(),
               |mut local, (w_a, coeff_a)| {
                   // Same logic as apply_cross_site's inner loop
                   // (weight, headroom, match, process_pair into &mut local)
                   local
               },
           )
           .reduce(
               || PauliSum::<T>::builder().n_qubits(n).build(),
               |mut a, b| {
                   for (w, c) in b.data().iter() {
                       a += (w.clone(), *c);
                   }
                   a
               },
           );
       for (w, c) in combined.data().iter() {
           *result += (w.clone(), *c);
       }
   }
   ```

2. Add dispatch in `apply()` (or in the cross-site call site within `rhs_into`):

   ```rust
   if let Some(cs) = &self.cross_site {
       let estimated_work = p.len() * cs.ops.len() * cs.ops.len();
       if estimated_work >= CROSS_SITE_PAR_THRESHOLD {
           apply_cross_site_par(cs, p, result, w_max);
       } else {
           apply_cross_site(cs, p, result, w_max);
       }
   }
   ```

   Set `CROSS_SITE_PAR_THRESHOLD` to a conservative value (e.g. 50_000) and tune via
   benchmarking.

3. Ensure the sequential `apply_cross_site` has zero rayon code in scope (it already
   does if `apply_cross_site_par` is in a separate function with `#[cold]`).

**Unit tests**:
- `parallel_fused_matches_sequential`: n=8 all-Ladder, dense rates. Force both paths
  (by temporarily overriding threshold). Assert identical results.

**Review checklist**:
- [ ] `apply_cross_site_par` is `#[cold] #[inline(never)]`.
- [ ] Sequential path has zero rayon in scope.
- [ ] `parallel_fused_matches_sequential` passes.
- [ ] **Performance**: benchmark n=6 (sequential should win) and n=10 (parallel may win).
      Record crossover. Adjust `CROSS_SITE_PAR_THRESHOLD` accordingly.

---

## Task 35 — Benchmarks, integration tests, superradiance example

**Goal**: Comprehensive validation and performance measurement.

**Steps**:

1. Add Criterion benchmarks in `benches/step.rs` (or new file `benches/cross_site.rs`):

   - `bench_rhs_fused_n6`: n=6 dense rates, all-Ladder, `CoefficientThreshold` strategy
     (w_max = MAX). Measure `rhs()` call.
   - `bench_rhs_fused_n6_wmax3`: same but `CombinedStrategy<MaxPauliWeight(3),
     CoefficientThreshold>`. Measure filter benefit.
   - `bench_rhs_fused_n10_wmax4`: n=10, w_max=4. Larger system.
   - `bench_solve_superradiance_n6`: full `solve()` call, n=6.

2. Integration test in `src/lindblad.rs::tests` or `src/solve.rs::tests`:

   `solve_superradiance_fused_matches_baseline`: n=6 superradiance, compare solve results
   against pre-Phase-5 reference values (hardcoded from a known-good run, or against the
   Generic expand path).

3. Update `examples/superradiance_flame.rs`:
   - Change strategy from `CoefficientThreshold` to
     `CombinedStrategy<MaxPauliWeight, CoefficientThreshold>`.
   - The fused kernel activates automatically (Ladder ops + Dense rates).
   - Add a comment explaining the weight filter.

4. Record all benchmark numbers in the commit message with before/after comparison.

**Review checklist**:
- [ ] Criterion benchmarks present and runnable.
- [ ] Integration test passes.
- [ ] superradiance_flame example compiles and runs.
- [ ] Benchmark numbers in commit message.

---

## Task 36 — Optimise `PauliWordTrait::weight()` using `count_ones` (ppvm-runtime)

**Goal**: Replace per-bit iteration with byte-level popcount.

**Steps**:

1. In `crates/ppvm-runtime/src/word/data.rs`, find the `weight()` implementation
   (line 109):

   Current:
   ```rust
   fn weight(&self) -> usize {
       (0..self.nqubits)
           .filter(|&i| self.xbits[i] || self.zbits[i])
           .count()
   }
   ```

   Replace with:
   ```rust
   fn weight(&self) -> usize {
       self.xbits.as_raw_slice().iter()
           .zip(self.zbits.as_raw_slice().iter())
           .map(|(x, z)| (x | z).count_ones() as usize)
           .sum()
   }
   ```

   `as_raw_slice()` on a `BitArray<[u8; N]>` returns `&[u8]`. The `count_ones()` method
   on `u8` uses hardware popcount when available. Verify that `bitvec::BitArray` provides
   `as_raw_slice()` — if not, use `self.xbits.as_raw_mut_slice()` or access the inner
   array directly via `.data` or similar.

   **Caution**: if the storage has padding bits beyond `nqubits` (e.g. 6 qubits in 1 byte
   = 8 bits, 2 padding), the padding bits must be zero for `count_ones` to be correct.
   Verify this is guaranteed by `PauliWord` construction. If not, mask the last byte:
   ```rust
   let mut total: usize = 0;
   let raw_x = self.xbits.as_raw_slice();
   let raw_z = self.zbits.as_raw_slice();
   for (x, z) in raw_x.iter().zip(raw_z.iter()) {
       total += (x | z).count_ones() as usize;
   }
   total
   ```
   If padding bits might be set, subtract the excess:
   ```rust
   let bits = raw_x.len() * 8;
   let excess = bits - self.nqubits;
   // But this only works if the padding is at the end of the last byte.
   // Safer: just mask the last byte.
   ```
   Investigate the `bitvec` storage layout before implementing.

2. Add or update tests in `crates/ppvm-runtime/tests/` or in the existing test module:

   - `weight_all_identity`: `PauliWord::new(10)` (all I) → weight 0.
   - `weight_all_x`: set all 10 qubits to X → weight 10.
   - `weight_mixed`: set qubits 0,3,7 to X,Y,Z, rest I → weight 3.
   - `weight_boundary`: n=8 (exactly 1 byte), n=9 (2 bytes with padding). Verify
     padding bits don't inflate the count.

**Review checklist**:
- [ ] Uses `count_ones()` on integer types, not per-bit iteration.
- [ ] Padding bits handled correctly (no false positives for partial bytes).
- [ ] All existing ppvm-runtime tests pass.
- [ ] `cargo test -p ppvm-runtime` clean.
- [ ] **Performance**: time `weight()` for n=6, n=20, n=100 before/after (use a
      micro-benchmark or `#[bench]`). Record in commit message.
