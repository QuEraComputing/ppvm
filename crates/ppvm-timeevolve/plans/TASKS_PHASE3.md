# Phase 3 Tasks

## Performance policy

Run `cargo bench bench_rhs` before and after every task and record both numbers in the
commit message. If a task introduces a regression (even a small one), document it explicitly:
note where it occurred and at which future task it is expected to be resolved. No regression
may be left unresolved at the end of the phase.

Benchmark fixture: `build_benchmark_lindblad` (n=5 dense all-Lower, 100 terms) in
`lindblad.rs::tests`. Use `--release` for all benchmarks.

## Task 21 — `LadderDirection`, `LadderOp`, `JumpOp`, `LadderOp::expand`

**Goal:** Introduce the public types needed for users to submit ladder operators. No changes
to `LindbladOp` yet — this task is purely additive.

**Steps:**
1. In `src/lindblad.rs`, add:
   - `pub enum LadderDirection { Raise, Lower }` with `pub fn flip(&self) -> Self` (swaps variants).
   - `pub struct LadderOp { pub qubit: usize, pub direction: LadderDirection }`.
   - `pub enum JumpOp<T: Config> { Generic(CollapseOp<T>), Ladder(LadderOp) }`.
2. Implement `LadderOp::expand<T: Config>(&self, n_qubits: usize) -> CollapseOp<T>`:
   - Create `word = T::PauliWordType::new(n_qubits)` (all-identity).
   - Build X word: `word.set_new(self.qubit, Pauli::X)`.
   - Build Y word: `word.set_new(self.qubit, Pauli::Y)`.
   - Phases: Lower → Y phase 1 (+i); Raise → Y phase 3 (-i). X always phase 0.
   - Push both terms (coeff 1.0 each) into a fresh `CollapseOp::new(n_qubits)`.
3. Re-export `LadderDirection`, `LadderOp`, `JumpOp` from `src/lib.rs`.

**Unit tests** (in `lindblad.rs::tests`):
- `ladder_direction_flip`: `Lower.flip() == Raise` and `Raise.flip() == Lower`.
- `ladder_op_expand_lower`: `LadderOp { qubit: 0, direction: Lower }.expand::<ByteF64<1>>(1)` produces a 2-term `CollapseOp` with X (phase 0) and Y (phase 1).
- `ladder_op_expand_raise`: same for Raise; Y term has phase 3.

**Review checklist:**
- [ ] `LadderDirection`, `LadderOp`, `JumpOp` exported from `lib.rs`.
- [ ] `expand` uses `T::PauliWordType::new` + `set_new` — no manual bit manipulation.
- [ ] No changes to `LindbladOp`, `LindbladTerm`, or `apply`.
- [ ] All existing tests pass.
- [ ] Benchmark: `bench_rhs` before/after; expect **no regression** (no hot-path changes).

---

## Task 22 — `LindbladTerm` enum + `LindbladOp::new` routing + on-site kernel

**Goal:** Change `LindbladTerm` to an enum, update `LindbladOp::new` to accept
`Vec<JumpOp<T>>`, and implement the fast on-site ladder kernel in `apply` (single-qubit
lookup for the `qi == qj` case).

**Steps:**
1. Change `LindbladTerm<T>` from a struct to an enum:
   ```rust
   pub(crate) enum LindbladTerm<T: Config> {
       Generic { left: PhasedPauliWord<...>, right: PhasedPauliWord<...>, weight: f64 },
       Ladder  { qi: usize, qj: usize,
                 left_dir: LadderDirection, right_dir: LadderDirection,
                 weight: f64 },
   }
   ```
   Update all existing field accesses in `apply`, `apply_par`, and tests to use
   `LindbladTerm::Generic { left, right, weight }` destructuring.

2. Update `LindbladOp::new` signature: `ops: Vec<CollapseOp<T>>` → `ops: Vec<JumpOp<T>>`.
   Routing for each `(i, j)` pair:
   - Both `Ladder(li)` and `Ladder(lj)`: push one `LindbladTerm::Ladder { qi: li.qubit, qj: lj.qubit, left_dir: li.direction.flip(), right_dir: lj.direction, weight: gamma_ij }`.
   - Otherwise: call `expand` on any `Ladder` variants to obtain a `CollapseOp`, then run the existing cross-product to produce `LindbladTerm::Generic` entries as before.
   - Add `impl<T: Config> From<CollapseOp<T>> for JumpOp<T>` so existing call sites can
     use `op.into()` instead of `JumpOp::Generic(op)`. Update all existing `LindbladOp::new`
     call sites in tests and the `build_benchmark_lindblad` fixture to pass
     `Vec<JumpOp<T>>` (using `.into()` or explicit wrapping). The existing test
     `two_term_op_x_plus_iy` asserts `lop.terms.len() == 4`; keep it valid by passing via
     `JumpOp::Generic`, and add a sibling assertion showing that the same op passed as
     `JumpOp::Ladder` yields `terms.len() == 1`.

3. In `apply`, add the match structure:
   ```rust
   match term {
       LindbladTerm::Generic { left, right, weight } => { /* existing code unchanged */ }
       LindbladTerm::Ladder { qi, qj, left_dir, right_dir, weight } if qi == qj => {
           for (w_a, coeff_a) in p.data().iter() {
               match w_a.get(*qi) {
                   Pauli::I => {}   // always zero
                   Pauli::Z => {    // both directions give +8γ·I and -8γ·Z
                       *result += (w_a.set_new(*qi, Pauli::I), 8.0 * weight * coeff_a);
                       *result += (w_a.clone(), -8.0 * weight * coeff_a);
                   }
                   Pauli::X | Pauli::Y => { /* sign from (left_dir, right_dir) */ }
               }
           }
       }
       LindbladTerm::Ladder { .. } => todo!(), // cross-site: Task 23
   }
   ```
   Derive the X/Y sign from Pauli algebra: for (Raise, Lower) — the standard lowering
   dissipator — X→+4γ and Y→-4γ; for (Lower, Raise) the signs flip. Derive all 4
   combinations and verify against `ladder_onsite_matches_generic`.

4. Update `apply_par` analogously inside the fold closure; cross-site case uses `todo!()`.

**Unit tests:**
- All existing tests pass (the `two_term_op_x_plus_iy` and `apply_lowering_op_lz` tests now
  exercise the `Ladder` path for the on-site case).
- `ladder_onsite_lower_all_paulis`: for a `JumpOp::Ladder(Lower, qubit=0)` with γ=1,
  verify L(I)=0, L(X)=+4X, L(Y)=-4Y, L(Z)=+8I-8Z against the previous `CollapseOp` result.
- `ladder_onsite_raise_all_paulis`: same for Raise; expect L(X)=-4X, L(Y)=+4Y, L(Z)=+8I-8Z.
- `ladder_onsite_matches_generic`: for n=1 and all 4 input Paulis, assert that the `Ladder`
  path and the `Generic` path (via `expand`) produce identical `PauliSum` outputs.

**Review checklist:**
- [ ] `LindbladTerm` is an enum; no struct definition remains.
- [ ] `LindbladOp::new` takes `Vec<JumpOp<T>>`.
- [ ] Ladder+Ladder pairs → one `LindbladTerm::Ladder`; all other pairs → `Generic`.
- [ ] `expand` is called at construction time for mixed pairs, not at apply time.
- [ ] `qi == qj` on-site case: no `comm_parity` call, no `MulAssign` — only `w_a.get()`, coefficient arithmetic, and `result +=`.
- [ ] Cross-site case has `todo!()` (not silently wrong).
- [ ] All existing tests pass.
- [ ] **Performance**: Run `bench_rhs` before and after. The enum dispatch adds a tiny
     match overhead to the `Generic` arm; if this causes a measurable regression on
     non-ladder benchmarks, record it in the commit message and note "to be resolved in
     Task 23 once the Ladder speedup offsets any overhead". For the all-ladder n=5 case,
     expect a clear speedup (on-site terms now O(1) vs O(n)).
- [ ] Benchmark numbers recorded in commit message.

---

## Task 23 — Cross-site kernel + `apply_par` update

**Goal:** Replace the `todo!()` with the two-qubit factorised kernel for the `qi != qj` case,
and update `apply_par` to handle both on-site and cross-site `Ladder` variants.

**Background:** When `qi != qj`, the left operator (direction `left_dir`) acts only on qubit
`qi` and the right operator (direction `right_dir`) acts only on qubit `qj`. The action is
the sum of 4 single-qubit sub-terms (left-X or left-Y) × (right-X or right-Y). Each
sub-term reduces to two single-qubit comm_parity checks + a single-qubit product lookup
+ `w_a.set_new_2(qi, out_qi, qj, out_qj)`.

**Steps:**
1. In `apply`, replace the `todo!()` cross-site branch with:
   - Extract the 2 sub-operators for left (X with phase 0, Y with phase from `left_dir`) and
     right (X with phase 0, Y with phase from `right_dir`).
   - For each of the 4 (left_sub, right_sub) combinations:
     a. Compute `p1`: does `left_sub.word` anticommute with `W[qi]`? Use `get_xbit`/`get_zbit`
        on the single qubit — no loop over all n qubits.
     b. Compute `p2`: does `W[qj]` anticommute with `right_sub.word`? Same.
     c. If `p1 + p2 > 0`: compute output Pauli at qi (`left_sub.word * W[qi]` → single-qubit
        multiply), output Pauli at qj (`W[qj] * right_sub.word`), combined phase, and
        `re_phase`. If non-zero, emit `w_a.set_new_2(qi, out_qi, qj, out_qj)` with coefficient.
2. Remove the `todo!()` from `apply_par`; apply the same pattern inside the fold closure.
3. Update the PAR_THRESHOLD check in `rhs_into` if needed (the term count is now lower
   for all-ladder systems; re-tune if benchmarks suggest a different crossover).

**Unit tests:**
- `ladder_crosssite_matches_generic_n2`: n=2, one `JumpOp::Ladder(Lower, qubit=0)` and one
  `JumpOp::Ladder(Lower, qubit=1)` with a 2×2 dense rate matrix. For all 16 two-qubit input
  Paulis (II, IX, IY, IZ, XI, …, ZZ), assert that the `Ladder` path output matches the
  `Generic` path (via `expand`).
- `parallel_matches_sequential_with_ladders`: adapt the existing `parallel_matches_sequential`
  fixture to use `JumpOp::Ladder` ops; assert sequential and parallel paths agree.

**Review checklist:**
- [ ] No loop over all n qubits inside the cross-site branch — only `get_xbit(qi)` / `get_zbit(qi)` calls.
- [ ] `set_new_2` used for cross-site output words.
- [ ] All 4 sub-term combinations handled (not just the contributing ones).
- [ ] `todo!()` removed from both `apply` and `apply_par`.
- [ ] `parallel_matches_sequential_with_ladders` passes.
- [ ] **Performance**: Run `bench_rhs` for the all-ladder n=5 case (sequential) and n=8 case
     (parallel). Cross-site kernel should deliver the bulk of the speedup. Verify that any
     regression recorded in Task 22 is now resolved. Record all numbers in commit message.
- [ ] If any non-ladder benchmark regressed in Task 22 due to enum dispatch overhead, confirm
     it has returned to baseline here.

---

## Task 24 — Update superradiance example to use `LadderOp`

**Goal:** Replace the manual `CollapseOp` X+iY construction in both examples with
`JumpOp::Ladder`, demonstrating the new API.

**Steps:**
1. In `examples/superradiance.rs` and `examples/probe_size.rs`, replace the manual
   `CollapseOp::new` + `op.push(ppw(...))` pattern with:
   ```rust
   ops.push(JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower }));
   ```
2. Update the `LindbladOp::new(ops, rates)` calls accordingly.
3. Verify both example outputs are unchanged.

**Unit tests:** None required — this is an example.

**Review checklist:**
- [ ] No manual `PhasedPauliWord` construction for ladder operators in the example.
- [ ] Example compiles and runs: `cargo run --example superradiance --release`.
- [ ] Output numerically matches the pre-Task-24 run.
- [ ] **Performance**: `bench_rhs` unchanged from Task 23 (example change is not in the hot path).
