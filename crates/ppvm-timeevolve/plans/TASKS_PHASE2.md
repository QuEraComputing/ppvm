# Phase 2 Tasks

## Task 16 — Commutator form in `LindbladOp::apply`

**Goal:** Replace the sandwich + anticommutator with a single product controlled by a
commutation-parity multiplicity, eliminating ~40% of MulAssign calls.

**Steps:**
1. In `src/lindblad.rs`, replace the inner loop body of `LindbladOp::apply` with the
   commutator form using `p1 + p2` multiplicity.
2. Remove `LindbladTerm::a_kl` (no longer needed after removing the anticommutator block).
3. Run benchmarks and record new mean times.

**Unit tests:**
- All existing `apply` and `rhs` tests must pass.
- `commutator_form_zero`: multiplicity=0 path produces no contribution.
- `commutator_form_double`: multiplicity=2 path gives coefficient 4×weight.

**Review checklist:**
- [ ] Anticommutator block removed.
- [ ] `multiplicity = p1 + p2`, not `p1 | p2`.
- [ ] No MulAssign when multiplicity=0.
- [ ] `LindbladTerm::a_kl` removed.
- [ ] All existing tests pass.
- [ ] Benchmark: `bench_rhs` before/after.

---

## Task 17 — Packed `comm_parity` using byte-level bit operations

**Goal:** Change `comm_parity` to operate byte-by-byte over packed storage instead of
per-qubit accessor loop.

**Steps:**
1. Change signature to take `&PauliWord<A, S>` directly.
2. Implement using raw byte fields.
3. Update call sites.

**Unit tests:** All existing tests pass; add n=20 and non-multiple-of-8 tests.

**Review checklist:**
- [ ] No new traits.
- [ ] Raw byte field used.
- [ ] Correct for NBYTES=1,2,3.
- [ ] Benchmark: `bench_rhs` before/after.

---

## Task 18 — Rayon parallelism over Lindblad terms

**Goal:** Parallelize outer `self.terms` loop using rayon fold/reduce.

**Steps:**
1. Add `rayon = "1"` to `[dependencies]`.
2. Rewrite `apply` with `par_iter().fold().reduce()`.
3. Benchmark and report speedup vs. core count.

**Unit tests:** All existing tests pass; `parallel_matches_sequential`.

**Review checklist:**
- [ ] `rayon` in `[dependencies]`.
- [ ] No write access to shared state in closures.
- [ ] Merge step correct.
- [ ] Benchmark: `bench_rhs` and `bench_solve` before/after.

---

## Task 19 — `Budget` truncation strategy

**Goal:** Add `Budget { target, min_threshold }` strategy capping `|P|` at target entries.

**Steps:**
1. Create `src/strategy.rs` with `Budget` implementing `Strategy`.
2. Re-export from `lib.rs`.
3. Update superradiance example.

**Unit tests:** `budget_limits_size`, `budget_matches_threshold`, `budget_accuracy`.

**Review checklist:**
- [ ] No changes to ppvm-runtime.
- [ ] When `|P| ≤ target`, identical to threshold-only truncation.
- [ ] Exported from `lib.rs`.
- [ ] Superradiance example updated.
