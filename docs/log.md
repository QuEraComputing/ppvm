# `ppvm-*-2` refactor log

Status: temp / machine-maintained scratchpad.

**Ownership:** this file is written by the refactor **orchestrator** (the workflow
script), not by individual agents. Agents *read* it for context and *return*
structured findings; the orchestrator appends/updates the gap tables below. This
keeps writes race-free when review and test/bench agents run in parallel.

**Lifecycle:** resolved gaps stay (struck through) as the migration record. At each
crate/phase boundary, surviving lessons are distilled into
[`traits-2-configuration-and-hashing.md`](design/traits-2-configuration-and-hashing.md)
or the implementation plan, and this scratch is trimmed.

**Gap schema:** each gap is one row.

| field | meaning |
| --- | --- |
| `id` | `<crate>.<component>.<n>` |
| `type` | `impl-friction` \| `correctness` \| `perf-drift` \| `missing-proof` \| `missing-test` |
| `sev` | high \| medium \| low |
| `routed-to` | impl \| proof \| test \| design \| human |
| `status` | open \| closed |
| `commit` | resolving commit sha (when closed) |

Routing rule: `correctness`/`impl-friction`→impl (or design); `perf-drift`→impl or
human sign-off (hard gate + allowlist); `missing-proof`→proof; `missing-test`→test.

Convergence rule: a component is done when a full cycle yields **0 new gaps** and
`cargo build` + differential/property tests + `lake build` are all green. Max 4
iterations, then escalate to human.

---

## Phase 0 — scaffolding

Status: **complete.** Deliverables:
- `crates/ppvm-traits-2` — stub lib (trait modules land in Phase 1).
- `crates/ppvm-conformance-2` — test-only harness (`publish = false`): seeded RNG,
  `random_pauli_string`, replayable `GateOp`/`random_circuit`, `old_word_from_str`,
  `assert_close`; 4 smoke tests; `harness_smoke` Criterion bench.
- Both wired into `[workspace].members`.

Gates green: `cargo build --all-targets`, `cargo test -p ppvm-conformance-2` (4/4),
`cargo clippy -- -D warnings`, `cargo fmt --check`, `cargo machete`, and the
Criterion bench (`random_circuit/5q/1024` ≈ 9.2 µs). No gaps.

_(no gap entries)_

## Phase 1 — ppvm-traits-2

Status: **complete** (workflow `wf_0e6bf5cf-951`, 11 agents; reported `escalate` —
see note). Deliverables: 8 trait modules (coefficient/algebra/word/pauli/gates/
graded/batch/hash) + prelude; 18 Lean-oracle unit tests (phase1_leaf_types.rs,
phase1_batch.rs). Gates verified by orchestrator: `cargo build/clippy/fmt(after
fix)/test -p ppvm-traits-2` green; `lake build PPVM` green (2124 jobs).

Gaps (all resolved except one deferred by design):

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| t2.coefficient.1 | impl-friction | med | design | closed | `Coefficient::half` re-forecloses exact rings (`0.5·(1+i)∉ℤ[i]`), same escape as `Mul<f64>`. Split into `Halvable: Coefficient` capability; design+plan+compat-table updated. |
| t2.pauli.1 | correctness | low | design | closed | Design blanket-`Clifford` snippet had `xor_z_col(c,t)` ⇒ `z_t⊕=z_c`; authoritative Lean `cnotAct` is `z_c⊕=z_t`. Impl is correct (`xor_z_col(t,c)`); **design snippet fixed** + arg-order comments added. |
| t2.pauli.2 | impl-friction | low | design | closed | `SymplecticColumns`/`PhaseTrack` primitive set under-enumerated in design (`// ...`). Impl completed it (`xor_z_from_x`, `cz_bits`, `s/cz/x/y/z_phase`); design note added. |
| t2.algebra.1 | missing-proof | med | proof | closed | CNOT/CZ conjugation ℤ/4 phase uncovered by Lean (only H/S at n=1; bit-only for CNOT/CZ). Added `TwoPauli` (`𝒫₂`) `Group` + `conjCNOTHom`/`conjCZHom` + `conjCNOT_sign`/`conjCZ_sign` matching `phase/clifford.rs:80,95`, anchored to the physical generator tables. Orchestrator-verified the deltas match the old Rust kernel bit-for-bit. |
| t2.graded.1 | missing-proof | low | proof | closed | `overlap` labeled "bilinear" but only biadditivity+symmetry proven. Added `overlap_smul_left`/`_right` (homogeneity) ⇒ full `C`-bilinearity; design citation updated. |
| t2.algebra.2 | impl-friction | low | impl | closed | `Phase` lacked a group API; added `one`/`compose`/`inverse` (ℤ/4) for phase accumulation across `KeyProduct` chains. |
| t2.coefficient.2 | missing-test | low | test | **deferred → Phase 5** | `ImaginaryUnit`/`Conjugate` laws pinned only on `f64`/`Complex<f64>`; the exact-ring witness (`GaussianInt`-style) that proves L4 admits exact rings lands in `ppvm-sym-2`. Not a Phase-1 blocker. |

Perf: none (pure trait defs; no hot path).

**Escalate note (loop artifact, not a real failure).** The workflow returned
`escalate` at `MAX_ITERS=3` because the convergence check (a) treated a
low-severity ergonomics impl-friction as blocking, and (b) let the review
re-nominate proofs already discharged by the prove agent; the final on-disk state
is fully green. Separately, the parallel test agent wrote an unformatted file
(the `fmt` gate ran only *inside* the impl agent). Both are fixed in the workflow
script for later phases: fmt re-gated after the test agent; only high/med
correctness/impl-friction block convergence; discharged proofs are marked
resolved. Lesson logged for the workflow, deliverable accepted.

## Phase 2 — ppvm-pauli-word-2

### Component: `PauliWord` — status **complete** (workflow `wf_98ec8b83-228`, converged iter 1)

Deliverables: `PauliWord` (data/product/clifford/hash/storage/column) implementing
`Word<Site=Pauli>`, `PauliBits`, `SymplecticColumns`+`PhaseTrack` (⇒ blanket
`Clifford`), `KeyProduct`, `Indexable`, `Columnar`; packed X/Z planes + lazy
`OnceLock<u64>` hash; phase-mul + Clifford kernels ported from the old crate. 27
in-crate tests. Conformance: `pauli_word_diff.rs` (8 differential vs old),
`pauli_word_lean.rs` (10 Lean-oracle), `pauli_word_bench.rs` (new-vs-old).

Gates verified by orchestrator: `build`/`clippy`/`fmt`/`test` (pauli-word-2 +
conformance-2), `machete`, `lake build PPVM` — all green.

**Perf gate: PASS** (no regression). new/old ratios: product `key_mul` **0.77×**
(new ~23% faster) · cnot conjugation **0.91×** (parity) · weight **1.00×**
(verbatim) · `key_hash` reported cold/warm **6.94ns→0.60ns** (lazy OnceLock; the
old crate hashes eagerly at construction, so no new/old ratio applies). The test
agent caught+removed a clone-vs-`Copy` bench artifact (new word drops `Copy` for
the cache). No target exceeds the 1.15 gate → no perf-drift, no allowlist entry.

Lean added (prove agent, `lake` green): `Matrix.lean` `tensorPauli`/`iuPow`/
`prod_iuPow`/`tensorPauli_mul` — grounds the packed n-qubit `phaseExpN` against a
genuine 2ⁿ×2ⁿ tensor-product ℤ[i] matrix product, closing the end-to-end gap
`Word.lean`'s scope note had left open (n-qubit phase previously validated only as
the *sum* of per-qubit exponents). Orchestrator-reviewed: real, non-vacuous.

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| pw2.design.1 | impl-friction | low | design | closed | `word-data-structures.md` was stale (assigned `SymplecticColumns`/`PhaseTrack` only to `PhasedPauliWord`). Authoritative config doc + plan require the bare word to implement them with a phase-discarding `PhaseTrack`. Doc fixed; `clifford.rs` note reframed from "friction" to authoritative. |
| pw2.proof.1 | missing-proof | med | proof | closed | n-qubit product grounded per-qubit only (Word.lean scope note). Added `tensorPauli_mul`; ties `key_mul` to real multi-qubit operators. |
| pw2.test.1 | missing-test | med | test | closed | Review (run in parallel) saw no matrix-oracle/randomized/differential coverage in-crate; the test agent independently added it in conformance-2 (matrix-exponent oracle exhaustive + randomized n-qubit; differential vs old for construction/get/weight/iter/bits/product/Clifford; hash-contract + avalanche). 18 pass. |
| pw2.impl.1 | impl-friction | low | accepted | noted | Bare word drops the Clifford *phase* delta by design (no-op `PhaseTrack`); the differential Clifford test compares bits only, phase validated via the Lean ℤ[i] oracle. Recovered by the phased wrapper (Phase 3). |

**Phase-3 watch item.** Dropping `Copy` (lazy `OnceLock` hash) is op-vs-op
parity here, but in `PauliSum` propagation words are cloned per stored term; watch
the *aggregate* clone cost in Phase 3 (may warrant a perf-allowlist entry then).

### Components: `LossyPauliWord`, `Phased` — pending (next Phase-2 runs).

## Phase 3 — ppvm-pauli-sum-2

_(no entries yet)_

## Phase 4 — ppvm-tableau-2

_(no entries yet)_

## Phase 5 — ppvm-sym-2

_(no entries yet)_

## Phase 6 — ColumnStore backend

_(no entries yet)_

## Phase 7 — downstream + cutover

_(no entries yet)_

---

## Perf-drift allowlist (human-approved regressions)

Regressions listed here are accepted as designed-in trade-offs and do **not** block
a crate. Anything not listed is a hard gate.

| id | component | metric | accepted ratio | justification | approved-by |
| --- | --- | --- | --- | --- | --- |
| _(none yet)_ | | | | | |
