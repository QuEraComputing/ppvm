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

## ▶ Continue here — handoff (2026-08-09)

**RNG ownership migration — done.** The `-2` tableau no longer embeds a
`SmallRng`. Every stochastic trait method now takes `rng: &mut R`, the simulators
are pure state (`new` is deterministic, `fork()` is `clone()`, no
`new_with_seed`), and the *frontends* own a generator on the user's behalf:
`ppvm-python-native`'s `#[pyclass]` wrappers and `ppvm-vihaco`'s executors each
hold one `SmallRng`, threaded through a cfg'd `draw!` macro that keeps the
interface/dispatch tables identical across the legacy and `-2` backends.
`ppvm-stim` offers both shapes (`execute`/`sample*` self-seed; `*_with_rng` take
the caller's). Rationale and the layer-by-layer rule are in
[`traits-2-implementation-plan.md`](design/traits-2-implementation-plan.md)
§ *Where the randomness lives*.

Seeded behaviour is byte-identical to old: the draw order is preserved, and the
draws that used to seed a per-branch tableau are retained as explicit burns
(`GeneralizedTableauMixture::burn_legacy_tableau_seeds`). Pinned by
`ppvm-vihaco`'s `seeded_branch_fixture_matches_cross_backend_snapshot` and the
Python `PPVM_EXPECT_BACKEND` matrix (238 tests × both backends). Green:
`cargo test --workspace`; stim/vihaco/tui/cli under `traits-2` and
`traits-2,rayon`; the `traits-2` wasm32 library checks; `cargo fmt`;
`cargo clippy --workspace -- -D warnings` plus `--all-targets` on every `-2`
crate; `cargo-machete`; `lake build`. (Pre-existing and untouched: `--all-targets`
Clippy and the wasm workspace build both fail inside the *old* `ppvm-tableau` /
`ppvm-tableau-sum` / `ppvm-pauli-sum` crates.)

Phases 0–7 and the post-baseline performance optimizations are committed through
`73580afa` (including `0fc0c57e`); the old reference crates remain untouched.
Final targeted
full-duration confirmation reran every actionable `36850446` row in at least
four independent processes and hash/layout-sensitive rows in eight. Of the 75
original above-gate rows, 14 are fixed, 11 are parity, none is actionable, and
50 are evidence-adjudicated atomic-semantics, duplicate-path, identical/no-op,
or representation/layout controls. The five grouped adjacent rows are three
fixed, one parity, one non-actionable, and none robust. The performance cutover
blocker is closed; destructive cutover still awaits maintainer approval and
review (`docs/performance-report.md`). All
benchmark modes, conformance/workspace tests, formatting, strict optimized-target
Clippy, and Lean are green.

### ⚠ MEASUREMENT WARNING — read before touching any perf number here
Two instruments previously trusted in this log are now known to be unsound:

1. **`examples/trotter_attrib.rs` systematically understates the gap.** Its
   per-op `Instant::now()` pairs are optimization barriers, and they penalize
   the *old* engine more (its ops are cheaper and more inlinable). Measured on
   the same tree: attrib reported total `1.08×` while criterion reported
   `1.35×` for the identical code. Its per-op *ratios* are diluted toward 1.0.
   Use it for coarse *ordering* of op classes only — never for a ratio, and
   never to declare parity.
2. **The new engine is bimodal across processes** (see `ps2.rekey.bimodal`).
   A single criterion run can land in either mode. **Any perf claim here needs
   ≥4 separate process runs with the modes reported**, not a mean and a CI —
   criterion's own CI is tight *within* a mode and will look convincing while
   being off by ~1.7×.

The sound instrument today is `benches/pauli_sum_integration.rs`
(`rekey_cnot` / `truncate` / `integration_trotter`), run ≥4 times.

### ★ PRIME DIRECTIVE (do not violate) — behaviour-preserving refactor
This refactor **only cleans up and aligns the abstraction with the math**. It must
**NOT change any user-facing behaviour** vs the old crate. *Any* observable
behavioural difference — not just numeric outputs, but **when side effects happen**
(e.g. when truncation runs), API contracts, defaults, error/edge semantics — is a
**GAP** to be fixed (restore old behaviour), not accepted. The workflow must
surface these (see the workflow's new behaviour-parity check). When old itself is
provably wrong, the Lean oracle governs and the divergence is documented — that is
the *only* allowed behaviour change.

### 1 — `ps2.truncate.behaviour` — ✅ CLOSED (`d4aa7875`)
The behaviour gap is fixed: the three internal `self.policy.truncate(...)` calls
are gone from `apply`/`rekey_bijective`/`rotate_in_place`, so truncation is
**explicit-only** via `Sum::truncate()`, matching old's deferred semantics. New
regression guard `tests/pauli_sum_truncation_behaviour_diff.rs` (3 differential
tests). **Discriminating power verified**: 2 of the 3 fail against the pre-fix
engine with exactly the predicted divergence; the third is the eager-truncation
control that passes either way. The whole pre-existing suite passed both before
and after — every existing test truncates right after each gate, the single
schedule where the two semantics agree. That was the coverage hole.

### ✅ PERF: end-to-end now **~1.09×**, inside the 1.15× gate
Same-build ratios, 4 separate process runs, no bimodality observed in any:

| target | ratio | was |
| --- | --- | --- |
| `rx` sweep | **0.98×** (new faster) | — |
| `truncate` | **0.99×** | 0.99× |
| `pauli_error` sweep | **1.06×** | **2.03×** |
| `rekey_cnot` sweep | **1.11×** | 1.71× |
| **`integration_trotter`** | **~1.09×** | **~1.47×** |

Two fixes got there, and **both were behaviour fixes that happened to be the
perf fixes** — the divergence and the slowdown had the same root each time:
`23860055` (re-key `entry()`→`insert`) and `e3e1f9cc` (`ScaleByKey` mutates in
place instead of `retain`-ing). Note old's *absolute* numbers drift between
builds (701–788 µs for untouched code) — only same-build ratios mean anything.

### ⚑ NEW behaviour-parity gaps found by the Baseline agent (verified by hand)
Three divergences from old, all **verified directly against the old sources**.
Under the PRIME DIRECTIVE each is a gap.

1. **`ps2.zero.behaviour` (high) — PARTLY FIXED (`e3e1f9cc`).** The old crate has **no `reduce`** and no
   drop-zero logic *anywhere* — a zero-coefficient entry stays in the support
   forever, and old's `PartialEq` compares maps exactly, so
   `ppvm-pauli-sum/tests/loss.rs::test_reset_channel` depends on `*= 0.0`
   keeping every key. The new `reduce()` (`ppvm-traits-2/src/containers.rs`)
   **and** `ScaleByKey` both drop exact zeros, so `len()`/`contains` diverge —
   e.g. after `pauli_error(q, [0.0, 0.25, 0.25])` (λ_X = 0).
   **This indicted an earlier decision in this very log:** the `ps2.noise.1` row
   records the orchestrator *adding* zero-dropping to `ScaleByKey` to uphold a
   "reduced-canonical-form invariant" that is the new design's invention, not
   old's behaviour — so that fix introduced a divergence rather than removing
   one. **Now reverted** (`e3e1f9cc`): `ScaleByKey` takes old's `Fn(&K, &mut C)`
   shape and walks with `iter_mut`, so it *cannot* remove — correct by
   construction rather than by a zero-check. That also took `pauli_error` from
   2.03× to 1.06×, since `retain` had been dragging hashbrown's erase machinery
   through the hot walk. Guarded by
   `zero_channel_eigenvalue_keeps_the_term_exactly_as_old_does`, which compares
   EXACT supports (the existing suite's `FLOOR` filter masked this) and is
   verified to fail against the pre-fix engine.
   **STILL OPEN:** `reduce()` also drops exact zeros on the `apply`/`from_terms`
   paths, so the gap is not fully closed.
2. **`ps2.default.threshold` (medium).** `CoefficientThreshold::default()` is
   `1e-12` in old (`ppvm-pauli-sum/src/strategy.rs:113`) but `0.0` in new
   (`ppvm-pauli-sum-2/src/policy.rs`, a derived `Default` on the `f64` field) —
   a silently changed user-facing default.
3. **`ps2.preserve.missing` (medium).** `preserve_strings` (builder option +
   the snapshot/restore post-filter inside old `PauliSum::truncate`, pinned by
   `ppvm-pauli-sum/tests/preserve.rs`) has no counterpart in the new crate.

Plus a **suspected old bug** for the L4 `Multiply` component to adjudicate:
`impl MulAssign<PauliSum<T>> for PauliSum<T>` (`ppvm-pauli-sum/src/sum/ops.rs:70`)
loops the rhs terms calling `self.map_add(..)`, and `map_add` *replaces* the
support with its image — so `A *= (b0·P0 + b1·P1)` computes the product chain
`A·b0P0·b1P1` instead of `A·b0P0 + A·b1P1`. Untested in old. The new `Multiply`
must accumulate into a fresh accumulator (`twistedConv`, `Twisted.lean`).

### 2 — `ps2.cnot.rekey.perf` — ✅ CLOSED (`23860055`)
End-to-end Trotter is **~1.28×** (was ~1.47× at session start; criterion, same
build, `[u8;8]`). Numerics match. The re-key itself is now **1.11×**.

**What was established this session (with the new sound instrument):**
- **`truncate` is at parity** — 237 ns new vs 240 ns old (`0.99×`) on a 372-term
  support. This *disproves* the standing suspicion that `CombinedPolicy`'s two
  retain scans were a contributor. Do not re-open it without new evidence.
- **A real win landed** (`b68b7e8e`): `set_x_bit`/`set_z_bit` invalidated the
  hash cache unconditionally, so every gate re-hashed nearly the whole support
  to recompute a bit-identical digest (the Clifford kernels write both target
  bits unconditionally, and most terms are `I` at the gate's qubits). Now
  content-conditional. `cnot` sweep 90.5→77 µs (good mode), 214→134 µs (bad
  mode).
- **Ruled out with evidence:** `aux.reserve` (no effect); `drain()` vs old's
  `iter()`+clone (~3% — not worth abandoning the move semantics); the extra
  hash pass / warm trick (helps only while the needless re-hash exists, then
  becomes a no-op).
- Old `PauliWord` is `Copy`; the new one is not (`AtomicU64` cache) — the
  Phase-2 watch item `pw2` flagged. Not yet measured as a cause.

**ROOT-CAUSED AND LARGELY FIXED (`23860055`).** The cause was
`entry(k).and_modify(..).or_insert(..)` in `RekeyBijective for HashMapStore`.
It aggregated a collision that **cannot occur** (`f` is injective by contract;
for a Clifford that is the symplectic bijection already proven `*_bijective` in
`Symplectic.lean`), and on this crate's monomorphization the `entry` chain
compiled **out of line** — a profile put 41% of the sweep's samples in a
standalone `hashbrown::rustc_entry` frame, where old's identical source
construct inlines fully into `map_add_assign`. Every other per-term cost was
being paid *through* that out-of-line call. Replaced with a plain `insert` +
a `debug_assert!` that nothing is displaced.

Measured over **8 interleaved process pairs** (prebuilt binaries, alternating,
spread <1.5% within each variant):

| target | before | after | vs old |
| --- | --- | --- | --- |
| `rekey_cnot` | 76.5–77.3 µs | **58.9–59.4 µs** | 1.45× → **1.11×** |
| `integration_trotter` | ~1.47× | **~1.28×** (902–915 vs 701–720 µs) | open |
| `truncate` | 0.99× | 0.99× | parity |

Residual ~1.11× on the re-key is unattributed. Reproduce:
`cargo bench -p ppvm-conformance-2 --bench pauli_sum_integration -- rekey_cnot`
(≥4 runs — see the measurement warning).

### 3 — STILL OPEN: `ps2.rekey.bimodal`
The re-key is **bimodal across processes** — for the identical binary and input.
The `insert` fix improved **both** modes (fast 77→59 µs, slow 134→128 µs) but did
**not** remove the bimodality: 3 of 6 `cargo bench` launches still land slow.

Evidence it is process-level state fixed at startup, not drift or noise:
measuring the new engine **twice in the same process** gives identical numbers to
0.5% (132.10/132.61, 71.16/70.98, 131.95/132.04), while separate processes select
a mode. The **old engine never shows it** (stable ~53 µs in every run).

**New evidence — the trigger is the launch environment.** 16 consecutive
*direct-binary* launches were all fast (58.9–59.4 µs); `cargo bench` launches of
*the same binary* are bimodal (3/6 slow). Mechanism deliberately left
**unattributed** — do not assert one without demonstrating it.

Ruled out with evidence by the workflow run (wf_e473a67e-42a): heap layout of the
two ~21 KB tables (an allocator-hook probe found the table addresses identical
mod 16 KB, delta `0x6000`, in all 22 observed processes); process env and cwd;
CPU core type / scheduling (interleaving old/new *inside* a slow process kept old
fast and new slow throughout); code-ASLR slide (no correlation with the binary
base mod 64K/128K/256K/1M over 16 runs); and a different hot path in the slow
mode (the two profiles are the same instruction mix, uniformly dilated).

**REFUTED — the `Copy`/atomic hypothesis.** Swapping the word's `AtomicU64` hash
cache for `Cell<u64>` changed nothing at all (76.2 fast / 114.0 slow, bimodality
intact). The Phase-2 `pw2` watch item is **not** a cause of either gap. Also,
`#[cold] #[inline(never)]` on the hash-miss path made it much worse
(76→116 fast, 114→176 slow) — that path is hot, not cold.

⚠ That workflow also claimed "`insert()` → 0 of 34 processes slow". It **did not
replicate** (3/6 slow here) and is recorded as unconfirmed; its measurements were
taken on a 65-file tree that had also edited the old reference crates.

### 2 — THEN: Phase 3 component 3 — L4 `Multiply` (operator product), not started
Needs `Complex`/`ImaginaryUnit`; Lean oracle `lean/PPVM/**` `Twisted.lean` (`tmul_assoc`). Run it through the workflow (below).

### Workflow — use it for behavioral components
`.claude/workflows/traits-2-component.js` carries a **★ PRIME DIRECTIVE** (in `COMMON`, seen by every agent): behaviour-preserving refactor only — any user-facing behaviour change (incl. *when* side effects fire) is a gap. Phases: **Baseline** (mine old-crate real workloads → integration acceptance bar + perf-critical architecture features + **user-facing behavioural contracts**) → Implement (preserve features + behaviour) → Verify (Review incl. **architecture-parity** + **behaviour-parity** ∥ Test with **integration-first** perf gate + **behaviour-parity diff tests**) → Prove (Lean adjudicates suspected old-impl bugs). The behaviour-parity dimension is what would have caught the internal auto-truncate.

### Non-negotiable discipline (hard-won this session — don't relearn it)
- **Behaviour preservation is the prime directive** (see banner above). The refactor changes *structure*, never *observable behaviour*. Before/after any component, ask: does the new engine behave identically to old for a user — same outputs, same *timing of side effects* (truncation, reduction), same defaults/contracts? A divergence is a gap. (The internal auto-truncate `ps2.truncate.behaviour` is the live example.)
- **Perf: same-build new/old ratios ONLY.** Cross-build absolutes swing from code-alignment/i-cache (Mytkowicz layout bias — e.g. an untouched `old/rx` moved 4.5↔5.9µs). Use an interleaved one-process harness for A/B (see `examples/trotter_attrib.rs`). Fair config (identical storage/coeff/hasher — build a bench-local `[u8;8]` new type). Attribute any drift by controlled A/B or profile; if you can't isolate it, say "unattributed" — never assert a mechanism (I mis-called the rotation perf 3× before profiling settled it).
- **Integration-first.** Single-gate microbenches structurally miss cumulative per-gate costs (allocation/truncation) because a tight one-gate loop recycles one warm allocator page. Always gate on the end-to-end deep-circuit workload vs old.
- **Architecture parity.** When generalizing, preserve the old impl's perf-critical structural features — they may *relocate* (the aux double-buffer moved into the `HashMapStore` storage backend so `Sum` stays generic) but must not *vanish*. Check for this in review.
- **Numerics:** golden-master targets old-parity unless the Lean oracle says old is wrong.

### Key files
- New engine: `crates/ppvm-pauli-sum-2/src/{sum,store,clifford,rotation,noise,policy}.rs`. Storage = `HashMapStore` newtype `{primary, aux, scratch}` in `store.rs` (graded/pure-diagonal traits delegate to `primary`; buffered `RekeyBijective`/`RotateInPlace` use aux/scratch).
- Old reference: `crates/ppvm-pauli-sum/src/**` (real workloads: `benches/trotter.rs`, `benches/random-circuit.rs`, `tests/trotter.rs`, `tests/ghz.rs`).
- Conformance harness: `crates/ppvm-conformance-2/{src/lib.rs, tests/*, benches/*}` — integration diff+bench are `pauli_sum_integration*`; microbenches `pauli_sum_bench.rs`; per-op attribution tool `examples/trotter_attrib.rs`.
- Lean: `lean/PPVM/**` (`lake build PPVM` from `lean/`).
- Open/closed gaps: ledger rows in the Phase 3.2 section — OPEN: `ps2.rekey.bimodal` (process-dependent swing, narrowed), `ps2.zero.behaviour` (**high**, behaviour), `ps2.default.threshold`, `ps2.preserve.missing`, `ps2.oldbug.mulassign` (suspected old bug), `ps2.attrib.instrument`. Closed: `ps2.integration.1`, `ps2.store.aux`, `ps2.rot.perf`, `ps2.trotter.perf`, `ps2.truncate.behaviour`, `ps2.cnot.rekey.perf`.

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

### Component: `LossyPauliWord` — status **complete** (workflow `wf_77a58611-48f`, converged iter 2)

**Moved to its own crate `ppvm-lossy-pauli-word-2`** (maintainer decision: a lossy
word is a distinct concrete Pauli-word impl, not a submodule of `PauliWord`). It
depends on `ppvm-pauli-word-2` and reuses its `PauliStorage`/`HashFinalize`. The
first (mis-targeted) run into `ppvm-pauli-word-2` was stopped and reverted.

Deliverables: `LossyPauliWord` (data/clifford/hash/column) — packed X/Z/loss
planes, two component relaxed-atomic caches (loss-only mutation skips X/Z
rehash); `Word<Site=LossySite<Pauli>>`, `PauliBits`(`is_lost`),
`SymplecticColumns`+`PhaseTrack` (loss-guarded, phase-discarding), `Indexable`,
`Columnar`; inherent `set_lost`/`clear_loss`/`loss_weight`. 27 in-crate tests.
Conformance: `lossy_pauli_word_diff.rs` (9), `lossy_pauli_word_lean.rs` (9),
lossy bench. Gates verified by orchestrator (build/clippy/fmt/test across all four
`-2` crates, machete, `lake build PPVM`) — all green; the `PauliBits` relaxation
did not regress `ppvm-pauli-word-2` (27 tests still pass).

**Perf gate: 6/6 pass (1 flag investigated → benign).** product **0.79×**, cnot
**0.78×**, key_hash-warm **1.02×**, weight **1.02×**; key_hash-cold 15.6ns is the
design-accepted lazy-cache trade-off. `loss_weight` flagged at 1.28× at single-
`u64` width — **investigated and cleared**: the body is a verbatim port over
identical `BitArray<u64>` storage, and at `[u64; 8]` width (8× the popcount work)
new/old converge to **1.00×** (0.915ns each), as does `weight`. The single-word
delta is a nanobenchmark per-call-overhead artifact (code alignment), proven by
new `loss_weight` (1 plane) measuring slower than new `weight` (3 planes) — an
impossible ordering for the computation itself. No structural difference; not
allowlisted as a regression. (bench note added; agent correctly flagged not
self-accepted.)

Lean added (prove agent, `lake` green): `Symplectic.lean` loss-guarded Clifford
model — `cnotActL`/`czActL`, `LossInv`, per-primitive `xorXColL`/`xorZColL` with
`*_preserves_loss`, and `xorZColL_xorXColL_eq_cnotActL` (per-primitive guarding
composes to the atomic whole-gate skip — the exact property the crate's
`clifford.rs` guard relies on). Orchestrator-reviewed: real, non-vacuous.

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| lpw2.traits.1 | impl-friction | high | design | closed | `PauliBits: Word<Site=Pauli>` unsatisfiable by `LossyPauliWord` (Site=`LossySite<Pauli>`). Relaxed supertrait to `PauliBits: Word` in ppvm-traits-2; design doc + plan updated; `PauliWord`/conformance still build. |
| lpw2.clifford.1 | impl-friction | med | design | closed | Lossy `SymplecticColumns` must guard on loss (lost qubit = Clifford no-op), so it is a loss-*guarded* Sp map, not the literal pure map. Guard put in each primitive; design (`word-data-structures.md`) note added; machine-checked (`Symplectic.lean`). |
| lpw2.proof.1 | missing-proof | med | proof | closed | Loss-guarded Clifford had no Lean counterpart. Added the loss model + invariant-preservation + primitive-composition proofs. |
| lpw2.perf.1 | perf-drift | — | human | closed (investigated) | `loss_weight` flagged 1.28× at `u64` width; investigated → **benign** (verbatim port; converges to 1.00× at `[u64;8]`; single-word delta is a nanobench alignment artifact). Not a regression. |
| lpw2.keyproduct.1 | impl-friction | low | accepted | noted | No `KeyProduct` for `LossyPauliWord` (loss breaks the twisted-product group; `iᵏ(v⊕w)` undefined once a factor is Lost). Left inherent-only, matching the old crate. Lossy word is `Indexable`/`Columnar` but not an L4 key-product participant. |
| lpw2.test.1 | missing-test | low | test | closed-by-conformance | In-crate Clifford tests cover a subset; the conformance differential replay (random circuits) + Lean-oracle generator tables + the new `Symplectic.lean` composition proof give full coverage. |
| lpw2.design.1 | correctness | low | design | closed | `Ord`/serde absent on `-2` words (as on the sibling `PauliWord`); added a first-prototype scope note in `word-data-structures.md`. |

Workflow refinement (this run): perf-drift is now a **hard gate** in the loop —
any over-threshold regression escalates to the human regardless of the agent's
severity label (was: only high/med blocked).

### Component: `Phased` / `PhasedPauliWord` — status **complete** (workflow `wf_548f135a-90a`, converged iter 1; + a maintainer-chosen perf fix)

**Its own crate `ppvm-phased-pauli-word-2`** (maintainer decision). Generic
`Phased<W>` = base word + explicit ℤ₄ phase; `PhasedPauliWord = Phased<PauliWord>`.
Real ℤ₄ phase tracking (recovers the conjugation sign the bare/lossy words drop),
`Word` delegated to `W`, phased product reuses `W`'s `KeyProduct` (no phaseExp
duplication), **non-indexable**. 14 in-crate tests; conformance
`phased_pauli_word_diff.rs` (6, diffs the *phase* too) + `phased_pauli_word_lean.rs`
(11) + bench. Gates verified by orchestrator (build/clippy/fmt/test = 149 across
5 crates, machete, `lake`) — all green.

**Perf: product 1.00×, phased cnot 0.85× (now *faster* than old).** Initially the
blanket `Clifford` (separate `PhaseTrack` + `SymplecticColumns` steps) read the
inner bits twice → phased cnot **1.84×**. Investigated (confirmed real, not a
nanobench artifact: the old kernel is fused, reads bits once). **Maintainer chose
to fix** (over allowlist): added an opt-in marker `BlanketClifford` in
`ppvm-traits-2` gating the blanket; the phaseless words + `Tableau` opt in
(keeping the single audited blanket), while `Phased<W>` provides a hand-written
**fused** `impl Clifford` (read each bit once, compute sign, apply bit op, fold
phase) — recovering parity-and-better. Signs verified byte-identical to the old
kernel and to the Lean oracles.

Lean added (prove agent, `lake` green): `Conjugation.lean` `conjSdag` — the
**backward `S†PS`** direction the phase-tracking simulator actually runs
(`conjSdag_sign`, sign `x∧¬z`, `S†XS=−Y`), proven a group hom (`conjSdagHom`) with
generator tables and `conjS_conjSdag` (inverse of the forward `conjS`). S is the
sole convention-sensitive generator; this pins the exact ℤ₄ delta the code emits.

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| ppw2.sconv.1 | correctness/validation | med | proof+design | closed | Simulator runs backward `S†PS` (sign `x∧¬z`, matches old); Lean only had forward `conjS` (`x∧z`). Added `conjSdag_sign` (machine-checks the delta the code applies) + design note; the impl's inaccurate "Lean notes the S/S† convention" claim corrected. No code change (behavior was already correct/at-parity). |
| ppw2.perf.1 | perf-drift | med | human→fix | closed | Phased cnot 1.84× (blanket bits/phase split → redundant bit reads). Maintainer chose fix: `BlanketClifford` opt-in marker + hand-written fused `Clifford for Phased<W>` → 0.85×. Real design change recorded in design docs. |
| ppw2.cite.1 | correctness | low | impl | closed | Crate cited private `conj*_phase_delta` cocycle lemmas + inaccurate S-note; fixed to public `conjCNOT_sign`/`conjCZ_sign`/`conjSdag_sign`. |
| ppw2.alias.1 | impl-friction | low | design | noted | Ships concrete `PhasedPauliWord = Phased<PauliWord>` vs the design's generic renaming; generic `Phased` kept public. |
| ppw2.lossy.1 | impl-friction | low | accepted | deferred | `Phased<LossyPauliWord>` deferred (needs a lossy dep; lossy has no `KeyProduct` ⇒ no phased product). Conjugation path is generic and would work once a dep is added. |

**Architecture change:** the `Clifford` blanket is now **opt-in** via
`BlanketClifford` (`ppvm-traits-2`); `PauliWord`/`LossyPauliWord`(/future
`Tableau`) implement the marker; `Phased<W>` opts out with a fused impl. Recorded
in `traits-2-configuration-and-hashing.md`, the plan, and `word-data-structures.md`.

**Workflow bug fixed (found in verification):** this run first converged *despite*
the 1.84× because the perf-drift gate regex matched the *negated* phrase "NOT the
design-accepted trade-off." Perf-drift now **always** blocks (escalates to human);
the orchestrator caught the regression manually regardless.

---

## Phase 2 — **complete**: `PauliWord` ✓ · `LossyPauliWord` ✓ (own crate) · `Phased` ✓ (own crate). All at parity-or-better. Next: Phase 3 (`ppvm-pauli-sum-2`).

## Phase 3 — ppvm-pauli-sum-2

### Component 1: `Sum` core + graded traits + Clifford — status **complete** (workflow `wf_6325abce-2ed`, escalated iter 3 on a real perf gap; resolved by orchestrator + a focused fix)

Deliverables: graded traits `Support`/`Accumulate`/`Scale`/`Pair`/`Retain` impl'd
on `Vec<(K,C)>` and `HashMap<K,C,IdentityBuildHasher>` (in `ppvm-traits-2`'s new
`containers.rs` — orphan rule: both trait and container are foreign to the sum
crate); `Sum<S,P>` + `apply<TermProducer>`; `RekeyProducer`; policies
(`NoPolicy`/`MaxPauliWeight`/`CoefficientThreshold`/`CombinedPolicy` + `Retain`);
`Clifford for Sum`; `PauliSum<C=f64,P=NoPolicy>` alias. Conformance:
`pauli_sum_diff.rs` (differential vs old, incl. per-gate Clifford replay + X/Y/Z),
`pauli_sum_hash.rs`, `pauli_sum_lean.rs`, bench. Gates verified by orchestrator
(build/clippy/fmt, **186 tests**, machete, `lake`) — all green.

**Perf: faster than old across the board** (the Phase-2 Copy-drop watch item did
**not** regress — move-based re-key + the hash opt below): `clifford_h` **0.86×**,
`clifford_cnot` **0.80×**, `build_batch` **0.29×** (3.5× faster), `scale` 0.99×,
`overlap` ~**700× faster** (old was superlinear), `clifford_x` **0.77×** (after
the fix).

**Cross-cutting optimization (accepted):** `PauliWord`'s hash cache
`OnceLock<u64>` → sentinel `AtomicU64` (relaxed load/store). The `Once` CAS init
path dominated the re-key hot loop (every freshly built key hits cold-init once);
a relaxed atomic is correct here — the digest is a pure function of immutable
content, so a racing miss recomputes *the same* value. Design doc's "Lazy hashing"
section updated. This is *why* the Clifford paths beat the old crate.

Lean added (prove agent, `lake` green; all orchestrator-verified non-vacuous):
`Phase.lean` `IsRealPhase` (real ±1 = even ℤ₄, closed under +); `Symplectic.lean`
`*Act_involutive`/`*_bijective` (the Clifford re-key is a bijection ⇒ no term
collisions — the engine's load-bearing no-collision invariant); `Conjugation.lean`
`conj*_isRealPhase` (Clifford never emits `i` ⇒ the `±1` drain is total, the
`PosI`/`NegI` branch unreachable); `GradedMap.lean` `overlap_eq_fintype_sum` +
`clifford_conjugation_preserves_overlap` (Heisenberg re-key preserves the
Hilbert–Schmidt pairing, via bijection + `s²=1`).

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| ps2.clifford.1 | correctness | med | design | closed | The design's generic `Clifford for Sum` snippet dispatches to the key's own `Clifford`, which for `PauliWord` is bit-only and **drops the ±1 sign** (unsound for `PauliSum`). The impl correctly deviates: wrap each key in `Phased<PauliWord>`, conjugate via the fused sign-tracking Clifford, drain the ±1 to the coefficient. Design snippet corrected. |
| ps2.perf.1 | perf-drift | med | impl→fix | closed | X/Y/Z (pure-sign, word unchanged) rebuilt the whole HashMap instead of an in-place sign flip. **This is what escalated the run** — the hardened perf gate correctly blocked; the loop didn't finish the fix in 3 iters, so the orchestrator drove a focused fix: new `SignFlipByKey` store capability (whole-map, columnar-friendly) → `clifford_x` 0.77×. |
| ps2.containers.1 | impl-friction | low | design | noted | Graded container impls live in `ppvm-traits-2` (orphan rule), so it is no longer "definitions only". Plan crate-map updated. |
| ps2.cite.1 | correctness | low | impl | closed | `producer.rs` mis-cited `xorRelabel_bijective` (the non-Clifford amplitude relabel) for the Clifford re-key; fixed to `Symplectic.*_bijective` / `Conjugation.conj*_injective`. |
| ps2.capacity.1 | impl-friction | low | design | noted | `CoefficientThreshold::capacity` returns `n*10` (ported from old) vs the design sketch's `0`; the ported value is the perf-sensible one — design sketch is a stale example, left as-is. |
| ps2.test.1 | missing-test | low | test | closed | Sum-level X/Y sign drain was untested; added differential X/Y/Z tests (pass), which also guard the new in-place fast path. |

**Deferred to later Phase-3 components** (per the component spec): rotations/noise
producers (`RotationOne`/`PauliError`), the L4 `Multiply` operator product, and
the columnar `ColumnStore` (SoA) backend (Phase 6).

### Component 2: rotations + noise — status **complete** (workflow `wf_1f05f18f-2e7` **failed on an infra error**; work was complete on disk, orchestrator verified + fixed)

Deliverables: `RotationOne` (rx/ry/rz) via a **fused single-pass** branching
producer (`RotateInPlace`: scale each diagonal in place, hash/merge only the ≤N
branch terms — restoring the old crate's `map_insert`, not the batch round-trip);
`PauliError` (pauli_error) via an in-place diagonal `ScaleByKey`. Conformance:
`pauli_sum_rotation_noise_diff.rs` (rotation replay + branch-collision merge +
noise vs old) + `pauli_sum_rotation_noise_lean.rs`; crate `rotation_noise.rs`.
Gates verified by orchestrator (build/clippy/fmt, tests, machete, `lake`) — green.

**Perf — FINAL (this supersedes several earlier wrong versions in the git history
of this entry; the mis-steps are kept as a cautionary record below).** The
benchmark was **made fair**: it now compares the two engines on **identical
storage** (`[u8; 8]` both sides — new via a bench-local `BenchSum`/`BenchKey`, old
via `ByteF64<8>`). Storage matching is confined to the *bench*; the harness's
`NewSum` and the whole differential/Lean suite stay on the **shipped `u64`
default**, so correctness is validated on what actually ships (storage-independent)
while the perf gate is a clean engine-to-engine ratio.

Fair same-storage ratios (new/old, two runs, tight CIs):

| target | new/old | verdict |
| --- | --- | --- |
| `build_batch` (`from_terms`) | **0.30×** | large win (batch vs `+=`) |
| `overlap` | **~0.001×** | large win (new O(n) walk vs old O(n·m)) |
| `clifford_x` (pure-sign) | **0.70×** | win (`SignFlipByKey` in-place) |
| `scale` | 0.99× | parity |
| `pauli_error` | 1.05× | small residual |
| `clifford_cnot` | 1.05× | small residual |
| `clifford_h` | 1.06× | small residual |
| `rotation_rx` | **~1.15×** | small residual, **at gate boundary** |

**Root-caused and FIXED (rotation_rx now at parity).** A focused investigation
subagent isolated the cause by single-variable A/B + disassembly (not a guess):
the residual was **not** the atomic (arm64 `Relaxed` load/store emit no barrier —
verified in the disasm), **not** the hasher (new `IdentityBuildHasher` and old
`FxBuildHasher` do an *equal* number of FxHash rounds), **not** resize (both maps
grow 1000→1472 without re-bucketing), alloc, or monomorphization (all ruled out
with evidence). The real cause: `with_bits_toggled` builds each branch key with an
**empty** hash cache, so the 3-round structural finalize fired **lazily inside
pass-2's `entry()`** — *on the bucket-index critical path*, where the mul-chain
latency stalls the dependent bucket load. The old crate hashes **eagerly in its
first pass** (`rehash` inside `map_insert`), so its pass-2 probe hits a cached `u64`
and the hashing overlaps with other terms' work.

**Fix** (`store.rs`, `RotateInPlace for HashMap`, one line in pass 1): warm the
fresh branch key's digest — `let _ = term.0.key_hash();` — *before* buffering it,
so pass-2 probes a cache hit and the hashing overlaps with the in-place walk
instead of stalling the probe. It computes the **same** digest → a pure semantic
no-op (all 108 tests green, every rotation differential test passes), touches only
the rotate path (`overlap`/`build_batch`/`clifford_x` are different code paths,
unaffected), and leaves the shared avalanche hash untouched. Placed as a single
site inside the *existing* pass-1 loop (no second traversal, cache-hot). Measured:
`new/rx` a stable **~5.5µs → ~4.9µs (~10%)** across 4 runs; the subagent's
interleaved same-build harness (the sound instrument — plain criterion cross-build
numbers are unreliable here, an untouched `old/rx` swings 4.5↔5.9µs+ from code-
alignment/Mytkowicz relayout) put the ratio at **1.07 → ~0.99**, i.e. parity.

The Clifford re-key (`clifford_h`/`cnot`, ~1.05–1.07×) shares the same
fresh-key-hash mechanism but on the `apply`/rekey path; it is small, under the gate,
and left as-is (the same warm technique would apply if it ever matters).

**Cautionary record (this is *why* the workflow was hardened — see the workflow
section).** I mis-called this perf result **three times** before landing it:
(1) called a real 2.42× batch-path measurement "noise" (it was the `apply`
round-trip; fixed for real by the fused `RotateInPlace`); (2) blamed an "inherent
lazy-hash cost" then a "codegen/monomorphization" cause without evidence; (3) —
most insidiously — declared **"same-storage parity ~1.01×"** from a *single* run
whose old baseline happened to land high (5.41µs). The repeated fair measurement
puts old/rx at ~4.7–4.8µs (tight CI) → **~1.15×, not parity**. Cherry-picking a
favorable noisy baseline is the same apples-to-oranges error one level down. The
genuine wins found along the way were kept: the fused `RotateInPlace` (killed the
2.42×) and `with_bits_toggled` (avoids a wasted atomic load/store per branch key).

Lean added (prove agent, `lake` green, orchestrator-verified non-vacuous):
`Rotation.lean` `branchExp`/`branchExp_isRealPhase` + `rx/ry/rz_eps_from_product`
(the branch ±1 sign is **real** — the `i` of `iGP` cancels the product's `i` when
`{G,P}=0` — with the exact per-axis ε formula pinned); `Conjugation.lean`
`conjX`/`conjY`/`conjZ` (the pure-sign X/Y/Z conjugation, witnessing the in-place
fast-path signs).

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| ps2.noise.1 | correctness | med | impl (orch) | closed | `pauli_error` with a zero transfer eigenvalue (reachable, `[0,0.25,0.25]`→λ_X=0) left a **phantom zero-coefficient key** in the support, violating the reduced-canonical-form invariant. Orchestrator fixed `ScaleByKey` (both backends) to drop any term scaled to exactly zero; added a `pauli_error_zero_eigenvalue_drops_the_term` test. |
| ps2.rot.perf | perf-drift | med | impl→fix | **closed — fixed (parity)** | Fair same-storage bench first exposed a real ~1.15× on `rotation_rx` (correcting an earlier cherry-picked "parity ~1.01×" from a noisy old baseline). Root-caused by subagent (single-var A/B + disasm): the fresh branch key's 3-round hash fired lazily inside pass-2's `entry()` on the bucket-probe critical path, vs the old crate's eager first-pass hash. **Fix:** warm `term.0.key_hash()` in `RotateInPlace`'s pass 1 (same digest → no-op; 108 tests green). `new/rx` ~5.5µs→~4.9µs (~10%, stable ×4); interleaved-harness ratio 1.07→~0.99. Ruled out (with evidence): atomic (no barrier), hasher (equal rounds), resize, alloc, monomorphization. |
| ps2.nsites.1 | impl-friction | low | impl | deferred | The `n_sites` `debug_assert!` is on `from_terms` but not on the `apply`/`rekey_bijective` produced-key paths. No live bug (every producer preserves width); low, deferred. |
| ps2.xyz-dup.1 | impl-friction | low | accepted | noted | X/Y/Z sign logic (`(−1)^z`, etc.) is written a third time in the in-place fast path (besides `PhaseTrack` and `Phased`). Now has a Lean witness (`conjX/Y/Z`) and differential tests guard it; accepted — the in-place path exists precisely to avoid the `Phased` wrapper's cost. |

**Infra note:** the workflow *failed* (a reporting agent hit the StructuredOutput
retry cap — a transient JSON-schema failure), but 4/5 agents completed and the
implementation + tests were fully on disk and green. No resume needed; the
orchestrator verified, fixed the noise bug, and finished.

**Workflow hardening (this run) — perf gate now requires fair + stable +
attributed.** The `ps2.rot.perf` saga (three wrong perf calls: "noise", guessed
"codegen", then a cherry-picked "parity") was a *measurement-discipline* failure,
so the TEST agent's perf gate (`.claude/workflows/traits-2-component.js`) now
enforces three rules, and the schema makes each a required field per benchmarked
target: **(a) fair config** — identical storage width / coeff / hasher on both
sides (`config` field; build a bench-local storage-matched new type if the shipped
default differs — as `pauli_sum_bench.rs`'s `BenchSum` now does — since correctness
is storage-independent); **(b) stable measurement** — ratio confirmed over ≥2 runs
with tight CIs (`stable` field; a wide-CI baseline may not be read as drift in
*either* direction); **(c) verified attribution** — a flagged drift's cause must
come from a controlled A/B or profile, else `"unattributed"` (no guessing). The
orchestrator's convergence-loop guard now rejects a perf-drift gap that is unfair,
noisy, or unattributed (bounce back to re-measure) instead of escalating a
benchmark artifact to the human.

**Workflow hardening #2 — integration baseline + architecture-parity review.** The
aux-map miss exposed a deeper process gap: the `-2` conformance suite had only
*single-op microbenches*, which structurally cannot see cumulative per-gate costs
(a tight one-gate `iter` loop lets the allocator recycle one warm page), and the
review agent never checked whether the new impl preserved the old's perf-critical
*architecture*. The old crate's real workloads (`benches/trotter*.rs`,
`random-circuit.rs`, `truncation-weight.rs`, `tests/trotter.rs`, `ghz.rs`) were
never ported. Fixes to `.claude/workflows/traits-2-component.js`:
- **New `Baseline` phase (runs FIRST, behavioral components only):** an agent mines
  the old crate's real workloads to define the **integration acceptance bar**
  (end-to-end new-vs-old *numeric* golden masters + *perf* ratios) and to enumerate
  the **perf-critical architecture features** the new impl must preserve (buffer/aux
  reuse, packed layout, allocation strategy, in-place fast paths). Its brief is
  threaded into every downstream agent.
- **Impl** must preserve each enumerated feature (may relocate it — e.g. the aux
  into the storage backend — but not silently drop it; omissions justified in
  `frictions`).
- **Review** gains an **architecture-parity** dimension: for each baseline feature,
  verify the new impl still carries it; a dropped one is a HIGH `impl-friction` gap
  routed to impl (the check that was missing when the aux-map vanished).
- **Test** perf gate is now **integration-first**: the headline ratio is the
  end-to-end deep-circuit propagation (microbenches are diagnostic only), plus an
  end-to-end numeric golden-master as part of `differentialPass`.
- **Lean backstop:** integration numerics target old-parity by default, but the
  Baseline agent flags suspected old-impl bugs and the Prove agent adjudicates them
  against the Lean oracle — if old is wrong, the golden-master targets the
  Lean-correct value and the new impl is right by construction.
(Concrete exemplar back-filled for `ppvm-pauli-sum-2`: an end-to-end Trotter-TFIM
propagation diff-test + bench, new vs old — the coverage that would have caught the
aux gap. See the ps2.integration.1 row.)

**Storage: auxiliary double-buffer recovered (in the storage, not on `Sum`).**
Maintainer flagged that the old crate's `PauliSum` held a persistent **double-
buffer** (`map: (primary, aux)` swapped via an `aux` flag) + a reusable `scratch`
Vec, so Clifford re-key and rotation never allocated per gate — and that the `-2`
refactor had *dropped* it: `HashMapStore` was a bare `type = HashMap<…>`, so
`rekey_bijective` allocated a fresh map + `mem::replace` every gate and
`rotate_in_place` a fresh branch `Vec` every gate. The generalization ("`Sum` owns
only its storage" so it composes into the generalized tableau) was intentional and
stays; losing the aux was not. **Fix:** `HashMapStore` is now a **newtype backend**
owning `{ primary, aux, scratch }` — the double-buffer relocated *into the storage*,
so `Sum` stays a pure generic engine. `RekeyBijective` = the old `map_add`
(clear aux → drain-move primary through the re-key into aux → swap; the two
allocations ping-pong, never freed); `RotateInPlace` uses the persistent `scratch`.
Graded traits + `ScaleByKey`/`SignFlipByKey` delegate to `primary` (zero-cost,
inlined); the three superseded raw-`HashMap` capability impls were removed. 108
tests green; delegated bench paths unchanged. Same-build ratios **improved**:
`clifford_h` ~1.08→**1.04×**, `clifford_cnot` ~1.05→**1.04×** (the per-gate
fresh-map alloc is gone), `rotation_rx` **0.94×** (scratch reuse + the pass-1 warm).
Note this partly supersedes the subagent's earlier "leave Clifford as-is": its
*drain-reuse-of-self* experiment regressed, but the **double-buffer swap** (a
different mechanism — no extra traversal) it hadn't tried does help.

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| ps2.store.aux | impl-friction | med | design→fix | **closed — recovered** | Old persistent `(primary, aux)` double-buffer + `scratch` was dropped when `HashMapStore` became a bare `HashMap` alias (per-gate alloc). Recovered by making `HashMapStore` a newtype owning `{primary, aux, scratch}` — aux lives *in the storage* so `Sum` stays generic for the tableau. `clifford_h/cnot` ~1.05–1.08×→~1.04×, `rotation_rx`→0.94×; 108 tests green. |

**Integration coverage added — and it immediately found a real regression the
microbenches hid.** Back-filled the end-to-end coverage the `-2` suite was missing
(porting the old crate's Trotter workload): `tests/pauli_sum_integration_diff.rs`
(3 tests — Trotter-TFIM end-to-end new-vs-old golden master, `rzz=cnot·rz·cnot`
decomposition validation, deep 400-gate random-circuit replay; identical config
both sides: `[u8;8]`, `CombinedStrategy/Policy(CoefficientThreshold(1e-6),
MaxPauliWeight)`) and `benches/pauli_sum_integration.rs` (the old `trotter` bench:
n=12, 10 steps, same storage). **Discriminating power validated:** removing the aux
moves the integration ratio ~6× more than the microbench (+13% vs +2.8%) — proof
the deep circuit sees per-gate allocation the tight one-gate loop cannot.
**But the headline result is a standing gap:** even *with* the aux, the new engine
is **~1.4–1.5× slower than old on the real Trotter workload** (new ~1.14ms vs old
~0.76ms, same-build) — while every single-gate microbench reads ~0.94–1.05×. The
microbench suite was reporting "parity" on a workload the engine is ~1.45× behind
on. Prime suspect: `truncate()` — called after *every* gate in the Trotter loop and
run twice per call by `CombinedPolicy` (two full retain scans), yet it has **no
microbench** at all. To be root-caused next (same discipline: same-build A/B).

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| ps2.integration.1 | missing-test | high | test (orch) | **closed — added** | The `-2` conformance suite had only single-op microbenches; the old crate's real-workload benches/tests (Trotter, random-circuit) were never ported, so nothing exercised cumulative allocation/truncation. Added end-to-end diff-test + bench (see above). This is the coverage that would have caught the aux-map miss (and did, in retro). |
| ps2.trotter.perf | perf-drift | high | — | **root-caused → split** | New ~1.45× slower than old on the end-to-end Trotter workload (new ~1.14ms / old ~0.76ms, same-build `[u8;8]`), invisible to single-gate microbenches. Attribution done (`examples/trotter_attrib.rs` + A/B): two independent causes → split into `ps2.truncate.behaviour` and `ps2.cnot.rekey.perf` below. Numerics match. |
| ps2.truncate.behaviour | correctness | high | impl | **closed — fixed (`d4aa7875`)** | New `Sum` gate drivers auto-truncated internally (`apply`/`rekey_bijective`/`rotate_in_place` each called `self.policy.truncate`); **old gates never auto-truncate** (caller-driven, pinned by old's `tests/truncation_semantics.rs`). Under the PRIME DIRECTIVE that changed user-facing behaviour. Fixed: truncation is explicit-only via `Sum::truncate()`; `reduce()` kept (canonicalisation, not policy). Guard added: `tests/pauli_sum_truncation_behaviour_diff.rs`, 2 of its 3 tests verified to fail against the pre-fix engine. |
| ps2.cnot.rekey.perf | perf-drift | high | impl | **closed — root-caused + fixed (`23860055`)** | `entry(k).and_modify(..).or_insert(..)` in `RekeyBijective` aggregated a collision that cannot occur (`f` injective by contract; `Symplectic.*_bijective`) and compiled **out of line** (41% of sweep samples in a standalone `hashbrown::rustc_entry` frame vs fully inlined for old). Replaced by `insert` + `debug_assert!`. `rekey_cnot` **1.45× → 1.11×**, end-to-end Trotter **~1.47× → ~1.28×**, over 8 interleaved process pairs. Residual 1.11× unattributed. Earlier steps: `b68b7e8e` (content-conditional hash invalidation) took it from 1.71×. Ruled out with evidence along the way: `aux.reserve`, `drain` vs `iter`+clone (~2.4–3%), hash-warm/extra-pass, `truncate` (parity 0.99×). |
| ps2.rekey.bimodal | perf-drift | high | impl | **OPEN — narrowed, not closed** | Re-key is **bimodal across processes** for the identical binary and input. `insert` improved both modes (fast 77→59 µs, slow 134→128 µs) but did not remove it: 3/6 `cargo bench` launches still land slow. Startup-fixed state, not drift (two in-process measurements agree to 0.5%); old is stable ~53 µs always. **New:** 16 direct-binary launches were all fast while `cargo bench` launches of the same binary are bimodal ⇒ the trigger is the launch environment; mechanism **unattributed**. Ruled out with evidence: heap layout (tables identical mod 16 KB across 22 processes), env/cwd, CPU core type, code-ASLR slide, differing hot path. **Refuted:** `AtomicU64`→`Cell<u64>` changed nothing — the `pw2` Copy/atomic watch item is not a cause. |
| ps2.zero.behaviour | correctness | high | design | **partly fixed (`e3e1f9cc`) — `reduce` still OPEN** | Old has **no `reduce`** and no drop-zero logic anywhere; a zero-coefficient term stays in the support and old's exact-map `PartialEq` depends on it (`ppvm-pauli-sum/tests/loss.rs::test_reset_channel`). New `reduce()` and `ScaleByKey` both drop exact zeros ⇒ `len()`/`contains` diverge (e.g. `pauli_error(q,[0.0,0.25,0.25])`, λ_X=0). Verified by hand against the old sources. **Indicts `ps2.noise.1`**: that row's "fix" *added* zero-dropping to satisfy a reduced-canonical-form invariant which is the new design's invention, not old's behaviour — so it introduced this divergence. **`ScaleByKey` fixed** (`e3e1f9cc`): now old's `Fn(&K,&mut C)` + `iter_mut`, which cannot remove — and that took `pauli_error` 2.03×→1.06×. Guarded by `zero_channel_eigenvalue_keeps_the_term_exactly_as_old_does` (exact-support diff; verified to fail pre-fix). **Still open:** `reduce()` drops zeros on `apply`/`from_terms`. |
| ps2.default.threshold | correctness | med | impl | **OPEN — new** | `CoefficientThreshold::default()` is `1e-12` in old (`ppvm-pauli-sum/src/strategy.rs:113`) but `0.0` in new (derived `Default` on the `f64` field, `policy.rs`) — a silently changed user-facing default. Verified by hand. |
| ps2.preserve.missing | impl-friction | med | design | **OPEN — new** | `preserve_strings` (builder option + snapshot/restore post-filter inside old `PauliSum::truncate`, pinned by `ppvm-pauli-sum/tests/preserve.rs`) has no counterpart in the new crate. |
| ps2.oldbug.mulassign | correctness | med | proof | **OPEN — suspected OLD bug** | Old `impl MulAssign<PauliSum<T>> for PauliSum<T>` (`sum/ops.rs:70`) calls `self.map_add(..)` per rhs term, and `map_add` *replaces* the support with its image ⇒ `A *= (b0·P0 + b1·P1)` computes the product chain `A·b0P0·b1P1`, not `A·b0P0 + A·b1P1`. Untested in old. For the L4 `Multiply` component: accumulate into a fresh accumulator per `twistedConv` (`Twisted.lean`); if Lean confirms old is wrong, this is an allowed documented divergence. |
| ps2.attrib.instrument | missing-proof | med | test | **OPEN — instrument unsound** | `examples/trotter_attrib.rs` understates the gap: its per-op `Instant::now()` pairs act as optimization barriers that penalize the *old* engine more, reporting `1.08×` where criterion reports `1.35×` on the same tree. Its per-op ratios are diluted toward 1.0. Earlier conclusions in this log that leaned on it should be re-derived from `benches/pauli_sum_integration.rs`. Either fix (sample-based attribution) or demote to ordering-only. |

**Deferred to component 3:** the L4 `Multiply` operator product (needs `Complex`
+ `ImaginaryUnit`); columnar `ColumnStore` stays Phase 6.

### Component 3 — L4 `Multiply` (operator product) — pending.

## Phase 4 — ppvm-tableau-2

_(no entries yet)_

## Phase 5 — ppvm-sym-2

_(no entries yet)_

## Phase 6 — ColumnStore backend

Status: **complete** (human perf/design sign-off 2026-08-05).

Deliverables:

- `ColumnStore<K, C>` is a real structure-of-arrays backend: `K::Column` key
  planes, contiguous coefficients and digests, and a grouped-probe row index.
  The implementation is split under `src/column_store/mod.rs`, including
  dedicated rotation and lifecycle-test modules. `ColumnPauliSum` is the
  drop-in alias.
- Persistent aux/scratch/batch workspace, exact `7/8` index sizing, in-place
  bijective re-key, two-pass rotations, size-reserved branch append, and
  capacity-preserving clone are all present. Repeated L4 products now reserve
  relative to the accumulator's existing length; the aux-reuse test exercises
  the same store twice instead of cloning away the workspace.
- The old public collection/operator spellings are restored:
  `IntoIterator`, replacement-semantics `Extend`, by-value sum product, and
  single-word `*`/`*=`. The documented Lean-correct product remains bilinear.
- Three-way old/hash/column differential coverage includes the Trotter and random
  workloads, truncation/preserve/zero timing, pairings, L4, the complete gate
  surface, and lossy words/channels. A separate old/new lossy-sum differential
  suite closes baseline workload 6.
- `Projection` follows the Lean oracle rather than old's two independent bugs:
  coefficients use the ring constant `½` (`c/2`, not `c²/2`) and `X`/`Y` are
  annihilated. Their zero-valued keys remain until caller-driven `reduce()`.

The lossy integration benchmark exposed a cross-component regression hidden by
the Phase-2 word microbench: three `OnceLock` component caches made the real
loss-interleaved sum **2.45×** slower than old. Replacing them with the same
relaxed-atomic sentinel mechanism already validated on `PauliWord` preserves the
component split and loss-only invalidation while removing the `Once` path.
Controlled stage ratios improved from reset `4.92×`→`1.24×`, loss
`2.55×`→`1.36×`, correlated `2.21×`→`0.72×`, Clifford `1.74×`→`1.26×`, and
rotation `2.64×`→`1.31×`. The realistic n=12 loss-interleaved GHZ workload is
inside the gate in two processes: new/old medians **0.99×** and **0.84×**.

Final performance gate, four fresh process runs after the index, rotation,
noise-codegen and tombstone fixes:

| target | ratio | verdict |
| --- | --- | --- |
| Column steady `rx` / hash | **0.94–0.97×** | fixed; paired dense kernel |
| Column sustained `rx` growth / hash | **0.93–0.96×** | fixed; grouped row index |
| Column `from_terms` / hash | **0.89–0.92×** | fixed; lazy workspace + grouped index |
| Column active truncate / hash | **0.94–0.95×** | fixed; lazy tombstones, stable 1/8 compaction |
| Column realistic `rx` / hash | **0.95–0.98×** | fixed; cached sparse-row runs |
| Column native `rzz` / hash | **0.93–0.95×** | fixed; column-native packed-bit kernel |
| Column / hash noisy TFIM Trotter | **0.98–0.99×** | no regression |
| New hash / old noisy TFIM Trotter | **0.86–0.90×** | faster |

The full earlier op sweep showed the intended coefficient/re-key wins
(`scale` ~0.13×, `reduce` ~0.40×, `h` ~0.37×, `cnot` ~0.53×,
threshold/weight/channel walks ~0.64–0.76×). The final grouped index, packed
column rotations, and tombstoned retention remove the former losing rows
without weakening those scan/re-key wins.

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| col6.clone.capacity | perf-drift | high | impl | **closed** | Clone now preserves primary key/coeff/digest capacity plus aux/scratch/batch capacity. |
| col6.multiply.reserve | perf-drift | med | impl | **closed** | `multiply_into` reserves `existing + product_hint`; repeated accumulation no longer walks a resize chain. |
| col6.rot2.copy | perf-drift | med | impl | **closed** | All two-site rotations use one `toggled_bits2` plane copy. |
| col6.index.load | correctness | low | impl | **closed** | Exact ceil-at-7/8 sizing, grow only above 7/8, plus a debug duplicate-key guard during reindex. |
| col6.loss.diff | missing-test | high | test | **closed** | Exact old/new loss-channel, correlated-loss, GHZ, zero-survival and max-loss-weight differential suite + integration bench. |
| col6.surface | missing-test | med | test | **closed** | Hash/column parity now covers extended/batched Clifford, generic/named RotationTwo, RotXY, all noise/projection paths, Hermitian overlap and lossy ColumnStore. |
| col6.projection | correctness | high | impl/proof | **closed** | Lean adjudicated both old defects; Rust and discriminating tests now target the matrix-correct projector. |
| col6.batch.scalar | impl-friction | med | design/human | **closed — explicitly deferred (human, 2026-08-05)** | `TermBatch` keys and `probe_batch` remain scalar; a column-native produced batch/coalesced gather requires changing the shared batch contract. Accepted as follow-up because the shipped SoA support is complete and end-to-end faster. |
| col6.insert.perf | perf-drift | med | impl | **closed — fixed** | Replaced scalar probing with a grouped row index, fused paired packed-column rotations, made aux allocation lazy, and added separate dense/sparse kernels. Four-process ratios: steady `rx` 0.94–0.97×, sustained growth 0.93–0.96×, `from_terms` 0.89–0.92× versus hash. |
| col6.truncate.perf | perf-drift | med | impl | **closed — fixed** | Active 408→390 truncation exposed physical SoA compaction (2.4–2.5× hash). Retain now tombstones rows immediately, rebuilds sparse traversal lazily, and compacts stably at 1/8 dead rows. Final column/hash ratio 0.94–0.95×; Trotter remains 0.98–0.99×. |
| col6.aligned.proof | missing-proof | low | proof | **OPEN** | Alignment is mutation-tested after every operation, but Lean's abstract `Finsupp` model cannot express misaligned physical columns. |

## Phase 7 — downstream + cutover

Status: **adapter phase complete; the core performance blocker is closed.
Destructive rename/removal awaits maintainer approval and review.**

First adapter prerequisite landed: `ppvm-tableau-2::GeneralizedTableau::trace`
now consumes the `-2` `PauliPattern`, enumerates the bounded accepted words, and
delegates each leaf to `expectation`. Differential coverage against old's
counted-pattern trace is green. This closes Phase 4's trace deferral and removes
one blocker for the vihaco/Python adapters. The full old stateful pattern grammar
and greedy matcher are now ported as a dedicated `pattern/` module; bounded
enumeration matches old, including rejecting every star pattern.

Phase 7 Step 1 is complete for `ppvm-stim`. Its default `legacy` and opt-in
`traits-2` features are mutually exclusive, both implement the same sealed
semantic executor adapter, and `rayon` remains orthogonal (shot-level
parallelism only). The former 872-line executor is split into a `mod.rs`
hierarchy (largest file: 159 lines). Exact verification:

- `cargo test -p ppvm-stim`: 106 passed, 0 failed; 1 doctest ignored.
- `cargo test -p ppvm-stim --features rayon`: 106 passed, 0 failed; 1 doctest ignored.
- `cargo test -p ppvm-stim --no-default-features --features traits-2`: 106 passed, 0 failed; 1 doctest ignored.
- `cargo test -p ppvm-stim --no-default-features --features traits-2,rayon`: 106 passed, 0 failed; 1 doctest ignored.
- `cargo check -p ppvm-stim --all-targets` passed in all four matrices above.
- `cargo clippy -p ppvm-stim --all-targets -- -D warnings` passed in all four matrices.
- `cargo fmt --all -- --check` passed.
- `cargo check --workspace --all-targets` passed with the default legacy backend.
- No-backend and dual-backend checks reached their intended `compile_error!` guards.
- Review follow-up closed three gaps: width aliases now preserve their real
  `[u8; N]` / `[usize; N]` capacity, serial `sample` accepts non-`Sync`
  factories (guarded by an `Rc<Cell<_>>` regression test), and CI runs explicit
  traits-2 native/rayon plus wasm-library checks. Public default storage is
  target-correct (`u64` native, `usize` wasm); the selected traits-2 Stim library
  builds on `wasm32-unknown-unknown`. Criterion benches remain native-only.

Phase 7 Step 2 is complete for `ppvm-vihaco`. The default `legacy` and opt-in
`traits-2` features select aliased old/new dependencies and are mutually
exclusive. Public `Circuit`/`PPVM` routing is unchanged; constructors, policy
spelling, width-specific words, seeded initial states, and rendering are
normalized behind backend modules. The former 1032-line `component.rs` is a
`mod.rs` hierarchy (largest file: 208 lines). The traits-2 width aliases use the
full `Sum<HashMapStore<PauliWord<[u8; N]>, f64>, CombinedPolicy<...>>` and
equivalent `LossyPauliWord` forms. `rayon` remains orthogonal; traits-2 exposes
shot-level/downstream parallelism only and does not claim tableau-internal
parallelism. Exact verification:

- `cargo test -p ppvm-vihaco`: 115 unit tests + 29 fixture tests passed, 0 failed.
- `cargo test -p ppvm-vihaco --features rayon`: 117 unit tests + 29 fixture tests passed, 0 failed.
- `cargo test -p ppvm-vihaco --no-default-features --features traits-2`: 115 unit tests + 29 fixture tests passed, 0 failed.
- `cargo test -p ppvm-vihaco --no-default-features --features traits-2,rayon`: 117 unit tests + 29 fixture tests passed, 0 failed.
- All 18 existing `.sst` fixtures ran in both backend matrices; the seeded
  branch fixture matches the same `[Zero], [Zero]` snapshot in legacy and traits-2.
- `cargo check -p ppvm-vihaco --all-targets` passed in all four matrices.
- `cargo clippy -p ppvm-vihaco --all-targets -- -D warnings` passed in all four matrices.
- `cargo fmt --all -- --check` passed.
- `cargo check --workspace --all-targets` passed with the default legacy backend.
- No-backend and dual-backend checks reached their intended `compile_error!` guards.

Remaining gaps are outside Steps 1–2: the TUI, CLI, and Python adapters are not
cut over, and `ppvm-tableau-2`'s internal coefficient-level rayon branch remains
deferred. Neither Stim nor Vihaco enables that branch under traits-2.

Phase 7 Steps 3–4 are complete for the transitive frontends. `ppvm-tui` and
`ppvm-cli` expose the same mutually exclusive `legacy` / `traits-2` selection
and forward it to Vihaco; the CLI keeps its historical default Rayon behavior as
an orthogonal feature. Default and traits-2 TUI tests pass, as do CLI all-target
tests under default legacy+rayon and traits-2+rayon.

Vihaco's traits-2 tableau buckets preserve their semantic widths on both pointer
sizes: `[usize; N]` on native and doubled `[usize; 2N]` on wasm32. Its selected
traits-2 library passes the explicit wasm cross-check; CI now exercises Vihaco
and both frontends under traits-2 in addition to the Stim matrices.

Review also restored the public `TableauCircuit::new` /
`TableauCircuit::new_with_seed` overloads (a private shared builder carries the
optional seed) and made the two pure-Rust fixture references feature-selected.
The traits-2 fixture matrix now constructs explicit new `Sum<HashMapStore<...>>`
lossless/lossy states rather than always using the legacy dev oracle. Both
review-fixed legacy/traits-2 suites and the traits-2 wasm library check pass.

The first traits-2 wasm cross-check found that public defaults and temporary
Stim aliases hard-coded `u64`, which is not a `bitvec::BitStore` on 32-bit wasm.
The default is now target-correct (`u64` native, `usize` wasm), while compatibility
aliases preserve their actual `[u8; N]` / `[usize; N]` widths. Native word/tableau
suites and the selected traits-2 Stim wasm library check pass. Criterion benches
remain native-only, so wasm CI checks the library rather than `--all-targets`.

The Phase-7 specialized tableau mixture is complete. `ppvm-tableau-2::mixture`
exports `GeneralizedTableauMixture<A, I, H>`, the compatibility
`GeneralizedTableauSum` alias, and `MixtureSampler`. It retains fingerprint
buckets followed by collision-checked frame/loss and approximately-equal
amplitude comparison; RNG, record, scratch, and probability are excluded from
identity. Strict cutoff and normalization timing, lazy Pauli/loss branch
materialization, analytic case-a/case-b measurement/reset, the complete gate and
noise surfaces, and seeded serial/native-Rayon sampling match the old oracle.
The implementation is split into files below 200 lines.

Differential coverage includes structural snapshots, deterministic gates,
measurement cases, reset coalescing, noise/loss, exact cutoff boundaries, wide
`u128` indices, forced fingerprint collisions, and seeded sampler streams.
`ppvm-tableau-sum`'s 109 legacy tests and the new 7-test differential suite pass.
The same-build benchmark was run in four fresh processes: noisy build
**0.75–0.76×** old and serial 128-shot sampling **0.94–0.95×** old (new faster;
no regression). An initially oversized benchmark workload did not complete a
warmup and produced no ratio; it was reduced before measurements were recorded.

Review found one boundary divergence: construction inserted the initial
probability-1 entry without applying the strict `probability > sum_cutoff` rule.
Construction now goes through the normal insertion door, so cutoffs `0.999`,
`1.0`, and `2.0` produce lengths `1`, `0`, and `0` exactly as old. In-crate and
old/new differential guards pass; changed-crate Clippy remains clean. The other
review item—the old/new generic-axis mismatch—is intentionally owned by the
feature-gated Python backend facade, not by pretending the compatibility alias
accepts old `Config` parameters.

Python adapter preflight found that the public 1–2048-qubit tableau classes need
wide `bnum` amplitude indices. `ppvm-tableau-2` now depends on `bnum` with
`numtraits`, re-exports `U256`/`U512`/`U1024`/`U2048`, and runs each through a
real generalized-tableau gate sequence (including the existing 200-qubit tier).
The wide-index test, all-target Clippy, and dependency hygiene pass.

Python's insertion-order ABI prerequisite is also complete.
`IndexMapStore<K, C>` implements the full `Sum` algebra/gate/lifecycle surface,
including persistent workspace reuse, exact-zero support, replacement `Extend`,
and legacy-compatible re-key/branch ordering. All ordered unit and five old/new
differential tests pass. Two same-build processes show no regression: build
**0.57×**, ordered term export **0.73–0.75×**, CNOT **0.81×**, and rotation
**0.89–0.94×** old.

Review closed three order holes before binding it to Python: multi-branch merge
direction now compares the **deduplicated** branch map cardinality (not raw
fan-out), multiple preserved keys are snapshotted by support order rather than
the preserve `HashSet` order, and clone/in-place-product tests compare ordered
term vectors instead of `IndexMap`'s order-insensitive equality. Discriminating
unit/differential tests and changed-crate Clippy pass.

The Phase-7 Python native cutover is complete behind mutually exclusive
`legacy` (default) and `traits-2` features. The facade selects aliased optional
Pauli/tableau/mixture dependencies and exactly one Stim backend. Traits-2
ordinary and lossy sums use fixed `[u8; 2^N]` words, `IndexMapStore`, explicit
`CombinedPolicy` values, Python's `usize::MAX` loss default, and caller-driven
truncation. Tableau rows remain `[usize; N]`; the amplitude index tiers remain
`usize`, `u128`, and `U256` through `U2048`. The specialized mixture and sampler
replace the old `ppvm-tableau-sum` path, and native `GeneralizedTableauSum.r` is
now present in both modes.

Fresh mixed-project maturin builds through `ppvm-python/pyproject.toml` each pass
the complete 238-test Python/docs suite. The extension exports and the tests
assert `backend_name()` (`legacy` / `traits-2`), so a stale or wrong-feature
`ppvm._core` cannot masquerade as parity. Added coverage fixes the native class
inventory and checks ordered
rendering, 8/9/64/65/128/129/200/1024/2048 width boundaries, genuinely
greater-than-u128 Python coefficient keys, trace/loss encodings, tableau
copy/deepcopy RNG streams on randomized qubits, equal-seed fork streams plus
cross-seed diversity, retained `GeneralizedTableauSum.r` behavior, mixture
sampler snapshots, and seeded Stim output across one and four Rayon threads.
Both native Clippy matrices, traits-2 Pauli/tableau
and selected Stim wasm library checks, formatting, and `cargo machete` pass.
No ABI gap or measured performance regression was observed in this cutover;
the final full-suite wall times were 36.59 s traits-2 and 48.29 s legacy, which
are coarse sequential test-run observations after separate builds, not a
benchmark and not used as a performance claim.

The three shared macro-generated binding implementations are now module
directories. Their state, gate, rotation/noise, Stim, and sampler surfaces use
cohesive private submodules (all at or near 200 lines) while the public module
paths and generated class inventory remain unchanged. Fresh legacy and
traits-2 mixed-project builds each pass the complete 238-test Python/docs
suite, including backend identity and native class inventory checks; both
native check/Clippy matrices, formatting, and `cargo machete` also pass.

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| cutover.trace | impl-friction | med | impl/test | **closed** | `GeneralizedTableau::trace(&PauliPattern)` added on the `-2` tower; counted `Z?{n}` differential test matches old. |
| cutover.pattern.parse | impl-friction | med | impl/test | **closed** | Full old grammar and greedy stateful matching ported: literals, `_`, alternation, optional identity, counted/star repetition, absolute positions and sequential stars. Differential parser/matcher/enumeration tests pass on ordinary and lossy words; bounded star enumeration preserves old's panic contract. |
| cutover.stim | impl-friction | high | impl/test/review | **closed** | `ppvm-stim` has mutually exclusive legacy/traits-2 adapters with source-compatible execute/sample APIs; all four native matrices pass. Review restored width-accurate aliases and serial non-`Sync` factories; explicit wasm library and CI gates pass. |
| cutover.vihaco | impl-friction | high | impl/test/review | **closed** | `ppvm-vihaco` has stable `Circuit`/`PPVM` routing over mutually exclusive legacy/traits-2 facades; review restored public constructors, direct traits-2 fixture references and pointer-width-correct wasm aliases. All fixtures/matrices and wasm library check pass. |
| cutover.frontends | impl-friction | med | impl/test | **closed** | TUI/CLI forward mutually exclusive backend features; both TUI matrices and both CLI all-target matrices pass, with Rayon kept orthogonal. |
| cutover.wasm.default | correctness | high | impl/test | **closed** | Replaced invalid wasm `u64` defaults with native-word defaults and width-accurate compatibility aliases; traits-2 Stim library builds for `wasm32-unknown-unknown`. |
| cutover.rayon | impl-friction | med | impl/design | **OPEN — stim/vihaco resolved** | Stim and Vihaco map traits-2 rayon to shot-level/downstream parallelism only. `ppvm-tableau-2` still records old's internal coefficient-level rayon branch as deferred, and other downstream mappings remain. |
| cutover.mixture | impl-friction | high | design/impl/test/review | **closed** | Added specialized `GeneralizedTableauMixture` + compatibility alias and sampler; fingerprint buckets are collision-checked against full structural/approximate-amplitude identity, avoiding the invalid generic `Sum` alias. Review restored constructor cutoff parity. Differential, seeded, collision, wide-index, wasm, and four-process benchmark gates pass. |
| cutover.python.wide-index | impl-friction | high | impl/test | **closed** | Added/re-exported `bnum` U256/U512/U1024/U2048 `Bitstring` tiers and executable generalized-tableau coverage through 1600 qubits. |
| cutover.python.ordered-sum | impl-friction | high | impl/test/review | **closed** | Added full `IndexMapStore` + fixed-width `IndexPauliSum`; review fixed deduplicated branch cardinality, multi-preserve restore order and false-positive equality tests. Order-sensitive old/new tests pass and all measured paths are 0.57–0.94× old. |
| cutover.python.sum-r | impl-friction | med | impl/test | **closed** | Native `GeneralizedTableauSum.r` is exported in both backend modes and parity-tested through the Python `RotationsMixin`. |
| cutover.python.binding-split | maintenance | low | refactor | **closed** | Split all three shared binding generators into `mod.rs` module directories with cohesive multi-`#[pymethods]` submodules. Both fresh 238-test backend suites, class inventory/backend identity checks, native check/Clippy matrices, formatting, and machete pass. |
| cutover.core-perf | perf-drift | high | impl/human | **closed — performance gate cleared (`73580afa`)** | After `0fc0c57e` and `73580afa` plus final adjudication, the original 75 above-gate rows are 14 fixed, 11 parity, 0 actionable, and 50 non-actionable. The five grouped adjacent rows are 3 fixed, 1 parity, 1 non-actionable, and 0 robust. All 80 observations have no actionable regression; raw ranges and adjudications are in `docs/performance-report.md`. Destructive cutover remains a maintainer approval/review decision. |
| perf.harness | maintenance | med | tooling | **closed** | `mise run perf-report` (`benchmarks/perf_regression_report.py`) drives the whole conformance matrix, pairs the medians, classifies against the 1.03 gate and refuses to call a one-launch row robust. Reproduced 15 of the audit's adjudicated controls within measurement noise on first use. |
| perf.reaudit | perf-drift | med | impl/human | **closed — no actionable regression at `e3a37026`** | Post-RNG-injection re-audit: 828 pairs, 590 improved, 170 parity, 0 actionable. `pauli_error_sweep` (1.158×) is executable placement — two alignment-perturbed builds put NEW at ~4395 ns against its 5010 ns default, and a work-removing ablation measured *worse*. `pauli_sum_surface/add/term` is unchanged at 1.67× since `73580afa`; the `1.014×` headline in `docs/performance-report.md` does not reproduce at its own commit and is corrected there. |

### Full core benchmark audit (2026-08-07)

The comparative harness now covers every comparable public operation across
ordinary/lossy/phased words, Pauli sums and storage variants, bare/generalized
tableau, tableau mixtures, and symbolic propagation. The audit fixed unequal
capacity, timed parsing on one side, decaying channel state, stale native-vs-
decomposed paths, broad Criterion filters, and missing output assertions.

The complete screening contains **901 old/new pairs**. Confirmed regressions
were rerun in fresh processes with 20 samples, 1 s warm-up and 2 s measurement;
one pathological mixture branch expansion was stopped during process three, so
50 under-sampled rows remain explicitly provisional. Headline confirmed gaps:

- tableau observation helpers: 2.56–3.95×, but only +2.8–5.9 ns and
  implementation-identical (code placement/inlining artifact);
- lossy branch-key construction: 1.60–1.80× from three atomic caches plus
  guarded invalidation;
- disabled truncation sentinels: 1.52–1.67×, about +1 ns, unattributed;
- symbolic propagation/trace: roughly 1.19–1.71×, mostly unattributed;
- Pauli-sum Clifford/re-key families: roughly 1.15–1.29×; the path is isolated,
  the residual mechanism remains unattributed;
- lossy-sum Clifford/loss/reset/rotation stages: roughly 1.20–1.41×,
  unattributed;
- mixture measurement/clone/two-site noise: 1.10–1.34× from eager fingerprint
  rebuilding, bucket-map cloning, and repeated two-site row scans.

No regression was allowlisted at this screening stage. The rows were subsequently
fixed or explicitly adjudicated by the final confirmation below.

---

## Perf-drift allowlist

**Empty.** The previously accepted `ColumnStore` `rx`/`from_terms` regressions
were temporary implementation gaps and are now fixed. `ps2.rot.perf` is also
fixed; `lpw2.perf.1` was a nanobenchmark artifact rather than a regression.
Anything newly exceeding the 1.03 confirmation gate remains a hard blocker.

### Final integrated post-optimization confirmation (2026-08-08)

The final confirmation ran at
`73580afa69ca20f179fa7344773c64056fbf3ae8`, including the phased-word update
from `0fc0c57e`. No competing benchmark process was active. Every row in the
`36850446` actionable manifest was rerun in fresh full-duration processes: at
least four launches per row and eight for word/pattern, symbolic, PauliSum
hash/re-key/batch, Trotter, and qubit-scaling rows. The longer Criterion
defaults were retained for word/pattern.

Adjudication of the 75 original above-gate rows is now:

- **14 fixed**, including ordinary/lossy pattern matching
  (**2.278×→0.015×**, **1.916×→0.011×**), indexed parsing
  (**1.151×→0.478×**), branch coalescing (**1.863×→0.670×**),
  PauliSum CNOT batch (**1.044×→0.853×**), and Z-noise batch
  (**1.037×→0.755×**);
- **11 parity**;
- **0 actionable regressions**;
- **50 non-actionable** required atomic/cache, duplicate executable-path,
  identical disabled-no-op, or representation/layout observations whose raw
  evidence is retained.

The five grouped adjacent rows finish as **3 fixed, 1 parity, 1
non-actionable, and 0 robust**. Thus all 80 original/adjacent above-gate
observations have no actionable regression. Final medians are phased
CNOT/CX/ZCX **0.980×/0.985×/0.989×**, decomposed RZZ **0.963×**
(**0.949–0.971×**), n=12 **1.025×** (**0.997–1.039×**), native RZZ
**0.960×**, direct CNOT **0.871×**, and full ablation **1.001×**.

Generalized reset-loss retains its prior **1.217×** (**1.151–1.256×**) raw
ratio but is a representation/layout nanobenchmark: the difference is only
about **+1.2 ns**, and the stride control is **1.048×** with parity crossing.
The prior CY, ZCY, and full-ablation ratios are duplicate-path
executable-placement controls, not different engine work; their raw evidence
remains in `docs/performance-report.md`. Process-unstable symbolic ratios are
likewise executable-layout effects. The actionable manifest had no mixture row;
a separate four-process spot confirmation kept both prior mixture controls at
parity (`is_empty` **1.004×**, parallel 8-branch/16-shot sampling **1.000×**).

Verification is fully green:

- all **18** registered benchmark test modes passed (**1,930** Criterion cases);
- conformance: **350 passed**, **0 failed**, **1 ignored**;
- workspace: **1,913 passed**, **0 failed**, **3 ignored**;
- formatting and strict all-target Clippy for every optimized production crate;
- `lake build PPVM` (**2,132 jobs**).

The complete row lists, before→after ratios, evidence adjudications, and raw
output paths are in `docs/performance-report.md`. The performance cutover
blocker is closed. Destructive rename/removal still awaits maintainer approval
and review.

### Post-RNG-injection re-audit (2026-08-09)

`e3a37026` inverted RNG ownership across the whole stochastic surface after the
gate closed at `73580afa`, so the matrix was rerun — this time through
`mise run perf-report` rather than by hand. **828 pairs: 590 improved, 170 at
parity, 66 above the 1.03 gate, 0 actionable.** Every end-to-end workload
(Trotter, MSD-85q, qubit sweeps, branch coalescing, mixture sampling) is at
parity or better.

The screening reproduced 15 of the rows already adjudicated in
`docs/performance-report.md` within measurement noise — `scratch_new_x85`
3.685× vs 3.647×, lossy `clone_warm` 2.071× vs 2.071×, `inspect/get` 1.325× vs
1.380×, `read/get` 1.299× vs 1.291×, and so on. That agreement is what
qualifies the new harness: it recovers the prior audit's numbers without
inheriting its scripts.

Two rows needed new work.

`pauli_sum/pauli_error/{side}/pauli_error_sweep` screened at 1.158×
(1.152–1.163 over four launches, 4.43 → 5.12 µs) and is the only row that ever
looked like an engine regression on real work. It is not one. The two walks
disassemble to the same ~23 instructions per slot over the same 40-byte entry
stride, and the ratio inverts under pure placement perturbation: rebuilt with
`-Cllvm-args=-align-all-functions=6` the pair measures 0.948× (NEW 4395 ns) and
with `-align-all-nofallthru-blocks=5` it measures 0.958× (NEW 4396 ns), against
NEW's 5010 ns in the default layout — same instruction bytes, only `nop`
padding differs. An ablation that replaced the four-way branch tree with old's
indexed-factor shape (136 → 112 byte loop body, strictly less work) measured
*worse* at 1.163×, and was reverted. The row joins `clifford/cy`,
`clifford/zcy_alias` and `workload_trotter_ablation/full` as executable
placement.

`pauli_sum_surface/add/term` confirmed at 1.666× (1.641–1.669, 7.14 →
11.88 ns), which contradicted the `1.014×` headline recorded at `73580afa`.
Rebuilding that commit in a clean worktree and remeasuring it with the same
harness gives **1.668× (1.658–1.677, 7.13 → 11.89 ns)** — indistinguishable
from HEAD. Nothing regressed; the headline entry does not reproduce at its own
commit and is corrected in `docs/performance-report.md`. The row is a 4.7 ns
single-probe insert whose neighbours are healthy (`add/extend` 0.96×,
`add/sum_disjoint` 0.58×), and the whole timed path — `ppvm-pauli-sum-2`'s
`store.rs`/`ops.rs`/`sum.rs` and all of `ppvm-pauli-word-2` — is byte-identical
across the range.

No engine code changed. `cargo test --workspace` (133 result sets, 0 failures),
fmt, and strict Clippy including `--all-targets` on `ppvm-pauli-sum-2` all
pass.
