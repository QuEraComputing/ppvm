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

## ▶ Continue here — handoff (2026-08-01)

Everything committed and green (whole workspace tests, clippy/fmt clean, `lake build PPVM` ok). Recent commits (newest first):
`b68b7e8e` content-conditional hash invalidation + cnot/truncate microbenches · `d4aa7875` truncation is explicit-only (behaviour gap closed) · `8a1dc79b` integration coverage + workflow Baseline phase · `4385ef3c` storage aux double-buffer recovered · `33f48fe5` Mytkowicz same-build perf rule.

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

### ⚑ NEW behaviour-parity gaps found by the Baseline agent (verified by hand)
Three divergences from old, all **verified directly against the old sources**,
none yet fixed. Under the PRIME DIRECTIVE each is a gap.

1. **`ps2.zero.behaviour` (high).** The old crate has **no `reduce`** and no
   drop-zero logic *anywhere* — a zero-coefficient entry stays in the support
   forever, and old's `PartialEq` compares maps exactly, so
   `ppvm-pauli-sum/tests/loss.rs::test_reset_channel` depends on `*= 0.0`
   keeping every key. The new `reduce()` (`ppvm-traits-2/src/containers.rs`)
   **and** `ScaleByKey` both drop exact zeros, so `len()`/`contains` diverge —
   e.g. after `pauli_error(q, [0.0, 0.25, 0.25])` (λ_X = 0).
   **Note this indicts an earlier decision in this very log:** the `ps2.noise.1`
   row records the orchestrator *adding* zero-dropping to `ScaleByKey` to uphold
   a "reduced-canonical-form invariant". That invariant is the new design's
   invention, not old's behaviour — so that fix introduced a divergence rather
   than removing one. Re-decide deliberately: either restore old's retain-zeros
   semantics, or get the Lean oracle to say old is wrong and document it.
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

### 2 — DO NEXT: `ps2.cnot.rekey.perf` — still OPEN, partly reduced
End-to-end Trotter is **~1.37–1.42×** (was ~1.47× at session start; criterion,
same build, `[u8;8]`). Numerics match. The gap is now concentrated almost
entirely in `cnot` / `RekeyBijective`.

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
planes, two component `OnceLock<u64>` caches (loss-only mutation skips X/Z
rehash); `Word<Site=LossySite<Pauli>>`, `PauliBits`(`is_lost`),
`SymplecticColumns`+`PhaseTrack` (loss-guarded, phase-discarding), `Indexable`,
`Columnar`; inherent `set_lost`/`clear_loss`/`loss_weight`. 27 in-crate tests.
Conformance: `lossy_pauli_word_diff.rs` (9), `lossy_pauli_word_lean.rs` (9),
lossy bench. Gates verified by orchestrator (build/clippy/fmt/test across all four
`-2` crates, machete, `lake build PPVM`) — all green; the `PauliBits` relaxation
did not regress `ppvm-pauli-word-2` (27 tests still pass).

**Perf gate: 6/6 pass (1 flag investigated → benign).** product **0.79×**, cnot
**0.78×**, key_hash-warm **1.02×**, weight **1.02×**; key_hash-cold 15.6ns is the
design-accepted lazy-OnceLock trade-off. `loss_weight` flagged at 1.28× at single-
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
| ps2.zero.behaviour | correctness | high | design | **OPEN — new** | Old has **no `reduce`** and no drop-zero logic anywhere; a zero-coefficient term stays in the support and old's exact-map `PartialEq` depends on it (`ppvm-pauli-sum/tests/loss.rs::test_reset_channel`). New `reduce()` and `ScaleByKey` both drop exact zeros ⇒ `len()`/`contains` diverge (e.g. `pauli_error(q,[0.0,0.25,0.25])`, λ_X=0). Verified by hand against the old sources. **Indicts `ps2.noise.1`**: that row's "fix" *added* zero-dropping to satisfy a reduced-canonical-form invariant which is the new design's invention, not old's behaviour — so it introduced this divergence. Re-decide: restore old semantics, or have Lean adjudicate old as wrong and document it. |
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

_(no entries yet)_

## Phase 7 — downstream + cutover

_(no entries yet)_

---

## Perf-drift allowlist (human-approved regressions)

Regressions listed here are accepted as designed-in trade-offs and do **not** block
a crate. Anything not listed is a hard gate.

| id | component | metric | accepted ratio | justification | approved-by |
| --- | --- | --- | --- | --- | --- |
| ps2.rot.perf | ppvm-pauli-sum-2 | `rotation_rx` | **FIXED — not allowlisted (now parity)** | Was a real ~1.15× (fair same-storage); root-caused (lazy branch-key hash firing on pass-2's probe critical path) and fixed by warming the digest in `RotateInPlace` pass 1 (semantic no-op, 108 tests green). `new/rx` ~5.5→~4.9µs; interleaved ratio ~0.99. Nothing to allowlist. |
| _(prior)_ | | | | `lpw2.perf.1` was a nanobench artifact (converges to 1.00× at wide width), not allowlisted. | |
