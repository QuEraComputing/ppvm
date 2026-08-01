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

**Perf — CORRECTED (an earlier version of this entry was wrong).** `pauli_error`
**0.89×** (faster). `rotation_rx` is a **real ~1.15–1.17× residual regression**
(new stable ~5.84µs vs old ~5.0µs), **over the 1.15× gate** — NOT "at parity" as
first written. The failed run's test-agent 2.42× was **not noise**: it was a real
measurement of the *batch/`apply`* path (see `store.rs` `RotateInPlace` friction
note, "the batch path benchmarked ~2.4× slower"), which the impl replaced with the
fused single-pass `RotateInPlace` → the ~1.15× residual. Attribution verified: the
`iter_batched_ref` per-iteration `clone` is symmetric (new-sum clone 1.36µs ≈
old-sum clone 1.41µs), so it is genuine rx cost — the branch path (per
anticommuting term: key `clone` + `set_bit` + cold `key_hash` + hashmap merge), a
modest partial materialization of the lazy-hash/`Copy`-drop word cost on the
branch-key-creating path.

**Optimized (maintainer chose "optimize").** Attribution (instrumented phase
timing, release): the obvious suspect — the per-call `branches` Vec allocation —
is **negligible** (67ns); cost splits ~evenly between pass-1 (branch-key build)
and pass-2 (hashmap hash+merge). The one cheap lever: the branch build did
`clone()` (which *loads* the source's cached digest) then `set_x_bit` (which
*stores* `HASH_UNCACHED`) — a wasted atomic load+store per branch key. Added
`PauliWord::with_bits_toggled` (copy planes, flip the bit, uncached hash) →
new/rx ~5.84µs → **~5.65µs**, ratio **~1.19× → ~1.16×**. The residual sits **at
the 1.15× gate boundary** and is inherent: the dominant remaining cost is the
cold `key_hash` finalization fold on each freshly-created branch key (pass-2),
which the old crate's `Copy`/eager-hash word pays a bit more cheaply — not
removable without changing the shared avalanche-tested hash.

**Pushed deeper (maintainer asked twice); the hash theory was DISPROVEN.** Ruled
out with evidence: per-call Vec alloc (67ns); a mid-merge resize (pre-reserve gave
no change); **and — critically — the hash is NOT the cause.** Head-to-head: the
new map uses `IdentityBuildHasher` (`key_hash()` passes straight through), while
the OLD map uses `FxBuildHasher` (an *extra* fxhash round on the cached `u64` per
lookup); and new's `with_bits_toggled` branch build is cheaper than old's
`clone()`+2×`set_bit`+eager `rehash()`. So new's hash+build path is **cheaper**,
not more expensive — my earlier "inherent lazy-hash cost" justification was wrong.
Also ruled out: truncate (`NoPolicy` no-op). **Root cause NOT pinned:** samply
would not symbolicate the binary. Most likely **codegen/monomorphization** (the
generic `Sum<S,P>` + RotateInPlace/closure/Policy layering optimizing slightly
worse than the old *concrete monolithic* `PauliSum<Config>` loop), possibly
compounded by the noisy old/rx baseline (4.9–5.5µs → the gap is ~10–16% depending
on the run). **Resolved: cheap win captured (`with_bits_toggled`); residual
~1.14× allowlisted** — a small, at-the-boundary, **non-algorithmic** difference
(the algorithmic path is demonstrably cheaper), exact source unpinned (codegen +
noise), NOT an inherent design cost. A future symbolicated profiling pass could
name it if it ever matters.

Lean added (prove agent, `lake` green, orchestrator-verified non-vacuous):
`Rotation.lean` `branchExp`/`branchExp_isRealPhase` + `rx/ry/rz_eps_from_product`
(the branch ±1 sign is **real** — the `i` of `iGP` cancels the product's `i` when
`{G,P}=0` — with the exact per-axis ε formula pinned); `Conjugation.lean`
`conjX`/`conjY`/`conjZ` (the pure-sign X/Y/Z conjugation, witnessing the in-place
fast-path signs).

| id | type | sev | routed | status | note |
| --- | --- | --- | --- | --- | --- |
| ps2.noise.1 | correctness | med | impl (orch) | closed | `pauli_error` with a zero transfer eigenvalue (reachable, `[0,0.25,0.25]`→λ_X=0) left a **phantom zero-coefficient key** in the support, violating the reduced-canonical-form invariant. Orchestrator fixed `ScaleByKey` (both backends) to drop any term scaled to exactly zero; added a `pauli_error_zero_eigenvalue_drops_the_term` test. |
| ps2.rot.perf | perf-drift | med | human | **allowlisted** | `rotation_rx` ~1.13–1.16× (branch path). The 2.42× was the batch path (real, fixed by `RotateInPlace` fusion). Optimized: `with_bits_toggled` (skip the wasted clone-hash-load+invalidate) → ~1.19×→~1.16×; alloc/resize ruled out empirically. Residual is inherent lazy-hash/`Copy`-drop cost on the key-creating path; the only bigger lever (Zobrist incremental hash) is a risky rewrite of the shared hash, not justified for a boundary case. Allowlisted. |
| ps2.nsites.1 | impl-friction | low | impl | deferred | The `n_sites` `debug_assert!` is on `from_terms` but not on the `apply`/`rekey_bijective` produced-key paths. No live bug (every producer preserves width); low, deferred. |
| ps2.xyz-dup.1 | impl-friction | low | accepted | noted | X/Y/Z sign logic (`(−1)^z`, etc.) is written a third time in the in-place fast path (besides `PhaseTrack` and `Phased`). Now has a Lean witness (`conjX/Y/Z`) and differential tests guard it; accepted — the in-place path exists precisely to avoid the `Phased` wrapper's cost. |

**Infra note:** the workflow *failed* (a reporting agent hit the StructuredOutput
retry cap — a transient JSON-schema failure), but 4/5 agents completed and the
implementation + tests were fully on disk and green. No resume needed; the
orchestrator verified, fixed the noise bug, and finished.

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
| ps2.rot.perf | ppvm-pauli-sum-2 | `rotation_rx` (branch fan-out) | ~1.14× (new ~5.65µs vs old ~4.9–5.5µs, noisy) | Deeply investigated. **NOT the hash** (new's `IdentityBuildHasher` pass-through is cheaper than old's `FxBuildHasher` extra round; `with_bits_toggled` cheaper than old clone+rehash), NOT alloc/resize/truncate. Cheap win captured (`with_bits_toggled`). Root cause unpinned (samply wouldn't symbolicate); most likely codegen/monomorphization of the generic `Sum` engine + measurement noise — a small, at-boundary, **non-algorithmic** difference, not an inherent design cost. | maintainer (optimize→accept after deep dive) |
| _(prior)_ | | | | `lpw2.perf.1` was a nanobench artifact (converges to 1.00× at wide width), not allowlisted. | |
