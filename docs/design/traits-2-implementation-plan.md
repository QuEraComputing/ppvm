# `ppvm-*-2`: implementation & migration plan

Status: plan

This plan turns the design in
[`traits-2-configuration-and-hashing.md`](traits-2-configuration-and-hashing.md)
(and the companion `word-data-structures.md` / `tableau-data-structure.md`) into a
sequence of new `ppvm-*-2` crates. The Lean development in [`../../lean`](../../lean)
is the mathematical specification: its machine-checked identities are reused as
test oracles.

## Guiding principles

1. **Additive, not in-place.** Every new crate is `ppvm-xxx-2`. No existing crate
   is modified until the final cutover. Old and new compile side by side in the
   same workspace.
2. **Lean is the spec.** The Lean theorems pin the semantics (phase cocycle,
   twisted-product associativity, symplectic isometries, rotation-as-2D-rotation,
   XOR-relabel bijection, truncation bounds, Pauli-channel eigenvalue). Each is
   reproduced as a Rust property test in the corresponding `-2` crate. Several
   Lean docstrings already cite the exact Rust source lines they refine
   (e.g. `crates/ppvm-pauli-word/src/phase/mul.rs:42`), which we carry forward.
3. **Differential testing against the old crates.** Behavioral equivalence with
   the current implementation is verified by a dedicated, test-only conformance
   crate that depends on *both* the old and new crates. Shipped `-2` crates never
   depend on old crates.
4. **Performance gate.** A crate is not "done" until Criterion benchmarks show its
   hot paths are at parity with or faster than the old crate on the same inputs.
5. **Cutover last.** Rename/replace happens only after all `-2` crates are green,
   at-parity, and downstream consumers have a shim. Until then the old crates
   remain the production path.

## Crate map (old → new)

| Old crate | New crate | Core content |
| --- | --- | --- |
| `ppvm-traits` (`Config`, `Coefficient`, `Word`/`PauliWordTrait`, `ACMap*`, `Strategy`, gate traits) | **`ppvm-traits-2`** | Trait *definitions* (+ the `BlanketClifford` blanket and, forced by the orphan rule, the graded-trait container impls in `containers.rs`): split `Coefficient`, `Angle`; `Word`/`Indexable`/`PauliBits`; `SymplecticColumns`/`PhaseTrack`/`StabilizerFrame`; gate/noise traits (`Clifford`/`CliffordExtensions` + the batched forms, `RotationOne`, `Reset`, `Measure`, `PauliError` and the channel family); graded map layers `Support`/`Accumulate`/`Scale`/`Pair`/`Multiply` (impl'd on `Vec`/`HashMap` here, since both trait and container are foreign to the sum crate); algebra capabilities `KeyProduct`/`ImaginaryUnit`/`Conjugate`; batch contract `Columnar`/`KeyColumn`/`KeyBatch`/`TermBatch`/`TermSink`/`TermProducer`; `IdentityHasher`; small concrete leaf types (`Pauli`, `Phase`, `LossySite`). |
| `ppvm-pauli-word` (word/) | **`ppvm-pauli-word-2`** | Concrete `PauliWord`; private packed X/Z planes + lazy relaxed-atomic sentinel hash cache + `HashFinalize`/`PauliStorage`/`PauliKeyColumn` (re-exported for the lossy/phased crates to reuse). Impls of `Word`/`Indexable`/`PauliBits`/`SymplecticColumns`/`PhaseTrack`/`BlanketClifford`/`KeyProduct`/`Columnar`. |
| `ppvm-pauli-word` (phase/) | **`ppvm-phased-pauli-word-2`** | The generic `Phased<W>` wrapper + `PhasedPauliWord = Phased<PauliWord>`; carries an explicit ℤ₄ phase and a hand-written **fused** `impl Clifford` (reads each inner X/Z bit once via `W: PauliBits`, computes the ℤ₄ sign, applies the bit update, folds in the sign) — *not* the blanket, so it does **not** implement `BlanketClifford`. Delegates `Word` to `W`, phased product via `W`'s `KeyProduct`. **Non-indexable.** |
| `ppvm-pauli-word` (loss/) | **`ppvm-lossy-pauli-word-2`** | `LossyPauliWord` — a *distinct* concrete Pauli-word impl in its own crate (adds a packed loss plane). Reuses `ppvm-pauli-word-2`'s packed-storage/hash infra. `Word<Site = LossySite<Pauli>>` + `PauliBits` (`is_lost`) + phase-discarding `SymplecticColumns`/`PhaseTrack` + `BlanketClifford` + `Indexable`; loss writes and `loss_weight()` inherent. |
| `ppvm-pauli-sum` | **`ppvm-pauli-sum-2`** | `Sum<S, P>`; graded traits `impl`'d on `Vec<(K,C)>` and `HashMap<K,C,IdentityBuildHasher>`; `Policy` + `NoPolicy`/`MaxPauliWeight`/`CoefficientThreshold`/`CombinedPolicy` + `Retain`; `TermProducer` impls (`RekeyProducer`, rotation/noise producers); gate/rotation/noise trait impls; `PauliSum`/`LossyPauliSum`/`FermionSum` aliases; `IdentityBuildHasher` + `HashMapStore` aliases. |
| `ppvm-tableau` | **`ppvm-tableau-2`** | `Tableau` (`Indexable` + `SymplecticColumns` + `PhaseTrack` + `StabilizerFrame` + `Clifford` + `Measure`); `GeneralizedTableau` (`frame: Tableau` + `amplitudes: Sum<Vec<(Bitstring, Complex<C>)>, CoefficientThreshold>`). |
| `ppvm-tableau-sum` | folds into **`ppvm-tableau-2::mixture`** as `GeneralizedTableauMixture` (`GeneralizedTableauSum` compatibility alias) | A specialized probability mixture over full generalized tableaux. Fingerprint buckets narrow candidates, then collision-checked structural comparison includes frame/loss and approximately-equal amplitudes. It cannot be a normal `HashMap` key because approximate amplitude equality is not an equivalence relation. |
| `ppvm-sym` | **`ppvm-sym-2`** | Exact/symbolic coefficient ring implementing the loosened bounds: `Coefficient` (no `Mul<f64>`), `ImaginaryUnit` (`i²=−one()`), `Conjugate`. Doubles as the concrete witness that L4 admits exact rings (mirrors the Lean `GaussianInt` instance). |
| `ppvm-stim`, `ppvm-vihaco`, `ppvm-cli`, `ppvm-python-native`, `ppvm-tui` | unchanged until cutover | Downstream consumers; adapters added in Phase 7. |
| (new, test-only) | **`ppvm-conformance-2`** | Not shipped. Depends on old + new. Houses cross-crate differential tests, shared random-circuit/random-word generators, and comparative benchmarks. |

### Dependency edges (new crates)

```
ppvm-traits-2  ──►  (nothing ppvm)
ppvm-pauli-word-2  ──►  ppvm-traits-2
ppvm-lossy-pauli-word-2  ──►  ppvm-traits-2, ppvm-pauli-word-2 (reuses packed-storage/hash infra)
ppvm-phased-pauli-word-2 ──►  ppvm-traits-2, ppvm-pauli-word-2 (wraps PauliWord; PhasedPauliWord alias)
ppvm-pauli-sum-2   ──►  ppvm-traits-2, ppvm-pauli-word-2
ppvm-tableau-2     ──►  ppvm-traits-2, ppvm-pauli-word-2, ppvm-pauli-sum-2
ppvm-sym-2         ──►  ppvm-traits-2, ppvm-pauli-sum-2
ppvm-conformance-2 (dev/test) ──►  every *-2 AND every old counterpart
```

This is the old graph with `Config` removed and `ppvm-tableau-sum` collapsed.

## Phase 0 — scaffolding and the equivalence harness

Goal: make "implement → prove against Lean → diff against old → benchmark" a
repeatable loop before writing real logic.

- Create `crates/ppvm-traits-2` (empty lib) and `crates/ppvm-conformance-2`
  (test-only: `publish = false`, all logic under `#[cfg(test)]` / `benches/`).
- Add both to `[workspace].members`.
- In `ppvm-conformance-2`, build the shared generators once:
  - random `PauliWord` of width `n` (paired old/new constructors from the same
    seed so both crates see identical operators);
  - random Clifford+rotation circuits (gate, qubit(s), angle) as data, replayable
    against any backend;
  - a `same_seed` RNG helper so old and new consume identical randomness.
- Establish the two test macros/utilities used everywhere below:
  - `assert_matches_old!(new_result, old_result)` for differential tests;
  - a Criterion group template that runs the *identical* workload on old and new
    and reports the ratio.
- CI: `cargo build --workspace --all-targets` already covers new crates;
  add `cargo test -p ppvm-conformance-2`. Benchmarks run on demand, not in the
  gate (but a regression check can be added later).

Done when: an empty `PauliWord`-vs-`PauliWord` identity diff test and a no-op
benchmark both run green.

## Phase 1 — `ppvm-traits-2` (the trait definitions)

Port the design's trait surface, one module per design section. No heavy logic;
the only concrete code is leaf types and a couple of blanket impls.

Modules:
- `coefficient.rs` — `Coefficient` (ring + `mul_sign`/`magnitude`, **no**
  `Mul<f64>`, **no** `half`), `Halvable: Coefficient` (`half`; the partial `0.5·x`
  the measurement projector needs, kept off `Coefficient` so exact rings qualify),
  `Angle<C>` (+ `impl Angle<f64> for f64`). Impl `Coefficient` + `Halvable` for
  `f64` and `Complex<f64>`.
- `algebra.rs` — `KeyProduct` (`key_mul -> (Self, Phase)`), `ImaginaryUnit`
  (`imaginary_unit()`, law `i·i == −one()`), `Conjugate` (`conj`). Impl
  `ImaginaryUnit`/`Conjugate` for `Complex<f64>`; `Conjugate` (identity) for `f64`.
  `Phase` enum (`{1, i, −1, −i}` ≅ ℤ/4).
- `word.rs` — `Word` (`Site`, `n_sites`, `get`, `weight`, `iter`), leaf types
  `Pauli`, `LossySite<S>`, `FermionSite`. `PauliBits: Word` (supertrait relaxed
  from `Word<Site = Pauli>` so `LossyPauliWord`, whose `Site` is
  `LossySite<Pauli>`, can implement it; propagation re-adds `Site = Pauli`).
- `pauli.rs` — `SymplecticColumns`, `PhaseTrack`, `StabilizerFrame`, the opt-in
  marker `BlanketClifford`; the blankets
  `impl<T: SymplecticColumns + PhaseTrack + BlanketClifford> Clifford for T` and
  the matching `CliffordExtensions` blanket (each extension gate derived as a
  product of audited generators, so no new phase primitives are needed)
  (the marker keeps the blankets coherence-legal alongside `Phased<W>`'s fused
  overrides).
- `gates.rs` — `Clifford`, `CliffordExtensions`, `CliffordBatch`,
  `CliffordExtensionsBatch`, `RotationOne<C, A = C>` (required `rotate_1`,
  defaulted `rx`/`ry`/`rz` + `*_many`), `RotationTwo<C, A = C>`,
  `RotXY<C, A = C>`, `CRx<C, A = C>`, `U3Gate<C, A = C>`, `TGate`,
  `Projection`, `Measure`, `Reset`,
  `PauliError<C>` (with the stim `x_error`/`y_error`/`z_error` + `*_many`
  defaults), and the channel family `PauliErrorAll`, `TwoQubitPauliError`,
  `Depolarizing`, `Depolarizing2`, `AmplitudeDamping`, `LossChannel`,
  `CorrelatedLossChannel`, `ResetLossChannel`, `AsymmetricLossChannel` — the old
  crate's behavioral surface, `<T: Config>` → `<C: Coefficient>`.
- `graded.rs` — `Support`, `Accumulate`, `Scale`, `Pair` (with both `overlap` and
  `hermitian_overlap where Coeff: Conjugate`), `Multiply`, `Retain`.
- `batch.rs` — `Columnar`, `KeyColumn`, `KeyBatch`, `TermBatch`, `TermSink`,
  `TermProducer`.
- `hash.rs` — `Indexable` (`key_hash`), `IdentityHasher`, `IdentityBuildHasher`.

Lean-oracle tests (unit): implement a trivial reference `Coefficient`/`Conjugate`
and check the stated laws (`i·i == −one()`, `conj∘conj == id`, `conj(i) == −i`)
against `ppvm-sym-2`'s exact ring in Phase 5 — for now, check `f64`/`Complex<f64>`.

Done when: crate compiles; trait method signatures match the design verbatim; the
`Clifford` blanket impl type-checks against a stub `SymplecticColumns + PhaseTrack`
that opts into `BlanketClifford`.

## Phase 2 — `ppvm-pauli-word-2` (the fundamental data structure)

Implement `PauliWord` first, then loss, then the phased wrapper. This is where the
bulk of the correctness-critical, Lean-validated logic lives.

Order:
1. **`PauliWord`** — private packed X/Z planes (const-generic width internally;
   `PauliStorage` is *not* re-exposed) + lazy relaxed-atomic sentinel hash cache.
   Impl `Word<Site = Pauli>`, `PauliBits`, `SymplecticColumns`, `PhaseTrack`,
   `BlanketClifford` (opt into the shared blanket `Clifford`),
   `KeyProduct` (the twisted product `v·w = iᵏ(v⊕w)`), `Indexable` (`key_hash`
   with the private `HashFinalize` fold), `Columnar`.
2. **`LossyPauliWord`** — adds the loss plane; `is_lost` override, `loss_weight`,
   inherent loss mutation; `Word<Site = LossySite<Pauli>>`.
3. **`Phased<W>` / `PhasedPauliWord`** — non-indexable phased wrapper with a
   hand-written *fused* `impl Clifford` (read-once bits + ℤ₄ sign; opts out of
   `BlanketClifford`).

Lean-oracle tests (these are the crate's backbone):
- **Phase product** — port `phaseExp`'s boolean `sign`/`imag` formulas exactly
  (Lean `Phase.lean` `phaseExp_eq_ref`, grounded in `Matrix.lean` `pauliMat_mul`);
  test `key_mul` against the ℤ[i] matrix reference for all single-qubit cases and
  randomized n-qubit words (n-qubit sum = `Word.lean` `phaseExpN_cocycle`).
- **Associativity / group laws** — `(u·v)·w == u·(v·w)`; `P·P == I` up to phase
  (`phaseExpN_self`); `P·Q == (−1)^{ω} Q·P` (`phaseExpN_sub_comm`).
- **Clifford conjugation** — `H`/`S`/`CNOT`/`CZ` bit rules; `HXH=Z`, `HYH=−Y`,
  `SXS†=Y`, … (Lean `Conjugation.lean`), and symplectic-form preservation
  `ω(gPg†, gQg†) == ω(P,Q)` (`Symplectic.lean` isometries).
- **Hash contract** — `Hash` writes exactly `key_hash()`; equal words ⇒ equal
  digest; a distribution property test for avalanche quality (the design's stated
  contract, not a type-level guarantee).

Differential tests vs old `ppvm-pauli-word`: for random words, assert equal
multiplication result+phase, equal Clifford results, equal weight/iter, and equal
*bucket distribution* (not necessarily equal raw digest — the fold may differ).

Benchmarks vs old: `mul`, single-gate Clifford conjugation, `key_hash` (cold and
cached), `weight`. Gate: parity or better. Watch the lazy-cache vs old eager
`hash_cache` + `Copy` trade-off (the design accepts losing `Copy` for correct lazy
caching — confirm the throughput cost is acceptable here).

Done when: all Lean-oracle tests pass, differential tests pass, benches ≥ parity.

## Phase 3 — `ppvm-pauli-sum-2` (the graded engine)

1. **Graded traits on containers** — `impl Support/Accumulate/Scale/Pair` for
   `Vec<(K,C)>` (linear scan, `K: Eq + Clone`) and for
   `HashMap<K,C,IdentityBuildHasher>` (`K: Indexable`). `Multiply` where
   `K: KeyProduct, C: ImaginaryUnit`. `Retain` for both.
2. **`Sum<S, P>`** — storage + policy + `n_sites`, the `apply<TP: TermProducer>`
   method (read → produce into `TermSink` → `accumulate_batch` → `reduce` →
   `policy.truncate`).
3. **Producers** — `RekeyProducer` (Clifford pushforward), rotation/noise
   branch producers.
4. **Gate/rotation/noise impls on `Sum`** — blanket `Clifford for Sum<S,P>` via
   per-key `Clifford`; `RotationOne`, `PauliError`, `Measure`.
5. **Policies** — `NoPolicy`, `MaxPauliWeight` (uses `Word::weight`),
   `CoefficientThreshold` (`magnitude() >= threshold`), `CombinedPolicy`.
6. **Aliases** — `IdentityBuildHasher`, `HashMapStore<K,C>`, `PauliSum`,
   `LossyPauliSum`.

Lean-oracle tests:
- **Module laws** — `accumulate` commutative/associative; `reduce` drops exactly
  zero coefficients (`GradedMap.lean` `reduce_structural`); `scale` distributes.
- **L4 multiply** — basis-monomial product `= key_mul` with the `iᵏ` folded onto
  the coefficient (`Twisted.lean` `tmul_assoc` — associativity of `multiply_into`
  over `Complex<f64>` and, in Phase 5, the exact ring); and right-multiplication
  by one word is the aggregation-free injective re-key `Sum::mul_word_assign`
  takes (`Twisted.lean` `twistedConv_single_right`,
  `twistedConv_single_right_apply`).
- **`overlap` / `hermitian_overlap`** — symmetry of `overlap` (`overlap_comm`),
  Pauli-basis orthonormality (`Noise.lean` `overlap_single_single`), and for
  complex amplitudes conjugate symmetry + PSD of `hermitian_overlap`
  (`GradedMap.lean` `hermitianOverlap_*`).
- **Rotation branch** — `c·P → cos·c·P + sin·c·P'` with `P' = iGP` distinct
  (`Rotation.lean` `anticommute_new_key`), norm preserved (`rot_norm_sq`), merge
  additive (`rot_rot`); and the whole-map two-pass `RotateInPlace` fast path
  equals the one-pass produce-and-accumulate batch for every walk order
  (`Rotation.lean` `accumulate_rotBatch`, with `eagerWalk_ne_twoPass` showing an
  interleaved single pass is observably different).
- **Truncation** — dropping-error `ℓ¹` bound (`Truncation.lean` `l1_bound`) as a
  runtime assertion in a randomized test; document the `≥` keep-rule vs the
  tableau's `>` (`cutoff_mismatch`).
- **Noise** — unital Pauli channel eigenvalue `λ_P = 1 − 2Σ_{anti} p_Q`
  (`Noise.lean` `pauli_channel_eigenvalue_omega`), and its contractivity
  `|λ_P| ≤ 1` on a sub-stochastic probability vector
  (`pauli_channel_eigenvalue_abs_le_one`, `l1_contractive`,
  `scaleByKey_support_subset`) — the licence for the channel fast path skipping
  truncation and the weight re-check.

Differential tests vs old `ppvm-pauli-sum`: replay identical random circuits on
old `PauliSum` and new `PauliSum`; after each gate assert equal support and
per-key coefficients (within tolerance for `f64`). Include rotation branching and
truncation with matched thresholds (accounting for the documented boundary rule).

Benchmarks vs old: Clifford propagation, rotation fan-out, `multiply_into`
(observable squaring), `overlap`, and a representative TFIM/Trotter sweep (the
existing `benchmarks/plot_*` scripts target these). Gate on parity.

## Phase 4 — `ppvm-tableau-2` and the mixture

- **`Tableau`** — `SymplecticColumns` + `PhaseTrack` (ℤ₂ sign + Aaronson–Gottesman
  `g`-rule) + `StabilizerFrame` + blanket `Clifford` + `Measure` (pure-Clifford
  algorithm) + `Indexable`.
- **`GeneralizedTableauMixture<A, I, H>`** — a specialized probability
  distribution over full `GeneralizedTableau` states, exported under the old
  `GeneralizedTableauSum` name for compatibility. It retains the old
  fingerprint-bucket + collision-check design: merge identity includes frame
  rows/phases, loss flags, and amplitude vectors equal within the coefficient
  threshold, while excluding RNG/record/scratch/probability state. It is
  deliberately **not** `Sum<HashMapStore<Tableau, _>, _>`: bare frames omit
  amplitudes/loss, and approximate amplitude equality is non-transitive, so it
  cannot lawfully be `Eq`/`Hash`.
- **`GeneralizedTableau`** — `frame: Tableau` + `amplitudes: Sum<Vec<(Bitstring,
  Complex<C>)>, CoefficientThreshold>`; Clifford updates the frame only,
  non-Clifford branches the amplitude `Sum`, `O(n²)` measurement is tableau-local.

Lean-oracle tests: symplectic-frame invariant preserved by every gate
(`Frame.lean` `isSymplecticFrame_*`) **and by the measurement projection**
(`Frame.lean` `isSymplecticFrame_projectFrame` — the justification for
`canonicalize` being a no-op); frame generators independent
(`frame_linearIndependent`) **and spanning** (`frame_surjective`), so
`compute_decomposition`'s anticommutation bitmasks are exactly the frame
coordinates (`Frame.lean` `frame_coordinate_expansion`, Yoder-2012 Lemma 5);
measurement dichotomy / deterministic⇔`X`-free
(`measurement_dichotomy`, `measure_deterministic_iff_xfree`); Clifford leaves
amplitudes fixed and the XOR relabel is a bijection (`Bitstring.lean`); case-a
measurement is the Born projector — the `ℤ/4` overlap sign table, `M² = I`,
`M† = M`, `P₀ + P₁ = I`, `P_b² = P_b`, `prob_1 = ⟨c, P₁ c⟩`, and
keep-`A`/transform-`B` `= 2·P_b` (`Projection.lean` `rustTerm_eq`,
`shiftOp_involutive`, `shiftOp_selfAdjoint`, `overlap_eq_inner`, `proj_add`,
`proj_idem`, `probOne_eq`, `projectRaw_eq_two_proj`); **case b** (`s = 0`, `Z`
already a stabilizer) is the projector with **factor 1** — phases are `ℤ/2`-valued
(`selfInverse_zero_phase_even`, the "measurement result cannot be imaginary"
assert), survivors are untouched (`proj_zero_apply`), a projection that drops
nothing is the identity and hence norm-preserving (`proj_zero_eq_self`, which is
why case b applies no magnitude filter and normalizes only when the support
shrank), and the crate's `retain` predicate is exactly that surviving set
(`proj_zero_eq_caseB_retain`); fused batch gates — a batched sign is a site
parity, `count_ones() & 1` (`Batch.lean` `two_mul_natCast`,
`seqApply_eq_batchApply` on `Nodup` sites, `czSeq_phase` on pairwise-disjoint
pairs, with `czSeq_phase_needs_disjoint` showing the disjointness precondition is
**necessary**, and `isSitewise_*` pinning each gate's sign predicate against the
audited `Conjugation.lean` tables); multi-site conjugation
cross-phase — the `(−1)^{popcount(z_running ∧ x_new)}` fold of
`compute_decomposition_word` is the genuine Pauli product (`Word.lean`
`phaseExpN_eq_canon`, `Canon.toG_mul`, `Canon.foldl_eq_prod`).

Differential tests vs old `ppvm-tableau` / `ppvm-tableau-sum`: identical Clifford
circuits ⇒ identical measurement outcome distributions (seeded RNG) and identical
stabilizer state; generalized-tableau expectation values within tolerance.

Benchmarks: Clifford throughput, measurement, generalized-tableau Trotter step.

## Phase 5 — `ppvm-sym-2` (exact rings, closing the L4 loop)

Provide an exact/symbolic coefficient ring implementing `Coefficient` (crucially
*without* `Mul<f64>`), `ImaginaryUnit`, and `Conjugate`. This is the runtime
counterpart of the Lean `GaussianInt` instance: it makes `PauliSum` and
`multiply_into` run over an exact ring with no floats, validating that the
loosened bounds actually compile and compute.

Tests: exact twisted-product associativity (no tolerance); `conj(i) == −i`
(`Matrix.lean` `star_iU`); exact expectation values on small symbolic circuits vs
hand computation. Differential vs old `ppvm-sym` where semantics overlap.

## Phase 6 — `ColumnStore` (SoA backend), optional/perf

The batch contract in Phases 1/3 is defined so this slots in with no signature
changes: a structure-of-arrays `ColumnStore` implementing the same
`Support/Accumulate/Scale/Pair` with vectorized `scale`, prefix-sum `reduce`, and
coalesced `probe_batch`. Requires `K: Columnar` (already impl'd in Phase 2).

Gate on: correctness via the *same* Lean-oracle and differential suites (a backend
swap must be observationally identical), then SIMD/prefetch benchmarks showing the
intended speedup over the `HashMap` backend on large support.

## Phase 7 — downstream adapters and cutover

1. Add adapters so `ppvm-stim`, `ppvm-vihaco`, `ppvm-cli`, `ppvm-python-native`,
   `ppvm-tui` can target `-2` types behind a feature flag; run their existing
   test suites against both.
2. When every `-2` crate is green + at-parity and downstream passes on `-2`:
   flip defaults, then **rename** `ppvm-xxx-2 → ppvm-xxx`, deleting the old
   crate in the same commit (git rename keeps history). Update paths in
   `Cargo.toml`. Drop `ppvm-tableau-sum` (absorbed) and `ppvm-conformance-2`.
3. Final workspace build + full test + benchmark regression check.

## Testing strategy (cross-cutting)

- **Lean-oracle property tests** live *in each `-2` crate* (they need no old code)
  and encode the machine-checked identities. Prefer exhaustive `decide`-style
  checks for the finite single-qubit cases (matching the Lean `by decide` proofs)
  and `proptest`/randomized checks for n-qubit / large-support cases.
- **Differential tests** live only in `ppvm-conformance-2`, driven by seeded
  random circuits replayed on old and new, asserting equal observable behavior
  (support, coefficients within tolerance, measurement distributions).
- Keep a short table in `ppvm-conformance-2`'s README mapping each Lean theorem →
  the Rust test that enforces it, so drift is visible.

## Benchmark strategy (cross-cutting)

- Mirror the old crates' `benches/` layout. Each comparative bench runs the
  identical workload on old and new and prints the ratio.
- Parity target: new ≥ old on the hot paths (`mul`, Clifford, `key_hash`,
  accumulate/reduce, overlap, measurement). Regressions must be justified by a
  correctness requirement the design already accepts (e.g. lazy caching costing
  `Copy`) and signed off explicitly, not silently.
- Reuse `benchmarks/plot_branch_coalesce.py` and `plot_tfim_sweep.py` for the
  end-to-end sweeps.

## Risks and open questions

- **Lazy hashing vs `Copy`.** The design drops `Copy` for a correct
  interior-mutable cache. Confirm in Phase 2 that the throughput cost on `mul`/hash
  is acceptable; if not, revisit the cache representation (the contract only fixes
  the `key_hash()` *value*, not the mechanism).
- **`Multiply` and commutativity.** The Lean `tmul_assoc` needs `[CommRing C]`;
  the `Coefficient` trait does not encode commutativity. Decide whether L4/`Multiply`
  should carry a commutativity marker bound or keep it as a documented obligation.
- **`ColumnStore` plane granularity/alignment** (design open question 2) — defer
  to Phase 6 benchmarks.
- **`probe_batch` return shape** (design open question 4) — coefficients vs opaque
  slot handles; decide during Phase 3 from the overlap/expectation call sites.
- **Batch reuse** (design open questions 1/3) — transient vs driver-owned
  `TermBatch`, and group size; benchmark decisions in Phase 3/6.

## Automated refactor workflow

The per-crate work is driven by a deterministic orchestrator (a `Workflow`
script) that loops four agent roles until convergence. Decisions locked with the
maintainer:

- **Orchestrator owns `docs/log.md`.** Agents *return structured findings*; the
  script writes the gap tables. Race-free under parallel agents; resumable.
- **Autonomous within a crate, human gate at boundaries.** The loop runs
  unattended until a component converges; the maintainer signs off at each
  crate boundary (especially perf) and at cutover.
- **Perf = hard gate + allowlist.** Any regression blocks the crate unless it is
  on the human-approved allowlist in `docs/log.md` with written justification.
- **Impl-first proofs.** Implement against the *existing* Lean spec; the proof
  agent only adds Lean for genuinely new algebraic invariants the review agent
  nominates, and updates the design-doc citations (keeping design ↔ lean ↔ rust
  a triangle).

### Granularity

Pipeline **per component** within a crate (`PauliWord` → `LossyPauliWord` →
`Phased`), not per whole crate — tight agent context, fine checkpoints.

### The cycle (one component)

```
iter 1: implement   | iter 2+: fix open gaps in log.md
   │
   ├─ [script gate] cargo build ─ fail ─► back to impl (with compiler output)
   │
   ├─ review ∥ test+bench           (parallel; both read the impl, RETURN findings)
   │     review  → consistency (design↔lean↔rust) + nominate missing-proof invariants
   │     test    → differential vs old crate where a twin exists, else Lean-oracle
   │               + property tests; RETURN correctness/perf/missing-test gaps
   │
   ├─ [orchestrator] aggregate findings → write log.md; route gaps by type
   │
   ├─ prove          (only if missing-proof gaps) → add Lean + update doc citations
   │                  → [script gate] lake build
   │
   └─ converged?  (0 new gaps AND cargo build + tests + lake build green)
         no  → loop (max 4 iters, else escalate to human)
         yes → git commit checkpoint; distill log.md survivors into design doc
```

Gap routing (written to `docs/log.md` by the orchestrator):
`correctness`/`impl-friction` → impl or design; `perf-drift` → impl or human
sign-off; `missing-proof` → proof; `missing-test` → test.

### Agent tool profiles

| role | reads | writes | runs | parallel? |
| --- | --- | --- | --- | --- |
| impl | design, lean, old crate | new crate src | `cargo build` | serial |
| review | design, lean, new src, log | (returns findings) | — | ∥ with test |
| test+bench | new+old crates, log | test/bench code | `cargo test`, `cargo bench` | ∥ with review |
| proof | log gaps, lean, new src | lean, design citations | `lake build` | serial (after aggregate) |

### Deterministic gates live in the script

`cargo build`, `cargo test` / differential, and `lake build` are authoritative
pass/fail run by the orchestrator; agents interpret their output to produce
fixes, but the branch decision is not an agent's judgment call.

### What deserves a Lean proof

Algebraic / semantic invariants only — products, cocycles, isometries, module
laws, error bounds, channel eigenvalues, relabel bijections. **Not** systems
mechanics (allocation, SIMD/plane layout, I/O, caching *mechanism*). The review
agent nominates; the proof agent discharges and cites.
