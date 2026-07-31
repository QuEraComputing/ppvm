export const meta = {
  name: 'traits-2-component',
  description:
    'Implement + review + test + (prove) one ppvm-*-2 component; loop until zero gaps + green gates. Returns structured findings; the caller writes docs/log.md and commits.',
  phases: [
    { title: 'Implement', detail: 'port from design+lean+old crate; build/clippy/fmt until green' },
    { title: 'Verify', detail: 'review (design↔lean↔rust consistency) ∥ test (reference-impl/property tests)' },
    { title: 'Prove', detail: 'add Lean for newly-nominated algebraic invariants; lake build; update doc citations' },
  ],
}

const REPO = '/Users/roger/Code/rust/ppvm/.worktrees/traits-2-impl'
const MAX_ITERS = 3

// ── Component under work this run (Phase 1). Generalize via `args` later. ──────
const C = {
  crate: 'ppvm-traits-2',
  crateDir: `${REPO}/crates/ppvm-traits-2`,
  // Trait-definition modules to produce, with their design content (plan §Phase 1).
  modules: [
    'coefficient.rs — `Coefficient` (ring + mul_sign/half/magnitude, NO Mul<f64>), `Angle<C>` (+ impl Angle<f64> for f64); impl Coefficient for f64 and num::Complex<f64>',
    'algebra.rs — `KeyProduct` (key_mul -> (Self, Phase)), `ImaginaryUnit` (imaginary_unit(); law i*i == -one()), `Conjugate` (conj); `Phase` enum (Z/4 ~ {1,i,-1,-i}); impls for Complex<f64> and f64',
    'word.rs — `Word` (Site, n_sites, get, weight, iter); leaf types `Pauli`, `LossySite<S>`, `FermionSite`, `FermionAction`; `PauliBits: Word<Site = Pauli>`',
    'pauli.rs — `SymplecticColumns`, `PhaseTrack`, `StabilizerFrame`; blanket `impl<T: SymplecticColumns + PhaseTrack> Clifford for T`',
    'gates.rs — `Clifford`, `RotationOne<C, A = C>`, `PauliError<C>`, `Measure`',
    'graded.rs — `Support`, `Accumulate`, `Scale`, `Pair` (overlap AND hermitian_overlap where Coeff: Conjugate), `Multiply`, `Retain`',
    'batch.rs — `Columnar`, `KeyColumn`, `KeyBatch`, `TermBatch`, `TermSink`, `TermProducer`',
    'hash.rs — `Indexable` (key_hash), `IdentityHasher`, `IdentityBuildHasher`',
  ],
  // No old *behavioral* twin for pure trait defs, so tests are Lean-oracle style
  // on the concrete leaf types, plus compile checks.
  hasOldTwin: false,
}

const IMPL_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['summary', 'filesWritten', 'buildPassed', 'clippyClean', 'fmtClean', 'frictions'],
  properties: {
    summary: { type: 'string' },
    filesWritten: { type: 'array', items: { type: 'string' } },
    buildPassed: { type: 'boolean' },
    clippyClean: { type: 'boolean' },
    fmtClean: { type: 'boolean' },
    frictions: {
      type: 'array',
      description: 'Design signatures that could not be implemented verbatim, and how they were resolved.',
      items: {
        type: 'object',
        additionalProperties: false,
        required: ['description', 'resolution', 'severity'],
        properties: {
          description: { type: 'string' },
          resolution: { type: 'string' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
        },
      },
    },
  },
}

const REVIEW_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['consistencyGaps', 'missingProofInvariants', 'verdict'],
  properties: {
    consistencyGaps: {
      type: 'array',
      items: {
        type: 'object',
        additionalProperties: false,
        required: ['type', 'severity', 'location', 'description', 'routedTo'],
        properties: {
          type: { type: 'string', enum: ['impl-friction', 'correctness', 'perf-drift', 'missing-test'] },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
          location: { type: 'string' },
          description: { type: 'string' },
          routedTo: { type: 'string', enum: ['impl', 'design', 'test', 'human'] },
        },
      },
    },
    missingProofInvariants: {
      type: 'array',
      description: 'Algebraic/semantic invariants a trait encodes that lack a Lean theorem and deserve one.',
      items: {
        type: 'object',
        additionalProperties: false,
        required: ['name', 'statement', 'whyNeeded'],
        properties: {
          name: { type: 'string' },
          statement: { type: 'string' },
          whyNeeded: { type: 'string' },
        },
      },
    },
    verdict: { type: 'string' },
  },
}

const TEST_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['testsAdded', 'testsPass', 'gaps', 'perfNote'],
  properties: {
    testsAdded: { type: 'array', items: { type: 'string' } },
    testsPass: { type: 'boolean' },
    gaps: {
      type: 'array',
      items: {
        type: 'object',
        additionalProperties: false,
        required: ['type', 'severity', 'description'],
        properties: {
          type: { type: 'string', enum: ['correctness', 'missing-test', 'impl-friction', 'perf-drift'] },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
          description: { type: 'string' },
        },
      },
    },
    perfNote: { type: 'string' },
  },
}

const PROVE_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['proofsAdded', 'designCitationsUpdated', 'lakeBuildPassed', 'newGapsRaised', 'notes'],
  properties: {
    proofsAdded: {
      type: 'array',
      items: {
        type: 'object',
        additionalProperties: false,
        required: ['file', 'theorem', 'statement'],
        properties: {
          file: { type: 'string' },
          theorem: { type: 'string' },
          statement: { type: 'string' },
        },
      },
    },
    designCitationsUpdated: { type: 'boolean' },
    lakeBuildPassed: { type: 'boolean' },
    newGapsRaised: { type: 'boolean' },
    notes: { type: 'string' },
  },
}

const COMMON = `
Repo root: ${REPO}
Run cargo from the repo root; run \`lake build PPVM\` from ${REPO}/lean.
Design doc: ${REPO}/docs/design/traits-2-configuration-and-hashing.md (the authoritative trait signatures).
Implementation plan: ${REPO}/docs/design/traits-2-implementation-plan.md (§Phase 1 lists the modules).
Lean spec (machine-checked semantics): ${REPO}/lean/PPVM/**.
Old reference crate (algorithm/detail): ${REPO}/crates/ppvm-traits.
Every .rs file MUST start with the two SPDX header lines used across the repo
(see any file under crates/ppvm-pauli-word/src). Do NOT edit docs/log.md — return
your findings as structured output; the orchestrator records them.`

const IMPL_PROMPT = `You are the IMPLEMENTATION agent for the ppvm-*-2 refactor, Phase 1: the
${C.crate} trait DEFINITIONS.
${COMMON}

Task: create these modules under ${C.crateDir}/src (and wire them in lib.rs),
porting the trait signatures VERBATIM from the design doc where they compile, and
matching the Lean-validated semantics:
${C.modules.map((m, i) => `  ${i + 1}. ${m}`).join('\n')}

Rules:
- This crate is trait DEFINITIONS + small leaf types ONLY. Do NOT implement the
  concrete PauliWord/Tableau (those are Phase 2/4). Provide the blanket
  \`impl<T: SymplecticColumns + PhaseTrack> Clifford for T\`.
- Add needed deps to Cargo.toml (e.g. \`num\`). Keep it minimal; no unused deps
  (\`cargo machete\` runs in the pre-commit gate).
- On each trait/method that encodes a machine-checked invariant, add a doc
  comment citing BOTH the design section AND the relevant Lean theorem
  (e.g. ImaginaryUnit law -> Matrix.lean \`iU_sq\`; KeyProduct twist -> Twisted.lean
  \`tmul_assoc\`; Conjugate -> Matrix.lean \`star_iU\`).
- Where a design signature CANNOT be implemented verbatim (e.g. the graded
  \`Accumulate::accumulate_batch\` takes \`TermBatch<Self::Key, Self::Coeff>\` but the
  batch types are bounded on \`Columnar\` while L0 \`Support::Key\` is only
  \`Eq + Clone\`): make the MINIMAL adjustment that compiles AND preserves the
  design intent, and record it in \`frictions\` (description + your resolution +
  severity). Do not silently diverge.
- Preserve the design's allocation shape (e.g. \`Word::iter\` returns
  \`impl Iterator\`; batch I/O is columnar; no \`&mut (K,C)\`/\`&mut [C]\` slot access).

Gates you MUST pass before returning (run them; iterate until green):
  cargo build -p ${C.crate} --all-targets
  cargo clippy -p ${C.crate} --all-targets -- -D warnings
  cargo fmt -p ${C.crate} -- --check   (run \`cargo fmt -p ${C.crate}\` to fix)
Return the schema. buildPassed/clippyClean/fmtClean must reflect the ACTUAL final
command exit status.`

const fixPrompt = (codeGaps) => `You are the IMPLEMENTATION agent continuing Phase 1 on ${C.crate}. A prior
review/test pass found gaps that route to code. Fix ONLY these; do not rewrite
working modules.
${COMMON}

Gaps to resolve (JSON):
${JSON.stringify(codeGaps, null, 2)}

Re-run and pass the same gates (build --all-targets, clippy -D warnings,
fmt --check). Record any new frictions. Return the schema.`

const REVIEW_PROMPT = `You are the REVIEW agent for Phase 1 of the ppvm-*-2 refactor. READ-ONLY: do not
edit files; return findings.
${COMMON}

Evaluate the freshly-written ${C.crate}/src against the design doc and the Lean
spec, on three axes:
1. Consistency — does every trait/method match the design's signature and the
   Lean-validated semantics? Is anything from the design's Phase-1 trait list
   MISSING or renamed? Did the impl agent's friction resolutions drift from the
   design intent (if so, route to 'design' or 'impl')?
2. Correctness (abstract algebra) — are the encoded contracts faithful (e.g.
   ImaginaryUnit law i·i=−1; Conjugate involution + conj(i)=−i; KeyProduct emits
   Phase; Pair split bilinear vs sesquilinear; graded layers L0–L4 bounds; the
   Clifford blanket over SymplecticColumns+PhaseTrack)?
3. Missing proofs — nominate any algebraic/semantic invariant this crate now
   encodes that is NOT yet covered by a Lean theorem and genuinely deserves one
   (algebra only — NOT systems mechanics like hashing mechanism or allocation).
   For Phase 1 most algebra is already covered; nominate sparingly and justify.

Cite exact file:line or trait/theorem names. Route each consistency gap
(impl/design/test/human). Return the schema.`

const TEST_PROMPT = `You are the TEST agent for Phase 1 of the ppvm-*-2 refactor.
${COMMON}

${C.crate} is pure trait defs + leaf types, so there is NO old behavioral twin to
diff. Instead add Lean-oracle-style unit tests (in ${C.crateDir}, e.g. a
\`#[cfg(test)]\` module or tests/) that pin the concrete leaf types and provided
impls:
- \`ImaginaryUnit for Complex<f64>\`: imaginary_unit()*imaginary_unit() == -one().
- \`Conjugate\`: conj(conj(x)) == x for Complex<f64>; conj is identity for f64;
  conj(i) == -i.
- \`Phase\` arithmetic behaves as ℤ/4 (compose/one/inverse) and matches the
  {1,i,-1,-i} interpretation.
- \`Angle<f64> for f64\`: sin_cos matches f64::sin_cos.
- \`IdentityHasher\`/\`IdentityBuildHasher\`: write_u64(n) then finish() == n
  (pass-through).
- A trivial reference \`Coefficient\` impl (e.g. over f64) compiles and satisfies
  the ring ops (smoke).
Also add a compile-only test that a stub type implementing SymplecticColumns +
PhaseTrack automatically gets \`Clifford\` (the blanket impl).

You ADD test files, so you MUST leave the crate green on ALL targets — run, in
order, and only set testsPass=true if every one passes:
  cargo fmt -p ${C.crate}                       (format YOUR new files)
  cargo fmt -p ${C.crate} -- --check            (must be clean)
  cargo clippy -p ${C.crate} --all-targets -- -D warnings
  cargo test -p ${C.crate}
(The fmt step matters: a prior run left an unformatted test file because only the
impl agent formatted.) Pure trait defs have no hot path, so do NOT add benchmarks
this phase — set perfNote to say so. Report tests added, pass/fail, and any gaps
(e.g. an invariant you could not test, or a correctness surprise). Return the
schema.`

const provePrompt = (noms) => `You are the PROOF agent for the ppvm-*-2 refactor. The review agent nominated
algebraic invariants the new ${C.crate} code encodes that may lack a Lean theorem.
${COMMON}

Nominated invariants (JSON):
${JSON.stringify(noms, null, 2)}

For each: check whether an existing theorem in ${REPO}/lean/PPVM already covers it
(many do). If genuinely uncovered AND it is a real algebraic/semantic invariant
(not systems mechanics), add a minimal Lean theorem in the appropriate file, and
update the design-doc citation so design↔lean↔rust stay a triangle. Skip
(with a note) anything already covered or not proof-worthy.

You MUST run \`lake build PPVM\` from ${REPO}/lean and it must pass. Return the
schema; lakeBuildPassed must reflect the actual exit status.`

// ── Convergence loop ──────────────────────────────────────────────────────────
const history = []
let converged = false
let pendingCodeGaps = null
let iter = 0

while (iter < MAX_ITERS && !converged) {
  iter++
  phase('Implement')
  const impl = await agent(iter === 1 ? IMPL_PROMPT : fixPrompt(pendingCodeGaps), {
    label: `impl#${iter}`,
    phase: 'Implement',
    schema: IMPL_SCHEMA,
  })
  if (!impl || !impl.buildPassed) {
    history.push({ iter, impl })
    return { component: C.crate, status: 'escalate', reason: 'build/gate failed', iter, history }
  }

  phase('Verify')
  const [review, test] = await parallel([
    () => agent(REVIEW_PROMPT, { label: `review#${iter}`, phase: 'Verify', effort: 'high', schema: REVIEW_SCHEMA }),
    () => agent(TEST_PROMPT, { label: `test#${iter}`, phase: 'Verify', schema: TEST_SCHEMA }),
  ])

  // Only high/medium correctness or impl-friction blocks convergence. Low-sev
  // ergonomics notes, and anything the reviewer routes to 'design'/'human'/'test'
  // or explicitly defers to a later phase, do NOT block (they are recorded, then
  // acted on by the orchestrator/human). This avoids the Phase-1 artifact where a
  // low-sev ergonomics note looped the crate to MAX_ITERS.
  const blocks = (g) =>
    (g.type === 'correctness' || g.type === 'impl-friction') &&
    (g.severity === 'high' || g.severity === 'medium') &&
    !/defer|phase\s*\d/i.test(g.description || '')
  const codeGaps = [
    ...((review && review.consistencyGaps) || []).filter((g) => g.routedTo === 'impl' && blocks(g)),
    ...((test && test.gaps) || []).filter(blocks),
  ]
  const proofNoms = (review && review.missingProofInvariants) || []

  // Discharge nominated proofs. Once the prove agent adds them and `lake build`
  // is green, they are RESOLVED for convergence — a later review must not
  // re-nominate the same invariant to keep the loop from spinning.
  let prove = null
  if (proofNoms.length) {
    phase('Prove')
    prove = await agent(provePrompt(proofNoms), { label: `prove#${iter}`, phase: 'Prove', effort: 'high', schema: PROVE_SCHEMA })
  }

  history.push({ iter, impl, review, test, prove })

  const testsOk = !test || test.testsPass !== false
  // Proof step is satisfied when there were no nominations, or the prover
  // discharged them with a green `lake build` (its own raised gaps are folded
  // into the next iteration's code fix only if they block).
  const proofOk = !proofNoms.length || (prove && prove.lakeBuildPassed)
  converged = codeGaps.length === 0 && testsOk && proofOk
  pendingCodeGaps = codeGaps
}

return {
  component: C.crate,
  status: converged ? 'converged' : 'escalate',
  iters: iter,
  history,
}
