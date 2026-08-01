export const meta = {
  name: 'traits-2-component',
  description:
    'Implement + review + test + (prove) one ppvm-*-2 component; loop until zero gaps + green gates. Behavioral components also differential-test vs the old crate and gate on perf. Component spec comes from `args` (Phase-1 traits default). Returns structured findings; the caller writes docs/log.md and commits.',
  phases: [
    { title: 'Implement', detail: 'port from design+lean+old crate; build/clippy/fmt until green' },
    { title: 'Verify', detail: 'review (design↔lean↔rust) ∥ test (Lean-oracle; +differential+bench if an old twin exists)' },
    { title: 'Prove', detail: 'add Lean for newly-nominated algebraic invariants; lake build; update doc citations' },
  ],
}

const REPO = '/Users/roger/Code/rust/ppvm/.worktrees/traits-2-impl'
const CONF = 'ppvm-conformance-2'
const MAX_ITERS = 3

// ── Phase-1 default component (used when no `args` spec is passed). ────────────
const DEFAULT = {
  crate: 'ppvm-traits-2',
  crateDir: `${REPO}/crates/ppvm-traits-2`,
  title: 'trait definitions',
  hasOldTwin: false,
  implGoal:
    'create the trait-definition modules (coefficient/algebra/word/pauli/gates/graded/batch/hash) porting the design signatures verbatim where they compile; trait defs + leaf types only, no concrete PauliWord/Tableau; provide the blanket `impl<T: SymplecticColumns + PhaseTrack> Clifford for T`.',
  leanOracles:
    'ImaginaryUnit i·i==−one; Conjugate involution + conj(i)=−i; Phase as ℤ/4; Angle<f64>; IdentityHasher pass-through; the Clifford blanket over a stub SymplecticColumns+PhaseTrack.',
}

// Component under work THIS run: `args` if provided, else the Phase-1 default.
const C = args && args.crate ? args : DEFAULT

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
    differentialPass: { type: 'boolean', description: 'new matches old on the diff surface (behavioral components)' },
    perfRatios: {
      type: 'array',
      description: 'new/old wall-clock ratio per benchmarked target (behavioral components). MUST be a FAIR, engine-to-engine comparison: identical algebraic config on both sides, and a ratio confirmed stable across ≥2 runs.',
      items: {
        type: 'object',
        additionalProperties: false,
        required: ['target', 'ratio', 'config', 'stable'],
        properties: {
          target: { type: 'string' },
          ratio: { type: 'number', description: 'new/old, median-vs-median' },
          config: {
            type: 'string',
            description: 'The algebraic config held IDENTICAL on both sides — storage width (e.g. [u8;8]), coefficient type, and comparable hasher. State it explicitly; a mismatch here invalidates the ratio (e.g. u64-new vs [u8;8]-old folds BitArray codegen into the ratio).',
          },
          stable: {
            type: 'boolean',
            description: 'true only if this is a SAME-BUILD new/old ratio (both benched in one binary — cross-build absolutes are unreliable from layout/Mytkowicz bias) that held across ≥2 runs with neither side on a wide CI. If a baseline is noisy or the ratio comes from comparing separate builds, set false and do NOT read a regression off it — rerun / widen measurement time / use an interleaved same-build harness first.',
          },
          attribution: {
            type: 'string',
            description: 'For any ratio flagged as drift: the VERIFIED cause from a controlled A/B (hold one variable) or a profile — NOT a guess. If not isolated, write "unattributed" rather than assert a mechanism.',
          },
          note: { type: 'string' },
        },
      },
    },
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
Design (authoritative signatures): ${REPO}/docs/design/traits-2-configuration-and-hashing.md
Concrete layouts: ${REPO}/docs/design/word-data-structures.md
Plan: ${REPO}/docs/design/traits-2-implementation-plan.md
Lean spec (machine-checked semantics): ${REPO}/lean/PPVM/**
Trait foundation (already built): ${REPO}/crates/ppvm-traits-2
${C.oldRef ? `Old reference crate (port its ALGORITHM to minimize perf drift): ${REPO}/${C.oldRef}` : `Old reference crate (algorithm/detail): ${REPO}/crates/ppvm-traits`}
Every .rs file MUST start with the two SPDX header lines used across the repo.
Do NOT edit docs/log.md — return findings as structured output; the orchestrator records them.`

const IMPL_PROMPT = `You are the IMPLEMENTATION agent for the ppvm-*-2 refactor. Component: ${C.crate}
(${C.title}).
${COMMON}

Task: ${C.implGoal}
${C.modules ? `Modules:\n${C.modules.map((m, i) => `  ${i + 1}. ${m}`).join('\n')}` : ''}

Rules:
${
  C.hasOldTwin
    ? `- Base the algorithm on the OLD crate (${C.oldRef}) to MINIMIZE performance drift; keep its packed layout and hot-path structure. Follow word-data-structures.md for the concrete layout. Backing fields (packed planes, hash cache) are PRIVATE; expose behavior only through the ppvm-traits-2 traits it implements.
- Use a LAZY structural hash (OnceLock<u64>) per the design (Copy is intentionally dropped); the digest value returned by Indexable::key_hash() must be the finalized, avalanche-quality digest.
- Add exactly the deps you use to Cargo.toml (add ppvm-traits-2; \`cargo machete\` gates unused deps).`
    : `- Trait definitions + small leaf types only; no heavy logic.`
}
- Cite the design section AND the relevant Lean theorem in the doc comment of any item that encodes a machine-checked invariant.
- Where a design signature can't be implemented verbatim, make the MINIMAL adjustment that compiles AND preserves intent, and record it in \`frictions\` (never silently diverge).

Gates you MUST pass before returning (iterate until green):
  cargo build -p ${C.crate} --all-targets
  cargo clippy -p ${C.crate} --all-targets -- -D warnings
  cargo fmt -p ${C.crate} -- --check   (run \`cargo fmt -p ${C.crate}\` to fix)
Return the schema; the boolean flags must reflect the ACTUAL final exit status.`

const fixPrompt = (codeGaps) => `You are the IMPLEMENTATION agent continuing on ${C.crate}. A prior review/test
pass found gaps routed to code. Fix ONLY these; do not rewrite working code.
${COMMON}

Gaps to resolve (JSON):
${JSON.stringify(codeGaps, null, 2)}

Re-run and pass the same gates (build --all-targets, clippy -D warnings,
fmt --check). If a gap is a perf regression, optimize toward parity with the old
crate; if you believe it is an unavoidable design-accepted trade-off, say so
precisely in \`frictions\` (the human owns the perf allowlist). Return the schema.`

const REVIEW_PROMPT = `You are the REVIEW agent (READ-ONLY: do not edit files; return findings) for
component ${C.crate}.
${COMMON}

Evaluate the freshly-written ${C.crate} against the design and the Lean spec:
1. Consistency — every trait impl matches the design contract and Lean-validated
   semantics; nothing from the design is missing; friction resolutions do not
   drift (route drift to 'design' or 'impl').
2. Correctness (abstract algebra) — the encoded contracts are faithful${
  C.hasOldTwin
    ? ' (Pauli phase product = phaseExp/pauliMat_mul; Clifford conjugation bit+phase rules incl. CNOT/CZ; symplectic-form preservation; hashing contract)'
    : ''
}.
3. Missing proofs — nominate ONLY genuinely-uncovered algebraic/semantic
   invariants that deserve a Lean theorem (not systems mechanics like hashing
   mechanism, allocation, or SIMD layout). Nominate sparingly and justify.

Cite exact file:line / trait / theorem names. Route each consistency gap. Return
the schema.`

const TEST_PROMPT = C.hasOldTwin
  ? `You are the TEST agent for component ${C.crate}. You have TWO jobs: differential
correctness vs the old crate, and a performance gate.
${COMMON}

Reuse the seeded generators in ${CONF} (\`seeded_rng\`, \`random_pauli_string\`,
\`random_circuit\`, \`GateOp\`). Add ${C.crate} as a dependency of ${CONF} (edit its
Cargo.toml) so tests/benches can see both crates.

1. DIFFERENTIAL correctness — in ${CONF}/tests/, for seeded random inputs, assert
   the NEW ${C.crate} matches the OLD (${C.oldRef}) on: ${C.diffSurface}.
   Compare OBSERVABLE algebra, NOT raw hash digests (the finalization fold differs
   by design). For hashing, test the CONTRACT instead: Hash writes exactly
   key_hash(); structurally-equal keys => equal digest; and an
   avalanche/low-collision distribution property test. Set differentialPass.
2. LEAN-ORACLE property tests — reproduce ${C.leanOracles} as Rust tests
   (exhaustive for finite single-qubit cases; randomized for n-qubit).
3. PERFORMANCE GATE — add Criterion benches in ${CONF}/benches/ comparing NEW vs
   OLD on: ${C.benchTargets}. Report the new/old ratio per target in \`perfRatios\`
   + \`perfNote\`. A perf ratio is only meaningful if the comparison is FAIR and
   the measurement is STABLE, and any drift you flag must be ATTRIBUTED, not
   guessed. Three MANDATORY rules — a ratio that violates any of them is invalid
   and must NOT be reported as a regression:

   (a) FAIR CONFIG (apples-to-apples). Hold the algebraic config IDENTICAL on both
       sides: same storage width (e.g. both \`[u8; 8]\` — do NOT bench the shipped
       \`u64\` default new-sum against an \`[u8; 8]\` old-sum; \`BitArray<u64>\` vs
       \`BitArray<[u8;8]>\` codegen differs a few percent and would fold into the
       ratio), same coefficient type, comparable hasher. If the shipped default
       differs from the old crate's, build a bench-LOCAL storage-matched new type
       (correctness is storage-independent and is covered by the differential
       suite on the shipped default). State the matched config in \`config\` for
       every target. A config mismatch INVALIDATES the number.

   (b) STABLE MEASUREMENT (no cherry-picking). Use enough measurement time for
       tight CIs, and confirm each ratio across ≥2 runs. If either side has a wide
       CI / noisy baseline, set \`stable:false\` and DO NOT read a regression off
       it — rerun or widen \`--measurement-time\` until it settles. Never quote a
       ratio (in EITHER direction — a flattering "parity" is as wrong as a scary
       "regression") from a single run whose baseline swings.
       SAME-BUILD RATIOS ONLY. The gate metric is the new/old ratio measured
       *within a single binary* (both benched under one code layout + thermal
       state). Treat ABSOLUTE cross-build numbers as unreliable: relinking for a
       code change relayouts the whole bench binary, so even an UNTOUCHED baseline
       swings from function alignment / i-cache effects (the Mytkowicz "producing
       wrong data" layout bias — here old/rx moved 4.5↔5.9µs+ with zero code
       change). So NEVER conclude from "old/rx was 4.8µs last build, 5.9µs this
       build → regression"; only same-build new-vs-old (where the layout bias
       cancels in the ratio) is sound. To A/B a fix, prefer an interleaved
       harness that measures both variants in ONE process (min-of-N), not two
       separate criterion builds.

   (c) VERIFIED ATTRIBUTION (no hand-waving). Only raise a 'perf-drift' gap for a
       ratio that is fair (a) AND stable (b). When you do, the cause in
       \`attribution\` must come from a CONTROLLED experiment (A/B holding ONE
       variable — storage, hasher, alloc, policy) or a profile, NOT a plausible
       story. If you cannot isolate it, write "unattributed" — do not assert a
       mechanism. Note whether it is a design-accepted trade-off (e.g. lazy hash
       cache costing the old \`Copy\` word on fresh-key paths). The HUMAN owns the
       allowlist, so do not silently accept it — but also do not manufacture a
       regression from an unfair or noisy bench.

Leave BOTH crates green (run, and set testsPass only if all pass):
  cargo fmt -p ${C.crate} -p ${CONF}   then   cargo fmt -p ${C.crate} -p ${CONF} -- --check
  cargo clippy -p ${C.crate} -p ${CONF} --all-targets -- -D warnings
  cargo test -p ${C.crate} -p ${CONF}
Return the schema.`
  : `You are the TEST agent for component ${C.crate} (pure trait defs — no old
behavioral twin to diff).
${COMMON}

Add Lean-oracle unit tests pinning the concrete leaf types / provided impls:
${C.leanOracles}
Leave the crate green (run, set testsPass only if all pass):
  cargo fmt -p ${C.crate}   then   cargo fmt -p ${C.crate} -- --check
  cargo clippy -p ${C.crate} --all-targets -- -D warnings
  cargo test -p ${C.crate}
No hot path here, so add NO benchmarks — set perfNote to say so. Report tests
added, pass/fail, and any gaps. Return the schema.`

const provePrompt = (noms) => `You are the PROOF agent. The review nominated algebraic invariants the new
${C.crate} code encodes that may lack a Lean theorem.
${COMMON}

Nominated invariants (JSON):
${JSON.stringify(noms, null, 2)}

For each: check whether an existing theorem in ${REPO}/lean/PPVM already covers it.
If genuinely uncovered AND a real algebraic/semantic invariant (not systems
mechanics), add a minimal Lean theorem in the right file and update the
design-doc citation (keep design↔lean↔rust a triangle). Skip (with a note)
anything already covered or not proof-worthy. You MUST run \`lake build PPVM\` from
${REPO}/lean and it must pass. Return the schema.`

// ── Convergence loop ──────────────────────────────────────────────────────────
// Perf-drift is a HARD gate: any regression over the bench threshold must reach
// the human (who owns the allowlist), so it blocks regardless of the agent's
// severity label — unless the agent marks it an already-accepted trade-off.
// Correctness/impl-friction block at high/medium; low-sev ergonomics and things
// routed to design/human/test or deferred to a later phase are recorded but do
// NOT spin the loop.
//
// VALIDITY GUARD (orchestrator, before escalating a perf-drift gap to the human):
// a perf-drift gap is only real if its `perfRatios` entry is FAIR and STABLE and
// ATTRIBUTED — i.e. `config` states an identical-on-both-sides algebraic config
// (matched storage width / coeff / hasher), `stable` is true (ratio held over ≥2
// runs, no wide-CI baseline), and `attribution` names a controlled-A/B or profiled
// cause (not a guess, not "unattributed"). If any of those fails, the number is a
// benchmark artifact, NOT a regression: bounce it back to the test agent to
// re-measure (fair config, more measurement time, isolate the cause) rather than
// escalate. This is the guard the ps2.rot.perf saga needed — the "1.15×" was
// u64-vs-[u8;8] apples-to-oranges, and the later "parity 1.01×" was a cherry-picked
// noisy old baseline; both would have been caught by requiring fair+stable+attributed.
const blocks = (g) =>
  g.type === 'perf-drift'
    ? // ALWAYS block: any perf regression the test agent raises as a gap needs
      // human sign-off (allowlist or fix). Do NOT try to self-clear from the
      // description — a negated mention ("NOT the design-accepted trade-off")
      // must not match. Genuinely design-accepted trade-offs (e.g. the lazy-hash
      // cold path) are explained in perfNote, never raised as a gap.
      true
    : (g.type === 'correctness' || g.type === 'impl-friction') &&
      (g.severity === 'high' || g.severity === 'medium') &&
      !/defer|phase\s*\d/i.test(g.description || '')

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

  const codeGaps = [
    ...((review && review.consistencyGaps) || []).filter((g) => g.routedTo === 'impl' && blocks(g)),
    ...((test && test.gaps) || []).filter(blocks),
  ]
  const proofNoms = (review && review.missingProofInvariants) || []

  let prove = null
  if (proofNoms.length) {
    phase('Prove')
    prove = await agent(provePrompt(proofNoms), { label: `prove#${iter}`, phase: 'Prove', effort: 'high', schema: PROVE_SCHEMA })
  }

  history.push({ iter, impl, review, test, prove })

  const testsOk = !test || test.testsPass !== false
  const diffOk = !C.hasOldTwin || !test || test.differentialPass !== false
  const proofOk = !proofNoms.length || (prove && prove.lakeBuildPassed)
  converged = codeGaps.length === 0 && testsOk && diffOk && proofOk
  pendingCodeGaps = codeGaps
}

return {
  component: C.crate,
  status: converged ? 'converged' : 'escalate',
  iters: iter,
  history,
}
