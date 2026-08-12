# Lean ↔ `ppvm-*-2` verification gap ledger

Round 1 (audit), 2026-08-10. Branch `codex/traits-2-impl`.

## What this file is

A durable, resumable ledger of the gaps between the Lean formalization in
`lean/` and the Rust it is supposed to verify. It exists because the campaign is
iterative: audit rounds open rows, close rounds write Lean theorems and
`*_lean.rs` oracle tests that retire them, and the campaign stops when two
consecutive audit rounds open nothing new. Every row therefore carries a stable
ID, an explicit status, and enough detail — the exact Rust site, the exact Lean
symbol, the claim that is unverified, why the current artifacts do not establish
it, and the closure that would — that a close-round agent can act on a row
without re-deriving the audit.

**Target.** The `ppvm-*-2` crates only: `ppvm-traits-2`, `ppvm-pauli-word-2`,
`ppvm-phased-pauli-word-2`, `ppvm-lossy-pauli-word-2`, `ppvm-pauli-sum-2`,
`ppvm-sym-2`, `ppvm-tableau-2`, and the oracle/differential suite in
`ppvm-conformance-2`. The legacy twins (`ppvm-pauli-word`, `ppvm-pauli-sum`,
`ppvm-tableau`, `ppvm-tableau-sum`, `ppvm-sym`, `ppvm-traits`) are read only to
adjudicate provenance; they are never the verification target.

**Verified means two legs.** A claim counts as verified only when both hold:

1. a Lean theorem in `lean/PPVM/**` states the claim over a model that is
   faithful to the shipped code (not to a legacy predecessor, not to a
   convenient abstraction that drops the load-bearing structure); and
2. a `*_lean.rs` oracle test in `crates/ppvm-conformance-2/tests/` pins the
   **real Rust type** to that theorem by name, and would **fail** under a
   plausible mutation of the code the theorem is about.

Leg 2's mutation criterion is the sharp one, and most `bridge` rows below fail
exactly there: the test runs, is green, cites a theorem, and is invariant under
the defect the theorem rules out.

**Agreement with legacy is not evidence of correctness.** The `*_diff.rs` tests
establish that a `-2` implementation reproduces its legacy predecessor. That is
evidence of a faithful port and of nothing else: where legacy is wrong, the diff
test is wrong with it, and where legacy and `-2` share a kernel byte for byte, a
diff test cannot detect a kernel defect by construction. This is why
`provenance` is a gap class of its own, and why several rows below are opened
against code that has extensive, green differential coverage.

**Gap classes.**

| class | meaning |
| --- | --- |
| `coverage` | The Rust does something no Lean definition or theorem models at all. |
| `fidelity` | A Lean model exists but does not match the shipped code (drifted citation, dropped side condition, different convention, wrong direction). |
| `strength` | The theorem is true but weaker than the claim credited to it (an identity restated, a hypothesis assumed, the interesting step inserted by hand). |
| `bridge` | Lean and Rust are both fine, but no `*_lean.rs` oracle pins one to the other, or the oracle that claims to is mutation-insensitive. |
| `provenance` | The only thing holding a behaviour in place is agreement with legacy. |

**Tiers.** Tier A is the mathematical kernel (Pauli algebra, conjugation,
measurement, channels, truncation bounds, symbolic ring laws). Tier B is
representation and engine machinery whose correctness is still a mathematical
claim (packed-bit layouts, store refinements, digests, batch contracts). Tier C
— `ppvm-cli`, `ppvm-tui`, `stim-parser`, `ppvm-stim`, `ppvm-python-native`,
`ppvm-vihaco`, `vihaco-circuit-isa` — is out of scope and was not read; see
[Out of scope](#out-of-scope).

**Severity** is the consequence if the claim is false, not the effort to close:
`high` = a wrong number reaches a user (an observable, a measurement outcome, a
sampled distribution) or a silently non-canonical state enters a hash table;
`medium` = a real defect class is undetectable by the current suite;
`low` = a documented convention or a citation is wrong, or the arithmetic is
uncovered but its consequence is contained.

**Statuses.** `open` (raised, not yet adjudicated or closed), `closed` (Lean +
oracle landed, with the commit in `Evidence`), `refuted` (a later round showed the
claim was already verified or the gap was miscategorized — see the
[appendix](#appendix-refuted-candidates), and do not re-litigate those), plus the
three adjudication statuses introduced in round 2:

| status | meaning |
| --- | --- |
| `adjudicated-defect` | Re-derived from the mathematics; a live defect confirmed. A fix is proposed in [Adjudications](#adjudications-round-2) with exact file and line, and it awaits the user's sign-off because it changes observable behaviour. **The row must not be closed by proving the current behaviour correct.** |
| `adjudicated-spec` | Re-derived; the code is right. The Lean statement and/or a docstring is wrong or overstated and must be corrected as part of closing the row. Any Lean theorem written for this row must state the *adjudicated* claim, not the one the row originally asserted. |
| `adjudicated-undecided` | Two independent derivations disagreed, or both were inconclusive. Must not be proved either way until the contest is resolved. No row carries this status after round 2. |

An adjudicated row is still open work: the status records *what the answer is*, not
that the Lean/oracle legs exist. Where an adjudication showed a row's
`Proposed closure` would have proved the wrong statement, that closure has been
edited in place and points at the adjudication subsection.

## Status dashboard

| status | count |
| --- | ---: |
| open | 56 |
| adjudicated-defect | 9 |
| adjudicated-spec | 18 |
| adjudicated-undecided | 0 |
| closed | 0 |
| refuted (appendix, not counted as open) | 16 |

All 83 round-1 rows are still live work; the 27 adjudicated ones have an
established answer (see [Adjudications](#adjudications-round-2)) but no Lean or
oracle leg yet. The `open`/`adjudicated` split is carried through every breakdown
below so that neither number is lost.

| class | open | adjudicated |
| --- | ---: | ---: |
| bridge | 19 | 7 |
| coverage | 15 | 10 |
| strength | 11 | 6 |
| fidelity | 9 | 3 |
| provenance | 2 | 1 |

| severity | open | adjudicated |
| --- | ---: | ---: |
| high | 15 | 12 |
| medium | 32 | 13 |
| low | 9 | 2 |

| tier | open | adjudicated |
| --- | ---: | ---: |
| A | 34 | 21 |
| B | 22 | 6 |

| sector | open | high | adjudicated | notes |
| --- | ---: | ---: | ---: | --- |
| clifford-conjugation | 5 | 0 | 0 | Generators fully grounded; extensions and the phaseless words unbridged. |
| graded-algebra-containers | 2 | 0 | 1 | Container L4 impls untested; `reduce`/`len` unmodelled. |
| hashing-digests | 3 | 0 | 1 | No Lean module exists for this sector at all. |
| lossy-word | 1 | 0 | 3 | Loss is an external parameter in Lean; canonicality unchecked after propagation. |
| measurement-branching | 1 | 1 | 6 | Every ℤ/4 sign in measurement is unpinned; `Projection.lean` has no bridge. |
| multiply-rotation | 5 | 1 | 4 | Product excellent; rotation direction/sign essentially unverified. |
| noise-observables | 9 | 5 | 3 | Eigenvalue never derived from a channel; three backends implement two channels. |
| products-and-channels (skeptic) | 0 | 0 | 1 | `asymmetric_loss_channel` untested at `p0 ≠ p1`. |
| sum-engine-stores | 6 | 2 | 0 | No abstraction function from any store to `K →₀ C`. |
| symbolic-coefficients | 4 | 1 | 2 | Ring is well modelled; the hom law is never pinned to a map-backed `Term`. |
| tableau-and-symbolic (skeptic) | 0 | 0 | 1 | The tableau's private g-rule is pinned to no phase model. |
| tableau-core | 4 | 1 | 4 | Bits verified; the sign column has no invariant and no ext-gate oracle. |
| truncation-policy-loss | 9 | 4 | 1 | Bounds proved for arbitrary drop sets, never tied to a policy parameter. |
| word-algebra | 6 | 0 | 0 | Phase kernel is the best-verified code in the repo; packed storage invisible. |
| word-and-clifford (skeptic) | 1 | 0 | 0 | `PauliSum`'s fourth copy of the extension signs is unbridged. |

Fifteen sector labels appear above; twelve are the audit sectors surveyed
end to end (see [Sector coverage](#sector-coverage)) and three suffixed
`(skeptic)` are cross-cutting rows discovered during adversarial verification
rather than by a sector sweep.

**Sectors with zero Lean coverage:** `hashing-digests` (confirmed by grep over
all 19 `.lean` files: no occurrence of `key_hash`, `Indexable`, `digest`, or
`IdentityHasher`). Two sub-areas are equally Lean-free without being whole
sectors: the entire `ppvm-tableau-2/src/mixture/**` subtree (branch weights,
sampler, state equality, fingerprint deltas, mixture noise and loss) and
`ppvm-tableau-2`'s trajectory noise samplers.

## Adjudications (round 2)

Round 2 took six clusters of suspected live defects out of the ledger and
re-derived the answer to each from the algebra, from scratch, twice: every
verdict below was independently re-derived by a second agent that wrote its own
probes and did not inherit the first's reasoning. This section is the record, so
that the next close round proves the **adjudicated** answer rather than cementing
whatever the code currently does. Where a row's `Proposed closure` would have
proved the wrong statement, that closure has been edited in place and points
here.

**Score.** Six units examined; six corroborated by independent re-derivation
(every second opinion came back `agree-with-corrections` — no unit is contested
at the verdict level); three units carry live defects (U1 correlated loss, U4
mixture weights/sampler, U6 unenforced preconditions). **U1's fix DIRECTION was
reversed after the round** — the paper draft entered evidence, contradicted the
in-repo documentation the unit had settled on, and was ruled normative, so the
defect is `ppvm-tableau-2` rather than `ppvm-pauli-sum-2`; the row count is
unchanged. See the ruling box under [U1](#u1-correlated-loss-convention). The disagreements that did arise are all
inside proposed *fixes*, and they are recorded as such under each unit's
**Contested inside the fix** heading rather than smoothed into a consensus that
does not exist. **Nothing has been applied to any crate**; every fix that changes
observable behaviour is flagged as needing the user's sign-off.

### U1 correlated-loss convention

Rows: **G-040** (`adjudicated-defect`), **G-035** (`adjudicated-defect`),
**G-043** (`adjudicated-defect`), **G-045** (`adjudicated-spec`).

> **CORRECTION + RULING (post-round — read this before acting on U1).**
> This unit was adjudicated **without the paper draft**
> (`../ppvm-paper/main.tex`), which is the definitions-of-record and which the
> unit's own conclusion depends on. With it in evidence, the verdict below is
> **withdrawn on the direction of the fix**. The unit's *derivation* stands unaltered and is
> correct — including its central finding that the two readings are the same CPTP
> family reparameterized by a factor of two, so **channel mathematics cannot
> pick one**. What it could not know is that the repository has *two*
> authoritative sources and they disagree:
>
> | Source | meaning of `p[1]` | agrees with |
> | --- | --- | --- |
> | Paper `main.tex:462`, `:523`, `:845` — "$p_{LQ}$ … the probability that both remain is $1-p_{LL}-2p_{LQ}$", $\mathcal{E}^\dagger[P_1P_2]=(1-2p_{LQ}-p_{LL})P_1P_2$, and the demo's $p_{LL}+2p_{LQ}=p$ | a *named* one lost; P(exactly one) $=2p_1$ | `ppvm-pauli-sum-2/src/loss.rs:137` (and legacy `ppvm-pauli-sum`), Lean `corrT` |
> | Shipped Python API `mixins.py:507` ("losing **exactly one** qubit … 50/50 random") + `test/generalized_tableau/test_loss.py:82`, `:173` | P(exactly one) $=p_1$ | `ppvm-tableau-2` mixture (`mixture/noise/loss.rs:95`) and trajectory (`noise.rs:365`) |
>
> So the in-repo evidence the unit weighed (trait doc, `mixins.py`, the Python
> suite) is real and verified — and it is contradicted by the paper.
>
> **RULED (2026-08-10, project lead): the paper is normative.** `p[1]` = $p_{LQ}$ =
> the probability that a *named* one of the pair is lost, so
> P(exactly one) $= 2p_1$ and the both-present survivor scales by $1-2p_1-p_0$.
> Therefore:
>
> * **`ppvm-pauli-sum-2/src/loss.rs` is CORRECT and does not change.** Legacy
>   `ppvm-pauli-sum` is correct. Lean's `corrT` (`Noise.lean:422`) is correct and
>   is transcribed as-is — G-034's closure is unblocked.
> * **`ppvm-tableau-2` is the defect on both backends** and is what gets fixed:
>   `mixture/noise/loss.rs:95` (survivor $1-p_0-p_1$ → $1-p_0-2p_1$; the two
>   single-loss branches $p_1/2$ → $p_1$) and `noise.rs:365` (the cumulative scan
>   over `p[..2]` must put $2p_1$, not $p_1$, on the exactly-one event).
> * **The documentation is unified on the paper's wording**, normatively in
>   `ppvm-traits-2/src/gates.rs:579` and cited from `ppvm-tableau-2/src/noise.rs:337`,
>   `ppvm-python/src/ppvm/mixins.py:507` and `paulisum.py:474` (the latter's
>   ambiguity is how the split survived undetected).
> * **Three shipped Python tests encode the wrong convention and must be
>   rewritten**: `test/generalized_tableau/test_loss.py:82`, `:90`, `:173`. Note
>   `:82`'s `p = [0, 1, 0]` is **inadmissible** under the ruling
>   ($p_0 + 2p_1 = 2 > 1$); the equivalent admissible witness for
>   "exactly one lost in every trial" is `p = [0, 0.5, 0]`.
> * **G-043's admissible region is $p_0 + 2p_1 \le 1$**, $p_0,p_1 \ge 0$,
>   $p_2 \in [0,1]$ — the original ledger hypothesis, restored.
>
> Two consequences the ruling resolves. (a) The paper's transport demo
> (§`sec:transport`, Pauli propagation) states the **exact** result
> $m(t)=(1-p)^{k(t)}$ from $[p_{LL},p_{LQ},p_{LN}]=[p/3,p/3,p/3]$, which requires
> $p_{LL}+2p_{LQ}=p$ — the pauli-sum reading, so that figure is consistent with
> the ruling and with the code that produced it. (b) The paper's MSD demo
> (§`sec:demos`, generalized tableau) sets $p_{LL}=p_{LN}=p_{LQ}=p_\mathrm{loss}/3$
> on a backend that reads `p[1]` as the total, which would make its per-two-qubit-gate
> loss $2p_\mathrm{loss}/3$ rather than $p_\mathrm{loss}$. If both figures were
> produced through the documented API, the two demos in the paper are using two
> different conventions — and under the ruling the MSD figure's effective per-CZ
> loss is $2p_\mathrm{loss}/3$, not $p_\mathrm{loss}$. **The MSD figure and its
> caption should be re-checked before publication**; fixing `ppvm-tableau-2` makes
> the code match the paper, but any number already produced by the old tableau
> path was generated under the other convention.
>
> Everything else in this unit — the trace-preservation derivation, the refutation
> of `loss.rs:106-111`'s "suspected old bug 4", the *existence* of the
> cross-backend factor-of-two split, and G-035's missing-guard finding — is
> unaffected and stands.

Verdict (as recorded by the unit, direction now withdrawn — see the correction
above): **code is wrong** — one backend answers a different question than the
documented API. Second opinion: **corroborated** (independent derivation, probe
numbers identical to the digit), with three corrections, all recorded below.

**Question.** Does `p[1]` in `correlated_loss_channel` mean P(exactly one atom
lost | both in the qubit subspace) — the reading `ppvm-tableau-2`'s mixture and
trajectory implement — or P(a *named* one lost), so that P(exactly one) = 2·p₁ —
the reading `ppvm-pauli-sum-2` and Lean's `corrT` implement?

**Derivation.** Per site the physical space is (2-dim qubit subspace) ⊕
span{|L⟩}. The observable basis is I (identity on the qubit block, zero on |L⟩),
X, Y, Z and L = |L⟩⟨L|, so the full-space identity splits as 𝟙 = I + L and trace
preservation for the Heisenberg map Λ\* is Λ\*(𝟙⊗…⊗𝟙) = 𝟙⊗…⊗𝟙.

*Single site.* Λ_loss(ρ) = |L⟩⟨L|·Tr ρ has adjoint Λ\*_loss(A) = ⟨L|A|L⟩·𝟙, so
Λ\*_loss kills I, X, Y, Z and sends L ↦ I + L. For Λ = (1−p)·id + p·Λ_loss:
Λ\*(k) = (1−p)k for k ≠ L and Λ\*(L) = L + p·I — an **unscaled survivor plus a
branch onto the loss-cleared key at weight p**, which is exactly
`ppvm-pauli-sum-2/src/loss.rs:82-95` and exactly Lean's `lossT`. Derived, not
assumed; the single-qubit loss arm is right on both backends.

*Two sites.* Decompose into the four orthogonal sectors S_qq, S_qL, S_Lq, S_LL (a
fixed decomposition, so a sector-wise channel is linear and CP). Free
parameters: a = P(both | both present), b = P(one *named* site | both present),
c = P(survivor lost | one already lost). Using Λ\*_{loss1}(x⊗y) = ⟨L|y|L⟩(x⊗𝟙),
Λ\*_{lossboth}(x⊗y) = ⟨L|x|L⟩⟨L|y|L⟩ 𝟙⊗𝟙 and Λ\*(A) = Σ_s Π_s Λ\*_s(A) Π_s:

- (x,y) with x,y ∈ {I,X,Y,Z}: every loss adjoint carries ⟨L|·|L⟩ = 0 and every
  non-qq projection vanishes ⇒ Λ\*(x⊗y) = (1 − a − 2b)·(x⊗y).
- (x,L), x ≠ L: S_qq contributes b·(x⊗I), S_qL contributes (1−c)(x⊗L), the other
  sectors 0 ⇒ Λ\*(x⊗L) = (1−c)(x⊗L) + b·(x⊗I). **The gain-branch weight on the
  one-already-lost row is b — the both-present single-loss weight for that site —
  not c.**
- (L,L): Λ\*(L⊗L) = (L⊗L) + c(I⊗L) + c(L⊗I) + a(I⊗I), survivor unscaled.

Column check: the I⊗I column sums (1−a−2b) + b + b + a = 1; the I⊗L and L⊗I
columns (1−c) + c = 1; the L⊗L column 1. **Trace preserving for every parameter
value, with no normalization hypothesis.**

Two consequences.

*(1) The "structural divergence" G-040 flags is not a second bug.*
pauli-sum-2's one-already-lost arm emitting a *gain* branch, and its (L,L) arm
emitting *recovery* branches at p₂ with an unscaled survivor, are precisely the
derived Heisenberg rows; the mixture's "lose the survivor at p₂" is the same
matrix element read as a Schrödinger column. T = Sᵗ on the loss sector:
S[qq→qL] = b ↔ T[(I,L)→(I,I)] = b, and S[qL→LL] = c ↔ T[(L,L)→(I,L)] = c. Both
directions are correct. Therefore `ppvm-pauli-sum-2/src/loss.rs:106-111`'s
"suspected old bug 4" — the arms "do not pair" and "the channel is not
trace-preserving unless `p[1] == p[2]`" — is **false on both counts**; both
agents verified Λ\*(𝟙) = 𝟙 exactly for the shipped arms at p = [0.07,0.11,0.19],
[0.5,0.02,0.25] and [0,1,0]. TP holds column-wise for every p.

*(2) The only real disagreement is b.* pauli-sum-2 takes b = p₁ (survivor
1−p₀−2p₁); the mixture and the trajectory take b = p₁/2 (survivor 1−p₀−p₁, two
branches at p₁/2; resp. a categorical draw over (p₀, p₀+p₁) then a fair coin).
Both are self-consistent CPTP families — the same family reparameterized by a
factor 2 — so channel mathematics alone cannot pick one. It is a
*parameterization convention*, and three in-repo arguments settle it:

- **Domain.** With b = p₁/2 the vector (p₀,p₁) is sub-stochastic over the
  disjoint named events {both, exactly one}, admissible on the standard simplex
  p₀,p₁ ≥ 0, p₀+p₁ ≤ 1 — the same convention G-035 is about, and the only region
  on which the trajectory's cumulative-sum scan over `p[..2]` is a valid
  categorical sampler. With b = p₁ the region is the non-standard p₀+2p₁ ≤ 1, so
  the API's own p₁ = 1 is inadmissible.
- **Documentation.** `ppvm-traits-2/src/gates.rs:577` — "p[1]: losing **either
  one** qubit when both are in the qubit subspace", i.e. the union of the two
  disjoint single-loss events; `ppvm-tableau-2/src/noise.rs:337` repeats it;
  `ppvm-python/src/ppvm/mixins.py:508` is unambiguous — "probability of losing
  **exactly one** qubit when both are in the qubit subspace (which qubit is lost
  is 50/50 random)". (`paulisum.py:474` says only "losing a single qubit" —
  ambiguous, not a counter-witness.)
- **Published contract.** `ppvm-python/test/generalized_tableau/test_loss.py:82`
  (p = [0,1,0] ⇒ "exactly one qubit lost in every trial"), `:90` (50/50 split),
  `:173` ("P(exactly one lost) should converge to p[1]" at 0.4).

**Correct answer.** `p[1]` = P(exactly one atom lost | both in the qubit
subspace), split 50/50 between the two qubits. Correct Heisenberg transfer:
both-present row × (1 − p₀ − p₁); one-already-lost row × (1 − p₂) plus a gain
branch onto the recovered key at **p₁/2**; both-lost row unscaled with branches
p₂, p₂, p₀. Admissible region p₀, p₁ ≥ 0, p₀ + p₁ ≤ 1, p₂ ∈ [0,1].

**Which implementation is right.** `ppvm-tableau-2`, both backends —
`mixture/noise/loss.rs:95` (survivor 1−p₀−p₁, two branches at p₁/2) and
`noise.rs:365` (cumulative sum over `p[..2]`, then a fair coin) — which is what
`ppvm-traits-2/src/gates.rs:577`, `ppvm-python/src/ppvm/mixins.py:508` and the
shipped Python suite document. `ppvm-pauli-sum-2/src/loss.rs:137,154,158` is
wrong by a factor of two on the single-loss mass, and Lean's `corrT`
(`Noise.lean:422`) transcribes that wrong convention, so
`correlatedLossChannel_trace_preserving` currently machine-checks the wrong
channel: it is TP, for the wrong parameterization. Both maps are individually
CPTP, so this is a parameterization defect, not an arithmetic one — but it is a
defect, because it makes one backend answer a different question than the
documented API.

**Live defect: YES — a wrong number, through `ppvm-pauli-sum-2`'s public
`CorrelatedLossChannel` impl.** At p = [0, 0.3, 0] every in-subspace observable
coefficient is multiplied by 0.4 where the documented convention (and both
tableau-2 backends) give 0.7 — a 43% under-report of e.g. ⟨ZZ⟩ after a
correlated-loss event. At p = [0, 1, 0] — the exact parameter the shipped Python
test documents as "exactly one qubit lost in every trial" — pauli-sum-2 returns
coefficient **−1.0**, a sign-flipped expectation value where the correct answer
for an in-subspace key is 0. The two backends disagree on the same circuit and
the disagreement grows with p₁. Secondary, mislabelling only, no wrong number:
`loss.rs:106-111`'s "suspected old bug 4" alleges two non-existent defects, and
`Noise.lean:438-460`'s prose presents `corrT` as the spec for all three
backends.

**Positivity (G-043, G-035).** CP ⟺ every event weight ≥ 0, since the map is then
a convex mixture of the CP maps id, loss0, loss1, lossboth, and a negative weight
makes it non-positive. So the **necessary and sufficient** admissibility
condition is p₀ ≥ 0, p₁ ≥ 0, p₂ ∈ [0,1], p₀ + p₁ ≤ 1 (under the corrected
convention; p₀ + 2p₁ ≤ 1 under the other — which is why G-043's proposed
hypothesis is wrong, see the closure edit there). Nothing states or enforces it:
`correlatedLossChannel_trace_preserving (p0 p1 p2 : ℝ)` has no hypotheses and is
true at (5, −3, 17); `loss.rs` has no assertion; the mixture has none; the
trajectory tests only `p[0] <= 0.0 && p[1] <= 0.0`. Measured: the mixture
silently truncates its negative survivor weight and renormalizes, turning an
invalid input into a plausible-looking distribution, and pauli-sum-2 returns a
coefficient of 2.0 (or −0.8) without complaint.

**G-045, asymmetric loss.** Independent derivation of the shipped `p_tot`: the
loss Kraus family K_L0 = √p0·|L⟩⟨0|, K_L1 = √p1·|L⟩⟨1| with survival
K₀ = √(1−p0)|0⟩⟨0| + √(1−p1)|1⟩⟨1| satisfies
K₀†K₀ + ΣK_Li†K_Li = (1−p0)|0⟩⟨0| + (1−p1)|1⟩⟨1| + p0|0⟩⟨0| + p1|1⟩⟨1| = 1, so
the family is CPTP, and the loss-branch trace is
p0⟨0|ρ|0⟩ + p1⟨1|ρ|1⟩ = p0(1+⟨Z⟩)/2 + p1(1−⟨Z⟩)/2 — verbatim
`ppvm-tableau-2/src/noise.rs:365`'s `p0*0.5*(1+z) + p1*0.5*(1-z)`. The sign
convention the first agent flagged as its one medium-confidence step was probed
by the second: `z_expectation(|0⟩) = +1`, `z_expectation(|1⟩) = −1`, and
`asymmetric_loss_channel(0, p0=1, p1=0)` fires on |0⟩ and not on |1⟩, so p0 is
the |0⟩-rate and the arithmetic is right. The shipped trajectory replaces K₀ by
the scalar √(1−p_tot)·1, dropping the survival back-action; that substitution is
exact iff p0 = p1 **or** the site is in a Z eigenstate, and its error is first
order in |p0−p1| (the conditional survivor's Bloch-z is biased by
≈ (p1−p0)/2·(1−⟨Z⟩²)/2). This is a documented, unavoidable approximation inside a
stabilizer simulator (the true K₀ is non-Clifford), so it is **not** a wrong
number relative to its own docstring: the gap is genuinely coverage (nothing in
the `-2` crates exercises p0 ≠ p1) plus a docstring that understates the
exactness region (Z eigenstates, not only p0 = p1).

**Proposed fix** — behaviour-changing on the pauli-sum-2 path; needs sign-off.

1. `crates/ppvm-pauli-sum-2/src/loss.rs:137` — replace
   `let both_present = C::one() - (p1.clone() + p1.clone()) - p0.clone();`
   with `let half_p1 = p1.half(); let both_present = C::one() - p1.clone() - p0.clone();`
2. `crates/ppvm-pauli-sum-2/src/loss.rs:154` and `:158` — branch weight
   `p1.clone()` → `half_p1.clone()` (the gain from the both-present sector into
   one *specific* single-loss key is p₁/2).
3. Bound change required by (1): add `C: ppvm_traits_2::Halvable` to the
   `CorrelatedLossChannel` impl at `loss.rs:120-126` (`Coefficient` has no `Div`;
   `Halvable::half` at `ppvm-traits-2/src/coefficient.rs:140` is the existing
   capability, implemented for `f64`/`Complex<f64>` and deliberately not for
   exact rings). This narrows the impl away from non-halvable coefficient rings —
   a 50/50 split is genuinely undefinable over ℤ[i], so that is the honest
   consequence; `proj.rs:68` already sets the precedent. Verified by the second
   agent: every `correlated_loss_channel` call site in the tree is `f64`
   (`ppvm-pauli-sum-2/tests/{loss,column_store_gates}.rs`,
   `ppvm-conformance-2/tests/{pauli_sum_loss_diff,pauli_sum_indexmap_diff}.rs`,
   the benches, `ppvm-stim/src/executor/adapter/traits_2.rs`, `ppvm-vihaco`), and
   `ppvm_sym_2::Term` has no `Halvable` impl, so the bound compiles today.
4. Docs, same file: `loss.rs:15` ("`correlated_loss_channel` to `c *= 1 − 2p₁ − p₀`")
   → `1 − p₁ − p₀`; delete the false "suspected old bug 4" paragraph at
   `loss.rs:102-115` and replace it with the derived transfer matrix.
5. `crates/ppvm-traits-2/src/gates.rs:577` — restate "losing either one" as
   "losing exactly one (either qubit, 50/50)", so the convention is stated once,
   normatively, where all three backends can cite it.
6. `lean/PPVM/Algebra/Noise.lean:422` `corrT` — both-present entry becomes
   `1 - p1 - p0`, the two one-already-lost gain entries become `p1 / 2`; the
   trace-preservation proof still goes through
   ((1−p₀−p₁) + p₁/2 + p₁/2 + p₀ = 1). The docstring at `:438-460` should cite
   `ppvm-tableau-2/src/mixture/noise/loss.rs` and `ppvm-pauli-sum-2/src/loss.rs`
   as the two pictures (T = Sᵗ) instead of presenting one arm ordering as "the"
   spec.
7. G-043/G-035, additive: `debug_assert!(p0 >= 0 && p1 >= 0 && p2 >= 0 && p0 + p1 <= 1.0 && p2 <= 1.0)`
   at all three correlated-loss entry points (`ppvm-pauli-sum-2/src/loss.rs:133`,
   `ppvm-tableau-2/src/mixture/noise/loss.rs:59`,
   `ppvm-tableau-2/src/noise.rs:352`) — note the second agent's caveat that on a
   coefficient-generic `C` this will not compile as written (no `PartialOrd`), so
   it must live in the f64 backends or behind the `Halvable`/f64 instantiations —
   plus a new Lean lemma proving `corrT`ᵗ is column-*stochastic* (entries ≥ 0,
   columns sum to 1) under `0 ≤ pᵢ`, `p₀ + p₁ ≤ 1`, which is what makes the map
   CP rather than merely TP, plus a zero-norm guard in
   `ppvm-tableau-2/src/mixture/data.rs:135 normalize_probabilities`.
8. G-045: no code change. Amend `ppvm-tableau-2/src/noise.rs:290-295` to say the
   approximation is exact iff `p0 == p1` *or* the target is in a Z eigenstate, and
   add the missing test at p0 ≠ p1 with ⟨Z⟩ ∈ {+1,−1,0}.

**Do NOT** "fix" this by doubling the tableau backends: that breaks the
documented Python contract (p₁ = 1 becomes inadmissible) and the three shipped
tests that pin P(exactly one) = p[1].

**Evidence.** Scratch probe `crates/ppvm-conformance-2/tests/zz_adj_loss.rs`
(deleted after the run), `cargo test -p ppvm-conformance-2 --test zz_adj_loss --
--nocapture --test-threads=1`. Verbatim:

```
running 4 tests
test probe_no_positivity_guard ... === p = [0.6, 0.6, 0.0]: p0 + p1 = 1.2 > 1 ===
  [pauli-sum-2] II coefficient = Some(-0.7999999999999999)
  [mixture] branches (lost0,lost1,weight): [(false, true, 0.25), (true, false, 0.25), (true, true, 0.5)]
  [mixture] survivor weight = -0, total lost-branch weight = 1
=== p = [5.0, -3.0, 17.0] (the ledger's p0 = 5) ===
  [pauli-sum-2] II coefficient = Some(2.0)
  [mixture] branches (lost0,lost1,weight): [(true, true, 1.0)]
  [mixture] survivor weight = -0, total lost-branch weight = 1
ok
test probe_one_already_lost_arm ... === p = [0, 0.3, 0.11], qubit1 already lost ===
  [pauli-sum-2 IL] terms (lost0,lost1,x0,z0)->coeff: [("(false, false, false, false)", 0.3), ("(false, true, false, false)", 0.89)]
  [pauli-sum-2 IL] sum of coefficients = 1.19
    IL -> Some(0.89),  II -> Some(0.3)
  [mixture after pre-loss] branches (lost0,lost1,weight): [(false, true, 1.0)]
  [mixture after pre-loss] survivor weight = -0, total lost-branch weight = 1
  [mixture after corr] branches (lost0,lost1,weight): [(false, true, 0.89), (true, true, 0.11)]
  [mixture after corr] survivor weight = -0, total lost-branch weight = 1
ok
test probe_p1_equals_one_is_the_documented_python_case ... === p = [0, 1, 0] ===
  [pauli-sum-2] ZZ coefficient = Some(-1.0)
  [pauli-sum-2] at p1=0.3: ZZ coefficient = Some(0.4) (mixture/trajectory imply 0.7)
  [mixture] branches (lost0,lost1,weight): [(false, true, 0.5), (true, false, 0.5)]
  [mixture] survivor weight = -0, total lost-branch weight = 1
ok
test probe_single_loss_mass_p1 ... === p = [0, 0.3, 0], both qubits present ===
  [pauli-sum-2 II] terms (lost0,lost1,x0,z0)->coeff: [("(false, false, false, false)", 0.4)]
  [pauli-sum-2 II] sum of coefficients = 0.4
  [pauli-sum-2] survivor scale = Some(0.4)  => implied P(exactly one lost) = 0.6
  [mixture] branches (lost0,lost1,weight): [(false, false, 0.7), (false, true, 0.15), (true, false, 0.15)]
  [mixture] survivor weight = 0.7, total lost-branch weight = 0.3
  [trajectory] over 200000 seeds: P(exactly one lost) = 0.29912, P(both) = 0
ok

test result: ok. 4 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 1.98s
```

Reading: (i) single-loss mass at p₁ = 0.3 is 0.6 on pauli-sum-2 vs 0.3 on the
mixture and 0.29912 (200k seeds) on the trajectory — **the ledger's proposed
cross-backend oracle fails exactly as it predicted, which is why that oracle
must not be landed as written**; (ii) the one-already-lost row on pauli-sum-2 is
(1−p₂) = 0.89 survivor + p₁ = 0.3 gain branch while the mixture's same matrix
element appears as 0.11 = p₂ of loss-of-survivor — the Heisenberg/Schrödinger
transpose, with the pauli-sum branch weight 0.3 where the corrected convention
wants 0.15 (the row sum 1.19 ≠ 1 is expected and harmless: TP is column-wise);
(iii) no positivity guard anywhere — p = [5, −3, 17] is accepted silently by both
backends, the mixture truncating its negative survivor and renormalizing to "both
lost with probability 1"; (iv) p = [0, 1, 0] yields a −1.0 coefficient on
pauli-sum-2.

**Second opinion — corroborated, with three corrections.** The second agent
re-derived the sector-pinched channel, the (x,L) row's gain weight being b rather
than p₂, TP column-wise for every p, both conventions being valid CPTP families,
and p₁ = exactly-one winning on documentation + domain + shipped tests; its probe
reproduced 0.4 vs 0.7, −1.0 at p₁ = 1, 0.29912 over 200k seeds, −0.8 and +2.0
with no positivity guard, and the mixture's truncate-and-renormalize at invalid
p. Corrections, none of which changes the verdict:

1. **The fix breaks the port's golden-master oracle, which the first write-up
   omits.** Legacy `crates/ppvm-pauli-sum/src/sum/noise.rs:241` is
   `*v *= 1.0 − p[1]*2.0 − p[0]` with the same p₁ gain branch at `:221`/`:233`, so
   `crates/ppvm-conformance-2/tests/pauli_sum_loss_diff.rs:92`
   (`correlated_loss_matches_all_four_loss_arms`, bit-exact at p = [0.07,0.11,0.19])
   and `crates/ppvm-conformance-2/tests/pauli_sum_indexmap_diff.rs:156`
   (`loss_and_correlated_loss_match_observable_order`, same p) **will fail**, as
   will the two pinned unit tests `crates/ppvm-pauli-sum-2/tests/loss.rs:224`
   (asserts `1 − 2·0.02 − 0.01`) and `:236` (asserts branch = p₁ = 0.02). Four
   tests, two of them the behaviour-preserving-port contract itself. The user must
   decide whether pauli-sum-2 may diverge from legacy or whether legacy is fixed
   in the same commit.
2. **Reachability is narrower in one direction and wider in another.**
   `crates/ppvm-python-native/Cargo.toml` binds the legacy crates only, so today
   the `-2` divergence is reachable from Rust (the `ppvm-pauli-sum-2` public impl,
   and `ppvm-stim/src/executor/adapter/traits_2.rs:108-116` which forwards `p`
   verbatim), not from Python. But the identical 2p₁ convention is live in shipped
   legacy `ppvm-pauli-sum`, so a Python user today already gets 2p₁ from
   `LossyPauliSum.correlated_loss_channel` and p₁ from `GeneralizedTableau`: the
   cross-backend disagreement is user-visible now, and the `-2` port inherits it.
3. **G-045's medium-confidence caveat is resolved in the first agent's favour** by
   the ⟨Z⟩ probe quoted above.

**Contested inside the fix.** Nothing in the ruling. The open judgement calls are
(a) legacy divergence vs fixing legacy in tandem (item 1 above) and (b) whether
narrowing `CorrelatedLossChannel` to `Halvable` coefficients is acceptable.

**Sign-off needed.** Items 1-3 (the p₁/2 arithmetic and the `Halvable` bound) and
item 7's asserts change observable behaviour / break four tests.

### U2 rotation direction and sign

Rows: **G-030**, **G-029**, **G-025**, **G-028** — all `adjudicated-spec`.
Verdict: **the Lean and the docstrings are wrong; the code is right.** Second
opinion: **corroborated** (independent derivation, independently written probe,
identical numbers), with four corrections — one of which makes the proposed fix
incomplete.

**Question.** For R = exp(−iθG/2) with {G,P} = 0, is the branch coefficient
+sinθ·(iGP) attached to the right conjugation direction, and is the shipped
`ppvm-pauli-sum-2` rotation therefore numerically correct or only mislabelled?

**Derivation.** G, P Hermitian Pauli words, G² = P² = I, GP = −PG. Since G² = I,
R = exp(−iθG/2) = cI − isG exactly, with c = cos(θ/2), s = sin(θ/2), R† = cI + isG.

R P R† = (cI − isG)P(cI + isG) = c²P + ics·PG − ics·GP + s²·GPG.
GPG = (−PG)G = −P·G² = −P, so the last term is −s²P; the diagonal is
(c²−s²)P = cosθ·P; the cross term is ics(PG − GP) = −2ics·GP = −sinθ·(iGP).
Hence

- **R P R† = cosθ·P − sinθ·(iGP)** and, by θ ↦ −θ,
- **R† P R = cosθ·P + sinθ·(iGP)**.

So the RHS `c·cosθ·P + c·sinθ·(iGP)` belongs to R†·(cP)·R, **not** to R·(cP)·R†.
G-030's algebra is confirmed. Explicit 2×2 anchors (Y = [[0,−i],[i,0]], so
iXZ = Y, iZX = −Y, iXY = −Z, iYX = +Z, iYZ = −X, iZY = +X): R Z R† with G = X is
[[cosθ, i sinθ],[−i sinθ, −cosθ]] = cosθ·Z − sinθ·Y ✓; R X R† with G = Z is
[[0, e^{−iθ}],[e^{iθ},0]] = cosθ·X + sinθ·Y = cosθ·X − sinθ·(iZX) ✓.

*What the shipped ε columns emit* (read from source, not from legacy or tests):
`rx` (`rotation.rs:225-231` → `store.rs:373-390`) fires iff the z-bit is set,
ε = −1 if x else +1, toggles x ⇒ Z ↦ +sinθ·Y, Y ↦ −sinθ·Z. `ry`
(`rotation.rs:234-256`) fires iff x ≠ z, ε = −1 if z else +1, toggles both ⇒
X ↦ +sinθ·Z, Z ↦ −sinθ·X. `rz` (`rotation.rs:259-276`) fires iff x is set,
ε = +1 if z else −1, toggles z ⇒ X ↦ −sinθ·Y, Y ↦ +sinθ·X. All six anticommuting
(G,P) pairs equal cosθ·P + sinθ·(iGP), i.e. **the code computes R† P R.** Two-site
by hand: `rzz` on X_aI_b gives −sin·Y_aZ_b and i(Z_aZ_b)(X_aI_b) = i(iY)_a⊗Z_b =
−Y_aZ_b ✓; `rxx` on Z_aI_b gives +sin·Y_aX_b and i(XZ)_a⊗X_b = +Y_aX_b ✓; `ryy`
on Z_aI_b gives −sin·X_aY_b and i(YZ)_a⊗Y_b = −X_aY_b ✓.

*Generic path and `levi_civita` (G-028).* `comm_2` is documented as
[Q,P]/2i = −i[G,P]/2 with G first; for {G,P} = 0, [G,P] = 2GP, so the coefficient
of iGP is −ε, and `rotate_2` multiplying sin by `−eps` is correct — which also
explains the +eps/−eps asymmetry against the native kernels. `levi_civita` is
documented as −i[P_i,P_j]/2 = ε·P_k and `rotate_1_branch` calls it as
(key, axis), i.e. i = P, j = G; then −i[P,G]/2 = +iGP for anticommuting pairs, so
+ε *is* the coefficient of iGP and `sin.mul_sign(eps)` is right. Hand check:
i = Z (0b10), j = X (0b01): rank(Z) = 2, rank(X) = 0, diff = (0−2) mod 3 = 1 ⇒
ε = +1, k = 0b11 = Y; and −i[Z,X]/2 = −i(2iY)/2 = +Y = iXZ ✓. **The key-first
argument order is load-bearing and currently correct**; swapping it negates every
ε.

*Why R†PR is the right answer here and not a free convention.* A Pauli sum in
this crate is an **observable**. Composing the per-gate maps in application order
gives O ↦ (R₁R₂…Rₙ)† O (R₁R₂…Rₙ), so feeding a circuit's gates in reverse order
yields U†OU, whose expectation on |ψ⟩ equals ⟨O⟩ on U|ψ⟩ — standard Heisenberg
back-propagation, which is exactly what `rotation.rs:451` documents ("a `Sum`
propagates observables backward"). `RotXY` corroborates it independently:
R(φ,θ) = exp(−iθ/2·(cosφ X + sinφ Y)) = Rz(φ)Rx(θ)Rz(−φ) because
e^{−iφZ/2} X e^{+iφZ/2} = cosφ X + sinφ Y; read right-to-left as a circuit, the
backward walk is exactly `rz(φ); rx(θ); rz(−φ)`, the shipped order
(`rotation.rs:474-476`), and at φ = π/2 the generator is Y so
`r(q,π/2,θ) == ry(q,θ)` — the documented consequence. A forward-conjugating
implementation with the same statement order would give `ry(−θ)`.

**Correct answer.** R P R† = cosθ·P − sinθ·(iGP); R† P R = cosθ·P + sinθ·(iGP).
The shipped kernels emit cosθ·P + sinθ·(iGP), i.e. R† P R — correct Heisenberg
(backward) observable propagation, matching `rotation.rs:451`. The statements that
attach that RHS to the LHS `e^{−iθG/2}·(cP)·e^{+iθG/2}` are the wrong ones.

**Which implementation is right.** All of them, and they agree. `HashMapStore`
(`store.rs:1421`), the columnar kernels (`column_store/rotations/rx.rs`, including
the fused closed-support 2×2 path) and the indexmap backend produce bit-identical
`rx` output over every word at n ≤ 3 and every qubit, and that output equals
R†·mat(A)·R to 0.0. The generic `rotate_2` (`comm_2`, `SIGN_NEG = 0x2840`) and the
native `rzz`/`rxx`/`ryy` agree with R†·mat·R exactly for all axis pairs at n = 2
and n = 3 with non-adjacent sites, and the `levi_civita` lossy fallback agrees
with the native single-site columns on all 4×4×4 cases. **No backend disagreement
exists in this unit**; the disagreement is between the code (backward) and the
prose (forward).

**Live defect: NO — mislabelling only.** No wrong number reaches a user of the
`-2` crates: every ε column, both two-site paths, all three stores and the lossy
fallback implement one and the same map O ↦ R†OR, which is the direction the
user-facing `RotXY` docstring (`rotation.rs:451`) and the behavioural contract
`r(q,π/2,θ) == ry(q,θ)` document. It is also coherent with the crate's Clifford
family: `Sum::s` maps X ↦ −Y = S†XS, identical to `rz(π/2)` (probe diff
6.1e-17), so mixing Clifford and rotation gates in one circuit stays coherent.
(`ppvm-tableau-2` uses cos(θ/2)/−i·sin(θ/2) on state amplitudes — Schrödinger
picture, a different object, no cross-crate conflict.) The residual risk is
documentation-driven and is exactly why this row is worth correcting *before*
proving: a reader who trusts `rotation.rs:9`, `store.rs:306-308` or
`Rotation.lean:17` and writes a Lean proof, a new backend, or a caller that feeds
gates in forward circuit order gets a sign-inverted answer, and the next round's
Lean would cement the false identity.

**Proposed fix** — prose only; no ε column, no key, no order changes; no
observable behaviour change.

1. `lean/PPVM/Instantiations/Rotation.lean:17` — replace
   `` `e^{-iθG/2} · (c · P) · e^{iθG/2} = c·cos θ · P + c·sin θ · (iGP)` `` with
   `` `R† · (c · P) · R = c·cos θ · P + c·sin θ · (iGP)`  with  `R = e^{-iθG/2}` ``
   (the second agent's wording, preferred over writing the exponent as
   `e^{+iθG/2}`, which would make the *gate* look like `e^{+iθG/2}` when it is
   not), plus one clause: "— the **backward (Heisenberg)** direction, which is
   what the crate propagates; the forward map `R P R†` is the same with `sin θ`
   negated." Same edit for the back-reference at `:31`. **No theorem changes**:
   `rx_/ry_/rz_/rzz_/rxx_/ryy_eps_from_product`, `branchExp`, the `mz`/`mx`/`my`
   matrices at `:596-610` and `rotAxis` (already "Rodrigues' formula with angle
   −θ") are all already the backward action and match the matrix computation term
   for term. Only the header prose is false.
2. `crates/ppvm-pauli-sum-2/src/rotation.rs:8-10` — replace "conjugates each
   stored term `(P, c)` to `c·cosθ·P + c·sinθ·(iGP)`" with "propagates each stored
   term `(P, c)` **backward** (Heisenberg, `R†·P·R` with `R = exp(−i·θ/2·G)` — the
   direction documented on [`RotXY`] below) to `c·cosθ·P + c·sinθ·(iGP)`; the
   forward map `R·P·R†` is the same with `θ ↦ −θ`."
3. **Added by the second opinion, and required for the fix to be complete:**
   `crates/ppvm-pauli-sum-2/src/store.rs:306-308`, the `RotateInPlace` docstring
   ("A single-qubit rotation `exp(−i·θ/2·G)` conjugates each term `(P, c)` to
   `(P, c·cosθ)` **plus** … a genuinely-new branch `(iGP, c·sinθ·ε)`"), makes the
   same undirected claim about the same generator. `RotateInPlace` is *precisely*
   the trait a new backend implements, so leaving it misleads exactly the reader
   this row protects. Add the same backward/`R†PR` clause. (The
   direction-free statements at `store.rs:345`, `ppvm-traits-2/src/gates.rs:229`,
   `ppvm-traits-2/src/word.rs:104`, `ppvm-tableau-2/src/gates.rs:15` and
   `lean/README.md:124` are fine as-is: they name no direction and inherit
   whatever the fixed headers say.)

Out of unit but the same defect class, flagged for the user:
`crates/ppvm-pauli-sum-2/src/clifford.rs:237` states `SXS† = −Y`, `SYS† = X`;
mathematically S X S† = diag(1,i)·X·diag(1,−i) = **+Y** and S Y S† = **−X** (both
agents verified independently; the probe distance from the shipped `s` to S†PS is
0.000 and to S P S† is 2.000). The code is right — it ships X ↦ −Y = S†XS,
consistent with `rz(π/2)` — and the docstring should read `S†XS = −Y`,
`S†YS = X`, `S†ZS = Z`.

For the coverage/bridge rows: G-025's proposed oracle values are correct as
written (`rz` on X gives {X: cosθ, Y: −sinθ}) and G-030's proposed
`mat(rx(θ)·A) == R†·mat(A)·R` oracle is the right one — both agents ran exactly it
and it passes. G-028's transcription task must pin the argument order explicitly:
`levi_civita(key, axis)` with `+ε` = coefficient of `iGP`.

**Evidence.** Probe `crates/ppvm-conformance-2/tests/zz_adj_rot.rs` (dense
2ⁿ×2ⁿ ℂ oracle built from Y = [[0,−i],[i,0]] and R = cos(θ/2)I − i sin(θ/2)G;
deleted after the run). `cargo test -p ppvm-conformance-2 --test zz_adj_rot --
--nocapture`, verbatim:

```
=== T0 site order ===
word 'XZI': x_bit(0)=true z_bit(0)=false x_bit(1)=false z_bit(1)=true to_string=XZI
=== T1 single-qubit grid, theta = 0.7 ===
cos(theta) = 0.764842187284, sin(theta) = 0.644217687238
rx(theta) on I: emitted [("I", 1.0)]  | diff vs R P Rd = 0.000e0 | diff vs Rd P R = 0.000e0
rx(theta) on X: emitted [("X", 1.0)]  | diff vs R P Rd = 0.000e0 | diff vs Rd P R = 0.000e0
rx(theta) on Y: emitted [("Y", 0.7648421872844885), ("Z", -0.644217687237691)]  | diff vs R P Rd = 1.288e0 | diff vs Rd P R = 1.110e-16
rx(theta) on Z: emitted [("Y", 0.644217687237691), ("Z", 0.7648421872844885)]  | diff vs R P Rd = 1.288e0 | diff vs Rd P R = 1.110e-16
ry(theta) on I: emitted [("I", 1.0)]  | diff vs R P Rd = 0.000e0 | diff vs Rd P R = 0.000e0
ry(theta) on X: emitted [("X", 0.7648421872844885), ("Z", 0.644217687237691)]  | diff vs R P Rd = 1.288e0 | diff vs Rd P R = 1.110e-16
ry(theta) on Y: emitted [("Y", 1.0)]  | diff vs R P Rd = 0.000e0 | diff vs Rd P R = 0.000e0
ry(theta) on Z: emitted [("X", -0.644217687237691), ("Z", 0.7648421872844885)]  | diff vs R P Rd = 1.288e0 | diff vs Rd P R = 1.110e-16
rz(theta) on I: emitted [("I", 1.0)]  | diff vs R P Rd = 0.000e0 | diff vs Rd P R = 0.000e0
rz(theta) on X: emitted [("X", 0.7648421872844885), ("Y", -0.644217687237691)]  | diff vs R P Rd = 1.288e0 | diff vs Rd P R = 1.110e-16
rz(theta) on Y: emitted [("X", 0.644217687237691), ("Y", 0.7648421872844885)]  | diff vs R P Rd = 1.288e0 | diff vs Rd P R = 1.110e-16
rz(theta) on Z: emitted [("Z", 1.0)]  | diff vs R P Rd = 0.000e0 | diff vs Rd P R = 0.000e0
=== T2 randomized single-qubit conjugation ===
worst |mat(result) - Rd M R| over all cases = 6.661e-16
best  |mat(result) - R M Rd| over all cases = 0.000e0
=== T3 two-site rotations (generic rotate_2 + native kernels) ===
worst |mat(rotate_2/native) - Rd M R| = 0.000e0
best  |mat(rotate_2/native) - R M Rd| = 0.000e0
=== T4 RotXY ===
key I: r(pi/2,th)=[("I", 1.0)]
         ry(+th)   =[("I", 1.0)]  diff=0.000e0
         ry(-th)   =[("I", 1.0)]  diff=0.000e0
key X: r(pi/2,th)=[("X", 0.6748757600712672), ("Y", 1.990811798769694e-17), ("Z", 0.7379313711099627)]
         ry(+th)   =[("X", 0.6748757600712672), ("Z", 0.7379313711099627)]  diff=1.991e-17
         ry(-th)   =[("X", 0.6748757600712672), ("Z", -0.7379313711099627)]  diff=1.476e0
key Y: r(pi/2,th)=[("X", 1.990811798769694e-17), ("Y", 1.0), ("Z", -4.518526458101167e-17)]
         ry(+th)   =[("Y", 1.0)]  diff=4.519e-17
         ry(-th)   =[("Y", 1.0)]  diff=4.519e-17
key Z: r(pi/2,th)=[("X", -0.7379313711099627), ("Y", 4.518526458101167e-17), ("Z", 0.6748757600712672)]
         ry(+th)   =[("X", -0.7379313711099627), ("Z", 0.6748757600712672)]  diff=4.519e-17
         ry(-th)   =[("X", 0.7379313711099627), ("Z", 0.6748757600712672)]  diff=1.476e0
worst |mat(r(phi,theta) applied) - Rd(phi,theta) P R(phi,theta)| = 1.241e-16
=== T5 rz(pi/2) vs Clifford s, and s vs S P Sd / Sd P S ===
key I: s -> [("I", 1.0)] | rz(pi/2) -> [("I", 1.0)] | diff(s, rz) = 0.000e0 | diff(s, S P Sd) = 0.000e0 | diff(s, Sd P S) = 0.000e0
key X: s -> [("Y", -1.0)] | rz(pi/2) -> [("X", 6.123233995736766e-17), ("Y", -1.0)] | diff(s, rz) = 6.123e-17 | diff(s, S P Sd) = 2.000e0 | diff(s, Sd P S) = 0.000e0
key Y: s -> [("X", 1.0)] | rz(pi/2) -> [("X", 1.0), ("Y", 6.123233995736766e-17)] | diff(s, rz) = 6.123e-17 | diff(s, S P Sd) = 2.000e0 | diff(s, Sd P S) = 0.000e0
key Z: s -> [("Z", 1.0)] | rz(pi/2) -> [("Z", 1.0)] | diff(s, rz) = 0.000e0 | diff(s, S P Sd) = 0.000e0 | diff(s, Sd P S) = 0.000e0
=== T6 columnar / indexmap backends ===
all three backends agree on rx; worst matrix diff (backward) = 0.000e0
=== T7 lossy rotate_2 fallback (levi_civita) vs rotate_1 ===
all 4x4x4 lossy fallbacks agree with the native single-site columns
test result: ok. 8 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.04s
```

Coverage: T1 = 3 axes × 4 keys at n = 1; T2 = 5 seeds × n ∈ {1,2,3} × 12 random
multi-term sums × 3 axes × random θ ∈ (−3.1, 3.1); T3 = n ∈ {2,3}, sites
(0,1)/(0,2), all 9 axis pairs × *every* Pauli word (16 and 64) through both the
generic `rotate_2` and the native kernels; T4 = φ ∈ {0, 0.3, π/2, 1.9, −0.7};
T6 = `HashMapStore` vs `ColumnStore` vs `IndexMapStore` over every word at n ≤ 3
and every qubit; T7 = `LossyPauliSum::rotate_2` with site 0 lost, all 4 × 4 × 4
cases. The 1.288e0 = 2·sin(0.7) figures are the smoking gun: **the forward model
is off by exactly twice the sinθ branch.**

The second agent's independently written probe
(`crates/ppvm-conformance-2/tests/zz_2nd_zz_adj_rot.rs`, also deleted) reproduced
these to the digit: 1.110e−16 backward vs 1.288e0 = 2·sin(0.7) forward on every
branching case; worst two-site diff 5.551e−17 with the native and generic paths
bit-equal; `r(0,π/2,θ)` vs `ry(0,θ)` 2.0e−17/4.5e−17 and vs `ry(0,−θ)`
1.476e0 = 2·sin(0.83); `s` distance 0.000 to S†PS and 2.000 to S P S†;
`ColumnPauliSum` (both the dense-support batch path and the singleton small-row
path) and `IndexPauliSum` equal to `HashMapStore` within 1e−14 over every word at
n ≤ 3, worst vs R†MR 8.882e−16.

**Second opinion — corroborated, with four corrections.**

1. **The fix is incomplete** — `store.rs:306-308` makes the same claim; folded
   into item 3 of the proposed fix above.
2. **Precision correction to the ledger's own wording, which G-030 repeats.**
   G-030 says the false identity appears "verbatim at `rotation.rs:7-9`". It does
   not: `rotation.rs` never writes the sandwich `e^{-iθG/2}·(cP)·e^{+iθG/2}`; it
   writes "A single-qubit rotation `exp(−i·θ/2·G)` … conjugates each stored term
   `(P, c)` to `c·cosθ·P + c·sinθ·(iGP)`", which is *ambiguous and misleading*,
   not a literally false equation. The only literally false equation in the repo
   is `Rotation.lean:17` (with `:31` referring back to it). The honest tally is
   **one false Lean equation plus three direction-ambiguous Rust docstrings**, not
   two verbatim-false statements — G-030's row text has been amended accordingly.
3. The Lean wording should be `R† · (c·P) · R` with `R = e^{-iθG/2}` rather than
   `e^{+iθG/2} · (c·P) · e^{-iθG/2}` (adopted above).
4. **The prose fix does not close G-025, G-029 or G-030 on the coverage axis.**
   Both agents ran the dense-ℂ matrix oracle and both deleted it. G-025's specific
   complaint — a global flip of any axis column passes every test in the repo,
   including all three copies of the `rx` sign in
   `column_store/rotations/rx.rs` — remains literally true after the prose fix;
   the second agent re-confirmed by construction that nothing in the suite pins
   the absolute ε. Closing these rows still requires landing the ~250-line oracle
   in `pauli_sum_rotation_noise_lean.rs`.

Two out-of-unit observations, offered as observations rather than findings: (i)
the `clifford.rs:237` docstring inversion above; (ii)
`crates/ppvm-vihaco/src/component/dispatch.rs:96-99` maps `T`/`TAdj` to
`rz(±π/8)`, whereas with `rz(θ) = exp(−iθZ/2)` the T gate is `rz(π/4)`. That
lives in the shared dispatch macro (identical for the `legacy` and `traits-2`
backends, so not introduced by this refactor) and is an angle-*magnitude* question
orthogonal to the direction question; it deserves its own row rather than being
folded into this unit.

**Contested inside the fix.** None. Both agents agree the fix is prose-only, safe,
and that "fixing" any sign here would inject a bug.

**Sign-off needed.** None — no observable behaviour changes.

### U3 lossy canonicality

Rows: **G-013**, **G-011**, **G-014**, **G-016** — all `adjudicated-spec`.
Verdict: **the code is right on every reachable path; two public branch builders
carry an undocumented precondition, and G-013's specific empirical claim is
false.** Second opinion: **corroborated** (independent probe reproduced identical
hash values and the identical mutation divergence string), with four corrections,
one of which widens the fix.

**Question.** Is the canonical lossy encoding (`lost ⇒ x = z = 0`) actually
preserved by the shipped `ppvm-lossy-pauli-word-2` / `ppvm-pauli-sum-2` code, is a
non-canonical word observationally distinguishable, and does the ledger's claim
that deleting the `xor_z_col` loss guard survives every test hold?

**Derivation.** A lossy site is an element of S = {I,X,Y,Z} ⊔ {L}, |S| = 5;
physically L is not a Pauli, because on the extended one-site space the identity
splits as 𝟙 = I ⊕ L. The crate encodes a site in three bits with
enc(I) = (0,0,0), enc(X) = (1,0,0), enc(Z) = (0,1,0), enc(Y) = (1,1,0),
enc(L) = (0,0,1), and decodes with dec(x,z,l) = L if l else Pauli(x,z)
(`Word::get`). `dec ∘ enc = id`, but dec is 2-to-1 on l = 1: the three triples
(1,0,1), (0,1,1), (1,1,1) decode to Lost and are **not** in enc's image. The
canonical set is therefore

  C = {(x,z,l) : l ⇒ x = z = 0},  5 of the 8 triples,

and enc : Sⁿ → Cⁿ is a bijection. Two consequences: (i) "compare/hash the raw
(x,z,l) blob" equals "compare/hash the logical word" **on C and only on C** —
`PartialEq` (`data.rs:541`) and `structural_hash_lossy` (`hash.rs:33`) read all
three planes, so canonicality is a genuine *precondition* of the hash-keyed sum,
not a representation detail; this is exactly G-011/G-014's uniqueness claim, true
with the hypothesis and false without it. (ii) Any mutator must map C → C, else it
produces a bit pattern that denotes no lossy word at all.

*What the gates must do.* The crate's model — "a gate touching a lost qubit does
not happen" — makes the correct output the input, bit-exact on all three planes.
Check the CNOT columns: the unguarded rules are x_t ⊕= x_c and z_c ⊕= z_t. With c
lost, canonicality gives x_c = z_c = 0 so the first column is automatically
harmless, but the second sets z_c := z_t, producing (0,1,1) at c — outside C.
Symmetrically with t lost, x_t ⊕= x_c breaks C. CZ writes z_a ⊕= x_b and
z_b ⊕= x_a and breaks C at whichever operand is lost. So **a guard on `xor_z_col`
is mathematically required, independently of the "skip" modelling choice**: no
semantics can be right that writes a Z bit underneath a loss bit, because that
triple denotes nothing. Guarding both operands in both primitives, as shipped, is
correct.

*Observability of a non-canonical word.* Take q lost with a stray x bit, (1,0,1).
Then `get(q)` = Lost, `Display` = 'L', `weight()`, `loss_weight()` and `is_lost`
are all blind — but `PartialEq` says unequal to the canonical twin, `key_hash()`
differs, `pauli_code(q)` returns 1 instead of 0, and — the sharp one —
`clear_loss(q)`/`loss_cleared(q)`, specified as "return the site to identity I"
and used as the branch key of the loss-recovery arm of `LossChannel`, yield
Present(X). So a non-canonical word both splits a hash-keyed sum into two entries
for one logical operator **and** turns into a wrong operator with a real
coefficient the moment the loss channel recovers the atom. The invariant is
load-bearing, not cosmetic.

*Injectivity (G-016).* No `SymplecticColumns` primitive writes `lbits` (only
`set_lost`/`clear_loss`/`set`/`with_lost`/`loss_cleared`/`set_*_bit` ever touch
it), so the loss mask is invariant and the gate map is sector-wise: on each fixed
loss mask it is the Sp(2n,2) element restricted to the present sub-block and the
identity on lost sites — a bijection — and distinct sectors have disjoint images
because the loss mask is part of the key. So the gate map **is** injective on
canonical lossy keys and the re-keying is collision-free: G-016's Rust-side claim
holds. The Lean simply cannot state it, because `Symplectic.lean:314` makes `lost`
an external `variable (lost : Fin n → Prop)` shared by both operands rather than a
component of the state. Note the Rust oracle is correctly scoped to that model:
`lossy_clifford_generators_preserve_symplectic_form` draws v and w with a *shared*
loss mask, which matters — with independent masks ω is genuinely not preserved
(v = "XI", w = "LZ", `cnot(0,1)`: x_v(1) flips 0 → 1, since a skipped w is not a
fixed point of the unguarded map). The test is sound and faithful to what it
states; what is missing is only the loss-plane-as-state extension.

**Correct answer.** C = {(x,z,l) : l ⇒ x = z = 0}; blob equality/hashing = logical
equality/hashing only on C; a gate touching a lost qubit must be the identity on
**all three** planes, so both `xor_x_col` and `xor_z_col` need the guard exactly as
shipped; every branch builder must map C → C or carry an explicit `¬lost`
precondition.

**Which implementation is right.** The shipped `-2` gate/channel code satisfies
all of this on every reachable path. The only holes are the public branch
builders, and there **two shipped implementations disagree on a lost site**:
`LossyPauliWord::toggled_bits2` (`data.rs:445`, the override) keeps the loss bit
and toggles X → (1,0,1), decoded as Lost; the trait-default `into_toggled_bits2`
(`ppvm-traits-2/src/word.rs:307`, *not* overridden by the lossy word) routes
through `set_x_bit`, which clears loss → (1,0,0), decoded as Present(X).
**Neither is right.** The override leaves C entirely (its output denotes no lossy
word); the default silently resurrects a lost atom as an X operator, which the
"a lost atom does not participate" model forbids. The correct specification is a
precondition `¬is_lost(i)` — the callers' `is_lost` guard already means "emit no
branch here". If a total function were wanted, the only choice consistent with the
gate model is loss-wins identity (return the word unchanged), which is neither of
the two. Both are currently unreachable: every in-crate call site
(`rotation.rs` 165/240/262/336/392/424, `store.rs` 380/396/1421/1443, `proj.rs`
76/93, `noise.rs` 316, `clifford.rs` 253/288/443/500 via `cy_toggles`,
`column_store/rotations/{mod,rx}.rs`) tests `is_lost` first, and
`PREFER_BORROWED_REKEY = true` for the lossy word keeps the owned path out of the
lossy Clifford entirely.

**Live defect: NO — mislabelling and an undocumented public precondition only.**
Word-level Clifford is canonicality-preserving on all 225 (2-qubit word,
generator) pairs and is a bit-exact no-op on the raw planes whenever it touches a
lost qubit; every sum-level gate/rotation/projection/noise/loss kernel and both
columnar rotation kernels guard on `is_lost` before calling a branch builder. IF a
downstream caller calls the public `PauliBits::toggled_bits` on a lost site the
consequences are real and measured: one logical operator occupies two entries in a
`LossyPauliSum` (`from_terms([LZ, LZ]) = [("LZ",1.0),("LZ",1.0)]`), and the
loss-recovery branch then emits the wrong operator with a real weight
(`loss_channel(0,0.25)` gives `("XZ",0.25)` where the canonical word gives
`("IZ",0.25)`).

**G-013's mutation claim is false.** Deleting the three-line
`if self.is_lost(ctrl) || self.is_lost(tgt) { return; }` from `xor_z_col`
(`clifford.rs:117`) does **not** survive the sector. The stray z bit is not
confined behind `Display`: `xor_z_col` is `z_ctrl ^= z_tgt`, so a later CNOT whose
*target* is the corrupted lost qubit XORs that stray bit onto a **present**
control's z bit, which `Display` shows.
`crates/ppvm-conformance-2/tests/lossy_pauli_word_diff.rs:248
clifford_replay_matches_old` fails at seed 1, n = 16 (output below). What survives
the mutation is `lossy_pauli_word_lean.rs` — measured by the second agent:
`clifford_leaves_loss_mask_invariant` reads only the loss plane, and
`lossy_clifford_generators_preserve_symplectic_form` uses a shared loss mask under
which ω is unchanged (x = 0 at the shared lost site on both operands, so the
corrupted z bit multiplies zero: 1 vs base 1 in every case), and the crate's own
table (`clifford.rs:186-199`) has no L-control-with-Z-target row. **So what
remains of G-013 is a low-severity coverage nit — no oracle asserts canonicality
on the BITS after propagation — not a high-severity bridge gap.**

**Proposed fix** — none of it changes release-mode observable behaviour.

1. `crates/ppvm-lossy-pauli-word-2/src/data.rs:421` (`fn toggled_bits`) — after
   the existing bounds `debug_assert!`, add
   `debug_assert!(!self.is_lost(i), "toggled_bits on lost site {i} breaks the canonical loss invariant");`
   and at `data.rs:445` (`fn toggled_bits2`) the same for both `i` and `j`.
2. **Widened by the second opinion:** the identical hole exists in
   `LossyPauliKeyColumn::toggled_bits`/`toggled_bits2`
   (`crates/ppvm-lossy-pauli-word-2/src/column.rs`, the two
   `LossyPauliWord::from_planes(..., BitArray::new(self.lplanes[row]), ...)`
   returns) and in the `KeyColumn` trait defaults
   (`crates/ppvm-traits-2/src/batch.rs:135`/`:145`, which delegate to
   `self.get(row).toggled_bits(...)` and inherit the lossy override). That is the
   surface the columnar rotation kernels use — i.e. the surface most likely to gain
   new gate kernels — so the guard rail must cover it or it stops at the scalar
   word.
3. `crates/ppvm-lossy-pauli-word-2/src/data.rs:38-44` — the type docstring asserts
   "Every mutator upholds it: …", which is **false** for the two branch builders.
   Name them as the exception with the precondition `¬is_lost(i)`.
4. `crates/ppvm-traits-2/src/word.rs:235` (and `:268`, `:307`) — document the
   precondition: "`i` (and `j`) must not be lost; the lossy override and this
   default differ on a lost site (loss-preserving vs loss-clearing), so callers
   must guard on `is_lost`."
5. Test closure, behaviour-preserving: in
   `crates/ppvm-conformance-2/tests/lossy_pauli_word_lean.rs`, call the existing
   `assert_site_invariants` on the POST-circuit word inside
   `clifford_leaves_loss_mask_invariant` (it currently compares only the loss
   plane), and add exhaustive 2-qubit coverage: for all 25 lossy words ×
   {h, s, cnot(0,1), cnot(1,0), cz, cy}, assert canonicality on the bits and that a
   gate touching a lost qubit is the identity on all three planes. That is the test
   that fails at gate 1 under the `xor_z_col` mutation, instead of relying on the
   leak-into-`Display` accident.

**Not proposed:** any change to the guards, to `into_toggled_bits2`, or to the loss
semantics. Making the two builders agree changes observable behaviour on an input
that is currently unreachable, and is the user's call.

**Evidence.** Probe `crates/ppvm-conformance-2/tests/zz_adj_lossy.rs` (unmutated
tree; deleted after the run), verbatim:

```
running 3 tests
canonical : [0: x=0 z=0 l=1][1: x=0 z=1 l=0] display=LZ
noncanon  : [0: x=1 z=0 l=1][1: x=0 z=1 l=0] display=LZ
canonical(bad) = false
toggled_bits2      -> [0: x=1 z=0 l=1][1: x=0 z=1 l=0] display=LZ canonical=false
into_toggled_bits2 -> [0: x=1 z=0 l=0][1: x=0 z=1 l=0] display=XZ canonical=true
eq=false same_hash=false
get(0): toggled=Lost into_toggled=Present(X)
display eq=true get eq=true weight 2==2 loss_weight 1==1 is_lost true==true
Eq: false | key_hash canon=0xf2703b327e0c3932 bad=0x8f3d55be65ef69de | pauli_code canon=0 bad=1
loss_cleared(0): canon -> IZ | bad -> XZ
clear_loss(0): canon -> IZ | bad -> XZ
test probe_3_toggled_vs_into_toggled_diverge_on_a_lost_site ... sum from_terms([LZ_canon, LZ_noncanon]) = [("LZ", 1.0), ("LZ", 1.0)]
ok
loss_channel(0, .25) on the noncanonical LZ term = [("XZ", 0.25), ("LZ", 1.0)]
loss_channel(0, .25) on the canonical    LZ term = [("LZ", 1.0), ("IZ", 0.25)]
reset_loss_channel(0) on noncanonical LZ = [("LZ", 0.0)]
test probe_2_toggled_bits_breaks_canonicality_and_is_observable ... ok
probe1: 225 (word, gate) pairs all canonical; LZ under cnot(0,1) -> [0: x=0 z=0 l=1][1: x=0 z=1 l=0]
test probe_1_word_level_clifford_keeps_canonical ... ok

test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
```

(probe_1 also asserted, for every lossy 2-qubit word and every generator touching
a lost qubit, bit-exact equality of all three planes with the input — 225 pairs,
no failure.)

Mutation probe (G-013's exact mutation, applied and then reverted with
`git checkout --`; both guards verified present at `clifford.rs:102` and `:117`,
tracked tree clean afterwards). `cargo test -p ppvm-lossy-pauli-word-2` →
`test result: ok. 28 passed` (the crate's own unit tests do not catch it, as the
ledger says). Conformance, verbatim:

```
     Running tests/lossy_pauli_word_diff.rs (target/debug/deps/lossy_pauli_word_diff-186aef8a20b1608b)
running 10 tests
...
test clifford_replay_matches_old ... FAILED
...
---- clifford_replay_matches_old stdout ----
thread 'clifford_replay_matches_old' (41215878) panicked at crates/ppvm-conformance-2/tests/lossy_pauli_word_diff.rs:248:13:
assertion `left == right` failed: H/S/CNOT YLYIIIZIYLLZXLXY seed 1
  left: "XLXXYZYXXLLXILXX"
 right: "XLXXYIYXXLLXILXX"
failures:
    clifford_replay_matches_old
test result: FAILED. 9 passed; 1 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.09s
```

The second agent, instructed not to mutate `crates/`, reproduced the same
divergence by modelling the mutated semantics outside `crates/` (validating the
model by checking that the *guarded* model reproduces the shipped word's `Display`
on all 48 replay circuits) and obtained the identical first divergence — seed 1,
n = 16, `YLYIIIZIYLLZXLXY` → mutated `XLXXYZYXXLLXILXX` vs correct
`XLXXYIYXXLLXILXX` — plus a hand case: "LZ", `cnot(0,1)` then `cnot(1,0)` →
mutated "LI" vs correct "LZ". Its own canonicality sweep was wider than the
first's: 275 (word, generator) pairs over the full 25-word alphabet ×
{h0, h1, s0, s_dag1, sqrt_x0, sqrt_y1, cnot(0,1), cnot(1,0), cz, cy(0,1), cy(1,0)},
all in C with the loss plane bit-identical, plus 48 random 200-gate circuits
(n = 2, 3, 5, 16 × 12 seeds) checked on the bits.

**Second opinion — corroborated, with four corrections.**

1. **Fix scope is incomplete** (`column.rs` and `batch.rs:135/145`) — folded into
   item 2 of the proposed fix above. This is the correction that matters.
2. **"No consumer" is imprecise.**
   `crates/ppvm-vihaco/src/component/backend/traits_2.rs:43-66` defines
   `Lossy64` … `Lossy2048` as
   `Sum<HashMapStore<LossyPauliWord<[u8; N]>, f64>, …>` under the `traits-2`
   feature, so there **is** a downstream consumer. It calls no bit-level API (a
   grep for `toggled_bits`/`set_*_bit`/`from_planes` across `ppvm-vihaco`,
   `ppvm-sym-2` and `ppvm-tableau-2` finds no lossy uses), so the conclusion is
   unchanged, but the exposure framing understates it.
3. **Math nit on "minimal".** Guarding both operands in both primitives is correct
   but not minimal: under canonical input a lost site reads as 0 in both planes, so
   each XOR column only needs the guard on the operand it *writes* (`xor_x_col` on
   tgt, `xor_z_col` on ctrl), and the H/S guards are entirely redundant for
   canonicality (swapping 0,0 and z ⊕= 0). Guarding both is the right conservative
   choice and matches "the whole gate skips" literally — but it means a mutation
   deleting the H/S guards or the read-side operand guards would survive every
   test **and be harmless**, which is worth knowing before anyone writes a
   mutation-coverage gate.
4. The first agent's "reasoned, not measured" sub-claim about the Lean oracle's
   blindness is now measured (ω unchanged, 1 vs base 1, for LZ/LX, LY/LZ, ZL/XL
   under `cnot(0,1)` and `cnot(1,0)`).

One near miss the second agent flagged for the debug_assert:
`crates/ppvm-conformance-2/benches/word_surface/lossy_branch.rs:18`/`:40` call
`toggled_bits`/`toggled_bits2` on a lossy word at SITE = 127 and SITE2 = 191;
those are safe only because the "IXYZL" cycle puts 'Y' at 127 (127 % 5 = 2) and
'X' at 191 (191 % 5 = 1). If the WIDTH/SITE constants ever become a multiple of 5,
any debug bench build panics. Pin the sites explicitly or comment the constraint.
The assert must also be validated with `cargo test --workspace --all-targets` in
debug, not by grep — neither agent ran the suite with the assert in place.

**Contested inside the fix.** None on the ruling. The judgement call the user owns
is that item 3 softens a currently-false absolute claim in the type docstring
("Every mutator upholds it"), which is a real contract change for downstream
callers, and that the fix documents rather than *enforces* the precondition in
release builds.

**Sign-off needed.** No release-mode behaviour changes. The docstring
weakening (item 3) is a published-contract change and should be acknowledged.

### U4 mixture weights and sampler

Rows: **G-018** (`adjudicated-defect`), **G-019** (`adjudicated-defect`),
**G-020** (`adjudicated-defect`), **G-071** (`adjudicated-spec`).
Verdict: **code is wrong** — three distinct wrong numbers reach any caller of the
public mixture API with `sum_cutoff > 0`. Second opinion: **corroborated**
(independent derivation and probes, numbers agreeing within Monte-Carlo noise),
with three corrections, one of which **inverts the sign of one proposed fix**.

**Question.** When `for_each_z_branch` discards a measurement child because its
conditional Born probability is `<= sum_cutoff`, must the mixture renormalize (and
must `MixtureSampler` draw from a distribution summing to 1), or is silently
leaving the weights sub-stochastic and letting the sampler's `.min(len-1)` clamp
dump the whole deficit on the smallest-weight entry defensible?

**Derivation.** (Kraus/instrument formalism only; legacy, Lean and the passing
tests were not used as an oracle.)

*What a mixture is.* ρ = Σ p_i ρ_i is a state iff p_i ≥ 0 and Σ p_i = 1
(Tr ρ = 1). A weight vector summing to S < 1 is not a density operator; it is the
output of a trace-decreasing map. `mixture/data.rs:19` documents the object as "A
classical probability distribution over complete generalized-tableau states", so
Σ = 1 is the stated invariant.

*Z measurement as an instrument.* Π_b = (I + (−1)^b Z_q)/2 with Π_b†Π_b = Π_b,
Π_bΠ_b' = δ_{bb'}Π_b and **Π₀ + Π₁ = I** (that completeness relation *is* trace
preservation). Unravelling ρ = Σ p_i |ψ_i⟩⟨ψ_i| gives children
weight(i,b) = p_i q_i(b) with q_i(b) = ⟨ψ_i|Π_b|ψ_i⟩ and state
Π_b ρ_i Π_b / q_i(b); since Σ_b q_i(b) = Tr(ρ_i(Π₀+Π₁)) = 1 we get
Σ_{i,b} weight(i,b) = Σ_i p_i = 1. **The correct unravelling maps a distribution
to a distribution: each child is parent × conditional Born probability, and the two
children of a parent sum to the parent** — G-018's claim, true of the exact map.
The reported marginals are P(b) = Σ_i p_i q_i(b), fully determined *before* any
retention decision, so a `measure()` that returns P(1) = 0 because it declined to
keep the b = 1 conditional state is reporting a **false number**, not a truncated
one.

*Truncation policies.* If a set D with mass m is discarded, the defensible options
are (a) renormalize, (b) keep the deficit explicit and report it, (c) refuse. For
(a): ρ − ρ' = Σ_D w ρ − (m/(1−m)) Σ_keep w ρ, so
‖ρ − ρ'‖₁ ≤ m + (m/(1−m))(1−m) = **2m**, i.e. trace distance ≤ m and
|⟨O⟩ − ⟨O⟩'| ≤ 2m for ‖O‖_∞ ≤ 1. (Both agents independently derived 2m; **G-071's
proposed 2m/(1−m) is valid but loose** — see the closure edit there.) Silently
leaving the weights sub-stochastic while continuing to *use* them as a distribution
is none of (a)/(b)/(c).

*The sampler's inverse-CDF lemma.* With c_i = Σ_{j≤i} p_j and U ~ U[0,1),
I = min{i : c_i > U} = `partition_point(|&b| b <= U)` satisfies
P(I = i) = c_i − c_{i−1} = p_i **only if c_{n−1} = 1**, so that the intervals tile
[0,1). If c_{n−1} = S < 1 then U ∈ [S,1) (probability 1−S) makes
`partition_point` return n, out of range, and `.min(n−1)` gives
P(I = n−1) = p_{n−1} + (1−S). The clamp is therefore mathematically equivalent to
sampling from **p + (1−S)·δ_{last}**. Because `sampler()` sorts entries in
*descending* weight before accumulating (`sampler.rs:106`), "last" is the
**smallest-weight** entry: relative error (1−S)/p_last, the worst possible
recipient. With `sum_cutoff = 1e-7` (the value this repo's own
`docs/notebooks/msd.py` uses) a retained entry can legitimately weigh ~1e-7 while
an accumulated deficit of ~1e-3 gives a factor ~1e4 over-sampling of that branch.

*Code check.* `crates/ppvm-tableau-2/src/mixture/measure.rs:50` and `:86` both read
`if p_other > self.sum_cutoff { …push branch… }` with **no `else`**. When the test
fails the child is never projected, never handed to `visit` (so its mass is absent
from `measure()`'s returned `(zero, one, lost)` triple), and never reaches
`insert_branches` (so the `dropped` flag that would trigger
`normalize_probabilities()` at `:102-104` stays false), while the parent's weight is
still multiplied by `p_likely` (`:62`/`:93`). `truncate()` renormalizes only
`if self.entries.len() != before`. Result: sub-stochastic weights **and** a false
reported Born probability. Per-call leak bound (second agent):
Σ_i p_i·p_other,i ≤ sum_cutoff·Σ_i p_i ≤ sum_cutoff, so mass ≥ (1−sum_cutoff)^k
after k measurements and the last entry's relative over-sampling is bounded by k —
bounded per call, unbounded in circuit depth. By contrast the noise paths
(`mixture/noise/loss.rs:41,98`; `pauli.rs:99,174`) drop only *inside* `insert_*`,
whose return value does trigger `normalize_probabilities()`, so the leak is
specific to `for_each_z_branch`'s pre-insert gate. Units mismatch, additionally:
`p_other` is a *conditional* probability while `sum_cutoff` is compared against
*absolute* weight in `truncate` (`data.rs:141`) and `insert_branches`
(`data.rs:186`); the mass actually removed is `other_probability = parent·p_other`.

*G-071 from the same algebra.* `truncate()` = drop `probability <= sum_cutoff`,
then renormalize iff something was removed. Given a stochastic input that is
exactly option (a) with the 2m bound, i.e. **correct** (the bound being unstated is
a documentation gap, not a numeric error). G-071's second complaint — that a
sub-cutoff branch colliding with a live entry is accumulated regardless
(`data.rs:186`) — is mathematically **correct, not a bug**: adding mass to an
existing entry is exact and lossless, and the cutoff properly applies to the
resulting entry's total, which the following `truncate()` re-tests. That sub-row is
a spec-characterization complaint, not a defect.

*G-020 axioms.* `structurally_equal` (`equality.rs:19`) = equal `is_lost` + equal
tableau rows + `amplitudes_equal`, where the latter is equal coefficient count and,
for every key of `left`, |left_k − right_k|² **<** c² with c =
`left.coefficient_threshold`.
**Reflexivity fails at c = 0** — the test is strict, |0|² < 0 is false — so at
`coefficient_threshold = 0` no two states are ever equal, not even a state and
itself, and dedup is dead (measured).
**Symmetry**: c and the iterated key set come from `left` only, but a left-only key
would need |left_k| < c, contradicting the amplitude pruner's invariant
`norm_sqr() > cutoff_sq` (`data.rs:752`/`:1241`), so within one mixture (shared
threshold) the asymmetry is inert; G-020's item (3) is **not** proven unsound.
**Transitivity fails unconditionally**: the relation is a per-index sup-norm ball,
so a, a + 0.6c, a + 1.2c with identical tableau rows/phases/loss (hence identical
`fingerprint`, which hashes only x/z words, phases and loss) land in the same
bucket and really are compared: A~B, B~C, A≁C. `insert_branches` takes the *first*
bucket match (`data.rs:177-185`/`:188`/`:194`), so which entries coalesce — and
hence the per-entry weights — is insertion-order dependent. `mixture/mod.rs:5-8`
already says this out loud. Merge error when collapsing a, b onto a:
≤ 2·min(p_a,p_b)·‖a−b‖₂ ≤ 2·min(p)·c·√D for stored support size D — bounded, but
D-dependent and unstated, so the tolerance is not an error budget in any norm the
readout uses. Item (4) (refusing to merge states equal up to a global phase) is
conservative: it costs entries and time, it does not corrupt weights (probe F:
three rounds of h+measure give 6 entries where the physical mixture has 2).

**Correct answer.** The mixture weights must sum to 1; `measure()` must report
P(b) = Σ_i p_i q_i(b). If children are discarded with mass m the only defensible
policies are renormalize by 1/(1−m) (trace distance ≤ m, |Δ⟨O⟩| ≤ 2m), keep the
deficit explicit and report it, or refuse. The sampler must draw from a cumulative
vector whose last element is exactly 1. `structurally_equal` must be an equivalence
relation to be a sound merge key; the strict `<` makes it non-reflexive at
threshold 0 and the sup-norm ball makes it non-transitive at any positive
threshold.

**Which implementation is right.** There is one mixture implementation. Within it,
the noise paths (`mixture/noise/loss.rs:41,98`, `pauli.rs:99,174`) are **right**;
`truncate()` (`data.rs:140-150`) is **right**; the measurement path
`for_each_z_branch` (`measure.rs:50`, `:86`) is **wrong** — its pre-insert gate
discards the child outside the flag-and-renormalize discipline the rest of the file
follows, and discards it before `visit`, so the mass is missing from the reported
triple as well as from the weights; and `MixtureSampler::choice`
(`sampler.rs:39-42`) is **wrong as a defensive layer** — instead of enforcing or
restoring the sum-to-1 precondition it converts a violated precondition into a
silently different distribution, with the descending sort making the
smallest-weight entry the recipient. (Separately: the same clamp turns an empty
mixture, reachable via `new(n, thr, 1.0)`, into an out-of-bounds panic at
`sampler.rs:52` rather than a diagnosable error.)

**Live defect: YES — three distinct wrong numbers** reach any caller of the public
`ppvm-tableau-2` mixture API with any `sum_cutoff > 0`:

1. `GeneralizedTableauMixture::measure` returns a **false analytic Born
   probability**: on a single-entry mixture prepared with RY so that the true
   (p0,p1) = (0.97, 0.03), at `sum_cutoff = 0.05` it returns
   `(0.96999999999999997, 0.0, 0.0)` — it reports p(1) = 0 when p(1) = 0.03, and
   the triple sums to 0.97.
2. The stored weights become sub-stochastic and decay multiplicatively with no
   bound and no renormalization: six such measurements give total mass
   `0.83297200492900003` = 0.97⁶.
3. `MixtureSampler::sample_shots` then draws from the wrong distribution, with the
   whole deficit on the smallest-weight branch. In an exchange-symmetric two-branch
   mixture reached entirely through the public gate/measure API
   (`h(0); measure(0); RY(1); measure(1)`) the two entries have **identical stored
   weight 0.485** yet 400k shots come out 0.485010 / 0.514990 — a +0.02998
   absolute, +6.18% relative excess on the last entry, and a physically impossible
   asymmetry between two exchange-symmetric branches. With an explicit deficit of
   0.1 over weights [0.5, 0.3, 0.1] the 0.1 entry is sampled at 0.200745 —
   **+100.7%** relative error on that branch's shot count.

G-020's non-reflexivity is live but is a **size/performance** defect, not a wrong
number: at `coefficient_threshold = 0` dedup never fires and the entry count
doubles per measurement round (2 → 4 → 8) instead of merging (2 → 4 → 6 at 1e-12);
weights and mass stay correct (mass = 1.0 in both runs).

**Non-monotonicity in `sum_cutoff`** (second agent, and the sharpest statement of
the bug): the *same* circuit is exactly right at `sum_cutoff = 0.02` — measure(1)
returns (0.97, 0.03), weights [0.5, 0.5], mass 1.0 — and wrong at 0.04/0.05. At
0.02, `p_other = 0.03` passes the pre-insert gate, the child is pushed,
`insert_branches` then drops it (absolute 0.015 ≤ 0.02) and sets `dropped`, so
`normalize_probabilities()` runs. That is evidence of an accident, not a policy,
and it localizes the bug precisely to the pre-insert conditional gate.

**There is no "deliberately sub-normalized" convention to appeal to.** The second
agent closed the first's own escape hatch: legacy
`ppvm-tableau-sum/src/data.rs:136` asserts
`debug_assert!(p_cum.last() >= 1 - sum_cutoff, "Normalization error in sum")`, and
running the legacy mixture through the same six measurements makes that assert
**fire** at mass 0.83297200492899992. Legacy carries the identical leaky gate
(`ppvm-tableau-sum/src/measure.rs:131`/`:221`, comment "intentionally stricter than
normal truncate"), so the `-2` crate is a faithful port — but the original authors'
own assertion says the sum must stay ≥ 1 − sum_cutoff, and the `-2` rewrite dropped
the guard, so today the violation is not even detectable in debug builds.

**Proposed fix.**

- **FIRST (the row-critical one)**, `crates/ppvm-tableau-2/src/mixture/measure.rs`,
  `for_each_z_branch`: bring the pre-insert drop under the flag-and-renormalize
  discipline the noise paths use — add `let mut dropped_below_cutoff = false;`
  next to `let mut keys_changed = false;` (line ~26), add
  `else { dropped_below_cutoff = true; }` to both `if p_other > self.sum_cutoff {`
  arms (`:50` and `:86`), and change `:104` to
  `if self.insert_branches(branches) | dropped_below_cutoff {` (bitwise `|`, so
  `insert_branches` still runs; `let inserted = …; if inserted || dropped` reads
  better). Weights sum to 1−m and `normalize_probabilities()` divides by that sum,
  preserving ratios — option (a), trace distance ≤ m. No division by zero:
  p_likely ≥ 0.5 always. Both callers (`measure()` at `:124` and `Reset::reset` at
  `:143`) are correct under it. **CHANGES OBSERVABLE BEHAVIOUR.**
- **SECOND** (independent, needed for the reported triple to be correct): count the
  dropped child's mass in the returned `(zero, one, lost)`. `visit` cannot simply be
  handed the un-projected fork, because `Reset`'s closure *mutates* the tableau it
  receives (`:143-150` applies `x(qubit)` on outcome `Some(true)`), so the callback
  must be split into a probability-reporting arm and a state-visiting arm.
  Observable: `measure()` starts returning 0.03 instead of 0.0.
- **THIRD**, `crates/ppvm-tableau-2/src/mixture/sampler.rs`: in `sampler()` (line
  ~110) add
  `debug_assert!((sum - 1.0).abs() <= 1e-9, "mixture weights are not normalized: {sum}")`
  — literally the legacy guard, and the best first move since it changes no release
  behaviour. Keep `.min(len-1)` at `:42` purely as a float-round-off guard, and make
  the empty case an explicit `assert!(!self.entries.is_empty())` rather than an
  out-of-bounds panic.
- **FOURTH (units)**: **DO NOT APPLY AS DESCRIBED — see Contested below.**
- **FIFTH**, `crates/ppvm-tableau-2/src/mixture/equality.rs:20`: change
  `delta.norm_sqr() < cutoff_sq` to `<=`. Restores reflexivity at every threshold
  ≥ 0 and degenerates to exact equality at 0. Callers are `structurally_equal`
  (`data.rs:188`/`:194`) and `structurally_equal_mutated` (`:229`/`:237`), i.e. the
  merge key for both the branch and the lazy/noise insertion paths. Observable only
  at threshold 0 and at exact-boundary ties, where it turns a dead dedup into a
  working exact-equality dedup: entry counts drop (8 → 6) and merged weights become
  one summed entry, with identical state and mass (1.0 both ways). Non-transitivity
  is **inherent to a tolerance ball** and should be closed by documenting the order
  dependence and the 2·min(p)·c·√D merge bound on `GeneralizedTableauMixture`, not
  by a code change.
- **Not a fix:** G-071's "insertion-time cutoff bypass" sub-row needs no change.

**Evidence.** `cargo test -p ppvm-conformance-2 --test zz_adj_mix -- --nocapture
--test-threads=1` (debug build; scratch file deleted). Verbatim:

```
running 5 tests
test probe_a_single_entry_subcutoff_branch_drop ... [A] entries before measure = 1
[A] mass before measure    = 1.00000000000000000
[A] measure(0) -> (zero, one, lost) = (0.96999999999999997, 0.00000000000000000, 0.00000000000000000)
[A] returned triple sum    = 0.96999999999999997
[A] true (zero, one)       = (0.97000000000000000, 0.03000000000000000)
[A] entries after measure  = 1
[A] mass after measure     = 0.96999999999999997
[A] mass deficit           = 0.03000000000000003
[A] sampler entry-weight sum = 0.96999999999999997
ok
test probe_b_repeated_measurement_multiplicative_mass_decay ... [B] after measure(0): triple=(0.970000000000,0.000000000000) sum=0.970000000000 entries=1 mass=0.96999999999999997
[B] after measure(1): triple=(0.940900000000,0.000000000000) sum=0.940900000000 entries=1 mass=0.94090000000000007
[B] after measure(2): triple=(0.912673000000,0.000000000000) sum=0.912673000000 entries=1 mass=0.91267300000000007
[B] after measure(3): triple=(0.885292810000,0.000000000000) sum=0.885292810000 entries=1 mass=0.88529281000000004
[B] after measure(4): triple=(0.858734025700,0.000000000000) sum=0.858734025700 entries=1 mass=0.85873402570000001
[B] after measure(5): triple=(0.832972004929,0.000000000000) sum=0.832972004929 entries=1 mass=0.83297200492900003
[B] 0.97^6 = 0.83297200492899992
ok
test probe_c_sampler_clamp_assigns_deficit_to_last_entry ... [C] measure(0) = (0.500000000000, 0.500000000000), entries=2
[C] measure(1) = (0.970000000000, 0.000000000000) sum=0.970000000000, entries=2
[C] weights = [0.485, 0.485]
[C] mass    = 0.96999999999999997
[C] sampler weights (desc) = [0.485, 0.485]  total = 0.96999999999999997
[C] shots = 400000
[C]   pattern 00: 194004 shots, empirical freq = 0.485010
[C]   pattern 10: 205996 shots, empirical freq = 0.514990
[C] stored weights normalized (the mathematically correct target) = [0.5, 0.5]
[C] stored weights as-is (the distribution the sampler claims)    = [0.485, 0.485]
ok
test probe_d_synthetic_three_entry_deficit ... [D] weights = [0.5, 0.3, 0.1], total = 0.90000000000000002, deficit = 0.09999999999999998
[D] sampler weights (desc) = [0.5, 0.3, 0.1]
[D]   pattern 00: 199683 shots, empirical freq = 0.499208
[D]   pattern 01: 80298 shots, empirical freq = 0.200745
[D]   pattern 10: 120019 shots, empirical freq = 0.300048
[D] as-stored weights   : 0.5 / 0.3 / 0.1
[D] renormalized target : 0.555556 / 0.333333 / 0.111111
[D] clamp prediction    : 0.5 / 0.3 / 0.1+0.1 = 0.2 for the LAST (smallest) entry
ok
test probe_e_truncate_renormalizes ... [E] before truncate: [0.6, 0.3, 0.1] mass=0.99999999999999989
[E] after  truncate: [0.6666666666666667, 0.33333333333333337] mass=1.00000000000000000
[E] sub-stochastic + nothing dropped: [0.6, 0.3] mass=0.89999999999999991
ok

test result: ok. 5 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.36s
```

```
running 1 test
test probe_f_reflexivity_at_zero_coefficient_threshold ... [F] thr=1e-12: after h+measure entries=2 weights=[0.5, 0.5]
[F] thr=1e-12: after h+measure again entries=4 weights=[0.25, 0.25, 0.25, 0.25]
[F] thr=1e-12: after 3rd round entries=6 mass=1.00000000000000000
[F] thr=0e0: after h+measure entries=2 weights=[0.5, 0.5]
[F] thr=0e0: after h+measure again entries=4 weights=[0.25, 0.25, 0.25, 0.25]
[F] thr=0e0: after 3rd round entries=8 mass=1.00000000000000000
ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 5 filtered out; finished in 0.00s
```

Reading: [A] the mass leak *and* a literally false reported p(1) = 0.0 where the
truth is 0.03; [B] exactly multiplicative decay (0.97⁶ to 14 digits); [C] the
headline, reached entirely through the public API — two entries with identical
stored weight 0.485 sampled at 0.485010 vs 0.514990; [D] the clamp exactly as the
inverse-CDF algebra predicts, with the other two entries untouched (0.499208 vs
0.5, 0.300048 vs 0.3), which rules out "the sampler renormalized"; [E] `truncate()`
renormalizes when it removes something and not otherwise, so whether the leak is
repaired depends on whether truncate happened to fire — an inconsistency, not a
policy; [F] non-reflexivity at threshold 0, mass still 1.0.

The second agent's independent probe (`zz_2nd_zz_adj_mix.rs`, deleted) reproduced
the identical `(0.96999999999999997, 0.0, 0.0)` triple and mass 0.97⁶ =
`0.83297200492899992`, got 400k shots at 0.486085 / 0.513915 on the
exchange-symmetric pair (+0.028915 absolute, +5.96% relative; algebra predicts
+0.03 / +6.19%, MC sd 0.0008) and 0.501185 / 0.300745 / **0.198070** on
[0.5, 0.3, 0.1], and reproduced the empty-mixture panic (`sum_cutoff = 1.0` leaves
0 entries; `sample()` panics "index out of bounds: the len is 0 but the index is 0"
at `sampler.rs:52`).

**Second opinion — corroborated, with three corrections.**

1. **The FOURTH proposed fix has the sign backwards, and alone it makes things
   worse.** The first agent proposed testing `other_probability > sum_cutoff`
   instead of `p_other > sum_cutoff` and wrote "Observable: yes — retains strictly
   more branches." Since `other_probability = parent · p_other ≤ p_other`, that
   test is strictly **harder** to pass: it retains strictly **fewer** branches. Its
   stated motivation is already satisfied (parent ≤ 1, so the absolute mass removed
   per drop is already ≤ sum_cutoff). Worse, probe 9 shows the current conditional
   gate sometimes accidentally *saves* the invariant (the `sum_cutoff = 0.02` case
   above), so switching to the absolute test moves exactly those drops back into
   the leaky pre-insert path: applied **without** FIRST it strictly increases the
   leak. It must be sequenced after FIRST, or dropped.
2. **The non-monotonicity framing** (right at 0.02, wrong at 0.04/0.05) is missing
   from the first write-up and is the sharpest statement of the bug; recorded above.
3. **The "maybe it's a deliberate sub-normalized convention" escape hatch is
   closed**: the legacy `debug_assert` fires on this exact path (recorded above).
   The first agent noted the assert exists and that `-2` dropped it, but did not
   test that the path violates it.

Minor: the first agent's "unbounded as p_last → sum_cutoff" is true only
cumulatively — per `measure()` call the leak is ≤ sum_cutoff, so mass ≥
(1−sum_cutoff)^k and the last entry's relative over-sampling is bounded by the
number of measurements k (unbounded in depth, not in one step); the 1e4 figure for
cutoff 1e-7 is consistent with that. The 2m bound, the "G-071 insertion sub-row is
not a defect" reading, and the whole G-020 analysis (including "symmetry is inert
given the pruner's |a| > c invariant", verified at `data.rs:752`/`:1241`) were
reached independently by both agents.

**Contested inside the fix.** Two items, neither resolved by consensus:

- **FOURTH (units).** Contested and **not** to be applied as written: the first
  agent claims it retains more branches, the second measured that it retains fewer
  and that applied alone it increases the leak. The second agent's argument
  (`other_probability ≤ p_other`) is arithmetic, so the safe reading is: sequence
  it after FIRST or drop it.
- **THIRD, second half (normalizing `cumulative` in `sampler()`).** The first agent
  offers "either divide every element by `sum` or add the `debug_assert`"; the
  second agent argues against dividing, because `MixtureSampler::entries` is a
  `pub` field that would then still show sub-stochastic weights while the sampler
  silently normalized, so the sampler and the inspectable weights would disagree
  and the source bug would be masked. Both agree the `debug_assert` half is right.
  **The user rules on the normalize half.**

Both agents also flag that FIRST deliberately breaks bug-for-bug parity with
legacy, which has the identical gate — and that this crate's design is explicitly
parity-driven (the "Compatibility schedule" comment at `sampler.rs:44`, the
`burn_legacy_*` seed burns, `mixture/data.rs:98-100`'s "bypassing this door
silently changed that boundary"). The differential suite
`crates/ppvm-conformance-2/tests/tableau_mixture_diff.rs` exercises the measurement
path only at `sum_cutoff = 1e-14` (lines 116, 131, 143, 194), where the dropped
children have essentially zero probability and the extra normalize is a division by
1.0, so it is *expected* to still pass — but the divergence should be recorded
rather than discovered later, and legacy should be fixed in tandem or the parity
claim narrowed.

**Sign-off needed.** FIRST, SECOND and FIFTH change observable behaviour (entry
weights, reported Born triples, shot distributions, entry counts) and deliberately
diverge from legacy; FOURTH is contested; the `sampler()` normalize half is
contested.

### U5 measurement sign

Rows: **G-021**, **G-022**, **G-023**, **G-062**, **G-058**, **G-059** — all
`adjudicated-spec`. Verdict: **the code is correct.** Every ℤ/4 sign in this unit
is the mathematically right one; all six rows are verification/bridge/coverage
gaps, not defects. Second opinion: **corroborated** (independent derivation and an
independently written probe; both mutation results reproduced verbatim), with three
corrections.

**Question.** Is the shipped Aaronson–Gottesman ℤ/4 measurement sign in
`ppvm-tableau-2` (deterministic branch outcome, case-a projection, row g-rule,
extension-gate phases) mathematically correct, or does an unpinned sign let a wrong
outcome reach a user?

**Derivation.** A row is (x, z, phase) read as R = i^phase ⊗_j P(x_j, z_j) with the
Hermitian one-site basis P(x,z) = i^{xz} X^x Z^z, i.e. P(0,0) = I, P(1,0) = X,
P(0,1) = Z, P(1,1) = Y (since XZ = −iY). The discriminators for that convention are
empirical, not aesthetic: under P = X^x Z^z every product phase is 2bc, always
*even*, yet the shipped kernel returns 3 for X·Z, so that reading is excluded; and
the −Y variant is excluded by direct comparison against genuine 2×2 matrices (done,
not argued — see probe part 1). It is also forced by the tableau's own `s`
(x' = x, z' = x⊕z, phase += 2(x∧z)) sending X ↦ +Y and Y ↦ −X.

*(1) The g-rule (G-058).* From P(a,b)·P(c,d) = i^{ab+cd} X^a Z^b X^c Z^d and
Z^b X^c = (−1)^{bc} X^c Z^b:

  P(a,b)·P(c,d) = i^g P(a⊕c, b⊕d),  **g = ab + cd + 2bc − (a⊕c)(b⊕d) (mod 4)**,

which is Aaronson–Gottesman's g (their g(1,1,1,0) = z₂ − x₂ = −1 ≡ 3 matches
Y·X = −iZ). The shipped kernel (`data.rs:311-319`) is
`sign = (a&b&c&!d)|(a&!b&!c&d)|(!a&b&c&d)`,
`imag = (a&!b&d)|(a&!c&d)|(!a&b&c)|(b&c&!d)`,
`phase += (2·pc(sign) + pc(imag)) % 4`, then `+= rhs.phase`. Exhaustively over all
16 (a,b,c,d) the shipped per-site value equals both the closed form **and** the
phase extracted from the genuine 2×2 complex matrix product (0 mismatches). Two
anchors: X·Z ⇒ 3 (= −iY ✓), Z·X ⇒ 1 (= +iY ✓), which also fixes the orientation:
`self.mul_assign(rhs)` computes lhs·rhs, not rhs·lhs. So the tableau's private
g-rule *is* the Pauli-group product and `+= rhs.phase` is the correct cocycle
bookkeeping. G-058 is a **bridge** gap (no Lean/oracle tie), not an error. The
second agent additionally pinned the orientation empirically (1038 failures when
compared as src·dst, 0 as dst·src) and checked the cross-word popcount fold at
n = 70 with `[u64;2]` storage against the per-site closed form: 480 checks, 0
mismatches.

*(2) The deterministic (case-b) sign (G-021).* Case b is exactly "no stabilizer has
x[q] = 1", i.e. ω(Z_q, s_i) = 0 for all i, so Z_q lies in the stabilizer span; in
the symplectic frame (ω(d_i,s_j) = δ_ij) the coefficient of s_i is
ω(Z_q, d_i) = x_{d_i}[q] and the destabilizer coefficients vanish. Hence
Z_q = ±∏_{i∈T} s_i with T = {i : d_i anticommutes with Z_q} — precisely
`get_deterministic_outcome`'s selection rule. Let M = ∏_{i∈T} s_i (order-free:
stabilizers commute). M has Z_q's bits, so M = i^p Z_q; each s_i fixes |ψ⟩, so
M|ψ⟩ = |ψ⟩ and Z_q|ψ⟩ = i^{−p}|ψ⟩ = (−1)^{p/2}|ψ⟩ for p ∈ {0,2}. **Eigenvalue −1 ⟺
p = 2 ⟺ outcome bit 1**, and the code returns `result.phase >= 2`. Correct under the
universal convention (outcome b ⟺ eigenvalue (−1)^b). Absolute anchor: on |0…0⟩ the
frame is d_i = X_i, s_i = Z_i with phase 0, T = {q}, M = Z_q, p = 0 ⇒ outcome 0;
after X_q the row Z_q carries phase 2 ⇒ outcome 1. Both right.
Generalized route: `compute_decomposition(q, Z)` starts p_word = Z_q and multiplies
the same selected stabilizers with `add_phase(8 − 2·stab.phase)` (= 0 on real rows,
i.e. it forms s_i^{−1}), giving Z_q·M^{−1} = i^{−p}·I, so `phase_decomp` = p and
`z_sign = (phase_decomp == 2)` is the same bit. In the amplitude basis
|j⟩ = ∏_l d_l^{j_l}|ψ₀⟩, commuting Z_q past the selected destabilizers contributes
(−1)^{⟨L,j⟩}, so the per-index eigenvalue exponent is `phase_decomp + 2⟨L,j⟩` —
exactly `compute_phase_with_mask_static`'s first term — hence
`z_overlap_re = ⟨ψ|Z_q|ψ⟩` and `prob_1 = (1 − ⟨Z_q⟩)/2` (Born), and the retain
predicate keeps indices with eigenvalue (−1)^outcome, i.e.
⟨L,j⟩ ⊕ (p == 2) = outcome, which is literally `(parity ^ outcome) == z_sign`.
Correct.

*(3) The random (case-a) 1/2 (G-022).* If stabilizer S anticommutes with Z_q then,
using S† = S, S² = I (row realness, see (4)) and S|ψ⟩ = |ψ⟩:
⟨Z_q⟩ = ⟨ψ|S Z_q S|ψ⟩ = −⟨ψ|Z_q|ψ⟩ ⇒ ⟨Z_q⟩ = 0 ⇒ **p1 = 1/2 exactly**. So the bare
frame's `rng.random::<bool>()` *is* the Born distribution, and the coefficient
path's `0.5 − 0.5·z_overlap_re` collapses to the same number on a pure stabilizer
state (support 1 ⇒ the shifted key never matches ⇒ `z_overlap_re` is exactly 0.0).
The two samplers agree in distribution; they differ only in RNG-stream consumption
(one bool vs one f64), which the crate documents.

*(4) Row Hermiticity (G-062).* Initial rows are X_i/Z_i with phase 0; all 29 phase
writes in `clifford.rs` are of the form `phase ^= (predicate) << 1`, i.e. deltas of
0 or 2, so evenness is gate-invariant. `mul_assign` on frame rows is only ever
applied to *commuting* pairs (stabilizer × stabilizer in the projection; the pivot
g_q into destabilizer d_i only for i ≠ q_idx, where ω(g_q, d_i) = δ = 0), and for
commuting Hermitian P, Q we have (PQ)† = QP = PQ. Hence phase stays in {0,2},
`odd_phase_destabilizer_mask()` is **identically 0** on a live frame, and the mask
term in `compute_phase_with_mask_static` is vacuous — a three-line theorem, not a
coincidence. So **G-023's "assumed hypothesis" hψ is derivable**, and the
projection's omission of the mask term is harmless. (The intermediate p_word inside
`compute_decomposition` legitimately *can* be imaginary in case a, which is why the
ℤ/4 phase exists at all.)

*(5) Extension-gate signs (G-059).* √P = e^{−iπP/4} (the choice making √Z = S,
matching the crate's `s`) gives R_x(π/2): X↦X, Y↦Z, Z↦−Y and R_y(π/2): X↦−Z, Y↦Y,
Z↦+X, with the daggers as transposes — exactly the table at `clifford.rs:426-434`,
and each bit/phase rule implements it (`sqrt_x`: x' = x⊕z, +2 iff z∧¬x ⇒ Z = (0,1)
↦ (1,1) with +2 = −Y ✓, Y = (1,1) ↦ (0,1) unflipped = +Z ✓).
CY = |0⟩⟨0|⊗I + |1⟩⟨1|⊗Y gives I⊗X ↦ Z⊗X and X⊗X ↦ −Y⊗Z, both re-derived by hand
via CY = (I⊗S)·CX·(I⊗S†) and matching the code's table and its
`xc & (xt⊕zt) & ¬(zc⊕zt)` phase predicate. The second agent resolved the gate
identification convention-free with `expectation`: every named tableau gate
implements the **forward** textbook unitary on the state — s = diag(1,i) (not S†),
sqrt_x = e^{−iπX/4} = R_x(+π/2), sqrt_y = R_y(+π/2), t = diag(1, e^{iπ/4}), the
`*_dag` gates their adjoints — so the tableau's ROW action is P ↦ U P U†, the
**transpose** of the `CliffordExtensions` doc table at
`ppvm-traits-2/src/gates.rs:95-107` (documented as U†PU, "backward Heisenberg").
625,856 full-4ⁿ-Pauli-set expectation checks over 400 random 14-gate circuits
(n = 2..4), 0 failures, worst |Δ| 1.110e−15.

*(6) Case-a post-measurement state (G-023).* Bare frame: after each measurement
every stabilizer row (with its i^phase) must fix the dense projected state — a
phase-exact test, not up-to-global-phase. Generalized: the first agent
reconstructed the dense state as Σ_j c_j (∏_{l: j_l=1} D_l)|ψ₀⟩ with |ψ₀⟩ obtained
by applying ∏_i (I+S_i)/2 to a fixed generic vector (well defined and order-free
because destabilizers pairwise commute; the global phase of ψ₀ cancels) and got
fidelity 1.000000000000 over 8400 comparisons after every step of 600 Clifford+T
circuits, including after case-a merges and the subsequent
`update_tableau_according_to_outcome` re-basing. The second agent avoided the basis
convention entirely by comparing the **full 4ⁿ Pauli expectation set** — which pins
the density matrix — and got 939,648 comparisons, 0 failures, worst |Δ| 7.494e−16.
So the merge phase `alpha + 2·⟨idx, destab⟩`, the omission of the odd-phase-mask
term there (justified by (4)), and the re-indexing induced by the frame update are
all correct.

**Correct answer.** Deterministic branch: outcome bit b = (accumulated ℤ/4 phase
== 2), i.e. b = 1 ⟺ Z_q eigenvalue −1 — what `data.rs:497` (`result.phase >= 2`)
and `measure.rs:336` (`z_sign = phase_decomp == 2`) already do. Row product phase =
g = ab + cd + 2bc − (a⊕c)(b⊕d) mod 4 = the shipped `(2·pc(sign) + pc(imag)) % 4`
kernel, exactly. Case-a probability = 1/2 exactly on any stabilizer state, so the
bare frame's fair coin and `0.5 − 0.5·z_overlap_re` agree. Row phases are provably
in {0,2}, so the odd-phase-destabilizer mask is identically 0. All six extension
gates match √P = e^{−iπP/4} / CY = ctrl-Y conjugation.

**Which implementation is right.** Both measurement backends, and they agree.
`Tableau::get_deterministic_outcome` computes phase(M) with M = ∏_{i∈T} s_i;
`GeneralizedTableau`'s `compute_decomposition` computes phase(Z_q·M^{−1}); since
M = i^p Z_q with Z_q² = I those are the same ℤ/4 quantity. The two case-a samplers
are the same distribution because prob_1 is provably exactly 0.5 on a pure
stabilizer state.

**Live defect: NO — no wrong number reaches a user.** Across the two agents,
25,077 + 15,340 independently adjudicated deterministic measurements (11,754 + 308
and 7,162 respectively with true eigenvalue −1 — the cases a +2 phase mutation
would corrupt) produced **zero** outcome mismatches against a dense state-vector
simulator, and 109,200 + 939,648 phase-exact post-measurement state checks passed
with residual 0 / fidelity 1 / worst |Δ| 7.5e−16. The six rows are
verification-strength, bridge and coverage gaps: the asserted invariants are true,
they are simply not pinned by Lean or by a mutation-sensitive oracle.

**Proposed fix.** **No change to any `crates/*/src` file.** Observable behaviour
must NOT change: the shipped signs are correct and any "fix" here would introduce a
bug. Three test/doc actions, in priority order:

1. `crates/ppvm-conformance-2/tests/tableau_lean.rs` — add the mutation-sensitive
   oracles the rows ask for, as *anchors* rather than re-derivations: (a) a full
   ℤ/4 g-rule test comparing a transcribed `mulWord`/`phaseExpN` against a real
   `Row` product, asserting the orientation `self·rhs` as well as the magnitude;
   (b) **the highest-value addition** — an absolute-sign anchor for the bare frame:
   `measure(q)` on `|0…0⟩` is `Some(false)` and on `X_q|0…0⟩` is `Some(true)`
   (both verified; both mutations below break it); (c) `assert!(phase % 2 == 0)` on
   all 2n rows inside `assert_symplectic_frame` (G-062); (d) extend the per-gate row
   table at `tableau_lean.rs:389`/`:439` with s_dag/sqrt_x(_dag)/sqrt_y(_dag)/cy —
   **transposing** the reference predicates relative to the `CliffordExtensions` doc
   table, or the new test will fail on correct code; (e) a case-a
   `prob_1 == 0.5` assertion for Clifford-only generalized tableaus (exact equality
   is legitimate: support 1 ⇒ `z_overlap_re` is exactly 0.0; a T-containing
   generalization would need a tolerance).
2. This ledger — correct two overstated premises, done in the row texts below.
3. Lean (`BranchPhase.lean`, `Frame.lean`) — the rows' proposed closures are the
   right ones and can now be written against a *verified* target: lift
   `phase_decomp` from the ℤ/2 shadow to full ℤ/4, add `IsRealFrame`, derive
   `oddPhaseMask = 0`.

**Evidence.** Probe 1, clean shipped code
(`cargo test -p ppvm-conformance-2 --test zz_adj_meas -- --nocapture`), verbatim
tail:

```
(x1z1)=(0,0) (x2z2)=(0,0) shipped_g=0 closed_form=0 matrix_ok=true
(x1z1)=(1,0) (x2z2)=(0,0) shipped_g=0 closed_form=0 matrix_ok=true
(x1z1)=(0,1) (x2z2)=(0,0) shipped_g=0 closed_form=0 matrix_ok=true
(x1z1)=(1,1) (x2z2)=(0,0) shipped_g=0 closed_form=0 matrix_ok=true
(x1z1)=(0,0) (x2z2)=(1,0) shipped_g=0 closed_form=0 matrix_ok=true
(x1z1)=(1,0) (x2z2)=(1,0) shipped_g=0 closed_form=0 matrix_ok=true
(x1z1)=(0,1) (x2z2)=(1,0) shipped_g=1 closed_form=1 matrix_ok=true
(x1z1)=(1,1) (x2z2)=(1,0) shipped_g=3 closed_form=3 matrix_ok=true
(x1z1)=(0,0) (x2z2)=(0,1) shipped_g=0 closed_form=0 matrix_ok=true
(x1z1)=(1,0) (x2z2)=(0,1) shipped_g=3 closed_form=3 matrix_ok=true
(x1z1)=(0,1) (x2z2)=(0,1) shipped_g=0 closed_form=0 matrix_ok=true
(x1z1)=(1,1) (x2z2)=(0,1) shipped_g=1 closed_form=1 matrix_ok=true
(x1z1)=(0,0) (x2z2)=(1,1) shipped_g=0 closed_form=0 matrix_ok=true
(x1z1)=(1,0) (x2z2)=(1,1) shipped_g=1 closed_form=1 matrix_ok=true
(x1z1)=(0,1) (x2z2)=(1,1) shipped_g=3 closed_form=3 matrix_ok=true
(x1z1)=(1,1) (x2z2)=(1,1) shipped_g=0 closed_form=0 matrix_ok=true
g-rule mismatches (of 16): 0
test g_rule_matches_real_matrices ... ok
--- case-a Born frequency vs dense p1 (T-containing states) ---
--- bare Tableau vs state vector ---
deterministic measurements checked: 24193  mismatches: 0  (of which truth=1: 11754)
random-branch measurements checked: 4607  prob != 1/2: 0  worst |p1-0.5| = 2.220e-16
case-a coin: 2291 ones of 4607 draws (freq 0.4973)
stabilizer-fixes-state checks: 100800  failures: 0  worst residual = 0.000e0
rows with odd (imaginary) phase: 0
test bare_tableau_vs_statevector ... ok
born-frequency cases: 73  worst |freq-p1|/sigma = 2.69
--- GeneralizedTableau vs state vector ---
state comparisons: 8400  failures: 0  worst fidelity: 1.000000000000
case-b (deterministic) outcomes checked: 884  mismatches: 0  (truth=1: 308)
case-a draws: 916  worst |p1-0.5| over Clifford-only runs: 3.886e-16
test generalized_tableau_vs_statevector ... ok

test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 3.10s
```

Circuits: n = 2..5 over x/y/z/h/s/s_dag/sqrt_x/sqrt_x_dag/sqrt_y/sqrt_y_dag/cnot/
cz/cy plus t for the generalized half, measurements interleaved, each followed by a
repeat measurement and an X-then-measure to force −1-eigenvalue deterministic
branches.

Mutation A, `p_word.add_phase(2)` at the end of `compute_decomposition` — the
existing Lean suite **does** catch it:

```
running 17 tests
test measurement_dichotomy_holds ... FAILED
...
thread 'measurement_dichotomy_holds' panicked at crates/ppvm-conformance-2/tests/tableau_lean.rs:700:13:
assertion `left == right` failed: seed 0: measure(0) is not idempotent
  left: false
 right: true
test result: FAILED. 16 passed; 1 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.65s
```

Mutation B alone, `result.add_phase(2)` in `get_deterministic_outcome` — the Lean
suite is blind, the crate's behaviour tests catch it, the probe catches every case:

```
test result: ok. 17 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.66s     [tableau_lean.rs]

DET MISMATCH q=1 p1=1 reported=false truth=true
DET MISMATCH q=1 p1=0 reported=true truth=false
DET MISMATCH q=0 p1=0 reported=true truth=false
DET MISMATCH q=0 p1=0 reported=true truth=false
DET MISMATCH q=0 p1=1 reported=false truth=true
deterministic measurements checked: 24193  mismatches: 24193  (of which truth=1: 11754)
random-branch measurements checked: 4607  prob != 1/2: 0  worst |p1-0.5| = 2.220e-16
stabilizer-fixes-state checks: 100800  failures: 0  worst residual = 0.000e0

---- frame_reset_returns_to_zero stdout ----
thread 'frame_reset_returns_to_zero' panicked at crates/ppvm-tableau-2/tests/behaviour.rs:649:5:
assertion `left == right` failed
  left: Some(true)
 right: Some(false)
---- frame_case_b_measurement_leaves_the_rng_untouched stdout ----
thread 'frame_case_b_measurement_leaves_the_rng_untouched' panicked at crates/ppvm-tableau-2/tests/behaviour.rs:627:9:
assertion `left == right` failed
  left: Some(true)
 right: Some(false)
test result: FAILED. 35 passed; 2 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.05s
```

Both mutations were reverted with `git checkout -- crates/ppvm-tableau-2/src/data.rs`
(`git diff --stat` empty for that crate) and the scratch file deleted. The second
agent reproduced both results verbatim.

Second-agent probe totals (own file, deleted): row-product checks 4400 with 0
failures as dst·src and 1038 as src·dst, 16/16 per-site combos; cross-word fold at
n = 70, 480 checks, 0 mismatches; bare frame 13,829 deterministic measurements
(truth = 1: 6,506) 0 mismatches, 4,171 random-branch with worst |p1 − 0.5| =
4.441e−16 and 2,084 ones of 4,171 draws, 63,000 phase-exact stabilizer-fixes-state
checks with worst residual 1.911e−15, 0 odd-phase rows; generalized 1,511
deterministic outcomes (truth = 1: 656) 0 mismatches, 2,800 ⟨Z⟩ comparisons worst
8.882e−16, 939,648 post-op full-Pauli state checks 0 failures worst 7.494e−16,
1,289 case-a branches (654 Clifford-only, worst |p1 − 0.5| = 3.331e−16); Born
frequencies over 4,000 shots/case, 26 cases, worst 2.69σ, no outlier > 5σ.

**Second opinion — corroborated, with three corrections.**

1. **The proposed g-rule oracle needs no new API surface.** The first agent wrote
   that a full ℤ/4 g-rule test "needs a `#[cfg(test)]`/`pub` row-multiply hook in
   `ppvm-tableau-2`". It does not: `StabilizerFrame::row_multiply` is a public trait
   method and `Tableau::rows()` is `pub`, and the second agent drove an exhaustive
   external matrix-grounded g-rule oracle from
   `crates/ppvm-conformance-2/tests/` with zero source changes. This matters because
   adding a hook would be a `src` change to a crate the fix claims not to touch.
   G-058's closure has been amended accordingly.
2. **Convention nit.** "P(1,1) = Y is the only convention under which phase ∈ {0,2}
   means Hermitian generator" is not the discriminator (P(1,1) = −Y also gives
   Hermitian generators at even phase). The real discriminators are the kernel
   returning an *odd* value (X·Z ⇒ 3), which kills X^x Z^z, and direct matrix
   comparison, which kills the −Y variant.
3. **Sharper statement of the gate convention**, folded into (5) above: the tableau
   implements the forward unitary (rows conjugate as P ↦ U P U†), the transpose of
   the `CliffordExtensions` doc table at `ppvm-traits-2/src/gates.rs:95-107`. The
   doc table is a real trap for anyone writing the G-059 oracle from it (R-14 covers
   this only in prose). Also: the first agent's claim that the extension gates were
   "exercised with phase-exact stabilizer comparisons, so a shipped-adjoint error
   would have appeared" holds only because the circuits left X/Y/Z eigenstates — a
   naive one-gate probe from |+⟩ cannot distinguish `sqrt_x` from `sqrt_x_dag` at
   all, so the discrimination comes from the random-circuit body (as it did for both
   agents).

**Contested inside the fix.** None. Both agents agree: no `src` change, and the
only judgement call is spending effort on Lean/oracle strength rather than on code.

**Sign-off needed.** None — nothing observable changes.

### U6 unenforced preconditions

Rows: **G-061** (`adjudicated-defect`), **G-060** (`adjudicated-defect`),
**G-054** (`adjudicated-defect`), **G-057** (`adjudicated-spec`),
**G-008** (`adjudicated-spec`). Verdict: **code is wrong** for three of the five.
Second opinion: **corroborated** (independent matrix-level derivation, own probes,
identical numbers), with three corrections — and two of the proposed fix items have
real problems.

**Question.** For each of the five Lean hypotheses with no Rust counterpart (Nodup
in the batch sweeps, pairwise-disjoint support in the fused CZ block, exact halving,
unbounded ℤ[i], additive associativity of a thresholded `Term`), is the precondition
reachable through a public API, and does violating it yield a wrong number, a panic,
or nothing?

**Ranking by reachability** (both agents, same order): **G-061** (a legal `.stim`
file today) > **G-060** (a `pub` fused-block call with offset < count or offset == 0)
> **G-054** (exact ring past |z| ≈ 3e9, release only) > **G-057** (needs
`set_min_eps` above default) > **G-008** (unreachable: the only call site is
`half` of 1).

**G-061, Nodup.** Applying gate G at qubit q twice is conjugation by G²: X² = Y² =
Z² = H² = I, so `g_many([q,q])` must be the **identity** map on every row; S² = Z,
(√X)² = X, (√Y)² = Y, so those must equal Z/X/Y-conjugation. `build_masks`
(`clifford.rs:594`) **ORs** each index into a per-word mask, so `[q,q]` contributes
one bit and the fused sweep applies G **once**. Probed: for all ten mask-based gates
(x, y, z, h, s, s_dag, sqrt_x, sqrt_x_dag, sqrt_y, sqrt_y_dag) `g_many(&[0,0])` is
bit-for-bit equal to a *single* `g(0)` and unequal to `g(0); g(0)`; for x/y/z/h the
per-index loop is the identity (mathematically correct) and the batch is not.
`cnot_many`/`cz_many`/`cy_many` are **not** affected — they iterate pairs and
re-read the bits each iteration, so they are correct even on repeated pairs.
Reachability is not hypothetical: `ppvm-stim/src/executor/gates.rs:36` lowers a Stim
single-qubit gate's whole target list with `tab.x_many(&qubits(targets))`,
`qubits()` (`executor/helpers.rs:17`) is a bare map with no dedup, and
`validate.rs` checks distinctness only for MPP products
(`check_mpp_distinct_qubits`), not for gate target lists. Stim's semantics for
`X 0 0` is "apply to each target in order" = identity. The second agent ran the real
public entry point `ppvm_stim::run_string` from a throwaway crate outside the repo:
`"X 0 0\nM 0"` → `[Some(true)]` while `"X 0\nX 0\nM 0"` → `[Some(false)]`, i.e. **a
flipped sampled bit from a legal `.stim` file** — and got identical results with
default features, because legacy `ppvm-tableau/src/gates/clifford.rs:365` ORs too.
So this is a pre-existing defect faithfully reproduced by the `-2` crate and **live
in the shipped default build as well**, not a `-2` regression.
`ppvm-vihaco`'s `CircuitMessage::QubitBatch` (`vihaco-circuit-isa/src/lib.rs:152`)
is likewise an unvalidated decoded operand dispatched to `x_many`
(`component/tableau.rs:86`).

**G-060, disjoint support.** CZ = diag(1,1,1,−1) and CZ·X_c·CZ = X_c Z_t, so
CZ(X⊗X)CZ = (X_cZ_t)(X_tZ_c) = (X_cZ_c)⊗(Z_tX_t) = −(XZ)⊗(XZ) = −(−iY)⊗(−iY) =
+Y⊗Y (verified entrywise). Then for U = CZ₁₂·CZ₀₁ on X₀X₁X₂: step 1 gives Y₀Y₁X₂;
CZ(X⊗Y)CZ = −Y⊗X so by symmetry CZ(Y⊗X)CZ = −X⊗Y, giving **−Y₀X₁Y₂**, i.e. in row
coordinates x = 111, z = 101, phase 2. The scalar per-pair loop reproduces exactly
that; the fused `cz_block_pairs(0,1,2)` produces x = 111, z = 111, phase 0 =
**+Y₀Y₁Y₂** — a different Pauli **and** a different sign. Two independent causes:
(a) the parity `xc & xt & (zc ^ zt)` reads pre-update z for all pairs at once, which
is Lean's `czSeq_phase_needs_disjoint`; (b) `z_delta = ((x>>offset)&mask_c) |
((x<<offset)&mask_t)` uses **OR**, so a z bit written by two overlapping pairs is
set once instead of XOR-cancelled. **(b) refutes G-060's own claim that "the bits
still come out right".** `cz_block(0,1,2)` on `GeneralizedTableau` reproduces it
(it forwards run = 2, offset = 1). The ledger's cited `cz_block_pairs(0,2,5)` also
diverges: on 200 randomly prepared 10-qubit tableaux the second agent found it
disagreeing with the per-pair loop on **158** (a single hand-picked state happened
to agree, which is why a randomized sweep matters). Degenerate offset == 0 (pairs
(q,q)) sets z ^= x over the whole block while the scalar `cz(q,q)` is a probed
no-op. The only guard is a `debug_assert_eq!` that all bits share one word; nothing
checks `offset >= count`. Both entry points are `pub`, and **`cz_block(0,1,n)` — CZ
on every adjacent pair, the most natural brickwork call a user would write — is
exactly the broken case**.

**G-008, exact halving.** `half` is `x/2.0`, and doubling a finite float never
rounds (same significand, exponent + 1; it can only overflow, and f64::MAX doubles
back to MAX under `h+h`). So `x.half() + x.half() == x` ⟺ **x/2 is representable**.
Writing x = m·2^−1074: for biased exponent ≥ 2 (|x| ≥ 2^−1021) halving just
decrements the exponent and is always exact; for the lowest normal binade
[2^−1022, 2^−1021) and for subnormals it needs (2^52+f)/2 resp. f/2 to be an
integer, i.e. the **raw bit pattern must be even**. Verified exhaustively over
4×2^17 bit patterns: every failure had an odd bit pattern *and* bits < 2<<52, and
there was no failure at or above 2^−1021. So the docstring's universally quantified
law is false — **and G-008's own proposed weakening ("exact wherever x is normal and
2·x is representable") is also false**: `f64::from_bits(0x10000000000001)` is normal
with 2x finite and fails. `Complex<f64>::half` divides componentwise and inherits
this. The only call sites in the `-2` tree are `C::one().half()` at
`ppvm-pauli-sum-2/src/proj.rs:74` and `:92` — x = 1 exactly, where the law holds —
so nothing depends on the false region.

**G-054, ℤ[i].** The positive claim is true and trivial: z + z = 1 + i needs
2·re = 1, impossible in ℤ (probed: no witness on a grid; the proof is parity of the
real part). The hazard is the representation: `norm_sq` (`exact.rs:72`) is
`re*re + im*im` in **i64** and `magnitude` is its sqrt, so both wrap once |z|² ≥
2^63, i.e. |z| ≳ 2^31.5 ≈ 3.04e9 — far below where the *value* stops being
representable. (1+i)^64 = (2i)^32 = 2^32 + 0i is computed exactly by 64 `Mul`s
(every intermediate product is by ±1), yet true norm_sq = 2^64 wraps to 0 in
release ⇒ `magnitude() = 0.0`; for re = 3037000500 the wrap is
−9223372036709301616 ⇒ `magnitude() = NaN`, and `NaN >= threshold` is false. Either
way `CoefficientThreshold::truncate`'s keep-rule `coeff.magnitude() >= threshold`
(`ppvm-pauli-sum-2/src/policy.rs:216`) **silently deletes a coefficient of true
modulus ≈ 4.3e9**. In debug both panic ("attempt to multiply with overflow" at
`exact.rs:73`), which is the only detection; release has none. The wrapping contract
*is* documented (`exact.rs:41-43`), but the laws on the same type are not
range-qualified: `norm_sq`'s doc says "exact in ℤ and strictly multiplicative" and
`magnitude`'s says it "satisfies every clause of the documented law … multiplicative
(|zw| = |z||w|)", both false past 2^31.5.

**G-057, additive associativity.** `Sum + Const` routes through `Sum::add_const`
(`term.rs:323`) which **drops** |c| < min_eps; `Sum + Sum` adds `s2.c0`
unconditionally (`add.rs:142`) and `Const + Const` adds unconditionally
(`add.rs:197`). With a = 1 + sin(x0) in Sum form and min_eps = 1e−3:
(a + 6e−4) + 6e−4 has c0 = 1.0 while a + (6e−4 + 6e−4) has c0 = 1.0012, and `Term`'s
`PartialEq` sees it. It is also non-commutative in effect, because truncation
parameters are inherited from the lhs only: a + b → c0 1.0, b + a → c0 1.0006.
**Exact associativity is unattainable for any nonzero thresholding rule** — that is
Lean's own `eps_drop_at_insert_ne_drop_at_end` — so the non-associativity is an
accepted consequence of thresholding, not a numerical bug; the defect is bounded by
(#dropped constants)·min_eps in c0. What *is* defective is the bookkeeping:
(i) `GradedMap.accumulate_assoc` is `add_assoc` on an untruncated `CMap`;
(ii) `add.rs:66-70` invokes that law to argue "the coefficient ring must be an
additive monoid … which forbids `x + c == x` for c ≠ 0" while the same file ships
`x + c == x` for 0 < |c| < min_eps in the Sum+Const arm; (iii) the arms are mutually
inconsistent, so `+` is not a function of the value a `Term` denotes but of which of
Const/One/Sum represents it (the additive analogue of `mulImpl_not_wellDefined`).
`min_eps` defaults to `f64::EPSILON` and no production path calls `set_min_eps`
(only the conformance harness and benches), so at the default the deviation is at
rounding-error scale.

**And the "holds only at eps = 0" escape is closed.** The second agent probed
min_eps = 0.0 with a = 1 + sin(x0), b = c = 1e−16: (a+b)+c has c0 = 1.0 while
a+(b+c) has c0 = 1.0000000000000002 — **f64 addition is itself non-associative, so
`Term`'s `+` is non-associative at *every* min_eps.** `accumulate_assoc` is a law of
an exact coefficient ring and is not a law of an f64-backed `Term` at any eps. This
directly invalidates part of G-057's proposed closure (see the closure edit there):
shipping the scope note "`accumulate_assoc` holds only at eps = 0" would install a
new false law where the old one was.

**Live defect: YES for G-061, G-060 and G-054.**

- **G-061:** a Stim circuit containing a repeated single-qubit target (`X 0 0`,
  legal Stim, semantically the identity) produces a flipped stabilizer sign and a
  flipped sampled measurement bit — measured through `run_string`: `Some(true)`
  where the truth is `false`.
- **G-060:** a caller of the `pub` fused CZ block with overlapping pairs gets the
  wrong Pauli in **both** bit planes and the wrong ℤ/4 sign (+Y₀Y₁Y₂ instead of
  −Y₀X₁Y₂), which then propagates into every measurement sign.
- **G-054:** in release, an exact ℤ[i] coefficient of modulus ≥ 2^31.5 reports
  magnitude 0.0 or NaN and is silently truncated away by the policy keep-rule; in
  debug it panics.
- **G-057:** with a user-set `min_eps`, (a+b)+c ≠ a+(b+c) in the constant term
  (1.0 vs 1.0012) — a truncation-scale discrepancy, not a gross error.
- **G-008:** none — mislabelling only.

**Proposed fix.**

1. **G-061**, `crates/ppvm-tableau-2/src/clifford.rs:594` `build_masks` — detect a
   duplicate index and **fall back to the per-index `Clifford` loop** (always
   correct by the trait contract) for every gate; apply the same to
   `y_many_skipping` (`:602`). Cost is one test-and-branch per index against a
   2n-row sweep. Independently,
   `crates/ppvm-stim/src/executor/gates.rs`/`helpers.rs` must expand or
   apply-per-target, since Stim's semantics is one application per target.
   **LOUD BEHAVIOUR CHANGE** (`X 0 0\nM 0` goes from wrong to right; measurement
   records for such circuits move), and it makes old and new disagree by design
   unless legacy is fixed in tandem. See **Contested** for the two variants that
   must NOT be used.
2. **G-060**, `crates/ppvm-tableau-2/src/data.rs:1483` — add
   `debug_assert!(offset >= count, "cz_block pairs must have pairwise-disjoint supports (Batch.lean::czSeq_phase_needs_disjoint)")`
   next to the existing same-word assert, and in `cz_block` (`data.rs:1621`) either
   assert `hi - lo >= count` or fall back to the per-pair `Clifford::cz` loop when
   `hi - lo < count`; document the precondition in both signatures. The predicate is
   exactly right: pairs (base+i, base+offset+i) collide iff i − j = offset (needs
   offset < count) **or** offset = 0, so `offset >= count` subsumes the separate
   offset == 0 case. No existing caller violates it — `tableau_lean.rs:557`
   (0,20,10), `tableau_behaviour_diff` (0,17,17) (offset == count is the boundary
   and is legitimately disjoint), the `tableau_diff` cases (0,32,17), (17,34,17),
   (60,70,8), (0,64,20), and `cz_block`'s per-segment calls (run ≤ count ≤ offset).
   The `debug_assert` half changes no release behaviour; the release-mode per-pair
   **fallback** is an observable change for overlapping calls (wrong → right) and
   must be called out.
3. **G-008**, docs only, `crates/ppvm-traits-2/src/coefficient.rs:141` — replace
   "Impls must be exact: `x.half() + x.half() == x`" with "Impls must be exact
   wherever `x/2` is representable; for `f64`/`Complex<f64>` that is every
   `|x| >= 2·f64::MIN_POSITIVE` (2^−1021), plus 0 and ±inf, and below that it holds
   iff the raw bit pattern is even — it fails for `f64::from_bits(1)`,
   `from_bits(3)` and the *normal* `from_bits(0x10000000000001)`." Optionally note
   at `proj.rs:74` that the law is only ever used at x = 1.
4. **G-054**, `crates/ppvm-sym-2/src/exact.rs:72`/`:151` — range-qualify both law
   claims (|z| < 2^31.5 for `norm_sq`/`magnitude`, and note that `Mul` wraps at the
   same scale whenever both operands are large). Minimal safe code change: leave
   `norm_sq` alone (documented range, debug panic) and compute `magnitude` as
   `(re as f64).hypot(im as f64)`, which cannot overflow, is correct to < 1 ulp for
   all i64, and removes the release-mode 0.0/NaN truncation path. Tolerance check
   done: the only assertions on `GaussianInt::magnitude` are `sym_lean.rs:259-282`
   (multiplicativity/subadditivity with 1e-9 slack, plus exact
   `magnitude() == 0.0 iff is_zero`), all of which `hypot` satisfies.
5. **G-057**, `crates/ppvm-sym-2/src/add.rs:66-70` — strike the
   `accumulate_comm`/`accumulate_assoc` citation as justification for the shipped
   `+`, and state instead that `accumulate_assoc` is a law of an exact ring, that
   `Term`'s `+` is non-associative at **every** min_eps (min_eps truncation above
   the default, f64 rounding at and below it), and that the defect is bounded by
   (#dropped constants)·min_eps plus rounding; note that Sum+Const truncates while
   Sum+Sum's c0 and Const+Const do not. Pin the 6e−4/6e−4 witness at min_eps = 1e−3
   in `tests/sym_lean.rs`.

**Evidence.** Probe `crates/ppvm-conformance-2/tests/zz_adj_pre.rs` (deleted),
`cargo test -p ppvm-conformance-2 --test zz_adj_pre -- --nocapture --test-threads=1`,
verbatim:

```
test g008_half_exactness ... G-008 smallest subnormal           x=bits(0x1) half+half==x: false (half=0e0, sum=0e0, x=5e-324) normal=false 2x_finite=true
G-008 3 ulp subnormal              x=bits(0x3) half+half==x: false (half=1e-323, sum=2e-323, x=1.5e-323) normal=false 2x_finite=true
G-008 2 ulp subnormal              x=bits(0x2) half+half==x: true  (half=5e-324, sum=1e-323, x=1e-323) normal=false 2x_finite=true
G-008 smallest normal              x=bits(0x10000000000000) half+half==x: true  (half=1.1125369292536007e-308, sum=2.2250738585072014e-308, x=2.2250738585072014e-308) normal=true 2x_finite=true
G-008 smallest normal + 1ulp       x=bits(0x10000000000001) half+half==x: false (half=1.1125369292536007e-308, sum=2.2250738585072014e-308, x=2.225073858507202e-308) normal=true 2x_finite=true
G-008 2*MIN_POSITIVE + 1ulp        x=bits(0x20000000000001) half+half==x: true  (half=2.225073858507202e-308, sum=4.450147717014404e-308, x=4.450147717014404e-308) normal=true 2x_finite=true
G-008 1.0                          x=bits(0x3ff0000000000000) half+half==x: true  (half=5e-1, sum=1e0, x=1e0) normal=true 2x_finite=true
G-008 f64::MAX                     x=bits(0x7fefffffffffffff) half+half==x: true  (half=8.988465674311579e307, sum=1.7976931348623157e308, x=1.7976931348623157e308) normal=true 2x_finite=false
G-008 Complex<f64>(2^-1074, 1.0) half+half==x: false
G-008 sampled 65599 bit patterns: 32780 failures; for x < 2^-1021 (bits < 2<<52) the law fails iff the bit pattern is ODD; for x >= 2^-1021 it never fails
ok
test g054_gaussian_int_overflow ... G-054 (1+i)^64 = 4294967296  (re=4294967296, im=0)
G-054 true norm_sq = 2^64 = 18446744073709551616; i64-wrapping norm_sq (release) = 0; magnitude would be 0
thread 'g054_gaussian_int_overflow' panicked at crates/ppvm-sym-2/src/exact.rs:73:9: attempt to multiply with overflow
G-054 debug-build z.norm_sq() panicked: true
thread 'g054_gaussian_int_overflow' panicked at crates/ppvm-sym-2/src/exact.rs:73:9: attempt to multiply with overflow
G-054 debug-build z.magnitude() panicked: true
G-054 z=3037000500: wrapping norm_sq = -9223372036709301616 (negative: true), sqrt = NaN
G-054 exists z with z+z == 1+i (|re|,|im| <= 4): false
ok
test g057_term_add_non_associative ... G-057 a = Sum(Sum { c0: 1.0, maps: Some(SumMaps { terms: {Prod { factors: [Factor { var: 0, sin: 1, cos: 0 }], sin_pow: 1, cos_pow: 0, phase: 0 }: 1.0}, aux: {} }) })
G-057 min_eps = 1e-3
      (a+b)+c c0 = 1.0
      a+(b+c) c0 = 1.0012
      equal: false
G-057 3x4e-4: left c0 = 1.0, right c0 = 1.0012, equal: false
G-057 a+b c0 = 1.0 vs b+a c0 = 1.0006 (min_eps inherited from lhs only), equal: false
ok
test g060_cz_block_overlapping_support ... G-060 prepared rows (x,z,phase):
    row 0: x=111 z=000 ph=0
    row 1: x=010 z=000 ph=0
    row 2: x=100 z=000 ph=0
    row 3: x=000 z=001 ph=0
    row 4: x=000 z=011 ph=0
    row 5: x=000 z=101 ph=0
G-060 cz_block_pairs(0,1,2) == cz(0,1);cz(1,2): false
    row 0: fused x=111 z=111 ph=0   |   seq x=111 z=101 ph=2
G-060 cz_block_pairs(0,0,3) == cz(q,q) loop: false
    row 0: fused x=111 z=111 ph=0   |   seq x=111 z=000 ph=0
    row 1: fused x=010 z=010 ph=0   |   seq x=010 z=000 ph=0
    row 2: fused x=100 z=100 ph=0   |   seq x=100 z=000 ph=0
G-060 cz_block(0,1,2) == cz(0,1);cz(1,2): false
    row 0: fused x=111 z=111 ph=0   |   seq x=111 z=101 ph=2
ok
test g061_duplicate_index_batch_vs_loop ... G-061 x            batch==loop: false  loop==identity: true   batch==identity: false
    row 3: batch (x=0b0,z=0b1,ph=2) vs loop (x=0b0,z=0b1,ph=0)
G-061 y            batch==loop: false  loop==identity: true   batch==identity: false
    row 0: batch (x=0b1,z=0b0,ph=2) vs loop (x=0b1,z=0b0,ph=0)
    row 3: batch (x=0b0,z=0b1,ph=2) vs loop (x=0b0,z=0b1,ph=0)
G-061 z            batch==loop: false  loop==identity: true   batch==identity: false
    row 0: batch (x=0b1,z=0b0,ph=2) vs loop (x=0b1,z=0b0,ph=0)
G-061 h            batch==loop: false  loop==identity: true   batch==identity: false
    row 0: batch (x=0b0,z=0b1,ph=0) vs loop (x=0b1,z=0b0,ph=0)
    row 3: batch (x=0b1,z=0b0,ph=0) vs loop (x=0b0,z=0b1,ph=0)
G-061 s            batch==loop: false  loop==identity: false  batch==identity: false
    row 0: batch (x=0b1,z=0b1,ph=0) vs loop (x=0b1,z=0b0,ph=2)
G-061 s_dag        batch==loop: false  loop==identity: false  batch==identity: false
    row 0: batch (x=0b1,z=0b1,ph=2) vs loop (x=0b1,z=0b0,ph=2)
G-061 sqrt_x       batch==loop: false  loop==identity: false  batch==identity: false
    row 3: batch (x=0b1,z=0b1,ph=2) vs loop (x=0b0,z=0b1,ph=2)
G-061 sqrt_y       batch==loop: false  loop==identity: false  batch==identity: false
    row 0: batch (x=0b0,z=0b1,ph=2) vs loop (x=0b1,z=0b0,ph=2)
    row 3: batch (x=0b1,z=0b0,ph=0) vs loop (x=0b0,z=0b1,ph=2)
G-061 cnot(0,1)x2  batch==loop: true   loop==identity: true   batch==identity: true
G-061 cz(0,1)x2    batch==loop: true   loop==identity: true   batch==identity: true
G-061 cz(0,0) is a no-op: true
ok
test g061_measurement_outcome_flips ... G-061 measure(0) after batched X 0 0 = Some(true); after two single X = Some(false) (truth: |0>, so 0)
ok

test result: ok. 6 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
```

(One earlier run failed an assertion of the *probe*, not of the code: the first
agent had assumed raw-bit parity tracked significand parity above the lowest
binade. `bits(0x20000000000001)` satisfies the law, which corrected the
characterization to the 2^−1021 boundary stated above.)

**Second opinion — corroborated, with three corrections.**

1. **G-057's "holds only at eps = 0" is false** — probed (1.0 vs
   1.0000000000000002 at min_eps = 0.0 with b = c = 1e−16). Recorded above and
   applied to G-057's closure.
2. **G-061's divergence is not confined to `--features traits-2`** — legacy ORs
   identically and `run_string` returns `[Some(true)]` on both backends, so this is
   a pre-existing defect live in the default build, which changes how the fix must
   be framed (fixing only the `-2` side makes old and new disagree by design).
3. **G-054's "faithful up to |re|,|im| < 2^63" is too generous** — true of the
   representation and of Add/Sub/Neg, but `Mul` computes `re*rhs.re` in i64, so two
   operands of magnitude ≳ 2^31.5 wrap as well; the (1+i)^64 example only escapes
   because every intermediate product is by ±1.

Additions, not disagreements: (a) `cz_block(0,1,n)` — a nearest-neighbour CZ chain,
the most natural use of the fused block — is precisely the broken offset < count
case; (b) `cz_block_pairs(0,2,5)`, the ledger's own witness, diverges on 158 of 200
random 10-qubit states.

**Contested inside the fix — two items where the first agent's proposal is unsafe
and the second agent's objection is decisive.**

- **G-061 option (a), "signal duplicates from `build_masks`".** `build_masks` cannot
  do that: its `Option` return means "nothing to do" and every caller does
  `let Some(...) = ... else { return }`, so signalling duplicates via `None` would
  silently turn the whole gate into a **no-op** — strictly worse than today. It needs
  a new return type (`Result`/enum) or the dedup must happen in each `*_many` entry
  point.
- **G-061 option (b), "build the mask with XOR instead of OR" — labelled
  *preferred* by the first agent — is wrong for four gate families.** `build_masks`
  is one shared helper for x/y/z/h/s/s_dag/sqrt_x(_dag)/sqrt_y(_dag); XOR would make
  `s_many(&[q,q])` skip q entirely (identity) when the truth is Z-conjugation, i.e.
  it replaces one wrong answer with a different wrong answer, since S² = Z ≠ I. The
  first agent does flag this, but it must not be labelled preferred. **The complete
  fix is the scalar fallback (item 1 above).**
- **`ppvm-stim/src/validate.rs` must NOT reject duplicate targets.** Stim accepts
  `X 0 0` and applies X per target, so rejecting it would break valid input;
  expand-or-apply-per-target is the only compatible option.
- **G-054's "i128 accumulation in `norm_sq`, then a saturating readout"** is not
  clean: `norm_sq` returns i64, so saturating swaps a debug panic + release wrap for
  a silently saturated "exact" integer norm — a new lie on a type whose whole point
  is exactness. Use `hypot` for `magnitude` instead (item 4 above).
- **G-035/G-043-style `debug_assert!` on a coefficient-generic `C`** (U1 item 7)
  will not compile as written — the same caveat applies to any generic assert here.

Both agents note that existing tests and benches all pass strictly increasing
indices and offset ≥ count, so nothing in-tree breaks on the asserts; only benches
call `cz_block*`, so no golden master moves for G-060.

**Sign-off needed.** G-061's fix and G-060's release-mode fallback change observable
behaviour (including sampled measurement records from legal `.stim` files) and make
`-2` diverge from legacy unless legacy is fixed in tandem; G-054's `hypot` change
alters results exactly where the current code already wraps.

## Ledger

`Rust cite` paths are relative to `crates/`; `Lean cite` paths to `lean/`.
IDs are assigned in a deterministic order — sector name, then class, then
severity high→low — and are never reused.

| ID | Sector | Class | Tier | Sev | Rust cite | Lean cite | Unverified claim | Status | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| G-001 | clifford-conjugation | bridge | A | medium | `ppvm-lossy-pauli-word-2/src/clifford.rs:59` | `PPVM/Pauli/Symplectic.lean` `sActL_cnotActL_sActL_eq_cyActL` | The blanket lossy `cy` cancels correctly when the control is lost and the target present. | open | — |
| G-002 | clifford-conjugation | bridge | A | medium | `ppvm-pauli-word-2/src/clifford.rs:125` | `PPVM/Pauli/Symplectic.lean` `cyAct` | Bare `PauliWord`'s own `cy`/`sqrt_*` bit maps are the Lean `Sp(2n,2)` maps. | open | — |
| G-003 | clifford-conjugation | bridge | A | low | `ppvm-phased-pauli-word-2/src/clifford.rs:66` | `PPVM/Pauli/Conjugation.lean` `conjX` | `Phased::x/y/z(q)` are literal group conjugation at site `q`. | open | — |
| G-004 | clifford-conjugation | coverage | A | medium | `ppvm-pauli-word-2/src/clifford.rs:102` | none | The `√X` transvection and `cyAct` are bijections on `Sp n`, so a Clifford re-key never collides. | open | — |
| G-005 | clifford-conjugation | strength | A | medium | `ppvm-phased-pauli-word-2/src/clifford.rs:134` | `PPVM/Pauli/Conjugation.lean` `conjCNOT_sign` | Two-site CNOT/CZ/CY signs are the correct conjugation phase inside 𝒫ₙ for n > 2. | open | — |
| G-006 | graded-algebra-containers | bridge | B | medium | `ppvm-traits-2/src/containers/hash_join.rs:181` | `PPVM/Algebra/Twisted.lean:213` `twistedConv` | The two shipped container `Multiply` impls compute the twisted convolution and agree. | open | — |
| G-007 | graded-algebra-containers | fidelity | B | low | `ppvm-traits-2/src/graded.rs:104` | `PPVM/Algebra/GradedMap.lean:261` `reduce_structural` | `reduce` drops exactly the zeros and `Support::len` is the canonical support size. | open | — |
| G-008 | graded-algebra-containers | fidelity | B | low | `ppvm-traits-2/src/coefficient.rs:141` | `PPVM/Instantiations/Projector.lean:25` | `Halvable::half` is exact: `x.half() + x.half() == x` for the shipped `f64`/`Complex<f64>`. | adjudicated-spec | [adj U6](#u6-unenforced-preconditions) |
| G-009 | hashing-digests | coverage | B | medium | `ppvm-pauli-word-2/src/hash.rs:46` | `PPVM/Pauli/Word.lean:45` (`Word`) | `key_hash`/`PartialEq`/`weight` depend only on the logical sites, not on padding. | open | — |
| G-010 | hashing-digests | coverage | B | medium | `ppvm-traits-2/src/hash.rs:28` | none | The three clauses of the `Indexable` pass-through contract hold for every shipped key. | open | — |
| G-011 | hashing-digests | coverage | B | medium | `ppvm-lossy-pauli-word-2/src/data.rs:433` | none | Logically equal lossy words are `Eq` and share a digest (requires canonical loss). | adjudicated-spec | [adj U3](#u3-lossy-canonicality) |
| G-012 | hashing-digests | provenance | B | medium | `ppvm-tableau-2/src/mixture/noise/pauli.rs:13` | none | The incremental fingerprint XOR delta equals the from-scratch fingerprint. | open | — |
| G-013 | lossy-word | bridge | A | high | `ppvm-lossy-pauli-word-2/src/clifford.rs:117` | `PPVM/Pauli/Symplectic.lean` `cnotActL_preserves_loss` | After any Clifford sequence every lost site still has `x = z = 0`. | adjudicated-spec | [adj U3](#u3-lossy-canonicality) |
| G-014 | lossy-word | coverage | B | medium | `ppvm-lossy-pauli-word-2/src/data.rs:421` | none | Logical lossy words biject with canonical `(x,z,loss)` triples and every mutator preserves that. | adjudicated-spec | [adj U3](#u3-lossy-canonicality) |
| G-015 | lossy-word | coverage | B | low | `ppvm-lossy-pauli-word-2/src/data.rs:288` | none | Counting `Lost` as non-identity is the right grading for a truncated `LossyPauliSum`. | open | — |
| G-016 | lossy-word | fidelity | A | high | `ppvm-lossy-pauli-word-2/src/clifford.rs:102` | `PPVM/Pauli/Symplectic.lean:314` `variable (lost …)` | The loss plane is word state, is never written by a gate, and the gate map is key-injective. | adjudicated-spec | [adj U3](#u3-lossy-canonicality) |
| G-017 | measurement-branching | bridge | A | high | `ppvm-tableau-2/src/measure.rs:407` | `PPVM/Tableau/Projection.lean:116` `rustTerm_eq` | The real `z_overlap_re`/`prob_1`/case-a merge/case-b retain compute what `Projection.lean` says. | open | — |
| G-018 | measurement-branching | coverage | A | high | `ppvm-tableau-2/src/mixture/measure.rs:86` | none | `for_each_z_branch` unravels a measurement into a correct classical mixture. | adjudicated-defect | [adj U4](#u4-mixture-weights-and-sampler) |
| G-019 | measurement-branching | coverage | A | high | `ppvm-tableau-2/src/mixture/sampler.rs:43` | none | `MixtureSampler::sample` draws shots from the Born distribution. | adjudicated-defect | [adj U4](#u4-mixture-weights-and-sampler) |
| G-020 | measurement-branching | coverage | B | medium | `ppvm-tableau-2/src/mixture/equality.rs:19` | none | `structurally_equal` is the right merge key, so merging mixture entries is exact. | adjudicated-defect | [adj U4](#u4-mixture-weights-and-sampler) |
| G-021 | measurement-branching | strength | A | high | `ppvm-tableau-2/src/data.rs:496` | `PPVM/Tableau/BranchPhase.lean:224` `FrameInvolution` | The reported deterministic outcome is the true eigenvalue, not its negation. | adjudicated-spec | [adj U5](#u5-measurement-sign) |
| G-022 | measurement-branching | strength | A | high | `ppvm-tableau-2/src/measure.rs:75` | `PPVM/Tableau/Projection.lean:225` `probOne_eq` | The anticommuting branch is a fair coin, and the two samplers agree. | adjudicated-spec | [adj U5](#u5-measurement-sign) |
| G-023 | measurement-branching | strength | A | high | `ppvm-tableau-2/src/measure.rs:449` | `PPVM/Tableau/Projection.lean:260` `projectRaw_eq_two_proj` | The case-a merge plus frame update leaves the normalized projection `P_b\|ψ⟩`. | adjudicated-spec | [adj U5](#u5-measurement-sign) |
| G-024 | multiply-rotation | bridge | A | high | `ppvm-pauli-sum-2/src/rotation.rs:353` | `PPVM/Instantiations/Rotation.lean:553` `comm2_generic_sign_eq_branchExp2` | The whole two-qubit rotation surface computes what the two-site Lean section proves. | open | — |
| G-025 | multiply-rotation | bridge | A | medium | `ppvm-pauli-sum-2/src/rotation.rs:266` | `PPVM/Instantiations/Rotation.lean:105` `rz_eps_from_product` | The absolute branch sign ε of each axis column is the Lean value. | adjudicated-spec | [adj U2](#u2-rotation-direction-and-sign) |
| G-026 | multiply-rotation | bridge | B | medium | `ppvm-pauli-sum-2/src/column_store/graded.rs:158` | `PPVM/Algebra/Twisted.lean:213` `twistedConv` | The columnar and indexmap products compute `twistedConv` and keep columns aligned. | open | — |
| G-027 | multiply-rotation | bridge | A | low | `ppvm-pauli-sum-2/src/multiply.rs:277` | `PPVM/Algebra/Twisted.lean:250` `twistedConv_single_right` | `mul_word_assign` equals the general single-term product. | open | — |
| G-028 | multiply-rotation | coverage | A | medium | `ppvm-pauli-sum-2/src/rotation.rs:75` | none | `levi_civita(i,j)` is `−i[P_i,P_j]/2 = ε·P_k` in the crate's argument order. | adjudicated-spec | [adj U2](#u2-rotation-direction-and-sign) |
| G-029 | multiply-rotation | provenance | A | medium | `ppvm-pauli-sum-2/src/rotation.rs:228` | `PPVM/Instantiations/Rotation.lean:90` `rx_eps_from_product` | The shipped rotation direction is correct (`rx(θ)` is not `rx(−θ)`). | adjudicated-spec | [adj U2](#u2-rotation-direction-and-sign) |
| G-030 | multiply-rotation | strength | A | high | `ppvm-pauli-sum-2/src/rotation.rs:9` | `PPVM/Instantiations/Rotation.lean:17,:75` `branchExp` | The rotation conjugation identity the ε columns are built on. | adjudicated-spec | [adj U2](#u2-rotation-direction-and-sign) |
| G-031 | multiply-rotation | strength | B | medium | `ppvm-pauli-sum-2/src/column_store/rotations/rx.rs:127` | `PPVM/Instantiations/Rotation.lean:388` `accumulate_rotBatch` | The columnar `rx` keeps the two-pass ordering `accumulate_rotBatch` licenses. | open | — |
| G-032 | multiply-rotation | strength | A | medium | `ppvm-pauli-sum-2/src/rotation.rs:248` | `PPVM/Instantiations/Rotation.lean:275` `rot_norm_sq` | The branch is the 2-D rotation `rot`, so norm preservation and angle additivity apply. | open | — |
| G-033 | noise-observables | bridge | A | high | `ppvm-pauli-sum-2/src/noise.rs:178` | `PPVM/Algebra/Noise.lean:235` `twoQubitPauliError_indices_anticommuting` | Each of the 16 arms of `two_qubit_pauli_error` scales by its own λ_P. | open | — |
| G-034 | noise-observables | bridge | A | high | `ppvm-pauli-sum-2/src/loss.rs:137` | `PPVM/Algebra/Noise.lean:450` `correlatedLossChannel_trace_preserving` | The shipped loss-channel arms realize `corrT`/`lossT`/`resetT`. | open | — |
| G-035 | noise-observables | bridge | A | medium | `ppvm-tableau-2/src/mixture/noise/pauli.rs:96` | `PPVM/Algebra/Noise.lean:150` `eigenvalue_abs_le_one_needs_substochastic` | Channel inputs are sub-stochastic, so the eigenvalue is contractive. | adjudicated-defect | [adj U1](#u1-correlated-loss-convention) |
| G-036 | noise-observables | coverage | A | high | `ppvm-pauli-sum-2/src/noise.rs:312` | none | `amplitude_damping(q, γ)` is the Heisenberg adjoint of the standard damping channel. | open | — |
| G-037 | noise-observables | coverage | A | high | `ppvm-tableau-2/src/mixture/noise/loss.rs:32` | none | The mixture's `loss_channel` is the same channel as the trajectory one. | open | — |
| G-038 | noise-observables | coverage | A | medium | `ppvm-tableau-2/src/noise.rs:74` | none | The trajectory samplers realize the channel probabilities exactly. | open | — |
| G-039 | noise-observables | coverage | B | low | `ppvm-pauli-sum-2/src/noise.rs:220` | none | `1 − 4p/3` and `1 − 16p/15` are the depolarizing transfer eigenvalues. | open | — |
| G-040 | noise-observables | fidelity | A | high | `ppvm-tableau-2/src/mixture/noise/loss.rs:95` | `PPVM/Algebra/Noise.lean:422` `corrT` | `corrT` is the spec for all three backends' correlated loss channel. | adjudicated-defect | [adj U1](#u1-correlated-loss-convention) — **RULED: paper normative**; fix `ppvm-tableau-2` both backends, `corrT` and pauli-sum-2 stand |
| G-041 | noise-observables | fidelity | B | medium | `ppvm-tableau-2/src/noise.rs:270` | `PPVM/Algebra/Noise.lean:289` `fire_nonstrict_fires_at_zero` | `loss_channel(q, 0.0)` is not the identity under the shipped convention. | open | — |
| G-042 | noise-observables | strength | A | high | `ppvm-pauli-sum-2/src/noise.rs:73` | `PPVM/Algebra/Noise.lean:57` `pauli_channel_eigenvalue_omega` | λ_P is the transfer eigenvalue of the unital Pauli channel. | open | — |
| G-043 | noise-observables | strength | A | medium | `ppvm-pauli-sum-2/src/loss.rs:82` | `PPVM/Algebra/Noise.lean:450` `correlatedLossChannel_trace_preserving` | The loss channels are CPTP, not merely trace-preserving linear maps. | adjudicated-defect | [adj U1](#u1-correlated-loss-convention) — region ruled: `p₀ + 2p₁ ≤ 1` |
| G-044 | noise-observables | strength | A | medium | `ppvm-pauli-sum-2/src/trace.rs:93` | `PPVM/Algebra/Noise.lean:490` `overlap_with_zero_xfree` | `⟨0ⁿ\|P\|0ⁿ⟩ = [P is X-free]`, so the zero-state read-out is right. | open | — |
| G-045 | products-and-channels (skeptic) | coverage | A | medium | `ppvm-tableau-2/src/noise.rs:317` | none | `asymmetric_loss_channel`'s state-dependent `p_tot` and its omitted back-action. | adjudicated-spec | [adj U1](#u1-correlated-loss-convention) |
| G-046 | sum-engine-stores | bridge | A | medium | `ppvm-pauli-sum-2/src/column_store/lifecycle.rs:151` | `PPVM/Algebra/GradedMap.lean:200` `pushforward_eq_reset_accumulate` | `apply_producer` computes the pushforward and the `reset` is not optional. | open | — |
| G-047 | sum-engine-stores | coverage | B | medium | `ppvm-pauli-sum-2/src/column_store/columns.rs:7` | `PPVM/Algebra/GradedMap.lean:51` (`CMap`) | Each store has a representation invariant and an abstraction map that commutes with every op. | open | — |
| G-048 | sum-engine-stores | coverage | B | low | `ppvm-pauli-sum-2/src/column_store/graded.rs:92` | none | `probe_batch(keys, out)` sets `out[i] = get(keys[i])` pointwise. | open | — |
| G-049 | sum-engine-stores | fidelity | A | high | `ppvm-pauli-sum-2/src/store.rs:910` | `PPVM/Algebra/GradedMap.lean:64` (`len`) | `Support::len` and `Sum`'s equality are the canonical zero-free Finsupp ones. | open | — |
| G-050 | sum-engine-stores | fidelity | A | high | `ppvm-pauli-sum-2/src/column_store/rotations/rx.rs:127` | `PPVM/Instantiations/Rotation.lean` `accumulate_rotBatch` | Every backend's rotation is the two-pass walk `twoPass`. | open | — |
| G-051 | sum-engine-stores | provenance | B | medium | `ppvm-pauli-sum-2/src/indexmap_store/branching.rs:23` | none | The ordered backend's term order and its dedup-cardinality merge rule are correct. | open | — |
| G-052 | symbolic-coefficients | bridge | A | medium | `ppvm-sym-2/src/mul.rs:227` | `PPVM/Instantiations/Symbolic.lean:467` `evalC_mul` | `eval`/`eval_complex` are multiplicative on the map-backed product surface. | open | — |
| G-053 | symbolic-coefficients | bridge | A | medium | `ppvm-sym-2/src/term.rs:376` | `PPVM/Instantiations/Symbolic.lean:173` `mulMono_drop_at_insert_eq_drop_at_end` | Dropping each monomial at insert equals truncating the full product. | open | — |
| G-054 | symbolic-coefficients | coverage | A | low | `ppvm-sym-2/src/exact.rs:22` | none | ℤ[i] admits no `half`, which is why `Halvable` was split off `Coefficient`. | adjudicated-defect | [adj U6](#u6-unenforced-preconditions) |
| G-055 | symbolic-coefficients | fidelity | A | high | `ppvm-sym-2/src/mul.rs:235` | `PPVM/Instantiations/Symbolic.lean:882` `mulImpl` | `mulImpl` is the product `ppvm-sym-2` actually implements. | open | — |
| G-056 | symbolic-coefficients | fidelity | A | medium | `ppvm-sym-2/src/eval.rs:138` | `PPVM/Instantiations/Symbolic.lean:260` `evalHom` | `evalHom` on the phase-free ring models `Term::eval`. | open | — |
| G-057 | symbolic-coefficients | strength | A | medium | `ppvm-sym-2/src/term.rs:323` | `PPVM/Algebra/GradedMap.lean:95` `accumulate_assoc` | `Term`'s `+` is an additive monoid operation. | adjudicated-spec | [adj U6](#u6-unenforced-preconditions) |
| G-058 | tableau-and-symbolic (skeptic) | bridge | A | medium | `ppvm-tableau-2/src/data.rs:318` | `PPVM/Pauli/Word.lean` `phaseExpN` | `Row::mul_assign`'s bits and ℤ/4 phase are `mulWord` and `phaseExpN`. | adjudicated-spec | [adj U5](#u5-measurement-sign) |
| G-059 | tableau-core | bridge | A | high | `ppvm-tableau-2/src/clifford.rs:491` | `PPVM/Tableau/Batch.lean:216` `isSitewise_extSqrtY` | The six tableau extension gates rewrite every row's ℤ/4 phase by the audited table. | adjudicated-spec | [adj U5](#u5-measurement-sign) |
| G-060 | tableau-core | bridge | B | medium | `ppvm-tableau-2/src/data.rs:1474` | `PPVM/Tableau/Batch.lean:260` `czSeq_phase` | The fused CZ-block parity phase equals the sequential per-pair loop. | adjudicated-defect | [adj U6](#u6-unenforced-preconditions) |
| G-061 | tableau-core | bridge | B | medium | `ppvm-tableau-2/src/clifford.rs:594` | `PPVM/Tableau/Batch.lean:91` `seqApply_eq_batchApply` | Every `*_many` fused sweep equals the per-index gate loop. | adjudicated-defect | [adj U6](#u6-unenforced-preconditions) |
| G-062 | tableau-core | coverage | A | high | `ppvm-tableau-2/src/data.rs:492` | `PPVM/Tableau/Frame.lean:53` `IsSymplecticFrame` | Every row of a live tableau carries a real phase, preserved by every operation. | adjudicated-spec | [adj U5](#u5-measurement-sign) |
| G-063 | tableau-core | coverage | A | medium | `ppvm-tableau-2/src/gates.rs:82` | `PPVM/Instantiations/Rotation.lean:654` `rotXY_heisenberg_order` | The tableau's `r(q,φ,θ)` and `u3` are the intended rotations in the forward order. | open | — |
| G-064 | tableau-core | coverage | A | low | `ppvm-tableau-2/src/gates.rs:28` | none | `t()`/`t_dag()` branch by the coefficient pair of `T = diag(1, e^{iπ/4})`. | open | — |
| G-065 | tableau-core | strength | A | high | `ppvm-tableau-2/src/data.rs:548` | `PPVM/Tableau/Frame.lean:296` `isSymplecticFrame_projectFrame` | The projection maps the frame of \|ψ⟩ to a frame of the projected state, sign included. | open | — |
| G-066 | tableau-core | strength | A | medium | `ppvm-tableau-2/src/clifford.rs:358` | `PPVM/Pauli/Conjugation.lean:85` `conjH_sign` | The four base conjugation tables are the signs of genuine `G·P·G†`. | open | — |
| G-067 | truncation-policy-loss | bridge | A | high | `ppvm-pauli-sum-2/src/policy.rs:168` | `PPVM/Algebra/Truncation.lean:44` `l1_bound` | The real Rust truncation incurs at most the ℓ¹ mass of the terms it dropped. | open | — |
| G-068 | truncation-policy-loss | bridge | B | medium | `ppvm-tableau-2/src/data.rs:1240` | `PPVM/Algebra/Truncation.lean:237` `cutoff_mismatch` | The two backends' `≥`-vs-`>` cutoff split is machine-checked. | open | — |
| G-069 | truncation-policy-loss | bridge | B | medium | `ppvm-pauli-sum-2/src/policy.rs:258` | `PPVM/Algebra/GradedMap.lean:600` `retain_seq_eq_retain_and` | `CombinedPolicy` is the conjunction, order-independent, and the sentinel skip is exact. | open | — |
| G-070 | truncation-policy-loss | bridge | B | medium | `ppvm-pauli-sum-2/src/sum.rs:488` | `PPVM/Algebra/Truncation.lean:315` `truncate_preserve_eq_widened_retain` | The snapshot/policy/restore composite equals one widened retain pass. | open | — |
| G-071 | truncation-policy-loss | coverage | A | high | `ppvm-tableau-2/src/mixture/data.rs:142` | none | The mixture's drop-then-renormalize truncation is a bounded approximation. | adjudicated-spec | [adj U4](#u4-mixture-weights-and-sampler) |
| G-072 | truncation-policy-loss | fidelity | A | high | `ppvm-tableau-2/src/data.rs:1240` | `PPVM/Algebra/Truncation.lean:219` `l2_bound` | `l2_bound` is the bound the stabilizer-tableau amplitude cutoff uses. | open | — |
| G-073 | truncation-policy-loss | fidelity | B | medium | `ppvm-lossy-pauli-word-2/src/data.rs:288` | `PPVM/Pauli.lean:73` `Pauli.weight` | The weight `MaxPauliWeight` thresholds is the count of non-identity Pauli slots. | open | — |
| G-074 | truncation-policy-loss | strength | A | high | `ppvm-pauli-sum-2/src/policy.rs:30` | `PPVM/Algebra/Truncation.lean:44` `l1_bound` | Each shipped policy has an error guarantee in its own parameter. | open | — |
| G-075 | truncation-policy-loss | strength | A | high | `ppvm-pauli-sum-2/src/sum.rs:629` | `PPVM/Algebra/Noise.lean:117` `l1_contractive` | The per-truncation ℓ¹ error telescopes over a noisy circuit. | open | — |
| G-076 | truncation-policy-loss | strength | B | medium | `ppvm-traits-2/src/graded.rs:99` | `PPVM/Algebra/GradedMap.lean:149` `accumulateTerms_perm` | A batch may be reordered and partitioned freely, exactly. | open | — |
| G-077 | word-algebra | bridge | B | medium | `ppvm-pauli-word-2/src/product.rs:58` | `PPVM/Pauli/Word.lean:55` `phaseExpN` | `key_mul` equals `phaseExpN` for words wider than one storage limb. | open | — |
| G-078 | word-algebra | coverage | B | medium | `ppvm-pauli-word-2/src/data.rs:38` | none | The canonical-unused-bits / tail-limb padding invariant. | open | — |
| G-079 | word-algebra | coverage | B | medium | `ppvm-pauli-word-2/src/data.rs:307` | `PPVM/Pauli.lean:73` `Pauli.weight` | `PauliWord::weight()` is the number of non-identity sites. | open | — |
| G-080 | word-algebra | coverage | B | low | `ppvm-pauli-word-2/src/data.rs:364` | none | `pauli_code(i)` equals `x_bit \| (z_bit << 1)` for every storage and every `i`. | open | — |
| G-081 | word-algebra | fidelity | A | low | `ppvm-traits-2/src/word.rs:47` | `PPVM/Pauli.lean:15`, `PPVM.Pauli.mul` | The Lean sector modules model the `-2` crates; the discriminant table is right. | open | — |
| G-082 | word-algebra | strength | A | medium | `ppvm-phased-pauli-word-2/src/product.rs:44` | `PPVM/Pauli/Phase.lean:241` `mul_one'` | `Phased::mul` is the group operation of the n-qubit phased Pauli group. | open | — |
| G-083 | word-and-clifford (skeptic) | bridge | A | medium | `ppvm-pauli-sum-2/src/clifford.rs:350` | `PPVM/Pauli/Conjugation.lean` `extSqrtX_sign` | `Sum`'s own extension-gate closures emit the Lean conjugation signs. | open | — |

## Gap details

Each section carries the claim that is unverified, why the current Lean and
oracle artifacts do not establish it, and the closure — the Lean theorem to
state *and* the oracle test that would pin the Rust to it. A close round should
be able to work a row from this section alone.

### G-001 — LossyPauliWord's blanket CY loss-cancellation has no `*_lean.rs` bridge, only a legacy diff

- Class `bridge` · Tier A · Severity medium · Sector clifford-conjugation · Status **open**
- Rust: `crates/ppvm-lossy-pauli-word-2/src/clifford.rs:59`
- Lean: `lean/PPVM/Pauli/Symplectic.lean` `sActL_cnotActL_sActL_eq_cyActL`

**Claim.** `LossyPauliWord`'s `cy` — emitted by the blanket as `s(t); cnot(c,t);
s(t); z(t)` with per-primitive guards that do NOT all agree — equals the atomic
whole-gate skip `cyActL`, preserves `lost ⇒ x=z=0`, and is the Sp isometry on
the present block.

**Why unverified.** This is the one place where the `-2` port genuinely changes
the algorithm (legacy had an atomic `if lost(c)||lost(t) { return }` in `fn cy`;
the blanket at `traits-2/src/pauli.rs:297` has three primitives whose guards
differ, so a lost control with a present target must *cancel* rather than skip).
Lean proves exactly this (`sActL_cnotActL_sActL_eq_cyActL`,
`cyActL_preserves_loss`, `cyActL_present_isometry`). The Rust side has no
bridge: `lossy_pauli_word_lean.rs` never calls `.cy(`, `.sqrt_x(`, `.s_dag(` —
its only gate test (`lossy_clifford_generators_preserve_symplectic_form`, lines
320-339) loops `gate 0..4` over h/s/cnot/cz. The crate's own module doc
(`clifford.rs:44-53`) names its evidence: "verified against the old reference
over the full 25-word lossy alphabet in
`ppvm-conformance-2::lossy_pauli_word_diff`" — i.e. legacy agreement, which is
evidence of a faithful port and nothing about correctness. Additionally
`clifford_leaves_loss_mask_invariant` (line 416) checks only the loss mask after
a circuit; `assert_site_invariants` (the `lost ⇒ x=z=0` check) is run only on
freshly parsed words, so `cnotActL_lost_target_stays_identity` is bridged by a
single hand-written 2-qubit case in the crate (`"XL"`, `clifford.rs:230`).

**Proposed closure.** Add to `lossy_pauli_word_lean.rs`: (a)
`lossy_cy_equals_atomic_whole_gate_skip` — for all 25 two-site lossy words plus
randomized n-qubit masks, assert `w.cy(c,t)` equals the reference `cyActL`
(identity if either operand lost, else the `cyAct` bit rule), covering the
lost-control/present-target cancellation; (b) extend the isometry loop to
cy/sqrt_x/sqrt_x_dag/sqrt_y/sqrt_y_dag (`cyActL_present_isometry`); (c)
re-assert `assert_site_invariants` after every gate of a random circuit that
includes cy, pinning `*ActL_preserves_loss`.

### G-002 — Bare PauliWord's hand-written CliffordExtensions (cy, sqrt_*) pinned only by the legacy diff

- Class `bridge` · Tier A · Severity medium · Sector clifford-conjugation · Status **open**
- Rust: `crates/ppvm-pauli-word-2/src/clifford.rs:125`
- Lean: `lean/PPVM/Pauli/Symplectic.lean` `cyAct` (and `Conjugation.lean` `conjCY_bits`)

**Claim.** `PauliWord`'s own atomic `cy` bit rule (`z_c ⊕= x_t ⊕ z_t`,
`x_t ⊕= x_c`, `z_t ⊕= x_c`) and its
`sqrt_x`/`sqrt_x_dag`/`sqrt_y`/`sqrt_y_dag`/`s_dag` aliases realize the same
Sp(2n,2) bit maps as Lean's `cyAct` / `extSqrtX_bits` / `extSqrtY_bits`.

**Why unverified.** `PauliWord` deliberately opts out of `BlanketClifford`
(`clifford.rs:8`) and hand-writes a third independent copy of every bit map (its
public gates at `clifford.rs:36-139` are themselves a fourth copy relative to
its own `SymplecticColumns` primitives at `clifford.rs:141-206`). The oracle
`pauli_word_lean.rs::crate_clifford_bit_map_matches_conjugation_oracle` (line
641) covers only h/s/cnot/cz, and `clifford_generators_preserve_symplectic_form`
(line 712) only `gate 0..4` = h/s/cnot/cz. The in-crate unit tests
(`clifford.rs:238-309`) likewise stop at h/s/cnot/cz/x/y/z — there is no `cy` or
`sqrt_*` test in the crate at all. The only coverage is
`pauli_word_diff.rs:243-276` against the legacy crate. So a transposition inside
`PauliWord::cy` (e.g. writing `zt ^ xc` into `x_t`, or reading `zc` where `zt`
is meant) would be caught by nothing but agreement with legacy, even though the
correct Lean model (`cyAct`, `conjCY_bits`) already exists.

**Proposed closure.** Extend
`crate_clifford_bit_map_matches_conjugation_oracle` to the extension gates:
exhaustively over the 16 two-site words assert `PauliWord::cy` matches `conj_cy`
(the Rust mirror of Lean `conjCY_bits`, bits only), and over the 4 single-site
words assert sqrt_x/sqrt_x_dag/sqrt_y/sqrt_y_dag/s_dag match
`extSqrtX_bits`/`extSqrtY_bits`/`conjS_bits`; add cy/sqrt_x to the isometry
loop. Also assert the two internal paths agree (public `w.h(q)` vs
`SymplecticColumns::swap_xz(q)`, etc.), which nothing currently does.

### G-003 — Phased word's pure-sign x/y/z conjugation has no oracle test (only an n=1 unit test)

- Class `bridge` · Tier A · Severity low · Sector clifford-conjugation · Status **open**
- Rust: `crates/ppvm-phased-pauli-word-2/src/clifford.rs:66`
- Lean: `lean/PPVM/Pauli/Conjugation.lean` `conjX`

**Claim.** `Phased::x/y/z(q)` fix the word and advance the phase by 2 exactly
when `z_q` / `x_q ⊕ z_q` / `x_q`, i.e. they are literal group conjugation by
X/Y/Z at site `q`.

**Why unverified.** The Lean side is the strongest in the file —
`conjX`/`conjY`/`conjZ` are stated as `mul (mul G p) G` over the bundled group
product, so the signs are *derived*, not asserted. But
`phased_pauli_word_lean.rs` never calls `.x(`, `.y(` or `.z(`: its
matrix-grounding tests cover h/s (n=1, n≤4), cnot/cz (n=2) and the six
extensions (n=1/n=2) only, and `named_conjugation_sign_theorems` covers
h/cnot/cz/s. The only coverage is the in-crate unit test
`single_qubit_gates_track_sign` (`clifford.rs:272-294`), which is n=1 — so a
site-indexing mistake (reading site q+1, or the wrong operand of the loss guard
at `clifford.rs:67`) is invisible, and the Lean theorem's cited consumer is the
`ppvm-pauli-sum-2` fast path rather than this kernel.

**Proposed closure.** Add `pauli_conjugation_signs_grounded_in_zi_matrices` to
`phased_pauli_word_lean.rs`: exhaustive n=1 plus randomized n=3,4 with a random
target q, asserting `G_q M G_q† == phased_mat(w)` after
`w.x(q)`/`w.y(q)`/`w.z(q)` using the existing `embed_single`, and cite
`conjX`/`conjY`/`conjZ`.

### G-004 — No Sp n model for the √X transvection; cyAct bijectivity missing despite the re-key no-collision claim

- Class `coverage` · Tier A · Severity medium · Sector clifford-conjugation · Status **open**
- Rust: `crates/ppvm-pauli-word-2/src/clifford.rs:102`
- Lean: none

**Claim.** Every Clifford gate's phase-stripped bit map on `Sp n` is a bijection
(so a Clifford re-key never collides two terms) and an ω-isometry — including
the `√X`-family transvection `x ⊕= z` and `cyAct`.

**Why unverified.** `Symplectic.lean` models exactly four bit maps at the
n-qubit level: hAct (swap), sAct (`z ⊕= x`), cnotAct, czAct — with isometry +
involutive + bijective for each — plus cyAct with an isometry only. The `√X`
family's map is the *other* transvection `x ⊕= z`
(`pauli-word-2/src/clifford.rs:104`, phased `clifford.rs:202`, tableau
`clifford.rs:469`); it has no `sqrtXAct` definition, hence no isometry and no
bijectivity, and `cyAct` has no `cyAct_involutive`/`cyAct_bijective`. That
matters concretely: `ppvm-pauli-sum-2/src/clifford.rs:33` asserts "A Clifford
re-key is a bijection, so colliding re-keyed terms never occur; reduce is a
no-op on the support size" for its whole gate surface, which includes
`sqrt_x`/`sqrt_x_dag`/`cy` (own re-key closures at `clifford.rs:350/365/447`),
while the Lean bijectivity theorems the §230 docstring points to cover only
h/s/cnot/cz. The bridge has the same hole:
`pauli_sum_lean.rs::clifford_rekey_is_support_preserving_bijection` (line 342)
cycles `i % 5` over h/s/cnot/cz/z only.

**Proposed closure.** Add to `Symplectic.lean`:
`sqrtXAct q v = update v q ((v q).1 + (v q).2, (v q).2)` with
`sqrtXAct_isometry`, `sqrtXAct_involutive`, `sqrtXAct_bijective`; and
`cyAct_involutive`/`cyAct_bijective` (CY is Hermitian and self-inverse, so this
follows from `sAct_cnotAct_sAct_eq_cyAct`). Extend
`clifford_rekey_is_support_preserving_bijection` and the `*_lean.rs` isometry
loops to sqrt_x/sqrt_x_dag/sqrt_y/sqrt_y_dag/cy.

### G-005 — Two-qubit conjugation signs proven and matrix-tested only at exactly n=2

- Class `strength` · Tier A · Severity medium · Sector clifford-conjugation · Status **open**
- Rust: `crates/ppvm-phased-pauli-word-2/src/clifford.rs:134`
- Lean: `lean/PPVM/Pauli/Conjugation.lean` `conjCNOT_sign`

**Claim.** On an n-qubit phased word, the CNOT/CZ/CY sign computed from just the
two touched sites is the correct conjugation phase inside 𝒫ₙ (the untouched
qubits contribute nothing to the delta).

**Why unverified.** Lean's phase model tops out at two qubits: `TwoPauli` is
literally `phase × (xc,zc,xt,zt)` (`Conjugation.lean:496-507`) and
`conjCNOT`/`conjCZ`/`conjCY` are automorphisms of that 64-element group only.
There is no 𝒫ₙ conjugation object anywhere in `lean/PPVM` (`Word.lean` /
`Matrix.lean` give the n-qubit *product* phase `phaseExpN` and its tensor-matrix
grounding, but no n-qubit gate action), so the locality/embedding step — obvious
by tensor factorization, but exactly the kind of step this discipline exists to
discharge — is unproven. The oracle mirrors the same ceiling:
`two_qubit_conjugation_signs_grounded_in_zi_matrices` and
`extension_conjugation_signs_grounded_in_zi_matrices` build only `kron` of two
sites, and the one n-qubit matrix test,
`n_qubit_conjugation_phase_grounded_in_zi_matrices` (line 534), uses
`embed_single` and exercises H and S only — there is no `embed_two`.
`named_conjugation_sign_theorems` is n=2. `pauli_sum_lean.rs`'s n>2 cnot/cz
usage only checks support-size bijection, not signs. So for n>2 the two-qubit
sign half of conjugation rests on inspection of the Rust.

**Proposed closure.** Lean: define the n-qubit phased word action
`conjCNOTN c t : PhasedWord n →* PhasedWord n` (or embed 𝒫₂ into 𝒫ₙ via a
MonoidHom on the two touched factors) and prove
`phase (conjCNOTN c t p) = phase p + cnotDelta (p c) (p t)` — the locality
theorem. Oracle: add `embed_two` to `phased_pauli_word_lean.rs` and extend
`n_qubit_conjugation_phase_grounded_in_zi_matrices` to cnot/cz/cy at randomized
(c,t) for n = 3,4, asserting `G M G†` against the tracked phase.

### G-006 — Shipped container `Multiply` impls (L4) are exercised by no test in the repo

- Class `bridge` · Tier B · Severity medium · Sector graded-algebra-containers · Status **open**
- Rust: `crates/ppvm-traits-2/src/containers/hash_join.rs:181`
- Lean: `lean/PPVM/Algebra/Twisted.lean:213` `twistedConv` / `twistedConv_single_single`

**Claim.** The two shipped `Multiply::multiply_into` impls compute the twisted
convolution `(A·B)[k] = Σ_{p·q=k} A[p]·B[q]·i^{β(p,q)}` of
`Twisted.twistedConv`, and agree with each other.

**Why unverified.** `rg 'multiply_into|Multiply' crates/ppvm-traits-2/{src,tests}`
returns only the two impls themselves (`hash_join.rs:181`,
`coordinate_list.rs:149`) plus a LOCAL newtype `ZiSum` at
`crates/ppvm-traits-2/tests/phase1_leaf_types.rs:1055` which re-implements the
coordinate-list algorithm by hand (via `acc.accumulate`, not the shipped inline
find/push) and tests that copy. That test's header at
`phase1_leaf_types.rs:996` asserts "The crate ships no L4 impl (`containers.rs`
stops at L3 + `Retain`)" — false, and doubly stale since `containers.rs` was
split by backend in `d92421b6`. Neither `containers/tests.rs` nor
`tests/phase1_containers.rs` imports `Multiply`. Nothing downstream uses these
impls either (`ppvm-pauli-sum-2` has its own
`Multiply for HashMapStore/ColumnStore/IndexMapStore`), so the entire workspace
test suite is invariant under mutating them: deleting the `phase.apply(...)`
fold at `hash_join.rs:185`, or swapping `and_modify/or_insert` for a
term-dropping `insert`, breaks nothing.

**Proposed closure.** Add a cross-backend L4 test to
`crates/ppvm-traits-2/tests/phase1_containers.rs` over the existing exact
`Zi`/`PauliKey` stub: for all 16 key pairs and BOTH backends assert
`multiply_into` equals `Twisted.twistedConv_single_single` (X·Z = −i·Y etc.),
assert accumulation into a non-empty `acc` (the L4 contract), assert no implicit
`reduce` (the cancelled zero key survives, per `twistedConv`), and assert the
`Vec` and `HashMap` results agree as multisets. The Lean side needs nothing new
for the Pauli key; see G-047 for the abstract key.

### G-007 — `reduce` and `Support::len` have no Lean model; `reduce_structural` is Mathlib's `mem_support_iff`

- Class `fidelity` · Tier B · Severity low · Sector graded-algebra-containers · Status **open**
- Rust: `crates/ppvm-traits-2/src/graded.rs:104`
- Lean: `lean/PPVM/Algebra/GradedMap.lean:261` `reduce_structural` (also `:64` `len`)

**Claim.** `graded.rs:87-88`: "`reduce` drops exactly the zero coefficients
(`reduce_structural`)"; `GradedMap.lean:64` docstring: "`Support::len` — the
size of the canonical (zero-free) support".

**Why unverified.** There is no Lean definition of `reduce` at all.
`reduce_structural : f k ≠ 0 ↔ k ∈ f.support` is `Finsupp.mem_support_iff.symm`
— a property of the MODEL, not of the algorithm, unfalsifiable by any mutation
of the Rust `reduce` (`hash_join.rs:72`, `coordinate_list.rs:73`). Worse, the
model provably cannot express the state `reduce` acts on: the design mandates
zeros persist until finalize (`graded.rs:104` "run **only** at finalize"),
`Scale` deliberately leaves them, and both containers report `len()` = stored
entries, whereas Lean's `len` = `f.support.card` = nonzero count. They disagree
on exactly the interesting states: `tests/phase1_containers.rs:152` asserts
`Support::len == 3` for a map whose Finsupp image has support card 1, and
`pauli_sum_lean.rs:397` asserts `len == 1` after `scale(&0.0)` where the model
gives 0. So `len`, `is_empty`, the cardinality of `iter`, and `reduce` itself
are all outside the model while being credited to it.

**Proposed closure.** Model the un-reduced backend explicitly —
`Store K C := List (K × C)` (or a `Finsupp` paired with a stored-key `Finset`)
with `toFinsupp` summing duplicates — define `reduceStore` as the zero-filter
and `lenStore` as stored length, then prove
`toFinsupp (reduceStore s) = toFinsupp s` (reduce is semantically the identity),
`lenStore (reduceStore s) = (toFinsupp s).support.card`, and idempotence.
Oracle: extend `tests/phase1_containers.rs:148` to both backends, asserting the
pre-reduce `len` equals the stored count while `get` of every zeroed key is
`Some(0)` rather than `None`, and that a second `reduce` is a no-op.

### G-008 — `Halvable::half`'s stated exactness law is false for the shipped `f64` / `Complex<f64>` impls

- Class `fidelity` · Tier B · Severity low · Sector graded-algebra-containers · Status **adjudicated-spec** (round 2 — [adj U6](#u6-unenforced-preconditions))
- Rust: `crates/ppvm-traits-2/src/coefficient.rs:141`
- Lean: `lean/PPVM/Instantiations/Projector.lean:25` (`projLin`, halving by the exact ring ½)

**Claim.** `coefficient.rs:141`: "Impls must be exact: `x.half() + x.half() ==
x`" — the law used to argue that exact rings are excluded from `Halvable` while
the float and `Complex<f64>` domains that back Phase-1 measurement satisfy it.

**Why unverified.** `impl Halvable for f64 { *self / 2.0 }`
(`coefficient.rs:237`) violates the law at the bottom of the exponent range,
verified numerically: for `x = f64::from_bits(1)` (2^-1074), `x/2` rounds to
`0.0` and `0.0 + 0.0 = 0.0 ≠ x`; for `x = f64::from_bits(3)`, `x/2` rounds UP to
2^-1073 so doubling gives 2^-1072 ≠ x. `Complex<f64>` (`coefficient.rs:271`)
inherits it componentwise. The Lean side halves by the ring constant ½ over ℝ
(`Projector.lean`'s `projLin`), where the law holds — so the model is exact
precisely where the shipped impl is not, and the trait's own justification for
splitting `Halvable` out ("an exact ring could only satisfy it with a lossy
integer `/2` for which `half(x) + half(x) != x`") applies verbatim to `f64`. No
test checks the law on any impl.

**Round-2 amendment.** The law is exact **iff `x/2` is representable**: for all
|x| ≥ 2·`f64::MIN_POSITIVE` = 2^−1021, for ±0 and ±inf, and below 2^−1021 iff the
raw bit pattern is even (verified exhaustively over 4×2^17 patterns — every failure
had an odd pattern and bits < 2<<52; none at or above 2^−1021). Reachability: none.
The only call sites in the `-2` tree are `C::one().half()` at
`ppvm-pauli-sum-2/src/proj.rs:74` and `:92`, i.e. x = 1, where the law holds. See
[adj U6](#u6-unenforced-preconditions).

**Proposed closure — INVALIDATED by round 2.** The proposed weakening ("exact
wherever `x` is normal and `2·x` is representable") is **itself false**:
`f64::from_bits(0x10000000000001)` = 2^−1022(1 + 2^−52) is normal with 2x finite
and fails the law, so an oracle asserting the law "on normals" would fail on
correct code. Replacement: document "exact wherever `x/2` is representable — every
`|x| >= 2·f64::MIN_POSITIVE` (2^−1021), plus 0 and ±inf; below that iff the raw bit
pattern is even", record `from_bits(1)`, `from_bits(3)` **and**
`from_bits(0x10000000000001)` as the witnesses, and pin the oracle in
`tests/phase1_leaf_types.rs` to the 2^−1021 boundary (not to "normal") for both
`f64` and `Complex<f64>`. If the strong law is wanted, state the projector kernel
in Lean over a field with an explicit representability hypothesis rather than over
ℝ.

### G-009 — Digest/Eq/weight fold the whole backing blob; padding-invariance is proved nowhere

- Class `coverage` · Tier B · Severity medium · Sector hashing-digests · Status **open**
- Rust: `crates/ppvm-pauli-word-2/src/hash.rs:46`
- Lean: `lean/PPVM/Pauli/Word.lean:45` (`abbrev Word` — models a word as `Fin n → Bool × Bool`, so unused capacity is inexpressible)

**Claim.** For a packed word of width `nqubits` in `A`-wide storage,
`key_hash()` (and `PartialEq`, and `Word::weight`) depend only on the logical
sites `0..nqubits`, i.e. the digest is invariant under any state of the bits at
positions `nqubits..8*size_of::<A>()`.

**Why unverified.** `structural_hash` hashes
`bytemuck::bytes_of(x)`/`bytes_of(z)` — the full `size_of::<A>()` blob, every
padding bit included — and `PartialEq` (`data.rs:537`) compares
`xbits.data == other.xbits.data` on the whole blob, and `Word::weight`
(`data.rs:307-337`) popcounts the whole blob. The design requires the opposite
(`word-data-structures.md:133` "Equality and hashing exclude unused capacity";
`traits-2-configuration-and-hashing.md:2449` "It excludes RNG, padding, cache
state"). The code satisfies that only via the *canonical-unused-bits invariant*
asserted in a docstring at `data.rs:36-42`; the enforcement is
`debug_assert!(i < self.nqubits)` in every setter, which compiles out in
release. The same invariant is load-bearing for `key_mul`'s sign/imag popcounts
over full machine words (`product.rs:58-70`) and for `Tableau::key_hash`, which
hashes every row's full planes (`tableau-2/data.rs:608-611`) with 6 qubits in
`[usize;2]` = 122 padding bits in the shipped test.
`word-data-structures.md:454` lists "tests proving unused high bits do not
affect identity" as required prototype validation; grep finds no such test in
any `-2` crate or in `ppvm-conformance-2`, and no Lean statement of the
invariant exists.

**Proposed closure.** Lean: model the packed word as `bits : Fin W → Bool × Bool`
plus `nqubits ≤ W` with `Canonical w := ∀ i ≥ nqubits, w i = (false,false)`, and
prove `digest`/`eq`/`weight` (defined as folds over all W slots, as the Rust
does) agree with the logical `Fin n` versions under `Canonical`, plus closure
lemmas `Canonical` under set_bit/toggle/xor/mask-kernels. Oracle: a
crate-internal test in `ppvm-pauli-word-2` (and `ppvm-tableau-2`) that, for
every mutator and Clifford/product kernel and every storage tier, asserts the
raw plane words above `nqubits` are still zero, plus a `*_lean.rs` test that a
narrow word in wide storage has the digest/weight of its logical content.

### G-010 — No Lean leg at all for the Indexable / pass-through hashing contract

- Class `coverage` · Tier B · Severity medium · Sector hashing-digests · Status **open**
- Rust: `crates/ppvm-traits-2/src/hash.rs:28`
- Lean: none

**Claim.** The three clauses of the `Indexable` contract — `Hash for K` is
exactly `write_u64(key_hash())`; `a == b ⇒ a.key_hash() == b.key_hash()`;
`KeyColumn::hash_into` reproduces `key_hash()` bit for bit — hold for every
shipped key, and clause 2 is what makes
`HashMap<K, C, IdentityBuildHasher>` sound at all.

**Why unverified.** `lean/PPVM/**` contains no module for hashing: grep for
`key_hash`, `Indexable`, `digest`, `IdentityHasher` over all 19 `.lean` files
returns nothing but design-doc URLs. Because `IdentityBuildHasher` makes
`finish() == key.key_hash()` verbatim (`hash.rs:39-52`), hashbrown's soundness
precondition reduces *exactly* to clause 2 — there is no second hasher to
launder a violation — yet the only justification anywhere is prose plus
hand-written property tests (`pauli_sum_hash.rs:46/64`, `tableau_hash.rs:71/91`,
`phase1_leaf_types.rs:1310`) whose docstrings cite the design doc, never a
theorem. None of the nine `*_lean.rs` oracles asserts anything about a digest,
so the sector has a bridge-only leg with nothing on the other end. Clause 3 is
the weakest: `hash_into` is exercised only by two crate unit tests over 4
hard-coded words at one storage width (`pauli-word-2/column.rs:341`,
`lossy/column.rs:372`), and it has no production call site at all (ColumnStore
fills its digest column from `key.key_hash()`), so the agreement it is credited
with is never load-bearing where it is tested.

**Proposed closure.** Lean: a small `PPVM/Hash.lean` defining a key as
(structural identity, digest = f identity) and proving
`keyHash_congr : a ≈ b → digest a = digest b`, `hash_eq_writeU64`, and
`hashInto_eq_keyHash` for the column model, plus the pass-through composition
`identityHasher ∘ hash = digest`. Oracle: rename/extend the contract tests into
`hash_lean.rs` citing those theorem names, add `hash_into` over seeded random
columns at every storage tier, and make ColumnStore's batch path actually call
`hash_into` so mutating it fails a test.

### G-011 — Lossy canonical-loss invariant is pinned only for parsed words, never for branch builders

- Class `coverage` · Tier B · Severity medium · Sector hashing-digests · Status **adjudicated-spec** (round 2 — [adj U3](#u3-lossy-canonicality))
- Rust: `crates/ppvm-lossy-pauli-word-2/src/data.rs:433`
- Lean: none

**Claim.** Two `LossyPauliWord`s with the same logical content (same loss mask,
same Pauli at every present site) are `Eq` and have the same `key_hash()` —
which requires the canonical-loss invariant `l[i] ⇒ x[i] = z[i] = 0`, since
`PartialEq` (`data.rs:541`) and `structural_hash_lossy` read the full X/Z planes
including lost slots.

**Why unverified.** `set_lost` (`data.rs:155`) and the
`set_x_bit`/`set_xz_bits`/`set_z_bit_pair` family (`data.rs:352-416`) all
maintain the invariant explicitly (clearing X/Z on loss, clearing loss on a
nonzero X/Z write). `PauliBits::toggled_bits`/`toggled_bits2` — the rotation
branch builders — do not: they copy the planes, toggle X/Z at site `i`, and pass
`self.lbits` through unchanged, so at a lost site they produce `(x=1, l=1)`.
`Word::get` reports `LossySite::Lost` for that word and `Display` prints 'L', so
it is logically identical to the canonical lost word, yet `Eq` says unequal and
the digest differs — the same operator would occupy two `PauliSum` entries.
Nothing local prevents it; safety comes only from `is_lost` guards at every call
site (`rotation.rs:165/240/262`, `loss.rs:85`, `proj.rs:76`). The oracle
`lossy_pauli_word_lean.rs:356` checks the invariant only on words built by
`From<&str>` and after Clifford circuits
(`clifford_leaves_loss_mask_invariant`, which asserts the loss mask, not the
X/Z-at-lost-site canonicality of branch output), and never states the digest
consequence. `word-data-structures.md:452` asks for "tests enforcing lost => X/Z
identity after every mutator"; the branch builders are the mutators it misses.

**Round-2 amendment.** The uniqueness claim is confirmed — blob equality/hashing
equals logical equality/hashing exactly on the canonical set
C = {(x,z,l) : l ⇒ x = z = 0} (5 of 8 triples) — and the `toggled_bits` hole is real
and measured: the non-canonical "LZ" has `key_hash` 0x8f3d55be65ef69de vs
0xf2703b327e0c3932, `from_terms([LZ, LZ])` keeps two entries for one operator, and
`loss_channel(0, 0.25)` on it emits ("XZ", 0.25) where the canonical word gives
("IZ", 0.25). It is unreachable today (every in-tree caller guards on `is_lost`), so
this is an API-contract defect, not a live one. The closure must also cover
`LossyPauliKeyColumn::toggled_bits`/`toggled_bits2` (`column.rs`) and the `KeyColumn`
defaults (`ppvm-traits-2/src/batch.rs:135`/`:145`). See
[adj U3](#u3-lossy-canonicality).

**Proposed closure.** Lean: define the lossy site model with `CanonicalLoss` and
prove `logicalEq ↔ blobEq` and `logicalEq → digest =` under it, plus closure of
`CanonicalLoss` under each mutator including `toggledBits`. Oracle: extend
`lossy_pauli_word_lean.rs` to apply `assert_site_invariants` to the output of
`toggled_bits`/`toggled_bits2`/`loss_cleared`/`set_*` at *lost* sites over
seeded random words, and assert `a.key_hash() == b.key_hash()` for any two words
with equal logical readout — which would fail today and force either a loss
guard in `toggled_bits` or a documented precondition.

### G-012 — Mixture fingerprint XOR deltas justified only by agreement with legacy `ppvm-tableau-sum`

- Class `provenance` · Tier B · Severity medium · Sector hashing-digests · Status **open**
- Rust: `crates/ppvm-tableau-2/src/mixture/noise/pauli.rs:13`
- Lean: none

**Claim.** The incrementally computed branch fingerprint equals the from-scratch
one: `fingerprint(apply_mutation(tab, m)) == fingerprint(tab) ^ delta(m)` for
every `Mutation` (Pauli/Pauli2/Loss/Loss2) — the homomorphism that makes
structurally identical branches land in the same bucket and merge.

**Why unverified.** `insert_lazy_branches`/`insert_branches`
(`mixture/data.rs:174-258`) merge a branch only if some existing entry shares
its supplied fingerprint, while `rebuild_buckets` (`data.rs:156-172`) recomputes
all fingerprints from scratch with `fingerprint()` whenever `dirty` — so the
incremental delta and the full recompute must agree exactly or merging becomes
dirty-flag dependent. The deltas are derived by hand (`pauli_deltas` XORs
`sign_mask(row)` over rows whose site is Y|Z for X, etc.; `loss_mask(qubit)` for
loss) against `phase_loss_hash`'s XOR-of-masks encoding, and the whole scheme
including `sign_mask`/`loss_mask` is a port of legacy
`ppvm-tableau-sum/src/storage/mod.rs:107` and `noise.rs:78/498`. Its only
cross-check is `tableau_mixture_diff.rs`, a differential test against that same
legacy scheme, which by construction cannot detect a shared error; the crate's
own `fingerprint_collision_still_checks_full_structure` (`mixture/tests.rs:20`)
only tests the false-positive direction. A delta that is wrong for one mutation
silently stops merging, and the duplicated entries are then each thresholded
against `sum_cutoff` and `truncate`d independently, perturbing the sampled
distribution rather than just costing time.

**Proposed closure.** Lean: model the fingerprint as a GF(2)^64-linear map on
(phase bits, loss bits) and prove
`fingerprint (pauliAct p q s) = fingerprint s + Σ_{rows anticommuting} signMask row`
and the loss analogue, i.e. that the hand-derived delta table is the image of
the mutation. Oracle: an independent test (not a legacy diff) that for random
mixtures and every `Mutation` variant asserts
`fingerprint(&{clone; apply_mutation}) == fingerprints[parent] ^ delta`, plus a
test that a mixture built with `dirty=true` (forcing full recompute) has the same
entry count as one built incrementally.

### G-013 — Canonical loss invariant never checked after propagation; all oracle checks are Display-blind

- Class `bridge` · Tier A · Severity high · Sector lossy-word · Status **adjudicated-spec** (round 2 — [adj U3](#u3-lossy-canonicality))
- Rust: `crates/ppvm-lossy-pauli-word-2/src/clifford.rs:117`
- Lean: `lean/PPVM/Pauli/Symplectic.lean` `cnotActL_preserves_loss` / `cnotActL_lost_target_stays_identity`

**Claim.** After any sequence of Clifford generators on a `LossyPauliWord`,
every lost site still has x = z = 0 (LossInv is preserved), so the physical
encoding of a logical lossy word stays canonical.

**Why unverified.** Lean proves this for the model. On the Rust side nothing
pins it. `lossy_pauli_word_lean.rs:374-398` (`assert_site_invariants`) is the
only place that reads `x_bit`/`z_bit` under a lost site, and it is called only
on freshly *parsed* words (lines 360, 368). The post-circuit test
`clifford_leaves_loss_mask_invariant` (line 416) compares only the loss plane.
Every Clifford comparison in `lossy_pauli_word_diff.rs` (lines 248, 273, 279,
319, 330) uses `to_string()`, and `Display` (`data.rs:519`) prints 'L' whenever
`lbits[q]` is set no matter what the X/Z bits hold — so those assertions are
structurally incapable of observing the invariant they are credited with. The
only real check is one hardcoded 2-qubit case, `clifford.rs:227`
`cnot_present_control_lost_target_preserves_invariant` (`"XL"`), which exercises
the `xor_x_col` guard only.

Concrete mutation that survives the whole sector: delete
`if self.is_lost(ctrl) || self.is_lost(tgt) { return }` from `xor_z_col`
(`clifford.rs:117`). On input `"LZ"`, `cnot(0,1)` then sets z_0 = 1 on a lost
site, breaking `lost ⇒ x=z=0`; the word still displays `"LZ"` but its
`key_hash`/`PartialEq` blob now differs from the canonical one. `clifford.rs`'s
own table (lines 186-199) contains LI/LX/LL/IL/XL but no L-control-with-Z-target
row, so it passes; `clifford_leaves_loss_mask_invariant` passes (lbits
untouched); `lossy_clifford_generators_preserve_symplectic_form` passes (the lost
site's x bits stay 0, so its ω contribution stays 0); both diff tests pass
(Display-blind).

**Round-2 amendment — this row's headline mutation claim is FALSE.** Deleting the
`xor_z_col` loss guard does **not** survive the sector: `xor_z_col` is
`z_ctrl ^= z_tgt`, so a later CNOT whose *target* is the corrupted lost qubit XORs
the stray bit onto a **present** control's z bit, which `Display` shows, and
`lossy_pauli_word_diff.rs:248 clifford_replay_matches_old` fails at seed 1, n = 16
(`XLXXYZYXXLLXILXX` vs `XLXXYIYXXLLXILXX`) — reproduced independently by both
agents. What survives is only `lossy_pauli_word_lean.rs` (its loss-mask check is
phase-blind and its ω check uses a shared loss mask, under which the corrupted bit
multiplies zero — measured 1 vs base 1). So the severity is a **low coverage nit**
(no oracle asserts canonicality on the bits after propagation), not a high-severity
bridge gap, and the shipped guards are mathematically required and present: writing
a Z bit under a loss bit produces a triple that denotes no lossy word at all. The
proposed closure remains the right one. See [adj U3](#u3-lossy-canonicality).

**Proposed closure.** Add to `lossy_pauli_word_lean.rs` a test citing
`cnotActL_preserves_loss` that, for every one of the 25 two-qubit lossy words
and for random n-qubit words under a 200-gate random circuit, asserts
`!w.x_bit(i) && !w.z_bit(i)` for every `i` with `is_lost(i)` — plus a bit-level
(not `to_string()`) comparison of the whole (x, z, loss) triple against a
transcription of `hActL/sActL/xorXColL;xorZColL/czActL`.

### G-014 — Loss-plane mutators and canonical-encoding uniqueness have no Lean model; `toggled_bits` can break the invariant

- Class `coverage` · Tier B · Severity medium · Sector lossy-word · Status **adjudicated-spec** (round 2 — [adj U3](#u3-lossy-canonicality))
- Rust: `crates/ppvm-lossy-pauli-word-2/src/data.rs:421`
- Lean: none

**Claim.** The logical lossy word `Fin n → LossySite Pauli` is in bijection with
canonical (x, z, loss) triples, every mutator maps canonical to canonical, and
therefore the raw-blob `PartialEq`/`key_hash` on the three planes is equivalent
to comparing the logical sites (`data.rs:36-46`).

**Why unverified.** Lean models only the *Clifford* action on loss, and only
with loss as an external parameter. Nothing models the operations that actually
write the loss plane — `set_lost` (`data.rs:155`), `clear_loss` (175),
`with_lost` (189), `set` (205), `loss_cleared` (492), and the loss-*clearing*
`set_x_bit`/`set_z_bit`/`set_xz_bits`/`set_xz_bits2`/`set_z_bit_pair` (353 ff.)
— nor the encoding-uniqueness claim that licenses hashing/equality over raw
blobs. Two concrete cracks. (a) `LossyPauliWord::toggled_bits` (`data.rs:421`)
copies `self.lbits` unchanged while toggling X/Z, so calling it on a lost site
produces `lost = 1, x = 1` — a non-canonical word that `get`/`Display` report as
`L` but that hashes and compares unequal to the canonical `L`; it is safe today
only because every caller happens to guard on `is_lost` first (`rotation.rs:165,
240, 262, 336, 392, 424`; `noise.rs:316`), a precondition documented nowhere on
`PauliBits::toggled_bits`. (b) The trait's *default* `toggled_bits`
(`ppvm-traits-2/src/word.rs:235`) routes through `set_x_bit(i, true)`, which
*clears* loss — so the default and the lossy override resolve the same
loss/Pauli conflict in opposite directions, and neither is stated as the spec.

**Round-2 amendment.** Both cracks are confirmed, and the adjudicated spec is that
**neither** resolution is right: the override leaves the canonical set entirely
(its output denotes no lossy word) and the trait default resurrects a lost atom as
X, which the "a lost atom does not participate" model forbids. The honest
specification is the precondition `¬is_lost(i)`, as this closure already proposes;
if a total function were ever wanted, the only choice consistent with the gate model
is loss-wins identity, which is neither shipped variant. Widen the fix to
`LossyPauliKeyColumn::toggled_bits`/`toggled_bits2` and to the `KeyColumn` defaults
at `ppvm-traits-2/src/batch.rs:135`/`:145`. Note the bench at
`ppvm-conformance-2/benches/word_surface/lossy_branch.rs:18`/`:40` calls these on a
lossy word at sites 127/191 and is safe only because 127 % 5 = 2 ('Y') and
191 % 5 = 1 ('X'): pin the sites before adding the `debug_assert`. See
[adj U3](#u3-lossy-canonicality).

**Proposed closure.** Define the lossy word in Lean as
`Fin n → LossySite Pauli` together with its canonical encoding
`enc : LossySite → (Bool × Bool × Bool)`, prove `enc` injective on the canonical
image and that each mutator (`setLost`, `clearLoss`, `setSite`, `lossCleared`)
maps canonical to canonical (hence `w₁ = w₂ ↔ blob w₁ = blob w₂`), and state
`toggledBits` with the explicit hypothesis `¬ lost i` (with `toggledBits` off a
lost site left undefined / proved outside the canonical set). Pin it with an
oracle test that builds each logical word by several distinct mutator paths and
asserts equal blobs and equal `key_hash`, plus a
`debug_assert!(!self.is_lost(i))` in `LossyPauliWord::toggled_bits`/
`toggled_bits2` and a documented precondition on the trait method.

### G-015 — `weight()` counting Lost as non-identity is justified only by agreement with legacy

- Class `coverage` · Tier B · Severity low · Sector lossy-word · Status **open**
- Rust: `crates/ppvm-lossy-pauli-word-2/src/data.rs:288`
- Lean: none

**Claim.** Lean defines no lossy weight function, so the `MaxPauliWeight`
grading actually used to threshold a `LossyPauliSum`
(weight = |{q : x_q | z_q | lost_q}|, lossy `data.rs:288`) is never instantiated
into `GradedMap.lean`'s abstract `w : K -> Nat` truncation theorems; the
L-as-non-identity convention is pinned only by in-crate literals
(`data.rs:627`) and the legacy diff.

**Why unverified.** `weight()` is a fused popcount of `x | z | l`. The only
assertion about it anywhere outside the crate is
`lossy_pauli_word_diff.rs:121` (`new.weight() == old.weight()`), i.e. agreement
with the legacy `ppvm-pauli-word` lossy `weight` the docstring says it was ported
"verbatim" from. `lossy_pauli_word_lean.rs` never asserts `weight()` at all
(only `loss_weight()`, line 401). No Lean file defines a lossy weight:
`GradedMap.lean`'s truncation theorems take an abstract `w : K → ℕ`
(`retain_weight_le_eq_self`) and say nothing about which sites count. This is
decision-relevant: it decides which terms a MaxPauliWeight-truncated lossy sum
keeps. Ironically the independent grounding already exists next door and is
unused — `Noise.lean:300-305` establishes that the per-site identity of the full
space is `𝟙 = I + L`, so `L` is *not* identity and counting it is right; nobody
connects the two.

**Proposed closure.** State in Lean, over a lossy site type,
`weight w = card {q | siteOf w q ≠ Present I}` and prove
`weight = presentWeight + lossWeight` plus the justification that `L ≠ 𝟙`-component
(reusing `Noise.lean`'s `unit1 = I + L`). Then add an oracle test in
`lossy_pauli_word_lean.rs` asserting
`w.weight() == (0..n).filter(|i| w.get(i) != LossySite::Present(Pauli::I)).count()`
for exhaustive n ≤ 2 and random n, citing that theorem, so the fact no longer
depends on the legacy diff.

### G-016 — Lean models loss as an external predicate, not as word state, so loss-plane invariance is unstatable

- Class `fidelity` · Tier A · Severity high · Sector lossy-word · Status **adjudicated-spec** (round 2 — [adj U3](#u3-lossy-canonicality))
- Rust: `crates/ppvm-lossy-pauli-word-2/src/clifford.rs:102`
- Lean: `lean/PPVM/Pauli/Symplectic.lean:314` `variable (lost : Fin n → Prop)` / `LossInv`
- Verification note: this is the one surviving row that returned **no skeptic verdict**; a close round should re-read the premise before spending budget on the (large) proposed Lean extension.

**Claim.** The loss plane is part of the word's state; the guard reads it from
the word, the Clifford generators never write it, and therefore (a) different
keys in the same `LossyPauliSum` get different guard decisions from one gate,
and (b) the gate map is injective on keys because it preserves each key's loss
mask.

**Why unverified.** In Lean the state is `Sp n = Fin n → ZMod 2 × ZMod 2` and
`lost` is a *separate, fixed* parameter, so the model literally cannot express
"the action does not change the loss plane" — that plane is not part of the
modeled state. In Rust the loss bits live in the word (`data.rs:78 lbits`) and
`is_lost` reads them (`data.rs:145`, `clifford.rs:102/117/127`). The
consequences that matter are precisely the ones outside the model: a gate
applied to a lossy sum is a *piecewise* map (a different Sp element per loss
sector), and its injectivity — needed for "no key collisions" when propagating a
`LossyPauliSum`, i.e. for correctness of the hashbrown re-keying — rests on
loss-mask invariance. The oracle file itself concedes the gap:
`lossy_pauli_word_lean.rs:19-25` lists loss invariance as "genuinely new facts
with no algebraic Lean counterpart, asserted as property tests", and
`clifford_leaves_loss_mask_invariant` (line 416) cites no Lean theorem.

**Round-2 amendment.** The Rust half of the claim is **established** by derivation:
no `SymplecticColumns` primitive writes `lbits`, so the loss mask is invariant and
the gate map is sector-wise (an Sp element on the present sub-block, the identity on
lost sites, with distinct sectors having disjoint images because the mask is part of
the key) — hence injective on canonical lossy keys, which is what `RekeyBijective`
needs. Measured over 275 (word, generator) pairs and 48 random 200-gate circuits at
the bit level. The gap is exactly the Lean fidelity one this row names: with `lost`
an external `variable`, loss-plane invariance and key-injectivity are unstatable.
Worth recording for the closure: with *independent* masks ω is genuinely not
preserved (v = "XI", w = "LZ", `cnot(0,1)`), so the existing oracle's shared-mask
scoping is correct and must be kept. See [adj U3](#u3-lossy-canonicality).

**Proposed closure.** Extend the Lean model to a lossy state
`LSp n := Fin n → (ZMod 2 × ZMod 2) × Bool` (or `Fin n → LossySite`), define the
guarded generators over it with the guard reading the state, and prove:
`lossMask (gAct v) = lossMask v` for each generator; `LossInv` preservation as a
corollary; and `Function.Injective (gAct)` on the lossy state (sector-wise
bijectivity). Then have `clifford_leaves_loss_mask_invariant` cite the new
`*_preserves_lossMask` theorem, and add an injectivity oracle: propagating a set
of distinct lossy keys through a random circuit yields a set of the same
cardinality.

### G-017 — No oracle test pins any measurement arithmetic to Projection.lean

- Class `bridge` · Tier A · Severity high · Sector measurement-branching · Status **open**
- Rust: `crates/ppvm-tableau-2/src/measure.rs:407`
- Lean: `lean/PPVM/Tableau/Projection.lean:116` `rustTerm_eq` (also `:225` `probOne_eq`, `:260` `projectRaw_eq_two_proj`, `:348` `proj_zero_eq_caseB_retain`)

**Claim.** The real Rust `z_overlap_re` / `prob_1` / case-a merge / case-b retain
compute what `Projection.lean` says they compute.

**Why unverified.** `Projection.lean` is the most load-bearing file in the sector
and it is entirely unbridged. Grepping all nine `*_lean.rs` oracle files for
`overlap`, `prob_1`, `probOne`, `z_overlap`, `projectRaw`, `caseB_retain`,
`proj_zero`, `Born` yields zero hits; `tableau_lean.rs` contains only frame-level
tests (symplecticity, dichotomy, linear independence, the mod-2 frame identity,
the XOR relabel, rotation additivity). Nothing evaluates `overlapRe`
independently and compares it with the value the crate accumulates.
Consequently a mutation of the four-way ℤ/4 dispatch — e.g. swapping the odd arms
to `1 => z_overlap_re -= im_w, 3 => += im_w` at `measure.rs:410`, which is
precisely the `Re(i^φ conj a · b)` vs `Re(conj(i^φ a) · b)` slip that
`rustTerm_eq`'s docstring warns about — passes every `*_lean.rs` test in the
repo. The odd arms are genuinely live: `FrameInvolution` forces
`parity4 pd = dot G L`, which is 1 for many case-a decompositions, so `phase` is
odd there. The same holds for `compute_overlap_case_a` (`measure.rs:671`), which
`expectation.rs` also uses.

**Proposed closure.** Add to `tableau_lean.rs`: (1) a transcription of
`overlapRe`/`rustTerm` over the live coefficient list plus
`compute_phase_with_mask_static`, asserted bit-for-bit against a
`z_overlap_re` the crate exposes (or against
`compute_overlap_case_a`/`_case_b`, which are already `pub`), on the Clifford+T
sweep; (2) `projectRaw_eq_two_proj` as an executable check: after a case-a
`project_case_a(outcome, ...)` on a `coefficient_threshold = 0` tableau, every
stored coefficient equals `2·P_b` of the pre-measurement vector up to the single
`normalize()` scale factor; (3) `proj_zero_eq_caseB_retain` as an executable
check of the survivor set of `project_case_b`. Each must fail under a sign
mutation of the corresponding arm.

### G-018 — Mixture branch weights and the sub-cutoff branch drop have no Lean model

- Class `coverage` · Tier A · Severity high · Sector measurement-branching · Status **adjudicated-defect** (round 2 — [adj U4](#u4-mixture-weights-and-sampler))
- Rust: `crates/ppvm-tableau-2/src/mixture/measure.rs:86`
- Lean: none

**Claim.** `for_each_z_branch` unravels a measurement into a correct classical
mixture: each child weight is parent × Born probability, the two children's
amplitude vectors are the normalized projections, and the retained mixture is a
probability distribution.

**Why unverified.** There is no Lean file for `GeneralizedTableauMixture` at all
(grep for mixture/sampler/probability across `lean/PPVM` finds only
`Frame.lean`'s prose remark that `TableauMixture = C[Tableau]` is the same
`C[K]`; `Algebra/Noise.lean`'s probability vectors are Pauli-channel weights, not
mixture weights), and the only tests are `tableau_mixture_diff.rs` — legacy
agreement, explicitly not verification. The unproven arithmetic is load-bearing
and demonstrably lossy: when `p_other <= sum_cutoff` (`measure.rs:86`, and `:50`
for case b) the second branch is never pushed into `branches`, so
`insert_branches` never sees it, never sets its `dropped` flag, and
`normalize_probabilities()` is not called; `truncate()` renormalizes only if it
*removed* an entry. Failure scenario: a 1-entry mixture with
`sum_cutoff = 0.05` measuring a qubit with `p_one = 0.97` leaves a single entry
of weight 0.97 and total mass 0.97, and `measure()` returns
`(0.03→dropped, 0.97, 0.0)` summing to <1; repeat over k measurements and the
mass decays multiplicatively with no bound and no renormalization. Whether that
is intended (a sub-stochastic "mass we chose not to track") or a bug is exactly
what a spec would settle, and there is none.

**Round-2 amendment.** The "whether that is intended … or a bug" question is
settled: it is a **bug**. Π₀ + Π₁ = I forces the two children to sum to the parent,
and `measure()`'s reported marginals P(b) = Σᵢ pᵢ qᵢ(b) are determined *before* any
retention decision — so returning p(1) = 0.0 when p(1) = 0.03 is a false number,
not a truncation. There is also no convention to appeal to: legacy
`ppvm-tableau-sum/src/data.rs:136` asserts
`p_cum.last() >= 1 - sum_cutoff` ("Normalization error in sum") and **fires** on
this path, and the `-2` rewrite dropped that guard. See
[adj U4](#u4-mixture-weights-and-sampler).

**Proposed closure — amended; one clause INVALIDATED by round 2.** The oracle's
alternative "(or equals the Lean-predicted sub-stochastic mass)" must be dropped:
asserting the leaked mass as the expected value would prove the bug. The oracle
must assert Σ weights == 1 after every `for_each_z_branch`, that each child weight
equals parent × the independently recomputed `prob_1`, **and** that the returned
`(zero, one, lost)` triple sums to 1 and reports the true Born marginals. The Lean
half stands as written (`measureStep` maps a distribution to two children summing
to the parent; the truncated/normalized result is again a distribution; the error
from dropping mass m is bounded — the tight elementary bound is trace distance ≤ m,
|Δ⟨O⟩| ≤ 2m, not the looser 2m/(1−m) quoted in G-071).

### G-019 — MixtureSampler's shot distribution is unproven and the clamp misassigns the missing mass

- Class `coverage` · Tier A · Severity high · Sector measurement-branching · Status **adjudicated-defect** (round 2 — [adj U4](#u4-mixture-weights-and-sampler))
- Rust: `crates/ppvm-tableau-2/src/mixture/sampler.rs:43`
- Lean: none

**Claim.** `MixtureSampler::sample` draws shots from the correct Born
distribution — entry i with probability p_i, then the per-entry measurement
distribution.

**Why unverified.** No Lean models the categorical inverse-CDF selection, and no
`*_lean.rs` test exercises the sampler (only `tableau_mixture_diff.rs:205`
compares shot streams against legacy, and `mixture/tests.rs` compares serial vs
parallel). The selection is `partition_point(|&bound| bound <= p)` over the
cumulative sums, clamped by `.min(self.entries.len() - 1)`. `sampler()` does
*not* renormalize, so combined with the mass leak in G-018 the cumulative
vector's last element is < 1. Failure scenario: entries sorted descending with
weights [0.6, 0.3] (total 0.9 after a sub-cutoff branch was dropped) give
`cumulative = [0.6, 0.9]`; every draw in [0.9, 1.0) — 10% of shots — falls past
the end and is clamped to index 1, so the *smallest*-weight entry is sampled
with probability 0.4 instead of 0.3, a 33% relative error on that branch's shot
counts. The `.min(...)` also silently makes the sampler total on an empty mixture
index 0 of an empty vector's saturating length. Nothing states, proves, or tests
the intended distribution.

**Round-2 amendment.** Confirmed as a **live defect**, and the failure scenario in
this row is measured rather than hypothetical. Through the public API
(`h(0); measure(0); RY(1); measure(1)`) two entries with **identical** stored weight
0.485 are sampled at 0.485010 / 0.514990 over 400k shots (+6.18% relative on the
last entry, breaking a physical exchange symmetry); with weights [0.5, 0.3, 0.1] the
0.1 entry is sampled at 0.200745 (+100.7%) while the other two are untouched, which
is the clamp's signature. The empty-mixture panic is real too
(`new(n, thr, 1.0)` → "index out of bounds: the len is 0 but the index is 0" at
`sampler.rs:52`). The closure is right as written; the cheapest first step is to
restore legacy's dropped guard `debug_assert!(p_cum.last() >= 1 - sum_cutoff)`
(`ppvm-tableau-sum/src/data.rs:136`), which fires on this path. See
[adj U4](#u4-mixture-weights-and-sampler).

**Proposed closure.** Lean: prove the inverse-CDF lemma — for
`cum i = Σ_{j≤i} p j` with `Σ p = 1`,
`P(partition_point(cum, U) = i) = p i` for `U ~ Uniform[0,1)` — and that it fails
(mass concentrates on the last index) when `Σ p < 1`, making the normalization
precondition explicit. Oracle: a `tableau_mixture_lean.rs` test asserting
`cumulative.last() == 1.0` within tolerance whenever a sampler is built, plus a
many-shot empirical check that per-entry shot frequencies match the entry weights
within a Chernoff bound on a mixture with known weights.

### G-020 — mixture/equality.rs state equality has no Lean model and is not an equivalence relation

- Class `coverage` · Tier B · Severity medium · Sector measurement-branching · Status **adjudicated-defect** (round 2 — [adj U4](#u4-mixture-weights-and-sampler))
- Rust: `crates/ppvm-tableau-2/src/mixture/equality.rs:19`
- Lean: none

**Claim.** `structurally_equal` is the right notion of "the same state up to
representation", so merging two mixture entries under it is exact
(`mixture/data.rs:19-22`: "Entries merge only after a full frame/loss comparison
and coefficient-wise approximate comparison").

**Why unverified.** No Lean file mentions this predicate, and it fails the
properties a merge key needs. (1) Not transitive: `amplitudes_equal` accepts when
every per-index `|Δ|² < cutoff_sq`, so with cutoff c, states A and B differing by
0.9c and B and C differing by 0.9c are each "equal" while A and C (1.8c apart) are
not — so which entries merge, and hence the resulting weights, depends on
insertion order (`insert_branches` takes the *first* bucket match,
`mixture/data.rs:188/194`). (2) Merging non-identical states is silently lossy:
two entries up to `coefficient_threshold` apart per coefficient are collapsed and
their probabilities added (`data.rs:198`), an error nothing bounds — unlike
`Truncation.lean`, which does bound the analogous coefficient-dropping error.
(3) Asymmetric: the tolerance is read from `left` only (`equality.rs:16`) and only
`left`'s keys are iterated, so equality is not symmetric when the two tableaus
carry different thresholds. (4) It is *representation*-structural, not physical:
entries differing by a global phase, or by a different generating set of the same
stabilizer group, never merge — so the predicate is neither sound nor complete
for state identity, and no statement says which of the two it is meant to be.
Related: `structurally_equal_mutated` re-derives Pauli conjugation signs through
a third implementation, `pauli_flip` (`equality.rs:144-151`, single-site
anticommutation), with no citation of `Conjugation.lean`'s
`conjX`/`conjY`/`conjZ` and no oracle pinning it.

**Round-2 amendment.** Two of the four sub-claims are adjudicated differently.
(1) Non-transitivity is **confirmed** (a per-index sup-norm ball; a, a+0.6c, a+1.2c
share a `fingerprint` bucket and are really compared, and `insert_branches` takes
the first match, so merges and weights are insertion-order dependent) and the merge
error is bounded by 2·min(p_a,p_b)·c·√D — D-dependent, so the tolerance is not an
error budget in the readout's norm. (3) The asymmetry is **inert**, not a defect: a
left-only key would need |left_k| < c, contradicting the amplitude pruner's
`norm_sqr() > cutoff_sq` invariant (`data.rs:752`/`:1241`), and within one mixture
the threshold is shared. (4) is conservative — it costs entries, not correctness.
**The row misses the live one: reflexivity fails at `coefficient_threshold = 0`**
because `equality.rs:20` tests `<` against `cutoff_sq = 0`, so dedup is dead and the
entry count doubles per measurement round (measured 2 → 4 → 8 at threshold 0 vs
2 → 4 → 6 at 1e−12, mass 1.0 in both). Fix: `<` → `<=`. Add reflexivity-at-0 to the
oracle; non-transitivity should be documented, not "fixed". See
[adj U4](#u4-mixture-weights-and-sampler).

**Proposed closure.** Lean: define the intended equivalence (equal loss mask,
equal stabilizer *group* with signs, equal amplitude vector) and prove it is an
equivalence relation and that merging under it preserves the distribution
exactly; separately prove a total-variation bound of the form
`≤ (#merges) · c · Σ|c_k|` for the approximate variant so the tolerance is a
documented error budget rather than an accident. Oracle: a
`tableau_mixture_lean.rs` test asserting reflexivity/symmetry/transitivity of
`structurally_equal` over a pool of generated tableaus (transitivity will fail at
the current tolerance, which is the point), plus a test pinning `pauli_flip`
exhaustively to `Conjugation.lean`'s `conjX`/`conjY`/`conjZ` sign deltas.

### G-021 — Deterministic-branch sign unverified: Lean pins phase_decomp only mod 2

- Class `strength` · Tier A · Severity high · Sector measurement-branching · Status **adjudicated-spec** (round 2 — [adj U5](#u5-measurement-sign))
- Rust: `crates/ppvm-tableau-2/src/data.rs:496`
- Lean: `lean/PPVM/Tableau/BranchPhase.lean:224` `FrameInvolution` / `:238` `selfInverse_branchPhase_iff`

**Claim.** The deterministic (case-b) outcome the crate reports is the true
eigenvalue of Z_q on the state: `get_deterministic_outcome`'s
`result.phase >= 2`, and equivalently case-b's `z_sign = phase_decomp == 2`
(`measure.rs:336`), pick the correct one of the two outcomes rather than its
negation.

**Why unverified.** This is the sign-bug hot spot and nothing in the development
constrains it. `Frame.lean` works in `Sp n = GF(2)^{2n}`, phase-stripped on
purpose (`Frame.lean:238-240`: "The outcome bit `b` is a *phase*, and `Sp n` is
the phase-stripped space, so `b` does not appear below"), so no theorem there
mentions a sign. `BranchPhase.lean`'s only constraint on `phase_decomp` is
`FrameInvolution`, whose first term is `parity4 pd` — a ℤ/2 reduction
(`parity4_eq_zero_iff pd ↔ pd = 0 ∨ pd = 2`); `pd` is otherwise a free variable
in every theorem (`frameOp_eq_shiftOp`, `probOne_eq_crate`). The oracle matches
that weakness exactly: `tableau_lean.rs:925` computes
`parity = (phase + (destab&stab).count_ones() + (stab&mask).count_ones()) % 2`
and asserts it is 0, and the case-b assertion at `:937` is
`phase == 0 || phase == 2` — which admits *both* signs. So the ℤ/4 phase
accumulated by the tableau's own private g-rule (`data.rs:299-320`,
`(2*sign_count + imag_count) % 4`, a second implementation distinct from
`ppvm-pauli-word-2`'s, and never pinned to `Phase.lean`'s `phaseExp` /
`Word.lean`'s `phaseExpN` by any `*_lean.rs` oracle — see G-058) is verified only
up to a factor of −1. Only the `*_diff.rs` agreement with legacy stands behind
the sign.

**Round-2 amendment.** The sign is **correct**: outcome b = (accumulated ℤ/4 phase
== 2) is the true Z_q eigenvalue, derived from Z_q = ±∏_{i∈T} s_i with T the
anticommuting destabilizer set (M = i^p Z_q, M|ψ⟩ = |ψ⟩ ⇒ eigenvalue (−1)^{p/2}).
25,077 + 15,340 deterministic measurements across the two agents, 12,062 + 7,162 of
them with true eigenvalue −1, zero mismatches against a dense simulator. One
premise here is **overstated**: a `+2` phase mutation is not invisible.
`result.add_phase(2)` in `get_deterministic_outcome` is caught by
`ppvm-tableau-2/tests/behaviour.rs:627`/`:649` (absolute |0⟩ anchors — genuine
verification, not legacy agreement), and `p_word.add_phase(2)` in
`compute_decomposition` is caught by `tableau_lean.rs:700`
(`measurement_dichotomy_holds`). The blind spot is real only for the `*_lean.rs`
suite in isolation (17/17 green under the bare-path mutation), which is what this
row's own Rust anchor is about. See [adj U5](#u5-measurement-sign).

**Proposed closure.** Lean: extend `BranchPhase.lean` past the ℤ/2 shadow —
define `pd` as the `phaseExpN` fold over the coordinate generators selected by
`frame_coordinate_expansion` and prove (i) that fold equals the phase
`compute_decomposition` accumulates in ℤ/4, and (ii) for `L = 0`,
`pd = 2 ↔ the state is a −1 eigenvector of Z_q`, i.e.
`get_deterministic_outcome` returns the eigenvalue. Oracle: add a
`tableau_lean.rs` test that recomputes `compute_decomposition`'s phase
independently from the Lean `phaseExpN` product over the rows selected by the two
masks and asserts **full ℤ/4 equality** (today only parity is checked), plus a
case-b test anchoring the sign absolutely (e.g. `measure(q)` on `|0..0⟩` and on
`X_q|0..0⟩` against the Lean projector's predicate). A `p_word.add_phase(2)`
mutation must fail these.

### G-022 — Random-branch probability is never shown to be the Aaronson–Gottesman 1/2

- Class `strength` · Tier A · Severity high · Sector measurement-branching · Status **adjudicated-spec** (round 2 — [adj U5](#u5-measurement-sign))
- Rust: `crates/ppvm-tableau-2/src/measure.rs:75`
- Lean: `lean/PPVM/Tableau/Projection.lean:225` `probOne_eq` (and `Frame.lean:185` `measurement_dichotomy`)

**Claim.** When a stabilizer anticommutes with Z_q, the outcome is a fair coin —
`rng.random::<bool>()` on the bare frame is the Born distribution, and it agrees
with the coefficient-aware path's `prob_1 = 0.5 - 0.5·z_overlap_re`
(`measure.rs:420`).

**Why unverified.** `Frame.lean`'s `measurement_dichotomy` is explicitly only the
shape of the case split (its own docstring: "Structurally this is just the
dichotomy that an 𝔽₂-valued function is either identically 0 or hits 1
somewhere"), and `measure_deterministic_iff_xfree` is stated only for
`identityFrame`, not for a general symplectic frame. Neither carries any
probability. `Projection.lean`'s `probOne_eq` is the Born rule for the
*generalized* tableau's abstract amplitude vector; nothing instantiates it at the
bare-frame state (support-1 `c`, `s ≠ 0`) to derive 1/2, and nothing states the
standard one-line AG argument
⟨ψ|M|ψ⟩ = ⟨ψ|MS|ψ⟩ = −⟨ψ|SM|ψ⟩ = −⟨ψ|M|ψ⟩ = 0. So the crate has two different
case-a samplers — an unconditional fair coin for `Tableau` and an
overlap-derived probability for `GeneralizedTableau` — with no theorem and no
test asserting they coincide. A mutation biasing `measure.rs:75` (e.g.
`rng.random::<f64>() < 0.4`) violates no Lean statement and fails no `*_lean.rs`
test.

**Round-2 amendment.** The AG argument is confirmed and gives p1 = **exactly**
1/2: ⟨Z_q⟩ = ⟨ψ|S Z_q S|ψ⟩ = −⟨Z_q⟩ = 0 for any anticommuting stabilizer S, using
S† = S and S² = I (row realness, G-062). So the bare frame's fair coin *is* the Born
distribution and the coefficient path's `0.5 − 0.5·z_overlap_re` collapses to the
same number on a pure stabilizer state (support 1 ⇒ `z_overlap_re` is exactly 0.0);
measured worst |p1 − 0.5| = 2.220e−16 over 4,607 case-a branches and 3.886e−16 over
916 generalized Clifford-only branches, with Born frequencies within 2.69σ over 73
T-containing cases. The closure's proposed exact-0.5 assertion is therefore
legitimate as *exact* equality for Clifford-only states. See
[adj U5](#u5-measurement-sign).

**Proposed closure.** Lean: state `overlapRe s φ c = 0` whenever `s ≠ 0` and `c`
has support 1 (immediate: every term has `c k · c (k+s) = 0`), then `probOne_eq`
gives `prob_1 = 1/2` — the bare-frame fair coin, derived. Better still, add the
frame-level AG argument as a corollary of `frame_coordinate_expansion`. Oracle:
for a Clifford-only (T-free) `GeneralizedTableau`, assert the case-a `prob_1` is
*exactly* 0.5 for every measurement on the sweep (the pure-frame path and the
coefficient path then provably agree), and add a chi-square/interval check on the
bare `Tableau::measure` coin over many seeded shots.

### G-023 — Case-a post-measurement state rests on the assumed hypothesis hψ (and an unmodeled re-basing)

- Class `strength` · Tier A · Severity high · Sector measurement-branching · Status **adjudicated-spec** (round 2 — [adj U5](#u5-measurement-sign))
- Rust: `crates/ppvm-tableau-2/src/measure.rs:449`
- Lean: `lean/PPVM/Tableau/Projection.lean:260` `projectRaw_eq_two_proj`

**Claim.** The keep-A / transform-B / merge arm followed by
`update_tableau_according_to_outcome` leaves the state equal to the normalized
projection P_b|ψ⟩ — i.e. the merged coefficients are the correct amplitudes *in
the new frame's basis*.

**Why unverified.** `projectRaw_eq_two_proj` is stated under
`hψ : ∀ k, iPow (ψ k) = sgn b * iPow (φ k)`, i.e. it assumes the very relation
between the projection's phase and the overlap's phase that is the content of the
claim; the file's own scope note (`Projection.lean:52-63`) records this and says
relating them "needs a Hilbert-space model of the frame that this development
does not have" and that the justification is that "old and the -2 crate agree
here verbatim" — agreement with legacy, which the ground rules say is evidence of
nothing. Concretely the Rust projection phase is `alpha + 2·⟨idx, destab⟩`
(`measure.rs:449-450`) while the overlap phase additionally folds
`2·popcount(idx ∧ stab ∧ odd_phase_mask)` (`data.rs:144-148`); they coincide up
to `sgn b` only because the mask is always 0 in a valid frame — a fact the oracle
*asserts* (`tableau_lean.rs:954`) but no Lean theorem supplies. The deeper half
is untouched: the amplitude basis is `|j⟩ = ∏_l d_l^{j_l}|ψ₀⟩`
(`BranchPhase.lean:28-42`), and `update_tableau_according_to_outcome` (called at
`measure.rs:507`, *after* the coefficients are written) replaces
`destab[q_idx]` and multiplies `g_q` into other destabilizers — so the very basis
the surviving indices are read in changes, and nothing models that re-indexing.
`relabelAmp` (`Bitstring.lean`) covers only XOR relabels, not a frame change.

**Round-2 amendment.** The composite is **correct** — measured over 8,400
statevector comparisons (fidelity 1.000000000000) and, avoiding any amplitude-basis
convention, 939,648 full-4ⁿ-Pauli-set comparisons after every step of 700
Clifford+T circuits (0 failures, worst |Δ| 7.494e−16), including immediately after
case-a merges and the `update_tableau_according_to_outcome` re-basing. And the
hypothesis this row calls assumed is **derivable in three lines**: all 29 gate phase
deltas are `^= b << 1`, and frame-row `mul_assign` only ever sees commuting
Hermitian pairs, so row phases stay in {0,2} and `oddPhaseMask` is identically 0
(G-062) — so `hψ` is a consequence, not an axiom. The closure stands; the Lean work
should derive `hψ` rather than assume it. See [adj U5](#u5-measurement-sign).

**Proposed closure.** Lean: give `Frame` a phase-carrying Hilbert-space semantics
(a state map `Frame → Amp`, or at minimum the three bullet points
`BranchPhase.lean:28-42` lists as axioms turned into definitions), then (i)
derive `hψ` from `oddPhaseMask = 0` — itself derivable from row Hermiticity,
`phase ∈ {0,2}` (G-062) — plus `alpha = pd + 2b`, and (ii) prove that
`projectFrame`'s new destabilizer set re-indexes the merged amplitude vector by
the identity, so `projectRaw` in the old basis *is* the state in the new basis.
Oracle: for n ≤ 4 with `coefficient_threshold = 0`, expand both the pre- and
post-measurement `(frame, coefficients)` pair to a dense 2^n statevector and
assert the post state equals the normalized `(I ± Z_q)/2` projection of the pre
state — this pins the composite (merge + frame update) that no current test
touches.

### G-024 — rzz/rxx/ryy/rotate_2/comm_2/RotXY::r have no `*_lean.rs` oracle at all

- Class `bridge` · Tier A · Severity high · Sector multiply-rotation · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/rotation.rs:353`
- Lean: `lean/PPVM/Instantiations/Rotation.lean:553` `comm2_generic_sign_eq_branchExp2` (also `comm2Coeff_eq_zero_iff`, `comm2_key_eq_mulBits2`, `rzz/rxx/ryy_*`, `rotXY_heisenberg_order`)

**Claim.** The shipped two-qubit rotation surface — the generic `rotate_2` over
`comm_2`'s `SIGN_NEG=0x2840` mask with its `sin.mul_sign(-eps)` convention, the
three native `+eps` kernels, and `RotXY::r`'s Heisenberg sub-rotation order —
computes what `Rotation.lean`'s two-site and RotXY sections prove.

**Why unverified.** `Rotation.lean:471-489` argues at length that the
off-diagonal axis pairs are "untested on both sides" and that "a +eps/−eps
asymmetry is exactly the kind of thing that is correct by coincidence on the
diagonal", then declares "This section closes it". It closes the Lean half only.
`rg` over `crates/` finds comm2Coeff / comm2Key / signNegMask / anti2 /
`*_eps_from_product` / rotXY_* referenced *only* in src docstrings, never in a
test. Neither `pauli_sum_rotation_noise_lean.rs` nor any other `*_lean.rs` calls
rzz, rxx, ryy, rotate_2 or r on a `Sum` — grep for
`\.rzz(|\.rxx(|\.ryy(|rotate_2(` across `crates/ppvm-conformance-2/tests/*_lean.rs`
returns nothing. What exists instead is precisely what the Lean docstring
dismisses: `crates/ppvm-pauli-sum-2/tests/gate_surface.rs:62-92`
`assert_matches_generic` (fast-path vs generic at the three diagonal axes only,
ported from legacy's `rxx_matches_generic`), `gate_surface.rs:97`
`rzz_explicit_values` ("old's rzz_explicit_values"), and legacy diffs in
`pauli_sum_integration_diff.rs:366-407`. So the `−eps` at `rotation.rs:353` and
the `0x2840` mask at `rotation.rs:123` are pinned by nothing except agreement
with legacy.

**Proposed closure.** Add a two-site block to
`pauli_sum_rotation_noise_lean.rs` that is exhaustive over the finite grid the
Lean theorems quantify over: all 16 axis pairs ([x,z]×[x,z]) × all 16 two-site
keys, asserting (a) whether the term branches == `anti2`
(`comm2Coeff_eq_zero_iff`), (b) the branch key == `mulBits2`
(`comm2_key_eq_mulBits2`), (c) the branch coefficient ==
`c·sinθ·realSign(branchExp2)` (`comm2_generic_sign_eq_branchExp2`) — 256 cases,
cheap, and it covers the off-diagonal axes for the first time. Then assert each
native kernel equals the generic path on the same grid, and add `r(q,φ,θ)` vs the
Rodrigues `rotAxis` of `rotXY_heisenberg_order` at a few (φ,θ) including φ=0 and
φ=π/2.

### G-025 — No oracle pins any branch sign ε; a global flip of any axis column passes every test

- Class `bridge` · Tier A · Severity medium · Sector multiply-rotation · Status **adjudicated-spec** (round 2 — [adj U2](#u2-rotation-direction-and-sign))
- Rust: `crates/ppvm-pauli-sum-2/src/rotation.rs:266`
- Lean: `lean/PPVM/Instantiations/Rotation.lean:105` `rz_eps_from_product` (also `rx_`/`ry_`/`rzz_`/`rxx_`/`ryy_eps_from_product`)

**Claim.** No `*_lean.rs` oracle asserts any branch coefficient, so none of
`rx_`/`ry_`/`rz_`/`rzz_`/`rxx_`/`ryy_eps_from_product` is bridged to the Rust;
the absolute sign is held only by in-crate ported tests
(`gate_surface.rs:97 + :144 + :204`), never by an oracle citing the Lean ε
theorems.

**Why unverified.** Every assertion in `pauli_sum_rotation_noise_lean.rs` is
invariant under negating an axis's ε column, i.e. under replacing `r{axis}(θ)` by
`r{axis}(−θ)`: `rotation_preserves_l2_norm` (`:80`) is a quadratic form;
`rotation_is_reversible` (`:106`) composes the mutated kernel with its own
mutated inverse; `rotation_angle_addition_and_trotter` (`:141`) likewise uses the
mutated kernel on both sides; `rotation_by_zero_is_identity` (`:267`) has sinθ=0;
`anticommuting_branch_produces_distinct_new_key_with_xor_bits` (`:217`) inspects
only the key set and support size at θ=π/2 and never reads a coefficient;
`branch_keys_have_correct_arity` (`:372`) checks arity only. So mutating
`rotation.rs:266` from `if z {1} else {-1}` to `if z {-1} else {1}` — or the
equivalent in ry (`:248`), in all three copies of the rx sign in
`column_store/rotations/rx.rs` (`:65, :117-125, :163-171`), or in
`store.rs:1421`'s HashMapStore `rotate_x` — leaves the entire `*_lean.rs` suite
green. The tests do incidentally pin the *relative* sign inside a branching pair
(norm preservation fails if ε(k) and ε(br k) are not opposite), but the absolute
sign, which is the whole content of the five `*_eps_from_product` theorems, is
unbridged.

**Round-2 amendment.** The complaint stands after adjudication — the absolute ε is
still pinned by nothing, and both agents re-confirmed that a global flip of any axis
column passes every test in the repo — but the *values* to assert are now settled and
this row's proposed grid is correct as written (`rz` on X gives
{X: cosθ, Y: −sinθ}; measured `rz(0.7)` on X = [("X", 0.7648421872844885),
("Y", −0.644217687237691)]). The oracle must compare against **R†·mat·R**, not
R·mat·R†. See [adj U2](#u2-rotation-direction-and-sign).

**Proposed closure.** The same matrix oracle as G-030, or a cheaper standalone:
for each axis and each of the 4 single-qubit keys, assert the exact branch
coefficient at θ (e.g. rz on X gives support {X: cosθ, Y: −sinθ} — derived from
the Lean ε value, not from legacy), exhaustively over the 3 axes × 4 keys grid,
and cite `rx_`/`ry_`/`rz_eps_from_product` by name. Extend the same grid to two
sites for rzz/rxx/ryy.

### G-026 — The operator product on the columnar and indexmap backends has no Lean-oracle bridge

- Class `bridge` · Tier B · Severity medium · Sector multiply-rotation · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/column_store/graded.rs:158`
- Lean: `lean/PPVM/Algebra/Twisted.lean:213` `twistedConv`

**Claim.** Every shipped backend's `multiply_into` / `multiply_in_place` computes
`twistedConv`, and (for the columnar SoA) preserves the key/coeff/digest/bucket
column alignment.

**Why unverified.** `pauli_sum_multiply_lean.rs` pins only
`PauliSum = Sum<HashMapStore<PauliWord, C>, P>` (`multiply.rs:48`/`lib.rs:110`).
`ColumnStore::multiply_into` (`column_store/graded.rs:158-179`, which reserves
via `reserve_for_live_len` and appends through `primary.add`) and
`IndexMapStore::multiply_into` / `multiply_in_place`
(`indexmap_store/algebra.rs:111`, `branching.rs:64`) are separate
implementations pinned only by diff tests. Sharper still:
`column_store_lean.rs:33` lists "the operator product" by name among the
mutations that must preserve column alignment and presents `assert_aligned` as
the executable form of that invariant — yet the file never calls
`multiply`/`multiply_in_place` at all (no occurrence of `multiply` in the file),
so the one hazard it was written to cover is uncovered for that operation. A
missed `reindex()` after a product's growth would leave iteration correct while
every subsequent `get`/`overlap`/branch-merge silently missed — exactly the
failure mode the file's own header describes.

**Proposed closure.** In `column_store_lean.rs` add a
`Sum<ColumnStore<PauliWord, Complex<f64>>, _>` product test: `assert_aligned`
after `multiply_into`, after `multiply_in_place` run twice through the same store
(aux reuse), and after `mul_word_assign`; and assert term-by-term agreement with
the hash-join `PauliSum` product and with the dense `Mat` oracle already written
in `pauli_sum_multiply_lean.rs` (factor `Mat` into `ppvm-conformance-2` so both
files share it). Same for the indexmap backend.

### G-027 — mul_word_assign's plain-insert re-key is not pinned to the product by any oracle

- Class `bridge` · Tier A · Severity low · Sector multiply-rotation · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/multiply.rs:277`
- Lean: `lean/PPVM/Algebra/Twisted.lean:250` `twistedConv_single_right` / `:264` `twistedConv_single_right_apply` / `:222` `mulWord_isRightCancellative`

**Claim.** No `*_lean.rs` oracle pins `mul_word_assign` to
`twistedConv_single_right` or to the dense matrix model; the only pin is a single
fixed n=3 in-crate unit test (`multiply.rs:539`).

**Why unverified.** `multiply.rs:267-275` credits three Lean theorems to this
path and `multiply.rs:256` notes the merge is a plain `insert` guarded only by a
`debug_assert` — so a violated injectivity precondition DROPS a term in release.
The bridge for it is an in-crate unit test (`multiply.rs:539`
`mul_word_assign_matches_the_general_product`) at one fixed n=3 support and one
fixed word, plus the by-value spelling test at `:551`. In
`pauli_sum_multiply_lean.rs`, `mul_word_assign` appears only as one step of the
buffer-interleaving sequence (`:615`), where it is compared against a replay of
*itself* in a fresh store — that detects buffer leakage, not a wrong re-key. So
the oracle suite never asserts `mul_word_assign == multiply(single-term)` nor
`mul_word_assign` against the matrix model, and the exhaustive single-qubit grids
in that file all go through `multiply`.

**Proposed closure.** Extend `pauli_sum_multiply_lean.rs`: exhaustive over the 16
single-qubit supports × 4 words, and randomized for n≤4, assert
`mul_word_assign(q)` equals `multiply(single-term q)` term-by-term AND equals
`mat(A)·mat(q)` under the dense oracle, plus `len()` invariance (the injectivity
consequence `twistedConv_single_right_apply` licenses).

### G-028 — levi_civita (the generic single-site ε/key table) has no Lean counterpart

- Class `coverage` · Tier A · Severity medium · Sector multiply-rotation · Status **adjudicated-spec** (round 2 — [adj U2](#u2-rotation-direction-and-sign))
- Rust: `crates/ppvm-pauli-sum-2/src/rotation.rs:75`
- Lean: none

**Claim.** `levi_civita(i, j)` returns (ε, k) with −i·[P_i,P_j]/2 = ε·P_k over the
2-bit code 00=I, 01=X, 10=Z, 11=Y, including the branch-free `rank` bit trick and
the mod-3 `diff` computation, and its ε agrees with the rx/ry/rz columns when
called as `levi_civita(key, axis)`.

**Why unverified.** Nothing in `lean/` transcribes it — `rg -i levi lean/` matches
only a prose mention at `Rotation.lean:115` saying the *two-site* kernels consult
no Levi-Civita table. `Rotation.lean`'s single-site facts are all stated as
branchExp(G,P), i.e. axis-first, whereas `rotate_1_branch` calls
`levi_civita(p_g, pauli_code(axis))` at `rotation.rs:171` with the *key* first
(the local is misleadingly named `p_g` but is built from `k.x_bit`/`k.z_bit` at
`:170`). Since [P,G] = −[G,P], the orientation is the one thing that must not be
wrong, and it is unmodeled. Verified by hand that it currently agrees (rank X=0,
Y=1, Z=2; ε(Z,X)=+1 matches rx's +1 on key Z; ε(Y,X)=−1 matches rx's −1 on key Y;
and −i[P,G]/2 = iGP for anticommuting pairs), so this is a coverage hole rather
than a live bug — but it is also unbridged: `rotate_1_branch` is reachable only
through the lossy fallbacks at `rotation.rs:337/340/393/396/424/427`, and no test
in `crates/ppvm-conformance-2/tests` (lean or diff) drives a two-qubit rotation on
a `LossyPauliSum` with a lost site, so nothing exercises the table at all in the
`-2` crates.

**Round-2 amendment.** The hand-checked agreement is confirmed by machine: the
lossy `rotate_2` fallback matches the native single-site columns on all 4×4×4 cases,
and `levi_civita(key, axis)` with **+ε = the coefficient of iGP** is the correct
orientation (−i[P,G]/2 = +iGP for anticommuting pairs; ε(Z,X) = +1 = rx's sign on
key Z). The transcription must pin that argument order explicitly — key first, axis
second — since swapping it negates every ε. See
[adj U2](#u2-rotation-direction-and-sign).

**Proposed closure.** Transcribe it: `def leviCivita (i j : Fin 4) : ℤ × Fin 4`
verbatim including `rank`/`diff`, then
`theorem leviCivita_eq_branchExp : ∀ i j, anti i j → (leviCivita i j).1 = realSign (branchExp (axisBits j) (keyBits i)) ∧ (leviCivita i j).2 = mulBits …`
by `decide` over all 16 pairs — which simultaneously fixes the argument-order
convention. Oracle: a `LossyPauliSum` test that `rotate_2`/`rxx`/`ryy` with one
site lost equals the corresponding `rotate_1` on the survivor, exhaustive over
the 4 axes × 4 survivor keys.

### G-029 — The shipped rotation direction is justified only by agreement with legacy

- Class `provenance` · Tier A · Severity medium · Sector multiply-rotation · Status **adjudicated-spec** (round 2 — [adj U2](#u2-rotation-direction-and-sign))
- Rust: `crates/ppvm-pauli-sum-2/src/rotation.rs:228`
- Lean: `lean/PPVM/Instantiations/Rotation.lean:90` `rx_eps_from_product`

**Claim.** The `-2` crate rotates observables in the correct direction
(equivalently: `rx(θ)` is not `rx(−θ)`).

**Why unverified.** Combining G-030 and G-025: the Lean side declares the
conjugation modeled-not-derived and encodes the sign by fiat in `branchExp`, and
no `*_lean.rs` assertion is sensitive to the sign. What is left holding the
direction in place is legacy: `column_store_diff.rs:325-330` freezes
`GOLDEN = 2.1610566562692544`, explicitly sourced as
`ppvm-pauli-sum/tests/trotter.rs:104`; `gate_surface.rs:97 rzz_explicit_values`
is labelled "old's rzz_explicit_values"; `gate_surface.rs:62
assert_matches_generic` is old's `*_matches_generic`; and `rotation.rs:190-193`
says the commute test, flipped bits and ε are "ported bit-for-bit from
ppvm-pauli-sum::sum::rot1". Per the audit's own rule, a passing diff against
legacy is evidence of a faithful port and of nothing about correctness — so if
legacy's rotation direction were inverted, every gate in this sector would be
inverted and every test in the repo would still pass. Note the in-crate
`r_is_heisenberg_ordered` / `rotXY_zero_eq_rx` / `rotXY_halfPi_eq_ry` checks are
all *relative* (r vs ry), so they cannot supply the missing absolute anchor.

**Round-2 amendment.** Adjudicated: the shipped direction is **correct** — it is
O ↦ R†OR, coherent with the crate's own Clifford family (`Sum::s` ships
X ↦ −Y = S†XS, distance 0.000 to S†PS and 2.000 to S P S†, and `rz(π/2)` matches `s`
to 6.1e−17) and with `r(q,π/2,θ) == ry(q,θ)` (2.0e−17, vs 1.476e0 = 2·sin(0.83) for
`ry(−θ)`). So legacy's direction is right and the golden master is not hiding an
inversion. The row survives only as the missing *absolute* anchor: the ℂ-matrix
oracle proposed here and in G-030 was run by both agents and passes, but neither
landed it. See [adj U2](#u2-rotation-direction-and-sign).

**Proposed closure.** Anchor the direction to physics once, independently of
legacy: the ℂ-matrix conjugation test proposed in G-030 (`mat` after `rx(θ)` ==
`R†·mat·R` with R built from cos/sin directly), plus a cross-check that
`rz(π/2)` on a `Sum` agrees with the crate's own Clifford S/S† re-key up to the
documented convention. Once either exists, the golden master becomes a regression
bar rather than the specification.

### G-030 — Rotation conjugation is modeled not derived — and the asserted identity is sign-wrong

- Class `strength` · Tier A · Severity high · Sector multiply-rotation · Status **adjudicated-spec** (round 2 — [adj U2](#u2-rotation-direction-and-sign))
- Rust: `crates/ppvm-pauli-sum-2/src/rotation.rs:9`
- Lean: `lean/PPVM/Instantiations/Rotation.lean:17` + `:75` `def branchExp` / `:81` `branchExp_isRealPhase`

**Claim.** For anticommuting single-qubit G and P, conjugating a stored term
(P,c) by the rotation exp(-i·θ/2·G) yields c·cosθ·P + c·sinθ·(iGP), so the ε
columns `rx/ry/rz_eps_from_product` and `rzz/rxx/ryy_eps_from_product` are the
physically correct branch prefactors.

**Why unverified.** `Rotation.lean:21` and `:31` state outright that the operator
identity is *modeled* by the 2-D `rot`, "**not** derived from operator algebra
here"; `lean/README.md:73` repeats it. The consequence is already visible: the
identity as written at `Rotation.lean:17` (and verbatim at `rotation.rs:7-9`) is
FALSE. Direct computation with R = cos(θ/2)I − i·sin(θ/2)G, using only G²=P²=I
and GP=−PG, gives R P R† = (c²−s²)P + i·cs·(PG−GP) = cosθ·P − sinθ·(iGP) — the
opposite sign of the stated RHS. Checked concretely at G=X, P=Z:
e^{-iθX/2} Z e^{+iθX/2} = cosθ·Z − sinθ·Y, whereas iGP = iXZ = +Y and the Rust
(`rotate_x`, sign = +1 when the x-bit is clear) emits +sinθ on key Y. So the
shipped code computes U†PU (the backward/Heisenberg propagation the `RotXY`
docstring at `rotation.rs:452` says a `Sum` performs), while both the Lean spec
docstring and `rotation.rs`'s module header claim it for UPU†. Because
`branchExp := 1 + phaseExp(G,P)` bakes in the +sin convention by fiat and every
downstream theorem is `by decide` over that definition, `lake build` proves the
whole tower without ever adjudicating the sign, and the 3x3 mz/mx/my matrices
underpinning `rotXY_heisenberg_order` are explicitly "read off the kernel"
(`Rotation.lean:582`), so they inherit the same unadjudicated orientation rather
than testing it.

**Round-2 amendment.** The algebra in this row is **confirmed**
(R P R† = cosθ·P − sinθ·(iGP), R† P R = cosθ·P + sinθ·(iGP)) and the shipped kernels
compute R† P R — the backward/Heisenberg direction, which is correct for an
observable and is what `rotation.rs:451` documents. Every ε column, both two-site
paths, all three stores and the lossy fallback agree with R†·mat·R (worst 6.661e−16;
the forward model is off by exactly 2·sinθ). Two precision corrections: the false
identity is **not** "verbatim at `rotation.rs:7-9`" — `rotation.rs` never writes the
sandwich, it writes a direction-ambiguous "conjugates … to", so the honest tally is
one literally false Lean equation (`Rotation.lean:17`, with `:31` referring back) plus
three direction-ambiguous Rust docstrings, the third being
`ppvm-pauli-sum-2/src/store.rs:306-308`'s `RotateInPlace` doc, which the fix must
also correct. The proposed oracle is right and passes as written. See
[adj U2](#u2-rotation-direction-and-sign).

**Proposed closure.** State the conjugation in Lean rather than asserting it. No
analysis is needed: in a Pauli algebra with G²=P²=I and GP=−PG, prove
`(c•1 - s•(i*G)) * P * (c•1 + s•(i*G)) = (c^2-s^2)•P - (2*c*s)•(i*G*P)` and
specialise c=cos(θ/2), s=sin(θ/2) to get `R P R† = cosθ·P − sinθ·(iGP)` and
`R† P R = cosθ·P + sinθ·(iGP)`; then define `branchExp` as the *conclusion* of
that lemma for the direction the crate ships, and fix the two docstrings. Oracle:
add to `pauli_sum_rotation_noise_lean.rs` the dense-matrix check already
available in `pauli_sum_multiply_lean.rs` — build R = cos(θ/2)I − i·sin(θ/2)G as
a 2ⁿ×2ⁿ ℂ matrix and assert `mat(rx(θ) applied to A) == R† · mat(A) · R` for
random n≤4 supports and random θ. That single test pins the direction, the
per-axis ε columns and the branch keys simultaneously.

### G-031 — Columnar rx fuses the two passes; accumulate_rotBatch does not license the paired 2x2 update

- Class `strength` · Tier B · Severity medium · Sector multiply-rotation · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/column_store/rotations/rx.rs:127`
- Lean: `lean/PPVM/Instantiations/Rotation.lean:388` `accumulate_rotBatch` / `:433` `eagerWalk_ne_twoPass`

**Claim.** The columnar `RotateInPlace` keeps the two-pass ordering (scale ALL
diagonals, then merge ALL branches), so `accumulate_rotBatch` is its correctness
licence.

**Why unverified.** `store.rs:355` states "Every implementer owes the two-pass
ordering" and `column_store/mod.rs:62-66` asserts the columnar backend keeps it,
citing `eagerWalk_ne_twoPass`. The shipped columnar rx does not: when
`closed_support` holds, `rotate_x_kernel` takes a *fused* path that, for a row i
whose x-toggled partner j is present, writes
`coeffs[i] = ci*cos + cj*sin*sign_j; coeffs[j] = cj*cos + ci*sin*sign_i`
(`rx.rs:127-129`, and again in the ≤512-row copy at `rx.rs:173-174`) — a branch
contribution merged into row i while rows after i are still unscaled. That is
exactly the interleaved shape `eagerWalk_ne_twoPass` exhibits as *wrong*. It
happens to be correct here, but only because of two facts the Lean never states:
on the anti set (z-bit set at the qubit) the branch map `toggle x` is an
involution, so the collision graph is a perfect matching with no chains; and
ε(br k) = −ε(k), so the fused 2x2 matrix is [[cos,−sin],[sin,cos]]. Neither the
involution property nor the ε antisymmetry appears anywhere in `Rotation.lean`
(`anticommute_new_key` only rules out br k = k), so the theorem credited with
licensing this backend does not cover the code it is cited on, and the doc claim
at `column_store/mod.rs:62` is false as written. No oracle covers it either:
`column_store_lean.rs` runs `s.rx` exactly once (`:207`) as an alignment probe at
≤8 sites, and `column_store/tests/rotations.rs` only diffs the dense kernel
against the sparse *same* kernel.

**Proposed closure.** State the side condition in Lean: `pairFuse` (the paired
in-place update) equals `twoPass` whenever `br` is an involution on the anti set,
`anti (br k)` holds, and `s (br k) = -s k`; discharge the Pauli instance with
`br = toggle x`, `branchExp G (mulBits G P) = branchExp G P + 2` (a `decide`).
Then either weaken the `column_store/mod.rs` claim to "two-pass or a fusion
satisfying pairFuse's hypotheses", or bridge it: in `column_store_lean.rs` build
a >512-row `ColumnPauliSum` with a closed x-partner support (the shape
`column_store/tests/rotations.rs:24` already constructs) and assert term-by-term
equality against the same rotation on a HashMapStore-backed `PauliSum`, which is
the backend whose two-pass shape `accumulate_rotBatch` actually models.

### G-032 — rot_norm_sq / rot_neg_rot / rot_rot are unconnected to the kernel's per-key ε

- Class `strength` · Tier A · Severity medium · Sector multiply-rotation · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/rotation.rs:248`
- Lean: `lean/PPVM/Instantiations/Rotation.lean:275` `rot_norm_sq` (and `:243` `def rot`, `:264` `rot_neg_rot`, `:284` `rot_rot`)

**Claim.** "The branch is a norm-preserving, angle-additive 2-D rotation on the
coefficient pair" (`rotation.rs:34-37`, `store.rs:344-347`,
`traits-2/src/gates.rs`, `sym-2/src/coeff.rs` all cite `rot_norm_sq`/`rot_rot`
for exactly this).

**Why unverified.** `rot` (`Rotation.lean:243`) is the fixed matrix
[[cos,−sin],[sin,cos]]: its off-diagonal signs are hard-coded, +sin below and
−sin above. The Rust kernel's off-diagonals are ε(k)·sin and ε(br k)·sin, where ε
is the per-key table (`if z {-1} else {1}` at `rotation.rs:248` for ry, etc.).
For `rot` to be the branch's coefficient map one needs ε(br k) = −ε(k) — i.e.
`branchExp G (mulBits G P) = branchExp G P + 2` — and one needs the branch map to
be an involution on the anti set so that the coefficient space really decomposes
into 2-D planes. Neither statement exists in the file; the ε section (`:57-107`)
and the `rot` section (`:235-286`) never reference each other, and `rot`'s
docstring at `:252` only relates it to a *pure* P (v = (1,0)), the case where the
partner coefficient is 0 and the antisymmetry is invisible. So the three theorems
the crate cites everywhere for norm preservation, invertibility and Trotter
merging are proved about a matrix that has not been shown to be the one the code
applies. (This is the Tier-A root of G-031; closing one closes both.)

**Proposed closure.** Add to `Rotation.lean`:
`eps_antisymm : ∀ gx gz x z, omega gx gz x z = 1 → branchExp gx gz (mulBits gx gz x z).1 (mulBits gx gz x z).2 = branchExp gx gz x z + 2`
(by `decide`), plus `br_involutive` and `anti (br k)` (both by `decide` from
`mulBits`), then a theorem that `twoPass` restricted to a {k, br k} pair *is*
`rot θ (A k, A (br k))` — which finally makes
`rot_norm_sq`/`rot_neg_rot`/`rot_rot` statements about the kernel. Oracle: the
norm test at `pauli_sum_rotation_noise_lean.rs:80` already exercises this
incidentally; make it deliberate by seeding supports that contain both members of
every branching pair and asserting the exact 2x2 image, citing the new pairing
theorem.

### G-033 — The 15-index-list theorem — the adjudication of suspected old bug 6 — is pinned to nothing

- Class `bridge` · Tier A · Severity high · Sector noise-observables · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/noise.rs:178`
- Lean: `lean/PPVM/Algebra/Noise.lean:235` `twoQubitPauliError_indices_anticommuting` (with `:241` `_length`, `:244` `_nodup`)

**Claim.** Each of the 16 match arms of `Sum::two_qubit_pauli_error` scales by
1 − 2·Σ_{Q anti P} p_Q for its own observed pair P, in the documented order
{IX,IY,IZ,XI,…,ZZ}.

**Why unverified.** Lean genuinely re-derives the lists (good, independent
`decide`), and `Noise.lean:157-173` explains why this matters: the lists were
hand-written in legacy with no derivation and legacy's tests only probe one-hot
vectors. But grep for `twoQubitPauliError_indices` across all of `crates/`
returns ZERO hits — no `*_lean.rs` oracle mentions it. The only Rust-side test of
this surface is `pauli_sum_gate_surface_diff.rs:116`, a legacy diff. So the Lean
result adjudicates legacy but never touches the `-2` code. The risk is concrete
and not hypothetical: the `-2` arms are keyed on the bit-packed code (2=Z, 3=Y)
while Lean's `code`/`qPair` use alphabet order (2=Y, 3=Z), so arm↔list
association is exactly where a transposition would hide. All 15 were checked by
hand and they agree — that check should be a test, not an audit note.

**Proposed closure.** Add to `pauli_sum_rotation_noise_lean.rs` an exhaustive n=2
test: for each of the 16 pairs P, build the pure sum {P: 1.0}, apply
`two_qubit_pauli_error` with a *linearly independent* p (e.g. p[i] = 3^-i), and
compare the resulting coefficient against 1 − 2·Σ_{Q: antiPair(P,Q)} p[Q]
computed in the test from an independent ω over the documented pair order
(16×16 = 256 cases, cheap). That would fail on any single transposed index.

### G-034 — Loss-channel trace-preservation theorems have no oracle; the Rust explicitly defers to them

- Class `bridge` · Tier A · Severity high · Sector noise-observables · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/loss.rs:137`
- Lean: `lean/PPVM/Algebra/Noise.lean:450` `correlatedLossChannel_trace_preserving` (also `:374` `lossChannel_trace_preserving`, `:401` `resetLossChannel_trace_preserving`)

**Claim.** The shipped arms of `Sum::correlated_loss_channel` / `loss_channel` /
`reset_loss_channel` realize the transfer matrices `corrT`/`lossT`/`resetT`, and
are therefore trace preserving in the Heisenberg picture (Λ*(I+L) = I+L).

**Why unverified.** `loss.rs:113-115` says in so many words: "No test in old
covers correlated loss with distinct p[0]/p[1]/p[2], so both are unpinned there;
reproducing them keeps the golden master exact and leaves the adjudication (a
CPTP / trace-preservation theorem) to Lean." Lean delivered the theorem and even
refutes the crate's suspicion — but grep for
`trace_preserving|corrT|lossT|resetT` across `crates/` returns zero hits. The
only test is `pauli_sum_loss_diff.rs:92
correlated_loss_matches_all_four_loss_arms`, a legacy diff at a single
probability triple. So the *entire* justification chain is: Lean proves a
hand-transcribed matrix is TP; a diff test proves `-2` equals legacy; nothing
proves `-2` equals the matrix. Mutating `both_present` from 1−2p₁−p₀ to 1−p₁−p₀
would break the diff test only because legacy happens to agree, not because any
oracle knows the Lean value.

**Round-2 note (this row was not itself adjudicated; status stays `open`).** Two
statements above are superseded by the U1 adjudication and must not be carried into
a closure. (a) The mutation named here — `both_present` from 1−2p₁−p₀ to 1−p₁−p₀ —
is not a defect but precisely the adjudicated **fix**: `p[1]` means P(exactly one
lost), so `1 − p₁ − p₀` is the correct scale and the shipped `1 − 2p₁ − p₀` is the
bug. (b) Trace preservation is not the discriminator here: it holds column-wise for
*both* conventions and for every p. See [adj U1](#u1-correlated-loss-convention).

**Proposed closure — the `corrT` half is INVALIDATED by round 2.** The 𝟙-is-fixed
test stands and is worth landing as written: build 𝟙 = Σ_{k ∈ {I,L}^n} 1·k on a
`LossyPauliWord`-keyed sum for n = 1, 2, apply `loss_channel` /
`reset_loss_channel` / `correlated_loss_channel` at distinct irrational-ish p, and
assert the result equals 𝟙 exactly on the {I,L} sector, citing the three theorem
names — but note it passes today under *both* conventions, so it cannot be the
oracle that pins `correlated_loss_channel`. The "per-arm coefficient check against
`corrT` evaluated in the test" must **wait for `corrT` to be corrected** (G-040
item 6: `1 - p1 - p0` and `p1 / 2`); transcribing today's `corrT` would
machine-check the 2p₁ bug into place. The arm that actually discriminates is the
single-loss mass: assert P(exactly one lost) = p₁ across all three backends.

### G-035 — The sub-stochastic precondition Lean calls load-bearing is unenforced and untested in -2

- Class `bridge` · Tier A · Severity medium · Sector noise-observables · Status **adjudicated-defect** (round 2 — [adj U1](#u1-correlated-loss-convention))
- Rust: `crates/ppvm-tableau-2/src/mixture/noise/pauli.rs:96`
- Lean: `lean/PPVM/Algebra/Noise.lean:150` `eigenvalue_abs_le_one_needs_substochastic` (with `:93` `pauli_channel_eigenvalue_abs_le_one`)

**Claim.** The channel inputs satisfy p_Q ≥ 0 and Σ_Q p_Q ≤ 1, so the eigenvalue
is contractive and the skipped truncation / skipped weight re-check in
`Sum::scale_by_key` (cited at `crates/ppvm-pauli-sum-2/src/sum.rs:626-631`) is
sound.

**Why unverified.** `Noise.lean:145-149` states the obligation explicitly — "the
Rust channel constructors owe the precondition" — and no `-2` code pays it.
`Sum::pauli_error` (`crates/ppvm-pauli-sum-2/src/noise.rs:67`) has no assertion
of any kind; `two_qubit_pauli_error` has none; the tableau path asserts only
per-element membership in [0,1]
(`debug_assert!(p.iter().all(is_probability))`,
`crates/ppvm-tableau-2/src/noise.rs:90`) which permits Σp = 3; the mixture
asserts nothing and then computes `self.entries[parent].1 *= 1.0 - total`, which
goes NEGATIVE for total > 1 and is subsequently silently deleted by `truncate`
(`retain(|(_, prob)| *prob > self.sum_cutoff)`,
`crates/ppvm-tableau-2/src/mixture/data.rs:145`) and re-normalized — so an
over-normalized vector produces a plausible-looking but wrong distribution rather
than an error. `normalize_probabilities` also divides by a sum with no zero check
(`data.rs:136-140`), yielding NaN weights if everything truncates away. No test
exercises Σp > 1 on any backend.

**Round-2 amendment.** Confirmed unenforced, and the correct condition is now
stated: for correlated loss, CP ⟺ p₀, p₁ ≥ 0, p₀ + p₁ ≤ 1, p₂ ∈ [0,1] (the
sub-stochastic reading of the *named disjoint events*, not p₀ + 2p₁ ≤ 1 — see
G-040). Measured: p = [0.6, 0.6, 0] gives pauli-sum-2 coefficient −0.8 while the
mixture truncates its −0.2 survivor and renormalizes to (0.25, 0.25, 0.5); p =
[5, −3, 17] gives +2.0 and "both lost w.p. 1", both silently. The `debug_assert`
must live in the f64 backends: on a coefficient-generic `C` it will not compile
(no `PartialOrd`). See [adj U1](#u1-correlated-loss-convention).

**Proposed closure.** Add `debug_assert!(px+py+pz <= 1.0)` (and the 15-ary / loss
analogues) to the `-2` channel entry points, and a zero-norm guard in
`normalize_probabilities`. Oracle: a `*_lean.rs` test citing
`eigenvalue_abs_le_one_needs_substochastic` that (i) asserts |λ_P| ≤ 1 over a
randomized sub-stochastic sweep for every backend, and (ii) asserts the
`debug_assert` fires (`should_panic`) on Σp > 1.

### G-036 — Amplitude damping (the only non-unital branching channel) has no Lean model at all

- Class `coverage` · Tier A · Severity high · Sector noise-observables · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/noise.rs:312`
- Lean: none

**Claim.** `amplitude_damping(q, γ)` implements the Heisenberg adjoint of the
standard amplitude-damping channel: X,Y ↦ √(1−γ)·(X,Y), Z ↦ (1−γ)Z + γI, I ↦ I,
and is unital (hence its Schrödinger dual is trace preserving).

**Why unverified.** grep over all of `lean/` for
`amplitude|damping|kraus|cptp` returns zero hits in any PPVM module. This is the
sector's only non-Pauli, non-diagonal channel, the only one that branches a key
(Z ↦ I), the only one that needs a square root, and the only one whose
correctness is a genuine Kraus computation (K₀ = diag(1,√(1−γ)),
K₁ = √γ|0⟩⟨1|). The shipped arithmetic was verified by hand to be right, but the
ONLY Rust-side evidence in the repo is
`pauli_sum_gate_surface_diff.rs:155 amplitude_damping_matches_old`, i.e. legacy
agreement — which by the repo's own rules is evidence of nothing about
correctness. The unital-channel eigenvalue formula that `Noise.lean` does prove
provably does NOT cover it (Λ* here is not diagonal: it moves Z onto I).

**Proposed closure.** Add a Lean section deriving the transfer matrix from Kraus
operators over ℂ or ℤ[i][√]: define K₀,K₁, prove Σ Kᵢ†Kᵢ = 1 (TP of the dual) and
compute Kᵢ†·σ·Kᵢ for σ ∈ {I,X,Y,Z}, obtaining exactly the shipped
(√(1−γ), (1−γ), γ) coefficients; state `amplitudeDamping_unital`. Oracle: a new
`pauli_sum_amplitude_damping_lean.rs` pinning n=1 exhaustive over {I,X,Y,Z} at
several γ (including 0 and 1) against those Lean-derived coefficients, plus the
accumulate-onto-existing-I case (contract 3(c)).

### G-037 — Mixture loss/erasure has no formal statement, and its Mutation::Loss never collapses the qubit

- Class `coverage` · Tier A · Severity high · Sector noise-observables · Status **open**
- Rust: `crates/ppvm-tableau-2/src/mixture/noise/loss.rs:32`
- Lean: none

**Claim.** `GeneralizedTableauMixture::loss_channel` implements the same loss
channel as the trajectory `GeneralizedTableau::loss_channel`, as a two-branch
convex mixture (survivor at 1−p, lost branch at p).

**Why unverified.** Asked directly: loss/erasure in the mixture path has *no*
formal statement anywhere in `lean/` (grep for `mixture|erasure` in `lean/PPVM`
returns only unrelated prose). And the two implementations are not the same
channel: the trajectory path (`crates/ppvm-tableau-2/src/noise.rs:250`
`lose_qubit`) measures Z, applies X on outcome 1, and only then sets `is_lost`,
i.e. it projects and resets to |0⟩; the mixture's branch mutation is
`Mutation::Loss { qubit } => tab.is_lost[qubit] = true`
(`crates/ppvm-tableau-2/src/mixture/equality.rs:64`) with no projection at all. A
mixture that is supposed to be the trajectory average must split a random Z
outcome into two branches at 1/2 each; it does not, so the lost qubit retains
coherence in the mixture and loses it in the trajectory. The only evidence offered
is `tableau_mixture_diff.rs:142 noise_loss_and_reset_loss_match`, legacy
agreement.

**Proposed closure.** State in Lean the mixture as a Schrödinger stochastic map on
`C[Tableau]` (`Frame.lean` already has the key type) and prove (i) Σ weights is
preserved by `loss_channel` absent truncation, (ii) the loss branch equals the
1/2·(outcome-0) + 1/2·(outcome-1) projection mixture when Z_q is random and the
single projection when deterministic — which is the statement that currently
fails. Oracle: a `tableau_mixture_lean.rs` comparing the mixture's post-loss
⟨X_q⟩ against the trajectory ensemble average over many seeds.

### G-038 — The tableau-2 trajectory samplers are never shown to realize the channel probabilities

- Class `coverage` · Tier A · Severity medium · Sector noise-observables · Status **open**
- Rust: `crates/ppvm-tableau-2/src/noise.rs:74`
- Lean: none

**Claim.** `depolarize_impl`'s comparison ladder (`p <= r` return, then `p > 3r` ⇒
X, `p > 1.5r` ⇒ Y, else Z) applies each of X,Y,Z with probability exactly p/3, and
`pauli_error_impl`'s cumulative scan applies X,Y,Z with probabilities exactly
p_X,p_Y,p_Z in the trait's documented order.

**Why unverified.** These are the *stochastic* realizations of the same channels
the Lean file is about, and Lean says nothing about them beyond three
trivialities (`fire_conventions_agree_off_diagonal`, `fire_strict_zero_noop`,
`fire_nonstrict_fires_at_zero`) which only compare `<` against `≤`. The
`p > r*3.0 / p > r*1.5` encoding is non-obvious (it is an inverse-CDF written as
a multiply on the *draw* rather than a divide on p, so its correctness depends on
r ≥ 0 and on the earlier `p <= r` rejection); nothing states that the induced law
is Uniform-driven Bernoulli(p/3)³. Same for the 15-way `scan`/`position` inverse
CDF at `crates/ppvm-tableau-2/src/noise.rs:129` and its `PAULI_PAIRS[i+1]`
re-indexing, where an off-by-one would silently permute the whole 15-Pauli
alphabet. Only `tableau_behaviour_diff.rs` covers this, and it tests RNG *draw
counts*, not the induced distribution.

**Proposed closure.** In Lean, model each sampler as a function r ↦ Pauli on
[0,1) and prove the preimage of each outcome is an interval of the stated length
(`Set.Ioo` measure or just interval endpoints), i.e.
`depolarize_preimage_X = Ico 0 (p/3)` etc., and that `PAULI_PAIRS[i+1]` is the
documented order. Oracle: a `tableau_noise_lean.rs` histogram test over a fixed
seed asserting each outcome's empirical frequency within a Chernoff bound of the
Lean-stated interval length, plus an exhaustive check that the sampler on a swept
grid of r reproduces the Lean interval boundaries exactly.

### G-039 — Depolarizing's 1−4p/3 and 1−16p/15 have no Lean corollary and no oracle test

- Class `coverage` · Tier B · Severity low · Sector noise-observables · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/noise.rs:220`
- Lean: none (nearest: `lean/PPVM/Algebra/Noise.lean:57` `pauli_channel_eigenvalue_omega`, `:241` `twoQubitPauliError_indices_length`)

**Claim.** `depolarize1`'s factor 1 − 4p/3 applied to every non-identity Pauli at
the site, and `depolarize2`'s 1 − 16p/15 applied to every term non-identity on
either endpoint, are the transfer eigenvalues of the uniform depolarizing
channels.

**Why unverified.** Both are instances of the proved formula (uniform p/3 gives
1 − 2·2·p/3; uniform p/15 with exactly 8 anticommuting partners gives
1 − 2·8·p/15, and the count 8 is even available as
`twoQubitPauliError_indices_length`) — but neither instantiation is stated
anywhere in Lean, and no `*_lean.rs` oracle touches either method:
`pauli_sum_gate_surface_diff.rs:129/141` are legacy diffs, and the mixture's
`depolarize1`/`depolarize2`
(`crates/ppvm-tableau-2/src/mixture/noise/pauli.rs:112, :193`) just forward p/3
and p/15 with no check either. The non-obvious part is the *predicate*:
`depolarize2` fires on "non-identity on either endpoint", which is only the right
support because λ_P is 1 − 16p/15 for all 15 non-identity pairs uniformly and 1
for II — a fact worth one corollary. `crates/ppvm-pauli-sum-2/src/noise.rs:258`
is the 16/15 site.

**Proposed closure.** Add two Lean corollaries:
`depolarize1_eigenvalue : λ_P = 1 - 4*p/3` for P ≠ I under the uniform vector,
and `depolarize2_eigenvalue : λ_P = 1 - 16*p/15` for every P ≠ II (using the
exactly-8 count), plus λ = 1 on the identity. Oracle: extend
`pauli_sum_rotation_noise_lean.rs` with an exhaustive n=1 (4 keys) and n=2 (16
keys) check of `depolarize1`/`depolarize2` coefficients against those closed
forms.

### G-040 — corrT models pauli-sum's 2·p₁ convention; tableau-2 and the mixture implement p₁ — different channels

- Class `fidelity` · Tier A · Severity high · Sector noise-observables · Status **adjudicated-defect** (round 2 — [adj U1](#u1-correlated-loss-convention))
- Rust: `crates/ppvm-tableau-2/src/mixture/noise/loss.rs:95`
- Lean: `lean/PPVM/Algebra/Noise.lean:422` `corrT` / `:450` `correlatedLossChannel_trace_preserving`

**Claim.** The three shipped `-2` backends implement two different correlated-loss
channels: pauli-sum-2 (`loss.rs:137`) puts total mass 2p₁ on the single-loss
events, while tableau-2's mixture (`mixture/noise/loss.rs:95`) and trajectory
(`noise.rs:365`) both put p₁; `corrT` models only the former, and no Lean
docstring or trait doc flags the split.

**Why unverified.** It does not. Read the both-present branch three ways. (a)
`ppvm-pauli-sum-2/src/loss.rs:137` scales the survivor by `1 − 2·p₁ − p₀`, i.e.
P(exactly one lost) = 2p₁ split p₁/p₁ — and `corrT`'s (I,L)/(L,I) columns each
carry p₁, so Lean models exactly this. (b)
`crates/ppvm-tableau-2/src/mixture/noise/loss.rs:95` scales by `1 − p₀ − p₁` with
the two single-loss branches at `p₁/2` each, i.e. P(exactly one lost) = p₁. (c)
`crates/ppvm-tableau-2/src/noise.rs:366` runs a cumulative sum over `p[..2]` and
then a fair coin, also giving P(exactly one lost) = p₁. Backends (b),(c) are
consistent with each other and inconsistent with (a) by a factor of two in the
trait's own documented "p[1]: losing either one".
`correlatedLossChannel_trace_preserving` is therefore a spec for at most one of
the three, and `Noise.lean` nowhere flags the split — it presents `corrT` as
"`correlated_loss_channel(p₀,p₁,p₂)` … arm for arm". Additionally the
one-already-lost arms differ structurally: (a) emits a *recovery* branch weighted
p₁ (Heisenberg gain column) while (b),(c) emit a *loss* of the survivor weighted
p₂.

**Round-2 amendment + RULING.** The convention question is **settled in favour of
the paper**: `p[1]` is the probability a *named* one of the pair is lost, so
P(exactly one) = 2p₁. `ppvm-pauli-sum-2` and `corrT` are **correct**;
`ppvm-tableau-2`'s mixture and trajectory are the defect. How it got contested: Round 2 initially concluded `p[1]` = P(exactly one
lost) on the strength of `ppvm-traits-2/src/gates.rs:577`,
`ppvm-python/src/ppvm/mixins.py:507` and three shipped Python tests; the paper
draft (`../ppvm-paper/main.tex:462`, `:523`, `:845`) says the opposite —
$p_{LQ}$ is the probability that a *named* one is lost, so P(exactly one)
$=2p_{LQ}$, which is what `ppvm-pauli-sum-2` and `corrT` implement. Both readings
are the same CPTP family reparameterized by a factor of two, so mathematics does
not decide; the repo simply ships two conventions and documents each in a
different place. **Change no code and no `corrT` until this is ruled on.** Full
detail, including the two paper demos that appear to use different conventions,
is in the correction box at [adj U1](#u1-correlated-loss-convention).

What round 2 *did* settle here, independent of the ruling:
the "one-already-lost arms differ structurally" half of this row is **refuted**:
pauli-sum-2's gain branch is the Heisenberg transpose (T = Sᵗ) of the mixture's
Schrödinger loss column, both are correct, and the channel is trace-preserving
column-wise for *every* p. See [adj U1](#u1-correlated-loss-convention).

**Proposed closure — RULED, direction fixed.** Round 2's replacement closure
(`corrT` → `1 - p1 - p0`, `loss.rs` → match) is **discarded**: it pointed the
wrong way. Round 1's option (ii) is likewise discarded — there is one channel, not
two to be modelled separately. The closure is:

1. State the convention **normatively, in one place** — `ppvm-traits-2`'s
   `CorrelatedLossChannel` doc (`gates.rs:579`) — in the paper's words: `p[1]` is
   the probability that a *named* one of the two is lost, so the chance of exactly
   one loss is `2·p[1]` and the both-present survivor scales by `1 − 2·p[1] − p[0]`.
   Make `ppvm-tableau-2/src/noise.rs:337`, `ppvm-python/src/ppvm/mixins.py:507` and
   `paulisum.py:474` all cite that one wording.
2. Fix `ppvm-tableau-2/src/mixture/noise/loss.rs:95` (survivor `1 − p0 − 2·p1`,
   single-loss branches at `p1` each) and `ppvm-tableau-2/src/noise.rs:365` (the
   categorical scan must place `2·p1` on the exactly-one event). Rewrite
   `ppvm-python/test/generalized_tableau/test_loss.py:82`, `:90`, `:173`, using
   `p = [0, 0.5, 0]` where `:82` used the now-inadmissible `p = [0, 1, 0]`.
   `ppvm-pauli-sum-2/src/loss.rs` and `Noise.lean:422 corrT` **do not change**.
3. Oracle: assert **all three** backends agree on P(exactly one lost) `= 2·p[1]`,
   and add the paper's exact end-to-end prediction $m(t)=(1-p)^{k(t)}$ from
   §`sec:transport` (with $[p_{LL},p_{LQ},p_{LN}]=[p/3,p/3,p/3]$) as an integration
   test — a closed-form number that discriminates the two conventions directly and
   would have caught this split on day one.

### G-041 — The firing-convention section's headline claim is false for the -2 code (p ≤ 0 is guarded first)

- Class `fidelity` · Tier B · Severity medium · Sector noise-observables · Status **open**
- Rust: `crates/ppvm-tableau-2/src/noise.rs:270`
- Lean: `lean/PPVM/Algebra/Noise.lean:289` `fire_nonstrict_fires_at_zero` (section at `:247-273`)

**Claim.** "So `loss_channel(q, 0.0)` is *not* the identity under the shipped
convention, while `depolarize1(q, 0.0)` is" (`Noise.lean:266-267`), presented as
"the whole observable content of the loss_channel convention divergence".

**Why unverified.** The shipped `-2` `loss_channel` reads
`if p <= 0.0 || self.is_lost[qubit] { return; }` *before*
`if p < rng.random::<f64>() { return; }`, so p = 0.0 returns early and IS the
identity; likewise `asymmetric_loss_channel` guards `p_tot <= 0.0 ||`
(`crates/ppvm-tableau-2/src/noise.rs:319`). The described divergence cannot be
exhibited by the code the section claims to be about. Compounding it: the section
cites `crates/ppvm-tableau/src/{noise.rs,tableau_like.rs}` (legacy paths; in `-2`
both live in `ppvm-tableau-2/src/noise.rs`), it attributes the
`if p <= r { return }` form to `pauli_error_impl`/`two_qubit_pauli_error_impl`
which do not contain that comparison at all (they use cumulative `> r` scans), and
the two "theorems" are `not_lt.mpr` and `le_refl 0` — a `strength` problem on top
of the fidelity one, since `(0:ℝ) ≤ 0` is billed as machine-checking a
behavioural divergence.

**Proposed closure.** Either delete the section or restate it truthfully: model
each channel's *guarded* firing predicate (`0 < p ∧ r < p` vs `0 < p ∧ r ≤ p`) and
prove they agree everywhere on [0,1) except {r = p}, so the shipped inconsistency
has no observable content at all; re-cite `ppvm-tableau-2/src/noise.rs` with the
correct comparison per method. Oracle: the existing zero-probability draw-count
tests in `tableau_behaviour_diff.rs` should move to a `*_lean.rs` and cite the
restated theorem.

### G-042 — λ_P is an arithmetic identity with hand-inserted signs, never derived from Q P Q† = (−1)^ω P

- Class `strength` · Tier A · Severity high · Sector noise-observables · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/noise.rs:73`
- Lean: `lean/PPVM/Algebra/Noise.lean:57` `pauli_channel_eigenvalue_omega` (and `:43` `pauli_channel_eigenvalue`)

**Claim.** The unital Pauli channel ρ ↦ Σ_Q p_Q QρQ† acts on the observable P by
multiplication by λ_P = Σ_Q p_Q(−1)^{ω(P,Q)}, so `pauli_error_factors`'s
[1−2(p_Y+p_Z), 1−2(p_X+p_Y), 1−2(p_X+p_Z)] are that channel's transfer
eigenvalues.

**Why unverified.** `pauli_channel_eigenvalue` proves only
`Σ_Q p_Q·(if anti P Q then -1 else 1) = 1 − 2·Σ_{anti} p_Q` — the
`(if anti then -1 else 1)` factor is *inserted by hand*, and the docstring
concedes it (`Noise.lean:55-56`: "the channel superoperator and its
diagonalization are not constructed here — only the eigenvalue's algebraic
form"). Nothing in `lean/` proves the physics step Q P Q† = (−1)^{ω(P,Q)}·P, nor
that Λ* is diagonal in the Pauli basis, nor that Λ is CP/TP. So the theorem is
true of ANY sign function called `anti`; it cannot distinguish the real Pauli
channel from one with a transposed commutation rule. Given the sector's whole
factor set (4p/3, 16p/15, the 15 index lists) hangs off this single identity, it
is the load-bearing unproved link. The tools exist and are unused:
`PPVM.PauliMatrix.tensorPauli_mul` gives the phase of a Pauli word product in
ℤ[i] and `PauliPhase.phaseExp_sub_comm` gives P·Q = (−1)^ω Q·P.

**Proposed closure.** In `Noise.lean` (or `Matrix.lean`) state and prove
`tensorPauli q * tensorPauli p * tensorPauli q = (-1)^(omega p q) • tensorPauli p`
from `tensorPauli_mul`, then define the superoperator
`Λ A = Σ_Q p_Q · (q·A·q)` on `CMap (Word n) GaussianInt` and prove
`Λ (single P c) = single P (λ_P * c)` with λ_P the existing expression. Oracle:
extend `pauli_sum_rotation_noise_lean.rs` to compute λ_P independently from the
matrix trace Tr(Λ*(P)·P)/2ⁿ over n=1,2 exhaustively and compare against the
shipped `pauli_error`/`pauli_error_many` output.

### G-043 — No positivity anywhere: the loss trace-preservation theorems hold for p₀ = 5 as well

- Class `strength` · Tier A · Severity medium · Sector noise-observables · Status **adjudicated-defect** (round 2 — [adj U1](#u1-correlated-loss-convention))
- Rust: `crates/ppvm-pauli-sum-2/src/loss.rs:82`
- Lean: `lean/PPVM/Algebra/Noise.lean:450` `correlatedLossChannel_trace_preserving`

**Claim.** The loss channels are valid CPTP maps (the sector's channels are
channels, not merely trace-preserving linear maps).

**Why unverified.** `correlatedLossChannel_trace_preserving` is advertised as
strong precisely because it needs "no normalization hypothesis at all"
(`Noise.lean:318-320, 438-439`). That is also its weakness: it is a
linear-algebra identity in ℝ that holds for p₀ = 5, p₁ = −3, p₂ = 17. Trace
preservation is one of the two CPTP conditions; complete positivity
(equivalently, that the arms come from a Kraus family / that the dual stochastic
matrix is entrywise nonnegative and column-summing to 1) is never stated for ANY
channel in the sector — grep for `kraus|cptp|positiv` over `lean/PPVM` finds
nothing relevant. So the file cannot distinguish the shipped loss channel from a
non-physical map with the same column sums, and the Rust side likewise never
validates the loss probabilities (`loss.rs:82/127` assert nothing).

**Round-2 amendment.** Adjudicated and confirmed: CP ⟺ every event weight ≥ 0
(the map is then a convex mixture of id, loss0, loss1, lossboth). With G-040 ruled in
favour of the paper, the necessary *and sufficient* admissibility condition is
**p₀ + 2p₁ ≤ 1**, p₀, p₁ ≥ 0, p₂ ∈ [0,1] — the original ledger hypothesis,
restored. (Round 2 had replaced it with p₀ + p₁ ≤ 1 under the reading the ruling
rejected.) Measured: p = [5, −3, 17] is accepted silently by both backends
(pauli-sum-2 returns coefficient +2.0; the mixture truncates its negative survivor
and renormalizes to "both lost with probability 1"). See
[adj U1](#u1-correlated-loss-convention).

**Proposed closure — original hypothesis restored by the G-040 ruling.** Use
`0 ≤ pᵢ`, `p₀ + 2p₁ ≤ 1`, `p₂ ≤ 1`. The rest of the closure stands: prove
the transpose of `corrT`/`lossT`/`resetT` is column-stochastic (entries ≥ 0, each
column summing to 1), i.e. the channel is a convex mixture of loss events, hence
CP; extend the proposed `pauli_sum_loss_lean.rs` (G-034) to assert every arm
coefficient is in [0,1] for admissible p; and add the `debug_assert`s at all three
entry points (`ppvm-pauli-sum-2/src/loss.rs:133`,
`ppvm-tableau-2/src/mixture/noise/loss.rs:59`, `ppvm-tableau-2/src/noise.rs:352`)
— noting that on a coefficient-generic `C` the assert will not compile as written
(no `PartialOrd`), so it must live in the f64 backends.

### G-044 — Zero-state read-out assumes ⟨0|P|0⟩ = [P is X-free] and is pinned by no oracle

- Class `strength` · Tier A · Severity medium · Sector noise-observables · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/trace.rs:93`
- Lean: `lean/PPVM/Algebra/Noise.lean:490` `overlap_with_zero_xfree` (and `:479` `overlap_with_zero`)

**Claim.** `Trace::trace(&PauliPattern::zero_state())` computes
⟨0ⁿ|O|0ⁿ⟩ = Σ_{P X-free} c_P.

**Why unverified.** The theorem is
`Σ_P c_P·(if diag P then 1 else 0) = Σ_{diag} c_P`, i.e. `Finset.sum_filter` — a
rearrangement. The entire physics content, that ⟨0ⁿ|P|0ⁿ⟩ is 1 exactly for X-free
P and 0 otherwise, is inserted as the indicator; the docstring itself calls it
"the *modeling assumption*" (`Noise.lean:486-489`), which is precisely the repo's
own definition of a strength gap. This is gratuitously assumed:
`PPVM.PauliMatrix.tensorPauli` gives the genuine ℤ[i] matrix, so ⟨0ⁿ|·|0ⁿ⟩ is the
(0,0) entry and the claim is a one-line computation there (the file already does
the harder `trace_tensorPauli_mul`). Separately, no `*_lean.rs` oracle pins the
Rust read-out at all: `zero_state()`/`trace` appear only in
`pauli_sum_loss_diff.rs:130`, `column_store_diff.rs:1361` and `sym_diff.rs` —
legacy diffs. Nor does anything state what the pattern does on a lossy word's L
site.

**Proposed closure.** Prove
`tensorPauli p 0 0 = if (∀ i, (p i).1 = false) then 1 else 0` in `Matrix.lean` and
derive `overlap_with_zero_xfree` from it (dropping the assumption), then extend to
the {I,X,Y,Z,L} alphabet for the lossy read-out. Oracle: a `*_lean.rs` test that,
for n ≤ 3, builds every Pauli word as a one-term sum, contracts against
`PauliPattern::zero_state()`, and compares with an independently computed
⟨0ⁿ|P|0ⁿ⟩ from explicit 2ⁿ×2ⁿ matrix products.

### G-045 — asymmetric_loss_channel's state-dependent p_tot has no Lean model and is tested only at p0 == p1

- Class `coverage` · Tier A · Severity medium · Sector products-and-channels (skeptic) · Status **adjudicated-spec** (round 2 — [adj U1](#u1-correlated-loss-convention))
- Rust: `crates/ppvm-tableau-2/src/noise.rs:317`
- Lean: none

**Claim.** `AsymmetricLossChannel::asymmetric_loss_channel(q, p0, p1)` loses `q`
with the correct total probability
`p_tot = p0·P(|0⟩) + p1·P(|1⟩) = p0·(1+⟨Z_q⟩)/2 + p1·(1−⟨Z_q⟩)/2`, and its
documented omission of the survival back-action
`K₀ = √(1−p0)|0⟩⟨0| + √(1−p1)|1⟩⟨1|` is an approximation whose error is bounded
(it is exact only at p0 == p1).

**Why unverified.** This is the sector's only *state-dependent* noise channel —
the loss rate is read out of the tableau via `z_expectation` at `noise.rs:316` and
mixed at `:317` — and it is the only channel in the `-2` crates whose own
docstring (`noise.rs:290-295`) declares itself an approximation of a non-Clifford
Kraus map, yet nothing anywhere characterizes or bounds that approximation.
`lean/` has no counterpart: grep for `asymmetric_loss` over `lean/PPVM/` hits
exactly one line, `Noise.lean:257`, and only in the firing-convention prose
(which G-041 confirms is wrong about this method — `noise.rs:319` guards
`p_tot <= 0.0` before the draw). `Noise.lean`'s Loss section models `lossT`
(`:364`), `resetT` (`:392`) and `corrT` (`:422`) and stops there. The Rust side is
no better: the ONLY test that calls it anywhere in the repo is
`tableau_behaviour_diff.rs:379-383
asymmetric_loss_pollutes_the_record_on_both_engines` at `(1.0, 1.0)` — a legacy
diff, at the one point of the parameter space where the whole z-dependence cancels
identically (p_tot = p0 = p1 = 1) and the channel saturates.
`ppvm-traits-2/tests/phase1_gate_surface.rs:804` calls it at (0.01, 0.02) but
through a hand-written stub impl at `:741`, not the real tableau. So the p0 ≠ p1
arithmetic that is the entire point of the channel — and the ⟨Z⟩ read-out it
depends on — is exercised by nothing, in any crate, at any parameter.

**Round-2 amendment.** The arithmetic is **correct**: the two-Kraus loss family
gives loss-branch trace p0(1+⟨Z⟩)/2 + p1(1−⟨Z⟩)/2, verbatim `noise.rs:365`, and the
⟨Z⟩ sign convention was probed (`z_expectation(|0⟩) = +1`, `z_expectation(|1⟩) = −1`,
and `asymmetric_loss_channel(p0=1, p1=0)` fires on |0⟩ and not on |1⟩, so p0 is the
|0⟩-rate). One sharpening for the closure: the dropped survival back-action makes the
shipped channel exact iff p0 = p1 **or the site is in a Z eigenstate**, with error
first order in |p0−p1| (the conditional survivor's Bloch-z biased by
≈ (p1−p0)/2·(1−⟨Z⟩²)/2) — so `noise.rs:290-295`'s "exact only at p0 == p1"
understates the exactness region. Coverage plus a docstring nit, not a wrong number.
See [adj U1](#u1-correlated-loss-convention).

**Proposed closure.** Lean: model the channel's total loss probability from the
two-Kraus loss family (K_L0 = √p0·|L⟩⟨0|, K_L1 = √p1·|L⟩⟨1|) and prove
`pLoss ρ = p0·(1+tr(Zρ))/2 + p1·(1−tr(Zρ))/2`, i.e. that the shipped
`p0*0.5*(1+z) + p1*0.5*(1-z)` is the trace of the loss branch; then state the
omitted survival back-action explicitly as
`K₀ = √(1−p0)|0⟩⟨0| + √(1−p1)|1⟩⟨1|` and prove `K₀†K₀ + Σ K_Li†K_Li = 1` so the
size of the dropped term is on the record (it vanishes iff p0 = p1). Oracle: a
`tableau_noise_lean.rs` test that, for a stabilizer state with known
⟨Z_q⟩ ∈ {+1, −1, 0} (|0⟩, |1⟩, and H|0⟩) and p0 ≠ p1, asserts the empirical loss
frequency over a seeded ensemble is within a Chernoff bound of the Lean-stated
p_tot — three points that pin the z-coefficient and would fail on either a
swapped p0/p1 or a dropped ⟨Z⟩ term, both of which the current suite accepts.

### G-046 — apply_producer is cited to two GradedMap theorems but pinned by no `*_lean.rs` oracle; two backends untested

- Class `bridge` · Tier A · Severity medium · Sector sum-engine-stores · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/column_store/lifecycle.rs:151`
- Lean: `lean/PPVM/Algebra/GradedMap.lean:200` `pushforward_eq_reset_accumulate` (and `:226` `merge_without_reset_ne_pushforward`)

**Claim.** `ApplyProducer::apply_producer` — produce the whole support, reset, then
accumulate — computes the pushforward `mapDomain φ ∘ mapRange g`, and the reset is
not optional.

**Why unverified.** `store.rs:578-584` names both theorems as the licence for the
reset, and `column_store/mod.rs:66-68` repeats the citation. But no file in
`crates/ppvm-conformance-2/tests` exercises `Sum::apply` at all (grepped for
`apply(`, `RekeyProducer`, `apply_producer` across all nine `*_lean.rs` and all
`*_diff.rs`). The only coverage is four hand-written examples in
`ppvm-pauli-sum-2/tests/engine.rs:179-226` — outside the oracle set, citing no
theorem, a single width, and only the HashMapStore backend.
`ColumnStore::apply_producer` (`column_store/lifecycle.rs:151`, which additionally
has to skip tombstoned rows and re-add through the accumulating `Columns::add`)
and `IndexMapStore::apply_producer` (`indexmap_store/lifecycle.rs:42`) are reached
by no test in the workspace: the whole producer path is latent because all four
hot gate families bypass it, so a mutation that dropped the `reset` in either
backend would be caught by nothing.

**Proposed closure.** Add to `pauli_sum_lean.rs` and `column_store_lean.rs` a
`pushforward_eq_reset_accumulate` test running `Sum::apply(RekeyProducer::new(..))`
on all three backends: exhaustive over the 4^n single- and two-qubit supports with
φ each single-qubit Clifford bit map (injective — check `pushforward_apply`
pointwise: the coefficient at φ(k) equals g(A k)), plus a deliberately
non-injective φ to pin the accumulating merge, plus the
`merge_without_reset_ne_pushforward` witness (assert that after `apply` the
pre-image key is ABSENT, which is what fails if the reset is removed). Seeded
proptest above n=2.

### G-047 — No abstraction function store → (K →₀ C) and no refinement theorems for any backend

- Class `coverage` · Tier B · Severity medium · Sector sum-engine-stores · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/column_store/columns.rs:7`
- Lean: `lean/PPVM/Algebra/GradedMap.lean:51` (`abbrev CMap` — the only object in the file; no store structure exists)

**Claim.** Each concrete store (ColumnStore's struct-of-arrays
`Columns{keys, coeffs, hashes, live, live_len, sparse_rows, live_runs, index}`,
HashMapStore's `{primary, aux, scratch, batch}`, IndexMapStore's ordered twin)
satisfies a representation invariant, and there is an abstraction map
α : store → (K →₀ C) such that every store operation commutes with α:
α(add k c s) = accumulateTerm (k,c) (α s), α(retain p s) = retain p (α s),
α(compact s) = α s, α(rekey f s) = pushforward f (α s), α(reset s) = 0.

**Why unverified.** `GradedMap.lean` is 725 lines about `Finsupp` and contains no
store type at all — grepped all of `lean/PPVM` for
`ColumnStore`/`tombstone`/`live`/`IndexMap`/`refine`; the only hits are prose in
docstrings ("the mathematical object every Sum backend refines",
`GradedMap.lean:50`). So none of the SoA machinery is modeled: the `live`
tombstone mask, the invariant `live.iter().filter(!=0).count() == live_len`
(`columns.rs:360 debug_assert_valid`), the `index : HashTable<(u32,u64)>` →
physical-row agreement, resurrection-by-append (`columns.rs:133-147` repoints a
stale bucket at a freshly pushed row while the dead row keeps its key and hash),
the 1/8-dead stable `compact()` + `reindex()` (`columns.rs:285-318`), the
lazily-rebuilt `sparse_rows`/`live_runs` cache, or the HashMapStore "aux/scratch/
batch are empty between operations" invariant that `Clone` and `PartialEq` rely on
(`store.rs:63-67`). The only artifact that even names the invariant is the ORACLE
test, whose own docstring concedes "the ONE invariant the layout genuinely
introduces, which no existing lemma states" (`column_store_lean.rs:22-23`). A
theorem about `Finsupp` cannot be violated by a missing `reindex()`; this is the
whole "spec ⇄ implementation" claim and it is absent.

**Proposed closure.** Add `lean/PPVM/Algebra/Store.lean`: a
`structure ColumnsModel K C` with `keys : List K`, `coeffs : List C`,
`live : List Bool` plus `Valid` (equal lengths; live keys pairwise distinct; index
modeled as a partial map row-lookup agreeing with the live rows), an
`abs : ColumnsModel → CMap K C` defined as the sum of live singletons, and
commuting lemmas `abs_add`, `abs_insert`, `abs_retain`, `abs_compact` (compaction
is α-invariant AND preserves live-row order), `abs_reindex`, `abs_rekey`, plus
`Valid` preservation for each. Bridge: extend `column_store_lean.rs` with a
debug-only accessor (or a `#[cfg(test)]` reflection hook) that exports (keys,
coeffs, live, index rows) and checks each commuting lemma exhaustively for n≤2 and
on seeded proptest above, citing the theorem names.

### G-048 — Pair::probe_batch has no Lean model and no test; ColumnStore trusts an unstated KeyBatch column invariant

- Class `coverage` · Tier B · Severity low · Sector sum-engine-stores · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/column_store/graded.rs:92`
- Lean: none

**Claim.** `Pair::probe_batch(keys, out)` sets `out[i] = get(keys[i])` for every i
— the "coalesced gather" the design advertises as the columnar read side.

**Why unverified.** `GradedMap.lean`'s L3 section mentions "the property the
batched probe relies on" but never defines a batched probe, so nothing states the
pointwise spec. On the Rust side `probe_batch` is implemented three times
(`store.rs:971`, `column_store/graded.rs:92`, `indexmap_store/algebra.rs:59`) and
appears in no test in `ppvm-pauli-sum-2/tests` or `ppvm-conformance-2/tests`
(grepped; the only test hits are in `ppvm-traits-2`'s own container tests). The
columnar impl also silently depends on an invariant nobody states or checks: it
uses the caller's precomputed digest column when
`hashes.len() == keys.keys().len()` and otherwise recomputes, so it is only
correct while `KeyBatch`'s parallel columns satisfy
`hashes()[i] == keys()[i].key_hash()` (asserted only in the doc comment of
`KeyBatch::fill_hashes`, `crates/ppvm-traits-2/src/batch.rs:292`). A stale digest
column makes every probe miss silently.

**Proposed closure.** State
`probeBatch (f : CMap K C) (ks : List K) : List (Option C) := ks.map (fun k => if k ∈ f.support then some (f k) else none)`
in `GradedMap.lean` with `probeBatch_get` (pointwise agreement with `get`) and
`probeBatch_perm`/`probeBatch_append` (gather/partition invariance). Oracle: add a
`probe_batch` test to `pauli_sum_lean.rs` and `column_store_lean.rs` comparing
`probe_batch` against per-key `get` for all three backends, exhaustively over the
4^n keys at n≤3 with and without `fill_hashes` called, with and without tombstones
present, citing `probeBatch_get`.

### G-049 — Lean models the zero-free canonical support; every -2 store keeps explicit zeros and len/eq observe them

- Class `fidelity` · Tier A · Severity high · Sector sum-engine-stores · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/store.rs:910`
- Lean: `lean/PPVM/Algebra/GradedMap.lean:64` (`def len := f.support.card`, docstring "`Support::len` — the size of the canonical (zero-free) support")

**Claim.** `Support::len` is the size of the canonical zero-free support, and
`Sum`'s equality/approx-equality are the Finsupp ones.

**Why unverified.** The Rust is the exact opposite by explicit design:
`Support::len` returns `primary.len()` / `live_len`, counting rows whose
coefficient is exactly zero, and `ops.rs:11-18` declares this a *user-facing
contract* (`state == (state2 *= 0.0)`), with every fast path documented as never
dropping a zero. Nothing in Lean models a non-canonical container, so no theorem
can state what `len()`, `PartialEq` or `AbsDiffEq` compute. The oracle makes the
mismatch executable: `column_store_lean.rs:526` asserts `zeroed.len() == 1` after
`scale(&0.0)`, while the Lean it cites in its own header gives `scale 0 f = 0`
with `len = 0`. `Truncation.lean:285-289` acknowledges the divergence in prose
("whether such a key is *listed* is that backend's reduce question") and then
never models it — so the acknowledgement removes the claim rather than proving it.

**Proposed closure.** Model the stored support as a pair
`(A : CMap K C, S : Finset K)` with `A.support ⊆ S` (S = the rows physically
present), define `lenStored := S.card`, `eqStored`, and restate
`reduce_structural` as `reduce (A,S) = (A, A.support)`, plus per-operation
`S`-propagation lemmas (accumulate/scale/scaleByKey/rotate add or keep keys; only
reduce/retain remove). Then a concrete corollary: `Sum::truncate` with preserve-set
P and a preserved key at coefficient 0 keeps the ROW (Rust) though the widened
retain has it off support — pin it with an oracle case
`n=1, preserve={X}, terms {X:0.0, Z:1.0}, CoefficientThreshold(0.5)`, asserting
`len()==2`, citing the new stored-support theorem.

### G-050 — ColumnStore's shipped rx kernel is a pair-fused interleaved walk, not the twoPass model it cites

- Class `fidelity` · Tier A · Severity high · Sector sum-engine-stores · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/column_store/rotations/rx.rs:127`
- Lean: `lean/PPVM/Instantiations/Rotation.lean` `accumulate_rotBatch` (and `eagerWalk_ne_twoPass`)

**Claim.** Every backend's rotation is the two-pass walk
`twoPass = diagPass + branchPass` — "scale ALL diagonals before merging ANY
branch" — which `accumulate_rotBatch` proves equals the one-pass produced batch;
`eagerWalk_ne_twoPass` shows an interleaved variant diverges.

**Why unverified.** `rotate_x_kernel` does not implement `twoPass`. When the
support is "closed" it walks rows, and for an anticommuting row `i` whose partner
row `j = toggle_x(i)` is present it writes BOTH final coefficients inside the walk
(`coeffs[i] = ci·cos + cj·sin·ε_j; coeffs[j] = cj·cos + ci·sin·ε_i`,
`rx.rs:127-129`, and again in `rotate_x_small_row` `rx.rs:173-174`), i.e. it
mutates a row it has not yet visited and resolves that row's branch merge before
its diagonal pass — precisely the interleaving `eagerWalk_ne_twoPass` exhibits as
wrong. It happens to be correct only because of four facts that appear nowhere in
Lean and nowhere in a Rust assertion: `br` (x-bit toggle) is an involution,
`anti(br k) = anti(k)` (z-bit and lostness are preserved), `br` is injective so two
unpaired branches never collide, and the `visited`/`j < i` guards make each pair
fire once. The `column_store/mod.rs:62-65` docstring asserts the opposite of the
code ("RotateInPlace keeps the two-pass ordering") while citing
`eagerWalk_ne_twoPass`; that is a citation over differing code. Additionally the
`original_rows > 512` monomorph (the `closed_support` probe + `visited` array path)
is reached by no test in the conformance suite — `column_store_lean.rs` uses 40-
and 64-term supports.

**Proposed closure.** State `pairedPass` in `Rotation.lean` (walk keys in any
order; on an anticommuting key whose `br`-image is in the support, write both
coefficients from the pre-walk values and mark both visited; else buffer the
branch) and prove `pairedPass = twoPass` from explicit hypotheses
`Function.Involutive br`, `∀ k, anti (br k) ↔ anti k`, `Function.Injective br`,
`s (br k) = ε-conjugate of s k` — plus a counterexample showing the theorem fails
if `anti (br k) ↔ anti k` is dropped. Oracle: add to `column_store_lean.rs` a test
that, for supports crossing 512 physical rows and with tombstones present,
`rx(q, θ)` via `ColumnStore` agrees bit-for-bit with the generic
`RotateInPlace::rotate_in_place` closure form on the same store, citing
`pairedPass_eq_twoPass`. (Shares its root with G-031/G-032.)

### G-051 — IndexMapStore's observable insertion order is justified only by a diff against the legacy indexmap backend

- Class `provenance` · Tier B · Severity medium · Sector sum-engine-stores · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/indexmap_store/branching.rs:23`
- Lean: none

**Claim.** The ordered backend's term order is correct: replacement keeps a key's
position, first insertion appends, retain preserves order, and the multi-branch
merge direction is chosen by comparing DEDUPLICATED map cardinalities rather than
the raw fan-out count.

**Why unverified.** `IndexMapStore` exists precisely because "term order is
observable" (`indexmap_store/mod.rs:5`), yet Lean has no ordered-support object
whatsoever — `CMap` is a `Finsupp`, and `accumulateTerms_perm` says order is
*irrelevant*, so the model cannot even express the property the backend is for.
The only conformance coverage is `pauli_sum_indexmap_diff.rs` (196 lines) which
compares against `ppvm_pauli_sum::config::indexmap` — legacy. The load-bearing
decision is reconstructed from legacy behaviour by the code's own admission:
"Legacy `consume` compares map cardinalities, not the raw fan-out count;
duplicate branch keys must not force branch-first ordering"
(`branching.rs:23-25`). That rule is not the same as HashMapStore's
(`if scratch.len() > primary.len()`, `store.rs:1551`), so for the same input the
two backends can pick opposite merge directions — the backends are demonstrably
NOT observationally equivalent on order, contradicting the "a backend swap is
observationally a no-op" claim at `column_store/mod.rs:49-51`.

**Proposed closure.** Add an ordered model to Lean
(`OrderedMap K C := List (K × C)` with distinct keys, an `abs` to `CMap`, and
`insertOrdered`/`addOrdered`/`retainOrdered`/`mergeSmallerFirst` definitions) and
prove: `abs` commutes with each op; `retainOrdered` is order-preserving;
`mergeSmallerFirst` with the deduplicated-cardinality rule yields a specified
order. Oracle: an `indexmap_store_lean.rs` that pins the exact term ORDER (not
just the key set) of `+=`, `extend`, `truncate`, `branch_in_place` and the
correlated loss channel against those definitions, exhaustively for n≤2, citing
them by name — replacing the legacy diff as the source of truth.

### G-052 — eval's ring-hom property is proved in Lean but never pinned to a map-backed Rust Term

- Class `bridge` · Tier A · Severity medium · Sector symbolic-coefficients · Status **open**
- Rust: `crates/ppvm-sym-2/src/mul.rs:227`
- Lean: `lean/PPVM/Instantiations/Symbolic.lean:467` `evalC_mul` (and `:270` `evalHom_mul`)

**Claim.** No `*_lean.rs` oracle states eval/eval_complex multiplicativity as a
law over the map-backed product surface (only two hardcoded instances at
`sym_lean.rs:300-307` and `335-363` exist), and `evalC_mul`/`evalHom_mul` are
`map_mul` on `AddMonoidAlgebra` rather than statements about the shipped
`mulImpl`.

**Why unverified.** `evalC_mul` / `evalHom_mul` are `map_mul` on genuine `AlgHom`s
and are cited in `eval.rs:225` and `coeff.rs:258` as the hom laws. But
`grep -rn 'evalC_mul|evalHom_mul' crates/` hits only docstrings: no `*_lean.rs`
test names them, and no test asserts multiplicativity of `eval` on a `Term`.
`sym_lean.rs`'s `symbolic_accumulation_laws_hold_denotationally` (line 634) checks
only additive comm/assoc, scalar scale_scale and the additive identity;
`symbolic_rotation_laws_hold` (line 691) does check eval of products, but only of
single atoms, i.e. only the `One × One` non-allocating arm. The whole map-backed
`Sum × Sum` / `Sum × One` product surface is covered only by `sym_diff.rs:225`
`ring_surface_matches_old_on_seeded_random_expressions`, which is
agreement-with-legacy — and legacy's `Sum × Sum` arm is where one of the two
confirmed bugs lived (`lean/README.md:179`, the signed `s2.c0 > min_eps` gate),
fixed by hand on this same branch, so the diff baseline is not independent
evidence.

**Proposed closure.** Add to `sym_lean.rs` a seeded property test over
`random_term` restricted to map-backed forms (force `+= 0.0` promotion) at
`max_sin = usize::MAX` / `min_eps = 0.0`: assert
`|(a.clone()*b.clone()).eval(&θ) - a.eval(&θ)*b.eval(&θ)| < tol`, the
`eval_complex` analogue, plus mul-associativity, distributivity over `+`, and
`a * Term::one()` denotational identity — naming `evalC_mul`/`evalHom_mul`. On the
Lean side add `evalC θ (multiply x y) = evalC θ x * evalC θ y` restated over
`mulImpl` at the untruncated bound so the theorem the test cites is about the
modelled implementation, not only about `AddMonoidAlgebra`'s `*`.

### G-053 — mulMono_drop_at_insert_eq_drop_at_end has no oracle test pinning the Rust to it

- Class `bridge` · Tier A · Severity medium · Sector symbolic-coefficients · Status **open**
- Rust: `crates/ppvm-sym-2/src/term.rs:376`
- Lean: `lean/PPVM/Instantiations/Symbolic.lean:173` `mulMono_drop_at_insert_eq_drop_at_end`

**Claim.** For the `max_sin` axis, dropping each produced monomial inside
`Sum::add_term` yields exactly the map obtained by forming the full product and
truncating afterwards — the exactness that makes drop-at-accumulate defensible
rather than an approximation.

**Why unverified.** This is the positive half of the truncation story and is cited
in `term.rs:337-347` and `mul.rs`, but no `*_lean.rs` test names it or checks it.
The only Rust-side evidence is `sym_diff.rs:1335 integration_sym_truncation_sweep`,
which asserts (a) old == new and (b) monotonicity of the retained monomial set in
k — a strictly weaker property that a buggy filter (e.g. testing `pow()` instead of
`sin_pow()`, or testing the multiplier's degree instead of the product's) would
still satisfy. So the theorem's Rust counterpart would survive mutation of the very
predicate it is about.

**Proposed closure.** Add a `sym_lean.rs` test: for seeded map-backed a,b at
`min_eps = 0.0`, compute `a.clone()*b.clone()` with `max_sin = usize::MAX`, filter
the resulting monomial table by `p.sin_pow() <= k` using `Term::inner()`/`iter()`,
and assert it equals (key-by-key, coefficient-by-coefficient) the table from the
same product computed with `max_sin = k` — naming
`mulMono_drop_at_insert_eq_drop_at_end`, and restricted to map-backed operands so
`mulImpl_not_wellDefined`'s caveat is respected.

### G-054 — The exact/inexact split's negative half (no halving in Z[i]) is prose only, and GaussianInt is i64-wrapping

- Class `coverage` · Tier A · Severity low · Sector symbolic-coefficients · Status **adjudicated-defect** (round 2 — [adj U6](#u6-unenforced-preconditions))
- Rust: `crates/ppvm-sym-2/src/exact.rs:22`
- Lean: none

**Claim.** `Halvable` must be split off `Coefficient` because ℤ[i] admits no
`half`: there is no x ∈ ℤ[i] with x + x = 1 + i — the claim `GaussianInt` exists to
witness (gap `t2.coefficient.1`).

**Why unverified.** The positive half is genuinely formalized: `Pauli/Matrix.lean`
grounds ℤ[i] in Mathlib's `GaussianInt` with `iU_sq`, `iU_pow_four`, `star_iU` and
twisted-product associativity over it, and `exact_ring.rs` / `sym_lean.rs` pin the
Rust exactly. The negative half — the actual reason the trait tower was re-cut — is
asserted only in prose at `exact.rs:21-23` and `lib.rs:71-78`;
`grep -rn 'Halvable|half' lean/PPVM` finds only `Projector.lean`'s
`oldStep_eq_half_iff`, which is about ℝ. Secondarily, the Rust ring is i64-backed:
`norm_sq` (`exact.rs:72`) and `Mul` wrap in release, so the exactness laws the
tests assert with zero tolerance hold only inside an unstated range, whereas the
Lean oracle is Mathlib's unbounded ℤ[i] (the tests use |re|,|im| <= 64 and one 2^40
case, never near the boundary).

**Round-2 amendment.** The positive claim is trivially true (2·re = 1 has no
solution in ℤ), so the Lean `no_half` theorem and its small-grid oracle stand. The
i64 half is **escalated to a live defect**: `norm_sq` and `magnitude` wrap at
|z| ≳ 2^31.5 ≈ 3.04e9, far below where the value stops being representable, so
(1+i)^64 = 2^32 — computed *exactly* by 64 `Mul`s — reports `magnitude() = 0.0` in
release, and z = 3037000500 wraps to −9223372036709301616 giving `magnitude() = NaN`;
either way `CoefficientThreshold::truncate`'s keep-rule
(`ppvm-pauli-sum-2/src/policy.rs:216`) **silently deletes** a coefficient of true
modulus ≈ 4.3e9 (debug panics at `exact.rs:73`). Two refinements: `Mul` itself wraps
at the same ~2^31.5 scale when both operands are large, so "faithful up to 2^63"
holds only for the representation and Add/Sub/Neg; and i128-with-saturating-readout
is the wrong repair (a saturated "exact" integer norm is a new lie) — compute
`magnitude` as `(re as f64).hypot(im as f64)`, which cannot overflow and satisfies
every existing assertion (`sym_lean.rs:259-282`). See
[adj U6](#u6-unenforced-preconditions).

**Proposed closure.** Add to `Pauli/Matrix.lean` (or a small `Exact.lean`):
`theorem no_half : ∀ z : GaussianInt, z + z ≠ ⟨1,1⟩` via parity of the real part,
i.e. `¬ ∃ z, 2*z = 1 + iU`, and state that ℤ[i] therefore satisfies the L4
capability set but not `Halvable`. Pin in `sym_lean.rs` with an exhaustive
small-grid assertion `∀ z in grid, z + z != GaussianInt::new(1,1)` naming that
theorem, plus either a documented range precondition or a checked-arithmetic
variant so the zero-tolerance ring laws are stated over the range where the i64
representation is faithful.

### G-055 — Lean's mulImpl omits the min_eps axis, including Sum x Sum's two cross-loop skips

- Class `fidelity` · Tier A · Severity high · Sector symbolic-coefficients · Status **open**
- Rust: `crates/ppvm-sym-2/src/mul.rs:235`
- Lean: `lean/PPVM/Instantiations/Symbolic.lean:882` `mulImpl` (and `:776` `epsClear_ne_retain_pointwise`)

**Claim.** `mulImpl k` is "the product ppvm-sym-2 actually implements"
(`Symbolic.lean:877`), so every truncation the shipped product performs is either
modelled by it or bounded by a companion theorem.

**Why unverified.** `mulImpl` filters on `sinDeg ≤ k` only; the coefficient axis
appears nowhere in it. The shipped `Sum × Sum` arm additionally skips an entire
cross loop when `s2.c0.abs() > self.min_eps` fails (`mul.rs:235`) and likewise for
`s1.c0` (`mul.rs:244`) — a whole-loop over-truncation of exactly the
`epsClear_ne_retain_pointwise` shape (a large stored coefficient can rescue a small
factor: min_eps=1e-6, s1={sin(x0) ↦ 1e3}, s2.c0=1e-7 drops a 1e-4 monomial that
`add_term`'s `|c| >= min_eps` rule keeps), and the guard's `>` disagrees with
`add_term`'s keep-rule at |c0| == min_eps exactly. `Sum::mul_term`'s per-term
min_eps and its `*self *= coeff` scalar arm (`mul.rs:112`, which applies no
truncation at all) are likewise unmodelled. `mul.rs`'s module doc justifies these
purely as "integration baseline, perf feature 6" — no Lean citation — while the
neighbouring `mul_term` clear got a full ℓ¹ treatment
(`epsClear_l1_eq`/`_lt`/`_error_lt`). This is the same code region as the confirmed
legacy bug, and the `-2` gate's only justification is that it matches the
hand-patched legacy.

**Proposed closure.** Extend `mulImpl` to `mulImpl (k : ℕ) (eps : ℝ)` with an
`epsRetained`-style filter on both the per-term rule and the two cross-loop skips,
then prove the analogue of `epsClear_l1_eq`/`_lt` for the cross-loop skip
(discarded mass = |c0|·ℓ¹(other table) < eps·ℓ¹) and a pointwise-inequality witness
like `epsClear_ne_retain_pointwise` for it. Pin in `sym_lean.rs`: build
s1={sin(x0) ↦ 1e3}, s2 with c0 = min_eps (boundary) and c0 = min_eps/10, and assert
the exact monomial set the guard drops versus what a per-monomial loop keeps.

### G-056 — Term::eval discards the i^k phase, but the Lean cites evalHom on a phase-free ring as its model

- Class `fidelity` · Tier A · Severity medium · Sector symbolic-coefficients · Status **open**
- Rust: `crates/ppvm-sym-2/src/eval.rs:138`
- Lean: `lean/PPVM/Instantiations/Symbolic.lean:260` `evalHom` (docstring: "`Term::eval` — evaluation at an angle vector θ")

**Claim.** `evalHom : SymRing →ₐ[ℝ] ℝ` models `Term::eval`, so the real read-out of
a propagated symbolic coefficient is the substitution homomorphism of the ring the
crate implements.

**Why unverified.** The ring the crate implements is `PhasedSymRing` (the file says
so itself at line 65: `Prod`'s phase byte is part of the hash key), yet `evalHom`
is defined on the phase-free `SymRing`, where the discrepancy cannot even be
stated. The shipped `Prod::eval` deliberately drops `i^k` (`eval.rs:131-136`), so
`Term::from(2.0).mul_i().eval(&[])` returns 2.0 for a coefficient whose value is 2i
— and this is the live read-out path: `Phase::apply` → `mul_i_pow` → `mul_phase`
phases coefficients on every operator product (X·Y = iZ), so a `trace`/`overlap`
fold that ends in `.eval()` reports a real number where the true value is
imaginary. The crate records this as `oldSuspectedBugs` #4 "preserved", and the
only pin is `sym_diff.rs:1194
divergence_real_eval_ignores_the_phase_on_both_crates` — agreement with legacy,
whose eval was `f64`-only. No Lean statement characterises the phase-forgetting map
or its defect.

**Proposed closure.** Define `forgetPhase : PhasedSymRing →ₐ[ℝ] SymRing` (the
augmentation along ℤ/4 → 1) and state `Term::eval = evalHom θ ∘ forgetPhase`, plus
the negative theorem `evalC θ x ≠ (evalHom θ (forgetPhase x) : ℂ)` for a phase-1
witness, i.e. the real read-out is not `Complex.re ∘ evalC`. Then add a
`sym_lean.rs` test asserting exactly that: for a phased `Term`, `eval()` equals the
phase-stripped value while `eval_complex()` equals `i^k` times it, and that
`eval() != eval_complex().re` — so the divergence is pinned to a Lean statement
rather than to old.

### G-057 — accumulate_assoc is credited as Term's additive-monoid law, but Term's `+` is non-associative under min_eps

- Class `strength` · Tier A · Severity medium · Sector symbolic-coefficients · Status **adjudicated-spec** (round 2 — [adj U6](#u6-unenforced-preconditions))
- Rust: `crates/ppvm-sym-2/src/term.rs:323`
- Lean: `lean/PPVM/Algebra/GradedMap.lean:95` `accumulate_assoc` (and `:92` `accumulate_comm`)

**Claim.** `Term` is an additive monoid, as machine-checked by
`GradedMap.accumulate_comm` / `accumulate_assoc` — the law `add.rs:68` invokes to
force divergence #1 ("x + c == x for c != 0 is indefensible") and `coeff.rs:30`
invokes for the value domain generally.

**Why unverified.** `accumulate_assoc` is `add_assoc` on an untruncated `CMap`; the
shipped `+` truncates. `Sum::add_const` (`term.rs:323`) drops |c| < min_eps while
the `Sum + Sum` arm (`add.rs:142`) adds c0 unconditionally, so with a = Sum-form
and min_eps = 1e-3: `(a + Term::from(6e-4)) + Term::from(6e-4)` drops both addends
(result c0 = 1.0) whereas `a + (Term::from(6e-4) + Term::from(6e-4))` folds the
constants first and keeps 1.2e-3 — the very `x + c == x` shape divergence #1 calls
indefensible, present by design at a different arm. The oracle test that claims to
pin these laws (`sym_lean.rs:634`) builds terms only at the default
min_eps = f64::EPSILON with coefficients in [-2,2], so the truncating branch never
fires and the test would survive removing `add_const`'s threshold entirely; the
truncated behaviour is pinned only diff-vs-legacy (`sym_diff.rs:595`).

**Round-2 amendment.** Adjudicated: the code is right and the *citation* is what
is wrong. Exact associativity is unattainable for any nonzero insert-time
thresholding rule (that is Lean's own `eps_drop_at_insert_ne_drop_at_end`), so the
non-associativity is an accepted consequence of thresholding, bounded by
(#dropped constants)·min_eps in c0; what is defective is that `add.rs:66-70` cites
`accumulate_assoc` to forbid `x + c == x` while the same file ships exactly that in
the Sum+Const arm, and that `+` is a function of which of Const/One/Sum represents
a value rather than of the value itself. Measured witnesses: at min_eps = 1e−3,
(a+6e−4)+6e−4 → c0 1.0 vs a+(6e−4+6e−4) → c0 1.0012; a+b → 1.0 vs b+a → 1.0006.
See [adj U6](#u6-unenforced-preconditions).

**Proposed closure — one clause INVALIDATED by round 2.** The alternative "an
explicit scope note that `accumulate_assoc` holds only at eps = 0" is **false and
must not be shipped**: at min_eps = 0.0 with a = 1 + sin(x0) and b = c = 1e−16,
(a+b)+c has c0 = 1.0 while a+(b+c) has c0 = 1.0000000000000002 — f64 addition is
itself non-associative, so `Term`'s `+` is non-associative at *every* min_eps, and
that scope note would install a new false law where the old one was. Replacement:
state in Lean the negative theorem for the truncated accumulation (`accumulateEps
eps` is not associative — explicit 6e−4 + 6e−4 witness at eps = 1e−3, mirroring
`eps_drop_at_insert_ne_drop_at_end`); say in `add.rs` that `accumulate_assoc` is a
law of an **exact** coefficient ring and is not a law of an f64-backed `Term` at any
eps, with the defect bounded by (#dropped constants)·min_eps plus f64 rounding; and
add the `sym_lean.rs` test pinning the witness and the bound.

### G-058 — The tableau's own g-rule row product is pinned to no Lean phase model

- Class `bridge` · Tier A · Severity medium · Sector tableau-and-symbolic (skeptic) · Status **adjudicated-spec** (round 2 — [adj U5](#u5-measurement-sign))
- Rust: `crates/ppvm-tableau-2/src/data.rs:318`
- Lean: `lean/PPVM/Pauli/Word.lean` `phaseExpN` (and `Pauli/Phase.lean` `phaseExp_eq_ref`, `Pauli/Matrix.lean` `pauliMat_mul`)

**Claim.** `Row::mul_assign`'s bit XOR and its ℤ/4 phase
`(2·popcount(sign) + popcount(imag)) % 4 + rhs.phase` compute the Pauli-group
product, i.e. they equal `mulWord` and `phaseExpN` site-by-site.

**Why unverified.** `ppvm-tableau-2` carries its own copy of the
Aaronson–Gottesman g-rule (`data.rs:299-320`; the sign/imag sum-of-products at
311-312, the fold at 318-319). It is the primitive under `row_multiply`
(`clifford.rs:256-261`), `get_deterministic_outcome` (`data.rs:485-489`),
`update_tableau_according_to_outcome` (`data.rs:536/540`) and
`compute_decomposition`, so every measurement sign flows through it. The
byte-identical formula in `ppvm-pauli-word-2` (`product.rs:64-73`) IS pinned to
Lean — `pauli_word_lean.rs:134/153` assert `phaseExp` against `pauliMat_mul` and
175-260 assert `phaseExpN_cocycle`/`_self`/`_sub_comm` — but no oracle pins the
tableau's copy: grepping `tableau_lean.rs` and `tableau_diff.rs` for
`mul_assign`/`phaseExp` returns nothing, `tableau_lean.rs`'s row tests cover only
per-gate conjugation, and `assert_symplectic_frame` (`tableau_lean.rs:490`) is
phase-blind. Its only stated provenance is the docstring's port citation to the
LEGACY path `ppvm-pauli-word/src/phase/mul.rs` (`data.rs:293`), and its only
behavioural check is `tableau_diff.rs`'s agreement with that legacy. A
`2*sign_count` → `sign_count` mutation would leave every `*_lean.rs` assertion
green (`tableau_lean.rs:923-927` checks the decomposition phase only mod 2), which
is the same blind spot G-021 describes one level up.

**Round-2 amendment.** The g-rule is **correct**: `Row::mul_assign` equals
g = ab + cd + 2bc − (a⊕c)(b⊕d) mod 4 plus `+= rhs.phase`, verified against genuine
2×2 complex matrix products at all 16 per-site patterns (0 mismatches), across word
boundaries at n = 70 (480 checks, 0 mismatches), and with the orientation pinned as
`self·rhs` (4,400 row-product checks: 0 failures as dst·src, 1,038 as src·dst). Two
corrections to the closure: no `pub(crate)`/`#[cfg(test)]` **hook is needed** —
`StabilizerFrame::row_multiply` is a public trait method and `Tableau::rows()` is
`pub`, and an exhaustive external matrix-grounded oracle was driven from
`crates/ppvm-conformance-2/tests/` with zero source changes; and the oracle should
assert the orientation as well as the magnitude, since orientation is the half a
magnitude-only test misses. See [adj U5](#u5-measurement-sign).

**Proposed closure.** Add a `tableau_lean.rs` test that transcribes
`mulWord`/`phaseExpN` from `Word.lean` and compares them against the real `Row`
product — exhaustively for n = 1, 2 over all (phase, x, z) pairs and seeded-random
for n = 6, 70 — asserting full ℤ/4 equality of the phase and equality of both bit
planes, naming `phaseExp_eq_ref` / `phaseExpN_cocycle` / `pauliMat_mul` as
`pauli_word_lean.rs` does for `ppvm-pauli-word-2`. This needs a
`pub(crate)`→test-visible hook (or a `#[cfg(test)]`/`pub` row-multiply accessor on
`Tableau`, e.g. via the existing `rows()` snapshot plus `row_multiply`), and it
retires the legacy citation at `data.rs:293` in favour of a Lean one.

### G-059 — No oracle pins the ℤ/4 sign of s_dag, sqrt_x(_dag), sqrt_y(_dag) or cy on any tableau row

- Class `bridge` · Tier A · Severity high · Sector tableau-core · Status **adjudicated-spec** (round 2 — [adj U5](#u5-measurement-sign))
- Rust: `crates/ppvm-tableau-2/src/clifford.rs:491`
- Lean: `lean/PPVM/Tableau/Batch.lean:216` `isSitewise_extSqrtY` (and `206`/`211`/`221`, `183`)

**Claim.** Each of the tableau's six extension gates rewrites every row's bits and
ℤ/4 phase by the audited conjugation table.

**Why unverified.** `tableau_lean.rs`'s row-level fidelity tests cover exactly five
single-qubit gates (h, s, x, y, z — lines 389-399) and two two-qubit gates (cnot,
cz — lines 439-446). The only test that touches
s_dag/sqrt_x/sqrt_x_dag/sqrt_y/sqrt_y_dag/cy/zcy is
`every_gate_preserves_the_symplectic_frame` (line 526), whose predicate is
`omega()` (line 480) computed from the x/z planes only — it never reads `row.2`
(the phase). So mutating `clifford.rs:507` from `xw & !zw` to `zw & !xw` (i.e.
shipping √Y† where the caller asked for √Y) leaves every assertion in every
`*_lean.rs` file green; the only thing that would fail is the legacy differential
suite (`tableau_diff.rs`) and `ppvm-tableau-2/tests/behaviour.rs`'s batch-vs-loop
macro, both of which are agreement-with-legacy evidence, not verification. Note
that the Lean names for these gates are adjoint-labelled relative to the tableau's
forward convention (see the appendix entry on `ext*` naming), so a systematic
adjoint flip in the extension set is exactly the mutation this sector cannot
currently detect.

**Round-2 amendment.** The six extension gates are **correct** (√P = e^{−iπP/4},
CY = ctrl-Y conjugation), pinned by 625,856 exact full-Pauli-set expectation
comparisons over 400 random 14-gate circuits, worst |Δ| 1.110e−15 — so the
adjoint-flip mutation this row describes is a real blind spot in the *oracle*, not a
shipped bug. Load-bearing correction for whoever writes the oracle: the tableau
implements the **forward** unitary, i.e. rows conjugate as P ↦ U P U†, which is the
**transpose** of the `CliffordExtensions` doc table at
`ppvm-traits-2/src/gates.rs:95-107` (documented as U†PU). Transcribing the reference
predicates from that table without transposing will fail on correct code. See
[adj U5](#u5-measurement-sign).

**Proposed closure.** Extend the (name, gate, reference) table at
`tableau_lean.rs:389` with s_dag/sqrt_x/sqrt_x_dag/sqrt_y/sqrt_y_dag transcribed
from the `isSitewise_*` predicates, and the two-qubit table at line 439 with cy
transcribed from `conjCY_bits`/`conjCY_sign`; the existing loop already asserts
bits AND phase per row over 10 seeded n=6 frames, so only the reference functions
are missing.

### G-060 — cz_block_pairs/cz_block do not enforce the disjoint-support precondition Lean proves is necessary

- Class `bridge` · Tier B · Severity medium · Sector tableau-core · Status **adjudicated-defect** (round 2 — [adj U6](#u6-unenforced-preconditions))
- Rust: `crates/ppvm-tableau-2/src/data.rs:1474`
- Lean: `lean/PPVM/Tableau/Batch.lean:260` `czSeq_phase` and `:291` `czSeq_phase_needs_disjoint`

**Claim.** The fused `count_ones() & 1` phase of a CZ block equals the sequential
per-pair CZ loop.

**Why unverified.** `Batch.lean` proves `czSeq_phase` only under
`P.Pairwise Disjoint2`, and `czSeq_phase_needs_disjoint` exhibits a concrete
counterexample (pairs (0,1),(1,2) on x=111,z=000: sequential flips the phase, the
batched parity counts zero). `data.rs:1462-1468` acknowledges this, but the only
guard in the code is a `debug_assert_eq!` that all bits share a word
(`data.rs:1483`); nothing checks `offset >= count`, and `cz_block`
(`data.rs:1621`) derives `offset = hi - lo` from caller-supplied bases without
checking it either. Both are `pub`. Failure scenario: `cz_block(0, 2, 5)` produces
pairs (0,2),(1,3),(2,4),(3,5),(4,6) whose supports overlap; the bits still come out
right (CZ only writes z and only reads x) but the ℤ/4 sign silently differs from
the per-pair `cz` loop. No test exercises overlap: every case in
`data.rs:1721-1788`, `tableau_diff.rs:373-382`
((0,17,17),(0,32,17),(17,34,17),(34,51,17),(51,68,17),(60,70,8),(0,64,20)) and
`tableau_lean.rs:557` (0,20,10) has offset ≥ count, so the existing
batched-vs-per-pair expansion assertions can never fire.

**Round-2 amendment — the "bits still come out right" claim in this row is
FALSE.** `z_delta = ((x>>offset)&mask_c) | ((x<<offset)&mask_t)` uses **OR**, so a z
bit written by two overlapping pairs is set once instead of XOR-cancelled: measured,
`cz_block_pairs(0,1,2)` on x = 111, z = 000 gives x = 111, z = 111, ph = 0
(+Y₀Y₁Y₂) where the per-pair loop gives x = 111, z = 101, ph = 2 (−Y₀X₁Y₂) — wrong
in both bit planes *and* in the sign. `cz_block(0,1,2)` reproduces it, and the
ledger's own `cz_block_pairs(0,2,5)` witness diverges on 158 of 200 random 10-qubit
states. Also: `offset >= count` is exactly the right predicate and **subsumes** the
degenerate offset == 0 case (pairs collide iff i − j = offset or offset = 0), which
currently XORs z ^= x over the block where the scalar `cz(q,q)` is a probed no-op.
Escalated to a live defect: `cz_block(0,1,n)` — CZ on every adjacent pair, the most
natural brickwork call — is precisely the broken case. No existing caller violates
the precondition, so the `debug_assert` breaks nothing; the release-mode per-pair
fallback is behaviour-changing (wrong → right) and needs sign-off. See
[adj U6](#u6-unenforced-preconditions).

**Proposed closure.** Add
`debug_assert!(offset >= count, "cz_block pairs must have disjoint supports (Batch.lean::czSeq_phase_needs_disjoint)")`
to `cz_block_pairs` and the equivalent `hi - lo >= count` check to `cz_block` (or
fall back to the per-pair loop), document the precondition in the signature, and
add an oracle case to `tableau_lean.rs` that (i) asserts the debug assert fires /
the fallback engages for offset < count and (ii) replays
`czSeq_phase_needs_disjoint`'s exact witness (x=111, z=000, pairs (0,1),(1,2))
against the per-pair loop.

### G-061 — build_masks collapses duplicate indices, silently violating seqApply_eq_batchApply's Nodup hypothesis

- Class `bridge` · Tier B · Severity medium · Sector tableau-core · Status **adjudicated-defect** (round 2 — [adj U6](#u6-unenforced-preconditions))
- Rust: `crates/ppvm-tableau-2/src/clifford.rs:594`
- Lean: `lean/PPVM/Tableau/Batch.lean:91` `seqApply_eq_batchApply`

**Claim.** Every `*_many` fused sweep equals the corresponding per-index gate loop,
which is the documented contract of `CliffordBatch` (its default bodies literally
loop, `ppvm-traits-2/src/gates.rs:139-176`).

**Why unverified.** `seqApply_eq_batchApply` carries `L.Nodup`, and
`Batch.lean:27-34` stresses that the independence side condition "belongs in the
record rather than in an unwritten assumption of a shipped public API". But
`build_masks` ORs each index into a per-word mask (`clifford.rs:595`), so a
repeated index contributes one bit, not two: `x_many(&[0,0])` flips the row phase
once whereas `x(0); x(0)` flips it twice (identity), and `h_many(&[0,0])` swaps the
x/z bits once whereas `h(0); h(0)` is the identity. The tableau's override
therefore disagrees with the trait's own default implementation on duplicate input,
with no assert and no doc note; `clifford.rs:581` mentions the Nodup hypothesis
only inside the `build_masks` doc-comment prose. The batch-vs-loop oracles
(`tableau_diff.rs:265-320`, `ppvm-tableau-2/tests/behaviour.rs:581-590`) all pass
strictly increasing index vectors, so the divergence is untested.

**Round-2 amendment.** Confirmed and escalated to a **live defect**: the correct
answer for `g_many([q,q])` is conjugation by G² (the identity for x/y/z/h, Z for s,
X for √X, Y for √Y — i.e. exactly what the per-index loop computes), the fused
sweep applies G once, and this is reachable from a legal `.stim` file today.
Measured through the public `ppvm_stim::run_string`: `"X 0 0\nM 0"` →
`[Some(true)]` while `"X 0\nX 0\nM 0"` → `[Some(false)]`, where the truth is
`false` — a flipped sampled bit, on **both** the `traits-2` and the default
(legacy) backends, since legacy `ppvm-tableau/src/gates/clifford.rs:365` ORs
identically. `cnot_many`/`cz_many`/`cy_many` are unaffected (they re-read bits per
iteration). See [adj U6](#u6-unenforced-preconditions).

**Proposed closure — amended by round 2; two variants ruled out.** (a) Signalling
duplicates out of `build_masks` via its `Option` is **unsafe**: `None` means
"nothing to do" and every caller does `let Some(..) = .. else { return }`, so that
would silently turn the whole gate into a no-op. (b) Building the mask with XOR
instead of OR is **wrong for four gate families**: `build_masks` is shared by
x/y/z/h/s/s_dag/sqrt_x(_dag)/sqrt_y(_dag), and XOR would make `s_many(&[q,q])` the
identity when the truth is Z-conjugation (S² = Z ≠ I). The complete fix is to detect
a duplicate and **fall back to the per-index `Clifford` loop** (correct by the trait
contract) for every gate, in `build_masks` and in the inline mask of
`y_many_skipping` (`clifford.rs:602`), plus expand-or-apply-per-target in
`ppvm-stim/src/executor/gates.rs` — note that *rejecting* `X 0 0` in
`validate.rs` would break valid Stim input. Document the precondition on
`CliffordBatch`. Oracle: a `tableau_lean.rs` case citing `seqApply_eq_batchApply`
comparing `g_many(&[q, q])` against `g(q); g(q)` for each of
x/y/z/h/s/s_dag/sqrt_* — asserting **equality with the loop** (or a panic), never
the shipped single-application semantics — plus an absolute anchor that
`X 0 0` on `|0…0⟩` measures 0. Behaviour-changing: needs sign-off, and it makes
`-2` diverge from legacy unless legacy is fixed in the same commit.

### G-062 — The frame's ℤ/4 sign column has no Lean invariant at all (row phases never proven real)

- Class `coverage` · Tier A · Severity high · Sector tableau-core · Status **adjudicated-spec** (round 2 — [adj U5](#u5-measurement-sign))
- Rust: `crates/ppvm-tableau-2/src/data.rs:492`
- Lean: `lean/PPVM/Tableau/Frame.lean:53` `IsSymplecticFrame` (none for the phase column)

**Claim.** Every row of a live tableau carries a real phase (ℤ/4 ∈ {0,2}), i.e. the
generators stay Hermitian, and that is preserved by every gate, by
`Row::mul_assign` on commuting rows, and by the measurement projection.

**Why unverified.** `Frame n` (`Frame.lean:43`) is a pair of maps into the
phase-stripped `Sp n`, and `Frame.lean:240` says outright "The outcome bit b is a
phase, and Sp n is the phase-stripped space, so b does not appear below". So no
Lean statement constrains the phase column. Yet the crate depends on its realness
in three places: `get_deterministic_outcome`'s
`debug_assert!(result.phase == 0 || result.phase == 2, "Measurement result cannot be imaginary!")`
(`data.rs:492`); `odd_phase_destabilizer_mask` (`data.rs:1134`) plus the mask term
in `compute_phase_with_mask_static` (`data.rs:146`), which `tableau_lean.rs:954`
asserts is *always empty* with the justification, in a comment, that "the rows of a
valid frame are Hermitian, so their ℤ/4 phases are even ... the mask term is
computed because old computes it (behaviour preservation)" — a load-bearing
invariant stated only as a test comment plus a runtime assert; and
`BranchPhase.lean`'s `FrameInvolution` hypothesis (M²=I), which every case-a/case-b
theorem in `Projection.lean` is stated under and which is checked empirically
(`tableau_lean.rs:914`) but never derived from the frame invariant. The missing
prerequisite is also absent from `Word.lean`: nothing says `phaseExpN p q` is even
when ω(p,q)=0, which is what makes `mul_assign` phase-real-preserving on the
commuting pairs the projection and `get_deterministic_outcome` multiply.

**Round-2 amendment.** The invariant is **true and provable in three lines**, so it
is a theorem rather than a coincidence: every one of the 29 phase writes in
`clifford.rs` has the form `phase ^= (predicate) << 1` (delta 0 or 2), and frame-row
`mul_assign` is only ever applied to commuting pairs — stabilizer × stabilizer in the
projection, and the pivot g_q into d_i only for i ≠ q_idx where ω(g_q, d_i) = 0 —
and for commuting Hermitian P, Q we have (PQ)† = QP = PQ. Hence phase ∈ {0,2},
`odd_phase_destabilizer_mask()` is identically 0 and the mask term is vacuous;
measured 0 odd-phase rows over every row of 1,200 + 48 circuits and 63,000+ row
inspections, with the runtime `debug_assert` never firing across 25,077 deterministic
measurements. The closure stands; item (iii)'s `assert!(phase % 2 == 0)` inside
`assert_symplectic_frame` is the cheap first step. See
[adj U5](#u5-measurement-sign).

**Proposed closure.** Extend `Frame` (or add `PhasedFrame`) with a ℤ/4 phase per
generator, define
`IsRealFrame T := ∀ i, IsRealPhase (T.destabPhase i) ∧ IsRealPhase (T.stabPhase i)`,
prove (i) `phaseExpN p q` is even when `omegaN p q = 0` (from `phaseExpN_sub_comm`),
(ii) each gate's delta is `if b then 2 else 0` so `IsRealFrame` is gate-invariant,
(iii) `IsRealFrame` is preserved by `rowUpdate`/`projectFrame` given the frame
relations, and (iv) `frameInvolution_zero_iff`'s hypothesis follows. Then have
`tableau_lean.rs` assert `phase % 2 == 0` on all 2n rows inside
`assert_symplectic_frame` (currently it only checks omega), which makes the
mask-is-empty assertion a consequence rather than a coincidence.

### G-063 — The tableau's `r` and `u3` have no Lean model; the only RotXY theorem is for the opposite order

- Class `coverage` · Tier A · Severity medium · Sector tableau-core · Status **open**
- Rust: `crates/ppvm-tableau-2/src/gates.rs:82`
- Lean: `lean/PPVM/Instantiations/Rotation.lean:654` `rotXY_heisenberg_order` (models the reverse order)

**Claim.** `GeneralizedTableau::r(q, φ, θ)` is the rotation by θ about the in-plane
axis cos φ·X + sin φ·Y (so r(q, π/2, θ) = ry(q, θ)), and
`u3(θ,φ,λ) = RZ(φ)·RY(θ)·RZ(λ)`.

**Why unverified.** `gates.rs:82` emits `rz(−φ); rx(θ); rz(φ)` — deliberately the
REVERSE of `ppvm-pauli-sum-2`'s order, per the docstring at `gates.rs:77-81` ("the
tableau runs in the Schrödinger picture, so the sub-rotations are applied in
forward order"). The only Lean statement about this family,
`rotXY_heisenberg_order`, proves
`mz(−φ)∘mx(θ)∘mz(φ) = rotAxis(cosφ,sinφ,0) θ` for the pauli-sum's order and
`Rotation.lean:585-590` explicitly flags that the forward order "composes to the
inverse rotation" and that a wrong-order implementation "passes every other
rotation test". Nothing states or proves the tableau's forward order is the correct
one for a state-space (as opposed to observable-space) simulator, and no Lean model
of `u3` exists at all. The only tests are tautological:
`ppvm-tableau-2/tests/behaviour.rs:770 r_matches_rz_rx_rz` re-lists the impl's own
three calls in the impl's own order, and its assertion is a measurement outcome on
a single qubit prepared from |0⟩ — where ry(θ) and ry(−θ) give identical outcome
distributions — so it cannot detect the very axis-sign error contract 10 names.

**Proposed closure.** State the tableau's order in Lean:
`theorem rotXY_schroedinger_order (φ θ) (v) : mz φ (mx θ (mz (-φ) v)) = rotAxis (cos φ, sin φ, 0) (-θ) v`
(or the state-space analogue with the sign that makes r(π/2,θ)=ry(θ)), plus a `u3`
composition theorem on the same Vec3 model. Oracle: add a `tableau_lean.rs` test
that `r(q, π/2, θ)` and `ry(q, θ)` agree on the full amplitude vector (indices and
coefficients), and `r(q, 0, θ)` vs `rx(q, θ)`, on a non-basis state (e.g. after
`h; t`) where the two signs are distinguishable — and the same for `u3` vs its
RZ·RY·RZ product with independently varied φ, λ.

### G-064 — The T-gate branch coefficients are pinned to nothing but the legacy constants

- Class `coverage` · Tier A · Severity low · Sector tableau-core · Status **open**
- Rust: `crates/ppvm-tableau-2/src/gates.rs:28`
- Lean: none

**Claim.** `t()`/`t_dag()` branch the amplitude vector by the coefficient pair of
T = diag(1, e^{iπ/4}), i.e. (cos π/8, −i sin π/8) up to an unobservable global
phase e^{iπ/8}.

**Why unverified.** The two hardcoded `Complex64` literals (`gates.rs:28-36`) are
the only specification of the T gate; no Lean statement mentions π/8, T, or the
pair. `gates.rs:14-18` credits `Rotation.lean`'s `rot_norm_sq` / `rot_rot`, but
those are theorems about a real (cos θ, sin θ) plane and constrain neither the
relative −i between the two branches nor the angle π/8 — any unit-norm pair
satisfies them. The only tests are `ppvm-tableau-2/tests/behaviour.rs:789` (t then
t_dag is the identity, which holds for any (c,s) with |c|²+|s|²=1 and a conjugate
partner) and the legacy differential suites; `tableau_lean.rs` uses `t()` only as a
state generator (lines 584, 888). So a wrong angle (e.g. π/16 constants) is caught
only by agreement with legacy.

**Proposed closure.** Add a Lean lemma that the crate's pair equals
`(Complex.exp (I*π/8) * cos (π/8), -I * Complex.exp (I*π/8) * sin (π/8))` and that
this is `exp(I*π/8) • Rz(π/4)`'s branch pair, i.e. `t = rz(π/4)` up to global phase
(with the mixture being classical, the global phase is provably unobservable).
Oracle: assert in `tableau_lean.rs` that `t(q)` and `rz(q, π/4)` produce identical
amplitude support and coefficients up to the single scalar e^{iπ/8}, and that `t`
applied 8 times is the identity on a branched state.

### G-065 — projectFrame/canonicalize proven only on bits — nothing says the represented state is preserved

- Class `strength` · Tier A · Severity high · Sector tableau-core · Status **open**
- Rust: `crates/ppvm-tableau-2/src/data.rs:548`
- Lean: `lean/PPVM/Tableau/Frame.lean:296` `isSymplecticFrame_projectFrame`

**Claim.** `update_tableau_according_to_outcome` maps the tableau representing |ψ⟩
to one representing the projected state (I+(−1)^b Z_q)|ψ⟩, hence
`StabilizerFrame::canonicalize` (`clifford.rs:284`) genuinely has nothing to do.

**Why unverified.** `isSymplecticFrame_projectFrame` is credited
(`clifford.rs:263-279` and `data.rs:499-508`) with justifying the `canonicalize`
no-op and "every downstream `compute_decomposition`", but it proves only that nine
ω-pairings still hold on the phase-stripped space. The Rust does strictly more: it
multiplies `g_q` into rows with the g-rule (`data.rs:535/538`, so signs move),
installs `destabilizers[q_idx] = g_q` and sets
`stab_q.phase = if outcome { 2 } else { 0 }` (`data.rs:548`) — the sign that
decides WHICH of the two ±Z_q eigenspaces is represented. Nothing in the Lean
development defines the state (stabilizer group / +1 eigenspace / density operator)
a frame represents, so "the projection restores the pairing" is proven in the weak
bit sense while the load-bearing sign assignment is unproven; `Projection.lean`'s
own scope note (lines 55-66) concedes the missing link and says "The old and the -2
crate agree here verbatim, so it is a specification gap", i.e. the sign leg of the
projection is currently justified by legacy agreement. The oracle is
correspondingly bit-only: `measurement_preserves_the_symplectic_frame`
(`tableau_lean.rs:568`) uses the phase-blind `omega()`.

**Proposed closure.** Add a Hilbert-space (or at least group-theoretic) model:
define `stabGroup T` from the phased rows and
`Represents T ψ := ∀ g ∈ stabGroup T, g • ψ = ψ`, then prove (a) each gate's row
update satisfies `Represents (gate T) (U ψ)`, and (b)
`Represents (projectFrame T q p, sign b) (P_b ψ / ‖·‖)` — at minimum for n=1,2 by
`decide` over the finite phased group with an explicit ℂ² / ℂ⁴ model, generalizing
the `Frame.lean` argument. Oracle: extend `assert_symplectic_frame` to also verify,
on the real engine, that every stabilizer row's signed Pauli fixes the
reconstructed state vector (feasible for n ≤ 8 by dense simulation), and that a
repeated measurement reproduces its own sign.

### G-066 — conjH/conjS/conjCNOT/conjCZ are fixed by asserted tables, not derived from matrices

- Class `strength` · Tier A · Severity medium · Sector tableau-core · Status **open**
- Rust: `crates/ppvm-tableau-2/src/clifford.rs:358`
- Lean: `lean/PPVM/Pauli/Conjugation.lean:85` `conjH_sign` (and `conjH_X/_Z/_Y:103-105`, `conjCNOT_Xc..:697-711`, `conjCZ_Xc..:713-724`)

**Claim.** The tableau's per-gate sign deltas ARE the signs of genuine conjugation
G·P·G† for G ∈ {H, S, CNOT, CZ}.

**Why unverified.** Every sign in the tableau reduces to these four maps, but each
is *defined* by its bit-and-sign rule (`Conjugation.lean:48/53/577/595`) and the
theorems credited with fixing it — `conjH_X`, `conjH_Y`, `conjCNOT_Xc`,
`conjCNOT_YcYt`, `conjCZ_Xt`, … — are `by decide` on that same definition, i.e.
tautologies restating the definition; the identification with G·P·G† comes from the
docstring calling them "the standard CNOT/CZ tableau tables" (line 690) and from
the hom-agrees-on-generators argument, whose premise (that true conjugation takes
those generator values) is itself the asserted table. Contrast conjX/conjY/conjZ
(lines 187-201) and `conjPauliZ_eq_conj` (line 271), which ARE derived as
conjugation inside the group, and `Pauli/Matrix.lean`, which grounds the product
phase in real ℤ[i] matrices. H and S are not Pauli-group elements and H is not
expressible over `GaussianInt` (it needs 1/√2), so the matrix leg that makes
`phaseExp` trustworthy is structurally absent for the conjugation layer —
`lean/README.md:126` already concedes "the operator conjugation is not derived" for
the rotation branch, and the same holds here. Effect: a sign convention error in
`conjH_sign` or `conjCNOT_sign` would be reproduced identically by the Lean model
and by the oracle (`tableau_lean.rs` transcribes the same predicate), and only the
legacy diff / stim conformance would notice.

**Proposed closure.** Extend `Pauli/Matrix.lean` over ℤ[i, 1/√2] (or ℚ(√2)[i], or
just ℂ with `Matrix (Fin 2) (Fin 2) ℂ` and `norm_num`) with explicit H, S, CNOT, CZ
matrices and prove
`hMat * pauliMat a b * hMat⁻¹ = iU^(conjH ⟨0,a,b⟩).phase • pauliMat (conjH …).x (conjH …).z`
for all 4 bit patterns (16 for the two-qubit gates) — turning
`conjH_sign`/`conjS_sign`/`conjCNOT_sign`/`conjCZ_sign` from definitions into
corollaries. Oracle: add a dense 2x2 / 4x4 complex matrix check in
`tableau_lean.rs` that the tableau's row update on a single-row frame agrees with
the literal matrix conjugation of the corresponding Pauli.

### G-067 — truncation_l1_bound oracle is a float tautology; it survives any mutation of Policy::truncate

- Class `bridge` · Tier A · Severity high · Sector truncation-policy-loss · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/policy.rs:168`
- Lean: `lean/PPVM/Algebra/Truncation.lean:44` `l1_bound`

**Claim.** The real Rust truncation incurs at most the ℓ¹ mass of the terms it
dropped (`pauli_sum_lean.rs:18` lists this as a reproduced law).

**Why unverified.** `pauli_sum_lean.rs:240-292` computes `dropped` from the Rust,
then `error = Σ c*e` and `l1 = Σ|c|` over that same list with
`e = sin(..).clamp(-1,1)`, and asserts `error.abs() <= l1 + TOL`. Since |e| ≤ 1 by
construction, that inequality is arithmetically true for *any* list — empty,
complete, or arbitrary. Mutating `MaxPauliWeight::truncate` to a no-op (dropped
empty → 0 ≤ 0) or to `map.retain(|_,_| false)` (drops all → the inequality still
holds) both pass the l1 assertion. The only Rust-sensitive assertion in the test is
the trailing `w.weight() <= cap` loop, which is one-sided: it cannot detect
over-dropping, and it does not pin the `<=`-vs-`<` boundary (`weight == cap` KEPT)
that `policy.rs`'s keep-rule turns on — that boundary is pinned only by
`pauli_sum_truncation_boundary_diff.rs`, i.e. by agreement with legacy. So the one
Tier-A truncation bound has, in effect, no bridge.

**Proposed closure.** Make the oracle pin the *predicate*: assert set equality
`kept_keys == {k ∈ full : k.weight() <= cap}` (and for `CoefficientThreshold`
`== {k : |c| >= t}`), including deliberately-placed `weight == cap` and `|c| == t`
boundary terms, and assert the error against a policy-parameter bound (see G-074)
rather than against a quantity recomputed from the same dropped list.

### G-068 — cutoff_mismatch's tableau side is never executed; the oracle assertion is a float identity

- Class `bridge` · Tier B · Severity medium · Sector truncation-policy-loss · Status **open**
- Rust: `crates/ppvm-tableau-2/src/data.rs:1240`
- Lean: `lean/PPVM/Algebra/Truncation.lean:237` `cutoff_mismatch`

**Claim.** `policy.rs:174-177` and `tableau-2/src/lib.rs:17,69` — the two backends
disagree at |c| == threshold (PauliSum keeps on `>=`, the tableau drops on strict
`>`), "machine-checked" by `cutoff_mismatch`.

**Why unverified.** `cutoff_mismatch` is `(t ≤ |t|) ≠ (t < |t|)` for `0 ≤ t` — a
statement about two ℝ relations that mentions neither backend's keep predicate; it
cannot detect a change to either. The oracle side is worse:
`pauli_sum_lean.rs:328-333` closes the mismatch test with
`assert_ne!(t <= t.abs(), t < t.abs())`, recomputed from a local f64 and true
unconditionally for t ≥ 0, and no tableau object is ever constructed.
`tableau_lean.rs` contains no occurrence of cutoff/threshold/trim at all, so the
strict-`>` rule at `data.rs:1240` / `:749` / `mixture/data.rs:145` is pinned only by
`tableau_behaviour_diff.rs`, i.e. by agreement with legacy. The PauliSum half
(`pauli_sum_lean.rs:317-326`, terms at exactly/above/below threshold at n=4) is
genuine and is the one good part.

**Proposed closure.** Restate the Lean side as two named predicates
(`keepSum t c := t ≤ N c`, `keepTab t c := t < N c`) with
`keepSum t t = true ∧ keepTab t t = false`, so mutating either spelling breaks the
proof. Then add a `tableau_lean.rs` test that seeds an amplitude with
`norm_sqr() == cutoff_sq` exactly (representable: cutoff = 0.5, amplitude 0.5)
through `trim`/`rotate_2` and asserts it is DROPPED, beside the existing PauliSum
assertion that the twin is KEPT.

### G-069 — CombinedPolicy's conjunction/commutation and the usize::MAX sentinel skip have no `*_lean.rs` oracle

- Class `bridge` · Tier B · Severity medium · Sector truncation-policy-loss · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/policy.rs:258`
- Lean: `lean/PPVM/Algebra/GradedMap.lean:600` `retain_seq_eq_retain_and` (also `retain_comm:612`, `retain_of_all_true:622`, `retain_weight_le_eq_self:639`)

**Claim.** `policy.rs:228-232` and `:158-164` — the two sequential retain passes
compute the conjunction of the keep-rules and the surviving set is independent of
pass order ("license a future backend to fuse or reorder the two walks"); and
`MaxPauliWeight`/`MaxLossWeight`'s `usize::MAX` early return is observationally
exact, zero-coefficient terms included.

**Why unverified.** These are the correct theorems for the code, but no `*_lean.rs`
pins the Rust to them. `column_store_lean.rs:165` does construct
`CombinedPolicy(CoefficientThreshold{0.5}, MaxPauliWeight(n))` and call
`truncate()`, yet its only assertion is `assert_aligned` (column/index
consistency) — it never checks *which* terms survived, so it passes whether the
composite computes the conjunction, the disjunction, or only P1.
Order-independence, the fused-order license, and the sentinel-skip exactness on
zero-coefficient terms are pinned only by
`pauli_sum_truncation_boundary_diff.rs` (old `CombinedStrategy` vs new
`CombinedPolicy`) and `pauli_sum_loss_diff.rs` (`MaxLossWeight(1)`), i.e. by
agreement with legacy — which by the audit's own rule is evidence of a faithful
port and nothing about correctness. `MaxLossWeight` in particular has zero
Lean-oracle coverage of any kind.

**Proposed closure.** Add to `pauli_sum_lean.rs`: (i)
`combined_policy_is_conjunction` — for random supports assert
kept == {k : |c| >= t} ∩ {k : weight(k) <= w}, and that
`CombinedPolicy(P1,P2).truncate` and `CombinedPolicy(P2,P1).truncate` yield equal
support sets (`retain_comm`); (ii) `sentinel_skip_is_identity` — a support
containing explicit zero-coefficient and maximal-weight terms is bit-for-bit
unchanged by `MaxPauliWeight(usize::MAX)` and `MaxLossWeight(usize::MAX)`, on both
HashMapStore and ColumnStore; (iii) the same conjunction test for a
`LossyPauliSum` under `MaxLossWeight(m)`.

### G-070 — truncate_preserve_eq_widened_retain — the best truncation theorem — has no `*_lean.rs` oracle

- Class `bridge` · Tier B · Severity medium · Sector truncation-policy-loss · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/sum.rs:488`
- Lean: `lean/PPVM/Algebra/Truncation.lean:315` `truncate_preserve_eq_widened_retain` (corollaries `truncatePreserve_apply_of_mem:336`, `truncatePreserve_empty:342`)

**Claim.** `sum.rs:478-486` — the snapshot/policy/restore composite equals one pass
with the widened keep-rule `keep ∨ k ∈ P`; a preserved key keeps *exactly* its
pre-truncate coefficient; the empty-keep-set fast path is exact.

**Why unverified.** The theorem is sound and independently motivated (it is what
makes the two invisible guards — the `storage.get(&key).is_none()` membership test
against `add_term`'s accumulation, and snapshot-before-policy — load-bearing rather
than incidental). But the only tests that touch `preserve` are
`pauli_sum_preserve_diff.rs` and `column_store_diff.rs:922`, both legacy diffs. So
the two mutations the theorem exists to catch — dropping the membership guard (a
survivor's coefficient DOUBLES) and moving the snapshot after the policy (a
preserved key is restored at its post-truncate residue, i.e. 0) — are caught today
only by legacy agreement. Additionally `sum.rs:474-477` asserts a restore
*ordering* property ("restored in their original relative order rather than the
keep-set hash table's iteration order"), observable through IndexMapStore's
insertion-ordered `iter()`/`Display`, which the Finsupp model cannot express and no
oracle checks.

**Proposed closure.** Add `truncate_preserve_is_widened_retain` to
`pauli_sum_lean.rs`: over random supports and both HashMapStore/IndexMapStore,
assert the post-truncate map equals the single-pass widened retain computed in the
test harness (`|c| >= t || preserve.contains(k)`) with exact bit equality on the
preserved keys' coefficients (catching the doubling), plus a case where a preserved
key survives the policy (guard must prevent doubling) and a case where a preserved
key is at a residue value (snapshot must be pre-policy). For the ordering claim,
assert IndexMapStore's `iter()` order of restored keys matches their original
relative order.

### G-071 — Mixture probability-cutoff truncation with renormalization has no Lean counterpart

- Class `coverage` · Tier A · Severity high · Sector truncation-policy-loss · Status **adjudicated-spec** (round 2 — [adj U4](#u4-mixture-weights-and-sampler))
- Rust: `crates/ppvm-tableau-2/src/mixture/data.rs:142`
- Lean: none

**Claim.** `GeneralizedTableauMixture::truncate` — drop entries with
`probability <= sum_cutoff`, then, if anything was dropped, **renormalize** the
remaining probabilities to sum to 1 — is a correct approximation with a bounded
error.

**Why unverified.** This is a sixth, structurally different retention strategy that
no reviewed Lean file touches
(`grep -i 'mixture|normaliz|probability' lean/**/*.lean` finds only `Noise.lean`'s
channel probabilities and `Projection.lean`'s amplitude normalize). Two features
put it outside every existing bound: (1) it *renormalizes*
(`mixture/data.rs:145-149` calls `normalize_probabilities`), so the incurred error
is NOT the dropped contribution that `l1_bound`/`l2_bound` bound — the surviving
entries are all rescaled by 1/(1−m), and the correct guarantee is a trace-distance
bound like 2m/(1−m), which is nowhere stated; (2) the cutoff is applied
asymmetrically at insertion — `mixture/data.rs:197-200` and `:244-248` apply
`probability > sum_cutoff` only to branches that fail the fingerprint join, so a
sub-cutoff branch that *collides* with an existing entry is accumulated regardless,
meaning the retained set is not characterizable by any predicate on the final
probabilities. Also note the strict `>` here versus PauliSum's `>=`, a third
instance of the boundary split that `cutoff_mismatch` only gestures at.
Measurement calls this (`mixture/measure.rs:105,119`) between reporting
pre-truncation mass and applying the cutoff, so it is on the observable path.

**Round-2 amendment.** `truncate` itself is **correct**: drop-then-renormalize is
the defensible policy, with trace distance ≤ m and |Δ⟨O⟩| ≤ 2m for m the dropped
mass; measured, it renormalizes exactly when it removes an entry
([0.6, 0.3, 0.1] → [0.6666666666666667, 0.33333333333333337], mass 1.0) and not
otherwise. Feature (2) of this row — the insertion-time cutoff "bypass" — is
**refuted**: accumulating a sub-cutoff branch onto a colliding live entry is exact
and lossless, and the cutoff properly applies to the merged total, which the
following `truncate()` re-tests. The real order-dependence of the retained set comes
from `structurally_equal`'s non-transitivity (G-020), not from the cutoff. The
genuine residue is that the bound is nowhere stated, and that the leak this row's
callers suffer lives in `for_each_z_branch` (G-018), not here. See
[adj U4](#u4-mixture-weights-and-sampler).

**Proposed closure — amended by round 2.** (i) stands. (ii) `2m/(1−m)` is valid but
loose: prove the tight elementary bound, trace distance ≤ m and
|⟨O⟩ − ⟨O⟩'| ≤ 2m for ‖O‖ ≤ 1 (from
‖ρ−ρ'‖₁ ≤ m + (m/(1−m))(1−m) = 2m). (iii) is **INVALIDATED** — do not prove the
insertion-time accumulation to be a defect; if an order-dependence witness is
wanted, state it against the non-transitive merge key (G-020). Oracle: keep the
`z_expectation`-before/after-truncate check against the proved bound, and keep the
"below-cutoff colliding branch is kept" test but label it as *characterizing exact
behaviour*, not as exhibiting a bug.

### G-072 — l2_bound is credited to the tableau path but bounds a linear functional; the tableau readout is sesquilinear

- Class `fidelity` · Tier A · Severity high · Sector truncation-policy-loss · Status **open**
- Rust: `crates/ppvm-tableau-2/src/data.rs:1240`
- Lean: `lean/PPVM/Algebra/Truncation.lean:219` `l2_bound` (and `l2_bound_normalized:226`)

**Claim.** "that is exactly the ℓ² truncation bound the stabilizer-tableau path
uses" (`Truncation.lean:217-218`) — i.e. the error the tableau's
`norm_sqr() > cutoff_sq` amplitude drop incurs is bounded by Cauchy–Schwarz on the
dropped coefficients.

**Why unverified.** `l2_bound` bounds `(Σ_{k∈D} c_k e_k)^2`, a squared *linear*
functional of the coefficients. But the generalized tableau's observable readout is
`⟨ψ|P|ψ⟩ = Σ_{α,β} ⟨α|P|β⟩ c_α* c_β`
(`crates/ppvm-tableau-2/src/expectation.rs:32-50`,
`compute_overlap_case_a`/`_b`), which is sesquilinear in the amplitude vector.
Deleting a set D of amplitudes therefore removes both the D×D block and the *cross*
block 2·Re Σ_{α∈D, β∉D} ⟨α|P|β⟩ c_α* c_β; `l2_bound` bounds neither. Second
mismatch: the tableau's inline cutoff (`data.rs:1240, 1308-1391, 1446`) and `trim`
(`data.rs:749`) drop amplitudes and do **not** renormalize, so the
post-truncation vector is sub-normalized and even the D×D contribution is measured
against a state whose norm has shifted — while the Lean docstring explicitly says
these bounds are "stated for … a normalized state" (`Truncation.lean:249-250`).
Finally there is no oracle at all: `tableau_lean.rs` contains no occurrence of
cutoff/threshold/trim, and `pauli_sum_lean.rs` never mentions l2.

**Proposed closure.** Either restate the tableau bound for the actual functional —
`|⟨ψ|P|ψ⟩ − ⟨ψ_D|P|ψ_D⟩| ≤ 2‖c_D‖₂‖c‖₂ + ‖c_D‖₂²` for ‖P‖ ≤ 1, derived from
Cauchy–Schwarz on both blocks — or re-scope `l2_bound`'s docstring to the PauliSum
linear-observable case and add the sesquilinear theorem separately. Oracle: a
`tableau_lean.rs` test that builds a magic-state superposition, records
`expectation(word)` before and after a `trim(cutoff)`, and asserts the difference
against the proved bound computed from the dropped amplitudes.

### G-073 — Lean's weight notion is a single-site Pauli count citing a legacy path; MaxPauliWeight on a lossy key counts Lost sites

- Class `fidelity` · Tier B · Severity medium · Sector truncation-policy-loss · Status **open**
- Rust: `crates/ppvm-lossy-pauli-word-2/src/data.rs:288`
- Lean: `lean/PPVM/Pauli.lean:73` `Pauli.weight` (used abstractly as `w : K → ℕ` in `GradedMap.lean:639` `retain_weight_le_eq_self`)

**Claim.** Lean docstring at `Pauli.lean:70-72`: "The n-qubit `weight()`
(`crates/ppvm-pauli-word/src/word/data.rs:150`) is the sum of these over all
slots" — i.e. the weight the `MaxPauliWeight` policy thresholds is the count of
non-identity Pauli slots.

**Why unverified.** Two problems. (1) The citation is to the legacy crate
`ppvm-pauli-word`, not to a `-2` target, and the `-2` code is not uniform:
`ppvm-pauli-word-2/src/data.rs:307` computes popcount(x|z), but
`ppvm-lossy-pauli-word-2/src/data.rs:288` computes popcount(x|z|**l**) — a Lost
site counts toward the Pauli weight (its own unit test at `data.rs:627` pins
`"XLIL".weight() == 3`). Since `MaxPauliWeight`'s bound is only
`W: Word + Indexable`, it is instantiable for `LossyPauliSum`, where its keep-rule
therefore drops on a quantity that is not "the sum of `Pauli.weight` over slots" —
the Lean sentence is false for that instantiation. (2) Nothing in Lean ever
*defines* the n-qubit weight function: `retain_weight_le_eq_self` abstracts it to
an arbitrary `w : K → ℕ`, so the only weight-policy theorems in the formalization
are the disable-sentinel identities, and no statement distinguishes `weight` from
`loss_weight` (the two policies' predicates, `policy.rs:168` vs `:321`).

**Proposed closure.** Define `wordWeight : (Fin n → Pauli) → ℕ` and
`lossyWeight`/`lossWeight` on a Lean lossy-site model (`Noise.lean` already has a
3-valued `Site` with `sL`), prove `lossyWeight = wordWeight_present + lossWeight`
and re-cite `Pauli.lean:71` to `crates/ppvm-pauli-word-2/src/data.rs:307` and
`crates/ppvm-lossy-pauli-word-2/src/data.rs:288` explicitly, noting the
Lost-counts-as-weight divergence. Oracle: extend `lossy_pauli_word_lean.rs` with an
exhaustive small-n check that `weight()` == the Lean `lossyWeight` and
`loss_weight()` == `lossWeight`, and a `LossyPauliSum` test showing
`MaxPauliWeight(w)` drops a word whose *present* Pauli weight is ≤ w but whose lost
count pushes it over.

### G-074 — l1_bound is proved for an arbitrary dropped set, never tied to any policy's keep-rule

- Class `strength` · Tier A · Severity high · Sector truncation-policy-loss · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/policy.rs:30`
- Lean: `lean/PPVM/Algebra/Truncation.lean:44` `l1_bound` (and `l1_bound_abv:82`)

**Claim.** "The truncation error a policy incurs is bounded" — i.e. each shipped
`Policy` has an error guarantee expressed in its own parameter (threshold t, weight
cap w, loss cap m).

**Why unverified.** `l1_bound`/`l1_bound_abv` quantify over an arbitrary `Finset D`
and bound |Σ_{k∈D} c_k e_k| ≤ Σ_{k∈D} N(c_k). That is the triangle inequality; it
holds for D chosen by any rule whatsoever, including "drop the largest
coefficients". Nothing in `Truncation.lean` ever instantiates D as the set the Rust
predicate rejects. There is no theorem
D ⊆ {k : N(c_k) < t} → Σ_{k∈D} N(c_k) ≤ D.card * t, which is the only statement
that would make `CoefficientThreshold`'s parameter mean anything. Worse, the claim
is attached to the *trait* doc (`policy.rs:30`, on `trait Policy`), so it is
credited to `MaxPauliWeight` (`policy.rs:168`) and `MaxLossWeight`
(`policy.rs:321`) too — and for those the dropped ℓ¹ mass is provably unbounded by
the parameter (a single weight-(w+1) term can carry coefficient 10^6). So the two
policies the headline workloads actually configure
(`CombinedPolicy(CoefficientThreshold(1e-6), MaxPauliWeight(w))` in
`pauli_sum_workload_diff.rs:367`) have no parameter-indexed error bound at all.

**Proposed closure.** State in `Truncation.lean`: (i)
`threshold_dropped_l1_le (t : ℝ) (D : Finset K) (hD : ∀ k ∈ D, N (c k) < t) : ∑ k ∈ D, N (c k) ≤ D.card * t`,
composed with `l1_bound_abv` into
`coefficientThreshold_error_le : N (Σ_{k∈D} c_k e_k) ≤ D.card * t`; (ii) an
explicit *counterexample* theorem
`maxPauliWeight_error_unbounded : ∀ M, ∃ c w, (∀ k ∈ D, weight k > cap) ∧ Σ_{k∈D} |c k| > M`,
so the weight-cap policies are documented as carrying no coefficient-mass
guarantee. Oracle: extend `truncation_l1_bound` to build the drop set from
`CoefficientThreshold` and assert
`error.abs() <= (dropped.len() as f64) * threshold`, which fails if the Rust
keep-rule is loosened.

### G-075 — No multi-step truncation bound; the cited telescoping fails for the rotation gate that drives truncation

- Class `strength` · Tier A · Severity high · Sector truncation-policy-loss · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/sum.rs:629`
- Lean: `lean/PPVM/Algebra/Noise.lean:117` `l1_contractive`

**Claim.** `sum.rs:629` — `l1_contractive` "is what makes `Truncation.l1_bound`
compose over a noisy circuit", i.e. the per-truncation ℓ¹ error telescopes over the
~1500 `truncate()` calls of the headline Trotter workload.

**Why unverified.** There is no theorem anywhere that states a multi-step bound: no
induction over a gate list, no `Σ_i ‖D_i‖₁` statement, nothing quantifying over a
sequence of propagate-then-truncate steps. `l1_contractive` (`Noise.lean:117`) is a
single-step, per-key statement `Σ|λ_k c_k| ≤ Σ|c_k|` for `|λ_k| ≤ 1` — it covers
only the *diagonal* `scale_by_key` channel. The gate that actually creates the
terms truncation exists to remove is the rotation, and it is ℓ¹-**expanding**:
`rotation.rs:175-176` leaves `c·cos` on the diagonal key and merges `c·sin·ε` onto
the branch key, so one rotation multiplies the ℓ¹ mass of an unmerged term by
|cosθ|+|sinθ| ∈ [1, √2]. Over R rotation gates the amplification is up to 2^{R/2},
and dropped mass that is measured before those gates is *not* an upper bound on the
final observable error. So the composition claim is asserted in a Rust docstring,
supported by a lemma about a different gate class, with the interesting hypothesis
(non-expansiveness of the gate mix) simply absent — and no oracle measures
accumulated error over a circuit.

**Proposed closure.** Prove
`l1_expansion_rotation : ‖rot θ A‖₁ ≤ (|cos θ| + |sin θ|)·‖A‖₁` and its tightness
witness, then a genuine telescoping theorem
`truncate_propagate_error : |⟨O⟩_exact − ⟨O⟩_trunc| ≤ Σ_i (Π_{j>i} κ_j)·‖D_i‖₁`
where κ_j is the step's ℓ¹ operator norm (1 for Clifford re-key by
`pushforward_apply`, 1 for a contractive channel, |cos|+|sin| for a rotation).
Oracle: a `pauli_sum_lean.rs` test running a fixed rotation circuit twice
(`NoPolicy` vs `CoefficientThreshold`), accumulating the per-truncate dropped ℓ¹
mass with the proved amplification factors, and asserting the final `overlap`
difference against that accumulated bound.

### G-076 — Batch order/partition invariance is proved over an exact AddCommMonoid; the shipped coefficient is f64

- Class `strength` · Tier B · Severity medium · Sector truncation-policy-loss · Status **open**
- Rust: `crates/ppvm-traits-2/src/graded.rs:99`
- Lean: `lean/PPVM/Algebra/GradedMap.lean:149` `accumulateTerms_perm` (and `accumulateTerms_add:159`)

**Claim.** `graded.rs:91-98` — "an impl is free to reorder [the batch] and to split
it across partitions/threads: order-invariance (`accumulateTerms_perm`) and
partition-invariance (`accumulateTerms_add`) are machine-checked".

**Why unverified.** Both theorems are stated for `CMap K C` with `C` an
`AddCommMonoid`, where `+` is exactly associative. The shipped coefficient rings
are `f64` and `Complex<f64>`, whose addition is neither associative nor (in the
presence of subnormals/cancellation) reorder-invariant, so the theorems do not
transfer to the license the docstring grants. This is not academic here: the
crate's own tests demand bit-exact agreement in places
(`column_store_lean.rs`'s `assert_aligned` compares `to_bits()`; the truncation
diff tests place boundary terms at |c| == τ and τ − 1 ulp), and the composition
with truncation is discontinuous — a reordered accumulation can move a merged
coefficient across the `CoefficientThreshold` keep-rule and change the *support*,
not just the last bits. The bridge is correspondingly soft:
`pauli_sum_lean.rs:73-105` reverses the term list and regroups it into three blocks
but asserts only `assert_maps_close(.., 1e-9)`, with no error model and no
truncation in the loop, so it would not detect a backend that genuinely reordered.
The same license is what `HashMapStore`'s data-dependent merge direction
(`store.rs:1551`, justified at `store.rs:461-463`) leans on.

**Proposed closure.** Either scope the Rust docstring's license to exact rings and
add a Lean witness `float_reassoc_ne : ∃ a b c : Float, (a+b)+c ≠ a+(b+c)` plus a
rounding-error statement `‖fold_l₁ B − fold_l₂ B‖ ≤ (|B|−1)·ε·‖B‖₁` for a
permutation of the same multiset, or restate `accumulateTerms_perm` over a
`LinearOrderedField`-with-error model. Oracle: a `pauli_sum_lean.rs` test that (i)
asserts bit-exact equality for an integer/rational coefficient instantiation under
reversal and block regrouping, (ii) for f64 asserts the reordered result within the
proved (|B|−1)·ε·ℓ¹ envelope, and (iii) runs `truncate()` after each ordering under
`CoefficientThreshold` and asserts the *support sets* agree — which is the property
the reorder license actually needs.

### G-077 — Every `*_lean.rs` oracle instantiates PauliWord<u64>, so multi-limb key_mul is unpinned

- Class `bridge` · Tier B · Severity medium · Sector word-algebra · Status **open**
- Rust: `crates/ppvm-pauli-word-2/src/product.rs:58`
- Lean: `lean/PPVM/Pauli/Word.lean:55` `PPVM.PauliWord.phaseExpN`

**Claim.** `key_mul`'s loop over all storage limbs accumulates sign/imag popcounts
across limbs before the single mod-4 reduction, so the emitted `Phase` equals
`phaseExpN` for words wider than one storage limb.

**Why unverified.** `pauli_word_lean.rs:33` is `type New = NewWord<u64>` and
`phased_pauli_word_lean.rs:34` uses the default `PhasedPauliWord`
(= `Phased<PauliWord<DefaultStorage = u64>>`). For storage `u64`,
`<u64 as BitView>::Store = u64`, so `as_raw_slice().len() == 1` and the loop body
at `product.rs:58-70` executes exactly once at every tested width (the oracles top
out at n = 60). The only other lean-oracle instantiation is `PauliWord<[u8; 8]>` in
`sym_lean.rs:49` and there n ≤ 6, i.e. limb 0 only. So mutating `product.rs:58` to
`for i in 0..1` leaves all nine `*_lean.rs` files green. Production keys are
multi-limb: `IndexPauliSum<N>` is `PauliWord<[u8; N]>`
(`ppvm-pauli-sum-2/src/lib.rs:117`) with `Store = u8`, i.e. N limbs of 8 bits, and
`pauli_sum_hash.rs:200-202` exercises `[u8; 8]`/`[u8; 16]`/`[u8; 32]` at 64/128/256
qubits — but only for hash distribution. Multi-limb `key_mul` is reached only by
`*_diff.rs` tests, and since the `-2` sign/imag/accumulate kernel is byte-identical
to legacy `ppvm-pauli-word/src/phase/mul.rs`, a diff test cannot detect a kernel
defect by construction. The residual risk is narrow (per-limb indexing, not the
cross-qubit accumulation, which n=60 in one u64 does exercise) but it is exactly the
part with no independent grounding.

**Proposed closure.** No new Lean needed — `phaseExpN` already quantifies over all
n. Add to `pauli_word_lean.rs` a second instantiation
`type Wide = NewWord<[u8; 16]>` and run
`phase_exp_equals_matrix_exponent_n_qubit_random` (n ≤ 5 for the matrix arm) plus
`phase_cocycle_associativity_random`, `square_is_plus_identity_random`,
`commutation_is_symplectic_form_random` at n in {9, 16, 17, 63, 64, 65, 128} so limb
boundaries at 8, 64 and the final limb are all crossed; mirror it in
`phased_pauli_word_lean.rs` with `Phased<PauliWord<[u8; 16]>>`.

### G-078 — Canonical-unused-bits / tail-limb padding invariant appears nowhere in Lean

- Class `coverage` · Tier B · Severity medium · Sector word-algebra · Status **open**
- Rust: `crates/ppvm-pauli-word-2/src/data.rs:38`
- Lean: none

**Claim.** For a word of width n stored in a blob of `8*size_of::<A>()` bits, every
bit at index >= n is permanently 0 in both planes; therefore (a) popcounting the
whole blob equals summing over the n logical qubits, (b) comparing/hashing the whole
blob equals comparing/hashing the logical bits, and (c) `key_mul`'s XOR preserves
the invariant.

**Why unverified.** `data.rs:36-42` names this the "canonical-unused-bits
invariant" and it is load-bearing in four places: `key_mul` popcounts
`2*8*size_of::<A>()` bits (`product.rs:58-70`), `weight()` popcounts the whole blob
(`data.rs:307-338`), `PartialEq` compares `xbits.data`/`zbits.data` in full
(`data.rs:534-541`), and `structural_hash` hashes `bytemuck::bytes_of` of the full
blob including padding (`hash.rs:46-47`). The Lean model has no notion of a blob at
all: `PPVM.PauliWord.Word n := Fin n -> Bool x Bool` (`Word.lean:47`) is exactly n
slots, so `phaseExpN` (`Word.lean:55`) can only model `key_mul` UNDER this
invariant, and the invariant is stated in no Lean definition or theorem — grep for
`unused|padding|tail|limb|nqubits|capacity` over `lean/PPVM` returns one unrelated
hit in `Tableau/Batch.lean:19`. There is also no Rust test for it: the lossy word
has `lost_site_has_canonical_bits_and_loss_is_exclusive`
(`lossy_pauli_word_lean.rs:356`), the bare word has no analogue, and no test in the
repo inspects `as_raw_slice`/`bytes_of` of a `PauliWord` to check high bits are
zero. (G-009 is the digest-side consequence of the same missing invariant; closing
one should close both.)

**Proposed closure.** Lean: add `PPVM.PauliWord.padWord : Word n -> Word m`
(m >= n, pad with `(false,false)`) plus (i)
`phaseExpN_pad : phaseExpN (padWord p) (padWord q) = phaseExpN p q`, which reduces
to `phaseExp false false false false = 0` (already available as
`phaseExp_id_left`), and (ii)
`mulWord_padded : IsPadded n p -> IsPadded n q -> IsPadded n (mulWord p q)`. Oracle:
in `pauli_word_lean.rs` add a `padding_bits_stay_zero` test over
`PauliWord<[u8; 16]>` at widths that are not multiples of the limb width
(n = 5, 37, 100, 127), asserting after `From<&str>`, every
`set_*`/`toggled_bits*` mutator, and `key_mul` that
`bytemuck::bytes_of(&w.raw_planes())` has all bits >= n clear — plus assert that
`key_mul` on such a word still matches the existing ℤ[i] matrix reference. Also
promote the width `debug_assert_eq!` at `product.rs:42` to a real `assert_eq!` or
document it as an unchecked precondition, since in release a width mismatch XORs the
wider operand's live bits into the narrower result's padding and violates (a)-(c)
silently.

### G-079 — weight()'s byte-cascade popcount is never tied to the Lean weight it feeds

- Class `coverage` · Tier B · Severity medium · Sector word-algebra · Status **open**
- Rust: `crates/ppvm-pauli-word-2/src/data.rs:307`
- Lean: `lean/PPVM/Pauli.lean:73` `PPVM.Pauli.weight`

**Claim.** `PauliWord::weight()` equals the number of non-identity sites, i.e. the
sum over i < n of `Pauli.weight (get i)` — the `w : K -> Nat` that the truncation
theorems are stated over.

**Why unverified.** The Rust is a fused popcount of `x | z` over the raw blob via an
8/4/2/1-byte cascade with `u64/u32/u16/u8::from_ne_bytes` (`data.rs:312-337`) —
endianness- and padding-sensitive, and structurally nothing like a per-site fold.
Lean defines only the SINGLE-site `Pauli.weight` and makes the n-qubit claim in
prose ("The n-qubit `weight()` ... is the sum of these over all slots"), citing the
LEGACY `crates/ppvm-pauli-word/src/word/data.rs:150`; there is no `weightN`
definition and no theorem. `lean/PPVM/Algebra/GradedMap.lean:627-643`
(`retain_weight_le_eq_self`, the `usize::MAX` sentinel) quantifies over an abstract
`w : K -> Nat`, so the truncation results never touch the concrete function. On the
oracle side, no `*_lean.rs` test asserts `w.weight()` against a per-site count:
`pauli_sum_lean.rs:288` only checks `w.weight() <= cap` after truncation (which
would still pass if `weight` systematically under-reported), and the only per-site
check in the sector is the in-crate unit test `weight_counts_nonidentity`
(`data.rs:675`) on a single 5-qubit literal. A `weight` that dropped the trailing
`if i < n` byte, or that read `x & z` instead of `x | z`, would pass every
Lean-oracle test.

**Proposed closure.** Lean: add
`PPVM.PauliWord.weightN (p : Word n) : Nat := sum over i of (if p i = (false,false) then 0 else 1)`
and prove `weightN p = card {i | (p i).1 || (p i).2}` (the popcount-of-`x|z` form),
plus `weightN_le : weightN p <= n`; then instantiate `GradedMap`'s
`retain_weight_le_eq_self` at `w := weightN` so the truncation theorems name the
concrete function. Oracle: in `pauli_word_lean.rs` add `weight_equals_site_count` —
seeded random words at n in {1,7,8,9,63,64,65,100} over both `PauliWord<u64>` and
`PauliWord<[u8; 16]>`, asserting
`w.weight() == (0..n).filter(|i| w.get(i) != Pauli::I).count()`.

### G-080 — pauli_code's endian-dependent byte extraction has no Lean model and no agreement test

- Class `coverage` · Tier B · Severity low · Sector word-algebra · Status **open**
- Rust: `crates/ppvm-pauli-word-2/src/data.rs:364`
- Lean: none

**Claim.** `PauliWord::pauli_code(i)` equals the trait default
`x_bit(i) | (z_bit(i) << 1)`, i.e. the packed 2-bit code `0=I, 1=X, 2=Z, 3=Y`, for
every backing storage `A` and every i < n.

**Why unverified.** The override replaces logical bit reads with a raw byte read:
`bytemuck::bytes_of(&self.xbits.data)[i >> 3] >> (i & 7) & 1` under
`cfg(target_endian = "little")`, falling back to the bit-indexed form on
big-endian. Its correctness is a statement about bitvec's `Lsb0` element ordering
composed with the little-endian in-memory byte order of the `Store` type — for
`Store = u64` it holds only because element e's byte b sits at offset 8e+b. Nothing
in Lean models the packed byte layout, so there is no theorem to pin it to; and
grepping `pauli_code` across `crates/ppvm-conformance-2/tests` and
`crates/ppvm-pauli-word-2` finds no assertion at all — no unit test in `data.rs`'s
test module, no oracle test. Yet it is the dispatch input for every two-qubit
Clifford re-key (`ppvm-pauli-sum-2/src/clifford.rs:84-85, 258-259, 273-274,
291-292, 306-307, 504-505`), the column-store site scan (`store.rs:289, 1229`),
noise (`noise.rs:90`) and pattern matching (`pattern/matches.rs:141,149`). A
swapped shift or a wrong byte index on a `[u64; N]` storage would be caught only
indirectly, by `*_diff.rs` agreement with a legacy crate that has the same fast
path.

**Proposed closure.** Lean: state the encoding once —
`PPVM.PauliWord.pauliCode (b : Bool x Bool) : Fin 4 := ...` with
`pauliCode_eq : pauliCode b = (if b.1 then 1 else 0) + 2 * (if b.2 then 1 else 0)`
and a lemma tying it to `PPVM.Pauli` (`I=0, X=1, Z=2, Y=3`), so the code is a named
spec rather than an implicit convention. Oracle: in `pauli_word_lean.rs` add
`pauli_code_matches_bit_reads` — for `PauliWord<u64>`, `PauliWord<[u8; 16]>` and
`PauliWord<[u64; 4]>` at n in {1,8,9,63,64,65,200}, assert for every i that
`w.pauli_code(i) == (w.x_bit(i) as u8) | ((w.z_bit(i) as u8) << 1)` and that it
agrees with the letter `w.get(i)`.

### G-081 — All four sector modules cite only the legacy crate; the Pauli discriminant table is now wrong

- Class `fidelity` · Tier A · Severity low · Sector word-algebra · Status **open**
- Rust: `crates/ppvm-traits-2/src/word.rs:47`
- Lean: `lean/PPVM/Pauli.lean:15` (module docstring table) and `PPVM.Pauli.mul` at `:62`

**Claim.** The Lean sector modules model the code under verification, the `-2`
crates.

**Why unverified.** `Phase.lean`'s module docstring and section headers cite
`crates/ppvm-pauli-word/src/phase/mul.rs:42` five times
(`:17, :65, :71, :77, :143, :214, :249`); `Word.lean:14-21` quotes legacy source
text verbatim (`self.add_phase(((2 * sign_count + imag_count) % 4) as u8);`) as the
kernel being modeled; `Pauli.lean` cites `ppvm-pauli-word/src/word/mul.rs`,
`ppvm-pauli-word/src/word/data.rs:150` and
`ppvm-traits/src/traits/word_trait.rs:113`. None of the four sector modules cites
`crates/ppvm-pauli-word-2/src/product.rs` or
`crates/ppvm-phased-pauli-word-2/src/product.rs`. The arithmetic was checked and
the sign/imag formulas are byte-identical, so this is not a math error — but two
concrete referents HAVE drifted. (1) `Pauli.lean:15-20` anchors the bit encoding to
"Rust discriminant (`ppvm-traits/src/char.rs`)" with
`I=0b00, X=0b01, Z=0b10, Y=0b11`; that legacy enum does declare `Z = 2, Y = 3`, but
`ppvm-traits-2/src/word.rs:47` declares the variants in the order `I, X, Y, Z`, so
`Pauli as u8` now gives `Y = 2, Z = 3` and the packed code survives only in
`PauliBits::pauli_code` (`word.rs:185`) plus a separate remap at
`ppvm-pauli-sum-2/src/rotation.rs:58`. The Lean table is a false statement about the
target crate. (2) In `-2` the product is split — the bare word emits a residual
`Phase` from `key_mul` (`product.rs:74`) and `Phased::mul` composes
`self.phase * rhs.phase * emitted` (phased `product.rs:48`) — so the
`MulAssign`/`add_phase` shape `Word.lean` quotes as "the Rust kernel" no longer
exists in the target. The direction of the citations also matters for provenance:
the reader cannot tell from the Lean that the `-2` code was ever read.

**Proposed closure.** Lean: retarget every Rust citation in `Pauli.lean`,
`Phase.lean`, `Word.lean` and `Matrix.lean` to the `-2` paths
(`crates/ppvm-pauli-word-2/src/product.rs:64` for sign/imag, `:73-74` for the mod-4
reduction and the emitted `Phase`,
`crates/ppvm-phased-pauli-word-2/src/product.rs:44-50` for the phased product,
`crates/ppvm-pauli-word-2/src/data.rs:295-301` for the (x,z) letter encoding); and
replace `Pauli.lean`'s discriminant column with the `PauliBits::pauli_code`
encoding (`ppvm-traits-2/src/word.rs:183-186`), noting explicitly that the `-2`
`Pauli` enum's discriminant is NOT the packed code. Oracle: add a
`pauli_letter_encoding` assertion in `pauli_word_lean.rs` pinning
`(x_bit, z_bit)` -> `Pauli` for all four patterns and `pauli_code` ->
`{I:0, X:1, Z:2, Y:3}`, so the table in Lean is test-anchored rather than prose.

### G-082 — n-qubit phased-product unit and inverse laws are credited to single-qubit theorems

- Class `strength` · Tier A · Severity medium · Sector word-algebra · Status **open**
- Rust: `crates/ppvm-phased-pauli-word-2/src/product.rs:44`
- Lean: `lean/PPVM/Pauli/Phase.lean:241` `PPVM.PauliPhase.PhasedPauli.mul_one'` (and `one_mul'`, `inv_mul_cancel'`, the `Group` instance at `:250`)

**Claim.** `Phased<PauliWord>::mul` is the group operation of the n-qubit phased
Pauli group P_n: `+I...I` is a two-sided unit and `(phi, w)^-1 = (-phi, w)` is a
two-sided inverse.

**Why unverified.** `Word.lean` lifts exactly three facts to general n —
`phaseExpN_cocycle` (`:64`), `phaseExpN_sub_comm` (`:73`), `phaseExpN_self`
(`:136`) — and there is no `phaseExpN_id_left`/`phaseExpN_id_right` and no n-qubit
`Group`/`Monoid` instance anywhere in the sector; the only bundled group is
`PhasedPauli` at n = 1 (`Phase.lean:202-257`), whose `phase : ZMod 4` field has no
counterpart in the bare `PauliWord` and whose `mul` is the single-qubit `phaseExp`.
Meanwhile `phased_pauli_word_lean.rs:194` (`phased_product_identity_laws`) and
`:212` (`phased_product_inverse_law`) cite `Phase.lean one_mul'/mul_one'` and
`inv_mul_cancel'` by name and run them at n in {1,2,3,5,16,60}. So for n > 1 the
unit and inverse laws are asserted by a Rust test whose named Lean referent does not
cover the instance under test. `Word.lean`'s `Canon` monoid (`:211-241`) is not a
substitute: it is the tableau's `i^phi X^x Z^z` normalization with the DIFFERENT
`crossPhase` cocycle, not `phaseExp`. The missing lemmas are one-liners
(`Finset.sum_eq_zero` + `phaseExp_id_left/right`), which is what makes this a
strength rather than a coverage gap.

**Proposed closure.** Lean: in `Word.lean` add
`phaseExpN_id_left : phaseExpN (fun _ => (false,false)) q = 0` and
`phaseExpN_id_right`, each by `Finset.sum_eq_zero` over
`phaseExp_id_left/right`; then define `PhasedWord n` (a `ZMod 4` phase plus
`Word n`) with `mul`/`one`/`inv` mirroring
`phased-pauli-word-2/src/product.rs:44-50` and install a genuine
`Group (PhasedWord n)` instance whose axioms come from `phaseExpN_cocycle`, the two
new identity lemmas, and `phaseExpN_self`. Oracle: retarget the docstring citations
of `phased_product_identity_laws` / `phased_product_inverse_law` to the new n-qubit
theorem names.

### G-083 — PauliSum's extension-gate sign kernels (s_dag/sqrt_x/sqrt_y/cy) reach no `*_lean.rs` oracle

- Class `bridge` · Tier A · Severity medium · Sector word-and-clifford (skeptic) · Status **open**
- Rust: `crates/ppvm-pauli-sum-2/src/clifford.rs:350`
- Lean: `lean/PPVM/Pauli/Conjugation.lean` `extSqrtX_sign` / `extSqrtY_sign` / `conjCY_sign`

**Claim.** `Sum`'s own fused re-key closures for
`s_dag`/`sqrt_x`/`sqrt_x_dag`/`sqrt_y`/`sqrt_y_dag`/`cy` emit the conjugation sign
of `Conjugation.lean`'s closed-form predicates (`extSqrtX_sign` = x∧z,
`extSqrtXdag_sign` = ¬x∧z, `extSqrtY_sign` = ¬x∧z, `extSqrtYdag_sign` = x∧¬z,
`conjCY_sign`) onto the coefficient.

**Why unverified.** `ppvm-pauli-sum-2/src/clifford.rs:318-517` is a FOURTH
independent transcription of these predicates (after `traits-2`'s blanket,
`phased-pauli-word-2`'s fused kernels and `pauli-word-2`'s own gates):
`clifford.rs:336` `s_dag`, `:350` `sqrt_x`, `:365` `sqrt_x_dag`, `:380` `sqrt_y`,
`:395` `sqrt_y_dag`, `:447` `cy` plus the shared `cy_toggles` at `:77-92`, each
hand-writing its own `pauli_code`-based bit reads and its own negate condition and
taking the plain-insert `rekey_bijective` fast path. grep for
`sqrt_x|sqrt_y|s_dag|\.cy(` over `pauli_sum_lean.rs`,
`pauli_sum_multiply_lean.rs`, `pauli_sum_rotation_noise_lean.rs` and `sym_lean.rs`
returns NOTHING — no lean oracle ever applies an extension gate to a `Sum`.
`clifford_rekey_is_support_preserving_bijection` (`pauli_sum_lean.rs:342`) cycles
i%5 over h/s/cnot/cz/z, and support-size preservation cannot see a sign anyway. The
only coverage is `pauli_sum_gate_surface_diff.rs` (legacy, `const N = 4`) against
`ppvm-pauli-sum/src/sum/clifford.rs:107-221`, from which this code is a direct port
— so a sign transposition between the `sqrt_y`/`sqrt_y_dag` pair (identical bit
maps, differing only in the predicate, exactly the case `Conjugation.lean:411`
flags as "the invisible bug this theorem rules out") is caught by nothing but port
agreement.

**Proposed closure.** No new Lean needed. Add
`sum_extension_gate_signs_match_conjugation_oracle` to `pauli_sum_lean.rs`: for
each of the 4 single-site and 16 two-site basis words, build a one-term `PauliSum`
with coefficient +1, apply
`s_dag`/`sqrt_x`/`sqrt_x_dag`/`sqrt_y`/`sqrt_y_dag`/`cy`, and assert the resulting
(key, coefficient sign) equals the Rust mirrors of
`extSqrtX_bits`/`extSqrtX_sign`, …, `conjCY_bits`/`conjCY_sign` already transcribed
in `phased_pauli_word_lean.rs`; then extend
`clifford_rekey_is_support_preserving_bijection` to those six gates so the bijection
claim at `clifford.rs:33` covers its whole gate surface.

## Sector coverage

Twelve sectors were surveyed end to end. The one-line verdicts:

| sector | verdict |
| --- | --- |
| word-algebra | Kernel solid and matrix-grounded; the packed-storage layer is invisible to Lean. |
| clifford-conjugation | Generators fully grounded both halves; extensions/phaseless words rest on legacy diffs. |
| lossy-word | Real loss theory in Lean, but loss is an external parameter and canonicality is unchecked post-propagation. |
| graded-algebra-containers | Batch/Retain/Twisted layers are genuine; L0/L2/`reduce` are Mathlib one-liners and the L4 impls are untested. |
| sum-engine-stores | Abstract algebra well modelled; no store refinement exists at all. |
| truncation-policy-loss | The policy *shape* is faithfully modelled; the error bounds are not tied to any policy parameter. |
| multiply-rotation | Product excellent and matrix-grounded; rotation direction and branch signs are effectively unverified. |
| noise-observables | Eigenvalue arithmetic and index lists solid; the physics step, amplitude damping and both mixture/trajectory paths are absent. |
| hashing-digests | Contract tests are good; there is no Lean leg whatsoever, and the digests fold padding the model cannot express. |
| tableau-core | Bit-level frame theory strong and well bridged; the ℤ/4 sign column has no invariant and no ext-gate oracle. |
| measurement-branching | Case-split shape and frame coordinates proven; every measurement sign, the Born 1/2, and all of `mixture/` are open. |
| symbolic-coefficients | Genuinely independent ring formalization; the hom law is never pinned to the map-backed product. |

### word-algebra

- Rust surveyed: `ppvm-pauli-word-2/src/{product,data,storage,hash}.rs`,
  `ppvm-phased-pauli-word-2/src/{product,data}.rs`,
  `ppvm-traits-2/src/{algebra,word}.rs` (legacy `ppvm-pauli-word/src/phase/mul.rs`
  and `ppvm-traits/src/char.rs` read only for the fidelity comparison).
- Lean: `PPVM/Pauli/{Phase,Word,Matrix}.lean`, `PPVM/Algebra/Twisted.lean`,
  `PPVM/Pauli.lean`.
- Oracles: `pauli_word_lean.rs`, `phased_pauli_word_lean.rs` (storage widths also
  observed in `sym_lean.rs`, `pauli_sum_multiply_lean.rs`, `column_store_lean.rs`,
  `pauli_sum_hash.rs`).

**Solidly verified.** The per-qubit phase kernel is the best-grounded thing in the
repo. `product.rs:64-65` computes the same `sign`/`imag` sum-of-products as legacy
(diffed byte for byte — identical), and `Phase.lean`'s `phaseExp_eq_ref` proves
those booleans equal an analytic `phaseRef`, which `Matrix.lean`'s `pauliMat_mul`
then pins to a genuine 2x2 ℤ[i] matrix product (exhaustive `decide` over all 16
patterns) — so this is NOT laundered from legacy behaviour. The n-qubit lift is
genuinely universal rather than witnessed: `phaseExpN` is a `Finset` sum over
`Fin n`, and `phaseExpN_cocycle`, `phaseExpN_self`, `phaseExpN_sub_comm`,
`mulWord_assoc`, `mulWord_right_injective` hold for all n, with `tensorPauli_mul` +
`prod_iuPow` closing phase-multiplicativity against honest 2ⁿ×2ⁿ Kronecker
matrices. `Twisted.lean` states the associativity obligation key-agnostically
(`gtmul_assoc`, `IsCocycle`), derives the Pauli case as an instance, and separately
proves right-cancellativity is independent of it. Both oracle files re-implement the
ℤ[i] matrix reference in Rust and run it on the REAL `key_mul` / `Phased::mul`
(exhaustive 16 single-qubit, exhaustive 256 phase×Pauli, seeded n≤5 vs Kronecker
matrices, seeded n≤60 for the cocycle/self/commutation laws).

### clifford-conjugation

- Rust surveyed: `ppvm-traits-2/src/{gates,pauli}.rs`,
  `ppvm-pauli-word-2/src/clifford.rs`,
  `ppvm-phased-pauli-word-2/src/clifford.rs`,
  `ppvm-lossy-pauli-word-2/src/clifford.rs`, `ppvm-tableau-2/src/clifford.rs`
  (convention comparison), `ppvm-pauli-sum-2/src/clifford.rs`,
  `ppvm-traits-2/tests/{phase1_gate_surface,phase1_leaf_types}.rs`.
- Lean: `PPVM/Pauli/Conjugation.lean` (884 lines), `PPVM/Pauli/Symplectic.lean`
  (591 lines), `PPVM/Pauli/Matrix.lean`, `PPVM/Tableau/Batch.lean` (`IsSitewise`).
- Oracles: `pauli_word_lean.rs`, `phased_pauli_word_lean.rs`,
  `lossy_pauli_word_lean.rs`, `tableau_lean.rs`, `pauli_sum_lean.rs`.

**Solidly verified.** The gate surface is 13 gates + 4 stim aliases and contains no
SWAP/ISWAP (grep-confirmed); `T`/`T†` live in the non-Clifford `TGate` trait. Both
halves — bits AND sign — of all seven `Clifford` generators are pinned on the real
`Phased` type against genuine ℤ[i] `G†PG` matrices, exhaustively at n=1 and n=2, and
all six `CliffordExtensions` gates likewise (√X/√Y as `exp(-iπP/4)` scaled by √2) —
real independent mathematics, not a legacy table. Lean gives genuine group
homomorphisms of 𝒫₁/𝒫₂ for conjH, conjS, conjSdag, conjPauliZ, conjCNOT, conjCZ,
conjST, conjSdagT and conjCY — the two-qubit ones proved structurally through a
2-cocycle-compatibility lemma rather than brute `decide` — plus
injectivity/involutivity, `conjS` order 4, and `conjS_conjSdag` fixing the
convention-sensitive backward `S` sign. The six extension gates are *derived* as
generator products with closed-form `_bits`/`_sign` theorems, dagger-inverse pairs,
√X²=X, √Y²=Y, and `IsRealPhase` closure. `Sp(2n,2)` isometry holds at general n with
involutivity/bijectivity for hAct/sAct/cnotAct/czAct, and the whole lossy-guard
theory is machine-checked. All 13 gates were hand-checked across the three word
`clifford.rs` files and agree with each other and with Lean.

### lossy-word

- Rust surveyed: `ppvm-lossy-pauli-word-2/src/{lib,data,clifford,column,hash}.rs`,
  `ppvm-traits-2/src/{word,pauli,batch}.rs`,
  `ppvm-phased-pauli-word-2/src/product.rs`,
  `ppvm-pauli-sum-2/src/{rotation,loss,noise,lib}.rs`,
  `ppvm-pauli-sum-2/src/column_store/rotations/{mod,rx}.rs`.
- Lean: `PPVM/Pauli/Symplectic.lean:294-591` (the loss section),
  `PPVM/Algebra/Noise.lean:290-455` (the lossy alphabet and channels),
  `PPVM/Algebra/GradedMap.lean` (abstract weight).
- Oracles: `lossy_pauli_word_lean.rs` (8 tests), `lossy_pauli_word_diff.rs`,
  `pauli_sum_loss_diff.rs`.

**Solidly verified.** The sector is not Lean-free: `Symplectic.lean:294-591` has a
substantial loss section — `LossInv`, the four loss-guarded generators, invariant
preservation, the critical `cnotActL_lost_target_stays_identity`, present-sub-block
isometries, `xorZColL_xorXColL_eq_cnotActL`, and the genuinely hard
`sActL_cnotActL_sActL_eq_cyActL` (lost control + present target must *cancel*, not
skip). That CY decomposition was checked against `ppvm-traits-2/src/pauli.rs:297`
plus the lossy S guard and the Lean is faithful. `Noise.lean:296-455` independently
grounds the lossy *alphabet* (`Fin 5`, `I` = qubit-subspace identity, `L` = loss
projector, `𝟙 = I + L`) and proves the three loss channels trace-preserving — real,
non-legacy mathematics. So "what the lossy alphabet denotes" IS stated, in the noise
file.

### graded-algebra-containers

- Rust surveyed: `ppvm-traits-2/src/{graded,algebra,coefficient,word,batch}.rs`,
  `ppvm-traits-2/src/containers/{mod,hash_join,coordinate_list,tests}.rs`,
  `ppvm-traits-2/tests/{phase1_containers,phase1_leaf_types}.rs`,
  `ppvm-pauli-sum-2/src/{store,multiply}.rs`, `ppvm-tableau-2/src/{lib,data}.rs`
  (`Amplitudes`).
- Lean: `PPVM/Algebra/GradedMap.lean`, `PPVM/Algebra/Twisted.lean`,
  `PPVM/Pauli/Matrix.lean`, `PPVM/Instantiations/Projector.lean`.
- Oracles: `pauli_sum_lean.rs`, `column_store_lean.rs`, `sym_lean.rs`
  (`ImaginaryUnit`/`Conjugate`/`Phase` sections), `ppvm-traits-2`'s own container
  and leaf-type tests.

**Solidly verified.** The batch layer of L1 is genuine content:
`accumulateTerms` is a fold over `Multiset (K × C)` with a `LeftCommutative`
instance, and `accumulateTerms_perm`/`_add`/`_singleton`/`batchMap` really do license
reordering, hash-partitioning and the scalar sugar;
`pushforward_eq_reset_accumulate` plus the negative
`merge_without_reset_ne_pushforward` are real, as is `accumulate_ne_overwrite`
(which adjudicates AGAINST legacy `HashMap::extend`). L3 is well done: `overlap`
biadditivity/symmetry/homogeneity, the sesquilinear `hermitianOverlap` block over a
`StarOrderedRing`, and `clifford_conjugation_preserves_overlap`. `Retain` is the
strongest section (`retain_seq_eq_retain_and`, `retain_comm`, `retain_key_add`,
`batchMap_filter_key`, plus the `truncMag_not_additive` witness). `Twisted.lean` is
model-citizen quality. On the Rust side the `Coefficient`/`ImaginaryUnit`/
`Conjugate`/`Phase`/`KeyProduct` laws are pinned by real tests rather than assumed:
`phase1_leaf_types.rs` checks `i·i = −1`, `i⁴ = 1`, `i² ≠ 1`, conjugation
involution/additivity/multiplicativity, the full ℤ/4 `Phase` group table against
both the packed exponent and the concrete {1,i,−1,−i} values, `magnitude`
subadditivity and multiplicativity, and re-derives `key_mul` against a hand-built
2×2 ℤ[i] model plus the 2-cocycle identity; `sym_lean.rs` re-runs the
`ImaginaryUnit`/`Conjugate` laws exhaustively over a `GaussianInt` box.

### sum-engine-stores

- Rust surveyed: `ppvm-pauli-sum-2/src/{store,sum,ops,policy}.rs` (all of
  `store.rs`'s 1791 lines), `column_store/{mod,columns,lifecycle,graded,gates}.rs`,
  `column_store/rotations/{mod,rx}.rs`,
  `indexmap_store/{mod,lifecycle,algebra,gates,branching}.rs`,
  `ppvm-traits-2/src/{batch,graded}.rs`, `ppvm-pauli-sum-2/tests/engine.rs`.
- Lean: `PPVM/Algebra/GradedMap.lean` (all 725 lines),
  `PPVM/Algebra/Truncation.lean`, `PPVM/Instantiations/Rotation.lean`
  (`twoPass` family).
- Oracles: `pauli_sum_lean.rs`, `column_store_lean.rs`, `column_store_diff.rs`,
  `pauli_sum_indexmap_diff.rs`.

**Solidly verified.** The abstract graded algebra. `GradedMap.lean` grounds L0–L4 in
genuine Mathlib structures rather than in what the code does, and earns real
content: order/partition invariance of `accumulate_batch`, the `reset` obligation for
`apply_producer` with its negative witness, `accumulate_ne_overwrite`, the whole
`Retain` block for `CombinedPolicy` and the `usize::MAX` sentinel,
`truncate_preserve_eq_widened_retain` with both preserve guards pinned,
`clifford_conjugation_preserves_overlap`, and the ℓ¹/ℓ² bounds with the
`l1_bound_needs_subadditive` / `l1_bound_seminorm_needs_zero` counterexamples.
`pauli_sum_lean.rs` and `column_store_lean.rs` bridge most of those to the real Rust
on two backends, and `column_store_lean.rs::assert_aligned` is a genuinely good
executable invariant over the SoA columns.

### truncation-policy-loss

- Rust surveyed: `ppvm-pauli-sum-2/src/{policy,loss}.rs`, `sum.rs:437-700`,
  `rotation.rs:143-260`, `lib.rs` (aliases), `column_store/graded.rs`,
  `column_store/columns.rs`, `ppvm-traits-2/src/{batch,graded}.rs`,
  `ppvm-lossy-pauli-word-2/src/data.rs:288`, `ppvm-pauli-word-2/src/data.rs:307`,
  `ppvm-tableau-2/src/mixture/data.rs:142-260`, `ppvm-tableau-2/src/data.rs:747-1400`,
  `ppvm-tableau-2/src/expectation.rs:35-80`.
- Lean: `PPVM/Algebra/Truncation.lean`, `PPVM/Algebra/GradedMap.lean`
  (retain/batch), `PPVM/Algebra/Noise.lean:63-143`, `PPVM/Pauli.lean:70-73`.
- Oracles: `pauli_sum_lean.rs` (`truncation_l1_bound`,
  `truncation_cutoff_mismatch`), `column_store_lean.rs`; everything else is a
  legacy diff (`pauli_sum_truncation_boundary_diff.rs`,
  `pauli_sum_truncation_behaviour_diff.rs`, `pauli_sum_preserve_diff.rs`,
  `pauli_sum_loss_diff.rs`, `pauli_sum_integration_diff.rs`,
  `pauli_sum_workload_diff.rs`).

**Solidly verified.** What the Rust implements is narrower than expected and, usefully,
entirely predicate-shaped: there is **no top-k / max-terms / sorted-prune strategy
anywhere in the `-2` crates**. The five policies all factor through the single
`Retain::retain(|k,c| bool)` capability and truncation is strictly caller-driven.
That shape *is* faithfully modelled: `GradedMap.retain` matches the Rust predicate
signature, `retain_seq_eq_retain_and`/`retain_comm` capture `CombinedPolicy`,
`retain_of_all_true`/`retain_weight_le_eq_self` justify the `usize::MAX` sentinel,
and `truncate_preserve_eq_widened_retain` is a genuinely good, independently
motivated theorem that collapses the three-step snapshot/policy/restore composite to
one widened retain and thereby pins both otherwise-invisible guards.
`l1_bound_abv`/`l1_bound_seminorm` plus their two counterexamples are real
mathematics about what law `Coefficient::magnitude` must satisfy, and they are used
to adjudicate `ppvm-sym-2`'s `magnitude` against the law rather than against legacy.
`cutoff_abs_iff_sq` correctly collapses the tableau's three cutoff spellings to two.

### multiply-rotation

- Rust surveyed: `ppvm-pauli-sum-2/src/{rotation,multiply,producer}.rs`,
  `column_store/rotations/{mod,rx}.rs`, `column_store/graded.rs`,
  `store.rs:348-1495`, `indexmap_store/{algebra,branching}.rs`,
  `ppvm-traits-2/src/containers/{hash_join,coordinate_list}.rs`,
  `ppvm-traits-2/src/coefficient.rs:153-207` (`Angle::sin_cos`).
- Lean: `PPVM/Algebra/Twisted.lean` (all 408 lines),
  `PPVM/Instantiations/Rotation.lean` (all 679 lines), `PPVM/Pauli/Matrix.lean`.
- Oracles: `pauli_sum_multiply_lean.rs` (648 lines),
  `pauli_sum_rotation_noise_lean.rs` (384 lines), `column_store_lean.rs` (one `rx`
  alignment probe); rotation is otherwise covered by `gate_surface.rs`,
  `column_store_diff.rs` and `pauli_sum_integration_diff.rs` — all legacy-derived.

**Solidly verified.** The L4 operator product. `Twisted.lean` grounds it
independently (Mathlib `Finsupp`/`CommRing`, the `phaseExp` 2-cocycle, `iPow_add`
from `i⁴=1` alone), proves the whole-map laws (`twistedConv_add_left`/`_right`,
`twistedConv_assoc`, `twistedConv_apply_id`, `twistedConv_single_right_apply`
spending `mulWord_right_injective`), and abstracts the obligation every `KeyProduct`
owes. `pauli_sum_multiply_lean.rs` bridges it against a genuine dense 2ⁿ×2ⁿ ℂ-matrix
oracle (`mat(A·B) == mat(A)·mat(B)`), exhaustive over all 16×16 single-qubit support
pairs and all 64 monomial triples, randomized to n=4, plus the trace tie and a
store-buffer interleaving test — and it deliberately asserts the *Lean* value where
legacy is wrong (old's non-bilinear `MulAssign` chain). No fidelity drift was found
in `multiply.rs`, `store.rs:1040/1074`, `column_store/graded.rs:158` or the two
container impls: all four are faithful transcriptions of `twistedConv`. Also worth
recording: `Angle::sin_cos` is an unconditional `f64::sin_cos`, so there are no
special-angle fast paths to verify, and the commuting case is a genuine no-op
(early `return None`, no zero term, no truncation).

### noise-observables

- Rust surveyed: `ppvm-pauli-sum-2/src/{noise,loss,trace,store}.rs`,
  `ppvm-tableau-2/src/noise.rs`, `ppvm-tableau-2/src/mixture/noise/{pauli,loss}.rs`,
  `ppvm-tableau-2/src/mixture/{equality,data}.rs`,
  `ppvm-traits-2/src/{coefficient,gates}.rs`.
- Lean: `PPVM/Algebra/Noise.lean` (all 20 declarations), `PPVM/Pauli/Matrix.lean`,
  `PPVM/Pauli/Symplectic.lean` (`omega`), `PPVM/Algebra/{Truncation,GradedMap}.lean`.
- Oracles: `pauli_sum_rotation_noise_lean.rs` (single-qubit λ_P only); everything
  else is a legacy diff (`pauli_sum_gate_surface_diff.rs`, `pauli_sum_loss_diff.rs`,
  `tableau_mixture_diff.rs`, `tableau_behaviour_diff.rs`).

**Solidly verified.** (1) The arithmetic collapse
λ_P = Σ_Q p_Q(−1)^{ω(P,Q)} = 1 − 2Σ_anti p_Q over an arbitrary finite key set and
tied to `Symplectic.omega`, with contractivity under sub-stochasticity,
`l1_contractive`, and `scaleByKey_support_subset` — and the single-qubit case IS
pinned to the real `PauliSum<f64>::pauli_error` exhaustively over P at n=1 with
randomized mixed p (`pauli_sum_rotation_noise_lean.rs:337`). (2) The fifteen
hand-written anticommuting index lists of `two_qubit_pauli_error` are genuinely
re-derived in Lean by `decide` (exactly-8, nodup, set-equal to the ω-anticommuting
set); all 15 arms of `noise.rs:178-192` were hand-checked against `oldIndices`/
`qPair` and the transcription is faithful despite the different site encodings.
(3) Heisenberg trace preservation Λ*(I+L)=I+L for the pauli-sum loss/reset/
correlated-loss transfer matrices — real independent mathematics that refutes the
crate's own "suspected old bug 4".

### hashing-digests

- Rust surveyed: `ppvm-traits-2/src/{hash,batch}.rs`,
  `ppvm-pauli-word-2/src/{hash,storage,data,column,product,clifford}.rs`,
  `ppvm-lossy-pauli-word-2/src/{hash,data,column}.rs`,
  `ppvm-pauli-sum-2/src/store.rs` and `column_store/{columns,graded}.rs`,
  `ppvm-tableau-2/src/data.rs` (`Tableau::key_hash`),
  `ppvm-tableau-2/src/mixture/{fingerprint,data,equality,noise/*}.rs`,
  `ppvm-sym-2/src/term.rs`.
- Lean: none. All 19 `.lean` files were grepped for `key_hash`, `Indexable`,
  `digest`, `IdentityHasher`: no hits outside design-doc URLs.
- Oracles: none. Coverage is hand-written contract tests (`pauli_sum_hash.rs`,
  `tableau_hash.rs`, `lossy_pauli_word_diff.rs:365-420`,
  `phase1_leaf_types.rs:1288-1326`, `phase1_batch.rs:81`) citing the design doc.

**Solidly verified.** Not by Lean, but the contract tests are good and worth keeping:
(1) `Hash for K` writes exactly one `write_u64(key_hash())` — pinned generically for
an arbitrary `Indexable` and for `PauliWord`/`LossyPauliWord`/`Tableau`, the last
with a recording hasher that also rejects extra `write` calls; (2) structurally
equal keys have equal digests over seeded random populations at n=1..16 and after
clone/mutate/replay; (3) `KeyBatch::fill_hashes` and `KeyColumn::hash_into`
reproduce `key_hash()` (thin: 4 fixed words per column type, one storage width);
(4) the low-bit avalanche and the `HashFinalize` storage-tier fold, with a genuinely
discriminating detector at every `[u8; N]` tier on a low-weight population.

### tableau-core

- Rust surveyed: `ppvm-tableau-2/src/clifford.rs` (all 1100 lines),
  `ppvm-tableau-2/src/data.rs` (Row/`mul_assign`, frame construction,
  `find_z_anticommuting_stabilizer`, `get_deterministic_outcome`,
  `update_tableau_according_to_outcome`, `compute_decomposition`, the `cz_block`
  family), `ppvm-tableau-2/src/gates.rs`, `ppvm-traits-2/src/{gates,pauli}.rs`,
  `ppvm-phased-pauli-word-2/src/clifford.rs` (convention cross-check),
  `ppvm-tableau-2/src/mixture/data.rs`.
- Lean: `PPVM/Tableau/Frame.lean`, `PPVM/Tableau/Batch.lean`,
  `PPVM/Pauli/{Symplectic,Conjugation,Phase,Word,Matrix}.lean`,
  `PPVM/Tableau/{BranchPhase,Projection}.lean` (indexes),
  `PPVM/Instantiations/Rotation.lean` (RotXY).
- Oracles: `tableau_lean.rs` (read in full); `tableau_diff.rs`,
  `tableau_behaviour_diff.rs`, `ppvm-tableau-2/tests/behaviour.rs` are legacy/batch
  agreement.

**Solidly verified.** (a) The tableau stores 2n rows — destabilizers in `0..n`,
stabilizers in `n..2n` — each with a ℤ/4 phase column, and the bit-level
symplectic-basis invariant is genuinely proven and genuinely bridged:
`IsSymplecticFrame` + `frame_linearIndependent` + `frame_surjective` +
`frame_coordinate_expansion` (a real Yoder-Lemma-5 spanning/coordinate result, not a
restatement of the code) are pinned by `identity_frame_is_symplectic`,
`every_gate_preserves_the_symplectic_frame` (n=70, 300 mixed steps, all 18 gate entry
points including the batch and `cz_block` kernels, asserted after EVERY step — so the
inductive-across-a-sequence question is answered empirically even though Lean has no
circuit fold), `frame_rows_are_linearly_independent`,
`measurement_preserves_the_symplectic_frame`, and
`msd_workload_preserves_the_symplectic_frame`. (b) The projection's bit action really
is `projectFrame` (`rowUpdate_eq_ite` matches the Rust's guarded `mul_assign(&g_q)`
and the pivot/overwrite pattern). (c) The g-rule in `Row::mul_assign` is grounded in
real ℤ[i] matrices via `phaseExp`→`phaseRef`→`pauliMat_mul` (though not bridged for
this crate — G-058). (d) `h`/`s`/`x`/`y`/`z`/`cnot`/`cz` row actions are pinned
bit-and-ℤ/4-exactly to the corresponding `conj*`. (e) `cy`'s bits and sign match
`conjCY` exactly, and `conjCY` is derived as `conjST∘conjCNOT∘conjSdagT`.

### measurement-branching

- Rust surveyed: `ppvm-tableau-2/src/measure.rs` (all 771 lines),
  `ppvm-tableau-2/src/data.rs:89-149`, `:299-320`, `:454-549`, `:1000-1167`,
  `ppvm-tableau-2/src/mixture/{measure,sampler,data,equality,gates,fingerprint}.rs`,
  `ppvm-pauli-sum-2/src/proj.rs` (cross-check).
- Lean: `PPVM/Tableau/Projection.lean`, `PPVM/Tableau/BranchPhase.lean`,
  `PPVM/Tableau/Frame.lean`, `PPVM/Instantiations/{Bitstring,Projector}.lean`,
  `PPVM/Pauli/{Phase,Word}.lean`.
- Oracles: `tableau_lean.rs` (996 lines; frame-level only — grep for
  `overlap`/`prob_1`/`probOne`/`projectRaw`/`Born`/`mixture`/`sampler` across all
  nine `*_lean.rs` returns zero hits).

**Solidly verified.** (a) The *shape* of the measurement case split —
`measurement_dichotomy` / `measure_deterministic_iff_xfree` with a real,
mutation-sensitive oracle (`measurement_dichotomy_holds`, which also checks
idempotence and that Z_q becomes a stabilizer). (b) That the one non-unitary frame
mutation keeps the 2n rows a symplectic basis (`isSymplecticFrame_projectFrame`),
bridged by `measurement_preserves_the_symplectic_frame` and the 85-qubit MSD sweep.
(c) `compute_decomposition`'s two bitmasks really are frame coordinates
(`frame_coordinate_expansion` + `frame_surjective`, done by counting, no
hand-waving). (d) The case-a/case-b ℤ/4 *model* arithmetic — `rustTerm_eq`,
`overlap_eq_inner`, `proj_add`/`proj_idem`, `probOne_eq`,
`proj_zero_eq_caseB_retain` — with the abstract `SelfInverse` hypothesis genuinely
discharged for the crate's own phase function by `BranchPhase.lean`
(`frameOp_eq_shiftOp`, `selfInverse_branchPhase_iff`), and that discharge pinned to
the real code by `decomposition_satisfies_the_lean_frame_identity`.
`BranchPhase.lean` in particular derives the phase formula rather than declaring it.

### symbolic-coefficients

- Rust surveyed: `ppvm-sym-2/src/{lib,term,mul,add,eval,coeff,exact}.rs`,
  `ppvm-sym-2/tests/exact_ring.rs` (legacy `ppvm-sym/src/mul.rs` for the patched
  cross-loop gate only).
- Lean: `PPVM/Instantiations/Symbolic.lean` (all 951 lines),
  `PPVM/Algebra/GradedMap.lean`, `PPVM/Pauli/Matrix.lean` (ℤ[i]),
  `PPVM/Algebra/Truncation.lean`, `PPVM/Instantiations/Projector.lean`.
- Oracles: `sym_lean.rs` (all 897 lines), `ppvm-sym-2/tests/exact_ring.rs`;
  `sym_diff.rs` is legacy agreement.

**Solidly verified.** The Lean file is genuine independent mathematics, not laundered
legacy: `SymRing = AddMonoidAlgebra ℝ (ℕ →₀ ℕ×ℕ)` and
`PhasedSymRing = AddMonoidAlgebra ℝ (Mono × ZMod 4)` are real Mathlib structures, and
`substHom`/`evalHom`/`evalC` are real `AlgHom`s built via
`AddMonoidAlgebra.lift`. The hardest facts are proved rather than assumed:
sine-degree additivity and the truncation *ideal* (`sinDeg_add`,
`truncIdeal_mul_right`), drop-at-insert = drop-at-end for `max_sin`, the `clear()`
soundness for the degree arm, the honest *negative* results
(`pythagorean_ne_one`, `symRot_norm_sq_ne_symbolically`, `evalC_not_injective`,
`epsClear_ne_retain_pointwise`, and especially `mulImpl_not_wellDefined`, which
refuses to over-claim `set_max_sin`), the ℓ¹ over-truncation bound for `mul_term`'s
eps arm, `phaseFold_eq_iSym_pow_mul` (which genuinely *forces* divergence #3 rather
than rationalising it), and `conjSym`/`evalC_conjSym`. `exact.rs`'s ℤ[i] is grounded
in Mathlib's `GaussianInt` in `Pauli/Matrix.lean` and pinned exhaustively and
exactly in `sym_lean.rs`.

## Out of scope

**Tier C crates — recorded, not omitted.** No file was read in `ppvm-cli`,
`ppvm-tui`, `stim-parser`, `ppvm-stim`, `ppvm-python-native`, `ppvm-vihaco`, or
`vihaco-circuit-isa`. They consume the `-2` crates (e.g.
`ppvm-vihaco/src/observable.rs` and `ppvm-python-native/src/interface.rs` use
`LossyPauliWord`) and they own the RNG on the user's behalf, but no mathematical
claim in this ledger depends on them. If a later round widens scope, the natural
first questions are which `Policy` each frontend configures and whether the
`PPVM_EXPECT_BACKEND` matrix would notice a `-2`-only defect.

**Legacy crates.** `ppvm-pauli-word`, `ppvm-pauli-sum`, `ppvm-tableau`,
`ppvm-tableau-sum`, `ppvm-sym`, `ppvm-traits` were read only to adjudicate Lean
citations and provenance. They are never verification targets.

**Auditor notes recorded here rather than filed as gaps.**

- *Not one of the five classes (Rust-internal doc/code drift).*
  `ppvm-pauli-word-2/src/hash.rs:73` still documents `Indexable::key_hash` as
  "lazily computed once and cached" though `3e43610f` made the digest an eager plain
  `u64` field; `5549c1f3` fixed the surrounding docs but missed this line.
  `ppvm-traits-2/tests/phase1_leaf_types.rs:996` claims "the crate ships no L4 impl
  (`containers.rs` stops at L3 + `Retain`)" — both halves are false since
  `d92421b6`. `ppvm-pauli-sum-2/src/lib.rs:51` says `Pair::probe_batch` "consumes
  the batch's precomputed digest column", but every shipped `probe_batch` outside the
  columnar one ignores `KeyBatch::hashes()`. `ppvm-traits-2/src/clifford.rs:26-28`
  claims the tableau's fused bodies mean "the *values* are the blanket's", which is
  false for sqrt_y/sqrt_y_dag (the blanket's h;z and z;h yield the adjoint sign);
  it holds for s_dag and sqrt_x/sqrt_x_dag only because S,Z commute and h;s;h is
  palindromic.
- *Statistical, not provable.* The avalanche contract's top-7 control-tag half
  (`ppvm-traits-2/src/hash.rs:14-16`, `ppvm-pauli-word-2/src/storage.rs:73-79`) is
  measured only in the low 12 bits by every test, and the `finalize_hash` width
  branch is unpinned in the direction it claims to protect. Hash avalanche quality is
  a statistical property of a concrete finalizer with no Lean statement available; see
  the appendix.
- *Design-doc edits.* `word-data-structures.md:452-454` asks for tests of two
  properties the shipped representations do not have (component digests on the lossy
  word; unused-high-bit identity tests that do not exist). `lean/README.md`'s "No
  open targets remain" is scoped to design-doc algebra claims and was never scoped to
  the loss, hashing, mixture, or trajectory-sampler surfaces.
- *Robustness, not formalization.* `GeneralizedTableauMixture::normalize_probabilities`
  (`mixture/data.rs:135-140`) divides by the weight sum with no zero guard, so a
  fully decayed mixture yields NaN weights and a NaN cumulative vector in the
  sampler. `ppvm-tableau-2`'s `rotate_2` with an identity axis reaches a
  `debug_assert_ne!(pauli, Pauli::I)` in `compute_decomposition` (panics in debug,
  computes the correct global-phase-only result in release).
  `ppvm-tableau-2/src/noise.rs:306 asymmetric_loss_channel` reproduces a legacy
  measurement-record contract violation (it does not pop the record entry its internal
  measure pushed, shifting every downstream `rec[-k]`); already in the crate's
  Deferrals.
- *Reproducibility artifacts, not claims.* `burn_legacy_tableau_seeds` in
  `mixture/noise/{pauli,loss}.rs` and `mixture/data.rs:105-117` /
  `sampler.rs:44-47` exist only to keep the RNG stream byte-compatible with legacy.
  The shot stream is therefore pinned to legacy's draw order; no Lean statement
  covers or needs that.
- *Deliberate, well-argued divergences that need no theorem.* `NoPolicy`'s
  `capacity()` clamp (`policy.rs:59-98`) versus legacy's unclamped `1 << (2n-1)` is
  allocation sizing with no observable semantics.
  `ppvm-pauli-word-2/src/data.rs:448-463 set_z_bit_pair` diverges from
  `into_toggled_bits2` at i == j (last-write-wins vs toggle-cancel) to match the
  trait default's two-scalar-set semantics; well covered by in-crate unit tests and
  no algebraic claim is credited to it. `Phased<LossyPauliWord>` deliberately gets
  Clifford but no product (`ppvm-phased-pauli-word-2/src/product.rs:20-23`), and
  `LossyPauliWord` implements no `KeyProduct` anywhere, so "is the lossy product
  associative" has no Rust referent to audit.
- *Cross-sector referrals recorded by the auditors.* `Noise.lean`'s loss transfer
  matrices are transcribed arm-for-arm from legacy `ppvm-pauli-sum/src/sum/noise.rs`
  (trace preservation is necessary but far from sufficient — a permutation of branch
  weights preserving column sums would satisfy every theorem there); that is what
  G-034 and G-043 exist for. `ppvm-sym-2`'s `Prod` canonicality (sorted, deduped, no
  all-zero factor, `sin_pow`/`cos_pow` matching the summed exponents) is enforced only
  by `debug_check` under `cfg(debug_assertions)` while Lean models the monomial as a
  `Finsupp` that cannot express a non-canonical representation — the same
  representation-refinement hole as G-047, one crate over.
  `ppvm-pauli-sum-2/src/proj.rs` is a *positive* worth recording: Lean adjudicated
  legacy's `c ↦ c²/2` and `_ => None` as wrong (`oldStep_eq_half_iff`,
  `oldProj_not_idem`, `twoProj_conj_X/Y` over honest ℤ[i] matrices) and the shipped
  `-2` code now deliberately diverges from legacy on Lean's authority; an oracle test
  pinning `proj.rs` to `projLin` still appears to be absent.

## Appendix: refuted candidates

Sixteen candidate gaps were raised and then refuted during adversarial
verification. They are recorded so later rounds do not re-litigate them. Do not
re-open a row here without new evidence that contradicts the stated reason.

| # | Candidate | Why refuted |
| --- | --- | --- |
| R-01 | `CliffordExtensions` declares one conjugation table; tableau-2 implements its adjoint (`ppvm-traits-2/src/gates.rs:95`) | The direction-polymorphic decomposition is already in Lean: `Conjugation.lean` defines BOTH orientations of every composite as generator products (`extSqrtX := conjH∘conjSdag∘conjH` at `:294`, `extSqrtXdag` at `:297`, `extSdag_eq_conjS` at `:332`) with the dagger-pair theorems at `:419-423`, and each implementor's oracle pins it to its own correct map (`tableau_lean.rs:395` `s`↦`conj_s` = x∧z; `phased_pauli_word_lean.rs:496` `s`↦`conjSdag` = x∧¬z). Residue is a Rust docstring scope nit. |
| R-02 | Extension `_sign`/`_bits` and `cnotDelta`/`czDelta` cite legacy files, not the `-2` kernels | Contradicted by the Lean text: `Conjugation.lean:354-366` says the legacy kernels ship the predicates "with no derivation" and that these are "now *derived* from the generator product rather than asserted"; `:479-490` justifies conjCNOT/conjCZ by MonoidHom plus the four generators, i.e. genuine `G·P·G†`. All six `-2` predicates were re-checked and agree, so the fidelity rule (legacy citation over *differing* code) does not apply. |
| R-03 | Lossy CY / extension gates pinned only by the legacy diff (`ppvm-lossy-pauli-word-2/src/clifford.rs:83`) | Duplicate of G-001 — same Lean symbol, same defect, same closure, differing only in the Rust anchor (the `xor_z_from_x` guard vs the `BlanketClifford` opt-in). |
| R-04 | `LossyPauliKeyColumn`'s per-qubit bit addressing has no test anywhere | False premise. `ppvm-pauli-sum-2/tests/column_store_gates.rs:131` declares `type ColumnLossy = Sum<ColumnStore<LossyPauliWord, f64>, NoPolicy>` and `lossy_channels_match_on_columnar_lossy_words` runs the loss channels against the HashMapStore `LossyPauliSum` term by term — an SoA-vs-AoS parity check (not a legacy diff) exercising the `plane_bit`/`is_lost` addressing. Coverage is thin (n=2, channels only) but the claim as written does not survive. |
| R-05 | `GradedMap`'s L4 section models the untwisted `AddMonoidAlgebra` | Misreads both sides. `GradedMap.lean:527-530` itself says `multiply_single` is "the phase-free core … (a 2-cocycle *twist* of this untwisted convolution)"; `graded.rs:190-197` cites `tmul_assoc` and `twistedConv` for the twisted laws in the same paragraph. The abstract-key obligation is stated at `Twisted.lean:87-131` (`gtmul`, `IsCocycle`, `gtmul_assoc`), and the only shipped `KeyProduct` impl is `PauliWord`. |
| R-06 | hash_join vs coordinate_list agreement is stated nowhere and pinned by one fixture | Not the repo's verification unit: both refinements are individually pinned to the one Lean model (`containers/tests.rs` and `phase1_containers.rs` carry the paired accumulate/reduce/scale/retain/`for_each_ref` tests), so observational agreement is a consequence up to float reassociation. The concrete holes it lists are already G-006, G-048 and G-076. |
| R-07 | `Vec<(K,C)>` admits duplicate-key states where its own `get` and `overlap` disagree | The arithmetic disagreement is real but the state is unreachable through the shipped API: `accumulate_batch`, `multiply_into`, `add_term` and `insert_term` all scan-then-update, and `RekeyBijective` is reachable only from `pub(crate) Sum::rekey_bijective` under proved-injective Clifford maps. What remains is API hygiene on a foreign type. |
| R-08 | Commutativity of the coefficient ring is a hidden hypothesis of every L2/L3/L4 theorem | The "discharged nowhere" half is false: `phase1_leaf_types.rs:593` already asserts `a*b == b*a` exhaustively over the ℤ[i] box with the comment "commutative ring (required by L4/tmul_assoc)", `sym_lean.rs:191-209` does the same for the exact `Term` ring, and `f64`/`Complex<f64>` are commutative by construction. What is left is an API-documentation obligation for hypothetical future impls. |
| R-09 | Batch order/partition invariance over `AddCommMonoid` vs the f64 backend (`store.rs:1551`) | Duplicate of G-076, which carries the same Lean citation and closure plus the support-discontinuity argument; the merge-direction choice at `store.rs:1551` is one instance of the license, not a second gap. |
| R-10 | `Policy::is_noop` short-circuits the preserve pipeline with no Lean obligation | Over-stated on both legs. `retain_of_all_true` is stated pointwise (so it covers stored zeros), `retain_le_top_eq_self` is exactly the `usize::MAX` sentinel, and the skipped composite follows in one step from `truncate_preserve_eq_widened_retain` plus `retain_of_all_true` on the widened rule. Every shipped `is_noop` answer is correct by inspection (`policy.rs:101, :149, :306, :249-251`). The oracle request duplicates G-069. |
| R-11 | "Lost endpoint ⇒ no noise at all" is punted by Lean but is not the marginal channel | The load-bearing assertion (that the correct restriction is the marginalized single-qubit channel) is a physics-modeling opinion, not a claim the repo makes: `noise.rs:194-197` states the opposite choice explicitly ("if just one atom is lost, then there is no well-defined noise channel on the other atom"). No theorem or docstring is credited with marginalization, so none of the five classes applies; the proposed closure is a semantics change needing sign-off. |
| R-12 | Top-7 control-tag half of the avalanche contract is asserted by no test | The observations hold (every measurement is low-bits-only) but there is no proof obligation to carry: avalanche quality and bit independence are statistical properties of a concrete finalizer, and the candidate's own closure concedes Lean cannot state it. Test-quality and docstring matter, not a gap class. |
| R-13 | Cited rationale "loss-only mutation must not rehash the X/Z planes" is false as shipped | Docs-consistency nit with no mathematical content: the module docstring the candidate cites (`ppvm-lossy-pauli-word-2/src/hash.rs:4-15`) already documents the shipped single-fold design and its real measured rationale. The stale request lives in `word-data-structures.md`, i.e. a design-doc edit. |
| R-14 | Lean's `ext*` maps are the ADJOINTS of the tableau extension gates they claim to model | The mathematical content is correct; only two sentences of prose are inverted. `Conjugation.lean`'s ext* section targets `ppvm-traits-2/src/pauli.rs`'s blanket (whose `s` is the phased word's backward `conjSdag`), and `Batch.lean:150-159` supplies the correct witness per tableau gate — each predicate matches the Rust bit for bit (sqrt_x↦`extSqrtXdag` z∧¬x = `clifford.rs:468`; sqrt_y↦`extSqrtYdag` x∧¬z = `:507`; …). HZ = R_y(π/2) was checked independently. The mutation risk it describes is G-059. |
| R-15 | Frame projection is proved invariant-preserving, not the correct post-measurement stabilizer group | Duplicate of G-065 — same Rust function and sign site (`data.rs:509-548`), same Lean symbol (`Frame.lean:280/:296`), same missing artefact (phase-carrying frame semantics), same proposed dense n ≤ 4 oracle. G-065 states it with the stronger `canonicalize` citation. |
| R-16 | `mulImpl` takes one global bound; the Rust inherits max_sin/min_eps from the LHS, so its product is non-commutative | Premise about the Lean is wrong: `mulImpl` is parametric in k, and both mixed Rust arms take the bound from the left operand (`mul.rs:271`, `:302`), so modelling both as `truncMulMono k …` is faithful call by call. The file already carries the stronger negative `mulImpl_not_wellDefined` and the Rust-side witness `truncation_parameters_come_from_the_left`. The min_eps half is G-055. |

## Round log

**Round 1 (audit) — 2026-08-10, branch `codex/traits-2-impl`, base commit
`d92421b6`.** 12 sectors surveyed end to end plus three cross-cutting skeptic
sweeps; 83 gaps opened (G-001…G-083), 16 candidates refuted (R-01…R-16). No gap
closed. Open by class: bridge 26, coverage 25, strength 17, fidelity 12,
provenance 3. Open by severity: high 27, medium 45, low 11. Sectors with zero Lean
coverage: `hashing-digests`. Sub-areas with zero Lean coverage: the whole
`ppvm-tableau-2/src/mixture/**` subtree and the tableau trajectory noise samplers.

Highest-leverage clusters for the first close round, by shared root:

1. **Measurement signs** — G-021, G-017, G-023, G-062, G-065, G-058 all fail
   because no Lean object carries the tableau's ℤ/4 phase. One phase-carrying
   `Frame` plus a full-ℤ/4 oracle retires six rows, three of them high.
2. **Rotation direction** — G-030 is the root of G-025, G-029, G-032, G-031 and
   G-050. One derived conjugation lemma plus one dense-matrix oracle retires six.
3. **Truncation parameters** — G-074 and G-067 are the same missing statement seen
   from Lean and from the oracle; G-075 extends it to the multi-step case.
4. **Packed-storage invariants** — G-078 and G-009 are one invariant; G-079, G-080
   and G-077 are its immediate neighbours in the same file.
5. **Mixture** — G-018, G-019, G-071, G-037, G-020 and G-012 are one absent Lean
   module (`PPVM/Tableau/Mixture.lean`) and one absent oracle file
   (`tableau_mixture_lean.rs`).

Later rounds append below.

**Round 2 (adjudication): six defect clusters examined, 6 corroborated by
independent re-derivation, 0 contested, 9 live defects found.** Each cluster's
answer was re-derived from the algebra twice, by two agents working independently
with their own probes; every second opinion returned `agree-with-corrections`, so
no verdict is contested *between the two agents*. 27 rows amended:
7 `adjudicated-defect` (G-035, G-018, G-019, G-020, G-061, G-060,
G-054), 18 `adjudicated-spec`. Seven `Proposed closure` entries were editing
targets because they would have proved the wrong statement — G-040 (would have
legitimized the 2p₁ convention), G-043 (`p₀ + 2p₁ ≤ 1`), G-008 (its own weakening
is false), G-018 ("or the sub-stochastic mass"), G-057 ("associative at eps = 0"),
G-071 (insertion-time bypass as a defect, and a loose bound), G-061 (XOR dedup /
"the documented behaviour"). Three ledger premises were **refuted**: G-013's
mutation does not survive the suite, G-060's fused block corrupts the bit planes
as well as the sign, and G-040's "structurally different one-already-lost arms"
are the Heisenberg transpose of each other. The three units whose code is right
(U2 rotation direction, U3 lossy canonicality, U5 measurement sign) each produced
prose/Lean corrections instead, and U5 is now positively verified against a dense
simulator (≈40k deterministic measurements, zero mismatches) rather than merely
unfalsified. Nothing was applied to any crate; see the
[sign-off list](#adjudications-round-2) per unit.

**Round 2 correction + ruling (post-round).** The round ran without the
paper draft (`../ppvm-paper/main.tex`) in evidence. With it, **U1's fix direction
is reversed**: the paper (`:462`, `:523`, and the demo arithmetic at `:845`) defines $p_{LQ}$ as the
probability a *named* one of the pair is lost, so P(exactly one) $=2p_{LQ}$ —
the `ppvm-pauli-sum-2`/`corrT`/legacy reading — while the shipped Python API
(`mixins.py:507`) and its test suite (`test_loss.py:82`, `:173`) document the
opposite. Both are the same CPTP family up to a factor of two, which the unit's
own derivation established, so this is a user-facing convention decision, not a
defect. **Ruled 2026-08-10 in favour of the paper**: `ppvm-pauli-sum-2`, legacy and
`corrT` are correct and do not change; `ppvm-tableau-2`'s mixture
(`mixture/noise/loss.rs:95`) and trajectory (`noise.rs:365`) are the defect, along
with the `mixins.py` docstring and three Python tests (one of which,
`test_loss.py:82`, uses a `p` that is inadmissible under the ruling). Live defect
count stays at 9; only the direction moved. Everything else in
U1 stands, including G-035's missing guard and the refutation of `loss.rs:106-111`'s
"suspected old bug 4". From round 3 on, `../ppvm-paper/main.tex` is a required
source for every agent, cited alongside the Rust — it is the definitions of
record and it settles several rows the audit could only mark ambiguous.
