# Log for 2026-07-17-drug-kossakowski

## 2026-07-17

Task: minimize the per-`pc_step` wall of the Kossakowski dissipator on the
**molecular dipolar-relaxation** workload (ZULF-NMR drug FID), which stresses
the path differently from the `2026-07-16-kossakowski` superradiance chain:

- operators are **2-local rank-2 tensors** with ~4 Pauli terms each (one
  channel per site pair × 5 spatial harmonics `m`), not single-site σ⁻;
- `K` is **block-diagonal** over the 5 components `m ∈ −2..2`, each block a
  dense `P×P` Gram matrix over the `P = C(N,2)` pairs;
- basis strings are dense/high-weight, so most candidate pairs take the
  **both-sided sandwich** (12-product) path — the arm the 07-16 ledger's
  final note flagged as remaining headroom ("both-sided pairs still pay the
  full 12-product path and dominate for heavy strings").

Harness: `crates/ppvm-lindblad/examples/drug_profile.rs` (a jittered N-spin
dipolar model with the block-diagonal rank-5 `Γ`, `pc_step_timed` phase
breakdown, Kossakowski path only — the eigenmode representation's O(N³) dense
blow-up is exactly what this path removes). Criterion bench
`benches/drug_dipolar.rs` carries the eigenmode-vs-Kossakowski comparison.
Metric: median `total/step`, `drug_profile 10 1500 5` (N=10, B=1500, 5 steps).
Same-session A/B for every keep/discard, per the 07-16 methodology.

Baseline (feature as of `4992d1e5`, i.e. the `2026-07-16-kossakowski` result):
  drug_profile n10 B1500  **7039 ms/step**
  phase profile: leak1 9%, expm1 30%, leak2 31%, expm2 30% — the dissipator
  action dominates *every* phase (expm passes are matrix-free = repeated `L*`),
  so the `KossakowskiPair` arm of `compute_action_terms` is ~100% of the step.

Iterations (each one commit; same-session A/B; kept unless noted):

- **it1 `2a24cd2e` — reuse `P_a·p` in the both-sided sandwich (keep).**
  Group the precompiled sandwich table by the left word so `P_a·p` is computed
  once per distinct `P_a` and reused across its `P_b` partners (32 → 20
  `pauli_mul` for 4-term ops). This is the 07-16 `it2` idea, which was
  noise-level *there* because σ⁻ ops are single-term (no `P_b` to reuse
  across); for the multi-term tensors here it is a real win.
  A/B: 7039 → **5765 ms/step = 18% (1.22×)**.

- **it2 `5a254d09` — fold conjugate pairs `(n,m)`+`(m,n)` (keep).**
  The action keeps only the real part, and `P_i·p·Q_j` / `Q_j·p·P_i` give the
  same word with conjugate phase, so a Hermitian-conjugate pair collapses to
  one upper-triangle (`n ≤ m`) entry using only the `(n,m)` products:
  both-sided sandwich doubles + takes Re; both-sided anticommutator uses
  `dd = 2·Re(K_nm A_n†A_m)`; one-sided uses `±½[F,p]` with the anti-Hermitian
  `F = −2i·Im(K_nm A_n†A_m)` (folded from the two conjugate one-sided
  contributions; sign derived and pinned by
  `kossakowski_pair_matches_eigenmode_jumps`, which passed first try).
  Halves the pair count (`P²` → `P(P+1)/2`) and drops complex arithmetic +
  per-pair imaginary-cancellation on the off-diagonal.
  A/B: 5765 → **2281 ms/step = 2.53×** (>2× because halving pairs also cuts
  candidate iteration and hashing). General for any Hermitian PSD `K`.

Cumulative baseline → it2: **7039 → 2281 ms/step = 3.08×**.

## Final note (2026-07-17)

DELIVERABLE SUMMARY.

Representation win vs the dense-jump ("eigenmode") form, measured on the real
**gemcitabine** 10-spin molecule (25 dense jumps ↔ block-diagonal rank-5 `Γ`,
B=2000, mean s/step over 8 steps, peak RSS; Python `pc_step_arr`):
  dense jumps       24.6 s/step   2239 MB
  Kossakowski (it2)  4.1 s/step   1720 MB
  → **6.0× faster, 1.3× less RAM** at N=10, and the gap grows with N (pair
  count `5·C(N,2)²` = O(N⁴), but the dense form additionally carries the
  ~2000-term jump `L†OL` intermediate that drives RSS). The eigenmode form at
  N=32 (fluticasone) does not complete a single step in the criterion budget —
  making the Kossakowski path the enabling representation for the 32-spin FID.

Correctness (unchanged across it1+it2): `kossakowski_pair_matches_eigenmode_jumps`
and the full lindblad suite green; 3-spin ZULF dense-vs-Kossakowski action
3.3e-15 and trajectory 3.2e-14 (also vs the stored fig10 `data.h5`);
gemcitabine (N=10) action 1.5e-15 over 90 random ≤5-local strings.
`cargo test --workspace` green; zero clippy warnings.

Remaining headroom (not pursued): the both-sided sandwich still rebuilds the
`local` hashmap per candidate pair; a support-bucketed candidate pass (group
pairs by shared `{a,b,c,d}` support) could cut redundant lookups for dense
strings. Left for a future campaign if the N=32 sweep needs it.
