# Log for 2026-07-16-kossakowski

## 2026-07-16

Task: minimize the per-`pc_step` wall time of the Kossakowski dissipator
path (criterion `pc_step_superradiance/kossakowski_n30`, B=4096, A=3B,
free-space superradiance chain d=0.1λ0, model of fig11).

Baseline (feature as merged, commit at branch point):
  eigenmode   n10 274 ms   n20 1.433 s   n30 4.349 s
  kossakowski n10  81 ms   n20 0.161 s   n30 0.250 s
  → representation win 3.4× / 8.9× / 17.4× (~0.55·N, the predicted
  N-fold scaling; recorded for the PR note).

Phase profile of the n30 Kossakowski step (pc_step_timed, basis 4096):
leakage1 15%, expm1 20%, leakage2 45%, expm2 19% — every phase is one
action pass over the (enriched) basis, so the target is the
KossakowskiPair arm of compute_action_terms.

Process note: iterations are implemented directly on this working branch
(small surgical diffs, one commit per attempt) rather than via worktree
subagents; keep/discard is decided on the criterion metric and recorded
here either way.

## Final note (2026-07-16)

DELIVERABLE SUMMARY.

Exact equivalence of representations (uncapped, same dt):
- random models (mixed σ⁻ + 2-term lincomb ops, real AND complex
  Hermitian PSD K, with Hamiltonian): max |Δcoeff| < 1e-12 over 3 steps.
- N=6 superradiance chain (fig11 model), T=1 full basis: Kossakowski vs
  eigenmode max |ΔR| < 1e-11; vs exact excitation-cascade ED 1e-4 target
  met; vs the stored fig11 data.h5 n6/exact reference: 1.11e-05.

Orbit-rep dissipative status: WORKED OUT OF THE BOX — no Rust fix
needed. The suspected stabilizer/identity-word miscounting does not
exist: under the 1/|G| projection convention of
canonicalize_basis_arr_complex (c_rep = coeff_word/|Stab|), the
phase-aware action is exact for any equivariant generator including
non-unital flows into the identity orbit. Pinned by: N=6 ring orbit-vs-
full equivalence (max |Δc| ~2e-15, both representations), R(t) vs ED
<1e-4, closed-form identity bookkeeping (c_I = −(1−e^{−Γ0 t}) = coeff_I/
|G| to 1e-9), truncated N=10 sanity (finite, admission-bounded, monotone
in rep budget). Corollary worth remembering: every orbit contributes
|G|·c_rep to a translation-invariant {I,Z} observable regardless of
stabilizer, so R = |G|·Σ c_rep over {I,Z} reps.

Measured scaling (criterion pc_step_superradiance, B=4096, A=3B, per
step; chain d=0.1λ0). Quiet-machine session (pre-it3 code):
  N=10: eigenmode 274 ms vs kossakowski  81 ms  ( 3.4×)
  N=20:           1.43 s vs             161 ms  ( 8.9×)
  N=30:           4.35 s vs             250 ms  (17.4×)
Final loaded session (it3 code; absolute walls inflated ~1.9× by
user load 10-14, RATIOS internally consistent):
  N=10:  3.8×   N=20: 10.0×   N=30: 19.0×
i.e. the predicted ~N-fold representation win, growing linearly in N.
Optimization loop: it1 keep (precompiled sandwich table, −9.8%),
it2 discard (left-product reuse, noise-level), it3 keep (one-sided
±½[D,p] commutator fast path, −32% same-session A/B). All decisions
from same-session A/Bs after the walls shifted ~+28% under load —
never compare across sessions on a loaded machine.

Remaining headroom (not pursued): both-sided pairs still pay the full
12-product path and dominate for heavy strings; the all-to-all
Hamiltonian is ~40% of the step and is representation-independent.
