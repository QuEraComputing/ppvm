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
