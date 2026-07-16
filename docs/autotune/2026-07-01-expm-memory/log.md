# Log for 2026-07-01-expm-memory

Goal: reduce the optimized expm's memory (it wins/ties wall for N>=10 after the
cache-action + tol-matched work in [[../2026-07-01-expm-pc-step]], but uses
~2-6x Trotter's RSS, growing with N).

Build reminder (learned the hard way last experiment): rebuild the native from
`ppvm-python/`, NOT `crates/ppvm-python-native/`, or `import ppvm._core` loads a
STALE .so:
  `cd ppvm-python && VIRTUAL_ENV=.venv uvx --from maturin maturin develop --release --uv`

## Profile first: where does the expm RSS go? (N=10, drop=1e-4, dim~262k, RSS 396MB)

MEMPROF (env-gated eprintln, since removed):
- cols (CSC action cache, nnz*16B for (u32,f64) + dim*24B vec-of-vecs hdrs): **66 MB**
- index (Word->u32): ~14 MB
- partial_ys (threads x dim dense matvec accumulators): **20 MB (~5%)**
- coeffs: 2 MB
- remainder (~294 MB): python/numpy/ppvm baseline ~150 MB + pc_step transient
  (leakage-enriched basis, coeffs_predict, expm work vectors, leak HashMaps).

=> The matvec accumulator `partial_ys` I originally suspected is only ~5% of
RSS. The expm-specific memory is dominated by the action cache `cols`, and
everything (cols, index, coeffs, partial_ys, the basis itself) scales with
`dim` = basis size. The real memory lever is therefore whatever bounds `dim`:
`drop_tol`, `max_basis`, `dt`, and the accuracy target.

## Iteration 1: atomic-scatter matvec (remove partial_ys) — DISCARD

Replaced the threads x dim dense accumulators with atomic f64 CAS into one
shared output (O(dim) memory, f64 preserved -> tests bit-exact, 9/9 pass).
Measured N=10 drop=1e-4 (correct build): wall 19.4-20.0s vs 17.8s baseline
(~12% SLOWER), RSS 363-393MB vs 396MB (negligible, since partial_ys was ~5%).
Bad trade (CAS-per-nnz cost > the small memory saved). Reverted; kept only a
code comment recording the negative result. No other speed-neutral code fix
found: `cols` (nnz*16) is inherent to the cache that gave the ~10x speedup, and
f32 coeffs would halve it but break the drop_tol=0 exact-reference tests.

## The real lever: knob characterization (N=10 ladder, T=2, exact-ED ref)

(1) max_basis sweep @ dt=0.05, drop=1e-5:
    2k->rel 3.6e-1, 5k->2.4e-1, 10k->1.2e-1, 50k->9.9e-2, unbounded(259k)->1.1e-4.
    Hard-capping bounds RSS (169->391MB) but costs accuracy STEEPLY here: the XY
    ladder operator genuinely spreads to ~250k terms.

(2) drop_tol sweep @ dt=0.05, unbounded: 3e-3->rel 6.8e-2 (2.6k terms, 282MB, 1.0s);
    1e-3->2.4e-2 (21k, 318MB); 3e-4->2.9e-3 (190k, 389MB); 3e-5->3.2e-4 (256k, 392MB).
    drop_tol keeps the magnitude-selected set and is the efficient knob;
    max_basis is a safety bound, not the primary control.

(3) expm dt sweep @ drop=1e-4 (KEY): dt=0.025->rel 3.0e-3/27.6s, 0.05->2.2e-3/15.4s,
    **0.1->7.4e-4/8.2s**, 0.2->2.3e-2/4.1s, 0.4->4.4e-1/0.9s. Optimal dt ~0.1 —
    MORE accurate AND ~2x faster than dt=0.05, because expm has no Trotter
    splitting error and fewer steps => fewer truncation events => less
    accumulated truncation loss. Earlier benchmarks under-tuned expm at dt=0.05.

(4) Trotter dt sweep @ mac=1e-5: 0.05->rel 9e-4/16.1s, 0.1->1.9e-3/8.2s,
    0.2->7.4e-3/3.9s; RSS flat ~173MB.

## Tuned conclusion

At N=10 with BOTH methods dt-tuned, wall is ~tied (~8s at rel~1e-3, dt=0.1);
expm uses ~2.3x Trotter RSS. Memory (not speed) is expm's cost, and it is only
reducible by accepting higher rel error (coarser drop / larger dt) — there is no
free code-level fix. At higher rel error (~1e-2) both are far cheaper and expm's
small-basis regime (drop=3e-3: 2.6k terms, 1s) is very competitive on wall.

## Tuned N=12 — CORRECTS the earlier "expm wins 0.6x" claim

The [[../2026-07-01-expm-pc-step]] scan reported expm 0.6x Trotter wall at N=12,
rel=1e-2. That was an ARTIFACT: it interpolated against a noisy Trotter point
(mac=3e-5 -> 491s while the TIGHTER mac=1e-5 -> 302s; non-monotonic = system-load
outlier). Re-measured cleanly, dt-tuned:

  expm    dt0.1 drop1e-3  rel 3.8e-2   24.8s   827MB   (41k terms)
  expm    dt0.1 drop3e-4  rel 2.1e-2  189.6s  1823MB  (518k)
  trotter dt0.1 mac1e-4   rel 1.6e-2   40.3s   263MB  (781k)
  trotter dt0.2 mac1e-4   rel 1.4e-2   34.2s   365MB  (1.32M)

At matched rel~1.5e-2, Trotter is ~34-40s / 263-365MB; expm (interp drop~4e-4)
is ~120s / ~1.4GB. => at N=12 Trotter wins BOTH wall (~3x) and memory (~5x).

Note expm's 827MB for only 41k final terms (drop=1e-3): the RSS is dominated by
the predictor-corrector LEAKAGE TRANSIENT (basis is enriched 2x before pruning,
and the action cache is built on the enriched basis), not the final basis. This
transient is inherent to the PC method and is the true memory driver — ~20x more
RSS/term than Trotter at N=12.

## Honest overall verdict (post-optimization, properly tuned)

The ~10.7x per-step win (pc-step experiment) moved expm from ~6-8x SLOWER than
Trotter to roughly COMPETITIVE, but it does NOT clearly beat Trotter:
- N=8:  Trotter wins (expm 3.3x wall, 1.4x RSS).
- N=10: wall ~tied, expm ~2.3x RSS.
- N=12: Trotter wins (~3x wall, ~5x RSS).
expm's structural weakness is memory (the PC leakage transient), growing with N.
Its genuine advantages: no Trotter splitting error (systematic, dt-controlled;
prefers dt~0.1) and a compact final basis at low accuracy. Measurement variance
on single runs at N>=12 is large — trust min-of-N microbenchmarks, not one-shot
walls.

## Wide (dt x drop_tol x max_basis) scan -> optimal expm params (N=10, min of 2)

Reran on a QUIET machine (the first pass was load-inflated). rel_err is
deterministic (bit-identical across both runs); only walls changed — these are
the clean walls.

dt x drop_tol grid (max_basis=inf) reveals a dt-ERROR FLOOR (the PC does only
2 leakage hops/step, so too-large dt can't spread the operator far enough):
- dt=0.05: as accurate as 0.1 but ~1.8x slower (more steps -> more truncation
  events). drop=3e-5 -> rel 3.2e-4 @ 10.4s.
- dt=0.10: OPTIMAL. best accuracy floor AND fastest to reach it:
  drop 3e-3->rel 3.3e-2/0.4s, 1e-3->1.25e-2/2.5s, 3e-4->4.0e-3/4.3s,
  1e-4->7.4e-4/4.7s, 3e-5->1.0e-3/6.9s.
- dt=0.15 -> accuracy floors at ~1.9e-2 (finer drop no longer helps);
  dt=0.20 -> ~2.2e-2; dt>=0.25 -> ~0.15. Above dt~0.1 you cannot reach small
  rel_err at ANY drop_tol. Only exceed dt~0.1 if low accuracy suffices (then
  dt=0.2-0.3 gives rel~3e-2 in ~0.3s).

Optimal (dt,drop) per accuracy target, max_basis=inf, MIN WALL (clean):
  rel~3e-2:  dt=0.2  drop=3e-3  -> 0.3s  (~8k terms)   [very cheap regime]
  rel~1e-2:  dt=0.1  drop=3e-4  -> 4.3s  (206k, rel 4.0e-3)
  rel~3e-3:  dt=0.1  drop=1e-4  -> 4.7s  (244k, rel 7.4e-4)
  rel<=1e-3: dt=0.1  drop=1e-4  -> 4.7s  (rel 7.4e-4)

max_basis interaction (base dt=0.1): capping BELOW the drop-tol-natural basis
size trades accuracy for RAM/wall. At the rel~1e-2 target, max_basis=200k
(vs inf) still meets it at 4.3s/310MB vs 5.3s/387MB — a mild win by trimming the
unneeded tail. For high accuracy the operator genuinely fills ~250k terms, so
max_basis must be >= that. => drop_tol is the primary knob; max_basis is a
memory safety cap that gives marginal savings at moderate accuracy.

RECOMMENDED expm defaults (this ladder/observable): drop_tol chosen for the
accuracy target, and dt chosen so the O(dt^3) per-step floor is below it:
- rel ~1e-2 : dt=0.1,   drop=1e-3   (~2.5s)
- rel ~1e-3 : dt=0.1,   drop=1e-4   (~4.7s)
- rel ~1e-4 : dt=0.05,  drop=3e-6   (~11s)
- rel ~1e-5 : dt=0.025, drop=1e-7   (~30s)
max_basis unbounded (or ~1.2x the saturated basis as a memory guard). Note the
optimal dt DECREASES with the accuracy target: the PC per-step error is ~O(dt^3)
(dt=0.1 floors ~1.1e-3, dt=0.05 ~1e-4, dt=0.025 <1e-5), so tighter targets need
smaller dt AND finer drop (the basis saturates ~261k at N=10, so fine drop is
nearly free once saturated — it only reduces per-step pruning loss).

## Tight-accuracy expm-vs-Trotter (N=10, exact-ED ref, clean, min-wall to reach)

Trotter error = O(dt^2) Strang splitting + truncation(min_abs_coeff); expm error
= O(dt^3) PC + truncation(drop_tol). Both basis-saturate ~261k terms, so it is
STEP COUNT that differs.

  target   expm (dt,drop)          wall   RSS     | Trotter (dt,mac)       wall   RSS
  1e-4     dt0.05 drop3e-6  9.9e-5  11.4s  393MB   | dt0.0125 mac1e-7 3.9e-5 37.2s 173MB
  1e-5     dt0.025 drop1e-7 9.7e-6  30.3s  422MB   | mac1e-7 FLOORS 3.9e-5  (>1e-8 needed)

=> At TIGHT accuracy the ranking INVERTS vs the moderate regime:
   - rel 1e-4: expm ~3.3x FASTER (Trotter's O(dt^2) forces dt=0.0125 = 4x more
     steps than expm's dt=0.05); Trotter still ~2.3x less RAM.
   - rel 1e-5: expm reaches it in 30s; Trotter's mac=1e-7 truncation floors at
     ~4e-5, needs mac<=1e-8 (basis -> near-full, cost >> expm).
Combined with the moderate-accuracy result (Trotter competitive/better on wall,
~2x less RAM for rel >= ~1e-3), the honest verdict:
   * moderate accuracy (rel >~ 1e-3): Trotter wins (wall ~tie, memory ~2x).
   * tight accuracy (rel <~ 1e-4):    expm wins WALL (higher-order in dt), and is
     the practical choice for rel ~1e-5+; Trotter keeps the ~2x memory edge.
This is the regime where the matrix-exponential method earns its keep.
