# Orbit-preserving Trotter vs expm comparison (2026-07-03)

Goal: compare Trotter vs expm in the ORBIT-PRESERVING (translation-symmetry /
momentum-k) representation, precision-matched at rel~1e-3, on an example with
~1M final basis. For expm, compare stream vs cache (leakage) with a max_basis cut.
Scan dt and drop_tol for both.

## Stage 0 (DONE): foundation verified vs exact ED (N=8 XY chain, k=1)
Both orbit-preserving methods reproduce the exact momentum-k autocorrelation
C_k(t) = <O_k, e^{tL}O_k>/<O_k,O_k>, O_k = sum_j e^{-2pi i k j/N} Z_j:
- expm orbit-rep (pc_step_orbit_rep, complex coeffs) == exact (to ~1e-4).
- Trotter symmetry-merged (real-pair PauliSum + momentum_merge) == exact
  (real part ~1e-4; tiny O(dt^2) equivariance drift in imag).
KEY API notes learned:
- expm: seed = k-mode Z (basis Z_j, coeffs e^{-i phi j}); project via
  canonicalize_basis_arr_complex(basis,coeffs,group,mom); evolve
  pc_step_orbit_rep(basis,coeffs,dt,max_basis,group,mom,drop,protected).
  TranslationGroup.chain_1d(N), mom=[k]. C_k via overlap with the seed reps.
- Trotter: keep O_k as a REAL PAIR (R=sum cos(phi j)Z_j, Im=sum -sin(phi j)Z_j);
  apply the SAME real gates to both; R.momentum_merge(Im,group,mom) folds to
  orbit-rep. GATE ANGLE = 2*c*dt (Rxx=e^{-i th/2 XX}, H coeff c=1) -- using dt
  (not 2dt) decays ~4x too slow. C_k via PauliSum.overlap:
  <A+iB,C+iD> = (<A,C>+<B,D>) + i(<A,D>-<B,C>), norm = <R,R>+<Im,Im> = N.
Harnesses: orbit_verify.py (expm), trotter_orbit_verify.py (Trotter).

## Remaining stages
1. DONE (858bbef7): orbit-rep streaming via PPVM_EXPM_STREAM (StreamOrbitOp).
2. DONE (2026-07-04): scaled to 1.40M orbit-rep peak at L=7 -- and got an EXACT
   ED reference there instead of a converged one (see below).
3. DONE (2026-07-04): dt x knob scan both methods at L=7; see the at-scale
   section below. NOTE: the L=5 Trotter numbers above used the legacy
   end-of-step truncation scheme; see the scheme-fix note below.

## L=5 (N=10) orbit k-resolved: Trotter vs expm, exact-ED referenced

Ladder L=5, k=2, T=2, momentum-k Z autocorrelation C_k(t). 2nd-order Strang
Trotter (palindrome, matches the verified real-space scheme) vs expm orbit-rep
(pc_step_orbit_rep). drop/mac=1e-4 (fine -> dt-floor-limited); dt scanned.

dt-convergence (drop=1e-4):
  expm:    dt0.1->2.0e-2  dt0.05->1.9e-3  dt0.025->7.8e-4   (~O(dt^3))
  trotter: dt0.1->1.15e-2 dt0.05->2.9e-3  dt0.025->1.26e-3  (O(dt^2), confirms 2nd order)

Matched rel~1e-3 (interp):
  expm:    ~47s, peak_basis 52k,  RSS ~262 MB
  trotter: ~44s, peak_basis 256k, RSS ~270 MB
=> WALL ~TIED (~1.05x, Trotter marginally faster). RSS ~TIED. But expm keeps
   ~5x FEWER terms (52k vs 256k) -- the k-resolved orbit reduction is real and
   large. It does NOT convert to a net wall/RAM win here because expm's per-term
   cost is ~higher (complex CSC action cache + doubly-enriched transient +
   complex coeffs), which cancels the 5x term-count advantage.
Caveat: at L=5 the RSS (~250-270MB) is baseline-dominated (python/numpy ~200MB),
which masks the per-term memory difference; at larger basis the expm-heavier-
per-term effect would be more visible in RSS.

Stream vs cache with a max_basis=15k CUT (from the first pass, dt0.1 drop1e-4):
  cache:  5.6s, 217 MB ; stream: 73.4s (~13x), 198 MB (~9% less), identical rel.
=> streaming is a LARGE-basis RAM tool; at 15k it costs ~13x wall for ~9% RAM
   (the cache is a small fraction of the ~200MB baseline there). Not worth it
   at this size.

## At scale: L=7 (N=14) ladder, k=3, T=2, EXACT-ED referenced (2026-07-04)

Design: instead of a self-converged reference at L=8+ (probes showed drop<=3e-4
grows unboundedly there: 904k terms at t=0.7, still x1.8/step), we chose L=7 --
big enough to be truncation-limited everywhere probed (orbit ceiling 4^14/7 =
38M, bases reached 5.5M), small enough that the EXACT reference is computable:
dim 2^14 dense ED. exact_ref_L7.py computes C_k(t) via the |U|^2 bilinear form
(same formula as the verified L=5 driver), with two speedups: H(XY) is REAL
symmetric -> dsyevd, and U-blocks via two real dgemms. 44 min one-time, saved
to exact_ref_L7_k3_T2.npz. Validated at L=5 against the known REF (5e-7).

### TROTTER SCHEME FIX (important)
orbit_bench.py ran the WHOLE Strang palindrome with truncate=False and one
truncate() per step ("legacy"). The VERIFIED real-space baseline
(demo/trotter_ladder.py) truncates per bond (rxx untruncated, each ryy
truncates). At N=10 the difference is mild (L=5, dt.05/mac1e-4: legacy rel
2.9e-3/44s, per-bond 3.7e-3/35s, same 256k peak). At N=14 it is fatal for the
legacy scheme: mac=3e-3 dt=0.1 ran >2.5h at 10.7GB RSS (killed; intra-palindrome
transient explodes over 84 untruncated gates) where per-bond takes 1.0s. All
L=7 Trotter rows below use per-bond ("trotter-pb", the honest baseline).
=> The L=5 "wall ~tied" conclusion above is an artifact of the weak legacy
   baseline. The at-scale picture replaces it.

### Scan (harness orbit_bench_L7.py; uncontended M-machine, 10 cores, 32GB)
method       dt    knob   rel_err  wall_s  RSS_mb    peak
expm-cache  0.1   3e-3   1.87e-2     79     965     23.5k
expm-cache  0.1   1e-3   2.21e-2    709    5334      212k
expm-cache  0.05  3e-3   7.41e-3    107     625     12.3k
expm-cache  0.05  1e-3   9.17e-3    961    3510      116k
expm-cache  0.05  3e-4   2.08e-3   5397    9240     1.40M   <- best rel achieved
expm-cache  0.025 1e-3   6.51e-3    951    1802     48.5k
trotter-pb  0.1   3e-3   3.23e-2      1     164     18.9k
trotter-pb  0.1   1e-3   1.45e-2      9     296      160k
trotter-pb  0.05  1e-3   3.01e-2      4     183     58.6k
trotter-pb  0.05  3e-4   4.13e-3    103     722      637k
trotter-pb  0.05  1e-4   3.31e-3   1412    4155     5.46M   <- trotter frontier
trotter-pb  0.025 3e-4   5.79e-3     58     308      208k
trotter-pb  0.025 1e-4   5.61e-3    761    1369     1.82M

### Findings
1. BOTH methods have a dt<->knob coupling: at fixed knob, SMALLER dt can be
   WORSE (trotter mac=3e-4: dt.05 4.1e-3 -> dt.025 5.8e-3; expm drop=1e-3:
   dt.05 9.2e-3 -> dt.025 6.5e-3 barely helps). Trotter: smaller rotations
   put more new terms under mac, error compounds PER GATE (84 gates/step).
   expm: tau_add = K*drop/dt stiffens the end-filter as dt shrinks. Tuning
   must move dt and knob together.
2. PRECISION FRONTIER: expm reaches rel 2.08e-3 (1.40M terms, 90 min, 9.2GB).
   Trotter saturates at ~3.3e-3: tightening mac 3e-4 -> 1e-4 bought x1.25 in
   rel for x8.6 terms (637k -> 5.46M) and x14 wall -- extrapolating to 2e-3
   needs >>20M terms, RAM-infeasible here. Below rel ~3e-3, expm is the only
   method that runs at all on this machine.
3. CROSSOVER at rel ~3e-3: looser than that, trotter-pb is far faster on wall
   (its per-gate PauliSum engine is much cheaper per term than expm's complex
   CSC/Krylov step: at rel~4e-3 trotter costs ~100s vs expm ~2-3000s interp).
   Tighter than that, trotter's cost curve goes vertical (finding 2) while
   expm still converges (dt^3 floor + per-STEP error accumulation).
4. TERM EFFICIENCY unchanged from L=5: at the frontier expm holds 3.9x fewer
   terms (1.40M vs 5.46M) at BETTER rel -- but expm's per-term RAM is ~5x
   (9.2GB vs 4.2GB at those sizes: complex coeffs + doubly-enriched transient
   + CSC cache), so the term advantage does not convert to an RSS win.
5. STREAM vs CACHE at mb=100k cut (dt.1/3e-4): identical rel 6.13e-2 (good),
   stream 649s vs cache 27s (24x) for 273 vs 288MB (-5%). Streaming's RAM
   savings only matter at multi-M bases where its wall (24x on 90min) is
   impractical -> cache is the right default on the orbit path; max_basis is
   the RAM lever but the hard cap costs a lot of accuracy (6.1e-2 capped at
   100k vs 2.2e-2 uncapped 212k at same dt) -- prefer picking drop to fit the
   RAM budget over capping (consistent with real-space scan_admission).
6. Growth probes (growth_probe.py): L=7 drop 1e-3 peaks ~212k at t=1.0 then
   declines; L=7 drop 3e-4 crosses 1M around t~1.0 (peak 1.40M full run);
   L=8 drop 1e-4 was 904k at t=0.7 still x1.8/step (unbounded within reach).

### Paper implication (momentum section)
Punchline shifts from "wall ~tied, 5x fewer terms" (L=5, weak baseline) to:
k-resolved CTPP owns the high-precision regime -- exact-in-dt (O(dt^3) two-hop
floor) + per-step truncation beats 2nd-order Strang whose per-gate truncation
error accumulation makes rel <~ 3e-3 unreachable at N=14 regardless of budget;
at loose precision gate-based propagation is the faster tool. Both statements
are exact-ED-referenced at N=14 with bases to 5.5M terms.

## k=1 (hydrodynamic mode) -- THE CANONICAL ROW (2026-07-04, user: "k=1 is enough")

Same L=7/T=2 setup, exact ED ref exact_ref_L7_k1_T2.npz (17 min uncontended;
C_k1 decays slowly through zero at t~1: 1, .886, .611, .307, .077, -.049, ...).
Driver takes PPVM_BENCH_K to select the reference.

method       dt    knob   rel_err  wall_s  RSS_mb    peak
expm-cache  0.1   3e-3   4.63e-3     30     681     19.0k
expm-cache  0.05  3e-3   1.62e-2     32     418      9.7k  <- coupling: worse than dt=.1
expm-cache  0.05  1e-3   1.81e-3    432    2828     84.8k  <- beats trotter's ceiling, 20x fewer terms
expm-cache  0.05  3e-4   5.46e-4   2752    9274     1.01M  <- SUB-1e-3, frontier
expm-cache  0.025 1e-3   6.02e-3    357    1224     35.9k
trotter-pb  0.1   3e-3   1.75e-2      1     162     17.5k
trotter-pb  0.1   1e-3   1.29e-2      5     261      142k
trotter-pb  0.05  1e-3   1.89e-2      3     179     55.5k
trotter-pb  0.05  3e-4   5.82e-3     44     706      567k
trotter-pb  0.025 3e-4   2.80e-3     20     274      190k
trotter-pb  0.025 1e-4   2.73e-3    280    1310     1.65M  <- ceiling: x8.7 terms for -2.5%

k=1 findings (sharper than k=3, same structure):
1. Trotter ceiling ~2.7e-3: mac 3e-4 -> 1e-4 bought 2.80 -> 2.73e-3 (x8.7 terms,
   x14 wall). Same saturation as k=3, now at a lower level (slow mode = smaller
   effective Trotter error), still a dead end.
2. expm goes SUB-1e-3: 5.46e-4 at 1.01M terms / 46 min / 9.3GB -- 5x more
   accurate than trotter's ceiling at ~0.6x the terms.
3. At trotter's ceiling precision (~2.7e-3): expm dt.05/drop1e-3 gives 1.81e-3
   at 432s vs trotter 280s -- wall comparable, terms 85k vs 1.65M (20x fewer),
   so at matched precision the orbit advantage DOES convert at k=1.
4. dt<->knob coupling confirmed on the expm side too: drop 3e-3, dt .1 -> .05
   WORSENS rel 4.6e-3 -> 1.6e-2 (tau_add = K*drop/dt stiffens the end-filter).
5. Paper row: k=1, exact-ED referenced, N=14 -- "CTPP reaches 5e-4 where
   2nd-order Strang saturates at 3e-3; at Strang's own best precision CTPP
   holds 20x fewer terms at equal wall."

## RAM comparison, best-of-each-method at matched precision (2026-07-04)

(All RSS from ru_maxrss in the worker; ~150-200 MB of it is python/numpy
baseline, which matters only for the small-basis cells.)

k=1, picking each method's CHEAPEST cell at a given rel:
  rel ~1.3-1.9e-2 (loose):   trotter .1/1e-3    13 e-3   4.5s   261 MB
                             expm    .05/3e-3   16 e-3    32s   418 MB
                             -> trotter ~7x wall, ~1.6x RAM cheaper
  rel ~5e-3 (mid):           trotter .05/3e-4  5.8e-3     44s   706 MB
                             expm    .1/3e-3   4.6e-3     30s   681 MB
                             -> genuinely TIED on both wall and RAM
                                (but dominated: see next row)
  rel ~2.8e-3 (trotter's     trotter .025/3e-4 2.8e-3     20s   274 MB
  ceiling):                  expm    .05/1e-3  1.8e-3    432s  2828 MB
                             -> trotter ~22x wall, ~10x RAM cheaper
  rel <= 1e-3:               expm    .05/3e-4  5.5e-4   46min  9274 MB
                             trotter: unreachable at any budget tried
                                (1e-4: 2.73e-3 at 280s/1310 MB, saturated)

k=3, same exercise: trotter's cheapest at its ~4e-3 class is .05/3e-4
(4.1e-3, 103s, 722 MB); expm needs .05/3e-4 (2.1e-3, 90min, 9240 MB) to beat
it on rel -- again trotter is cheaper in wall AND RAM everywhere it can reach.

CORRECTION to k=1 finding 3 above: "20x fewer terms at comparable wall" compared
expm against trotter's SATURATED 1e-4 cell (280s/1310 MB). Trotter's efficient
cell at the same rel (.025/3e-4: 20s/274 MB) is ~22x faster and ~10x lighter
than expm's 1.81e-3 cell. Fewer terms (85k vs 190k) does NOT convert to less
RAM: expm's per-term footprint (complex coeffs + CSC cache + enriched
transient) eats the entire term advantage.

HONEST HEADLINE (both k): within the precision range Trotter can reach at all,
gate-based propagation is cheaper in BOTH wall and RAM (up to ~20x/~10x at its
ceiling). CTPP's momentum-space value is exclusively the extended precision
range (5x lower error at k=1: 5.5e-4 vs 2.7e-3 ceiling, at 46min/9.3GB) plus
what gate methods cannot do at all (Lindbladians, exact-in-dt). The paper
must state the trade this way and not claim a term-count win as a memory win.

## Harness migration to xy-experiments + dt x drop x K mid-precision scan (2026-07-04)

User preference: use the clean xy-experiments harnesses (k_pec_run.py /
k_trotter_run.py + main_k_{pec,xy}_ladder.py) instead of orbit_bench_L7.py.
Differences that matter:
- TROTTER: k_trotter_run evolves the momentum-MERGED pair (momentum_merge
  every step). orbit_bench evolved the full real-space pair (merge only at
  readout). Merging redistributes weight onto reps (coefficient = sum over L
  translates), so the SAME min_abs_coeff is effectively ~L x looser: better
  rel per knob, but far bigger live basis. Neither variant dominates -- see
  below. Both validated vs exact ED at L=5 (pec 1.86e-3 == ledger; trot
  3.20e-3 vs 2.9e-3 unmerged).
- PEC: identical pc_step_orbit_rep kernel (k_pec_run passes no protected_arr;
  negligible). Records C_k every step; supports dephasing jumps.
- Fixed stale `from ppvm_python_native import ...` (module no longer
  installed; symbols live in ppvm._core) -- committed in xy-experiments
  (a044a27). Peak RSS measured externally via /usr/bin/time -l (the h5
  peak_rss_mb attr is psutil END-rss, not peak).
Driver: scan_xy_mid.py (this folder). L=7, k=1, T=2, exact-ED referenced.

Scan (K = PPVM_K_LEAKAGE, tau_add = K*drop/dt; K=0 no end-filter):
method     dt     drop   K   rel_err  wall_s  peakRSS  peak_reps
trot(m)   0.1    1e-3    -   8.12e-3    70      910MB    906k
trot(m)   0.1    3e-4    -   7.42e-3   718     7.0GB    6.90M   <- dt-floor, wasted
trot(m)   0.05   3e-3    -   1.49e-2     7      171MB     61k
trot(m)   0.05   1e-3    -   7.60e-3    65      586MB    522k
trot(m)   0.05   3e-4    -   2.55e-3   728     3.6GB    4.75M   <- below unmerged ceiling!
trot(m)   0.025  3e-3    -   3.09e-2     4      133MB     22k   <- coupling
trot(m)   0.025  1e-3    -   7.51e-3    38      247MB    189k
trot(m)   0.0125 1e-3    -   1.86e-2    21      168MB     62k   <- coupling
pec       0.1    3e-3    0   5.10e-3    32      647MB     19k
pec       0.1    3e-3    1   7.26e-3     8      316MB     40k
pec       0.1    3e-3    5   2.36e-2     1      167MB      4.5k <- over-filtered
pec       0.05   3e-3    1   9.45e-3     7      249MB     21k
pec       0.05   3e-3    5   7.88e-3     1      136MB      1.3k <- absurdly cheap ~8e-3
pec       0.1    1e-3    0   3.40e-3   341     3.8GB     159k
pec       0.1    1e-3    1   3.17e-3    93     1.5GB     354k
pec       0.1    1e-3    5   1.01e-2     8      273MB     43k
pec       0.05   1e-3    1   1.30e-3    82      808MB    211k   <- STAR CELL
pec       0.05   1e-3    5   7.26e-3     4      219MB     11k
pec       0.025  1e-3    1   5.59e-3    45      329MB     70k
(NOTE trot(m) 0.1/1e-4 was killed pre-emptively: extrapolated ~30M reps /
 ~30GB. merged 0.025/3e-4 not run -- gap, likely dominated anyway.)

FINDINGS:
1. K=1 IS A STRICT PARETO WIN for pec at drop=1e-3: vs K=0 at dt=0.05
   (old-harness K0: 1.81e-3/432s/2.8GB) K=1 gives 1.30e-3 / 82s / 808MB --
   better rel, 5.3x wall, 3.5x RAM. Same at dt=0.1 (3.17e-3/93s/1.5GB vs
   3.40e-3/341s/3.8GB). K=5 over-filters at drop<=1e-3 but is a bargain
   dial at loose drop (7.9e-3 at 1s/136MB!). K=1 default recommended for
   momentum runs; the old K-leakage guidance ("K~1 helps 1.1-2x") strongly
   UNDERSTATED the k=1 benefit.
2. Merged vs unmerged Trotter: NEITHER dominates. Merged reaches 2.55e-3
   (unmerged saturated at 2.7e-3) but at 728s/3.6GB/4.75M reps; unmerged
   0.025/3e-4 hit 2.80e-3 at 19.5s/274MB. Effective-knob loosening from the
   merge is the cause. "Trotter" below = best variant per point.
3. dt<->knob coupling is now confirmed in BOTH methods and BOTH trotter
   variants (trot(m) 0.025->0.0125 at 1e-3: 7.5e-3 -> 1.9e-2).

MID-ACCURACY (~5e-3) OPTIMAL POINTS -- the head-to-head requested:
  trotter (unmerged, dt=.025, mac=3e-4): rel 2.80e-3, 19.5s, 274 MB
  pec (dt=.1, drop=3e-3, K=0):           rel 5.10e-3, 32 s,  647 MB
  pec (dt=.025, drop=1e-3, K=1):         rel 5.59e-3, 45 s,  329 MB
  -> AT the mid class Trotter's best point is still the cheapest (1.6-2.3x
     wall, 1.2-2.4x RAM) and overshoots accuracy.
  BUT one small step up in budget flips it: pec .05/1e-3/K1 = 1.30e-3 at
  82s/808MB -- 2.1x better rel than Trotter's best-ever (2.55e-3 at
  728s/3.6GB), at 8.9x less wall and 4.5x less RAM than that cell.
REVISED CLASS BOUNDARY (k=1): Trotter owns rel >~ 2.7e-3 (cheap, saturating);
pec+K1 owns rel <~ 2.5e-3 outright (wall AND RAM). The pec frontier costs
from the earlier k=1 section (K=0) are superseded by ~5x cheaper K=1 cells.

## Merged vs unmerged Trotter: knob-calibrated dissection (2026-07-04, user challenge)

User objection: merging is just a re-representation, the Lx storage win should
be unconditional. Test: merged cells at tau' = L*tau (folded coefficients are
~Lx larger, so same tau ==> ~Lx weaker truncation). dt=0.025, k=1:

  merged  7e-3 (=unmerged 1e-3):  2.23e-2   1.4s  113MB    4.0k reps
  unmerged 1e-3:                  3.27e-2   1.8s  161MB   18.6k terms
  -> LOOSE REGIME: user is right. Merged wins on rel AND storage (~L x
     compression) AND RAM at calibrated knobs.

  merged  2.1e-3 (=unmerged 3e-4): 2.75e-2   8.7s  145MB   43.8k reps
  merged  1e-3:                    7.51e-3    38s  247MB    189k reps
  merged  5e-4:                    3.49e-3   166s  620MB    717k reps
  unmerged 3e-4:                   2.80e-3  19.5s  274MB    190k terms
  -> TIGHT REGIME: merged converges but needs ~3.8x more STORED reps (each
     rep = 7 translates, so ~26x more represented weight) for slightly worse
     rel; unmerged wins ~8x wall / ~2.3x RAM at the ~3e-3 class. Merged rel
     is NON-MONOTONIC in drop (2.2e-2 -> 2.75e-2 -> 7.5e-3 -> 3.5e-3), a
     coherent-error signature.

Interpretation (mechanism hypothesis, not fully dissected): in merged form
every per-gate truncation decision on the mid-layer transient drops ALL L
translates coherently (error amplitude ~L per decision); in unmerged form
drops are translate-level with quasi-random phases (partial cancellation),
and the readout sector-projection filters out-of-sector truncation noise for
free. Loose truncation: few decisions -> storage win dominates. Tight: the
coherent amplification dominates.

NOTE the pec orbit-rep path is the *proper* re-representation the objection
assumes: it never unfolds (whole evolution in canonical rep space, truncation
once per step on the invariant |c_orbit|), so it gets the Lx compression
without the per-gate coherent-drop pathology. Paper angle: symmetry-compressed
truncation is natural in the CTPP orbit basis; grafting it onto gate-based
propagation via per-step fold/unfold injects coherent truncation error.
