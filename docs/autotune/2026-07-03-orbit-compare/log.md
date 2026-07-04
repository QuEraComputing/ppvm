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

## Extended K=5 grid (2026-07-04, user: "K=5 got very close to Trotter")

pec, k=1, all dt x drop with K=5 (tau_add = 5*drop/dt):
  dt     drop   rel_err  wall_s  peakRSS  peak_reps
  0.025  3e-3   1.04e-1    0.8    121MB      384    <- over-filtered (drop/dt big)
  0.025  1e-3   5.92e-2    1.9    156MB     3.0k    <- over-filtered
  0.0125 1e-3   3.82e-2    1.4    131MB     0.8k    <- over-filtered
  0.1    3e-4   2.86e-3   69.0    1.4GB     537k
  0.05   3e-4   2.78e-3   35.2    524MB     134k   <- AT Trotter's best rel
  0.025  3e-4   4.82e-3   17.6    277MB      34k   <- mid class, Trotter-like cost
  0.0125 3e-4   6.18e-3    9.5    196MB     8.4k
  0.1    1e-4   3.89e-3    869    5.9GB    4.75M   <- dominated (dt floor + huge basis)
  0.05   1e-4   1.04e-3    500    2.1GB    1.30M
  0.025  1e-4   1.22e-3    225    923MB     325k

FINDINGS:
1. USER HUNCH CONFIRMED: pec 0.05/3e-4/K5 = 2.78e-3 at 35s/524MB sits at
   Trotter's best-ever rel (2.80e-3 at 19.5s/274MB) within 1.8x wall / 1.9x
   RAM. With K=0 the same rel class cost 432s/2.8GB -- the end-filter buys
   ~12x wall / ~5x RAM at Trotter-ceiling precision.
2. K=5 over-filters iff drop/dt is large (tau_add >~ 0.1); rule of thumb:
   pick K*drop/dt ~ 0.01-0.03.
3. Combined pec Pareto frontier (k=1) now TRACKS Trotter down to its ceiling
   and extends beyond:
     7.9e-3: 1s/136MB (3e-3,K5)   | 4.8e-3: 18s/277MB (0.025/3e-4,K5)
     2.8e-3: 35s/524MB (0.05/3e-4,K5) | 1.3e-3: 82s/808MB (0.05/1e-3,K1)
     1.0e-3: 500s/2.1GB (0.05/1e-4,K5) | 5.5e-4: 46min/9.3GB (0.05/3e-4,K0)
   Trotter's remaining edge at its own best point is <2x wall/RAM; below
   ~2.5e-3 pec owns everything.
REVISED HEADLINE: with the K end-filter tuned (K*drop/dt ~ 0.01-0.03), CTPP
matches gate-based propagation within <2x across Trotter's entire reachable
range and is the only method below ~2.5e-3.

## CORRECTION (2026-07-04, user): "within 2x" is NOT competitive

The bar for "competitive" at Trotter's frontier point (2.80e-3 at 19.5s/274MB)
is ~20s wall. pec's best there (0.05/3e-4/K5: 2.78e-3, 35s, 524MB) MISSES it
(1.8x wall, 1.9x RAM). Retract the "tracks within <2x = competitive" framing.
Standing summary: at rel >~ 2.5e-3 Trotter remains the cheaper method; pec is
the only method below ~2.5e-3, and K-tuning cut its cost there by ~5-12x vs
K=0. Untested knobs if the ~20s target is pursued later: interior K (3-7) at
0.05/3e-4, drop between 3e-4 and 1e-3 with K rescaled to keep K*drop/dt ~
0.01-0.03, dt ~ 0.03-0.04.

## METRIC CHANGE: median pointwise relative error (2026-07-04, user directive)

New standard metric: median over the 11 sample times of |C(t)-REF(t)|/|REF(t)|
(pointwise rel first, then aggregate; NOT norm-ratio, NOT avg-abs/avg-ref).
Rationale: L2 norm-ratio is dominated by early times where |REF| is O(1);
the mean of pointwise rel is dominated by the zero crossing at t=1.8
(REF=-0.0026); the median measures typical relative accuracy across the
decay. Adopted after user review.

All cells with stored curves, recomputed (k=1; wall/RSS from original runs;
old-harness cells re-run ONLY if wall < 30s per user rule):

  method          K   drop    dt    median_rel   wall    RSS
  pec             1   1e-3   0.05    1.02e-3      81s   808MB  <- best overall
  pec             5   1e-4   0.05    1.36e-3     499s   2.1GB
  pec             5   1e-4   0.025   2.29e-3     224s   923MB
  pec             1   3e-3   0.05    5.39e-3     6.1s   249MB  <- mid winner
  pec             5   1e-3   0.05    7.21e-3     3.5s   219MB  <- loose winner
  pec             5   3e-4   0.0125  6.94e-3     8.8s   196MB
  pec             5   3e-4   0.05    2.23e-2      34s   524MB  <- ex-"star cell"!
  pec             0   3e-3   0.1     1.65e-2      32s   647MB
  trot(merged)    -   5e-4   0.025   5.42e-3     165s   620MB  <- trot best
  trot(merged)    -   1e-3   0.05    5.68e-3      64s   586MB
  trot(merged)    -   3e-4   0.05    5.86e-3     728s   3.6GB
  trot(unmerged)  -   3e-4   0.025   1.16e-2      21s   274MB
  trot(unmerged)  -   1e-3   0.1     5.37e-2      7s    261MB
  (trot(unmerged) 0.025/1e-4 and pec K0 0.05/3e-4 NOT re-measured: >30s rule.)

FINDINGS UNDER THE MEDIAN METRIC:
1. THE STORY FLIPS: pec wins EVERY precision class, by large margins.
   ~7e-3: pec 3.5s/219MB vs trot 21s/274MB (unmerged, 1.16e-2).
   ~5.5e-3: pec 6.1s/249MB vs trot(merged) 64s/586MB -- 10x wall, 2.4x RAM.
   ~1e-3: pec 81s/808MB; NO measured Trotter cell below 5.4e-3.
   Trotter saturates at median ~5.4e-3 among measured cells: its tight-knob
   cells improve early-time (L2) accuracy but not late-time relative accuracy.
2. METRIC SENSITIVITY IS LARGE. The L2 "star cell" (pec 3e-4/0.05/K5:
   L2 2.78e-3) has median 2.23e-2 -- excellent early, mediocre late.
   Conversely pec K1 1e-3/0.05 is uniform (L2 1.30e-3, median 1.02e-3).
   Same for trot(unmerged) 0.025/3e-4: L2 2.80e-3 but median 1.16e-2.
   The L2-based conclusion "Trotter cheaper at rel >~2.5e-3" was an artifact
   of early-time weighting; for k-resolved DECAY-RATE physics (late-time
   relative accuracy), the median verdict is the physically relevant one.
3. Best-known settings under median: pec dt=0.05 with (drop=3e-3,K=1) for
   fast/mid, (drop=1e-3,K=1) for precision. dt=0.05 is consistently the pec
   sweet spot at T=2.

## Merged vs unmerged Trotter, revisited under the MEDIAN metric (2026-07-04)

Calibrated pairs (tau' = L*tau = 7*tau), dt=0.025, median pointwise rel:
  merged 7e-3:  1.35e-1  0.7s   4.0k reps | unmerged 1e-3: 5.94e-2  1.8s  18.6k
  merged 2e-3:  3.33e-2  8.0s  43.8k reps | unmerged 3e-4: 1.16e-2 21.2s   190k
-> UNMERGED is ~2.3-2.9x BETTER rel at BOTH calibrated pairs. The earlier L2
   result "merged wins the loose pair" was an early-time-weighting artifact;
   under the median (late-time-sensitive) metric the coherent-drop penalty of
   merged evolution shows at every knob setting. Unmerged costs ~2.6x wall and
   ~4.3x stored terms at the tight pair (its per-translate storage), roughly
   ~1.5x RSS.

Frontier (best measured, any knob):
  merged:   5.42e-3 @ 165s/620MB (0.025/5e-4);  5.68e-3 @ 64s/586MB (0.05/1e-3)
  unmerged: 1.16e-2 @ 21s/274MB (0.025/3e-4)   [tighter cells >30s: unmeasured]
-> merged reaches ~2x lower median than unmerged's best MEASURED point, but
   only because it was run at tighter effective knobs (tau'=5e-4 = unmerged
   mac 7e-5, whose unmerged run would be a multi-M-term/,>400s cell). At
   matched effective truncation unmerged wins rel; at matched rel ~1.2-1.5e-2
   they cost about the same (unmerged 1.16e-2@21s/274MB vs merged
   1.46e-2@38s/247MB).

Merged saturation: mac 1e-3 -> 3e-4 at dt=0.05 gives 5.68e-3 -> 5.86e-3 (11x
wall for nothing), and dt=0.025/5e-4 gives 5.42e-3 -- a dt-INDEPENDENT floor
~5.4e-3, consistent with coherent whole-orbit truncation noise rather than
Trotter dt error.

Both variants remain dominated by pec under the median metric: pec K1
3e-3/0.05 = 5.39e-3 @ 6.1s/249MB matches merged Trotter's best at ~27x less
wall and ~2.5x less RAM.

## CORRECTION: merged vs unmerged at FIXED PRECISION, 1% median target (2026-07-04)

User (correctly) rejected the calibrated-knob verdict: compare at fixed
precision, then wall+RAM. Target: median rel ~1e-2 ("good to go"). Filled the
two missing unmerged cells (curves re-measured, canonical uncontended
wall/RSS from the original runs):
  unmerged 0.05/3e-4:  median 1.75e-2   44s   706MB
  unmerged 0.025/1e-4: median 8.06e-3  280s  1310MB

BEST-OF-VARIANT AT median <= 1e-2 (k=1):
  pec (1e-3/0.05/K5):        7.21e-3     3.5s   219MB   <- overall winner
  trot MERGED (1e-3/0.05):   5.68e-3    63.9s   586MB
  trot UNMERGED (1e-4/0.025):8.06e-3   279.7s  1310MB
-> AT THE 1% TARGET MERGED BEATS UNMERGED: 4.4x wall, 2.2x RAM. The
   calibrated-knob analysis (unmerged 2-3x better rel at same effective tau)
   measured truncation efficiency per threshold -- NOT the practical
   criterion. On the Pareto frontier the picture is:
     >= 1.5e-2: unmerged cheapest (0.05/1e-3: 1.44e-2 at 3.0s/179MB)
     ~1e-2 and below: merged wins among Trotter variants (crossover ~1.2e-2)
   Both variants saturate, but merged saturates LOWER (~5.4e-3) and reaches
   its level far cheaper than unmerged approaches its (~8e-3: x13 wall for
   -30% from 1.16e-2 -> 8.06e-3).
-> So the Lx compression of merged evolution IS the practically better
   scheme at the paper's working precision -- the user's original
   re-representation intuition holds at fixed precision. The coherent-drop
   effect is real (it sets merged's ~5.4e-3 floor and its loose-regime
   disadvantage) but does not decide the 1% comparison.
pec context at 1%: 18x faster and 2.7x lighter than the best Trotter variant.
BENCHMARK TARGET going forward: median pointwise rel ~1e-2 (user directive).

## Per-STEP truncation for merged Trotter: NEGATIVE result (2026-07-04)

Test of "would truncating only once per step (post-merge, canonical reps) fix
merged Trotter's coherent-noise floor?" — k_trotter_run gained a
per_gate_trunc kwarg (main_k_xy_ladder: --gate_trunc 0). Cell: L=7, k=1,
dt=0.05, drop=1e-3, T=2.
RESULT: the run DIED without completing (~1.5h+, peak RSS 16.8 GB, no h5
written) — the unpruned 84-gate intra-step transient explodes, the same
end-of-step-truncation blowup seen unmerged at N=14, only ~L x smaller and
still fatal. Per-step truncation is NOT a viable fix for merged Trotter at
this size; the per-gate/per-bond truncation that causes the coherent-noise
floor is also what keeps the transient bounded. (pec avoids the dilemma
structurally: its generator action never unfolds the orbit basis, so there
is no intra-step transient to prune.)

## REAL-SPACE MSD comparison, L=7, gamma=0, T=2 (2026-07-04)

Observable: MSD(t) from the localized center-rung seed (Z_{j0,0}+Z_{j0,1})/2,
via main_realspace_ladder.py (--mode trotter|adaptive) + msd.py. Metric:
median pointwise rel of MSD(t) over t=0.2..2.0 (t=0 excluded, MSD=0) vs the
sector-reduced exact reference (exact_msd_L7.py: magnetization-block eigh +
|U_m|^2 bilinear -- 18 s total, ~80x faster than the naive full-space |U|^2;
validated at L=5 to machine precision). K = PPVM_K_LEAKAGE reaches pc_step_arr
directly. Driver: scan_realspace_msd.py.

mode      dt     knob    K  median_rel  wall_s  peakRSS  peak_terms
trotter   0.1    1e-3    -   2.23e-2     0.3     112MB     10.1k
trotter   0.05   1e-3    -   3.22e-2     0.2     111MB      3.7k
trotter   0.05   3e-4    -   1.34e-2     2.4     120MB     39.2k
trotter   0.025  3e-4    -   1.34e-2     1.5     113MB     12.9k
trotter   0.05   1e-4    -   5.88e-3    25.8     190MB      337k
trotter   0.025  1e-4    -   7.91e-3    13.7     134MB      111k
trotter   0.025  3e-5    -   2.07e-3   209.4     415MB     1.20M
adaptive  0.1    1e-3    0   5.56e-3     8.9    1135MB     40.2k
adaptive  0.1    1e-3    1   5.19e-3     2.2     324MB     77.8k
adaptive  0.05   1e-3    1   8.87e-3     1.9     314MB     42.0k
adaptive  0.05   1e-3    5   2.70e-2     0.1     154MB      2.8k
adaptive  0.05   3e-3    1   2.89e-2     0.2     171MB      4.7k
adaptive  0.05   3e-3    5   1.72e-1     0.0     116MB       314
adaptive  0.05   3e-4    1   3.97e-3    43.3    1472MB      491k
adaptive  0.05   3e-4    5   8.90e-3     1.6     313MB     34.3k
adaptive  0.025  1e-3    1   1.57e-2     1.9     226MB     16.8k
adaptive  0.1    3e-4    1   2.32e-3    50.5    2977MB      838k
adaptive  0.1    1e-3    5   2.95e-2     0.3     191MB     10.4k

FIXED-PRECISION VERDICT (1% median target):
  pec/adaptive (0.1/1e-3/K1):  5.19e-3   2.2s   324MB
  trotter      (0.025/1e-4):   7.91e-3  13.7s   134MB
  -> pec ~6x FASTER; trotter ~2.4x LIGHTER. Split verdict, unlike momentum
     space (where pec won both axes).
Deeper class (~2e-3):
  pec (0.1/3e-4/K1):    2.32e-3   50.5s  3.0GB
  trotter (0.025/3e-5): 2.07e-3  209.4s  415MB
  -> same split, sharper: pec 4x wall advantage, trotter 7x RAM advantage.
FINDINGS:
1. Real space reverses the RAM story: without the orbit compression, pec's
   per-term footprint (CSC cache + transient) makes it RAM-heavy while
   per-bond Trotter stays lean. Matches the standing handoff note ("expm's
   real-space weakness is memory"). Wall: pec wins every class (4-6x).
2. K=1 again a strict Pareto win over K=0 (0.1/1e-3: 4x wall, 3.5x RAM,
   better rel). K=5 over-filters below drop 3e-4 at these dt.
3. pec's dt optimum here is dt=0.1 (coupling: at fixed drop/K, dt 0.05 and
   0.025 are WORSE - tau_add and per-step pruning tighten with 1/dt).
4. Both methods are dramatically cheaper in real space than momentum space
   at L=7 (seconds, not minutes): the local seed spreads only ~half the ring
   by T=2, whereas the k-mode is delocalized from t=0.

## Real-space MSD: K=5 completion + streaming test (2026-07-04)

adaptive  0.1    3e-4    5   4.58e-3     2.7s    396MB    124k  <- new 1%-class best
adaptive  0.1    1e-4    5   1.87e-3    33.5s   2917MB   1.25M  <- new deep frontier
adaptive  0.05   1e-4    5   6.47e-3    12.2s    839MB    325k
adaptive  0.025  3e-4    5   3.43e-2     0.7s    168MB    9.6k  <- over-filtered
adaptive  0.1    1e-3   K1+stream  5.19e-3   8.5s   344MB  <- identical rel, 3.9x wall, NO RAM gain
adaptive  0.1    3e-4   K1+stream  2.32e-3 245.4s  2574MB  <- identical rel, 4.9x wall, -13% RAM

1. K=5 at dt=0.1 is the real-space sweet spot (tau_add = 5*drop/0.1 lands in
   the ~5e-3..1.5e-2 zone): 0.1/1e-4/K5 dominates the previous deep cell
   (1.87e-3 vs 2.32e-3, 33.5s vs 50.5s, same RAM class).
2. STREAMING DOES NOT RESCUE THE RAM AXIS here: -0..13% RSS for 4-5x wall.
   At these basis shapes RSS is dominated by the leakage transient (basis
   arrays/maps), not the CSC action cache. The earlier "1.8x less RAM"
   characterization applied to larger cache-dominated shapes; for the paper,
   pec's real-space RAM cost is structural at this scale, not a mode choice.
   (All scan cells use the DEFAULT cached-CSC path = the in-basis generator
   built once per step as a sparse matrix and reused across expm matvecs;
   PPVM_EXPM_STREAM=1 recomputes the action per matvec instead.)
Updated 1% verdict: pec 0.1/3e-4/K5 = 4.58e-3 @ 2.7s/396MB vs trotter
0.025/1e-4 = 7.91e-3 @ 13.7s/134MB -> pec ~5x faster, trotter ~3x lighter.
Deep ~2e-3: pec 1.87e-3 @ 33.5s/2.9GB vs trotter 2.07e-3 @ 209s/415MB ->
pec ~6x faster, trotter ~7x lighter. The split verdict stands, sharpened.

## Real-space 1% cell refinement: dt=0.1, drop 4e-4..7e-4 (2026-07-04, user suggestion)

adaptive  0.1  4e-4  K5   5.10e-3    1.4s   324MB   68.7k
adaptive  0.1  5e-4  K5   4.06e-3    0.9s   296MB   43.8k
adaptive  0.1  7e-4  K5   6.00e-3    0.5s   221MB   23.3k   <- RSS-minimal <=7.9e-3
adaptive  0.1  4e-4  K1   3.43e-3   16.1s  1499MB    472k
adaptive  0.1  5e-4  K1   3.20e-3   10.0s   883MB    303k

UPDATED 1% VERDICT (real space): pec 0.1/7e-4/K5 = 6.00e-3 @ 0.5s/221MB vs
trotter 0.025/1e-4 = 7.91e-3 @ 13.7s/134MB -> pec 27x FASTER at 1.65x RSS
(and ~110-150MB of both is python baseline, so basis-attributable RSS is
pec ~80MB vs trotter ~10MB). The interior-drop sweep turned the "split
verdict" into "pec wins wall massively, RAM near-parity" at the 1% target.
The deep-precision split (pec 6x wall / trotter 7x RAM at ~2e-3) stands.
Also: rel is non-monotonic in drop within the 4-6e-3 band (coherent
cancellations) - pick cells by measured rel, not knob interpolation.

## K sweep 1..10 at the real-space sweet spot (dt=0.1, drop=5e-4) (2026-07-04)

K   tau_add  median_rel  wall_s  RSS     peak
1   5e-3     3.20e-3      9.6    881MB   303k
2   1e-2     4.64e-3      5.2    569MB   202k
3   1.5e-2   4.71e-3      2.4    358MB   105k
4   2e-2     5.09e-3      1.4    325MB    64k
5   2.5e-2   4.06e-3      0.9    317MB    44k
6   3e-2     5.00e-3      0.6    268MB    33k
7   3.5e-2   6.07e-3      0.5    237MB    24k
8   4e-2     1.97e-2      0.3    193MB    18k   <- CLIFF
9   4.5e-2   2.03e-2      0.3    208MB    16k
10  5e-2     2.99e-2      0.2    182MB    12k

Shape: a long PLATEAU (rel 3.2-6.1e-3, non-monotonic wiggles) from K=1..7
while wall drops 19x and RSS 3.7x, then a sharp accuracy CLIFF between
tau_add 3.5e-2 and 4e-2 (K=7 -> 8). Practical rule: increase K until just
before the cliff (tau_add <~ 3.5e-2 here); accuracy is insensitive across
the plateau, so K is nearly a free 10-20x wall / 3-4x RAM dial. Refines the
earlier "K*drop/dt ~ 0.01-0.03" guidance: plateau extends to ~0.035.

## 2D scan: K (1..10) x drop_tol at dt=0.1, real-space MSD (2026-07-04)

median_rel (rows=drop, cols=K1..K10; dt=0.1, so tau_add = K*drop*10):
  2e-3: 1.8e-2 2.5e-2 2.3e-2 2.5e-2 3.5e-2 4.2e-2 4.2e-2 7.0e-2 7.2e-2 1.3e-1
  1e-3: 5.2e-3 5.0e-3 5.2e-3 |1.8e-2 2.9e-2 2.9e-2 2.0e-2 2.5e-2 4.1e-2 3.0e-2
  5e-4: 3.2e-3 4.6e-3 4.7e-3 5.1e-3 4.1e-3 5.0e-3 6.1e-3 |2.0e-2 2.0e-2 3.0e-2
  2e-4:   killed 3.0e-3 3.8e-3 4.0e-3 4.7e-3 4.6e-3 5.0e-3 5.8e-3 4.8e-3 5.1e-3
  (| marks the cliff; 2e-4 has no cliff through K=10; 2e-4/K1 aborted -- cell
   >200s, user: abort 1e-4-class cells. 1e-4 row: only K5 = 1.9e-3 measured.)
wall_s ranges: 2e-3 row 0.0-0.5s; 1e-3 row 0.1-2.0s; 5e-4 row 0.2-9.6s;
2e-4 row 1.6-201s (K2). RSS: 120MB (2e-3/K10) to 4.1GB (2e-4/K1-2).

FINDINGS:
1. THE CLIFF IS A PURE tau_add THRESHOLD: 1e-3 cliffs at K=4, 5e-4 at K=8 --
   both tau_add = K*drop/dt = 4e-2; 2e-4 through K=10 only reaches 2e-2 and
   shows no cliff. K and drop enter the end-filter ONLY through tau_add;
   the safe zone is tau_add <~ 3.5e-2 regardless of how it is composed.
2. PLATEAU ACCURACY IS SET BY drop ALONE (the per-step prune): ~5e-3 at
   1e-3, ~3-6e-3 at 5e-4, ~3-5e-3 at 2e-4 -- improving only weakly below
   5e-4 (deepest: 1e-4/K5 = 1.9e-3).
3. RECIPE: pick drop for the target accuracy, then set K ~ 0.03*dt/drop.
   Multiple (drop,K) combos on the same tau_add diagonal give near-identical
   cost -- a degeneracy ridge (e.g. ~5e-3 at ~0.5s/230-270MB via 1e-3/K3,
   7e-4/K5, or 5e-4/K7).

## tau_add scaling study: cliff position vs dt and drop (2026-07-04)

36 cells, real-space MSD, dt in {0.2, 0.05, 0.025} x drop in {1e-3, 5e-4},
fractional K chosen so tau_add = K*drop/dt sweeps {.01,.02,.03,.04,.06,.09}
identically in every row (plus the dt=0.1 rows from the 2D scan).

RESULT: the admission cliff sits at tau_add ~ 0.035-0.04 for EVERY
(dt, drop) combination with dt <= 0.1 -- rows at dt=0.05 and 0.025, both
drops, all depart between tau*=0.03 and 0.06, exactly where dt=0.1 cliffed.
At dt=0.2 the O(dt) plateau (~1.5e-2) masks the cliff (noisy row, bump at
.04 still visible). tau_add is an ABSOLUTE RATE THRESHOLD, independent of
both dt and drop_tol.

=> K is the wrong parameterization: the correct invariant is tau_add itself,
   and K should be DERIVED as K = tau_add * dt / drop_tol. Recommended
   default tau_add ~ 0.02-0.03 (one octave below the cliff). Caveat: the
   cliff value is observable/scale-dependent (momentum k=1 tolerated up to
   ~0.1; here J=1, seed normalized to 1) -- expect tau_add* to scale with
   the Hamiltonian coupling / seed norm, so expose it as a direct parameter
   rather than hardcoding.
Secondary: plateau rel at fixed drop WORSENS as dt shrinks (dt=0.05: 8.8e-3
-> dt=0.025: 1.1e-2 at drop 1e-3; dt=0.1: 5.2e-3) -- the per-step prune
frequency effect; dt~0.1 remains the real-space sweet spot at T=2.

## CORRECTION: "plateau" was wrong; measured sensitivities (2026-07-04)

Local log-log slopes over all measured adaptive pairs (dt<=0.1, real-space
MSD; pairs differ in exactly one knob):
  d(ln rel)/d(ln tau_add), tau_add<0.035:  median +0.24  IQR [+0.05,+0.37]  n=68
  d(ln rel)/d(ln tau_add), tau_add>0.035:  median +0.94  IQR [+0.55,+1.54]  n=152
  d(ln rel)/d(ln drop) at fixed tau_add, below: median +0.03 IQR [-0.01,+0.25] n=32
  d(ln rel)/d(ln drop) at fixed tau_add, above: median -0.02                 n=17
=> 1. There is NO true plateau: rel rises ~ tau_add^0.24 below the cliff,
      steepening to ~ tau_add^1 (with a sharp step) above ~0.04. "Slow rise
      then steepening" is the correct description.
   2. USER OBSERVATION CONFIRMED: in the sampled regime the error depends
      FAR more strongly on tau_add than on drop_tol -- at fixed tau_add,
      drop has near-zero effect (median slope +0.03) except at the very
      lowest tau_add where it re-emerges (~+0.3, e.g. 1e-4/tau=.005:
      1.9e-3 vs 5e-4/tau=.005: 3.2e-3).
   3. Model: rel ~ f(drop) + g(tau_add), f ~ drop^~0.3, g ~ tau_add^~0.24
      below / ^~1 above 0.04; whichever dominates sets the sensitivity. My
      earlier "drop sets the plateau accuracy" conflated the two knobs
      (drop was varied at fixed K, dragging tau_add along).
   4. Practical recipe (revised): lower tau_add and drop TOGETHER so
      f ~ g (balanced); tau_add is the primary accuracy knob in the
      commonly-sampled regime, drop the secondary.

## Addendum: the efficient locus is ~fixed K at fixed dt (2026-07-04, user)

Balance f(drop)~drop^0.3 against g(tau_add)~tau_add^0.24 => the efficient
descent path is tau_add ~ drop^1.25, i.e. K = tau_add*dt/drop ~ drop^0.25*dt:
at FIXED dt this is nearly constant K (data: good cells at K in [1,5] across
drop 1e-3..1e-4 at dt=0.1). So fixed-K tuning accidentally tracks the
efficient path within one dt; its failure modes are (a) dt changes rescale
the filter silently and (b) the hard cliff is absolute in tau_add (~0.04),
so K_cliff = 0.04*dt/drop moves 16x across our grid. The dimensionless knob
that is actually dt-invariant is K' = tau_add/drop (= K/dt); the cliff and
the efficient path are both naturally expressed in (tau_add, drop).

## dt=0.025 verification of the tau_add / K' picture (2026-07-04)

drop   K      tau_add  K'=tau/drop  median_rel  wall_s  RSS
2e-4   0.625  0.005    25           2.27e-3      34.0   881MB
2e-4   1.25   0.01     50           3.34e-3      23.5   882MB
2e-4   2.5    0.02     100          8.84e-3       5.5   327MB
2e-4   7.5    0.06     300          3.13e-2       0.7   177MB  <- cliffed (>0.04) OK
1e-4   1.25   0.005    50           8.79e-4     151.6   3.0GB  <- first sub-1e-3 real-space cell
1e-4   2.5    0.01     100          2.87e-3      36.6   837MB
1e-4   5      0.02     200          8.84e-3       8.4   357MB

ALL CLAIMS CONFIRMED AT dt=0.025:
1. Cliff at fixed tau_add: 2e-4 row cliffs by tau=0.06, same as everywhere.
2. tau_add dominates at moderate tau: at tau=0.02, drop 2e-4 and 1e-4 give
   IDENTICAL rel (8.84e-3 both); drop re-emerges at low tau (tau=0.005:
   2.27e-3 vs 8.79e-4 - f-dominated).
3. Efficient window is dt-INVARIANT in K' = tau_add/drop: Pareto cells at
   K' 25-50, same as dt=0.1's ~10-50. In raw K the window shifted 4x down
   (K 0.6-1.25 vs 1-5), exactly the dt ratio.
4. Cleanest single-pair demonstration of the K failure mode: SAME (K=5,
   drop=1e-4) gives 1.9e-3 at dt=0.1 but 8.8e-3 at dt=0.025 - the fixed-K
   filter silently quadrupled tau_add.
Bonus: 1e-4/K1.25 = 8.79e-4 is the deepest real-space cell so far;
2e-4/K0.625 (2.27e-3, 34s, 881MB) Pareto-improves the old deep cell
(1.87e-3, 33.5s, 2.9GB) on RAM by 3.3x at comparable rel/wall.

## max_basis as the primary knob (2026-07-04, user hypothesis - CONFIRMED)

1. COLLAPSE: below the cliff, rel ~ peak_basis^(-0.2..-0.5) with only
   x1.1-1.35 residual scatter at fixed dt, across ALL (drop, tau_add)
   compositions. Basis size is the first-order accuracy predictor and the
   ~exact cost predictor (RSS and wall track peak). Second order: at equal
   peak, low-tau_add compositions are slightly better.
2. RANK CAP BEATS THRESHOLDS AT MATCHED SIZE (dt=0.025, drop=1e-4, K=1.25):
     M600k:  9.21e-4  127.9s   658MB   peak 547k
     M300k:  5.73e-4   58.3s   455MB   peak 289k  <- vs threshold cells at
                                          ~300k: 2.3-3.3e-3 (4-5x better!)
     M100k:  7.22e-4   18.6s   265MB   peak 98.6k <- NEW FRONTIER
     M30k:   9.90e-3    6.0s   218MB   peak 29.9k <- cap below natural
                                          weight distribution: back on the
                                          size-rel curve
     uncapped (1.3M):  8.79e-4  151.6s  3.0GB
   M100k beats the uncapped run on EVERY axis (8x wall, 11x RAM, better
   rel). keep-top-N is the optimal size-B selection; thresholds are not.
   The cap also bounds the leakage transient => much lower RSS per kept term.
   (M600k > M300k rel: non-monotone, cancellation noise ~1.5x - don't
   over-read small rel differences.)
3. ONCE THE CAP BINDS, THRESHOLDS ARE SECONDARY: drop 1e-5 vs 1e-4 at M100k:
   1.10e-3 vs 7.22e-4, same wall/RSS. Keep mild thresholds only to cheapen
   the pre-cap transient.
REVISED RECIPE (real space): set max_basis to the RAM/wall budget (the ONE
primary knob), keep tau_add ~ 0.005-0.02 and drop ~ 1e-4 as mild transient
hygiene. New best cells: 7.2e-4 @ 18.6s/265MB (M100k), 5.7e-4 @ 58s/455MB
(M300k) - both far beyond the old threshold-only frontier.

## Head-to-head: tau_add vs max_basis as the size-setter, matched peak (2026-07-04)

dt=0.025, drop=1e-4 unless noted; peak ~100k and ~300k classes:

  ~100k:  tau-only (K5, no cap):        8.84e-3    8.4s   357MB   82k
          cap-only (M100k, K0):         2.65e-3   18.8s   281MB   99k
          cap-only, drop 1e-5 (pure):   9.79e-4   22.0s   266MB  100k
          cap + tau (M100k, K1.25):     7.22e-4   18.6s   265MB   99k
  ~300k:  tau-only (K2.5, no cap):      2.87e-3   36.6s   837MB  325k
          cap-only (M300k, K0):         1.40e-3   63.7s   459MB  289k
          cap + tau (M300k, K1.25):     5.73e-4   58.3s   455MB  289k

VERDICT AT MATCHED SIZE:
1. ACCURACY: cap > tau_add (2-3.3x better rel as the size-setter), and cap
   composition + weak thresholds is best (7.2e-4 / 5.7e-4).
2. WALL: tau_add-only is ~2x cheaper than cap-only at the same size (K=0
   admits the full transient before the cap acts) - tau_add's remaining
   value is transient hygiene, not selection.
3. RSS: cap wins (bounds the transient): 459 vs 837MB at 300k.
4. SUBTLE: with the cap binding, drop=1e-4 HURTS accuracy vs drop=1e-5
   (2.65e-3 vs 9.79e-4 at K0/M100k): the per-step prune pre-deletes strings
   the rank cap would have ranked into the top-N. With a binding cap, keep
   thresholds MINIMAL (tiny drop, mild tau_add for wall only).
FINAL KNOB HIERARCHY: max_basis = the accuracy/cost dial; tau_add = wall
hygiene (~2x); drop = keep tiny once the cap binds.

## Framing note: fixed SIZE, not fixed SET (2026-07-04)

The cap-primary mode is "adaptive basis of constant rank": leakage still
proposes candidates every step (the expansion machinery is the engine), and
the rank cap selects the top-B of old+new on equal footing - membership
churns while size stays fixed. A truly frozen basis (static Galerkin, no
expansion) would fail once the operator front leaves the initial support.
Paper framing: max_basis is the Pauli-basis analog of the MPS bond dimension
chi - fixed-chi TDVP practice, with drop_tol as the (inferior) singular-
value-cutoff analog. Convergence protocol: rerun at 2B, compare.

## Follow-up: "cap + tau wins at 300k" was the drop confound (2026-07-04)

Completed 2x2 at M300k (dt=0.025):
                    drop=1e-4              drop=1e-5
  K=0 (no tau):     1.40e-3  63.7s 459MB   5.53e-4  74.4s 430MB  <- pure cap
  K=1.25:           5.73e-4  58.3s 455MB   9.73e-4  69.5s 429MB
(and M100k: K0/1e-5: 9.79e-4 22s; K1.25/1e-5: 1.10e-3 23s; K1.25/1e-4: 7.22e-4 19s)
-> With the drop=1e-4 confound removed, PURE CAP ties the best cap+tau cell
   (5.53e-4 vs 5.73e-4). All sane-threshold compositions at matched cap land
   in the same ~5-11e-4 band (the ~2x cancellation-noise band); tau_add's
   rel effect is noise-level, its wall effect ~10-25%. Hierarchy unchanged:
   cap primary; tau_add = modest wall trim; avoid mid-range drop (1e-4-ish)
   with a binding cap UNLESS paired with tau_add (the bad cell is
   specifically K0 + drop=1e-4: prune deletes candidates the cap would keep,
   and nothing pre-filters the transient either).

## CONCLUDING SUMMARY: truncation schemes (2026-07-04, user-approved framing)

How the threshold scheme works: two rules at different pipeline points.
The admission filter (add_tol = tau_add = K*drop_tol/dt, PPVM_K_LEAKAGE)
gates the LEAKAGE stage - newly proposed strings enter only if their
coefficient rate exceeds tau_add, keeping the doubly-enriched transient
lean. The pruning filter (drop_tol) deletes retained strings whose
coefficients fall below it after each step. The knobs decouple: error rises
slowly (~tau_add^0.24) up to a sharp cliff at tau_add ~ 0.04 that is
INDEPENDENT of dt and drop_tol (tau_add, not K, is the natural parameter);
drop_tol controls the retained tail and hence basis size. Both are indirect
ways of setting one quantity - error and cost collapse (~30%) onto a single
function of peak basis size.

Why max_basis supersedes it: at MATCHED basis size the rank cap (keep-top-B
over retained+proposed, membership still churning every step - the fixed-
bond-dimension TDVP analog) yields 2-5x lower error than any threshold
combination, and bounds the transient RAM thresholds never touch. Best cap
cells beat best threshold cells on every axis simultaneously. Thresholds on
top of the cap trim wall by only ~10-25%, their accuracy effect is within
the cancellation-noise band, and mid-range drop_tol with a binding cap can
actively hurt. RECOMMENDATION: run with max_basis alone (thresholds ~0);
B is the single convergence dial, verified by re-running at 2B.

## Correction to the MPS analogy (2026-07-04, user)

"Fixed-chi TDVP" was imprecise: the cap-primary scheme grows the basis
freely from the seed until it reaches B, then saturates (peak == cap) with
churning membership. The right analogy is two-site TDVP/TEBD with a MAXIMUM
bond dimension and negligible singular-value cutoff (B = chi_max), not
fixed-chi single-site TDVP (a fixed manifold, no growth). Notes updated
(sec:rank-cap).
