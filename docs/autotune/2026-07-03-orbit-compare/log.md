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

## Real-space MSD: Trotter vs CAP-PRIMARY pec (2026-07-04, final comparison)

pec grid, pure cap (drop=1e-5, K=0), dt x B:
  dt     B      median_rel  wall_s  RSS
  0.1    10k    7.65e-3      0.3    163MB
  0.1    30k    2.49e-2      0.8    218MB   <- cancellation outlier
  0.1    100k   5.23e-3      3.1    250MB
  0.1    300k   2.66e-3     10.5    428MB
  0.05   10k    4.68e-2      0.6    179MB
  0.05   30k    1.22e-2      1.7    223MB
  0.05   100k   1.50e-3      5.6    283MB   <- sweet spot
  0.05   300k   8.58e-4     18.7    428MB
  0.025  10k    1.37e-2      1.0    170MB
  0.025  30k    1.02e-2      3.1    233MB
  0.025  100k   9.79e-4     22.0    266MB
  0.025  300k   5.53e-4     74.4    430MB
Cap mode moves pec's dt optimum from 0.1 (threshold mode) to 0.05.

FIXED-PRECISION HEAD-TO-HEAD vs Trotter (best measured per class):
  ~1e-2:   pec 0.1/M10k    7.7e-3   0.3s  163MB | trot 0.025/1e-4  7.9e-3  13.7s 134MB
           -> pec 46x faster at ~RAM parity (both near python baseline)
  ~2e-3:   pec 0.05/M100k  1.50e-3  5.6s  283MB | trot 0.025/3e-5  2.07e-3 209s  415MB
           -> pec 37x faster AND 1.5x lighter, better rel
  <1e-3:   pec only: 8.6e-4 @ 18.7s/428MB; 5.5e-4 @ 74s/430MB
THE SPLIT VERDICT DISSOLVES: with max_basis as the knob, pec dominates
Trotter in real space on wall at every precision (37-46x) and reaches RAM
parity-or-better; Trotter's former RAM edge was an artifact of comparing
against threshold-tuned pec (whose transient the thresholds never bounded).

## CORRECTION: "pure cap" cells were drop=1e-5, and that matters (2026-07-04, user)

True drop=0 cells (K=0), vs the mislabeled "pure cap" (drop=1e-5) cells:
  dt=0.025 M100k: drop0 5.90e-3 10.3s | 1e-5: 9.79e-4 22.0s   (6x worse)
  dt=0.025 M300k: drop0 1.50e-3 34.5s | 1e-5: 5.53e-4 74.4s   (2.7x)
  dt=0.05  M100k: drop0 2.64e-3  5.3s | 1e-5: 1.50e-3  5.6s   (1.8x)
  dt=0.05  M300k: drop0 1.35e-3 17.6s | 1e-5: 8.58e-4 18.7s   (1.6x)
Consistent across 4 cells -> real effect, not cancellation noise. With
drop=0 the peak sits exactly at B (no clearing); with 1e-5 slightly under.
Hypothesized mechanism: a tiny prune clears stagnant near-zero weight from
the retained set, freeing rank slots for freshly proposed strings; with
drop=0 a zombie tail occupies the bottom of the top-B. drop=0 is somewhat
faster (up to 2x at dt=0.025) but strictly Pareto-dominated on accuracy.
REVISED RECOMMENDATION: cap + a pruning threshold FAR BELOW the working
coefficient scale (e.g. 1e-5 here) - "near zero", not exactly zero.
Cap-primary comparison table (previous section) used the 1e-5 cells and is
unaffected; only the "pure cap" label was wrong.

## drop sweep at fixed cap: 0 / 1e-10 / 1e-6 / 1e-5 / 1e-4 (2026-07-04)

dt=0.025 M100k: 0: 5.90e-3 | 1e-10: 5.22e-3 | 1e-6: 6.72e-4 | 1e-5: 9.79e-4 | 1e-4: 2.65e-3
dt=0.025 M300k: 0: 1.50e-3 | 1e-10: 1.85e-3 | 1e-6: 2.57e-3 | 1e-5: 5.53e-4 | 1e-4: 1.40e-3
dt=0.05  M100k: 0: 2.64e-3 | 1e-10: 4.53e-3 | 1e-6: 2.48e-3 | 1e-5: 1.50e-3
REVISION of the previous section's claim:
- ROBUST: 1e-10 == 0 in every respect (peak pinned at B; similar rel): the
  beneficial prune must act at a scale comparable to the bottom-of-cap
  coefficients, not at numerical dust. drop in [1e-6,1e-5] is never worse
  than 0 beyond noise and sometimes much better.
- CLEAN at M100k/dt=0.025: {0,1e-10} ~5-6e-3 vs {1e-6,1e-5} ~0.7-1e-3
  (5-8x, beyond noise). MUDDY at M300k: non-monotonic (1e-6 worse than 0,
  1e-5 best), spread ~4.7x ~ noise x2 - the "consistent 1.6-6x" claim was
  overconfident there.
- Recommendation unchanged in practice (small nonzero prune ~1e-6..1e-5),
  but the strength of the effect is cap-size- and cell-dependent; the
  zombie-clearing mechanism is supported at B=100k and unproven at B=300k.

## MECHANISM RESOLVED: the cap acts at ADMISSION; drop=0 freezes the basis (2026-07-04)

Code (lib.rs): add_leakage_capped admits new strings only into
room = max_basis - len(basis), top leakage magnitudes first; retained
strings are NEVER displaced (cap_basis is a no-op in the normal flow since
admission never overfills). Therefore:
- drop ~ 0 (incl. 1e-10): basis fills to exactly B early, room -> 0,
  NO admission ever again -> FROZEN static basis for the rest of the run.
  Explains: peak == B exactly, lower wall (enrichment idle), and the ~9x
  accuracy loss (5.9e-3 vs 6.7e-4 at B=100k) - static-Galerkin error.
- drop in [1e-6, 1e-5]: the prune clears the low-weight tail each step,
  opening room; admission refills with the best candidates -> membership
  churns at fixed size. This is the actual algorithm behind all good cap
  cells. The prune is the CHURN VALVE, not "hygiene".
- drop=1e-4: valve too wide, deletes strings still carrying weight.
This supersedes the "zombie-clearing" hypothesis and explains the entire
drop-sweep table, including the wall differences. Notes sec:rank-cap
corrected (the previous "top-B of union / displacement" description was
factually wrong about the implementation).
POSSIBLE CODE IMPROVEMENT (future): true top-B-of-union displacement
(rank competition between retained and proposed) would remove the need for
the pruning threshold entirely - worth implementing and testing.

## Top-B-of-union displacement scheme: implemented, NEGATIVE result (2026-07-04)

Implemented admit_basis (A) on pc_step/pc_step_arr: enrichment may grow the
working set to A >= B; the final cap_basis (now live) keeps the top-B by
evolved |coeff| over the whole union - genuine rank displacement, drop_tol
no longer needed for turnover. Backward compatible (admit_basis=None = old
behaviour; regression cell M100k/1e-5 reproduces 9.79e-4 exactly). Tests pass.

Results (dt=0.025, drop=0):
  B=100k: A=1.25B: 3.74e-3 31s | A=2B: 2.96e-3 46s | A=4B: 1.38e-3 89s
          vs valve (cap + drop 1e-5): 9.79e-4 23s/288MB  <- VALVE WINS
  B=300k: A=2B: 5.51e-4 157s/583MB vs valve 5.53e-4 74s/430MB  <- tie, 2x wall
  dt=0.05 B=100k: A=2B: 1.99e-3 21s vs valve 1.50e-3 5.6s
Monotone improvement with A (converging toward the valve result from above).

MECHANISM HYPOTHESIS (boundary-layer cycling): the union scheme deletes
A-B strings EVERY step, most of them front-layer strings that get
re-admitted next step, re-accumulate a small coefficient, and are deleted
again - each cycle discards accumulated weight (repeated coherent error
injection) and wastes work. The valve scheme's turnover (~1e2/step at these
knobs) makes near-permanent swaps instead. Lesson: LOW-churn evolution
accumulates less truncation error than aggressive per-step re-selection;
the 2TDVP analogy breaks because Pauli-dictionary "Schmidt vectors" can't
rotate - discrete swap-in/out has a per-swap error cost that SVD rotation
does not. The valve (cap + small drop) stays the recommended scheme.

## Late-time / max-rel analysis: valve vs displacement (2026-07-04, user suspicion CONFIRMED)

max pointwise rel over t=0.2..2.0 (same cells as previous section):
  0.025 100k valve 1e-5:   median 9.8e-4  MAX 3.80e-2 @t=2.0  (profile: 1e-5-ish early, then 9e-3,2e-2,3e-2,4e-2)
  0.025 100k disp A=1.25B: median 3.7e-3  MAX 1.13e-2 @t=2.0
  0.025 100k disp A=2B:    median 3.0e-3  MAX 8.1e-3  @t=1.6
  0.025 100k disp A=4B:    median 1.4e-3  MAX 6.3e-3  @t=1.4  (rel(T)=1.0e-3)
  0.025 100k frozen:       median 5.9e-3  MAX 5.3e-1  @t=2.0  (catastrophic)
  0.025 300k valve:        median 5.5e-4  MAX 1.05e-2 @t=2.0
  0.025 300k disp A=2B:    median 5.5e-4  MAX 1.30e-2 @t=2.0  (tie)
  0.05  100k valve:        median 1.5e-3  MAX 3.04e-2 @t=1.8
  0.05  100k disp A=2B:    median 2.0e-3  MAX 1.20e-2 @t=1.8  (2.5x better)

READING: the valve's low churn (~1e2 swaps/step) is TOO SLOW to track the
operator front at late times - it is quasi-frozen on the front's timescale,
and its error diverges toward t=T (mini version of the frozen catastrophe;
bigger B delays it - B=300k valve stays ok). The displacement scheme
refreshes the basis every step and holds a flat late-time error (max 3-6x
better than the valve at B=100k), at the cost of a small early/mid plateau
(~2e-3, the cycling injection) and ~2x wall.
VERDICT IS METRIC-DEPENDENT:
  median rel (whole-curve typical) -> valve wins;
  max rel / late-time (what matters for transport-coefficient extraction
  from late-time MSD slopes) -> displacement wins at B=100k, ties at 300k.
Neither scheme dominates. For D(gamma) production runs (late-time slopes),
prefer displacement (admit_basis ~ 2-4x B, drop=0) or valve with generous B.

## T=10 comparison, every-step reference (2026-07-05)

Dense exact reference: exact_msd_L7_T10.npz (M=400, every dt=0.025 step,
12 min sector-reduced ED). Driver now takes SCAN_T env + max_rel column.
dt=0.025, metrics over ALL 400 points (t=0.025..10):

  scheme                med_rel   max_rel   wall    RSS
  frozen  B=100k        3.17e-1   5.98e-1    53s   238MB   <- permanent ~0.3-0.5 plateau
  valve   B=100k        8.67e-3   4.09e-2    66s   279MB
  disp A=2B B=100k      3.88e-3   1.44e-2   177s   341MB   <- best B=100k
  disp A=4B B=100k      4.67e-3   1.55e-2   340s   442MB
  valve   B=300k        5.95e-3   1.53e-2   302s   460MB
  disp A=2B B=300k      5.71e-3   1.44e-2   669s   625MB

Time-resolved picture (rel_vs_time_T10.png): the valve B=100k error PEAKS
~4e-2 at t~2 (the quasi-freeze episode seen in the T=2 study) and then
RECOVERS to the common ~5e-3-1e-2 saturation band - once the operator
equilibrates the slow churn suffices again; it does not keep diverging.
All adaptive schemes converge to that band; disp A=2B/B=100k has the lowest
long-time curve (~2-6e-3). Frozen never recovers.
T=10 VERDICT (B=100k): displacement wins BOTH median (2.2x) and max (2.8x)
at 2.7x wall; disp A=2B/B=100k even beats valve/B=300k on both error metrics
at 0.6x wall and 0.74x RSS. At B=300k valve==disp within noise, valve 2.2x
faster. A=4B adds nothing over A=2B at T=10.
Practical: for long-horizon runs at tight memory, displacement A=2B is the
best configuration measured; the valve's weakness is specifically the
transient regime around the equilibration time.

## Frozen-basis B scan at T=10 (2026-07-05)

  frozen B=100k: median 3.17e-1  max 5.98e-1   53s  238MB
  frozen B=200k: median 2.10e-1  max 3.30e-1  111s  296MB
  frozen B=300k: median 3.87e-2  max 4.69e-1  180s  397MB
Raw curves (msd_raw_T10.png): B=100k undershoots the plateau (~2.4-2.9),
B=200k lands ~3.0-3.3, B=300k OVERSHOOTS through the wrap transient (peak
5.6 at t~2.3!) then wanders around ~3.7-4.1 with persistent oscillations
(max rel still 0.47; long-time median 3.9e-2 = 5-10x above any adaptive
scheme at a THIRD of the basis budget of valve/300k... equal budget).
Frozen error is NOT monotone-decreasing in B in the transient (300k's max
~ 100k's), and no frozen size equilibrates to the correct plateau cleanly.
Adaptivity - even the valve's ~1e2 swaps/step - is qualitatively load-
bearing; basis size cannot substitute for it.

## Valve drop=1e-6 at T=10 + combined figure (2026-07-05)

  valve 1e-6 B=100k: median 1.39e-2  max 6.96e-2   59s  271MB
  valve 1e-6 B=300k: median 7.67e-3  max 3.24e-2  207s  432MB
Both WORSE than valve 1e-5 at T=10 (8.67e-3/4.09e-2 and 5.95e-3/1.53e-2):
smaller drop = fewer swaps/step = deeper quasi-freeze through the wrap
transient. (At T=2 1e-6 had looked better - horizon-dependent ranking.)
The zoom panel (msd_combined_T10.png, 3.7-4.1) shows valve 1e-6 drifting
3.71-3.95 around the plateau out to t=10 while disp A=2B and valve 1e-5
track the exact curve's fine wiggles at the ~1e-2 level.
Combined 3-panel figure (raw / plateau zoom / rel err): msd_combined_T10.png.

## Threshold-only runs added to the T=10 comparison (2026-07-05)

  thresholds tau=0.02 (1e-4/K5):   median 2.48e-2  max 4.19e-2   69s   519MB  peak 142k
  thresholds tau=0.01 (1e-4/K2.5): median 2.15e-2  max 3.29e-2  376s  1014MB  peak 475k
(T=2 peaks were 82k/325k: uncapped bases DRIFT 1.5-1.7x by T=10.)
vs cap schemes at SMALLER bases: valve 1e-5/100k 8.7e-3/4.1e-2 @66s/279MB;
disp A=2B/100k 3.9e-3/1.4e-2 @177s/341MB -> thresholds are 3-6x worse on
median at 2-3x the RSS.
Time-resolved (msd_combined_T10.png): thresholds fail at BOTH ends -
(a) EARLY: rel ~1e-2 already at t<0.5 (the absolute-rate admission filter
bites hardest during the initial fast spread, where cap schemes are at
1e-5); (b) LATE: both threshold curves FLATLINE at a wrong plateau
(4.07/4.09 vs exact ~3.96-3.98) from t~4.5 on - the truncated dynamics
reaches a stagnant fixed point ~2-3% high, error frozen at 2-3e-2.
FIGURE msd_combined_T10.png is now the complete truncation-scheme taxonomy:
frozen (catastrophic), thresholds (uncontrolled size + early/late bias),
valve (transient episode, drop-sensitive), displacement (uniform tracking).

## Is the T=10 saturation a dt artefact? (2026-07-05, user question)

dt scan at fixed configs (T=10, B=100k):
  disp A=2B: dt=0.05: 6.20e-3/1.54e-2 | 0.025: 3.88e-3/1.44e-2 | 0.0125: 1.86e-2/3.06e-2
  valve 1e-5: 0.025: 8.67e-3/4.09e-2 | 0.0125: 4.09e-3/2.76e-2
  valve 1e-5 B=300k, dt=0.0125: 5.38e-3/1.54e-2 (435s/459MB)
ANSWER: partially, and scheme-dependently.
1. VALVE was dt-limited at 0.025: halving dt halves its median (8.7->4.1e-3)
   and shrinks its transient peak. Its churn/time also doubles (confounded
   with dt error - drop is per-step).
2. DISPLACEMENT is ANTI-dt-limited: its injection is per-STEP (A-B strings
   cycled each step), so error/time ~ 1/dt; interior optimum dt ~ 0.025.
   The universal band at dt=0.025 was the two mechanisms crossing.
3. Once dt-unblocked, B still does NOT improve the valve median
   (100k: 4.09e-3 vs 300k: 5.38e-3 at dt=0.0125) though max improves 1.8x.
   A persistent post-transient offset ~4-6e-3 dominates the T=10 median:
   damage done during the wrap transient persists in the equilibrated state,
   and none of (B, A, drop, dt) measured so far pushes the median below
   ~4e-3. OPEN: dt=0.00625 valve run would test whether the residual is
   still dt-limited (confounded with churn/time doubling); alternatively
   the transient-episode error itself must be attacked (displacement
   at its optimal dt already does: max 1.44e-2).
Best T=10 cells so far: valve 1e-5/100k/dt=0.0125 (4.09e-3 med, 131s/276MB)
and disp A=2B/100k/dt=0.025 (3.88e-3 med, max 1.44e-2, 177s/341MB).

## Displacement B=300k dt scan (2026-07-06)

  disp A=2B B=300k: dt=0.05: 4.70e-3/1.38e-2 280s | 0.025: 5.71e-3/1.44e-2 669s | 0.0125: 8.61e-3/1.78e-2 1048s
Anti-dt trend confirmed at 300k but much WEAKER than at 100k (x1.8 over the
dt range vs x5): the per-step cycling injection scales with the
boundary/basis ratio, so bigger B damps it. Displacement dt-optimum shifts
to LARGER dt at larger B (best 300k cell: dt=0.05, 4.70e-3 at 280s).
B-scaling remains weak in all schemes at T=10 (transient scar dominates).
Figure rel_vs_time_dtscan_T10.png updated with B=300k (dashed).

## Late-time floor DECOMPOSED: revivals + scar (2026-07-06, user question)

Exact plateau t>=4: mean 3.9738, range [3.9405, 4.0021]. A CONSTANT
prediction has median rel 1.84e-3, max 8.5e-3 - the size of the coherent
finite-size revival signal itself.
Scheme late-time errors: 4.0-5.6e-3 vs exact; 3.0-6.0e-3 vs the SMOOTHED
exact; correlation with the exact wiggles: -0.27..+0.28 ~ ZERO.
=> The late-time floor = (a) the revival component (~2e-3 median / ~8e-3
max): coherent recurrences living in high-weight scrambled strings that NO
truncated basis retains at any feasible budget - knob-independent by
nature; + (b) a residual plateau bias ~3e-3 (the transient scar proper),
improvable only by surviving the wrap transient better.
CONSEQUENCES: (i) the apparent "non-convergence" at late times is mostly a
metric artifact - we were grading truncated methods on reproducing
revivals; (ii) for transport physics (slopes/plateaus, thermodynamic
limit) the relevant late-time accuracy is ~3e-3 and the revival component
should be excluded (smoothed reference) or pushed out by larger L;
(iii) this cleanly documents the method's fundamental scope: hydrodynamic
/ dissipative components YES, coherent recurrences NO - worth a sentence
in the paper's discussion.

## Revival-carrier hypothesis TESTED with exact numerics at N=10 (2026-07-06)

Setup (revival_test_N10.py): L=5 ladder, O(0)=center-rung Z pair, exact
evolution; full 4^10 Pauli decomposition (per-qubit transform, round-trip
verified to 3e-16, Parseval exact).
Coefficient structure at t0=2 (CONFIRMS the dust picture):
  260,830 strings nonzero; top-100 hold 21% of the norm; 90% needs
  B=100k (=10% of all strings); kept mean Pauli weight 4.2-7.4 vs dust 7.5.
SINGLE cut at t0=2, then exact evolution (REFUTES the strong one-shot form):
  B=1k:   future wiggle corr +0.43, amplitude 20%
  B=10k:  corr +0.85, amplitude 43%
  B=100k: corr +0.99, amplitude 94%   <- one cut at 10%-budget keeps revivals
REPEATED cut every dt=0.1 (100 cuts), exact evolution between cuts
(IDENTIFIES the real mechanism):
  B=10k:  corr -0.14, amplitude 5%
  B=100k: corr -0.09, amplitude 22% (uncorrelated -> noise, not signal)
  and the late-time mean shifts to the exactly-uniform value (2.000/1.993
  vs exact 1.9386) - the flat, slightly-off plateau seen in production.
CONCLUSION: truncation acts as WEAK CONTINUOUS DEPHASING on the coherent
recurrence channel. Per-cut amplitude retention is high (~0.94-0.99 at
these budgets) but compounds over hundreds of steps to ~zero; the same
budget that survives one cut erases revivals under repetition. This also
explains part of the plateau bias (truncation pushes the profile toward
exactly-uniform, losing the coherent depression of the exact plateau).
Paper angle: implicit truncation-dephasing is the flip side of DAOE's
explicit artificial dissipation; hydrodynamics is immune (carried by few
large-coefficient strings re-selected every step), recurrences are not.

## Trotter + truncation at T=10: same dephasing, plus its own wiggles (2026-07-06)

  trotter mac=1e-4 (peak 111k, STATIONARY == its T=2 peak): med 8.60e-3 max 3.08e-2 37s/134MB
  trotter mac=3e-5 (peak 1.20M): med 3.88e-3 max 1.68e-2 623s/415MB
Late-window (t>=4) wiggle analysis vs exact (exact wiggle std 0.0141):
  scheme                    late mean   corr    amp ratio
  Trotter 1e-4 (111k)        4.0042    +0.60      1.92   <- OVERSHOOT: own coherent wiggles
  Trotter 3e-5 (1.2M)        3.9932    -0.22      0.59
  valve B=100k               4.0058    +0.46      0.67
  disp A=2B B=100k           3.9889    +0.21      0.81
  valve B=300k               4.0067    +0.30      0.65
  disp A=2B B=300k           3.9984    +0.24      0.25
(exact late mean 3.9738 - ALL schemes sit 0.4-0.8% high.)
ANSWER: yes - Trotter+truncation shows the same effective dephasing of the
true recurrences (no scheme tracks the exact dip structure, e.g. the t=5.7
excursion to 3.945; at tight mac Trotter flattens exactly like the pec
schemes). TWIST: at loose mac Trotter additionally GENERATES its own
coherent oscillations (amplitude 1.9x the exact wiggles, partially phase-
aligned corr 0.6) - spurious revival-like artifacts from the Trotter+
per-gate-truncation dynamics, absent in the exact-in-dt pec schemes.
REFINEMENT of the earlier corr~0 claim: correlations across schemes span
-0.2..+0.6 - the N=14 late-time signal contains a partially retrievable
low-weight coherent component on top of the erased scrambled part; the
dephasing strength ranks with discarded weight per event (disp most
aggressive: amp 0.25; valve intermediate; loose Trotter least, but noisy).
Figure: msd_trotter_overlay_T10.png.

## Comparison restricted to t <= MSD maximum (t=1.225) (2026-07-06, user)

Exact peak MSD 4.5080 at t=1.225. Median/max rel over t in [0.2, 1.225]
(walls/RSS are full-T=10 costs; the up-to-peak portion is ~12% of steps):
  valve 1e-6 B=100k    3.5e-5 / 1.6e-3    59s   271MB   <- best & cheapest pec
  disp A=2B B=100k     3.6e-5 / 4.7e-3   177s   341MB
  disp A=4B B=100k     4.5e-5 / 5.0e-3   340s   442MB
  valve 1e-6 B=300k    6.4e-5 / 5.1e-3   207s   432MB
  disp A=2B B=300k     6.7e-5 / 1.3e-3   669s   625MB
  frozen B=300k        8.2e-5 / 2.1e-3   180s   397MB
  frozen B=200k        9.6e-5 / 7.6e-3   111s   296MB
  valve 1e-5 B=300k    1.1e-4 / 2.4e-3   302s   460MB
  frozen B=100k        1.5e-4 / 8.1e-3    53s   238MB
  valve 1e-5 B=100k    1.9e-4 / 1.7e-3    66s   279MB
  thresholds tau=.01   1.1e-3 / 4.1e-3   376s  1014MB
  Trotter mac=3e-5     1.4e-3 / 3.5e-3   623s   415MB
  thresholds tau=.02   2.0e-3 / 1.2e-2    69s   519MB
  Trotter mac=1e-4     2.7e-3 / 1.0e-2    37s   134MB
FINDINGS:
1. In the ballistic/transport window ALL cap-based pec variants (valve,
   displacement, even frozen!) collapse to median 3.5e-5..1.9e-4 - the
   scheme distinctions that consumed the T=10 analysis are IRRELEVANT
   before the wrap; even freezing only starts to bite at the window edge.
2. pec is 10-80x more accurate than Trotter/thresholds here: Trotter is
   floored by O(dt^2) splitting error during the fast ballistic dynamics
   (1.4e-3 even at 1.2M strings); thresholds by early-time tau_add bias.
3. USER'S RED-HERRING CALL VALIDATED: the late-time scheme distinctions
   were driven by finite-size recurrence physics outside the method's
   scope. In the physically extractable window the story is simple:
   cap-based CTPP at modest B dominates everything, and the cheapest cell
   (valve 1e-6, B=100k: 3.5e-5 at 59s/271MB full-run cost) wins outright.
For production D extraction at large L (window entirely pre-wrap), this
is the regime that matters.

## L=21 (N=42) MSD to the plateau: rough -> converged (2026-07-06)

No exact reference; convergence by B/dt-doubling + cross-method agreement.
Driver: SCAN_L env. Physics: NO overshoot peak at L=21 - the MSD rises
monotonically to the uniform plateau 36.67 (reaches 36.47 by t=15; the
L=7-style coherent overshoot is a small-L effect). The 'peak' window is
the full rise t in [0, ~15].
ROUGH PASS pitfalls (documented): B<=300k valve cells scatter WILDLY at
L=21 (non-monotone in B; pec100k collapsed to MSD~20 at t=6); loose
trotter (mac 1e-3, 10k strings) underestimates MSD by 40% - the ballistic
front carries the dj^2 weight and is exactly what truncation cuts.
CONVERGED (displacement A=2B, dt=0.1): B=150k vs 300k median 1.7e-3
(max 7e-3) over the FULL window; dt 0.1 vs 0.05 median 1.3e-3. Trotter
ladder approaches the same curve from below (mac 3e-4 -> 1e-4: 28.4 ->
30.9 at t=6 vs disp 31.6).
FINAL TABLE (rel vs disp300k over [0.2,15]; wall/RSS for T=15 dt=0.1):
  disp A=2B B=300k  (reference)             482s   666MB
  disp A=2B B=150k   1.7e-3 / 7.3e-3        145s   420MB  <- converged, cheap
  valve 1e-6 B=300k  3.5e-3 / 1.2e-1(!)     156s   443MB  <- erratic late excursions
  Trotter mac=1e-4   2.5e-2 / 3.8e-2         75s   344MB  <- 2.5% LOW (front loss),
     saturates at 35.09 vs 36.47; mac-convergence is x8 basis / x10 wall per
     x2 error - reaching 2e-3 would need ~1e-5..3e-6 mac = 10-100M strings,
     INFEASIBLE. At N=42 CTPP(displacement) converges where Trotter cannot.
Figure: msd_L21_T15.png.

## L=21 Trotter dt x mac ladder, t<=6, vs converged displacement (2026-07-06)

  run                     med_rel   max_rel  MSD(6)  wall    rss~   peak
  trot 0.1/3e-4           5.3e-2    1.0e-1   28.37    11s   143MB   102k
  trot 0.1/1e-4           1.4e-2    2.3e-2   30.91   151s   344MB   854k
  trot 0.05/3e-4          6.4e-2    1.5e-1   30.67     8s   120MB    40k
  trot 0.05/1e-4          3.8e-2    7.4e-2   30.16    67s   187MB   343k
  trot 0.05/3e-5          7.4e-3    1.4e-2   31.36  1175s   1.1GB  3.67M  <- best trotter
  trot 0.025/1e-4         7.9e-2    1.4e-1   28.45    62s   140MB   115k
  disp B=150k (check)     2.8e-3    7.3e-3   31.56    78s   375MB   150k
  (reference: disp B=300k; 0.025/3e-5 + 0.025/1e-5 cells pending/partial)
FINDINGS:
1. At fixed mac, SMALLER dt is WORSE for Trotter (0.025/1e-4 = 7.9e-2 vs
   0.1/1e-4 = 1.4e-2; basis shrinks 854k->115k from per-gate culling) -
   the joint dt<->mac coupling at L=21 scale. Only mac-tightening converges.
2. Best Trotter rung (0.05/3e-5, 3.67M strings, 20 min, 1.1GB) is still
   2.6x LESS accurate than the CHEAP displacement run (150k strings, 78s,
   375MB) - a 15x wall / 3x RAM / 24x string handicap AND worse error.
3. Extrapolating the ladder (x10 wall per ~x2 error): matching disp-150k's
   2.8e-3 needs mac ~1e-5..3e-6 at dt~0.05 = 10-40M strings / hours /
   several GB. Trotter mac-convergence at N=42 is impractical; the
   displacement-CTPP advantage at scale is decisive on every axis.
Figure: msd_L21_trotter_ladder.png.

## L=41 (N=82) vs external Pauli-propagation reference (2026-07-06)

User-provided reference figure: MSD vs time, truncation ladder 2^-13..2^-18
(PauliPropagation-style, RvKP XX ladder; their time axis = 4x ours). Ran
disp A=2B at L=41, T=5 (our units), dt=0.1:
  B=150k: 100s / 393MB ;  B=300k: 209s / 650MB
  B-convergence: median 8.5e-3, max 1.8e-2 (150k vs 300k) - slightly looser
  than at L=21 (front support is 2x wider), fine at digitization accuracy.
Comparison at digitized brown (2^-18) points (t_theirs/4):
  t=2.0: 12.76 vs 12.6 (+1.3%) | t=2.5: 16.43 vs 16.5 (-0.4%)
  t=3.0: 20.14 vs 21.0 (-4.1%) | t=3.75: 26.11 vs 26.0 (+0.4%)
  t=4.0: 27.60 vs 28.0 (-1.4%) | t=5.0: 35.25 vs 35.0 (+0.7%)
  (t<1.5 points off by 6-11% - digitization error dominates there, small
   values on a coarse plot.)
=> Agreement within ~1-4% everywhere = digitization accuracy. EXTERNAL
CROSS-VALIDATION of the 82-qubit displacement-CTPP result against an
independent code/method, at 209s/650MB on a laptop. The L=21 curve peels
off at t>3 exactly as finite size predicts (fronts wrap); L=41 tracks the
unconstrained growth to t=5. Time-convention factor 4 confirmed by the
match itself. Figure: msd_L41_vs_reference.png.

## L=41 external comparison FINALIZED with accurately digitized data (2026-07-06)

History: my eyeballed digitization was wrong; the user's accurate data
first appeared mismatched (ours ~25% low) - a free one-parameter time-
rescale fit found lambda=1.252, and the user then identified the cause:
their exported x-axis needed *5/4 (fit matched the correction exactly,
a good pipeline sanity check). Corrected comparison (62 ref points):
  disp B=300k vs 2^-18: median 5.5e-3  max 1.7e-2 ; endpoint (t=4.95):
    34.86 vs 34.63 (+0.7%)
  disp B=150k vs 2^-18: median 1.5e-2 (slightly low at late t)
  disp B=300k vs 2^-17: median 1.1e-2, growing to 3-7e-2 at t>4 - the
    reference's own truncation ladder bending down, as expected.
=> 82-qubit displacement-CTPP (209s/650MB laptop) agrees with the external
Pauli-propagation 2^-18 curve at the ~0.5% level across the full window,
and the deviation vs their looser 2^-17 curve has the right sign/shape.
Strong independent cross-validation. reference_digitized_L41.py holds the
data; figure msd_L41_vs_reference.png.

## Reference identified + their resources (2026-07-06)

The external reference = Begusic & Chan, PRX Quantum 6, 020302 (2025)
"Real-Time Operator Evolution in 2D and 3D via Sparse Pauli Dynamics",
Fig. 2(c): XX-ladder [their Eq. (17), H = (1/4)*ours -> t_theirs = 4*t_ours,
confirming our conversion], L=41, n=82, q_j=(Z_j1+Z_j2)/2 - IDENTICAL setup.
Their D ~ 0.94 (their units) from regression t in [10,20] (= ours [2.5,5]).
THEIR COST for the 2^-18 curve (p.4): "the two points in Fig. 2(d) at fixed
delta/dt = 2^-18/0.02 take around 84 h (dt=0.01) and 43 h (dt=0.02) to
simulate on six cores" (Xeon Platinum 8352Y 2.2GHz). RAM for Fig 2 not
stated; their stated memory elsewhere: 2D delta=2^-23: 8.5e9 Pauli ops,
">1 TB", 36h/16 cores; 3D: "memory budget of about 1.5 TB" at <1e9 ops
(~125 B/op storage). Their ladder 2^-18 run's N is not given (bound
N <= c0/delta^2 ~ 3.4e10; plausibly 1e8-1e9 ops -> tens of GB, estimate).
OURS at matched accuracy (0.5% median agreement): disp A=2B B=300k,
dt=0.1 (our units; 20x coarser step than theirs - exact-in-dt makes this
possible), 209 s / 650 MB on a laptop.
=> wall ratio ~740x (43h vs 209s); core-hours ~445x (258 vs 0.58);
RAM: 650 MB vs (unstated, plausibly 10s of GB). Their delta/dt coupling
(their Sec II: error is a function of delta/dt) forces small dt AND large
N; CTPP's dt-freedom is the structural advantage on display.

## L=21 Trotter deep rungs completed (2026-07-06, background cells landed)

  trot 0.025/3e-5: median 1.51e-2 max 2.7e-2  MSD(6)=31.18  1385s  350MB   1.25M strings
  trot 0.025/1e-5: median 2.86e-3 max 8.9e-3  MSD(6)=31.35  14845s 2.0GB  11.08M strings
The 1e-5 rung FINALLY matches disp B=150k accuracy (2.8e-3): the
extrapolated cost was right. MATCHED-ACCURACY HEADLINE AT L=21 (t<=6):
  CTPP-displacement: 78 s, 150k strings, ~0.4 GB
  Trotter (2nd-order): 4.1 h, 11.1M strings, 2.0 GB
  -> 190x wall, 74x strings, ~5x RAM at equal accuracy.
Also note trot 0.025/3e-5 (1.5e-2) is WORSE than 0.05/3e-5 (7.4e-3):
the dt-mac coupling again - halving dt at fixed mac degrades.

## A = K*B scan at L=41, B=300k, T=5 (2026-07-07)

  K=1:   29s  480MB  med 5.9e-2 vs K=4  MSD(5)=20.50  <- FROZEN (no headroom
         with drop=0: admission stops at fill; curve bends over at t~2.5)
  K=1.25 50s  519MB  1.8e-2   33.51
  K=1.5  67s  548MB  9.4e-3   34.78
  K=2   209s  650MB  3.5e-3   35.25   (yesterday's run; wall possibly inflated)
  K=3   128s  806MB  4.4e-4   35.42   <- converged (vs K=4)
  K=4   168s 1041MB  anchor   35.46
Agreement vs Begusic-Chan 2^-18 saturates ~5e-3 for K>=2 (digitization/their
truncation limited). RSS ~ linear in A; wall non-monotone (K=3 cheaper than
K=2). RECIPE UPDATE for large L: A=3B preferred (10x tighter than 2B at ~no
wall cost, +25% RAM); A=2B when RAM-bound; A=B is the frozen limit - never.
Figure: msd_L41_Kscan.png.

## A = K*B scan at L=41, B=150k (2026-07-07)

  K=1:   13s 292MB  med 3.3e-1 vs K=4  MSD(5)= 9.48  <- frozen, collapses even
         earlier than at 300k (fills sooner)
  K=1.25 22s 324MB  3.5e-2  33.50
  K=1.5  27s 327MB  2.3e-2  33.48
  K=2   100s 393MB  7.3e-3  34.62   (2026-07-06 wall; possibly inflated)
  K=3    52s 487MB  1.3e-3  35.38   <- converged in K
  K=4    66s 580MB  anchor  35.55
Residual B-bias (150k/K=4 vs 300k/K=4): median 2.3e-3, max 5.3e-3 - the
irreducible B=150k error once A is converged; the K and B errors are
roughly additive. Same structure as B=300k: frozen at K=1, converged at
K=3, wall non-monotone (K=3 cheaper than K=2). BEST CHEAP CELL at L=41:
B=150k/A=3B = 52s / 487MB / med ~2-3e-3 total.
Figure: msd_L41_Kscan_B150k.png.

## Repo cleanup (2026-07-07)

1. admit_basis (displacement) added to pc_step_orbit_rep + binding + python
   wrapper - the winning scheme is now available in momentum space (the k=1
   comparison should be re-run cap-primary). Verified: L=5 momentum sanity
   reproduces 1.86e-3 exactly; displacement smoke holds cap under admit=2B.
2. PPVM_EXPM_STREAM streaming paths REMOVED (MfOp, per_col_norms,
   StreamOrbitOp, per_col_orbit_stream; ~200 lines): measured useless in
   both spaces (real: 0-13% RAM for 4-5x wall; momentum: -5% for 24x).
3. rk4_step / rk4_step_arr marked deprecated (drop-only truncation, weakest
   scheme; kept for --mode rk4).
4. xy-experiments: README rewritten to the current script set (7 stale
   entries removed), trotter_ladder.py tracked, lockfile committed.
pc_step_complex intentionally KEPT (not deprecated) - see handoff.
Real-space regression after streaming removal: M100k/1e-5 T=10 cell
reproduces 8.67e-3/4.09e-2 bit-identically.

## Full K scan to K=10, both B (2026-07-07)

  B=150k: K=1.25..10: 3.9e-2, 2.7e-2, 1.2e-2, 1.7e-3, 4.0e-3, 4.7e-3,
          1.7e-3, 4.0e-3 (vs 300k/K10 anchor); plateaus at the B-bias floor
          (~2.3e-3) from K=3 on, fluctuating within the ~2x noise band.
          Walls 22->203s, RSS 324MB->1.24GB (linear in A).
  B=300k: 1.8e-2, 9.4e-3, 2.9e-3, 1.3e-3, 5.9e-4, 1.7e-3, 2.8e-4, anchor;
          post-K=3 values are sub-noise scatter (note K=5 bump in BOTH
          families - cancellation noise). Walls 50->705s, RSS 519MB->2.4GB.
VERDICT: K=3 is simultaneously the convergence onset and the cost optimum;
beyond it cost grows linearly in A for noise-level changes. Improve via B,
not A. Figure msd_L41_Kscan_full.png.

## pc_step_complex removed (2026-07-07, user request)

Production fn + inner + expm_step_complex + PyO3 binding + python wrapper
deleted. Its remaining role (the full-space complex reference bridge in the
orbit-rep equivalence tests) is preserved by a test-local helper
`pc_step_complex_full` (untruncated two-hop enrichment + exact in-basis
complex exponential via expm_apply_mf_cxvec). The two sector-check feature
tests were removed with the feature; the k=0 real/complex equivalence test
and the orbit-vs-full projection test both pass. Building blocks kept:
leakage_complex, compute_action_sum_complex, expm_apply_mf_cxvec.
L=5 momentum exact-ED sanity still 1.86e-3 bit-identical.

## rk4_step removed (2026-07-07)

Fully removed (was only doc-deprecated): rk4_step + rk4_step_inner,
PyO3 binding, rk4_step_arr wrapper, and the --mode rk4 harness branch in
xy-experiments/main_realspace_ladder.py. compute_action_sum (the L*-apply
primitive it used) is retained. Tests 7/7; harness smoke reproduces the
recorded dt=0.05/M100k cell bit-identically. The crate now exposes exactly
two step functions: pc_step (real) and pc_step_orbit_rep (momentum), with
the unified (max_basis, admit_basis, drop_tol) truncation API.

## K=2 walls re-measured on idle machine (2026-07-07)

  B=150k/K=2: 37.8s/390MB (was 100s - contended on 07-06)
  B=300k/K=2: 131.2s/650MB (was 209s - contended)
Wall is now monotone-ish linear in A across both families; the earlier
"K=3 faster than K=2" was contention. Cost statement: K=3 ~ same wall as
K=2, converged -> default A=3B unchanged. All other same-day walls stand.

## tau_add promoted to a first-class argument; PPVM_K_LEAKAGE removed (2026-07-07)

Both step functions now take `tau_add: Option<f64>` directly (the natural,
dt/drop-independent parameterization established by the cliff study);
k_leakage() and the env var are gone. Binding + python wrappers updated;
harnesses (main_realspace_ladder --tau_add, k_pec_run/main_k_pec_*
tau_add + admit_basis args); scan_realspace_msd converts K tokens to
--tau_add = K*drop/dt at launch, preserving old cell semantics exactly.
pc_step_timed runs without the filter (documented). REGRESSION: the
historical K=1 cell (0.1/1e-3, L=7 T=2) reproduces median 5.19e-3 /
peak 77,820 bit-identically through the new path; L=5 momentum sanity
1.86e-3 unchanged; tests 7/7.
Final public truncation API, both spaces:
  (max_basis, admit_basis=None, drop_tol, tau_add=None) + protected.
Zero env-var numerics knobs remain.

## B=600k baseline at L=41 (2026-07-07)

  600k/K=2: 245s/1.16GB  med 1.57e-3 vs 600k/K3  MSD(5)=35.230
  600k/K=3: 462s/1.55GB  (new anchor)            MSD(5)=35.445
Rebased ladder (vs 600k/K3, t in [0.2,5]):
  150k/K3: 1.13e-3 (52s/487MB)   MSD(5)=35.380
  300k/K3: 5.88e-4 (128s/806MB)  MSD(5)=35.420
  300k/K7: 1.13e-3 (494s/1.9GB)  - within its B-bias band, K>3 buys nothing
FINDINGS:
1. B-convergence is cleanly first-order: bias halves per B-doubling
   (1.13e-3 -> 5.9e-4, ratio 1.9). Extrapolated residual of the 600k anchor
   ~3e-4; Richardson limit MSD(5) ~ 35.47(2).
2. 600k/K=2 (1.57e-3) is WORSE than 300k/K=3 (5.9e-4) at 2x the cost:
   under-admission (K=2) costs more than B-doubling gains. A=3B confirmed
   at all three B.
3. vs Begusic-Chan 2^-18: anchor median 6.3e-3 - unchanged within
   digitization scatter; that comparison is digitization-limited ~5-6e-3.

## Fixed-A (fixed RAM) K/B split test, A=1.2M, L=41 (2026-07-07)

  B=300k K=4: 1041MB  median 4.75e-4
  B=400k K=3: 1164MB  median 2.69e-3     <- OUTLIER (should be best, isn't)
  B=600k K=2: 1162MB  median 1.57e-3
HONEST RESULT: NON-MONOTONE in both B and K -> we are in the cancellation-
NOISE FLOOR at this accuracy (~5e-4..3e-3, vs an anchor with ~3e-4 residual).
Single-cell median-rel differences below ~3x are NOT meaningful here. So the
fixed-A test does NOT resolve the (B,K) split, and my earlier "increase B at
K=3 traces the Pareto frontier / K=4 on the frontier" claims OVER-READ this
noise. What IS robust across all data:
  - K<3 measurably worse (under-admission): K=1 frozen, K=1.25-1.5 ~1-4e-2,
    K=2 ~3e-3-1.6e-2 (B-dependent). Need K>=3.
  - B-ladder AT K=3: 150k 1.1e-3 -> 300k 5.9e-4 -> 600k anchor: monotone,
    ~first-order, the one clean accuracy trend.
  - K>=3 at fixed B, or (B,K) split at fixed A: within noise, unresolved.
RAM ARGUMENT (theory, not from noisy cells): RAM ~ A = K*B. Beyond the K~3
knee, extra admission is provably neutral in the large-A limit (admitted-
then-truncated strings can't change the kept set), while B keeps improving.
So spend RAM on B, keep K just at the knee (~3). The cells can't PROVE the
fixed-RAM optimum at this accuracy, but the K-neutrality argument + the clean
B-ladder both point the same way. To resolve empirically would need seed/
observable averaging to beat down the ~2-3x cancellation noise.

## Branch triage (2026-07-07)

The two dangling worktree-agent branches were inspected before deletion —
both turned out to be FULLY TRIAGED already, with proper merge-then-revert
bookkeeping in the campaign and load-robust verdicts in their ledgers:
- parallel map_insert (3d5f0409): cherry-picked as 15f96554, final verdict
  DISCARD (1.1x wall for +50% RAM on a quiet machine; see
  2026-07-02-trotter-ladder for the triple-measurement saga + process
  lesson), reverted as 9ee8b3ba.
- leakage2 action-cache reuse (191e92aa): cherry-picked as 2411240f, verdict
  DISCARD (2.3x RSS - load-independent - for a wall gain bounded by the 17%
  leakage share; see 2026-07-02-expm-ladder), reverted as 42a5fa9c. Holds
  a fortiori under displacement (A=2-3B working sets).
Both branch refs deleted; content preserved in campaign history + ledgers.
Remaining branch state: autotune/ladder-tuning = the campaign superset
(contains expm-pc-step, expm-memory, symmetry-merging); local main is 50
behind origin/main -> upstream sync is the outstanding git task before PR.

## Branch triage part 2: perf/mimalloc-allocator + lindblad-shim (2026-07-07)

perf/mimalloc-allocator (June 17, 1 commit): PR #129 CLOSED unmerged; the
campaign independently carries the measured version (53f7ff29: mimalloc +
chunked leakage -> ~50% peak RSS cut) and current code has mimalloc in
Cargo.toml. Superseded -> local ref DELETED (remote cleanup optional).

lindblad-shim (June 17, 5 unique commits vs campaign): PR #98 ("Adaptive
Pauli-Lindbladian shim + adaptive-evolution demo") is still OPEN; the
campaign branched off mid-review and is now 162 commits ahead - the de
facto successor. Its unique commits: (a) b07f75de MAINTAINER-REQUESTED API
changes: PcStepConfig struct bundling the step knobs (replacing the long
positional list + clippy allow) and an error.rs extraction - the campaign
never adopted this and has since grown MORE positional args (max_basis,
admit_basis, drop_tol, tau_add, num_threads). The same review feedback
will return on the campaign PR; (b) f14edafc drops mimalloc - CONFLICTS
with the campaign's measured keep (campaign wins, it has data);
(c) expm-engine removals done differently/superseded by mf_expm evolution.
DISPOSITION: keep the branch while PR #98 is open. TODO before the campaign
PR: port the PcStepConfig + error.rs pattern onto the current API (even
more justified now with 5 knobs), then close #98 as superseded.

## Why momentum is SLOW: orbit-rep per-term canonicalization tax (2026-07-07)

User asked why momentum-space MSD is so slow when the orbit basis is |G|x
SMALLER than real space. Microbenchmark (matched basis B=20-30k, one pc_step):
  real-space L=11:        5.6 us/term  (flat in L)
  momentum L=6  |G|=12:  30.2 us/term
  momentum L=11 |G|=22:  66.4 us/term  (13.5x real space)
  momentum L=16 |G|=32: 142.5 us/term
Per-term cost grows ~|G|*N: the culprit is group.canonicalize_with_shift(q)
in build_orbit_rep_cols (orbit_rep.rs:95) and leakage_orbit_rep (:381) -
called for EVERY action term of EVERY rep, every build (2x/step), scanning
all |G|=2L translations (each an O(N) word compare) to map each generated
term back to its canonical rep. Real space has no canonicalization.

IMPLICATION (the real answer): for a FULL MSD you need all |G| momentum
modes. Converging one k-mode needs ~B_real/|G| terms (the compression), so
  total momentum work = |G| modes x (B_real/|G|) terms x c_orbit
                      = B_real x c_orbit  ~=  13.5x  x (real-space work).
=> reconstructing a real-space quantity (MSD) via all momentum modes is
NET SLOWER than doing it in real space, by the per-term tax (13.5x at L=11,
growing with L). Momentum space is the WRONG tool for full MSD.
Momentum WINS for: (a) MEMORY - |G|x fewer terms/mode, modes done serially;
(b) SINGLE k-resolved quantities (one transport channel, one RP resonance,
D(k) at one k) - there you evolve B_real/|G| terms at 13.5x/term = net
~1.6x FASTER than a full real-space run AND |G|x less RAM. That is the
method's actual niche, not real-space-observable reconstruction.

OPTIMIZATION OPPORTUNITY: the canonicalization map (rep->canonical-rep,shift)
is FIXED across steps and across the 2 expm calls for the stable part of the
basis; it is recomputed every build. Caching it (keyed by rep Word) would
drop the orbit per-term cost toward the real-space ~6 us, reviving the
compression as a genuine WALL advantage for single-mode runs. Also
canonicalize_with_shift itself is O(|G|*N); minimal-rotation tricks could
make it ~O(N). Filed for a future optimization session.
CONSEQUENCE FOR THE MOMENTUM MSD COMPARISON: not pursued to completion -
full-MSD-in-momentum is the wrong benchmark. The right momentum benchmark
is single-k D(k) vs Trotter, where the compression pays off.

## D extraction protocol + dt=0.2 B-convergence (2026-07-08)

Begusic-Chan D extraction (from their process.ipynb, confirmed):
  D = linregress(t[-51:], MSD[-51:]).slope/2, i.e. a single linear fit over
  t_theirs in [10,20] = t_ours [2.5,5], slope/2; then delta/dt->0 extrapolation
  across thresholds (Fig 2d). IDENTICAL to our scalar-D protocol; our recomputed
  per-threshold D match theirs exactly (2^-17/-18/-19 -> 3.287/3.515/3.661 ours).
  Our windowed D(t) (fig06) is the LOCAL-slope version of the same estimator.

dt=0.2 B-convergence (L=41, A=3B unless noted; [2.5,5] regression D):
  dt=0.1 K=3: D flat ~3.79 across B=50k-600k (converged). MSD(5)=35.25.
  dt=0.2 K=3: D = 3.43/3.74/3.58/3.75 (B=50/100/150/300k) - SCATTERED, LOW,
              not converged. MSD(5)~34.2.
  dt=0.2 K=5: D = 3.81 (B=150k), 3.82 (B=300k) - RECOVERS to ~3.79-3.81.
  dt=0.2 K=10 B=150k: 3.79.
FINDING: the admission knee K scales with dt. At larger dt each exp(dt L*)
step spreads the operator more, so A=3B headroom is insufficient (admission-
limited) - need K>=5 at dt=0.2 vs K=3 at dt=0.1. Once K is large enough,
dt=0.2 converges to the same D=3.79-3.81. NB dt=0.1 remains the sweet spot
for MSD(5) absolute value (35.25 vs 34.6 at dt=0.2/K5); the SLOPE (D) converges
at both once K is adequate. Recipe refinement: K_knee grows with dt; use
K~3 at dt<=0.1, K~5 at dt=0.2.

## D vs B for different dt (L=41, A=3B unless noted) (2026-07-08)

  dt=0.05 K3: D = 3.63/3.65/3.73/3.78 (B=50/100/150/300k) - climbs to 3.79,
              needs larger B (more truncation events at small dt).
  dt=0.1  K3: D ~ 3.79 flat for all B (sweet spot; B=100k=3.91 is noise).
  dt=0.2  K3: 3.43/3.74/3.58/3.75 - admission-limited, scattered.
  dt=0.2  K5: 3.81/3.82 - recovers.
COMPLETE PICTURE: D converges to ~3.79 for all dt once the ADEQUATE knob is
used, but the binding knob differs: small dt is B-limited (truncation-event
accumulation -> need larger B), large dt is A/K-limited (per-step spreading
-> need larger K). dt=0.1/K=3 is the joint sweet spot (converged at the
smallest B and K). Figure: /tmp/D_vs_B_dt.png (diagnostic; not yet a repo fig).

## "Just push to larger B?" - answered at B=600k (2026-07-08)

At B=600k, K=3:  dt=0.05 D=3.787, dt=0.1 D=3.796, dt=0.2 D=3.684.
=> Larger B CONVERGES dt=0.05 and dt=0.1 to ~3.79 (dt=0.05 was B-limited,
   climbs 3.63->3.79 over B=50k->600k). But dt=0.2 STAYS at 3.68 even at
   B=600k - larger B does NOT fix it because it is ADMISSION/K-limited, not
   B-limited. dt=0.2 needs K=5 (->3.81), independent of B.
Clean rule: small dt is B-limited (push B); large dt is K-limited (push K).
dt<=0.1 with K=3 + adequate B is the efficient converged regime.

## B=1M points added (2026-07-08)

B=1M, K=3:  dt=0.1 D=3.783 (MSD5 35.37);  dt=0.05 D=3.752 (MSD5 35.18).
vs B=600k:  dt=0.1 3.796;  dt=0.05 3.787.
=> No systematic climb from 600k->1M; both sit in 3.75-3.80, i.e. the
plateau is confirmed (flat within the ~+-0.03 extraction/cancellation-noise
floor). The converged D at dt<=0.1 is 3.76-3.79, matching Begusic-Chan 3.76.
Extraction noise (not B) is now the limiting uncertainty at large B.

## The B=100k D-outlier explained (2026-07-08)

The B=100k/dt=0.1/K=3 point gives D=3.91 (vs ~3.79 neighbors). CAUSE: its
MSD OVERSHOOTS the converged curve across [2.5,5] by a growing margin
(+0.25 at t=2.5 -> +0.84 at t=5; MSD(5)=36.2 vs 35.4 converged). Growing
offset => steeper slope => inflated D. NOT under-resolution (that gives LOW
D, cf. 50k/75k). This is DETERMINISTIC non-monotone truncation bias in B
(the method has no randomness; the top-B kept set changes discontinuously
with B and the coherent retained weight can over/undershoot) - the same
non-monotonicity Begusic-Chan report for SPD vs delta. CORRECTION: earlier
"cancellation-noise floor" mislabels this as statistical; it is deterministic
non-monotone convergence. Implication unchanged: use the median/band and the
D(t) collapse, never a single D(B) point (B=100k would mislead by ~3%).

## B=3M, K=3 (overnight, 2026-07-08) - all dt converge to D~3.79

  dt=0.2 : D[2.5,5]=3.791  MSD(5)=34.63  (58 min, 8.1 GB, 3.0M strings)
  dt=0.1 : D[2.5,5]=3.774  MSD(5)=35.34  (71 min, 9.0 GB)
  dt=0.05: D[2.5,5]=3.794  MSD(5)=35.43  (137 min, 8.8 GB)
All three land at D=3.77-3.79 = converged, matching Begusic-Chan 3.76 (~1%).

REVISION of the earlier dt=0.2 conclusion: I claimed "dt=0.2 is K/admission-
limited; larger B does NOT fix it, needs K=5." The B=3M/K=3 point (3.79)
shows larger B DOES converge dt=0.2 at K=3 - just slowly/noisily (the
B<=600k/K=3 series scattered 3.43-3.75; it was converging in B all along,
not stuck). Correct statement: dt=0.2/K=3 converges in B but needs MUCH
larger B (~3M) than dt=0.1 (~300k); K=5 is the cheaper route (converged by
B=150k). So "small dt B-limited / large dt K-limited" is too binary: large
dt is BOTH slower-in-B AND helped by K; either knob reaches 3.79.
dt=0.05: B=1M gave 3.752, B=3M gives 3.794 - confirms it was mildly
B-limited and is now converged (climbed up as expected). Plateau confirmed
at 3x the previous largest basis.

## B=6M, K=3 (overnight, 2026-07-09) - plateau confirmed at 6x

  dt=0.1 : D[2.5,5]=3.765  MSD(5)=35.32  (171 min, 15.6 GB, 6M strings)
  dt=0.05: D[2.5,5]=3.762  MSD(5)=35.19  (306 min, 15.6 GB)
B-ladder (D[2.5,5]): dt=0.1: 1M 3.783 / 3M 3.774 / 6M 3.765; dt=0.05: 1M
3.752 / 3M 3.794 / 6M 3.762. All within +-0.03 noise, centered ~3.76 =
Begusic-Chan value exactly. Plateau confirmed at 6x the earlier largest B
(2x the B=3M overnight run). Converged D = 3.76-3.77. Peak RSS 15.6 GB
(A=18M working set); this is near the practical single-node ceiling on 34 GB.
