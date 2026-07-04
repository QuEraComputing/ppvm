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
