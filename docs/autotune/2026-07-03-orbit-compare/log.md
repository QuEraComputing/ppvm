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
1. Add orbit-rep STREAMING (matrix-free) to orbit_rep.rs (PPVM_EXPM_STREAM does
   NOT reach it -- it has its own build_orbit_rep_cols/OrbitRepCscOp). Needed for
   the expm stream-vs-cache comparison in the orbit-preserving case.
2. Scale to ~1M orbit-rep basis: large-N chain (orbit-rep basis ~ full/|G|), long
   enough T. No exact ED at that scale -> converged reference (tightest drop/dt).
3. Scan dt x drop_tol for BOTH methods; match rel~1e-3; compare wall + PEAK RSS.
   For expm: stream vs cache (leakage), both with a max_basis cut.

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
