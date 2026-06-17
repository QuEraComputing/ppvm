#!/usr/bin/env julia
# Goal (1): PauliPropagation.jl benchmark — MSD of the 2-leg XY ladder under
# real-space 2nd-order (Strang) Trotter, NO symmetry merging. Mirrors the ppvm
# run in main_realspace_ladder.py (same bonds, Strang structure, localized Z
# seed, |coeff| truncation) so the MSD curves overlay and the runtime/memory
# can be compared head-to-head.
#
# Usage:  julia -t 10 xy_ladder_msd.jl [L] [dt] [nsteps] [min_abs_coeff] [max_weight] [out.csv]
#   defaults: L=41 dt=0.02 nsteps=100 min_abs_coeff=1e-6 max_weight=0(=off) out=data/g1_pp.csv
#
# Output: CSV with columns  t,msd  (one row per step).  Wall time is reported
# by the NERSC `/usr/bin/time -v` wrapper; we also print it.
#
# Verified locally (PauliPropagation.jl, Julia 1.12) against ppvm at L=5,
# dt=0.02: MSD(t) curves agree to ~1e-9. If the NERSC PauliPropagation version
# differs and the API names (PauliSum / PauliString / PauliRotation /
# propagate / getcoeff) or the `theta` angle convention change, re-run the
# small-L cross-check against main_realspace_ladder.py and adjust.
using PauliPropagation
using Printf

arg(i, d) = length(ARGS) >= i ? ARGS[i] : d
L         = parse(Int,     arg(1, "41"))
dt        = parse(Float64, arg(2, "0.02"))
nsteps    = parse(Int,     arg(3, "100"))
mincoeff  = parse(Float64, arg(4, "1e-6"))
maxweight = parse(Int,     arg(5, "0"))          # 0 => no weight cap
outfile   = arg(6, "data/g1_pp.csv")

nlegs = 2
nq = nlegs * L
j0 = div(L, 2)
site(j, a) = j + a * L + 1                        # j∈0:L-1, a∈0:1  -> 1-indexed qubit
chain_coord(q) = (q - 1) % L                      # chain position of qubit q

# Brick-wall bond layering (matches ppvm): vertex-disjoint layers so adjacent
# gates within a sweep commute. Degree-3 ladder -> up to 4 layers: even leg-
# bonds, odd leg-bonds, the wrap seam (odd ring isn't 2-edge-colourable), rungs.
legbond(j, a) = (site(j, a), site(mod(j + 1, L), a))
bonds = Tuple{Int,Int}[]
for a in 0:1, j in 0:2:L-2; push!(bonds, legbond(j, a)); end   # even leg bonds
for a in 0:1, j in 1:2:L-2; push!(bonds, legbond(j, a)); end   # odd leg bonds
for a in 0:1; push!(bonds, legbond(L - 1, a)); end             # wrap seam
for j in 0:L-1; push!(bonds, (site(j, 0), site(j, 1))); end    # rungs

# Second-order (O(dt^3)) per-bond Strang: forward sweep over the bonds then a
# reversed sweep, each bond rxx then ryy. rxx and ryy commute on a bond, so the
# reversed bond order makes the step a palindrome -> O(dt^3). Truncation MATCHES
# ppvm: rxx never drops (min_abs_coeff=0; it creates the hopping intermediate),
# ryy drops at |coeff|<D right after -- keeping the Z-hopping flux, hence total
# Z, well conserved under truncation.
theta = dt
function trotter_step!(ps, D)
    for (a, b) in vcat(bonds, reverse(bonds))
        ps = propagate([PauliRotation([:X, :X], [a, b])], ps, [theta]; min_abs_coeff = 0.0)
        ps = maxweight > 0 ?
            propagate([PauliRotation([:Y, :Y], [a, b])], ps, [theta]; min_abs_coeff = D, max_weight = maxweight) :
            propagate([PauliRotation([:Y, :Y], [a, b])], ps, [theta]; min_abs_coeff = D)
    end
    return ps
end

# localized seed observable O(0) = 0.5 Z_{(j0,0)} + 0.5 Z_{(j0,1)}  (Σ a_q = 1)
psum = PauliSum(nq)
add!(psum, PauliString(nq, :Z, site(j0, 0), 0.5))
add!(psum, PauliString(nq, :Z, site(j0, 1), 0.5))

# single-Z target strings for the profile readout
ztargets = [PauliString(nq, :Z, q, 1.0) for q in 1:nq]
# min-image chain displacement from j0, per qubit
dj = [let d = mod(chain_coord(q) - j0, L); d > div(L, 2) ? d - L : d end for q in 1:nq]

function msd(ps)
    num = 0.0; den = 0.0
    for q in 1:nq
        c = getcoeff(ps, ztargets[q])             # coefficient of Z_q in the sum
        num += dj[q]^2 * c
        den += c
    end
    return num / den
end

ts  = Float64[]
ms  = Float64[]
push!(ts, 0.0); push!(ms, msd(psum))
t0 = time()
for s in 1:nsteps
    global psum = trotter_step!(psum, mincoeff)
    push!(ts, s * dt); push!(ms, msd(psum))
end
wall = time() - t0

mkpath(dirname(outfile))
open(outfile, "w") do io
    println(io, "t,msd")
    for (t, m) in zip(ts, ms); @printf(io, "%.10g,%.10g\n", t, m); end
end
@printf("done -> %s  (L=%d, nq=%d, %d steps, wall %.1fs, MSD(T)=%.6g, basis=%d)\n",
        outfile, L, nq, nsteps, wall, ms[end], length(psum))
