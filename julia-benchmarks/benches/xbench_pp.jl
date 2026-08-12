# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
#
# The PauliPropagation.jl side of the cross-library Pauli-propagation
# benchmark. See `benchmarks/cross-library/README.md` for the shared circuit
# definitions, the parameter contract, and the CSV schema — every runner reads
# the same environment variables and prints the same columns.
#
#   MODEL=tfim QUBITS=8,16,24 STEPS=10 DT=0.1 JCOUP=1.0 HFIELD=1.0 ATOL=1e-6 \
#     julia --project=@. -t1 benches/xbench_pp.jl

using PauliPropagation
using Printf

const MODEL = get(ENV, "MODEL", "tfim")
const STEPS = parse(Int, get(ENV, "STEPS", "10"))
const DT = parse(Float64, get(ENV, "DT", "0.1"))
const JCOUP = parse(Float64, get(ENV, "JCOUP", "1.0"))
const HFIELD = parse(Float64, get(ENV, "HFIELD", "1.0"))
const ATOL = parse(Float64, get(ENV, "ATOL", "1e-6"))
const ITERS = parse(Int, get(ENV, "ITERS", "3"))

const THETA_BOND = 2 * JCOUP * DT
const THETA_SITE = 2 * HFIELD * DT

"""The seed observable: `Σ_i Z_i` for TFIM, `Z_1` for Heisenberg."""
function seed_operator(n::Int)
    ps = PauliSum(n)
    if MODEL == "tfim"
        for i in 1:n
            add!(ps, PauliString(n, [:Z], [i]))
        end
    else
        add!(ps, PauliString(n, [:Z], [1]))
    end
    return ps
end

"""One first-order Trotter step, gates in the order the shared spec fixes."""
function trotter_step!(state, n::Int)
    if MODEL == "tfim"
        for i in 1:n
            state = propagate(PauliRotation([:X], [i], THETA_SITE), state; min_abs_coeff = ATOL)
        end
        for i in 1:(n - 1)
            state = propagate(
                PauliRotation([:Z, :Z], [i, i + 1], THETA_BOND), state; min_abs_coeff = ATOL
            )
        end
    else
        for i in 1:(n - 1)
            for axes in ([:X, :X], [:Y, :Y], [:Z, :Z])
                state = propagate(
                    PauliRotation(axes, [i, i + 1], THETA_BOND), state; min_abs_coeff = ATOL
                )
            end
        end
        for i in 1:n
            state = propagate(PauliRotation([:Z], [i], THETA_SITE), state; min_abs_coeff = ATOL)
        end
    end
    return state
end

function run_model(n::Int)
    state = seed_operator(n)
    for _ in 1:STEPS
        state = trotter_step!(state, n)
    end
    return state
end

"""`⟨0…0|O|0…0⟩` for TFIM; the `Z_1` autocorrelator for Heisenberg."""
function readout(state, n::Int)
    if MODEL == "tfim"
        return real(overlapwithzero(state))
    else
        c = getcoeff(state, PauliString(n, [:Z], [1]))
        return real(c)
    end
end

"""
Print the propagated support as `word coefficient`, largest first, site 1
leftmost — the same format the Rust and Python runners emit under `DUMP=1` so
the driver can diff all four term-for-term.
"""
function dump_support(n::Int)
    state = run_model(n)
    terms = Dict{String,Float64}()
    for (ps, c) in state
        w = join([['I', 'X', 'Y', 'Z'][getpauli(ps, i) + 1] for i in 1:n])
        terms[w] = get(terms, w, 0.0) + real(c)
    end
    println("# $(length(terms)) terms")
    for w in sort(collect(keys(terms)); by = w -> (-abs(terms[w]), w))
        println("$w $(Printf.@sprintf("%+.12e", terms[w]))")
    end
end

function main()
    qubits = parse.(Int, split(get(ENV, "QUBITS", "8,12,16,20,24,28,32"), ","))
    if haskey(ENV, "DUMP")
        dump_support(first(qubits))
        return
    end
    println(
        stderr,
        "PauliPropagation.jl $MODEL: steps=$STEPS dt=$DT J=$JCOUP h=$HFIELD atol=$ATOL " *
        "iters=$ITERS threads=$(Threads.nthreads())",
    )
    println("model,library,qubits,steps,dt,atol,time_s,terms,observable")
    for n in qubits
        # Warm up so the reported time excludes JIT for this width's type
        # specialization.
        st = run_model(n)
        nterms = length(st)
        obs = readout(st, n)
        best = Inf
        for _ in 1:ITERS
            best = min(best, @elapsed run_model(n))
        end
        println("$MODEL,pauli-propagation-jl,$n,$STEPS,$DT,$ATOL,$(round(best, sigdigits=7)),$nterms,$obs")
        println(stderr, "  n=$n  $(round(best, digits=4))s  $nterms terms  obs=$obs")
        flush(stdout)
    end
end

main()
