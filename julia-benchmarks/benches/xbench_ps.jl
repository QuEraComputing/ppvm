# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
#
# The PauliStrings.jl side of the cross-library Pauli-propagation benchmark.
# See `benchmarks/cross-library/README.md` for the shared circuit definitions,
# the parameter contract, and the CSV schema.
#
# PauliStrings.jl's own front door for this is `evolve(H, O, tspan;
# method=Trotter())`, which builds the gate list from the Hamiltonian's internal
# string order. We build the `TrotterGate` list by hand instead, so the gate
# sequence is the one the shared spec fixes and the propagated operator is
# comparable term-for-term with the other three engines.
#
#   MODEL=tfim QUBITS=8,16,24 STEPS=10 DT=0.1 JCOUP=1.0 HFIELD=1.0 ATOL=1e-6 \
#     julia --project=@. -t1 benches/xbench_ps.jl

using PauliStrings
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

"""The `PauliString` generator for a one- or two-site Pauli term."""
function generator(n::Int, ops::Vector{String}, sites::Vector{Int})
    o = Operator(n)
    if length(sites) == 1
        o += ops[1], sites[1]
    else
        o += ops[1], sites[1], ops[2], sites[2]
    end
    return o.strings[1]
end

"""
The Trotter gate list in *application* order, per the shared spec.

`trotter_step!` consumes gates in matrix-multiply order and applies them in
reverse, so the caller reverses this before handing it over.
"""
function gate_list(n::Int)
    gates = PauliStrings.TrotterGate[]
    if MODEL == "tfim"
        for i in 1:n
            push!(gates, PauliStrings.TrotterGate(generator(n, ["X"], [i]), THETA_SITE))
        end
        for i in 1:(n - 1)
            push!(
                gates,
                PauliStrings.TrotterGate(generator(n, ["Z", "Z"], [i, i + 1]), THETA_BOND),
            )
        end
    else
        for i in 1:(n - 1)
            for op in ("X", "Y", "Z")
                push!(
                    gates,
                    PauliStrings.TrotterGate(
                        generator(n, [op, op], [i, i + 1]), THETA_BOND
                    ),
                )
            end
        end
        for i in 1:n
            push!(gates, PauliStrings.TrotterGate(generator(n, ["Z"], [i]), THETA_SITE))
        end
    end
    return gates
end

"""The seed observable: `Σ_i Z_i` for TFIM, `Z_1` for Heisenberg."""
function seed_operator(n::Int)
    o = Operator(n)
    if MODEL == "tfim"
        for i in 1:n
            o += "Z", i
        end
    else
        o += "Z", 1
    end
    return o
end

function run_model(n::Int, gates)
    o = seed_operator(n)
    truncation(x) = cutoff(x, ATOL)
    for _ in 1:STEPS
        PauliStrings.trotter_step!(o, gates; truncation = truncation, truncate_every = 1)
    end
    return o
end

"""`⟨0…0|O|0…0⟩` for TFIM; the `Z_1` autocorrelator for Heisenberg."""
function readout(o::Operator, n::Int)
    if MODEL == "tfim"
        return real(expect(o, "0"^n))
    else
        z1 = seed_operator(n)
        return real(trace_product(z1, o; scale = 1))
    end
end

"""
Print the propagated support as `word coefficient`, largest first, site 1
leftmost — the same format the other three runners emit under `DUMP=1`.

PauliStrings.jl carries `im^{#Y}` inside the stored coefficient (its `Matrix`
convention), so `op_to_strings` is the readout that puts the coefficients on the
same footing as the other engines' real ones.
"""
function dump_support(n::Int, gates)
    o = run_model(n, gates)
    coeffs, strings = op_to_strings(o)
    terms = Dict{String,Float64}()
    for (c, s) in zip(coeffs, strings)
        w = replace(s, '1' => 'I')
        terms[w] = get(terms, w, 0.0) + real(c)
    end
    println("# $(length(terms)) terms")
    for w in sort(collect(keys(terms)); by = w -> (-abs(terms[w]), w))
        @printf("%s %+.12e\n", w, terms[w])
    end
end

function main()
    qubits = parse.(Int, split(get(ENV, "QUBITS", "8,12,16,20,24,28,32"), ","))
    if haskey(ENV, "DUMP")
        n = first(qubits)
        dump_support(n, reverse(gate_list(n)))
        return
    end
    println(
        stderr,
        "PauliStrings.jl $MODEL: steps=$STEPS dt=$DT J=$JCOUP h=$HFIELD atol=$ATOL " *
        "iters=$ITERS threads=$(Threads.nthreads())",
    )
    println("model,library,qubits,steps,dt,atol,time_s,terms,observable")
    for n in qubits
        # `gate_list` is circuit construction, not propagation — built once and
        # excluded from the timed region, as in every other runner.
        gates = reverse(gate_list(n))
        # Warm up so the reported time excludes JIT.
        o = run_model(n, gates)
        nterms = length(o)
        obs = readout(o, n)
        best = Inf
        for _ in 1:ITERS
            best = min(best, @elapsed run_model(n, gates))
        end
        println("$MODEL,pauli-strings-jl,$n,$STEPS,$DT,$ATOL,$(round(best, sigdigits=7)),$nterms,$obs")
        println(stderr, "  n=$n  $(round(best, digits=4))s  $nterms terms  obs=$obs")
        flush(stdout)
    end
end

main()
