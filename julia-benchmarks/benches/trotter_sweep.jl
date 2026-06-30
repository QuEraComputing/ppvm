# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
#
# TFIM Trotter runtime vs qubit count in PauliPropagation.jl, mirroring
# crates/ppvm-runtime/examples/trotter_qubit_sweep.rs so the curves line up
# with the ppvm fxhash/gxhash sweep.
#
# Circuit: initial observable sum_i Z_i; per Trotter step, RX(theta_x) + a
# depolarizing channel on every qubit, then RZZ(theta_zz) + depolarizing on
# every bond. Truncation min_abs_coeff = 1e-6.
#
# Params come from the environment so they match the Rust run exactly:
#   J=1.0 STEPS=20 QUBITS="8,16,..." ITERS=2 \
#     julia --project=@. -t1 benches/trotter_sweep.jl > pp.csv

using PauliPropagation

const H = 1.0
const DT = 0.1 / H
const J = parse(Float64, get(ENV, "J", "1.0")) * H
const STEPS = parse(Int, get(ENV, "STEPS", "20"))
const MIN_ABS_COEFF = parse(Float64, get(ENV, "MIN_ABS_COEFF", "1e-6"))
const P_DEPOL = 1e-4

const THETA_X = DT * H
const THETA_ZZ = DT * J

function build_state(n::Int)
    ps = PauliSum(n)
    for i in 1:n
        add!(ps, PauliString(n, [:Z], [i]))
    end
    return ps
end

function run_trotter(n::Int)
    state = build_state(n)
    for _ in 1:STEPS
        for i in 1:n
            rx = PauliRotation([:X], [i], THETA_X)
            state = propagate(rx, state; min_abs_coeff = MIN_ABS_COEFF)
            state = propagate(DepolarizingNoise(i, P_DEPOL), state; min_abs_coeff = MIN_ABS_COEFF)
        end
        for i in 1:(n - 1)
            rzz = PauliRotation([:Z, :Z], [i, i + 1], THETA_ZZ)
            state = propagate(rzz, state; min_abs_coeff = MIN_ABS_COEFF)
            state = propagate(DepolarizingNoise(i, P_DEPOL), state; min_abs_coeff = MIN_ABS_COEFF)
            state = propagate(DepolarizingNoise(i + 1, P_DEPOL), state; min_abs_coeff = MIN_ABS_COEFF)
        end
    end
    return state
end

function main()
    qubits = parse.(Int, split(get(ENV, "QUBITS", "8,16,24,32,40,48,56,64,80,96,112,122"), ","))
    iters = parse(Int, get(ENV, "ITERS", "2"))

    println(stderr, "PP.jl params: steps=$STEPS theta_x=$THETA_X theta_zz=$THETA_ZZ min_abs_coeff=$MIN_ABS_COEFF threads=$(Threads.nthreads())")
    println("qubits,hasher,bytes,time_s,terms")
    for n in qubits
        # Warm up (JIT-compiles the term-type specialization for this tier).
        st = run_trotter(n)
        nterms = length(st.terms)
        best = Inf
        for _ in 1:iters
            t = @elapsed run_trotter(n)
            best = min(best, t)
        end
        println("$n,pauli_propagation_jl,0,$(round(best, sigdigits = 6)),$nterms")
        println(stderr, "n=$n  PP.jl  $(round(best, digits = 4))s  ($nterms terms)")
        flush(stdout)
    end
end

main()
