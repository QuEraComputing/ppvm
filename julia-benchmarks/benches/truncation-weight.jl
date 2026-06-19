# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
#
# Mirror of crates/ppvm-runtime/benches/truncation-weight.rs.
# Measures PauliPropagation.jl's truncate! over the same parameter grid:
#   - 3 term-weight profiles (3, 50, 120)
#   - 4 weight cutoffs (10, 100, 1000, Inf)
#   - 1000 terms on 128 qubits
#
# PP.jl's truncate! always builds a single closure combining min_abs_coeff
# and max_weight checks, then walks the dictionary once. To approximate the
# Rust "single-strategy" benches we set the other threshold to a sentinel
# value that never trips.

using PauliPropagation
using BenchmarkTools

const N_QUBITS = 128
const N_TERMS = 1_000
const COEFF_EPS = 1e-12

const PROFILES = [
    ("w3", 3),
    ("w50", 50),
    ("w120", 120),
]

const CUTOFFS = [
    ("10", 10.0),
    ("100", 100.0),
    ("1000", 1000.0),
    ("inf", Inf),
]

const PAULI_SYMS = [:X, :Y, :Z]

# Build N_TERMS PauliStrings each with exactly target_weight non-identity
# slots (clamped to N_QUBITS), spread by stride. Pauli choices depend on k
# so terms are mostly distinct. Matches make_terms() in the Rust mirror.
function make_psum(target_weight)
    weight = min(target_weight, N_QUBITS)
    stride = max(div(N_QUBITS, max(weight, 1)), 1)
    psum = PauliSum(N_QUBITS)
    for k in 0:(N_TERMS - 1)
        positions = Int[]
        symbols = Symbol[]
        for w in 0:(weight - 1)
            push!(positions, mod(k + w * stride, N_QUBITS) + 1)  # 1-indexed
            push!(symbols, PAULI_SYMS[mod(k * 31 + w, 3) + 1])
        end
        coeff = 1.0 / (k + 1)
        add!(psum, PauliString(N_QUBITS, symbols, positions, coeff))
    end
    return psum
end

function bench_max_weight_alone()
    println("\n=== max-weight-only (min_abs_coeff = 0, never trips) ===")
    results = Dict{String, Any}()
    for (profile, weight) in PROFILES
        psum = make_psum(weight)
        println("[$profile] psum has $(length(psum)) terms")
        for (cutoff_name, cutoff) in CUTOFFS
            id = "$profile/cut-$cutoff_name"
            t = @benchmark PauliPropagation.truncate!(state;
                    min_abs_coeff = 0.0, max_weight = $cutoff) setup = (state = deepcopy($psum)) samples=200
            results[id] = t
            println("  $id: median = $(median(t).time) ns")
        end
    end
    return results
end

function bench_coeff_threshold_alone()
    println("\n=== coeff-threshold-only (max_weight = Inf, never trips) ===")
    results = Dict{String, Any}()
    for (profile, weight) in PROFILES
        psum = make_psum(weight)
        id = profile
        t = @benchmark PauliPropagation.truncate!(state;
                min_abs_coeff = $COEFF_EPS, max_weight = Inf) setup = (state = deepcopy($psum)) samples=200
        results[id] = t
        println("  $id: median = $(median(t).time) ns")
    end
    return results
end

function bench_combined()
    println("\n=== combined (both checks active, matches Rust CombinedStrategy) ===")
    results = Dict{String, Any}()
    for (profile, weight) in PROFILES
        psum = make_psum(weight)
        for (cutoff_name, cutoff) in CUTOFFS
            id = "$profile/cut-$cutoff_name"
            t = @benchmark PauliPropagation.truncate!(state;
                    min_abs_coeff = $COEFF_EPS, max_weight = $cutoff) setup = (state = deepcopy($psum)) samples=200
            results[id] = t
            println("  $id: median = $(median(t).time) ns")
        end
    end
    return results
end

function bench_clone_baseline()
    println("\n=== clone (deepcopy) baseline ===")
    results = Dict{String, Any}()
    for (profile, weight) in PROFILES
        psum = make_psum(weight)
        t = @benchmark deepcopy($psum) samples=200
        results[profile] = t
        println("  $profile: median = $(median(t).time) ns")
    end
    return results
end

println("Using $(Threads.nthreads()) threads")
println("N_QUBITS = $N_QUBITS, N_TERMS = $N_TERMS")

bench_max_weight_alone()
bench_coeff_threshold_alone()
bench_combined()
bench_clone_baseline()
