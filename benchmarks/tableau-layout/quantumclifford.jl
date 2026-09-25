# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
# Run through the tableau-layout Python harness, or pass fixture_path samples min_ms.
using QuantumClifford
using QuantumClifford: Tableau, mul_right!
using LinearAlgebra: Adjoint

const FNV_OFFSET = UInt64(14695981039346656037)
const FNV_PRIME = UInt64(1099511628211)
const TIMING_SINK = Ref(UInt64(0))
const TIMING_STATE = Ref{Any}(nothing)
const OPERATIONS = (:comm, :mul, :h, :s, :cnot, :circuit)

function read_fixture(path, ::Type{W}, layout) where {W<:Unsigned}
    lines = readlines(path)
    n = parse(Int, first(lines))
    @assert length(lines) == 2n + 1
    t = zero(Tableau{Vector{UInt8},Matrix{W}}, 2n, n)
    for (row, line) in enumerate(lines[2:end])
        phase, xs, zs = split(line)
        @assert length(xs) == length(zs) == n
        t.phases[row] = parse(UInt8, phase)
        for q in 1:n
            t[row, q] = (xs[q] == '1', zs[q] == '1')
        end
    end
    layout == "fastrow" ? fastrow(t) : fastcolumn(t)
end

# Base.copy of an Adjoint can materialize a Matrix. Copy its parent explicitly.
copy_layout(xzs::Adjoint) = copy(parent(xzs))'
copy_layout(xzs::Matrix) = copy(xzs)
clone(t::Tableau) = Tableau(copy(t.phases), t.nqubits, copy_layout(t.xzs))

Base.@noinline function sweep!(t, ::Val{:comm})
    value = UInt64(0)
    for row in eachindex(t)
        value += comm(t, row, mod1(row + 1, length(t)))
    end
    value
end

Base.@noinline function sweep!(t, ::Val{:mul})
    for row in eachindex(t)
        mul_right!(t, row, mod1(row + 1, length(t)))
    end
    UInt64(0)
end

Base.@noinline function sweep!(t, ::Val{OP}) where {OP}
    s = Stabilizer(t)
    n = nqubits(t)
    for q in 1:n
        OP in (:h, :circuit) && apply!(s, sHadamard(q))
        OP in (:s, :circuit) && apply!(s, sPhase(q))
        OP in (:cnot, :circuit) && apply!(s, sCNOT(q, mod1(q + 1, n)))
    end
    UInt64(0)
end

function run!(t, operation, iterations)
    value = UInt64(0)
    for _ in 1:iterations
        value += sweep!(t, operation)
    end
    value
end

function checksum(t, ::Val{:comm})
    h = FNV_OFFSET
    for row in eachindex(t)
        h = (h ⊻ UInt64(comm(t, row, mod1(row + 1, length(t))))) * FNV_PRIME
    end
    h
end

function checksum(t, operation)
    sweep!(t, operation)
    h = FNV_OFFSET
    for row in eachindex(t)
        h = (h ⊻ UInt64(t.phases[row])) * FNV_PRIME
        for q in 1:nqubits(t)
            x, z = t[row, q]
            h = (h ⊻ (UInt64(x) | (UInt64(z) << 1))) * FNV_PRIME
        end
    end
    h
end

function time_batch(t, operation, iterations)
    work = clone(t)
    start = time_ns()
    value = run!(work, operation, iterations)
    elapsed = time_ns() - start
    TIMING_SINK[] ⊻= value
    TIMING_STATE[] = work
    elapsed
end

function benchmark(t, layout, operation, samples, minimum_ns)
    # Compile the whole timed path before calibration. Fixture copies and checksums
    # are outside the timer; every sample starts from the same tableau.
    digest = checksum(clone(t), operation)
    time_batch(t, operation, 1)
    iterations = 1
    while time_batch(t, operation, iterations) < minimum_ns
        iterations *= 2
    end
    n = nqubits(t)
    op = typeof(operation).parameters[1]
    count = op in (:comm, :mul) ? 2n : op == :circuit ? 3n : n
    for sample in 1:samples
        elapsed = time_batch(t, operation, iterations)
        println(join(("quantumclifford", layout, 8sizeof(eltype(t.xzs)), n,
                      op, iterations, sample, elapsed, count,
                      string(digest; base=16, pad=16)), ','))
        flush(stdout)
    end
end

function main(args)
    length(args) == 3 || error("usage: quantumclifford.jl fixture_or_directory samples min_ms")
    fixture, sample_arg, min_arg = args
    samples = parse(Int, sample_arg)
    minimum_ns = parse(Float64, min_arg) * 1e6
    @assert samples > 0 && minimum_ns > 0
    paths = isdir(fixture) ? sort(filter(p -> endswith(p, ".txt"), readdir(fixture; join=true))) : [fixture]
    layouts = split(get(ENV, "QC_LAYOUTS", "fastrow,fastcolumn"), ',')
    bits = parse.(Int, split(get(ENV, "QC_WORD_BITS", "8,16,32,64,128"), ','))
    words = Dict(8=>UInt8, 16=>UInt16, 32=>UInt32, 64=>UInt64, 128=>UInt128)
    println("implementation,layout,word_bits,n,operation,iterations,sample,elapsed_ns,operations_per_iteration,checksum")
    for path in paths, width in bits, layout in layouts
        @assert layout in ("fastrow", "fastcolumn")
        t = read_fixture(path, words[width], layout)
        for op in OPERATIONS
            try
                benchmark(t, layout, Val(op), samples, minimum_ns)
            catch exception
                # SIMD.jl does not define UInt128 vector lanes. All other failures
                # are errors in the run and must fail the harness.
                width == 128 && layout == "fastrow" && op == :mul && exception isa TypeError || rethrow()
                println(stderr, "UNSUPPORTED quantumclifford layout=$layout word_bits=$width n=$(nqubits(t)) operation=$op: ",
                        sprint(showerror, exception))
                flush(stderr)
            end
        end
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
