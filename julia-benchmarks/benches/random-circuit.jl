using BenchmarkTools
using PauliPropagation


function layer!(circuit, n)
    for i in 1:n
        rz = PauliRotation([:Z], [i], 1.1)
        ry = PauliRotation([:Y], [i], 2.1)

        push!(circuit, rz)
        push!(circuit, ry)
        push!(circuit, rz)
    end
end

function entangle!(circuit, n)
    for i in 1:n
        push!(circuit,
            CliffordGate(:CNOT, [i, mod1(i + 1, n)])
        )
    end
end

function random_circuit(n, depth)
    circuit = Gate[]

    # NOTE: inverse order since it propagation is done with reverse(circuit)
    for _ in 1:depth
        entangle!(circuit, n)
        layer!(circuit, n)
    end

    return circuit
end


function run_circuit(circuit, state)
    propagate!(circuit, state)
    overlapwithzero(state)
end

function run_benchmark(n, depth)
    s = PauliSum(n)
    add!(s, PauliString(n, [:Z, :Z], [1, 2]))

    circuit = random_circuit(n, depth)

    return @benchmark run_circuit($circuit, state) setup = (state = copy($s))
end

result = run_benchmark(4, 2)
