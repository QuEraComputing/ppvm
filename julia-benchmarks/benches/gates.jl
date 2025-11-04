using BenchmarkTools
using PauliPropagation


function build_initial_state(n_qubits)
    initial_state = PauliSum(n_qubits)
    add!(initial_state, PauliString(n_qubits, [:Z, :Z], [1, 2]))

    for _ in 1:2
        for i in 1:n_qubits
            rz = PauliRotation([:Z], [i], 1.1)
            ry = PauliRotation([:Y], [i], 2.1)
            initial_state = propagate(rz, initial_state; min_abs_coeff=0)
            initial_state = propagate(ry, initial_state; min_abs_coeff=0)
            initial_state = propagate(rz, initial_state; min_abs_coeff=0)
        end

        for i in 1:n_qubits
            j = i + 1 - (i == n_qubits) * n_qubits
            cnot = CliffordGate(:CNOT, [i, j])
            initial_state = propagate(cnot, initial_state; min_abs_coeff=0)
        end
    end

    return initial_state
end

function benchmark_suite()
    # parameters
    n_qubits = 12

    # setup initial state
    initial_state = build_initial_state(n_qubits)

    # temporary state used in PP.jl
    tmp_state = similar(initial_state)
    
    println("Using $(Threads.nthreads()) threads")
    println("Initial state has $(length(initial_state.terms)) terms")

    group = BenchmarkGroup()

    # define the set of gates we want to benchmark
    x = CliffordGate(:X, [1])
    y = CliffordGate(:Y, [1])
    z = CliffordGate(:Z, [1])
    h = CliffordGate(:H, [1])
    
    cnot = CliffordGate(:CNOT, [1, 2])
    cz = CliffordGate(:CZ, [1, 2])

    rx = PauliRotation([:X], [1], 0.5)
    ry = PauliRotation([:Y], [1], 0.5)
    rz = PauliRotation([:Z], [1], 0.5)

    rxx = PauliRotation([:X, :X], [1, 2], 0.5)
    ryy = PauliRotation([:Y, :Y], [1, 2], 0.5)
    rzz = PauliRotation([:Z, :Z], [1, 2], 0.5)

    # collect the gate applications into the benchmark group
    # group["x"] = @benchmarkable applytoall!($x, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    # group["y"] = @benchmarkable applytoall!($y, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    # group["z"] = @benchmarkable applytoall!($z, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    # group["h"] = @benchmarkable applytoall!($h, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    
    # group["cnot"] = @benchmarkable applytoall!($cnot, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    # group["cz"] = @benchmarkable applytoall!($cz, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))

    # group["rx"] = @benchmarkable applytoall!($rx, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    # group["ry"] = @benchmarkable applytoall!($ry, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    # group["rz"] = @benchmarkable applytoall!($rz, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))

    # group["rxx"] = @benchmarkable applytoall!($rxx, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    # group["ryy"] = @benchmarkable applytoall!($ryy, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    # group["rzz"] = @benchmarkable applytoall!($rzz, nothing, state, cache_state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))


    group["x"] = @benchmarkable propagate!($x, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    group["y"] = @benchmarkable propagate!($y, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    group["z"] = @benchmarkable propagate!($z, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    group["h"] = @benchmarkable propagate!($h, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    
    group["cnot"] = @benchmarkable propagate!($cnot, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    group["cz"] = @benchmarkable propagate!($cz, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))

    group["rx"] = @benchmarkable propagate!($rx, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    group["ry"] = @benchmarkable propagate!($ry, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    group["rz"] = @benchmarkable propagate!($rz, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))

    group["rxx"] = @benchmarkable propagate!($rxx, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    group["ryy"] = @benchmarkable propagate!($ryy, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))
    group["rzz"] = @benchmarkable propagate!($rzz, state) setup = (state = copy($initial_state); cache_state = copy($tmp_state))

    return group
end


function run_benchmark()
    group = benchmark_suite()
    tune!(group)
    return run(group)
end

results = run_benchmark()

