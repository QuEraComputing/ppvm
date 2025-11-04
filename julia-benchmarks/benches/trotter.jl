using PauliPropagation
using BenchmarkTools


function trotter_circuit(n, total_time, dt, interaction_strength, external_field, noise_param)

    circuit = Gate[]

    steps = Int(total_time / dt)
    theta_zz = dt * interaction_strength
    theta_x = dt * external_field

    
    for _ in 1:steps

        for i in 1:n
            rx = PauliRotation([:X], [i], theta_x)
            push!(circuit, rx)
            push!(circuit, DepolarizingNoise(i, noise_param))
        end

        for i in 1:n-1
            rzz = PauliRotation([:Z, :Z], [i, i+1], theta_zz)
            push!(circuit, rzz)

            push!(circuit, DepolarizingNoise(i, noise_param))
            push!(circuit, DepolarizingNoise(i + 1, noise_param))
        end



    end

    return circuit
end

function run_benchmark()
    # set parameters here
    n_qubits = 12
    h = 1.0
    dt = 0.1 / h
    time = 1.0 / h
    j = 1.0 / 8.0 * h
    noise_param = 1e-4


    initial_state = PauliSum(n_qubits)

    for i = 1:n_qubits
        add!(initial_state, PauliString(n_qubits, [:Z], [i]))
    end

    println("Using $(Threads.nthreads()) threads")



    circuit = trotter_circuit(n_qubits, time, dt, j, h, noise_param)

    return @benchmark propagate!($circuit, state; min_abs_coeff = 1e-6) setup = (state = copy($initial_state))
end



bench = run_benchmark()