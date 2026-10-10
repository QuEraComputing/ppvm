# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
using PauliPropagation
using Printf
using Statistics

function circuit(n, depth, shared)
    gates = Gate[]
    ids = Int[]
    for _ in 1:depth
        for q in 1:n
            push!(gates, PauliRotation([:X], [q]))
            push!(ids, shared ? 1 : length(gates))
        end
        for q in 1:n-1
            push!(gates, PauliRotation([:Z,:Z], [q,q+1]))
            push!(ids, shared ? 2 : length(gates))
        end
    end
    gates, ids
end
angles(count, sample) = [0.07 + mod(i*17 + sample*11,53)/100 for i in 0:count-1]
function build(n, gates)
    seed = PauliString(n, [:Z], [1])
    surrogate = propagate(gates, wrapcoefficients(seed, NodePathProperties); max_weight=Inf,max_freq=Inf,max_sins=Inf)
    zerofilter!(surrogate)
end
function numeric(n,gates,values)
    propagate(gates,PauliString(n,[:Z],[1]),values; min_abs_coeff=0.0,max_weight=Inf,max_freq=Inf,max_sins=Inf)
end
word(w,n) = join([['I','X','Y','Z'][getpauli(w,q)+1] for q in 1:n])
function main()
    n=parse(Int,get(ENV,"QUBITS","3")); depth=parse(Int,get(ENV,"STEPS","3"))
    mode=get(ENV,"MODE","shared"); @assert mode in ("shared","independent")
    gates,ids=circuit(n,depth,mode=="shared")
    count=mode=="shared" ? 2 : length(gates)
    seed=PauliString(n,[:Z],[1])
    full=propagate(gates,wrapcoefficients(seed,NodePathProperties);max_weight=Inf,max_freq=Inf,max_sins=Inf)
    terms=length(full)
    max_error=0.0
    for sample in 0:4
        values=angles(count,sample)[ids]
        evaluate!(full,values)
        reference=numeric(n,gates,values)
        a=Dict(word(w,n)=>PauliPropagation.tonumber(c) for (w,c) in full)
        b=Dict(word(w,n)=>real(c) for (w,c) in reference)
        for w in union(keys(a),keys(b))
            max_error=max(max_error,abs(get(a,w,0.0)-get(b,w,0.0)))
        end
        if haskey(ENV,"VERIFY")
            for w in sort(collect(keys(a)))
                @printf("C,%d,%s,%.17e\n",sample,w,a[w])
            end
        end
    end
    @assert max_error<1e-10 "surrogate/numeric mismatch: $max_error"
    haskey(ENV,"VERIFY") && return
    trials=parse(Int,get(ENV,"TRIALS","7")); batch=parse(Int,get(ENV,"BATCH","100"))
    @assert trials>0 && batch>0
    assignments=[angles(count,sample)[ids] for sample in 0:batch-1]
    surrogate=build(n,gates)
    # Warm up construction, evaluation, and numeric propagation before timing.
    evaluate!(surrogate,assignments[1]); overlapwithzero(surrogate)
    build(n,gates); overlapwithzero(numeric(n,gates,assignments[1]));
    for trial in 0:trials-1
        GC.gc()
        start=time_ns(); fresh=build(n,gates); build_s=(time_ns()-start)/1e9
        start=time_ns()
        for values in assignments
            evaluate!(surrogate,values)
            overlapwithzero(surrogate)
        end
        eval_s=(time_ns()-start)/1e9/batch
        start=time_ns()
        for values in assignments
            overlapwithzero(numeric(n,gates,values))
        end
        numeric_s=(time_ns()-start)/1e9/batch
        evaluate!(surrogate,assignments[1])
        @printf("paulipropagation-jl,%d,%d,%s,%d,%.12e,%.12e,%.12e,%d,%.17e,%.3e\n",n,depth,mode,trial,build_s,eval_s,numeric_s,terms,overlapwithzero(surrogate),max_error)
    end
end
main()
