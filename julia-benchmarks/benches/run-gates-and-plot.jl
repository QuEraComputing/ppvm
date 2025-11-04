include("gates.jl")

using BenchmarkPlots, StatsPlots
plot(results)
savefig("julia-gates-benchmarks.svg")
