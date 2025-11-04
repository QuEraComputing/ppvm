# Julia benchmarks for comparison

Start Julia in the local project with

```bash
# cd julia-benchmarks
julia --project=@. -t 4
```

Then, run benchmarks with e.g.

```julia
julia> include("benches/trotter.jl")
```
