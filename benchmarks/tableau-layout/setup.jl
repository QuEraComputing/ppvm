# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
# Use a small environment instead of instantiating QC's docs/test workspace.
using Pkg
using TOML
length(ARGS) == 2 || error("usage: setup.jl quantumclifford_checkout environment_directory")
checkout, environment = abspath.(ARGS)
# Replay the recorded dependency versions while relocating the two path packages.
recorded = joinpath(@__DIR__, "results")
if isfile(joinpath(recorded, "Manifest.toml")) && !isfile(joinpath(environment, "Manifest.toml"))
    mkpath(environment)
    cp(joinpath(recorded, "Project.toml"), joinpath(environment, "Project.toml"); force=true)
    manifest = TOML.parsefile(joinpath(recorded, "Manifest.toml"))
    manifest["deps"]["QuantumClifford"][1]["path"] = checkout
    manifest["deps"]["QECCore"][1]["path"] = joinpath(checkout, "lib", "QECCore")
    open(joinpath(environment, "Manifest.toml"), "w") do output
        TOML.print(output, manifest)
    end
end
Pkg.activate(environment)
Pkg.develop([PackageSpec(path=checkout), PackageSpec(path=joinpath(checkout, "lib", "QECCore"))]; preserve=Pkg.PRESERVE_ALL)
Pkg.instantiate()
