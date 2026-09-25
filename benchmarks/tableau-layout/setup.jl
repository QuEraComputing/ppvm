# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
# Use a small environment instead of instantiating QC's docs/test workspace.
using Pkg
length(ARGS) == 2 || error("usage: setup.jl quantumclifford_checkout environment_directory")
checkout, environment = abspath.(ARGS)
Pkg.activate(environment)
Pkg.develop([PackageSpec(path=checkout), PackageSpec(path=joinpath(checkout, "lib", "QECCore"))])
Pkg.instantiate()
