# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""End-to-end wall-time scaling of the pure-Rust predictor-corrector step
with rayon thread count.

The entire ``pc_step`` body (both leakage calls, the generator build, and
both ``expm_multiply`` calls) is wrapped in a rayon pool of the requested
size — leakage and generator parallelise over basis elements, the matrix
exponential parallelises over SpMV rows. So the speedup numbers reflect
overall PC throughput, not just SpMV.

Usage::

    python lindblad_pc_scaling.py --L 51 --steps 20 --max-cores 16

For HPC: ``--max-cores`` is the upper bound of the thread sweep, so on a
32-core node you'd pass ``--max-cores 32`` (or any subset).
"""

from __future__ import annotations

import argparse
import importlib.util as _ilu
import sys as _sys
import time
from pathlib import Path as _Path
from statistics import median

import numpy as np

# Bypass the `ppvm` package __init__ to avoid a pre-existing stim/bloqade
# version mismatch in the dev env. Pull the bare submodule path directly.
_spec = _ilu.spec_from_file_location(
    "_ppvm_lindblad",
    _Path(__file__).resolve().parent.parent / "src" / "ppvm" / "lindblad.py",
)
_mod = _ilu.module_from_spec(_spec)
_sys.modules["_ppvm_lindblad"] = _mod
_spec.loader.exec_module(_mod)
Lindbladian = _mod.Lindbladian


def build_nn_xy_dephasing(L: int, J: float, gamma: float):
    h_terms = []
    for i in range(L - 1):
        a, b = i, i + 1
        xs = ["I"] * L
        xs[a] = xs[b] = "X"
        ys = ["I"] * L
        ys[a] = ys[b] = "Y"
        h_terms += [("".join(xs), J), ("".join(ys), J)]
    jump_terms = [("I" * j + "Z" + "I" * (L - j - 1), gamma) for j in range(L)]
    return h_terms, jump_terms


def build_long_range_xy_dephasing(L: int, J: float, alpha: float, gamma: float):
    """All-to-all XY with 1/r^α couplings (matches the lindblad_adaptive
    demo). Operator weight grows much faster than NN because every bond
    activates every step, giving a basis size that meaningfully exercises
    parallel scaling."""
    pairs = [
        (a, b, 1.0 / min(b - a, L - b + a) ** alpha)
        for a in range(L)
        for b in range(a + 1, L)
    ]
    kac = sum(j for _, _, j in pairs) / L
    pairs = [(a, b, j / kac) for a, b, j in pairs]
    h_terms = []
    for a, b, j in pairs:
        for q in "XY":
            term = ["I"] * L
            term[a] = term[b] = q
            h_terms.append(("".join(term), J * j))
    jump_terms = [("I" * j + "Z" + "I" * (L - j - 1), gamma) for j in range(L)]
    return h_terms, jump_terms


def run_pc_steps(L_op, L, site0, dt, n_steps, tau_add, num_threads):
    z_strings = ["I" * j + "Z" + "I" * (L - j - 1) for j in range(L)]
    basis = [z_strings[site0]]
    coeffs = np.array([1.0])
    protected = [z_strings[site0]]
    times = []
    for _ in range(n_steps):
        t0 = time.perf_counter()
        basis, coeffs = L_op.pc_step(
            basis,
            coeffs,
            dt,
            tau_add,
            protected=protected,
            parallel_threshold=0,
            num_threads=num_threads,
        )
        times.append(time.perf_counter() - t0)
    return times, len(basis)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--L", type=int, default=51, help="number of qubits")
    parser.add_argument("--steps", type=int, default=20, help="PC steps per measurement")
    parser.add_argument("--dt", type=float, default=0.05, help="time step")
    parser.add_argument("--J", type=float, default=1.0, help="XY exchange")
    parser.add_argument("--gamma", type=float, default=1.0, help="Z dephasing rate")
    parser.add_argument(
        "--model",
        choices=("nn", "long-range"),
        default="long-range",
        help="NN XY or all-to-all 1/r^alpha XY (default: long-range)",
    )
    parser.add_argument(
        "--alpha", type=float, default=1.0, help="power-law exponent for long-range model"
    )
    parser.add_argument(
        "--tau", type=float, default=1e-8, help="leakage threshold for basis expansion"
    )
    parser.add_argument(
        "--max-cores", type=int, default=4, help="upper bound of the thread sweep"
    )
    parser.add_argument(
        "--warmup-steps",
        type=int,
        default=4,
        help="PC steps to discard before timing (lets the basis settle)",
    )
    args = parser.parse_args()

    L = args.L
    site0 = L // 2
    if args.model == "nn":
        h_terms, jump_terms = build_nn_xy_dephasing(L, args.J, args.gamma)
    else:
        h_terms, jump_terms = build_long_range_xy_dephasing(
            L, args.J, args.alpha, args.gamma
        )
    L_op = Lindbladian(L, h_terms, jump_terms)

    # Warmup: pre-build pools at each thread count and amortise one-time setup.
    print(f"Warmup (L={L}, {args.warmup_steps} steps per pool size)...")
    for n in range(1, args.max_cores + 1):
        run_pc_steps(L_op, L, site0, args.dt, args.warmup_steps, args.tau, n)
    L_op.clear_cache()

    print()
    model_str = (
        "NN XY (OBC)" if args.model == "nn" else f"long-range XY alpha={args.alpha}"
    )
    print(
        f"{model_str} + Z-dephasing  L={L}  J={args.J}  gamma={args.gamma}  "
        f"dt={args.dt}  tau={args.tau:g}"
    )
    print(
        f"{args.steps} PC steps per measurement, parallel_threshold=0 "
        f"(force parallel SpMV); rayon pool wraps the entire PC step "
        f"(leakage + generator + expm)."
    )
    print()
    header = (
        f"{'threads':>8s}  {'first-step (ms)':>16s}  {'steady (ms)':>12s}  "
        f"{'speedup':>9s}  {'|basis|':>8s}"
    )
    print(header)
    print("-" * len(header))

    baseline = None
    for n in range(1, args.max_cores + 1):
        L_op.clear_cache()
        times, basis_size = run_pc_steps(
            L_op, L, site0, args.dt, args.steps, args.tau, n
        )
        first = times[0] * 1000.0
        steady = median(times[1:]) * 1000.0
        if baseline is None:
            baseline = steady
        speedup = baseline / steady
        print(
            f"{n:>8d}  {first:>16.1f}  {steady:>12.2f}  {speedup:>8.2f}x  "
            f"{basis_size:>8d}"
        )


if __name__ == "__main__":
    main()
