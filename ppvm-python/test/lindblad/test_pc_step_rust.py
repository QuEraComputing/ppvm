# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Pure-Rust :meth:`Lindbladian.pc_step` (Al-Mohy & Higham expm + parallel
SpMV): agrees with the numpy-eigendecomp PC reference at FP precision and
shows the same cubic dt-scaling against the bilinear reference. Plus a
sanity check that a length-1 real lincomb routes to the Hermitian fast path.
"""

from __future__ import annotations

from itertools import pairwise

import numpy as np

from ppvm import Lindbladian

from ._helpers import (
    adaptive_z_correlator_pc,
    bilinear_nn_xy_z_dephasing_obc,
    nn_xy_z_dephasing_obc,
    random_pauli_str,
    xy_dephasing,
)


def _adaptive_z_correlator_pc_rust(L_op, L, site0, dt, n_steps, tau_add):
    """Same as :func:`_helpers.adaptive_z_correlator_pc` but the per-step PC
    work (leakage expansion, predictor expm, second-hop expansion, corrector
    expm) all runs in Rust through :meth:`Lindbladian.pc_step`."""
    z_strings = ["I" * j + "Z" + "I" * (L - j - 1) for j in range(L)]
    basis = [z_strings[site0]]
    coeffs = np.array([1.0])
    protected = [z_strings[site0]]

    corr = np.zeros((n_steps + 1, L))
    corr[0, site0] = 1.0

    for step in range(n_steps):
        basis, coeffs = L_op.pc_step(basis, coeffs, dt, tau_add, protected=protected)
        index = {s: i for i, s in enumerate(basis)}
        for j in range(L):
            if z_strings[j] in index:
                corr[step + 1, j] = coeffs[index[z_strings[j]]]
    return corr


def test_pc_step_rust_matches_python_pc():
    """The pure-Rust PC step agrees with the numpy-eigendecomp PC reference
    at FP precision.

    Pins the Rust matrix exponential (Al-Mohy & Higham) against an
    independent reference (numpy ``eig``) under the exact same
    basis-expansion schedule."""
    L = 4
    J = 1.0
    gamma = 1.0
    site0 = L // 2
    dt = 0.01
    n_steps = 5
    tau_add = 1e-12

    h_terms, jump_terms = nn_xy_z_dephasing_obc(L, J, gamma)
    L_op = Lindbladian(L, h_terms, jump_terms)

    rust = _adaptive_z_correlator_pc_rust(L_op, L, site0, dt, n_steps, tau_add)
    py_ref = adaptive_z_correlator_pc(L_op, L, site0, dt, n_steps, tau_add)

    diff = float(np.max(np.abs(rust - py_ref)))
    assert diff < 1e-10, f"Rust PC differs from numpy-eigendecomp PC by {diff:.3e}"


def test_pc_step_rust_dt_scaling_is_cubic():
    """End-to-end: the Rust-only PC step matches the bilinear reference with
    cubic dt-scaling, confirming the Rust matrix exponential is not the
    accuracy bottleneck."""
    L = 4
    J = 1.0
    gamma = 1.0
    site0 = L // 2
    T = 0.05
    tau_add = 1e-12

    h_terms, jump_terms = nn_xy_z_dephasing_obc(L, J, gamma)
    L_op = Lindbladian(L, h_terms, jump_terms)

    err = []
    for dt in (0.01, 0.005, 0.0025):
        n_steps = round(T / dt)
        times = np.arange(n_steps + 1) * dt
        exact = bilinear_nn_xy_z_dephasing_obc(L, J, gamma, times, site0)
        rust = _adaptive_z_correlator_pc_rust(L_op, L, site0, dt, n_steps, tau_add)
        err.append(float(np.max(np.abs(rust[-1] - exact[-1]))))

    # Halving dt should drop the error by ≥5× (cubic gives 8×).
    for prev, curr in pairwise(err):
        assert curr < prev / 5, f"Rust PC dt-halving ratio < 5: errors {err}"
    assert err[-1] < 1e-7, f"Rust PC tight-dt error too large: {err[-1]:.3e}"


def test_lincomb_single_term_matches_hermitian_fast_path():
    """A length-1 real lincomb should route to the Hermitian fast path and
    produce numerically identical results to passing the string directly."""
    L = 4
    h_terms, jump_simple = xy_dephasing(L, alpha=1.0, gamma=0.3)
    L_simple = Lindbladian(L, h_terms, jump_simple)
    # Same operator, expressed as a length-1 complex lincomb.
    jump_lincomb = [([(s, 1.0 + 0.0j)], g) for s, g in jump_simple]
    L_lincomb = Lindbladian(L, h_terms, jump_lincomb)

    rng = np.random.default_rng(99)
    for _ in range(5):
        p = random_pauli_str(rng, L)
        got_a = L_simple.action(p)
        got_b = L_lincomb.action(p)
        for k in set(got_a) | set(got_b):
            assert got_a.get(k, 0.0) == got_b.get(k, 0.0), (
                f"lincomb fast path mismatch at p={p!r}: {got_a} vs {got_b}"
            )
