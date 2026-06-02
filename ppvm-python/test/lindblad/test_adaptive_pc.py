# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""End-to-end convergence of the adaptive predictor-corrector evolution
(numpy eigendecomposition reference) against the closed bilinear NN-XY +
Z-dephasing solution.

For NN interactions the JW bilinears stay closed under the adjoint
Lindbladian, so the spin correlator obeys a tractable L²×L² ODE
(:func:`bilinear_nn_xy_z_dephasing_obc`). OBC keeps this exact; PBC would
introduce a parity-twist that only matches up to 1/L corrections.
"""

from __future__ import annotations

from itertools import pairwise

import numpy as np

from ppvm import Lindbladian

from ._helpers import (
    adaptive_z_correlator,
    adaptive_z_correlator_pc,
    bilinear_nn_xy_z_dephasing_obc,
    nn_xy_z_dephasing_obc,
)


def test_adaptive_converges_to_nn_xy_z_dephasing_bilinear():
    """Halving dt drives the adaptive shim toward the closed bilinear solution.

    Single-hop has local truncation O(dt²) per step → global error O(T·dt).
    """
    L = 4
    J = 1.0
    gamma = 1.0
    site0 = L // 2
    T = 0.05
    tau_add = 1e-12  # tight enough that integrator (T·dt) dominates

    h_terms, jump_terms = nn_xy_z_dephasing_obc(L, J, gamma)
    L_op = Lindbladian(L, h_terms, jump_terms)

    errors = []
    final_corr = None
    for dt in (0.01, 0.005, 0.0025):
        n_steps = round(T / dt)
        times = np.arange(n_steps + 1) * dt
        exact = bilinear_nn_xy_z_dephasing_obc(L, J, gamma, times, site0)
        L_op.clear_cache()
        shim = adaptive_z_correlator(L_op, L, site0, dt, n_steps, tau_add)
        # Endpoint comparison only — independent of which step counts we use.
        errors.append(float(np.max(np.abs(shim[-1] - exact[-1]))))
        if dt == 0.0025:
            final_corr = (shim[-1], exact[-1])

    # Halving dt should roughly halve the error; allow 2× slack.
    assert errors[1] < 0.8 * errors[0], f"dt-halving did not help: {errors}"
    assert errors[2] < 0.8 * errors[1], f"dt-halving did not help: {errors}"

    # Integrator floor at dt = 0.0025 is ~T·dt = 1.25e-4; expect <1e-3.
    assert errors[-1] < 1e-3, (
        f"shim vs bilinear at smallest dt: max abs error = {errors[-1]:.3g}; "
        f"shim={final_corr[0]}, exact={final_corr[1]}"
    )


def test_predictor_corrector_lifts_dt_scaling_to_cubic():
    """The predictor-corrector basis expansion lifts the single-hop scheme's
    local O(dt²) truncation to O(dt³). PC error is also strictly smaller at
    every dt we test.
    """
    L = 4
    J = 1.0
    gamma = 1.0
    site0 = L // 2
    T = 0.05
    tau_add = 1e-12

    h_terms, jump_terms = nn_xy_z_dephasing_obc(L, J, gamma)
    L_op = Lindbladian(L, h_terms, jump_terms)

    err_single = []
    err_pc = []
    for dt in (0.01, 0.005, 0.0025):
        n_steps = round(T / dt)
        times = np.arange(n_steps + 1) * dt
        exact = bilinear_nn_xy_z_dephasing_obc(L, J, gamma, times, site0)
        L_op.clear_cache()
        single = adaptive_z_correlator(L_op, L, site0, dt, n_steps, tau_add)
        L_op.clear_cache()
        pc = adaptive_z_correlator_pc(L_op, L, site0, dt, n_steps, tau_add)
        err_single.append(float(np.max(np.abs(single[-1] - exact[-1]))))
        err_pc.append(float(np.max(np.abs(pc[-1] - exact[-1]))))

    # PC strictly more accurate than single-hop at every dt (by ~100× in this
    # regime). Threshold loose enough to absorb expm_multiply tolerance noise.
    for s, p, dt in zip(err_single, err_pc, (0.01, 0.005, 0.0025)):
        assert p < s / 50, (
            f"PC ({p:.3e}) not meaningfully better than single-hop ({s:.3e}) at dt={dt}"
        )

    # dt-scaling one order steeper: halving dt should drop the error ~8×
    # (dt³ vs single-hop's ~4×). Require >5× per halving with safety margin.
    for prev, curr in pairwise(err_pc):
        assert curr < prev / 5, f"PC dt-halving ratio < 5: errors {err_pc}"

    # Smallest-dt PC error should sit at FP noise of the bilinear reference.
    assert err_pc[-1] < 1e-7, f"PC at smallest dt: error = {err_pc[-1]:.3e}"
