# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Dissipative generators on the orbit-rep path.

`pc_step_orbit_rep` evolves canonical translation-orbit representatives
with complex coefficients under the convention set by the momentum
projector `canonicalize_basis_arr_complex` (1/|G| prefactor):

    c_rep = coeff_word / |Stab(word)|.

Under this convention the phase-aware action is exact for ANY equivariant
generator, including jump and Kossakowski dissipators: transitions between
orbits of different sizes (e.g. the non-unital Z → I flow of σ⁻ decay,
where I is stabilized by the whole group) carry the correct weight because
the stabilizer factors on both sides cancel against the projector
normalization. A corollary used below: every orbit contributes
``|G| · c_rep`` to a translation-invariant {I,Z} sum, independent of its
stabilizer, so the emission rate in rep space is ``R = |G| · Σ c_rep``
over {I,Z} reps.

Model: ring of N emitters (positions on a circle → circulant J and Γ,
exact C_N symmetry), collective σ⁻ decay, momentum sector k = 0.
"""

import numpy as np
import pytest

from ppvm import Lindbladian
from ppvm._core import TranslationGroup, canonicalize_basis_arr_complex
from ppvm.lindblad import _basis_to_codes, _codes_to_basis, sigma_minus

from ._dissipative_refs import (
    couplings,
    eigenmode_jumps,
    exact_rate_sector,
    hamiltonian_terms,
    rate_observable,
    ring_positions,
)

BIG = 10_000_000


def ring_model(n, d_over_lam=0.1):
    J, Gam = couplings(ring_positions(n, d_over_lam))
    assert np.allclose(Gam, np.roll(np.roll(Gam, 1, 0), 1, 1)), "Gamma not circulant"
    assert np.linalg.eigvalsh(Gam).min() > -1e-12, "Gamma not PSD"
    ops = [sigma_minus(j, n) for j in range(n)]
    return J, Gam, hamiltonian_terms(n, J), ops, rate_observable(n, Gam)


def to_rep(basis, coeff, group, mom):
    b, c = canonicalize_basis_arr_complex(basis, np.asarray(coeff, dtype=np.complex128), group, mom)
    return dict(zip(_codes_to_basis(b), c))


@pytest.mark.parametrize("representation", ["kossakowski", "eigenmode"])
def test_orbit_matches_full_basis(representation):
    """N=6 ring, uncapped: orbit-rep evolution equals the full-basis
    real-space evolution projected to rep space, coefficient by
    coefficient."""
    n, dt, steps = 6, 0.02, 4
    _, Gam, h_terms, ops, obs = ring_model(n)
    if representation == "kossakowski":
        lind = Lindbladian(n, h_terms, kossakowski=(ops, Gam))
    else:
        lind = Lindbladian(n, h_terms, eigenmode_jumps(ops, Gam))
    group = TranslationGroup.chain_1d(n)
    mom = np.array([0], dtype=np.int32)

    strings = list(obs)
    basis0 = _basis_to_codes(strings, n)
    coeff0 = np.array([obs[s] for s in strings])

    b, c = basis0.copy(), coeff0.copy()
    for _ in range(steps):
        b, c = lind.pc_step_arr(b, c, dt, max_basis=BIG, drop_tol=0.0)
    full = to_rep(b, c, group, mom)

    br, cr = canonicalize_basis_arr_complex(basis0, coeff0.astype(np.complex128), group, mom)
    for _ in range(steps):
        br, cr = lind.pc_step_orbit_rep(
            br, cr, dt, max_basis=BIG, group=group, momentum=mom, drop_tol=0.0
        )
    orbit = dict(zip(_codes_to_basis(br), cr))

    assert set(full) == set(orbit)
    max_dev = max(abs(full[s] - orbit[s]) for s in full)
    assert max_dev < 1e-12, f"orbit vs full-basis: max |dc| = {max_dev:.2e}"


def orbit_rate_trace(lind, obs, n, dt, steps, group, mom, max_basis=BIG, admit=None):
    """R(t) from the orbit-rep evolution: |G| * sum of {I,Z}-rep coeffs."""
    strings = list(obs)
    basis0 = _basis_to_codes(strings, n)
    coeff0 = np.array([obs[s] for s in strings])
    br, cr = canonicalize_basis_arr_complex(basis0, coeff0.astype(np.complex128), group, mom)
    R = np.zeros(steps + 1)
    peak = 0
    for k in range(steps + 1):
        iz = np.all((br == 0) | (br == 2), axis=1)
        R[k] = n * cr[iz].sum().real
        peak = max(peak, len(cr))
        if k == steps:
            break
        br, cr = lind.pc_step_orbit_rep(
            br,
            cr,
            dt,
            max_basis=max_basis,
            group=group,
            momentum=mom,
            drop_tol=0.0,
            admit_basis=admit,
        )
    return R, peak


def test_orbit_rate_vs_exact_cascade():
    """N=6 ring, T=1, full rep basis: R(t) traced down from the orbit-rep
    evolution matches the excitation-cascade ED to < 1e-4."""
    n, dt, T = 6, 0.01, 1.0
    steps = round(T / dt)
    J, Gam, h_terms, ops, obs = ring_model(n)
    lind = Lindbladian(n, h_terms, kossakowski=(ops, Gam))
    group = TranslationGroup.chain_1d(n)
    mom = np.array([0], dtype=np.int32)

    R_orbit, _ = orbit_rate_trace(lind, obs, n, dt, steps, group, mom)
    _, R_exact = exact_rate_sector(n, J, Gam, T_run=T, dt_out=dt)
    err = np.abs(R_orbit - R_exact).max()
    assert err < 1e-4, f"orbit-rep R(t) vs exact cascade: max |dR| = {err:.2e}"


def test_orbit_truncated_sanity():
    """N=10 ring, genuinely truncated: no NaNs or blowup, and the error
    against the (matched-capacity) real-space run improves monotonically
    with the rep budget."""
    n, dt, steps = 10, 0.02, 20
    _, Gam, h_terms, ops, obs = ring_model(n)
    lind = Lindbladian(n, h_terms, kossakowski=(ops, Gam))
    group = TranslationGroup.chain_1d(n)
    mom = np.array([0], dtype=np.int32)

    # Real-space reference at matched effective capacity B_full = n * B_reps.
    b = _basis_to_codes(list(obs), n)
    c = np.array([obs[s] for s in obs])
    R_full = np.zeros(steps + 1)
    for k in range(steps + 1):
        iz = np.all((b == 0) | (b == 2), axis=1)
        R_full[k] = c[iz].sum()
        if k == steps:
            break
        b, c = lind.pc_step_arr(
            b, c, dt, max_basis=n * 2048, drop_tol=0.0, admit_basis=3 * n * 2048
        )

    errs = {}
    for b_reps in (512, 2048):
        R, peak = orbit_rate_trace(
            lind,
            obs,
            n,
            dt,
            steps,
            group,
            mom,
            max_basis=b_reps,
            admit=3 * b_reps,
        )
        assert np.all(np.isfinite(R)), f"non-finite R(t) at B_reps={b_reps}"
        assert np.abs(R).max() < 5 * n, f"R(t) blowup at B_reps={b_reps}"
        assert peak <= 3 * b_reps, "admission bound violated"
        errs[b_reps] = np.abs(R - R_full).max()

    assert errs[2048] <= errs[512], (
        f"no improvement with rep budget: err(2048)={errs[2048]:.3e} > err(512)={errs[512]:.3e}"
    )
    # At matched capacity the two representations should agree closely.
    assert errs[2048] < 0.05 * np.abs(R_full).max(), (
        f"orbit-rep tracks real-space poorly: {errs[2048]:.3e}"
    )


def test_identity_bookkeeping_closed_form():
    """Uniform single-site decay (K = Γ0·1): from O = Σ_j Z_j the exact
    solution is coeff_{Z_j}(t) = e^{-Γ0 t} per site and
    coeff_I(t) = -n(1 - e^{-Γ0 t}) (each site pours into the identity).
    In rep space c_Z is the member coefficient (trivial stabilizer) and
    c_I = coeff_I / |G| = -(1 - e^{-Γ0 t}) (the identity is stabilized by
    the whole group) — this pins the non-unital bookkeeping."""
    n, dt, steps = 6, 0.01, 40
    ops = [sigma_minus(j, n) for j in range(n)]
    lind = Lindbladian(n, [], kossakowski=(ops, np.eye(n)))
    group = TranslationGroup.chain_1d(n)
    mom = np.array([0], dtype=np.int32)

    # O = Σ_j Z_j → one Z rep with c = 1 (k = 0 eigenstate);
    # canonicalize_first rewrites the row to the lex-min representative.
    strings = ["Z" + "I" * (n - 1)]
    br = _basis_to_codes(strings, n)
    cr = np.array([1.0 + 0.0j])
    for _ in range(steps):
        br, cr = lind.pc_step_orbit_rep(
            br,
            cr,
            dt,
            max_basis=BIG,
            group=group,
            momentum=mom,
            drop_tol=0.0,
            canonicalize_first=True,
        )
    out = dict(zip(_codes_to_basis(br), cr))

    t = steps * dt
    (z_rep,) = [s for s in out if s.count("Z") == 1 and set(s) <= {"I", "Z"}]
    c_z = out[z_rep]
    c_i = out["I" * n]
    assert abs(c_z - np.exp(-t)) < 1e-9, f"c_Z = {c_z} vs {np.exp(-t)}"
    assert abs(c_i - (-(1 - np.exp(-t)))) < 1e-9, (
        f"c_I = {c_i} vs coeff_I/|G| = {-(1 - np.exp(-t))}"
    )
