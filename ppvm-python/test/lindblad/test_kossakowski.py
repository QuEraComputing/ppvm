# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Kossakowski-form dissipator: exact equivalence with the eigenmode-jump
representation, physics regression against the exact excitation-cascade
reference, and leakage coverage.
"""
import pathlib

import numpy as np
import pytest

from ppvm import Lindbladian
from ppvm.lindblad import _basis_to_codes, _codes_to_basis, sigma_minus

from ._dissipative_refs import (
    G0,
    chain_positions,
    couplings,
    eigenmode_jumps,
    exact_rate_sector,
    hamiltonian_terms,
    rate_observable,
)

BIG = 10_000_000  # uncapped max_basis

FIG11_H5 = pathlib.Path(
    "/Users/alexschuckert/dev/26_ppvm/CTPP Figures/fig11_superradiant_burst/data.h5"
)


def run_steps(lind, obs, n, dt, steps):
    """Evolve {string: coeff} a few uncapped pc steps; return final dict."""
    strings = list(obs)
    basis = _basis_to_codes(strings, n)
    coeff = np.array([obs[s] for s in strings], dtype=np.float64)
    for _ in range(steps):
        basis, coeff = lind.pc_step_arr(basis, coeff, dt, max_basis=BIG, drop_tol=0.0)
    return dict(zip(_codes_to_basis(basis), coeff))


def random_psd(rng, m, complex_k=False):
    a = rng.standard_normal((m, m))
    if complex_k:
        a = a + 1j * rng.standard_normal((m, m))
    return a @ a.conj().T / m


def random_lincomb(rng, n):
    """A random 2-term Pauli lincomb with complex coefficients."""
    ops = "IXYZ"
    out = []
    for _ in range(2):
        s = "".join(rng.choice(list(ops)) for _ in range(n))
        if s == "I" * n:
            s = "X" + s[1:]
        c = complex(rng.standard_normal(), rng.standard_normal())
        out.append((s, c))
    return out


@pytest.mark.parametrize("complex_k", [False, True])
def test_random_model_equivalence(complex_k):
    """Eigenmode jumps built from K and kossakowski=(ops, K) produce the
    same evolution to near machine precision (same dt, uncapped basis)."""
    rng = np.random.default_rng(7 if complex_k else 3)
    n, dt, steps = 4, 0.02, 3
    # ops: all single-site sigma^- plus one random 2-term lincomb
    ops = [sigma_minus(j, n) for j in range(n)] + [random_lincomb(rng, n)]
    K = random_psd(rng, len(ops), complex_k)
    h_terms = [("XX" + "I" * (n - 2), 0.9), ("I" + "ZZ" + "I" * (n - 3), -0.4)]
    obs = {"Z" + "I" * (n - 1): 1.0, "IXY" + "I" * (n - 3): 0.3}

    out_k = run_steps(Lindbladian(n, h_terms, kossakowski=(ops, K)), obs, n, dt, steps)
    out_j = run_steps(Lindbladian(n, h_terms, eigenmode_jumps(ops, K)), obs, n, dt, steps)

    assert set(out_k) == set(out_j)
    max_dev = max(abs(out_k[s] - out_j[s]) for s in out_k)
    assert max_dev < 1e-12, f"representations diverged: max |dc| = {max_dev:.2e}"


def test_kossakowski_coexists_with_jump_terms():
    """kossakowski= and jump_terms may both contribute."""
    n = 2
    ops = [sigma_minus(j, n) for j in range(n)]
    K = [[1.0, 0.5], [0.5, 1.0]]
    both = Lindbladian(n, [], [("ZI", 0.3)], kossakowski=(ops, K))
    only_k = Lindbladian(n, [], kossakowski=(ops, K))
    only_j = Lindbladian(n, [], [("ZI", 0.3)])
    a_both = both.action("XI")
    a_sum = {}
    for d in (only_k.action("XI"), only_j.action("XI")):
        for s, c in d.items():
            a_sum[s] = a_sum.get(s, 0.0) + c
    for s in set(a_both) | set(a_sum):
        assert abs(a_both.get(s, 0.0) - a_sum.get(s, 0.0)) < 1e-13


def superradiance_chain(n, d_over_lam=0.1):
    J, Gam = couplings(chain_positions(n, d_over_lam))
    ops = [sigma_minus(j, n) for j in range(n)]
    return J, Gam, hamiltonian_terms(n, J), ops, rate_observable(n, Gam)


def rate_trace(lind, obs, n, dt, steps):
    """R(t) on the fully inverted state = sum of {I,Z}-string coefficients."""
    strings = list(obs)
    basis = _basis_to_codes(strings, n)
    coeff = np.array([obs[s] for s in strings], dtype=np.float64)
    R = np.zeros(steps + 1)
    for k in range(steps + 1):
        iz = np.all((basis == 0) | (basis == 2), axis=1)  # codes: I=0, Z=2
        R[k] = coeff[iz].sum()
        if k == steps:
            break
        basis, coeff = lind.pc_step_arr(basis, coeff, dt, max_basis=BIG, drop_tol=0.0)
    return R


def test_superradiance_physics_regression():
    """N=6 subwavelength chain, full basis, T=1: the Kossakowski path matches
    the eigenmode path to ~1e-12 and the exact cascade reference to < 1e-4."""
    n, dt, T = 6, 0.01, 1.0
    steps = round(T / dt)
    J, Gam, h_terms, ops, obs = superradiance_chain(n)

    R_k = rate_trace(Lindbladian(n, h_terms, kossakowski=(ops, Gam)), obs, n, dt, steps)
    R_j = rate_trace(Lindbladian(n, h_terms, eigenmode_jumps(ops, Gam)), obs, n, dt, steps)
    assert np.abs(R_k - R_j).max() < 1e-11, (
        f"kossakowski vs eigenmode R(t): {np.abs(R_k - R_j).max():.2e}"
    )

    _, R_exact = exact_rate_sector(n, J, Gam, T_run=T, dt_out=dt)
    err = np.abs(R_k - R_exact).max()
    assert err < 1e-4, f"kossakowski vs exact cascade: max |dR| = {err:.2e}"


@pytest.mark.skipif(not FIG11_H5.exists(), reason="fig11 data.h5 not present")
def test_superradiance_vs_fig11_reference():
    """Cross-check R(t) against the stored exact reference of the
    superradiant-burst study (same model, N=6, first 1/Gamma_0)."""
    h5py = pytest.importorskip("h5py")
    n, dt, T = 6, 0.01, 1.0
    steps = round(T / dt)
    _, Gam, h_terms, ops, obs = superradiance_chain(n)
    R_k = rate_trace(Lindbladian(n, h_terms, kossakowski=(ops, Gam)), obs, n, dt, steps)
    with h5py.File(FIG11_H5, "r") as h5:
        R_ref = h5["n6/exact"][: steps + 1]
    err_ref = np.abs(R_k - R_ref).max()
    assert err_ref < 1e-4, f"vs fig11 data.h5 n6/exact: {err_ref:.2e}"


def test_leakage_covers_kossakowski_terms():
    """Leakage of a Z string under a pure-Kossakowski dissipator is nonzero
    and identical to the eigenmode-jump leakage (admission sees the same
    physics)."""
    n = 3
    ops = [sigma_minus(j, n) for j in range(n)]
    _, Gam = couplings(chain_positions(n, 0.1))
    lk = Lindbladian(n, [], kossakowski=(ops, Gam))
    lj = Lindbladian(n, [], eigenmode_jumps(ops, Gam))

    basis = ["ZII"]
    coeffs = np.array([1.0])
    leak_k = lk.leakage(basis, coeffs)
    leak_j = lj.leakage(basis, coeffs)
    assert leak_k, "pure-Kossakowski dissipator produced empty leakage"
    assert set(leak_k) == set(leak_j)
    for s in leak_k:
        assert abs(leak_k[s] - leak_j[s]) < 1e-12
