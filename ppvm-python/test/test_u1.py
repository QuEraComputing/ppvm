# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Tests for U(1)-symmetric (z-magnetization-conserving) gate helpers."""

import math

import pytest

from ppvm import PauliSum


def _terms_dict(state: PauliSum) -> dict[str, float]:
    return dict(state.terms)


def _approx_equal(a: dict[str, float], b: dict[str, float], tol: float = 1e-10) -> bool:
    keys = set(a) | set(b)
    return all(abs(a.get(k, 0.0) - b.get(k, 0.0)) < tol for k in keys)


def test_exchange_exists_and_runs():
    """``PauliSum.exchange`` is callable and propagates a simple observable."""
    ps = PauliSum.new(2, "IZ")
    ps.exchange(0, 1, 0.3)
    # IZ is in the (1,1) sector, so the exchange must produce a non-trivial spread.
    assert len(ps) >= 1


def test_xyzz_exists_and_runs():
    """``PauliSum.xyzz`` is callable and propagates a simple observable."""
    ps = PauliSum.new(2, "IZ")
    ps.xyzz(0, 1, theta_xy=0.2, theta_zz=0.05)
    assert len(ps) >= 1


def test_exchange_matches_rxx_then_ryy():
    """``exchange(θ)`` is equivalent to ``rxx(θ)`` then ``ryy(θ)``."""
    terms = ["IZ", "ZI", "ZZ", "XY", "YX", "IX", "XI", "YZ"]
    angles = [0.0, 0.3, -0.7, math.pi / 3]
    for term in terms:
        for theta in angles:
            fused = PauliSum.new(2, term)
            composed = PauliSum.new(2, term)
            fused.exchange(0, 1, theta)
            composed.rxx(0, 1, theta)
            composed.ryy(0, 1, theta)
            assert _approx_equal(_terms_dict(fused), _terms_dict(composed)), (
                f"exchange disagrees on {term} at θ={theta}: "
                f"{_terms_dict(fused)} vs {_terms_dict(composed)}"
            )


def test_xyzz_matches_exchange_then_rzz():
    """``xyzz`` is equivalent to ``exchange`` followed by ``rzz``."""
    terms = ["IZ", "ZI", "ZZ", "XY", "YX", "XX", "YY"]
    pairs = [(0.0, 0.0), (0.3, 0.1), (-0.4, 0.7), (math.pi / 5, -math.pi / 7)]
    for term in terms:
        for theta_xy, theta_zz in pairs:
            combined = PauliSum.new(2, term)
            stepwise = PauliSum.new(2, term)
            combined.xyzz(0, 1, theta_xy=theta_xy, theta_zz=theta_zz)
            stepwise.exchange(0, 1, theta_xy)
            stepwise.rzz(0, 1, theta_zz)
            assert _approx_equal(_terms_dict(combined), _terms_dict(stepwise)), (
                f"xyzz disagrees on {term} at θ_xy={theta_xy}, θ_zz={theta_zz}: "
                f"{_terms_dict(combined)} vs {_terms_dict(stepwise)}"
            )


def test_exchange_zero_angle_is_identity():
    """``exchange(a, b, 0)`` leaves any observable unchanged."""
    for term in ["IZ", "ZI", "XY", "YX", "II", "ZZ"]:
        ps = PauliSum.new(2, term)
        before = _terms_dict(ps)
        ps.exchange(0, 1, 0.0)
        assert _approx_equal(_terms_dict(ps), before), (
            f"exchange(0) altered {term}: got {_terms_dict(ps)}"
        )


def test_exchange_preserves_total_z():
    """``Σ_i Z_i`` commutes with `XX + YY`; exchange propagation is a no-op."""
    ps = PauliSum.new(3, [("Z0", 1.0), ("Z1", 1.0), ("Z2", 1.0)])
    before = _terms_dict(ps)
    ps.exchange(0, 1, 0.41)
    ps.exchange(1, 2, -1.2)
    assert _approx_equal(_terms_dict(ps), before), (
        f"exchange perturbed total Z: {_terms_dict(ps)} vs {before}"
    )


def test_xyzz_preserves_total_z():
    """`XX + YY + ZZ` also commutes with `Σ_i Z_i`."""
    ps = PauliSum.new(3, [("Z0", 1.0), ("Z1", 1.0), ("Z2", 1.0)])
    before = _terms_dict(ps)
    ps.xyzz(0, 1, theta_xy=0.5, theta_zz=0.3)
    ps.xyzz(1, 2, theta_xy=-0.2, theta_zz=0.7)
    assert _approx_equal(_terms_dict(ps), before)


def test_apply_u1_trotter_step_matches_manual_application():
    """The high-level helper must reproduce the documented manual gate order."""
    edges = [(0, 1), (1, 2), (2, 3)]
    theta_xy = 0.21
    theta_zz = 0.07
    fields_z = [0.1, 0.0, -0.1, 0.05]

    auto = PauliSum.new(4, [("Z0", 1.0), ("Z1", -1.0), ("Z2", 0.5)])
    manual = PauliSum.new(4, [("Z0", 1.0), ("Z1", -1.0), ("Z2", 0.5)])

    auto.apply_u1_trotter_step(
        edges=edges,
        theta_xy=theta_xy,
        theta_zz=theta_zz,
        fields_z=fields_z,
    )
    for i, j in edges:
        manual.xyzz(i, j, theta_xy=theta_xy, theta_zz=theta_zz)
    for site, h in enumerate(fields_z):
        if h != 0.0:
            manual.rz(site, h)

    assert _approx_equal(_terms_dict(auto), _terms_dict(manual)), (
        f"u1 trotter helper disagrees with manual application:\n"
        f"auto={_terms_dict(auto)}\nmanual={_terms_dict(manual)}"
    )


def test_apply_u1_trotter_step_with_per_edge_couplings():
    """Per-edge `(i, j, J)` triples must scale each edge's angles by `J`."""
    base_xy = 0.3
    couplings = [1.0, 0.5, -2.0]
    edges = [(0, 1, couplings[0]), (1, 2, couplings[1]), (2, 3, couplings[2])]

    auto = PauliSum.new(4, "Z0")
    manual = PauliSum.new(4, "Z0")

    auto.apply_u1_trotter_step(edges=edges, theta_xy=base_xy)
    for i, j, J in edges:
        manual.exchange(i, j, base_xy * J)

    assert _approx_equal(_terms_dict(auto), _terms_dict(manual))


def test_apply_u1_trotter_step_per_edge_angles():
    """`theta_xy` may be a list of per-edge angles."""
    edges = [(0, 1), (1, 2)]
    angles = [0.4, -0.25]

    auto = PauliSum.new(3, "Z0")
    manual = PauliSum.new(3, "Z0")

    auto.apply_u1_trotter_step(edges=edges, theta_xy=angles)
    for k, (i, j) in enumerate(edges):
        manual.exchange(i, j, angles[k])

    assert _approx_equal(_terms_dict(auto), _terms_dict(manual))


def test_apply_u1_trotter_step_preserves_total_z():
    """A full Trotter slice still preserves any observable in the total-Z sector."""
    ps = PauliSum.new(4, [("Z0", 1.0), ("Z1", 1.0), ("Z2", 1.0), ("Z3", 1.0)])
    before = _terms_dict(ps)
    ps.apply_u1_trotter_step(
        edges=[(0, 1), (1, 2), (2, 3)],
        theta_xy=0.3,
        theta_zz=0.15,
        fields_z=[0.1, -0.2, 0.05, 0.0],
    )
    # `\\sum Z_k` is invariant under the full Hamiltonian, including the
    # diagonal Z field (which is itself in the total-Z algebra).
    assert _approx_equal(_terms_dict(ps), before), (
        f"trotter step perturbed total Z: {_terms_dict(ps)} vs {before}"
    )


def test_apply_u1_trotter_step_rejects_bad_lengths():
    """The helper validates per-edge / per-site list lengths."""
    ps = PauliSum.new(3, "Z0")
    with pytest.raises(ValueError, match="per-edge"):
        ps.apply_u1_trotter_step(edges=[(0, 1), (1, 2)], theta_xy=[0.1])
    with pytest.raises(ValueError, match="fields_z"):
        ps.apply_u1_trotter_step(edges=[(0, 1)], theta_xy=0.1, fields_z=[0.1, 0.0])


def test_xy_dynamics_against_xxz_evolution():
    """Two consecutive `xyzz` calls on the same edge with `θ_xy = 0`, `θ_zz ≠ 0`
    must commute with everything in the Z-diagonal sector — a direct check
    that the ZZ piece of `xyzz` matches a pure rzz."""
    ps_xyzz = PauliSum.new(2, "ZZ")
    ps_rzz = PauliSum.new(2, "ZZ")
    ps_xyzz.xyzz(0, 1, theta_xy=0.0, theta_zz=0.42)
    ps_rzz.rzz(0, 1, 0.42)
    assert _approx_equal(_terms_dict(ps_xyzz), _terms_dict(ps_rzz))


# =============================================================================
# Truncation-robust conservation hardening
# =============================================================================
#
# `apply_u1_trotter_step` propagates {I, Z}-polynomial observables exactly
# in the conserved sector, modulo per-gate floating-point ε. These tests
# pin that guarantee against aggressive truncation policies: as long as
# the cutoff (`min_abs_coeff` or `max_pauli_weight`) leaves headroom above
# the per-gate ε and below the conserved-coefficient magnitude (~1),
# conservation is preserved within the `1e-10` tolerance used by
# `_approx_equal`.


def test_apply_u1_trotter_step_total_z_under_aggressive_coefficient_truncation():
    """`Σ Z_i` survives many Trotter sweeps with `min_abs_coeff = 0.5`.

    The conserved coefficients are O(1) — well above the cutoff. The
    transient cross terms cancel to O(ε) ≪ 0.5 before truncation runs,
    so the cutoff drops only the ε residues and leaves `Σ Z_i` intact.
    """
    n = 5
    edges = [(i, i + 1) for i in range(n - 1)]
    ps = PauliSum.new(
        n,
        [(f"Z{i}", 1.0) for i in range(n)],
        min_abs_coeff=0.5,
    )
    expected = _terms_dict(ps)
    for _ in range(4):
        ps.apply_u1_trotter_step(
            edges=edges,
            theta_xy=0.3,
            theta_zz=0.1,
            fields_z=[0.05] * n,
        )
    assert _approx_equal(_terms_dict(ps), expected), (
        f"total Z drifted under aggressive coefficient truncation:\n"
        f"got={_terms_dict(ps)}\nwant={expected}"
    )


def test_apply_u1_trotter_step_total_z_under_max_pauli_weight_one():
    """`Σ Z_i` survives Trotter sweeps with `max_pauli_weight = 1`.

    Cross terms have weight 2 and are dropped by the discrete weight
    check, independent of any floating-point sensitivity in the cutoff.
    """
    n = 5
    edges = [(i, i + 1) for i in range(n - 1)]
    ps = PauliSum.new(
        n,
        [(f"Z{i}", 1.0) for i in range(n)],
        max_pauli_weight=1,
    )
    expected = _terms_dict(ps)
    for _ in range(3):
        ps.apply_u1_trotter_step(
            edges=edges,
            theta_xy=0.25,
            theta_zz=0.07,
            fields_z=[0.04] * n,
        )
    assert _approx_equal(_terms_dict(ps), expected), (
        f"total Z drifted under MaxPauliWeight=1 truncation:\n"
        f"got={_terms_dict(ps)}\nwant={expected}"
    )


def test_full_zz_correlator_under_aggressive_truncation():
    """`Σ_{i<j} Z_iZ_j` (the sum over *all* pairs, not just nearest
    neighbours) commutes with `XX+YY+ZZ` on every edge, so the full
    correlator survives Trotter sweeps under aggressive truncation. The
    nearest-neighbour-only sum is not conserved, so the full pair
    decomposition is what locks down `xyzz`'s ZZ semantics."""
    n = 4
    # All C(n, 2) = 6 unordered pairs: (0,1),(0,2),(0,3),(1,2),(1,3),(2,3).
    all_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    terms = [(f"Z{i}Z{j}", 1.0) for i, j in all_pairs]
    ps = PauliSum.new(n, terms, min_abs_coeff=0.5, max_pauli_weight=2)
    expected = _terms_dict(ps)
    edges = [(i, i + 1) for i in range(n - 1)]
    for _ in range(2):
        ps.apply_u1_trotter_step(
            edges=edges,
            theta_xy=0.21,
            theta_zz=0.07,
        )
    assert _approx_equal(_terms_dict(ps), expected), (
        f"full ZZ correlator drifted under aggressive truncation:\n"
        f"got={_terms_dict(ps)}\nwant={expected}"
    )
