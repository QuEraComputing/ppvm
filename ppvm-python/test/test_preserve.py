# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Tests for observable-aware (preserve-set) truncation."""

import math

import pytest

from ppvm import PauliSum, preserve_single_z

# =============================================================================
# Helper: `preserve_single_z` returns the right strings.
# =============================================================================


def test_preserve_single_z_helper():
    assert preserve_single_z(1) == ["Z"]
    assert preserve_single_z(2) == ["ZI", "IZ"]
    assert preserve_single_z(4) == ["ZIII", "IZII", "IIZI", "IIIZ"]


# =============================================================================
# Plumbing: when preserve_strings is set, the field round-trips and `new()`
# accepts it.
# =============================================================================


def test_preserve_strings_round_trip():
    ps = PauliSum.new(
        3,
        "Z1",
        preserve_strings=preserve_single_z(3),
        preserve_threshold=1e-4,
        preserve_weight_lambda=0.0,
    )
    assert list(ps.preserve_strings) == ["ZII", "IZI", "IIZ"]
    assert ps.preserve_threshold == 1e-4
    assert ps.preserve_weight_lambda == 0.0


def test_preserve_strings_length_validated():
    with pytest.raises(ValueError, match="length n_qubits"):
        PauliSum.new(
            3,
            "Z1",
            preserve_strings=["ZI"],  # wrong length
            preserve_threshold=0.1,
        )


# =============================================================================
# Behavior: a preserved string with tiny coefficient survives truncation.
# =============================================================================


def test_preserved_string_survives_below_threshold():
    """A single-Z string with coefficient well below the cutoff must
    survive truncation. Without preserve, it would be dropped."""
    # Start with a tiny-coefficient Z0 plus a normal-coefficient XII; the
    # preserve-aware truncate should keep both — Z0 because it's preserved,
    # XII because its magnitude (0.5) is above the 1e-3 threshold.
    ps = PauliSum.new(
        3,
        [("Z0", 1e-8), ("X0", 0.5)],
        preserve_strings=preserve_single_z(3),
        preserve_threshold=1e-3,
    )
    # The XII below-threshold-non-preserved version: also include something
    # that should be dropped to make sure the policy does drop non-preserved
    # below-threshold things.
    ps2 = PauliSum.new(
        3,
        [("Z0", 1e-8), ("X0", 1e-6), ("X1", 0.5)],
        preserve_strings=preserve_single_z(3),
        preserve_threshold=1e-3,
    )
    # Manually trigger the truncation via a no-op-on-this-state gate that
    # auto-truncates: apply rx(addr0, 0.0) — does nothing arithmetically
    # but triggers the auto-truncate.
    ps.rx(0, 0.0)
    ps2.rx(0, 0.0)
    kept_ps = {t for t, _ in ps.terms}
    kept_ps2 = {t for t, _ in ps2.terms}
    assert "ZII" in kept_ps, "preserved tiny Z0 must survive (ps)"
    assert "XII" in kept_ps, "above-threshold XII must survive (ps)"
    assert "ZII" in kept_ps2, "preserved tiny Z0 must survive (ps2)"
    assert "IXI" in kept_ps2, "above-threshold X1 must survive (ps2)"
    assert "XII" not in kept_ps2, "below-threshold non-preserved X0 must be dropped"


# =============================================================================
# Transport diagnostic: <Σ_j Z_j(t) Z_i(0)>.sum() is conserved exactly with
# preserve, but drifts without.
# =============================================================================


def test_total_z_conservation_with_preserve_vs_without():
    """Reproduces a tiny-scale version of the user's main.py setup:
    propagate a localized Z_i under XY exchange + Z-dephasing, with
    aggressive truncation, and check that single-Z preserve substantially
    reduces the drift in `result.sum(axis=1)` versus the same run with
    plain `CoefficientThreshold` truncation.

    The drift is not zero in either case (small θ ⇒ direct loss
    dominates and preserve closes most of the gap; but some indirect loss
    via dropped off-diagonal terms always remains). The test asserts
    only that preserve reduces drift by an order of magnitude or more.
    """
    L = 8
    i = L // 2
    threshold = 0.02
    gamma = 0.5
    dt = 0.1
    steps = 6
    noise = (1 - math.exp(-gamma * dt)) / 2

    # Small angle so far-site Z_j coefficients drift toward the cutoff
    # — the regime where preserve actually helps.
    edges = [(a, a + 1, 0.08) for a in range(L - 1)]
    z_observables = [PauliSum.new(L, f"Z{j}") for j in range(L)]

    def evolve(ps):
        sums = []
        for _ in range(steps + 1):
            sums.append(sum(ps.overlap(zz) for zz in z_observables))
            for q in range(L):
                ps.pauli_error(q, [0.0, 0.0, noise])
            # XY exchange on each edge: rxx then ryy (both commute and
            # together preserve total Z magnetization).
            for a, b, th in reversed(edges):
                ps.rxx(a, b, th)
                ps.ryy(a, b, th)
            for q in range(L):
                ps.pauli_error(q, [0.0, 0.0, noise])
        return sums

    # Plain truncation.
    ps_plain = PauliSum.new(L, f"Z{i}", min_abs_coeff=threshold, max_pauli_weight=L)
    drift_plain = evolve(ps_plain)

    # Preserve-aware truncation with single-Z exempt.
    ps_pres = PauliSum.new(
        L,
        f"Z{i}",
        min_abs_coeff=threshold,
        max_pauli_weight=L,
        preserve_strings=preserve_single_z(L),
        preserve_threshold=threshold,
    )
    drift_pres = evolve(ps_pres)

    plain_drift = 1.0 - drift_plain[-1]
    pres_drift = abs(1.0 - drift_pres[-1])
    assert plain_drift > 1e-3, (
        f"sanity check: plain truncation should drift here "
        f"(got drift={plain_drift:.2e}); try lowering threshold"
    )
    assert pres_drift < plain_drift / 10, (
        f"preserve should reduce drift by >= 10x; "
        f"got preserve_drift={pres_drift:.2e}, plain_drift={plain_drift:.2e}"
    )


# =============================================================================
# Virtual DAOE: weight_lambda > 0 makes high-weight terms get truncated more
# aggressively, even when the base threshold would have kept them.
# =============================================================================


def test_weight_lambda_drops_high_weight_more_aggressively():
    """A weight-3 term at coefficient 0.05 should survive a uniform 0.01
    threshold, but should be dropped by a `λ = 1.0` weight-biased
    threshold (whose effective cutoff at weight 3 is 0.01·e^3 ≈ 0.2)."""
    # Build a PauliSum with several terms of varying weight.
    ps = PauliSum.new(
        4,
        [("X0", 0.05), ("X0X1", 0.05), ("X0X1X2", 0.05)],
        preserve_strings=[],  # no preserve; use only weight-biased cutoff
        preserve_threshold=0.01,
        preserve_weight_lambda=1.0,
    )
    # We can't pass preserve_strings=[] currently — the dataclass treats
    # falsy as "no preserve". So instead pass a dummy preserve set with
    # something we don't have, plus set the threshold. Re-do via a
    # weight-only style:
    ps = PauliSum.new(
        4,
        [("X0", 0.05), ("X0X1", 0.05), ("X0X1X2", 0.05)],
        preserve_strings=["IIII"],  # never-matching keep-set; pure weighted cutoff
        preserve_threshold=0.01,
        preserve_weight_lambda=1.0,
    )
    # Trigger truncate (no-op gate).
    ps.rx(3, 0.0)
    kept = {t for t, _ in ps.terms}
    # Weight 1: cutoff = 0.01·e ≈ 0.0272. 0.05 > cutoff → survives.
    assert "XIII" in kept, "weight-1 above weighted cutoff should survive"
    # Weight 2: cutoff = 0.01·e^2 ≈ 0.0739. 0.05 < cutoff → dropped.
    assert "XXII" not in kept, "weight-2 below weighted cutoff should be dropped"
    # Weight 3: cutoff = 0.01·e^3 ≈ 0.2008. 0.05 < cutoff → dropped.
    assert "XXXI" not in kept, "weight-3 below weighted cutoff should be dropped"


# =============================================================================
# Default behavior: when preserve_strings is None, the existing strategy
# (CoefficientThreshold + MaxPauliWeight) is used unchanged.
# =============================================================================


def test_no_preserve_uses_existing_strategy():
    """Without preserve_strings, behaviour is identical to before this
    change. Below-threshold strings are dropped uniformly."""
    ps = PauliSum.new(2, [("ZI", 0.5), ("XI", 1e-8)], min_abs_coeff=1e-3)
    ps.rx(1, 0.0)  # no-op gate triggers truncate
    kept = {t for t, _ in ps.terms}
    assert "ZI" in kept
    assert "XI" not in kept
