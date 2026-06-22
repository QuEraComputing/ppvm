# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Tests for observable-aware (preserve-set) truncation."""

import math

import pytest

from ppvm import PauliSum


def _single_z(n_qubits: int) -> list[str]:
    return ["".join("Z" if i == j else "I" for i in range(n_qubits)) for j in range(n_qubits)]


# =============================================================================
# Plumbing: when preserve_strings is set, the field round-trips and `new()`
# accepts it.
# =============================================================================


def test_preserve_strings_round_trip():
    ps = PauliSum.new(
        3,
        "Z1",
        preserve_strings=_single_z(3),
    )
    assert ps.preserve_strings is not None
    assert list(ps.preserve_strings) == ["ZII", "IZI", "IIZ"]


def test_preserve_strings_length_validated():
    with pytest.raises(ValueError, match="length n_qubits"):
        PauliSum.new(
            3,
            "Z1",
            preserve_strings=["ZI"],  # wrong length
        )


# =============================================================================
# Behavior: a preserved string with tiny coefficient survives truncation,
# regardless of which truncation strategy is active. The preserve mechanism
# is a post-filter that re-inserts dropped preserved keys after the strategy
# runs — it composes with any strategy.
# =============================================================================


def test_preserved_string_survives_coefficient_truncation():
    """`min_abs_coeff` would drop a tiny-coefficient single-Z, but the
    preserve mechanism puts it back."""
    ps = PauliSum.new(
        3,
        [("Z0", 1e-8), ("X0", 0.5), ("X1", 1e-8)],
        min_abs_coeff=1e-3,
        preserve_strings=_single_z(3),
    )
    # Trigger auto-truncate via a no-op gate.
    ps.rx(0, 0.0)
    kept = {t for t, _ in ps.terms}
    assert "ZII" in kept, "preserved tiny Z0 must survive"
    assert "XII" in kept, "above-threshold XII must survive"
    assert "IXI" not in kept, "below-threshold non-preserved IXI must be dropped"


def test_preserved_string_survives_weight_truncation():
    """`max_pauli_weight` would drop a high-weight string; if that string
    happens to be in the preserve set, it's restored."""
    # Build a 4-qubit sum where one term has weight 3 (X0X1X2 → "XXXI")
    # and we cap max_pauli_weight at 2. Without preserve, the weight-3
    # term is dropped. With preserve including it, it survives.
    ps = PauliSum.new(
        4,
        [("Z0", 1.0), ("X0X1X2", 0.7)],
        max_pauli_weight=2,
        preserve_strings=["XXXI"],
    )
    ps.rx(0, 0.0)  # no-op triggers truncate
    kept = {t for t, _ in ps.terms}
    assert "ZIII" in kept, "weight-1 Z0 must survive"
    assert "XXXI" in kept, "weight-3 XXXI is in preserve set and must survive"


def test_preserved_string_survives_combined_truncation():
    """The strategy combines coefficient *and* weight cuts; preserve still
    works orthogonally on top."""
    ps = PauliSum.new(
        4,
        [("Z0", 1e-8), ("X0X1X2", 1e-8), ("Y0", 0.5)],
        min_abs_coeff=1e-3,
        max_pauli_weight=2,
        preserve_strings=["ZIII", "XXXI"],  # one tiny-coef, one high-weight
    )
    ps.rx(0, 0.0)
    kept = {t for t, _ in ps.terms}
    assert "ZIII" in kept, "preserved tiny Z0 must survive coefficient cut"
    assert "XXXI" in kept, "preserved weight-3 XXXI must survive weight cut"
    assert "YIII" in kept, "above-threshold YIII must survive"


# =============================================================================
# Transport diagnostic: <Σ_j Z_j(t) Z_i(0)>.sum() is conserved exactly with
# preserve, but drifts without.
# =============================================================================


def test_total_z_conservation_with_preserve_vs_without():
    """Propagate a localized Z_i under XY exchange + Z-dephasing with
    aggressive truncation. Single-Z preserve substantially reduces the
    drift in `result.sum(axis=1)` versus the same run with plain
    `min_abs_coeff` truncation."""
    L = 8
    i = L // 2
    threshold = 0.02
    gamma = 0.5
    dt = 0.1
    steps = 6
    noise = (1 - math.exp(-gamma * dt)) / 2

    edges = [(a, a + 1, 0.08) for a in range(L - 1)]
    z_observables = [PauliSum.new(L, f"Z{j}") for j in range(L)]

    def evolve(ps):
        sums = []
        for _ in range(steps + 1):
            sums.append(sum(ps.overlap(zz) for zz in z_observables))
            for q in range(L):
                ps.pauli_error(q, [0.0, 0.0, noise])
            for a, b, th in reversed(edges):
                ps.rxx(a, b, th)
                ps.ryy(a, b, th)
            for q in range(L):
                ps.pauli_error(q, [0.0, 0.0, noise])
        return sums

    # Plain truncation.
    ps_plain = PauliSum.new(L, f"Z{i}", min_abs_coeff=threshold, max_pauli_weight=L)
    drift_plain = evolve(ps_plain)

    # Same strategy + preserve-set on top.
    ps_pres = PauliSum.new(
        L,
        f"Z{i}",
        min_abs_coeff=threshold,
        max_pauli_weight=L,
        preserve_strings=_single_z(L),
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
