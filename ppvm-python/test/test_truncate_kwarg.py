# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0

"""Tests for the per-gate ``truncate: bool`` kwarg and explicit
``PauliSum.truncate()``.

The kwarg lets a user defer truncation across a sequence of commuting
gates so that small intermediate coefficients survive long enough to
contribute to a downstream gate's output. Default ``truncate=True``
preserves the historical behaviour.
"""

import math

from ppvm import PauliSum


def _total_z(ps: PauliSum, n: int) -> float:
    """Sum of expectation values of every single-Z observable."""
    return sum(ps.overlap(PauliSum.new(n, f"Z{j}")) for j in range(n))


# =============================================================================
# Backward compatibility: gates without the kwarg behave exactly as before.
# =============================================================================


def test_default_truncate_kwarg_is_true():
    """Calling gates with no ``truncate=`` kwarg matches the historical
    behaviour (truncate after every gate)."""
    n = 3
    ps = PauliSum.new(n, "Z1", min_abs_coeff=1e-12)
    ps.rxx(0, 1, 0.5)
    ps.ryy(0, 1, 0.5)
    ps_explicit = PauliSum.new(n, "Z1", min_abs_coeff=1e-12)
    ps_explicit.rxx(0, 1, 0.5, truncate=True)
    ps_explicit.ryy(0, 1, 0.5, truncate=True)
    assert dict(ps.terms) == dict(ps_explicit.terms)


# =============================================================================
# `truncate=False` defers the cut; `ps.truncate()` then applies it.
# =============================================================================


def test_truncate_false_then_explicit_matches_implicit():
    """A pair of gates with ``truncate=False`` followed by an explicit
    ``ps.truncate()`` yields the same final state as truncating between
    the gates — when the threshold is small enough that nothing actually
    gets dropped at either intermediate point."""
    n = 3
    ps_a = PauliSum.new(n, "Z1", min_abs_coeff=1e-12)
    ps_a.rxx(0, 1, 0.5, truncate=False)
    ps_a.ryy(0, 1, 0.5, truncate=False)
    ps_a.truncate()

    ps_b = PauliSum.new(n, "Z1", min_abs_coeff=1e-12)
    ps_b.rxx(0, 1, 0.5)
    ps_b.ryy(0, 1, 0.5)
    assert dict(ps_a.terms) == dict(ps_b.terms)


def test_truncate_false_keeps_intermediate_terms_alive():
    """``truncate=False`` leaves the post-gate state untruncated; the
    plain call (default ``truncate=True``) immediately runs the strategy
    and may drop terms. Checked here mid-sequence (before any explicit
    final truncate) so the difference is visible."""
    n = 3
    threshold = 0.5  # sits between sin(0.4)≈0.39 and cos(0.4)≈0.92

    ps_plain = PauliSum.new(n, "Z1", min_abs_coeff=threshold)
    ps_plain.rxx(0, 1, 0.4)
    # The w=2 cross term `XYI` has coefficient sin(0.4) ≈ 0.39 < 0.5 and
    # is dropped by the immediate strategy run; only `IZI` survives.
    assert dict(ps_plain.terms) == {"IZI": math.cos(0.4)}

    ps_def = PauliSum.new(n, "Z1", min_abs_coeff=threshold)
    ps_def.rxx(0, 1, 0.4, truncate=False)
    # The same call with truncate deferred leaves both `IZI` and the
    # w=2 cross term in place — they would only get dropped on a
    # subsequent explicit `truncate()`.
    keys = {k for k, _ in ps_def.terms}
    assert "IZI" in keys and len(keys) > 1, (
        "deferred-truncate state should contain more than just IZI"
    )


# =============================================================================
# Headline: `rxx(truncate=False) + ryy(truncate=False) + truncate()` is
# semantically equivalent to the (now-deleted) `exchange` gate, which
# called rxx and ryy back-to-back in Rust and let the PyO3 wrapper
# truncate once at the end.
# =============================================================================


def test_truncate_false_pair_reproduces_old_exchange():
    """The new idiom replaces the bespoke `exchange` gate. Build it both
    ways: (1) rxx + ryy with intermediate truncate disabled, (2) what
    `exchange` did internally (which is exactly that). Identical state."""
    n = 4
    threshold = 0.5

    # Path 1: rxx(truncate=False) + ryy(truncate=False) + truncate()
    ps1 = PauliSum.new(n, "Z2", min_abs_coeff=threshold)
    ps1.rxx(1, 2, 0.4, truncate=False)
    ps1.ryy(1, 2, 0.4, truncate=False)
    ps1.truncate()

    # Path 2: the same primitives at the PyO3 level — what exchange did
    # in Rust was rxx then ryy with the single PyO3-level truncate at the
    # end. We replicate it here via the kwarg.
    ps2 = PauliSum.new(n, "Z2", min_abs_coeff=threshold)
    ps2.rxx(1, 2, 0.4, truncate=False)
    ps2.ryy(1, 2, 0.4)  # default truncate=True ⇒ single truncate at end

    assert dict(ps1.terms) == dict(ps2.terms)


# =============================================================================
# Sanity: with a loose enough threshold, deferred truncate preserves
# total Σ Z exactly through an XY-exchange chain (because no Pauli is
# ever near the threshold).
# =============================================================================


def test_deferred_truncate_preserves_total_z_with_loose_threshold():
    """Apply rxx+ryy on each edge with the truncate deferred until the
    end of the pair. With a non-aggressive threshold, total Σ Z is
    conserved to FP precision."""
    n = 4
    ps = PauliSum.new(n, "Z2", min_abs_coeff=1e-12)
    for a in range(n - 1):
        ps.rxx(a, a + 1, 0.3, truncate=False)
        ps.ryy(a, a + 1, 0.3, truncate=False)
        ps.truncate()
    drift = abs(1.0 - _total_z(ps, n))
    assert drift < 1e-10, f"Σ Z conservation broken: drift={drift:.3e}"


# =============================================================================
# `truncate()` is callable on its own.
# =============================================================================


def test_truncate_method_drops_below_threshold_keys():
    """Explicit ``ps.truncate()`` drops keys below the configured cutoff.
    (The constructor does not auto-truncate, so this is the only way to
    apply the strategy to the initial state.)"""
    n = 2
    ps = PauliSum.new(n, [("ZI", 0.5), ("XI", 1e-8)], min_abs_coeff=1e-3)
    before = dict(ps.terms)
    assert "XI" in before, "constructor leaves below-threshold entries in place"
    ps.truncate()
    after = dict(ps.terms)
    assert "XI" not in after, "truncate() should drop XI below 1e-3"
    assert after["ZI"] == 0.5


# =============================================================================
# Sanity: noise channels also take the kwarg.
# =============================================================================


def test_noise_channels_accept_truncate_kwarg():
    """`pauli_error` accepts the kwarg; defaults match the previous
    behaviour."""
    n = 3
    p_z = (1 - math.exp(-2 * 0.5 * 0.1)) / 2
    ps = PauliSum.new(n, "Z1", min_abs_coeff=1e-10)
    ps.pauli_error(0, [0.0, 0.0, p_z])
    ps.pauli_error(1, [0.0, 0.0, p_z], truncate=False)
    ps.pauli_error(2, [0.0, 0.0, p_z], truncate=False)
    ps.truncate()
    ps_ref = PauliSum.new(n, "Z1", min_abs_coeff=1e-10)
    for q in range(n):
        ps_ref.pauli_error(q, [0.0, 0.0, p_z])
    assert dict(ps.terms) == dict(ps_ref.terms)
