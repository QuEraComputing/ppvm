"""Unit tests for the regen-stim core helpers.

Tests the math (per_bit_sigma, within_tolerance) and a smoke test for
run_stim. ppvm-side smoke (run_ppvm, write_distribution_fixture) lives
in the integration tests but exercising it requires the ppvm wheel built.
"""

import math

from regen_stim import core


def test_per_bit_sigma_interior_probabilities():
    sigmas = core.per_bit_sigma([0.5, 0.1, 0.9], num_shots=100)
    assert math.isclose(sigmas[0], math.sqrt(0.5 * 0.5 / 100))
    assert math.isclose(sigmas[1], math.sqrt(0.1 * 0.9 / 100))
    assert math.isclose(sigmas[2], math.sqrt(0.9 * 0.1 / 100))


def test_per_bit_sigma_boundary_probabilities():
    sigmas = core.per_bit_sigma([0.0, 1.0], num_shots=100)
    assert math.isclose(sigmas[0], math.sqrt(1.0 / 100))
    assert math.isclose(sigmas[1], math.sqrt(1.0 / 100))


def test_within_tolerance_tight_match():
    assert core.within_tolerance([0.5], [0.5], test_num_shots=1024, tolerance_sigma=5.0)


def test_within_tolerance_drift():
    # 50% off at p=0.5 with N=1024: sigma = sqrt(0.25/1024) ≈ 0.0156.
    # 5*sigma ≈ 0.078. Drift 0.10 is outside tolerance.
    assert not core.within_tolerance(
        [0.6], [0.5], test_num_shots=1024, tolerance_sigma=5.0
    )


def test_within_tolerance_length_mismatch_is_failure():
    assert not core.within_tolerance(
        [0.5, 0.5], [0.5], test_num_shots=1024, tolerance_sigma=5.0
    )


def test_run_stim_smoke_x_then_m():
    ref = core.run_stim("X 0\nM 0\n", num_shots=64)
    assert ref.bit_means == [1.0]
    assert ref.num_shots == 64
    assert ref.stim_version


def test_run_stim_smoke_h_then_m_is_random():
    ref = core.run_stim("H 0\nM 0\n", num_shots=10_000, seed=42)
    assert 0.4 < ref.bit_means[0] < 0.6


def test_max_qubit_in_source_indexes_correctly():
    assert core._max_qubit_in_source("M 0\n") == 0
    assert core._max_qubit_in_source("CX 0 7\nM 0 7\n") == 7
