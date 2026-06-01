"""Statistical and structural tests for `TableauSumSampler`.

The sampler draws shots from a `GeneralizedTableauSum` (a probability-weighted
collection of stabilizer branches). The reference is the pure
`GeneralizedTableau` trajectory simulator: re-prepare a fresh tableau per shot
and measure every qubit. In the many-shot limit both must yield the same joint
distribution over per-qubit `MeasurementResult` outcomes.

Each statistical test draws N shots from both backends with fixed seeds and
asserts the total variation distance (TVD) between the empirical distributions
is below a finite-sample threshold, matching the Rust `sampler_vs_pure.rs`
suite. Tolerances and shot counts are taken to be comfortably above ~5 sigma.
"""

import math

from ppvm import GeneralizedTableau, MeasurementResult
from ppvm.generalized_tableau_sum import GeneralizedTableauSum

# Deterministic seeds so the suite never flakes.
SEED_SUM = 0xC0FFEE
SEED_PURE = 0xDEADBEEF


def _frequencies(shots):
    """Empirical distribution over shots, keyed by the per-qubit outcome tuple."""
    n = len(shots)
    counts: dict[tuple[MeasurementResult, ...], float] = {}
    for shot in shots:
        key = tuple(shot)
        counts[key] = counts.get(key, 0.0) + 1.0 / n
    return counts


def _tvd(a, b):
    keys = set(a) | set(b)
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in keys)


def _sum_shots(n_qubits, circuit, shots, seed=SEED_SUM):
    tab = GeneralizedTableauSum(n_qubits, seed=seed)
    circuit(tab)
    return tab.sampler().sample_shots(shots)


def _pure_shots(n_qubits, circuit, shots, seed=SEED_PURE):
    """Reference: one fresh single-tableau trajectory per shot."""
    out = []
    for i in range(shots):
        tab = GeneralizedTableau(n_qubits, seed=seed + i)
        circuit(tab)
        out.append([tab.measure(q) for q in range(n_qubits)])
    return out


def _assert_distributions_match(n_qubits, circuit, shots, tol, label):
    sum_shots = _sum_shots(n_qubits, circuit, shots)
    pure_shots = _pure_shots(n_qubits, circuit, shots)
    d = _tvd(_frequencies(sum_shots), _frequencies(pure_shots))
    assert d < tol, f"[{label}] TVD = {d:.4f} >= tol {tol}"


# ---------------------------------------------------------------------------
# Structural / shape
# ---------------------------------------------------------------------------


def test_sample_shape_and_type():
    tab = GeneralizedTableauSum(3, seed=1)
    tab.h(0)
    sampler = tab.sampler()
    shot = sampler.sample()
    assert len(shot) == 3
    assert all(isinstance(r, MeasurementResult) for r in shot)


def test_sample_shots_shape_and_type():
    tab = GeneralizedTableauSum(3, seed=1)
    tab.h(0)
    sampler = tab.sampler()
    shots = sampler.sample_shots(20)
    assert len(shots) == 20
    assert all(len(s) == 3 for s in shots)
    assert all(isinstance(r, MeasurementResult) for s in shots for r in s)


def test_raw_sample_returns_plain_ints():
    tab = GeneralizedTableauSum(2, seed=1)
    tab.h(0)
    sampler = tab.sampler()
    raw = sampler.raw_sample()
    assert len(raw) == 2
    assert all(isinstance(r, int) and not isinstance(r, MeasurementResult) for r in raw)


def test_raw_shots_shape():
    tab = GeneralizedTableauSum(2, seed=1)
    tab.h(0)
    sampler = tab.sampler()
    raw = sampler.raw_shots(10)
    assert len(raw) == 10
    assert all(len(s) == 2 for s in raw)


def test_sample_shots_zero_is_empty():
    tab = GeneralizedTableauSum(2, seed=1)
    sampler = tab.sampler()
    assert sampler.sample_shots(0) == []


# ---------------------------------------------------------------------------
# Determinism / reproducibility
# ---------------------------------------------------------------------------


def test_same_seed_reproducible():
    # Identical full runs (same tableau seed, one sampler each) must agree.
    def run():
        tab = GeneralizedTableauSum(3, seed=12345)
        tab.h(0)
        tab.cnot(0, 1)
        tab.cnot(0, 2)
        return tab.sampler().sample_shots(200)

    assert run() == run()


# ---------------------------------------------------------------------------
# Deterministic circuits (every shot identical) -- also exercise the parallel
# path: N=200 > 4*n_threads, so sample_shots dispatches to rayon. Correct
# results here confirm the parallel branch produces sound samples.
# ---------------------------------------------------------------------------


def test_x_gate_all_ones_parallel_path():
    tab = GeneralizedTableauSum(1, seed=1)
    tab.x(0)
    shots = tab.sampler().sample_shots(200)
    assert all(s[0] == MeasurementResult.ONE for s in shots)


def test_zero_state_all_zeros():
    tab = GeneralizedTableauSum(2, seed=1)
    shots = tab.sampler().sample_shots(200)
    assert all(s == [MeasurementResult.ZERO, MeasurementResult.ZERO] for s in shots)


def test_rx_pi_flips_to_one():
    tab = GeneralizedTableauSum(1, seed=1)
    tab.rx(0, math.pi)
    shots = tab.sampler().sample_shots(200)
    assert all(s[0] == MeasurementResult.ONE for s in shots)


def test_loss_certain_marks_lost():
    tab = GeneralizedTableauSum(1, seed=1)
    tab.x(0)
    tab.loss_channel(0, 1.0)
    shots = tab.sampler().sample_shots(200)
    assert all(s[0] == MeasurementResult.LOST for s in shots)


# ---------------------------------------------------------------------------
# Sampler vs pure trajectory: noiseless Clifford
# ---------------------------------------------------------------------------


def test_bell_state_correlated():
    def circuit(t):
        t.h(0)
        t.cnot(0, 1)

    # Outcomes are correlated (00 / 11 only), each ~0.5.
    sum_shots = _sum_shots(2, circuit, 8000)
    assert all(s[0] == s[1] for s in sum_shots)
    _assert_distributions_match(2, circuit, 8000, 0.04, "bell_state_correlated")


def test_ghz_state_correlated():
    def circuit(t):
        t.h(0)
        t.cnot(0, 1)
        t.cnot(0, 2)

    sum_shots = _sum_shots(3, circuit, 8000)
    assert all(s[0] == s[1] == s[2] for s in sum_shots)
    _assert_distributions_match(3, circuit, 8000, 0.04, "ghz_state_correlated")


def test_rx_half_pi_unbiased():
    def circuit(t):
        t.rx(0, math.pi / 2)

    _assert_distributions_match(1, circuit, 8000, 0.04, "rx_half_pi_unbiased")


# ---------------------------------------------------------------------------
# Sampler vs pure trajectory: noise channels
# ---------------------------------------------------------------------------


def test_depolarize_on_ground_state():
    p = 0.6

    def circuit(t):
        t.depolarize(0, p)

    sum_shots = _sum_shots(1, circuit, 8000)
    ones = sum(1 for s in sum_shots if s[0] == MeasurementResult.ONE) / len(sum_shots)
    expected = 2.0 * p / 3.0
    assert abs(ones - expected) < 0.04, f"P(1)={ones:.4f}, expected {expected:.4f}"
    _assert_distributions_match(1, circuit, 8000, 0.04, "depolarize_on_ground_state")


def test_loss_channel_after_hadamard():
    # Three outcomes: ZERO, ONE, LOST.
    def circuit(t):
        t.h(0)
        t.loss_channel(0, 0.3)

    _assert_distributions_match(1, circuit, 8000, 0.04, "loss_channel_after_hadamard")


def test_bell_pair_with_loss_on_q0():
    def circuit(t):
        t.h(0)
        t.cnot(0, 1)
        t.loss_channel(0, 0.3)

    _assert_distributions_match(2, circuit, 8000, 0.05, "bell_pair_with_loss_on_q0")


def test_pauli_error_nonuniform():
    p = [0.15, 0.25, 0.35]

    def circuit(t):
        t.pauli_error(0, p)

    sum_shots = _sum_shots(1, circuit, 8000)
    ones = sum(1 for s in sum_shots if s[0] == MeasurementResult.ONE) / len(sum_shots)
    expected = p[0] + p[1]  # X and Y flip the Z-basis outcome
    assert abs(ones - expected) < 0.04, f"P(1)={ones:.4f}, expected {expected:.4f}"
    _assert_distributions_match(1, circuit, 8000, 0.04, "pauli_error_nonuniform")


def test_bell_pair_with_depolarize_breaks_correlation():
    p = 0.3

    def circuit(t):
        t.h(0)
        t.cnot(0, 1)
        t.depolarize(0, p)

    sum_shots = _sum_shots(2, circuit, 8000)
    same = sum(1 for s in sum_shots if s[0] == s[1]) / len(sum_shots)
    expected = 1.0 - 2.0 * p / 3.0
    assert abs(same - expected) < 0.04, f"P(same)={same:.4f}, expected {expected:.4f}"
    _assert_distributions_match(2, circuit, 8000, 0.05, "bell_pair_depolarize")


def test_ghz_three_qubits_with_per_qubit_noise():
    def circuit(t):
        t.h(0)
        t.cnot(0, 1)
        t.cnot(0, 2)
        for q in range(3):
            t.depolarize(q, 0.1)
            t.loss_channel(q, 0.05)

    _assert_distributions_match(3, circuit, 8000, 0.08, "ghz_per_qubit_noise")


# ---------------------------------------------------------------------------
# Sampler vs pure trajectory: non-Clifford
# ---------------------------------------------------------------------------


def test_t_gate_distribution():
    # H, T, H on |0> gives non-trivial Z-basis probabilities; both backends
    # must agree.
    def circuit(t):
        t.h(0)
        t.t(0)
        t.h(0)

    _assert_distributions_match(1, circuit, 8000, 0.04, "t_gate_distribution")


def test_rotation_sequence_with_depolarize():
    def circuit(t):
        t.ry(0, 0.41 * math.pi)
        t.rz(0, 0.23 * math.pi)
        t.ry(0, 0.17 * math.pi)
        t.depolarize(0, 0.12)

    _assert_distributions_match(1, circuit, 8000, 0.04, "rotation_seq_depolarize")


# ---------------------------------------------------------------------------
# Reset
# ---------------------------------------------------------------------------


def test_reset_after_hadamard_collapses_to_zero():
    def circuit(t):
        t.h(0)
        t.reset(0)

    shots = _sum_shots(1, circuit, 2000)
    assert all(s[0] == MeasurementResult.ZERO for s in shots)


def test_reset_bell_pair_decorrelates():
    def circuit(t):
        t.h(0)
        t.cnot(0, 1)
        t.reset(0)

    sum_shots = _sum_shots(2, circuit, 8000)
    assert all(s[0] == MeasurementResult.ZERO for s in sum_shots)
    _assert_distributions_match(2, circuit, 8000, 0.04, "reset_bell_pair")
