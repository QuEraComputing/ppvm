# Mirror of test_magic_states.py, but evolving a `GeneralizedTableauSum` and
# drawing all shots from its `Sampler` in one call, instead of re-running a
# single `GeneralizedTableau` trajectory per shot. The noise channels branch the
# sum once at build time; sampling then draws from the resulting mixture. The
# infidelity statistics must match the trajectory simulator.
#
# From tsim example
# https://bloqade.quera.com/latest/digital/examples/tsim/magic_state_distillation/

import numpy as np

from ppvm.generalized_tableau_sum import GeneralizedTableauSum

theta = -np.arccos(np.sqrt(1 / 3))  # Distillation angle
p = 0.05  # Noise probability


def test_simple_infidelity():
    tab = GeneralizedTableauSum(n_qubits=1, seed=0)

    # Prepare magic state SH
    tab.reset(0)
    tab.rx(0, theta=theta)
    tab.t_dag(0)

    # Add noise
    tab.depolarize1(0, p=p)

    # Undo prep without noise
    tab.t(0)
    tab.rx(0, theta=-theta)

    n_shots = 500_000
    shots = tab.sampler().sample_shots(n_shots)
    count_one = sum(shot[0] for shot in shots)

    infidelity = float(count_one) / n_shots
    print(f"Infidelity: {infidelity}")

    # NOTE: result from the example at the time of writing
    assert np.isclose(infidelity, 0.03355, atol=1e-3)


# NOTE: helper for below -- builds the circuit on the sum (no final measurement).
def build_distillation_sqrt(tab: GeneralizedTableauSum, q: list[int], noise: bool):
    # initial state - prepare noisy magic states on all 5 qubits
    for i in q:
        tab.reset(i)

    for qi in q:
        tab.rx(qi, theta=theta)
    for qi in q:
        tab.t_dag(qi)
    if noise:
        for qi in q:
            tab.depolarize1(qi, p=p)

    # distillation circuit
    for qi in [q[0], q[1], q[4]]:
        tab.sqrt_x(qi)
    for a, b in zip([q[0], q[2]], [q[1], q[3]]):
        tab.cz(a, b)
    for qi in [q[0], q[3]]:
        tab.sqrt_y(qi)
    for a, b in zip([q[0], q[3]], [q[2], q[4]]):
        tab.cz(a, b)
    tab.sqrt_x_dag(q[0])
    for a, b in zip([q[0], q[1]], [q[4], q[3]]):
        tab.cz(a, b)
    for qi in q:
        tab.sqrt_x_dag(qi)

    # undo magic state preparation on first qubit to measure infidelity
    tab.t(q[0])
    tab.rx(q[0], theta=-theta)


def build_distillation_rot(tab: GeneralizedTableauSum, q: list[int], noise: bool):
    # initial state - prepare noisy magic states on all 5 qubits
    for i in q:
        tab.reset(i)

    for qi in q:
        tab.rx(qi, theta=theta)
    for qi in q:
        tab.t_dag(qi)
    if noise:
        for qi in q:
            tab.depolarize1(qi, p=p)

    # distillation circuit
    for qi in [q[0], q[1], q[4]]:
        # tab.sqrt_x(qi)
        tab.rx(qi, theta=np.pi / 2)
    for a, b in zip([q[0], q[2]], [q[1], q[3]]):
        tab.cz(a, b)
    for qi in [q[0], q[3]]:
        # tab.sqrt_y(qi)
        tab.ry(qi, theta=np.pi / 2)
    for a, b in zip([q[0], q[3]], [q[2], q[4]]):
        tab.cz(a, b)
    # tab.sqrt_x_dag(q[0])
    tab.rx(q[0], theta=-np.pi / 2)
    for a, b in zip([q[0], q[1]], [q[4], q[3]]):
        tab.cz(a, b)
    for qi in q:
        # tab.sqrt_x_dag(qi)
        tab.rx(qi, theta=-np.pi / 2)

    # undo magic state preparation on first qubit to measure infidelity
    tab.t(q[0])
    tab.rx(q[0], theta=-theta)


def _accepted_infidelity(shots: list[list]) -> float:
    """Post-select on the syndrome [1, 0, 1, 1] and return P(qubit 0 == 1)."""
    count_one = 0
    count_accepted = 0
    for result in shots:
        if result[1:] != [True, False, True, True]:
            # wrong syndromes, don't select this sample
            continue
        count_accepted += 1
        # sel should now be in 0 since we undid the magic state before measuring
        count_one += result[0]
    return float(count_one) / count_accepted


def test_distillation_infidelity_sqrt():
    tab = GeneralizedTableauSum(5, seed=0)
    q = list(range(5))
    build_distillation_sqrt(tab, q, True)

    n_shots = 50_000
    shots = tab.sampler().sample_shots(n_shots)
    infidelity = _accepted_infidelity(shots)
    print(f"Sqrt Infidelity: {infidelity}")

    # NOTE: result from the example at the time of writing
    assert np.isclose(infidelity, 0.00683, atol=5 * 1e-3)


def test_distillation_infidelity_rot():
    tab = GeneralizedTableauSum(5, seed=0)
    q = list(range(5))
    build_distillation_rot(tab, q, True)

    n_shots = 50_000
    shots = tab.sampler().sample_shots(n_shots)
    infidelity = _accepted_infidelity(shots)
    print(f"Rot Infidelity: {infidelity}")

    # NOTE: result from the example at the time of writing
    assert np.isclose(infidelity, 0.00683, atol=1e-3)


def test_distillation_infidelity_sqrt_noiseless():
    tab = GeneralizedTableauSum(5, seed=0)
    q = list(range(5))
    build_distillation_sqrt(tab, q, False)

    shots = tab.sampler().sample_shots(10_000)
    infidelity = _accepted_infidelity(shots)

    assert np.isclose(infidelity, 0.0)


def test_distillation_infidelity_rot_noiseless():
    tab = GeneralizedTableauSum(5, seed=0)
    q = list(range(5))
    build_distillation_rot(tab, q, False)

    shots = tab.sampler().sample_shots(10_000)
    infidelity = _accepted_infidelity(shots)

    assert np.isclose(infidelity, 0.0)


def test_single_qubit_magic_state_noiseless():
    def build_single_qubit(tab: GeneralizedTableauSum, q: list[int]):
        for i in q:
            tab.reset(i)
        for qi in q:
            tab.rx(qi, theta=theta)
        for qi in q:
            tab.t_dag(qi)

        # distillation circuit
        for qi in [q[0]]:
            tab.sqrt_x(qi)
        for qi in [q[0]]:
            tab.sqrt_y(qi)
        tab.sqrt_x_dag(q[0])
        for qi in q:
            tab.sqrt_x_dag(qi)

        # undo magic state preparation on first qubit to measure infidelity
        tab.t(q[0])
        tab.rx(q[0], theta=-theta)

    def build_single_qubit_rots(tab: GeneralizedTableauSum, q: list[int]):
        for i in q:
            tab.reset(i)
        for qi in q:
            tab.rx(qi, theta=theta)
        for qi in q:
            tab.t_dag(qi)

        # distillation circuit
        for qi in [q[0]]:
            tab.rx(qi, theta=np.pi / 2)
        for qi in [q[0]]:
            tab.ry(qi, theta=np.pi / 2)
        tab.rx(q[0], theta=-np.pi / 2)
        for qi in q:
            tab.rx(qi, theta=-np.pi / 2)

        # undo magic state preparation on first qubit to measure infidelity
        tab.t(q[0])
        tab.rx(q[0], theta=-theta)

    tab = GeneralizedTableauSum(n_qubits=1, seed=0)
    build_single_qubit_rots(tab, [0])
    shots = tab.sampler().sample_shots(100)
    assert all(not shot[0] for shot in shots)

    tab = GeneralizedTableauSum(n_qubits=1, seed=0)
    build_single_qubit(tab, [0])
    shots = tab.sampler().sample_shots(100)
    assert all(not shot[0] for shot in shots)
