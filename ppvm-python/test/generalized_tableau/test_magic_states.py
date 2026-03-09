# From tsim example
# https://bloqade.quera.com/latest/digital/examples/tsim/magic_state_distillation/

import time

import numpy as np

from ppvm import GeneralizedTableau

theta = -np.arccos(np.sqrt(1 / 3))  # Distillation angle
p = 0.05  # Noise probability


def test_simple_infidelity():
    tab = GeneralizedTableau(n_qubits=1)

    def run_shot(tab: GeneralizedTableau):

        # Prepare magic state SH
        tab.reset(0)
        tab.rx(0, theta)
        tab.t_adj(0)

        # Add noise
        tab.depolarize(0, p)

        # Undo prep without noise
        tab.t(0)
        tab.rx(0, -theta)
        return tab.measure(0)

    # NOTE: need a lot of shots to ensure convergence, otherwise test becomes flaky
    n_shots = 500_000
    count_one = 0
    total_time = 0.0
    for _ in range(n_shots):
        start = time.time()
        result = run_shot(tab)
        total_time += time.time() - start
        count_one += result

    infidelity = float(count_one) / n_shots
    print(f"Total time for {n_shots} shots: {total_time} s")
    print(f"Average time per shot: {total_time / n_shots * 1e3} ms")
    print(f"Infidelity: {float(count_one) / n_shots}")

    # NOTE: result from the example at the time of writing
    assert np.isclose(
        infidelity, 0.03355, atol=1e-3
    )  # could increase tolerance for fewer shots


# NOTE: helper for below
def distillation_sqrt(tab: GeneralizedTableau, q: list[int], noise: bool):
    # initial state - prepare noisy magic states on all 5 qubits
    for i in q:
        tab.reset(i)

    for qi in q:
        tab.rx(qi, theta)
    for qi in q:
        tab.t_adj(qi)
    if noise:
        for qi in q:
            tab.depolarize(qi, p)

    # distillation circuit
    for qi in [q[0], q[1], q[4]]:
        tab.sqrt_x(qi)
    for a, b in zip([q[0], q[2]], [q[1], q[3]]):
        tab.cz(a, b)
    for qi in [q[0], q[3]]:
        tab.sqrt_y(qi)
    for a, b in zip([q[0], q[3]], [q[2], q[4]]):
        tab.cz(a, b)
    tab.sqrt_x_adj(q[0])
    for a, b in zip([q[0], q[1]], [q[4], q[3]]):
        tab.cz(a, b)
    for qi in q:
        tab.sqrt_x_adj(qi)

    # undo magic state preparation on first qubit to measure infidelity
    tab.t(q[0])
    tab.rx(q[0], -theta)

    return [tab.measure(qi) for qi in q]


def distillation_rot(tab: GeneralizedTableau, q: list[int], noise: bool):
    # initial state - prepare noisy magic states on all 5 qubits
    for i in q:
        tab.reset(i)

    for qi in q:
        tab.rx(qi, theta)
    for qi in q:
        tab.t_adj(qi)
    if noise:
        for qi in q:
            tab.depolarize(qi, p)

    # distillation circuit
    for qi in [q[0], q[1], q[4]]:
        # tab.sqrt_x(qi)
        tab.rx(qi, np.pi / 2)
    for a, b in zip([q[0], q[2]], [q[1], q[3]]):
        tab.cz(a, b)
    for qi in [q[0], q[3]]:
        # tab.sqrt_y(qi)
        tab.ry(qi, np.pi / 2)
    for a, b in zip([q[0], q[3]], [q[2], q[4]]):
        tab.cz(a, b)
    # tab.sqrt_x_adj(q[0])
    tab.rx(q[0], -np.pi / 2)
    for a, b in zip([q[0], q[1]], [q[4], q[3]]):
        tab.cz(a, b)
    for qi in q:
        # tab.sqrt_x_adj(qi)
        tab.rx(qi, -np.pi / 2)

    # undo magic state preparation on first qubit to measure infidelity
    tab.t(q[0])
    tab.rx(q[0], -theta)

    return [tab.measure(qi) for qi in q]


def test_distillation_infidelity_sqrt():
    tab = GeneralizedTableau(5)
    q = list(range(5))

    n_shots = 50_000
    count_one = 0
    count_accepted = 0
    total_time = 0.0
    for _ in range(n_shots):
        start = time.time()
        result = distillation_sqrt(tab, q, True)
        total_time += time.time() - start
        if result[1:] != [True, False, True, True]:
            # wrong syndromes, don't select this sample
            continue

        count_accepted += 1
        sel = result[0]

        # sel should now be in 0 since we undid the magic state before measuring
        count_one += sel

    infidelity = float(count_one) / count_accepted
    print(f"Total time for {n_shots} shots: {total_time} s")
    print(f"Average time per shot: {total_time / n_shots * 1e3} ms")
    print(f"Sqrt Infidelity: {infidelity}")

    # NOTE: result from the example at the time of writing
    assert np.isclose(
        infidelity, 0.00683, atol=5 * 1e-3
    )  # could increase tolerance for fewer shots


def test_distillation_infidelity_rot():
    tab = GeneralizedTableau(5)
    q = list(range(5))

    n_shots = 50_000
    count_one = 0
    count_accepted = 0
    total_time = 0.0
    for _ in range(n_shots):
        start = time.time()
        result = distillation_rot(tab, q, True)
        total_time += time.time() - start
        if result[1:] != [True, False, True, True]:
            # wrong syndromes, don't select this sample
            continue

        count_accepted += 1
        sel = result[0]

        # sel should now be in 0 since we undid the magic state before measuring
        count_one += sel

    infidelity = float(count_one) / count_accepted
    print(f"Total time for {n_shots} shots: {total_time} s")
    print(f"Average time per shot: {total_time / n_shots * 1e3} ms")
    print(f"Rot Infidelity: {infidelity}")

    # NOTE: result from the example at the time of writing
    assert np.isclose(
        infidelity, 0.00683, atol=5 * 1e-3
    )  # could increase tolerance for fewer shots


def test_distillation_infidelity_sqrt_noiseless():
    tab = GeneralizedTableau(5)
    q = list(range(5))

    n_shots = 10_000
    count_one = 0
    count_accepted = 0
    total_time = 0.0
    for _ in range(n_shots):
        start = time.time()
        result = distillation_sqrt(tab, q, False)
        total_time += time.time() - start
        if result[1:] != [True, False, True, True]:
            # wrong syndromes, don't select this sample
            continue

        count_accepted += 1
        sel = result[0]

        # sel should now be in 0 since we undid the magic state before measuring
        count_one += sel

    infidelity = float(count_one) / count_accepted

    assert np.isclose(infidelity, 0.0)


def test_distillation_infidelity_rot_noiseless():
    tab = GeneralizedTableau(5)
    q = list(range(5))

    n_shots = 10_000
    count_one = 0
    count_accepted = 0
    total_time = 0.0
    for _ in range(n_shots):
        start = time.time()
        result = distillation_rot(tab, q, False)
        total_time += time.time() - start
        if result[1:] != [True, False, True, True]:
            # wrong syndromes, don't select this sample
            continue

        count_accepted += 1
        sel = result[0]

        # sel should now be in 0 since we undid the magic state before measuring
        count_one += sel

    infidelity = float(count_one) / count_accepted

    assert np.isclose(infidelity, 0.0)


def test_single_qubit_magic_state_noiseless():
    def distillation_single_qubit(tab: GeneralizedTableau, q: list[int]):
        # initial state - prepare noisy magic states on all 5 qubits
        for i in q:
            tab.reset(i)

        for qi in q:
            tab.rx(qi, theta)
        for qi in q:
            tab.t_adj(qi)

        # distillation circuit
        for qi in [q[0]]:
            tab.sqrt_x(qi)
        for qi in [q[0]]:
            tab.sqrt_y(qi)
        tab.sqrt_x_adj(q[0])
        for qi in q:
            tab.sqrt_x_adj(qi)

        # undo magic state preparation on first qubit to measure infidelity
        tab.t(q[0])
        tab.rx(q[0], -theta)

        return [tab.measure(qi) for qi in q]

    def distillation_single_qubit_rots(tab: GeneralizedTableau, q: list[int]):
        # initial state - prepare noisy magic states on all 5 qubits
        for i in q:
            tab.reset(i)

        for qi in q:
            tab.rx(qi, theta)
        for qi in q:
            tab.t_adj(qi)

        # distillation circuit
        for qi in [q[0]]:
            tab.rx(qi, np.pi / 2)
        for qi in [q[0]]:
            tab.ry(qi, np.pi / 2)
        tab.rx(q[0], -np.pi / 2)
        for qi in q:
            tab.rx(qi, -np.pi / 2)

        # undo magic state preparation on first qubit to measure infidelity
        tab.t(q[0])
        tab.rx(q[0], -theta)

        return [tab.measure(qi) for qi in q]

    for i in range(100):
        tab = GeneralizedTableau(n_qubits=1)
        result = distillation_single_qubit_rots(tab, [0])
        assert not result[0]

    for i in range(100):
        tab = GeneralizedTableau(n_qubits=1)
        result = distillation_single_qubit(tab, [0])
        assert not result[0]


test_single_qubit_magic_state_noiseless()
test_distillation_infidelity_rot()
test_distillation_infidelity_sqrt()
