"""Noisy 85-qubit magic-state-distillation circuit on the summed tableau.

Python port of `crates/ppvm-tableau-sum/examples/msd-noisy.rs`. Unlike the
noiseless `demo/msd.py` (which runs a single `GeneralizedTableau` trajectory),
this evolves a `GeneralizedTableauSum` -- a probability-weighted collection of
stabilizer branches -- applying a loss channel and a depolarizing channel after
every gate. Each noise op branches the sum; the `sum_cutoff` prunes branches
whose weight falls below the threshold. Finally a `Sampler` draws shots from the
resulting mixture in parallel.

Circuit from Rafael:
https://www.notion.so/Simulating-85-qubit-MSD-circuit-using-stabilizer-rank-decomposition-and-pyzx-288f86eeff3c802fb262ef1cfa69dfae
"""

import time

from ppvm.generalized_tableau_sum import GeneralizedTableauSum

QUBITS_PER_CODE_BLOCK = 17


def noise(tab: GeneralizedTableauSum, q: int, p_loss: float, p_depolarize: float) -> None:
    """Apply the per-gate noise: a loss channel then a depolarizing channel."""
    tab.loss_channel(q, p_loss)
    tab.depolarize1(q, p=p_depolarize)


def encode(
    tab: GeneralizedTableauSum,
    qubits: list[int],
    p_loss: float,
    p_depolarize: float,
) -> None:
    if len(qubits) != QUBITS_PER_CODE_BLOCK:
        return

    for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16]:
        tab.sqrt_y(qubits[i])
        noise(tab, qubits[i], p_loss, p_depolarize)
    print(f"Branches: {len(tab)}")

    for i, j in [[1, 3], [7, 10], [12, 14], [13, 16]]:
        tab.cz(qubits[i], qubits[j])
        noise(tab, qubits[i], p_loss, p_depolarize)
        noise(tab, qubits[j], p_loss, p_depolarize)
    for i in [7, 16]:
        tab.sqrt_y_dag(qubits[i])
        noise(tab, qubits[i], p_loss, p_depolarize)
    for i, j in [[4, 7], [8, 10], [11, 14], [15, 16]]:
        tab.cz(qubits[i], qubits[j])
        noise(tab, qubits[i], p_loss, p_depolarize)
        noise(tab, qubits[j], p_loss, p_depolarize)
    for i in [4, 10, 14, 16]:
        tab.sqrt_y_dag(qubits[i])
        noise(tab, qubits[i], p_loss, p_depolarize)
    for i, j in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]]:
        tab.cz(qubits[i], qubits[j])
        noise(tab, qubits[i], p_loss, p_depolarize)
        noise(tab, qubits[j], p_loss, p_depolarize)
    for i in [3, 6, 9, 10, 12, 13]:
        tab.sqrt_y(qubits[i])
        noise(tab, qubits[i], p_loss, p_depolarize)
    for i, j in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]]:
        tab.cz(qubits[i], qubits[j])
        noise(tab, qubits[i], p_loss, p_depolarize)
        noise(tab, qubits[j], p_loss, p_depolarize)
    for i in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14]:
        tab.sqrt_y(qubits[i])
        noise(tab, qubits[i], p_loss, p_depolarize)
    for i, j in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]]:
        tab.cz(qubits[i], qubits[j])
        noise(tab, qubits[i], p_loss, p_depolarize)
        noise(tab, qubits[j], p_loss, p_depolarize)
    for i in [0, 2, 5, 6, 8, 10, 12]:
        tab.sqrt_y_dag(qubits[i])
        noise(tab, qubits[i], p_loss, p_depolarize)


def main() -> None:
    n_qubits = QUBITS_PER_CODE_BLOCK * 5  # 85
    p_loss = 1e-4
    p_depolarize = 1e-4
    sum_cutoff = 1e-7
    n_shots = 1000

    build_start = time.perf_counter()

    tab = GeneralizedTableauSum(n_qubits, min_abs_coeff=1e-10, sum_cutoff=sum_cutoff)
    qubit_addrs = list(range(n_qubits))

    # Split qubits into 5 code blocks of QUBITS_PER_CODE_BLOCK each.
    ql = [
        qubit_addrs[i * QUBITS_PER_CODE_BLOCK : (i + 1) * QUBITS_PER_CODE_BLOCK]
        for i in range(5)
    ]

    # Phase 1: encoding (H + T on the encoding qubit, then encode the block).
    for q in ql:
        encoding_qubit = q[7]

        tab.h(encoding_qubit)
        noise(tab, encoding_qubit, p_loss, p_depolarize)

        tab.t(encoding_qubit)
        noise(tab, encoding_qubit, p_loss, p_depolarize)

        encode(tab, q, p_loss, p_depolarize)

    print(f"Branches: {len(tab)}")

    # Phase 2: middle gates (sqrt_x / cz / sqrt_y / sqrt_x_dag layers).
    for i in [0, 1, 4]:
        for q in ql[i]:
            tab.sqrt_x(q)
            noise(tab, q, p_loss, p_depolarize)

    print(f"Branches: {len(tab)}")

    for control, target in zip(ql[0], ql[1]):
        tab.cz(control, target)
        noise(tab, control, p_loss, p_depolarize)
        noise(tab, target, p_loss, p_depolarize)
    for control, target in zip(ql[2], ql[3]):
        tab.cz(control, target)
        noise(tab, control, p_loss, p_depolarize)
        noise(tab, target, p_loss, p_depolarize)
    for q in ql[0]:
        tab.sqrt_y(q)
        noise(tab, q, p_loss, p_depolarize)
    for q in ql[3]:
        tab.sqrt_y(q)
        noise(tab, q, p_loss, p_depolarize)
    for control, target in zip(ql[0], ql[2]):
        tab.cz(control, target)
        noise(tab, control, p_loss, p_depolarize)
        noise(tab, target, p_loss, p_depolarize)
    for control, target in zip(ql[3], ql[4]):
        tab.cz(control, target)
        noise(tab, control, p_loss, p_depolarize)
        noise(tab, target, p_loss, p_depolarize)
    for q in ql[0]:
        tab.sqrt_x_dag(q)
        noise(tab, q, p_loss, p_depolarize)
    for control, target in zip(ql[0], ql[4]):
        tab.cz(control, target)
        noise(tab, control, p_loss, p_depolarize)
        noise(tab, target, p_loss, p_depolarize)
    for control, target in zip(ql[1], ql[3]):
        tab.cz(control, target)
        noise(tab, control, p_loss, p_depolarize)
        noise(tab, target, p_loss, p_depolarize)
    for block in ql:
        for q in block:
            tab.sqrt_x_dag(q)
            noise(tab, q, p_loss, p_depolarize)

    print(f"Branches: {len(tab)}")
    build_ms = (time.perf_counter() - build_start) * 1e3
    print(f"Build time: {build_ms:.0f} ms")

    sampler = tab.sampler()
    sample_start = time.perf_counter()
    sampler.sample_shots(n_shots)
    # sampler.raw_shots(1000)
    sample_us = (time.perf_counter() - sample_start) * 1e6
    print(f"Time to {n_shots} samples: {sample_us:.0f} us")
    print(f"Time per shot: {sample_us / n_shots * 1e3:.0f} ns")


if __name__ == "__main__":
    main()
