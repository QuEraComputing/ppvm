import time

from ppvm import GeneralizedTableau

# from Rafael:  https://www.notion.so/Simulating-85-qubit-MSD-circuit-using-stabilizer-rank-decomposition-and-pyzx-288f86eeff3c802fb262ef1cfa69dfae?source=copy_link#28df86eeff3c80bfa087ed15bcf49b77

QUBITS_PER_CODE_BLOCK = 17


def encode(tab: GeneralizedTableau, qubits: list[int]) -> None:
    if len(qubits) not in (7, 17):
        raise ValueError(f"Unsupported number of qubits {len(qubits)}")

    for q in qubits:
        tab.reset(q)

    if len(qubits) == 7:
        for idx, q in enumerate(qubits):
            if idx == 6:
                continue
            tab.sqrt_y_adj(q)

        tab.cz(qubits[1], qubits[2])
        tab.cz(qubits[3], qubits[4])
        tab.cz(qubits[5], qubits[6])

        tab.sqrt_y(qubits[6])

        tab.cz(qubits[0], qubits[3])
        tab.cz(qubits[2], qubits[5])
        tab.cz(qubits[4], qubits[6])

        for idx, q in enumerate(qubits):
            if idx < 2:
                continue
            tab.sqrt_y(q)

        tab.cz(qubits[0], qubits[1])
        tab.cz(qubits[2], qubits[3])
        tab.cz(qubits[4], qubits[5])

        tab.sqrt_y(qubits[1])
        tab.sqrt_y(qubits[2])
        tab.sqrt_y(qubits[4])
        return

    # len == 17
    for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16]:
        tab.sqrt_y(qubits[i])

    for i, j in [[1, 3], [7, 10], [12, 14], [13, 16]]:
        tab.cz(qubits[i], qubits[j])
    for i in [7, 16]:
        tab.sqrt_y_adj(qubits[i])
    for i, j in [[4, 7], [8, 10], [11, 14], [15, 16]]:
        tab.cz(qubits[i], qubits[j])
    for i in [4, 10, 14, 16]:
        tab.sqrt_y_adj(qubits[i])
    for i, j in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]]:
        tab.cz(qubits[i], qubits[j])
    for i in [3, 6, 9, 10, 12, 13]:
        tab.sqrt_y(qubits[i])
    for i, j in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]]:
        tab.cz(qubits[i], qubits[j])
    for i in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14]:
        tab.sqrt_y(qubits[i])
    for i, j in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]]:
        tab.cz(qubits[i], qubits[j])
    for i in [0, 2, 5, 6, 8, 10, 12]:
        tab.sqrt_y_adj(qubits[i])


def main():
    n_qubits = QUBITS_PER_CODE_BLOCK * 5
    # print(f"Number of qubits: {n_qubits}")

    tab = GeneralizedTableau(n_qubits)
    qubit_addrs = list(range(n_qubits))

    # split qubits into 5 groups of QUBITS_PER_CODE_BLOCK
    ql = [
        qubit_addrs[i * QUBITS_PER_CODE_BLOCK : (i + 1) * QUBITS_PER_CODE_BLOCK]
        for i in range(5)
    ]
    assert len(ql) == 5

    t_gate_counter = 0
    for q in ql:
        encoding_qubit = q[6] if len(q) == 7 else q[7]
        tab.h(encoding_qubit)
        tab.t(encoding_qubit)
        encode(tab, q)
        t_gate_counter += 1

    for i in [0, 1, 4]:
        for q in ql[i]:
            tab.sqrt_x(q)

    for control, target in zip(ql[0], ql[1]):
        tab.cz(control, target)

    for control, target in zip(ql[2], ql[3]):
        tab.cz(control, target)

    for q in ql[0]:
        tab.sqrt_y(q)

    for q in ql[3]:
        tab.sqrt_y(q)

    for control, target in zip(ql[0], ql[2]):
        tab.cz(control, target)

    for control, target in zip(ql[3], ql[4]):
        tab.cz(control, target)

    for q in ql[0]:
        tab.sqrt_x_adj(q)

    for control, target in zip(ql[0], ql[4]):
        tab.cz(control, target)

    for control, target in zip(ql[1], ql[3]):
        tab.cz(control, target)

    for i in range(5):
        for q in ql[i]:
            tab.sqrt_x_adj(q)

    # print(f"# T gates: {t_gate_counter}")
    # print(f"2 ^ t: {2**t_gate_counter}")

    bit_string = "".join("1" if tab.measure(i) else "0" for i in range(n_qubits))
    return bit_string
    # print(bit_string)


if __name__ == "__main__":
    avg_time = 0.0
    total_time = 0.0
    n_shots = 1000
    for _ in range(n_shots):
        start = time.perf_counter()
        main()
        elapsed = time.perf_counter() - start
        total_time += elapsed
        avg_time += elapsed / n_shots

    print(f"Average per shot: {avg_time * 1e3} ms")
