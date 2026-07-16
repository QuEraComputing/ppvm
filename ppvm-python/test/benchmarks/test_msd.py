"""pytest-benchmark mirror of crates/ppvm-tableau/benches/tableau-msd-fused.rs.

Builds the 85-qubit magic-state-distillation circuit (5 code blocks of 17
qubits) on a `GeneralizedTableau`, timing the circuit construction. Every gate
is *splatted* -- applied to a collection of qubits in a single call, the Python
equivalent of the Rust ``*_many`` / ``cz_block_pairs`` fused methods:
single-qubit gates broadcast over each target, two-qubit gates consume
consecutive pairs.

Two arms, mirroring the Rust bench's two ``bench_function`` calls:

* ``test_msd_fused`` -- circuit construction only (``msd_func_fused::<false>``).
* ``test_msd_fused_measure`` -- a full shot: construction plus a single
  ``measure_many`` readout of all 85 qubits (``msd_func_fused::<true>``). The
  difference between the two arms isolates the measurement cost.

Circuit from Rafael:
https://www.notion.so/Simulating-85-qubit-MSD-circuit-using-stabilizer-rank-decomposition-and-pyzx-288f86eeff3c802fb262ef1cfa69dfae

Run the timed benchmark with:

    uv run --project ppvm-python --group dev pytest ppvm-python/test/benchmarks/test_msd.py --benchmark-enable
Without ``--benchmark-enable`` it runs once as a smoke test (see ``addopts`` in
pyproject.toml).
"""

import pytest

from ppvm import GeneralizedTableau

QUBITS_PER_CODE_BLOCK = 17
N_BLOCKS = 5
N_QUBITS = QUBITS_PER_CODE_BLOCK * N_BLOCKS  # 85


def _at(qubits: list[int], idxs: list[int]) -> list[int]:
    """Map block-local indices to absolute qubit addresses."""
    return [qubits[i] for i in idxs]


def _pairs(qubits: list[int], index_pairs: list[tuple[int, int]]) -> list[int]:
    """Flatten block-local index pairs into a consecutive (a, b, ...) cz target list."""
    return [qubits[i] for pair in index_pairs for i in pair]


def encode(tab: GeneralizedTableau, qubits: list[int]) -> None:
    if len(qubits) not in (7, 17):
        raise ValueError(f"Unsupported number of qubits {len(qubits)}")

    if len(qubits) == 7:
        tab.sqrt_y_dag(_at(qubits, [0, 1, 2, 3, 4, 5]))
        tab.cz(_pairs(qubits, [(1, 2), (3, 4), (5, 6)]))
        tab.sqrt_y(qubits[6])
        tab.cz(_pairs(qubits, [(0, 3), (2, 5), (4, 6)]))
        tab.sqrt_y(_at(qubits, [2, 3, 4, 5, 6]))
        tab.cz(_pairs(qubits, [(0, 1), (2, 3), (4, 5)]))
        tab.sqrt_y(_at(qubits, [1, 2, 4]))
        return

    # len == 17
    tab.sqrt_y(_at(qubits, [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16]))
    tab.cz(_pairs(qubits, [(1, 3), (7, 10), (12, 14), (13, 16)]))
    tab.sqrt_y_dag(_at(qubits, [7, 16]))
    tab.cz(_pairs(qubits, [(4, 7), (8, 10), (11, 14), (15, 16)]))
    tab.sqrt_y_dag(_at(qubits, [4, 10, 14, 16]))
    tab.cz(_pairs(qubits, [(2, 4), (6, 8), (7, 9), (10, 13), (14, 16)]))
    tab.sqrt_y(_at(qubits, [3, 6, 9, 10, 12, 13]))
    tab.cz(_pairs(qubits, [(0, 2), (3, 6), (5, 8), (10, 12), (11, 13)]))
    tab.sqrt_y(_at(qubits, [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14]))
    tab.cz(_pairs(qubits, [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (12, 15)]))
    tab.sqrt_y_dag(_at(qubits, [0, 2, 5, 6, 8, 10, 12]))


def build_msd() -> GeneralizedTableau:
    """Construct the full splatted MSD circuit and return the tableau."""
    tab = GeneralizedTableau(N_QUBITS)
    qubit_addrs = list(range(N_QUBITS))

    # Split qubits into N_BLOCKS code blocks of QUBITS_PER_CODE_BLOCK each.
    ql = [
        qubit_addrs[i * QUBITS_PER_CODE_BLOCK : (i + 1) * QUBITS_PER_CODE_BLOCK]
        for i in range(N_BLOCKS)
    ]

    # Encoding: H + T on each block's encoding qubit, then encode the block.
    for block in ql:
        encoding_qubit = block[6] if len(block) == 7 else block[7]
        tab.h(encoding_qubit)
        tab.t(encoding_qubit)
        encode(tab, block)

    # Middle gates: sqrt_x / cz / sqrt_y / sqrt_x_dag layers, all splatted. The
    # cross-block CZ layers entangle two contiguous registers (constant offset),
    # so they use the word-fused cz_block instead of a per-pair cz -- this is the
    # Python analogue of the Rust bench's cz_block_pairs / _cross_word calls.
    block_len = QUBITS_PER_CODE_BLOCK
    tab.sqrt_x(ql[0])
    tab.sqrt_x(ql[1])
    tab.sqrt_x(ql[4])
    tab.cz_block(ql[0][0], ql[1][0], block_len)
    tab.cz_block(ql[2][0], ql[3][0], block_len)
    tab.sqrt_y(ql[0])
    tab.sqrt_y(ql[3])
    tab.cz_block(ql[0][0], ql[2][0], block_len)
    tab.cz_block(ql[3][0], ql[4][0], block_len)
    tab.sqrt_x_dag(ql[0])
    tab.cz_block(ql[0][0], ql[4][0], block_len)
    tab.cz_block(ql[1][0], ql[3][0], block_len)
    for block in ql:
        tab.sqrt_x_dag(block)

    return tab


def build_and_measure() -> list:
    """Full shot: construct the circuit, then read out all qubits at once."""
    tab = build_msd()
    return tab.measure_many(range(N_QUBITS))


@pytest.mark.benchmark(group="msd")
def test_msd_fused(benchmark):
    # Construct a fresh tableau and apply the whole splatted circuit each round,
    # mirroring Rust's iter_batched_ref(|| {}, |_| msd_func_fused::<false>()).
    benchmark(build_msd)


@pytest.mark.benchmark(group="msd")
def test_msd_fused_measure(benchmark):
    # Construction plus a single measure_many readout of all 85 qubits, mirroring
    # Rust's msd_func_fused::<true>(). One FFI call for the whole readout.
    benchmark(build_and_measure)
