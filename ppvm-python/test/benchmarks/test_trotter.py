"""pytest-benchmark mirror of crates/ppvm-runtime/benches/trotter.rs.

Same circuit, parameters, truncation and initial observable as the Rust
trotter benchmark, run through the Python ``PauliSum``. Python hard-picks the
IndexMap + FxHash storage variant, so this corresponds to exactly the
``ByteF64FxIndexMap`` Rust config -- a single backend, by design.

Run the timed benchmark with:

    uv run pytest test/benchmarks/test_trotter.py --benchmark-enable

Without ``--benchmark-enable`` it runs once as a smoke test (see ``addopts`` in
pyproject.toml).
"""

import copy

import pytest

from ppvm import PauliSum

# parameters (match trotter.rs)
N_QUBITS = 12
H = 1.0
DT = 0.1 / H
TIME = 1.0 / H
J = 1.0 / 8.0 * H
MIN_ABS_COEFF = 1e-6
# pauli_error([p/4]*3) damps each X/Y/Z coeff by (1-p): the depolarizing
# channel for p = 1e-4, matching DepolarizingNoise(1e-4) in PP.jl.
NOISE = [1e-4 / 4.0] * 3

ROUNDS = 200


def build_state():
    # initial observable: sum_i Z_i
    terms = [f"Z{i}" for i in range(N_QUBITS)]
    return PauliSum.new(N_QUBITS, terms, min_abs_coeff=MIN_ABS_COEFF)


def trotter(state):
    steps = int(TIME / DT)
    theta_zz = DT * J
    theta_x = DT * H
    # gates default to truncate=True, i.e. truncate after every gate
    # application, to stay consistent with the Rust benchmark / PP.jl.
    for _ in range(steps):
        for i in range(N_QUBITS):
            state.rx(i, theta_x)
            state.pauli_error(i, NOISE)
        for i in range(N_QUBITS - 1):
            state.rzz(i, i + 1, theta_zz)
            state.pauli_error(i, NOISE)
            state.pauli_error(i + 1, NOISE)


@pytest.mark.benchmark(group="trotter")
def test_trotter(benchmark):
    state = build_state()
    # clone a fresh state each round (not timed), mirroring Rust's
    # iter_batched_ref(|| state.clone(), ...).
    benchmark.pedantic(
        trotter,
        setup=lambda: ((copy.copy(state),), {}),
        iterations=1,
        rounds=ROUNDS,
    )
