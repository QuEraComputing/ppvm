import time

from bloqade.decoders.dialects.annotate.types import MeasurementResultValue

from bloqade import squin
from ppvm import GeneralizedTableauSimulator

ZERO = MeasurementResultValue.Zero
ONE = MeasurementResultValue.One
LOST = MeasurementResultValue.Lost


def _run_kernel(kernel, n_qubits: int, *, seed: int | None = None):
    options = {"min_abs_coeff": 1e-8}
    if seed is not None:
        options["seed"] = seed

    sim = GeneralizedTableauSimulator(n_qubits, options=options)
    return list(sim.task(kernel).run())


def test_basic_execution():
    @squin.kernel
    def main():
        q = squin.qalloc(85)

        squin.h(q[0])
        squin.t(q[0])

        for i in range(1, 85):
            squin.cx(q[0], q[i])

        squin.broadcast.t(q[1:5])

        m = squin.measure(q[80])

        if m == MeasurementResultValue.One:
            squin.x(q[81])

        return squin.broadcast.measure(q)

    options = {"min_abs_coeff": 1e-8}
    sim = GeneralizedTableauSimulator(85, options=options)
    task = sim.task(main)

    start = time.time()
    result = task.run()
    print(result)
    print(f"Runtime: {(time.time() - start) * 1e3} ms")

    assert result[81] == ZERO


def test_basis_state_gates_return_expected_measurements():
    @squin.kernel
    def main():
        q = squin.qalloc(6)

        squin.x(q[1])
        squin.y(q[2])

        squin.h(q[3])
        squin.z(q[3])
        squin.h(q[3])

        squin.rx(3.141592653589793, q[4])
        squin.ry(3.141592653589793, q[5])

        return squin.broadcast.measure(q)

    assert _run_kernel(main, 6) == [ZERO, ONE, ONE, ONE, ONE, ONE]


def test_bell_pair_measurements_are_correlated_after_collapse():
    @squin.kernel
    def main():
        q = squin.qalloc(4)

        squin.h(q[0])
        squin.cx(q[0], q[1])

        first = squin.measure(q[0])
        second = squin.measure(q[1])

        if first == MeasurementResultValue.One:
            squin.x(q[2])

        if second == MeasurementResultValue.One:
            squin.x(q[3])

        return squin.broadcast.measure(q)

    outcomes = {tuple(_run_kernel(main, 4, seed=seed)) for seed in range(16)}

    assert outcomes == {
        (ZERO, ZERO, ZERO, ZERO),
        (ONE, ONE, ONE, ONE),
    }


def test_reset_loss_and_deterministic_pauli_channels():
    @squin.kernel
    def main():
        q = squin.qalloc(5)

        squin.x(q[0])
        squin.reset(q[0])

        squin.x(q[1])
        measured_one = squin.measure(q[1])
        if measured_one == MeasurementResultValue.One:
            squin.reset(q[1])

        squin.qubit_loss(1.0, q[2])

        squin.bit_flip(1.0, q[3])

        squin.h(q[4])
        squin.single_qubit_pauli_channel(0.0, 0.0, 1.0, q[4])
        squin.h(q[4])

        return squin.broadcast.measure(q)

    assert _run_kernel(main, 5, seed=0) == [ZERO, ZERO, LOST, ONE, ONE]
