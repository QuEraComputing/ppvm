import math

from bloqade.decoders.dialects.annotate.types import MeasurementResultValue
from kirin.dialects import ilist

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

    # The CX chain entangles q[0..84]; measuring q[80] collapses every
    # qubit to the same value. The conditional `x(q[81])` then flips q[81]
    # to ZERO when q[80] is ONE, so q[81] is always ZERO and every other
    # qubit agrees with q[80].
    for seed in range(8):
        result = _run_kernel(main, 85, seed=seed)
        assert result[81] == ZERO
        bulk = [result[i] for i in range(85) if i != 81]
        assert all(m == bulk[0] for m in bulk)


def test_basis_state_gates_return_expected_measurements():
    @squin.kernel
    def main():
        q = squin.qalloc(6)

        squin.x(q[1])
        squin.y(q[2])

        squin.h(q[3])
        squin.z(q[3])
        squin.h(q[3])

        squin.rx(math.pi, q[4])
        squin.ry(math.pi, q[5])

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


def test_single_qubit_clifford_gates_with_adjoints():
    @squin.kernel
    def main():
        q = squin.qalloc(8)

        squin.s(q[0])

        squin.h(q[1])
        squin.s(q[1])
        squin.s(q[1])
        squin.h(q[1])

        squin.h(q[2])
        squin.s(q[2])
        squin.s_adj(q[2])
        squin.h(q[2])

        squin.sqrt_x(q[3])
        squin.sqrt_x(q[3])

        squin.sqrt_x(q[4])
        squin.sqrt_x_adj(q[4])

        squin.sqrt_y(q[5])
        squin.sqrt_y(q[5])

        squin.sqrt_y(q[6])
        squin.sqrt_y_adj(q[6])

        squin.h(q[7])
        squin.t(q[7])
        squin.t_adj(q[7])
        squin.h(q[7])

        return squin.broadcast.measure(q)

    assert _run_kernel(main, 8) == [ZERO, ONE, ZERO, ONE, ZERO, ONE, ZERO, ZERO]


def test_rz_leaves_z_eigenstates_unchanged():
    @squin.kernel
    def main():
        q = squin.qalloc(3)

        squin.rz(0.0, q[0])
        squin.rz(0.5, q[1])

        squin.x(q[2])
        squin.rz(math.pi, q[2])

        return squin.broadcast.measure(q)

    assert _run_kernel(main, 3) == [ZERO, ZERO, ONE]


def test_controlled_y_and_z_gates():
    @squin.kernel
    def main():
        q = squin.qalloc(6)

        squin.x(q[0])
        squin.cy(q[0], q[1])

        squin.x(q[2])
        squin.cz(q[2], q[3])

        squin.x(q[4])
        squin.x(q[5])
        squin.cz(q[4], q[5])

        return squin.broadcast.measure(q)

    assert _run_kernel(main, 6) == [ONE, ONE, ONE, ZERO, ONE, ONE]


def test_u3_and_phased_xz_gates():
    @squin.kernel
    def main():
        q = squin.qalloc(6)

        squin.u3(0.0, 0.0, 0.0, q[0])
        squin.u3(math.pi, 0.0, 0.0, q[1])
        squin.u3(0.0, math.pi, 0.0, q[2])

        squin.phased_xz(0.0, 0.0, 0.0, q[3])
        squin.phased_xz(math.pi, 0.0, 0.0, q[4])
        squin.phased_xz(0.0, math.pi, 0.0, q[5])

        return squin.broadcast.measure(q)

    assert _run_kernel(main, 6) == [ZERO, ONE, ZERO, ZERO, ONE, ZERO]


def test_depolarizing_channels_at_zero_probability_preserve_state():
    @squin.kernel
    def main():
        q = squin.qalloc(3)

        squin.depolarize(0.0, q[0])
        squin.depolarize2(0.0, q[1], q[2])

        return squin.broadcast.measure(q)

    assert _run_kernel(main, 3) == [ZERO, ZERO, ZERO]


def test_two_qubit_pauli_channel_applies_selected_pauli():
    @squin.kernel
    def main():
        q = squin.qalloc(4)

        zero_probs = ilist.IList(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )
        squin.two_qubit_pauli_channel(zero_probs, q[0], q[1])

        # Index 3 corresponds to the XI Pauli product in the order
        # {IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ}, so the
        # control qubit gets flipped while the target qubit is untouched.
        xi_probs = ilist.IList(
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )
        squin.two_qubit_pauli_channel(xi_probs, q[2], q[3])

        return squin.broadcast.measure(q)

    assert _run_kernel(main, 4) == [ZERO, ZERO, ONE, ZERO]


def test_correlated_qubit_loss_loses_both_qubits_in_pair():
    @squin.kernel
    def main():
        q = squin.qalloc(4)

        squin.correlated_qubit_loss(1.0, [q[0], q[1]])
        squin.correlated_qubit_loss(0.0, [q[2], q[3]])

        return squin.broadcast.measure(q)

    assert _run_kernel(main, 4) == [LOST, LOST, ZERO, ZERO]
