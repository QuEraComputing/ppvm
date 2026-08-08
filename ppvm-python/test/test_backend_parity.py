"""Backend-independent ABI contracts for the native extension."""

import copy
import json
import os
import subprocess
import sys

import pytest

from ppvm import GeneralizedTableau, PauliSum, StimProgram, _core, sample_stim
from ppvm.generalized_tableau_sum import GeneralizedTableauSum
from ppvm.paulisum import LossyPauliSum


def test_compiled_backend_matches_fresh_build_request():
    expected = os.environ.get("PPVM_EXPECT_BACKEND")
    backend_name = vars(_core)["backend_name"]
    assert backend_name() in {"legacy", "traits-2"}
    if expected is not None:
        assert backend_name() == expected


def test_native_class_inventory_is_stable():
    expected = {"StimProgram", "backend_name"}
    expected |= {f"PauliSumIndexMapFxHash{i}" for i in range(16)}
    expected |= {f"PauliSumLossIndexMapFxHash{i}" for i in range(16)}
    expected |= {f"GeneralizedTableau{i}" for i in range(1, 33)}
    expected |= {f"GeneralizedTableauSum{i}" for i in range(1, 33)}
    expected |= {f"TableauSumSampler{i}" for i in range(1, 33)}
    actual = {name for name in dir(_core) if not name.startswith("_")}
    assert actual == expected


@pytest.mark.parametrize("width", [8, 9, 64, 65, 128, 129, 200, 1024, 2048])
def test_width_boundaries_preserve_terms_and_rendering(width):
    first = "X" + "I" * (width - 1)
    last = "I" * (width - 1) + "Z"
    state = PauliSum([first, last], width, [2.0, -3.0])
    assert state.terms == [(first, 2.0), (last, -3.0)]
    assert str(state._interface) == f"2.000 * {first} + -3.000 * {last}"


def test_order_trace_preserve_and_loss_encoding():
    state = PauliSum(["XI", "IZ", "XI"], 2, [1.0, 2.0, 3.0])
    assert state.terms == [("XI", 4.0), ("IZ", 2.0)]
    assert state.trace("X0") == 4.0

    kept = PauliSum(["XI"], 2, [1e-12], min_abs_coeff=1.0, preserve_strings=["XI"])
    kept.truncate()
    assert kept.terms == [("XI", 1e-12)]

    lossy = LossyPauliSum(["LI"], 2)
    assert lossy.terms == [("LI", 1.0)]


def test_wide_coefficients_copy_fork_and_mixture_streams():
    tab = GeneralizedTableau(200, seed=7)
    randomized = range(192, 200)
    for qubit in randomized:
        tab.h(qubit)
        tab.t(qubit)
    coefficients = tab.coefficients()
    assert all(isinstance(index, int) for index in coefficients)
    assert max(coefficients) >= 1 << 128

    left = copy.copy(tab)
    right = copy.deepcopy(tab)
    left_stream = [left.measure(qubit) for qubit in randomized]
    right_stream = [right.measure(qubit) for qubit in randomized]
    original_stream = [tab.measure(qubit) for qubit in randomized]
    assert left_stream == right_stream == original_stream

    fork_source = GeneralizedTableau(8, seed=7)
    for qubit in range(8):
        fork_source.h(qubit)
        fork_source.t(qubit)
    fork_a = fork_source.fork(seed=11)
    fork_b = fork_source.fork(seed=11)
    assert [fork_a.measure(q) for q in range(8)] == [fork_b.measure(q) for q in range(8)]
    seeded_streams = {
        tuple(fork_source.fork(seed=seed).measure(q) for q in range(8)) for seed in range(8)
    }
    assert len(seeded_streams) > 1

    rotated = GeneralizedTableauSum(1, seed=13)
    rotated.r(0, 0.3, 0.7)
    rotated_values = {shot[0].value for shot in rotated.sampler().sample_shots(128)}
    assert rotated_values == {0, 1}

    mixture = GeneralizedTableauSum(2, seed=13)
    sampler = mixture.sampler()
    mixture.x(0)
    assert sampler.raw_shots(16) == [b"\x00\x00"] * 16


def test_stim_seed_is_independent_of_rayon_thread_count():
    code = """
import json
from ppvm import StimProgram, sample_stim
p = StimProgram.parse("H 0\\nM 0")
print(json.dumps(sample_stim(p, n_qubits=1, num_shots=128, seed=17)))
"""
    outputs = []
    for threads in ("1", "4"):
        env = os.environ.copy()
        env["RAYON_NUM_THREADS"] = threads
        result = subprocess.run(
            [sys.executable, "-c", code],
            check=True,
            capture_output=True,
            text=True,
            env=env,
        )
        outputs.append(json.loads(result.stdout))
    assert outputs[0] == outputs[1]


def test_measurement_and_trace_encoding():
    tab = GeneralizedTableau(1, seed=0)
    assert tab.trace("Z0") == 1.0
    tab.loss_channel(0, 1.0)
    assert tab.measure(0).value == 2
    assert sample_stim(StimProgram.parse("M 0"), 1, num_shots=1, seed=0) == [[0]]
