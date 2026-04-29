"""generated/codes/: stim gen sweeps over surface, repetition, and color codes."""

from __future__ import annotations

import stim

from . import core

# Harness-side tableau cap. Circuits whose max qubit index ≥ this are skipped
# because the cargo harness uses a fixed-size 64-qubit tableau (storage `<8>`
# in stim_corpus.rs). Bumping this requires the harness's storage parameter
# bumped in lockstep and every committed mean regenerated.
HARNESS_MAX_QUBITS = 64

CODES_TASKS: dict[str, list[str]] = {
    "surface_code": [
        "unrotated_memory_x",
        "unrotated_memory_z",
        "rotated_memory_x",
        "rotated_memory_z",
    ],
    "repetition_code": ["memory"],
    "color_code": ["memory_xyz"],
}

DISTANCES = [3, 5, 7]
ROUNDS = [1, 3, 5]
NOISE_VALUES: list[float | None] = [None, 0.001, 0.01]

PHASE1_SUPPORTED = {
    # Resets / single-qubit Cliffords / two-qubit Cliffords.
    "R", "RZ", "I", "X", "Y", "Z", "H", "H_XZ",
    "S", "S_DAG", "SQRT_Z", "SQRT_Z_DAG", "SQRT_X", "SQRT_X_DAG",
    "SQRT_Y", "SQRT_Y_DAG",
    "CX", "ZCX", "CNOT", "CY", "ZCY", "CZ", "ZCZ",
    # Measurements.
    "M", "MZ", "MR",
    # Noise.
    "DEPOLARIZE1", "DEPOLARIZE2", "PAULI_CHANNEL_1", "PAULI_CHANNEL_2",
    "X_ERROR", "Y_ERROR", "Z_ERROR", "I_ERROR",
    # Annotations.
    "DETECTOR", "OBSERVABLE_INCLUDE", "TICK", "QUBIT_COORDS", "SHIFT_COORDS",
    "REPEAT",
}

def _first_unsupported_in(circuit: stim.Circuit) -> str | None:
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            inner = _first_unsupported_in(inst.body_copy())
            if inner is not None:
                return inner
            continue
        if inst.name not in PHASE1_SUPPORTED:
            return inst.name
    return None


def first_unsupported_instruction(source: str) -> str | None:
    """Return the first phase-1-unsupported instruction name, or None.

    Walks the parsed circuit in source order, descending into REPEAT bodies.
    Why: ppvm's normalizer recurses into REPEAT in lexical order, so a regex
    over the raw text would miss instructions inside indented REPEAT blocks
    and mis-classify the awaiting_phase2_instruction field.
    """
    return _first_unsupported_in(stim.Circuit(source))


def gen_circuit(code: str, task: str, distance: int, rounds: int, noise: float | None) -> str:
    code_task = f"{code}:{task}"
    kwargs: dict = dict(distance=distance, rounds=rounds)
    if noise is not None:
        kwargs.update(
            after_clifford_depolarization=noise,
            before_round_data_depolarization=noise,
            before_measure_flip_probability=noise,
            after_reset_flip_probability=noise,
        )
    circuit = stim.Circuit.generated(code_task, **kwargs)
    return str(circuit) + "\n"


def fixture_name(code: str, task: str, distance: int, rounds: int, noise: float | None) -> str:
    noise_str = "noiseless" if noise is None else f"p{noise:g}".replace("0.", "")
    return f"{code}_{task}_d{distance}_r{rounds}_{noise_str}"


def shot_count_for(distance: int) -> int:
    return 64 if distance >= 7 else 256


def run() -> int:
    paths = core.CorpusPaths.default()
    failures: list[str] = []
    written: list[str] = []
    for code, tasks in CODES_TASKS.items():
        for task in tasks:
            for distance in DISTANCES:
                for rounds in ROUNDS:
                    # Stim's color_code generator requires rounds >= 2.
                    if code == "color_code" and rounds < 2:
                        continue
                    for noise in NOISE_VALUES:
                        try:
                            src = gen_circuit(code, task, distance, rounds, noise)
                        except Exception as e:
                            failures.append(
                                f"stim gen failed for {code}/{task}/d{distance}/r{rounds}: {e}"
                            )
                            continue
                        if core._max_qubit_in_source(src) + 1 > HARNESS_MAX_QUBITS:
                            continue
                        name = fixture_name(code, task, distance, rounds, noise)
                        unsupported_instr = first_unsupported_instruction(src)
                        if unsupported_instr is not None:
                            meta = core.FixtureMeta(
                                name=name,
                                category="unsupported",
                                source=src,
                                test_num_shots=0,
                            )
                            try:
                                core.write_unsupported_fixture(
                                    meta, paths, awaiting_phase2_instruction=unsupported_instr
                                )
                                written.append(f"unsupported/{name}")
                            except Exception as e:
                                failures.append(f"unsupported/{name}: {e}")
                            continue
                        test_shots = shot_count_for(distance)
                        meta = core.FixtureMeta(
                            name=name,
                            category="generated/codes",
                            source=src,
                            test_num_shots=test_shots,
                        )
                        try:
                            core.write_distribution_fixture(meta, paths)
                            written.append(f"generated/codes/{name}")
                        except Exception as e:
                            failures.append(f"generated/codes/{name}: {e}")
    print(f"regen-stim codes: wrote {len(written)} fixtures")
    if failures:
        print("regen-stim codes: failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0
