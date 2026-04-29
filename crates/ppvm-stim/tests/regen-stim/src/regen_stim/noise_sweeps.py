"""generated/noise_sweeps/: per-channel parameter sweeps."""

from __future__ import annotations

from . import core

NUM_SHOTS = 4096

SINGLE_QUBIT_PROBS = [0.001, 0.01, 0.1, 0.5]
TWO_QUBIT_PROBS = [0.001, 0.01, 0.1, 0.5]
SMALL_BIG_PROBS = [0.001, 0.01, 0.5]
READOUT_PROBS = [0.001, 0.01, 0.5]

PAULI1_PARAM_SETS: list[tuple[float, float, float]] = [
    (0.05, 0.0, 0.0),
    (0.0, 0.0, 0.1),
    (0.05, 0.05, 0.05),
]

# 15-arg Pauli2 channel order: IX IY IZ XI XX XY XZ YI YX YY YZ ZI ZX ZY ZZ
PAULI2_PARAM_SETS: list[tuple[float, ...]] = [
    (0.01,) + (0.0,) * 14,
    (0.0,) * 7 + (0.01,) + (0.0,) * 7,
    tuple([1.0 / 16.0] * 15),
]


def _prob_suffix(p: float) -> str:
    return f"p{p:g}".replace("0.", "")


def fixture_source_per_qubit_channel(
    channel: str, args: str, n_qubits: int, repeat: int = 5
) -> str:
    qs = " ".join(str(i) for i in range(n_qubits))
    body = "\n".join(f"{channel}{args} {qs}" for _ in range(repeat))
    return f"H {qs}\n{body}\nM {qs}\n"


def fixture_source_pair_channel(
    channel: str, args: str, n_qubits: int, repeat: int = 5
) -> str:
    assert n_qubits % 2 == 0
    qs = " ".join(str(i) for i in range(n_qubits))
    body = "\n".join(f"{channel}{args} {qs}" for _ in range(repeat))
    return f"H {qs}\n{body}\nM {qs}\n"


def fixture_source_readout_noise(measure: str, p: float, n_qubits: int) -> str:
    qs = " ".join(str(i) for i in range(n_qubits))
    return f"X {qs}\n{measure}({p}) {qs}\n"


def _emit(name: str, src: str, paths: core.CorpusPaths, failures: list[str]) -> int:
    meta = core.FixtureMeta(
        name=name,
        category="generated/noise_sweeps",
        source=src,
        test_num_shots=NUM_SHOTS,
    )
    try:
        core.write_distribution_fixture(meta, paths)
        return 1
    except Exception as e:
        failures.append(f"{name}: {e}")
        return 0


def run() -> int:
    paths = core.CorpusPaths.default()
    failures: list[str] = []
    written = 0

    # DEPOLARIZE1: 4 probs × 2 qubit counts = 8.
    for n in (1, 4):
        for p in SINGLE_QUBIT_PROBS:
            src = fixture_source_per_qubit_channel("DEPOLARIZE1", f"({p})", n)
            written += _emit(f"depolarize1_n{n}_{_prob_suffix(p)}", src, paths, failures)

    # DEPOLARIZE2: 4 probs × 2 even qubit counts = 8.
    for n in (2, 4):
        for p in TWO_QUBIT_PROBS:
            src = fixture_source_pair_channel("DEPOLARIZE2", f"({p})", n)
            written += _emit(f"depolarize2_n{n}_{_prob_suffix(p)}", src, paths, failures)

    # PAULI_CHANNEL_1: 3 sets × 2 qubit counts = 6.
    for n in (1, 4):
        for i, params in enumerate(PAULI1_PARAM_SETS):
            args = "(" + ", ".join(str(x) for x in params) + ")"
            src = fixture_source_per_qubit_channel("PAULI_CHANNEL_1", args, n)
            written += _emit(f"pauli_channel_1_set{i}_n{n}", src, paths, failures)

    # PAULI_CHANNEL_2: 3 sets × 2 even qubit counts = 6.
    for n in (2, 4):
        for i, params in enumerate(PAULI2_PARAM_SETS):
            assert len(params) == 15
            args = "(" + ", ".join(str(x) for x in params) + ")"
            src = fixture_source_pair_channel("PAULI_CHANNEL_2", args, n)
            written += _emit(f"pauli_channel_2_set{i}_n{n}", src, paths, failures)

    # X_ERROR / Y_ERROR / Z_ERROR: 3 channels × 3 probs × 1 qubit = 9.
    for ch in ("X_ERROR", "Y_ERROR", "Z_ERROR"):
        for p in SMALL_BIG_PROBS:
            src = fixture_source_per_qubit_channel(ch, f"({p})", 1)
            written += _emit(f"{ch.lower()}_{_prob_suffix(p)}", src, paths, failures)

    # M(p) / MR(p): 2 measurements × 3 probs = 6.
    for measure in ("M", "MR"):
        for p in READOUT_PROBS:
            src = fixture_source_readout_noise(measure, p, n_qubits=1)
            written += _emit(
                f"{measure.lower()}_readout_{_prob_suffix(p)}", src, paths, failures
            )

    print(f"regen-stim noise-sweeps: wrote {written} fixtures")
    if failures:
        print("regen-stim noise-sweeps: failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0
