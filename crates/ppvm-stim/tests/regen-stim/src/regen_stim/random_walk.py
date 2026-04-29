"""generated/random/: random-walk sequences of supported instructions."""

from __future__ import annotations

import random as pyrand

from . import core

CLIFFORD_GATES_1Q = ["H", "S", "S_DAG", "X", "Y", "Z", "SQRT_X", "SQRT_Y"]
CLIFFORD_GATES_2Q = ["CX", "CY", "CZ"]
NOISE_GATES_1Q: list[tuple[str, list[float]]] = [
    ("DEPOLARIZE1", [0.001, 0.01, 0.1]),
    ("X_ERROR", [0.001, 0.01, 0.1]),
    ("Y_ERROR", [0.001, 0.01, 0.1]),
    ("Z_ERROR", [0.001, 0.01, 0.1]),
]

# Sub-sampled grid (full 3 × 4 × 3 × 8 = 288 is too many; spec asks ~30–40).
REGIMES = ["clifford-only", "+noise", "+readout"]
QUBIT_COUNTS = [2, 4, 8]
LENGTHS = [10, 50]
SEEDS = [0, 1]


def _emit_1q_clifford(n_qubits: int, lines: list[str]) -> None:
    gate = pyrand.choice(CLIFFORD_GATES_1Q)
    q = pyrand.randrange(n_qubits)
    lines.append(f"{gate} {q}")


def _emit_2q_clifford(n_qubits: int, lines: list[str]) -> None:
    gate = pyrand.choice(CLIFFORD_GATES_2Q)
    a, b = pyrand.sample(range(n_qubits), 2)
    lines.append(f"{gate} {a} {b}")


def _emit_noise(n_qubits: int, lines: list[str]) -> None:
    gate, probs = pyrand.choice(NOISE_GATES_1Q)
    p = pyrand.choice(probs)
    q = pyrand.randrange(n_qubits)
    lines.append(f"{gate}({p}) {q}")


def _emit_readout(n_qubits: int, lines: list[str]) -> None:
    p = pyrand.choice([0.001, 0.01, 0.1])
    q = pyrand.randrange(n_qubits)
    lines.append(f"M({p}) {q}")


def gen_program(n_qubits: int, n_instructions: int, regime: str, seed: int) -> str:
    pyrand.seed(seed)
    qs = " ".join(str(i) for i in range(n_qubits))
    lines = [f"R {qs}"]
    for _ in range(n_instructions):
        roll = pyrand.random()
        if regime == "clifford-only":
            if roll < 0.6 or n_qubits < 2:
                _emit_1q_clifford(n_qubits, lines)
            else:
                _emit_2q_clifford(n_qubits, lines)
        elif regime == "+noise":
            if roll < 0.30:
                _emit_noise(n_qubits, lines)
            elif roll < 0.7 or n_qubits < 2:
                _emit_1q_clifford(n_qubits, lines)
            else:
                _emit_2q_clifford(n_qubits, lines)
        elif regime == "+readout":
            if roll < 0.10:
                _emit_readout(n_qubits, lines)
            elif roll < 0.30:
                _emit_noise(n_qubits, lines)
            elif roll < 0.7 or n_qubits < 2:
                _emit_1q_clifford(n_qubits, lines)
            else:
                _emit_2q_clifford(n_qubits, lines)
        else:
            raise ValueError(f"unknown regime {regime}")
    lines.append(f"M {qs}")
    return "\n".join(lines) + "\n"


def run() -> int:
    paths = core.CorpusPaths.default()
    failures: list[str] = []
    written = 0
    for regime in REGIMES:
        for n in QUBIT_COUNTS:
            for length in LENGTHS:
                for seed in SEEDS:
                    src = gen_program(n, length, regime, seed)
                    regime_slug = regime.replace("+", "plus_").replace("-", "_")
                    name = f"{regime_slug}_n{n}_len{length}_s{seed}"
                    meta = core.FixtureMeta(
                        name=name,
                        category="generated/random",
                        source=src,
                        test_num_shots=128,
                    )
                    try:
                        core.write_distribution_fixture(meta, paths)
                        written += 1
                    except Exception as e:
                        failures.append(f"{name}: {e}")
    print(f"regen-stim random: wrote {written} fixtures")
    if failures:
        print("regen-stim random: failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0
