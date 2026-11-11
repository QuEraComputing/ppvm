"""Stim-generated Clifford memory circuits, as both simulators run them: noise
removed, and the two instructions ppvm's harness lacks rewritten.

- `C_XYZ` is `H` then `SQRT_X_DAG` (checked with `stim.Tableau`).
- `MX` and `MY` (only the color code's final data measurements) are `H`, `M`
  and `S_DAG`, `H`, `M`: they map X or Y to Z first; the state afterwards is
  not used.
"""

import re

import stim

# (family, Stim generator task, distances); rounds = distance.
FAMILIES = [
    ("surface", "surface_code:rotated_memory_z", [3, 5, 7, 9, 11, 15, 19, 23, 27, 31]),
    ("repetition", "repetition_code:memory", [3, 9, 25, 75, 225, 675]),
    ("color", "color_code:memory_xyz", [3, 5, 7, 9, 13, 17, 21, 25, 31, 37, 43]),
]


def generate(task: str, distance: int) -> str:
    """The circuit text for `task` at `distance`, noise-free and rewritten."""
    noise = 0.001  # only so the generator emits the full structure; removed below
    circuit = stim.Circuit.generated(
        task,
        distance=distance,
        rounds=distance,
        after_clifford_depolarization=noise,
        before_round_data_depolarization=noise,
        before_measure_flip_probability=noise,
        after_reset_flip_probability=noise,
    ).without_noise()
    return rewrite(str(circuit))


def rewrite(text: str) -> str:
    def c_xyz(m: re.Match) -> str:
        targets = m.group(1)
        return f"H {targets}\nSQRT_X_DAG {targets}"

    def mx(m: re.Match) -> str:
        targets = m.group(1)
        return f"H {targets}\nM {targets}"

    def my(m: re.Match) -> str:
        targets = m.group(1)
        return f"S_DAG {targets}\nH {targets}\nM {targets}"

    text = re.sub(r"^\s*C_XYZ ([\d ]+)$", c_xyz, text, flags=re.MULTILINE)
    text = re.sub(r"^\s*MX ([\d ]+)$", mx, text, flags=re.MULTILINE)
    return re.sub(r"^\s*MY ([\d ]+)$", my, text, flags=re.MULTILINE)


def num_qubits(text: str) -> int:
    return stim.Circuit(text).num_qubits
