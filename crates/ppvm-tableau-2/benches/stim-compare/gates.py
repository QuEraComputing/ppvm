"""Per-gate cost over size: a brickwork `CX` layer and a full `H` layer, on
ppvm-tableau-2's `run_stim` and on official Stim. Each circuit is timed at two
repeat counts; the difference removes per-shot overhead (simulator
construction, the Python call).

    uv run --no-project --with stim python gates.py [--extra NAME=PATH]...
"""

import argparse
import statistics

import common

SIZES = [32, 64, 128, 274, 512, 1024, 2048]
REPEATS = (50, 150)


def layer(gate: str, n: int) -> tuple[str, int]:
    """(one layer, gates in it)."""
    if gate == "CX":
        pairs = [(q, q + 1) for q in range(0, n - 1, 2)] + [
            (q, q + 1) for q in range(1, n - 1, 2)
        ]
        even, odd = pairs[: (n // 2)], pairs[(n // 2) :]
        lines = [
            f"CX {' '.join(f'{a} {b}' for a, b in even)}",
            f"CX {' '.join(f'{a} {b}' for a, b in odd)}",
        ]
        return "\n".join(lines), len(pairs)
    return f"H {' '.join(map(str, range(n)))}", n


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extra", action="append", default=[])
    parser.add_argument("--rounds", type=int, default=3)
    args = parser.parse_args()
    common.build_ppvm()
    tools = common.tools(args.extra)
    (common.OUT / "gates").mkdir(parents=True, exist_ok=True)
    for gate in ("CX", "H"):
        for n in SIZES:
            body, count = layer(gate, n)
            paths = []
            for reps in REPEATS:
                path = common.OUT / "gates" / f"{gate}_{n}_{reps}.stim"
                path.write_text(f"REPEAT {reps} {{\n{body}\n}}\n")
                paths.append(path)
            cells = []
            for name, cmd in tools.items():
                shots = common.calibrated_shots(cmd, paths[1], 0.2)
                deltas = []
                for _ in range(args.rounds):
                    lo, hi = (common.time_per_shot(cmd, p, shots) for p in paths)
                    deltas.append((hi - lo) * 1e3 / ((REPEATS[1] - REPEATS[0]) * count))
                cells.append(f"{name} {statistics.median(deltas):6.2f}")
            print(f"{gate:2} n={n:<5} ns/gate: " + "  ".join(cells), flush=True)


if __name__ == "__main__":
    main()
