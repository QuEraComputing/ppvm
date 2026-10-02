"""Per-gate cost: `REPEAT 200` of a gate layer on n qubits, ppvm vs official Stim.

`CX` is a brickwork (pairs (0,1),(2,3),... then (1,2),(3,4),...), `H` hits
every qubit. Prints ns per gate. Usage: `gates.py [shots]`.
"""

import pathlib
import re
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
PYTHON = "/Users/david/git/ppvm/ppvm-python/.venv/bin/python"
SIZES = [32, 64, 128, 274, 512, 1024, 2048]
REPS = 200


def layers(gate: str, n: int) -> tuple[str, int]:
    """(circuit, gates per shot)."""
    if gate == "CX":
        even = list(range(0, n - n % 2))
        odd = list(range(1, n - 1 - (n - 1) % 2 + 1))[: (n - 1) // 2 * 2]
        body = f"CX {' '.join(map(str, even))}\nCX {' '.join(map(str, odd))}\n"
        count = (len(even) + len(odd)) // 2
    else:
        body = f"H {' '.join(map(str, range(n)))}\n"
        count = n
    return f"REPEAT {REPS} {{\n{body}}}\n", REPS * count


def execute_us(out: str) -> float:
    m = re.search(r"execute median\s+([\d.]+)(ns|µs|ms|s)\b", out)
    return float(m.group(1)) * {"ns": 1e-3, "µs": 1.0, "ms": 1e3, "s": 1e6}[m.group(2)]


def main() -> None:
    shots = sys.argv[1] if len(sys.argv) > 1 else "50"
    out_dir = HERE / "circuits/gates"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"{'circuit':10s} {'gates':>7s} {'ppvm ns/gate':>13s} {'stim ns/gate':>13s} {'ratio':>6s}")
    for gate in ["CX", "H"]:
        for n in SIZES:
            text, count = layers(gate, n)
            path = out_dir / f"{gate}_n{n}.stim"
            path.write_text(text)
            ppvm = execute_us(subprocess.run(
                [str(HERE / "target/release/stim-compare"), "file", str(path), shots, "ppvm", str(n)],
                capture_output=True, text=True, check=True).stdout)
            stim = execute_us(subprocess.run(
                [PYTHON, str(HERE / "stim_bench.py"), "file", str(path), shots],
                capture_output=True, text=True, check=True).stdout)
            print(f"{gate}_n{n:<6d} {count:7d} {1e3 * ppvm / count:13.2f} {1e3 * stim / count:13.2f} {ppvm / stim:6.2f}")


if __name__ == "__main__":
    main()
