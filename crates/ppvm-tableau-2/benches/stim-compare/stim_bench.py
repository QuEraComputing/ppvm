"""Time official Stim's TableauSimulator on a circuit, the way ppvm's
`run_stim` runs it: a fresh simulator per shot, one `do_circuit` per shot.

Usage: `stim_bench.py <circuit.stim> <shots>`. Prints `<x> µs/shot`.
The per-shot Python call adds about a microsecond; it matters only for the
smallest circuits.
"""

import pathlib
import sys
import time

import stim


def main() -> None:
    circuit = stim.Circuit(pathlib.Path(sys.argv[1]).read_text())
    shots = int(sys.argv[2])
    n = circuit.num_qubits
    for _ in range(2):  # warm-up
        stim.TableauSimulator().do_circuit(circuit)
    start = time.perf_counter()
    for _ in range(shots):
        sim = stim.TableauSimulator()
        sim.set_num_qubits(n)
        sim.do_circuit(circuit)
    elapsed = time.perf_counter() - start
    print(f"{elapsed / shots * 1e6:.2f} µs/shot (stim {stim.__version__}, {n} qubits)")


if __name__ == "__main__":
    main()
