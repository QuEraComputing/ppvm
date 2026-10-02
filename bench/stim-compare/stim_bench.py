"""Time Stim's TableauSimulator (official PyPI build) on the same surface_d30
variants as the Rust harness, one shot per `do_circuit` call.

Usage: `stim_bench.py [shots] [variant]`, or `stim_bench.py file <path> [shots]`.
"""

import pathlib
import statistics
import sys
import time

import stim

SRC = (
    pathlib.Path(__file__).parent / "../../crates/ppvm-stim/examples/surface_d30.stim"
).read_text()
WARMUP = 2

MEASURE_OPS = {"M", "MR", "R", "DETECTOR", "OBSERVABLE_INCLUDE"}
NOISE_OPS = {"DEPOLARIZE1", "DEPOLARIZE2", "X_ERROR"}

# Same (name, stripped ops, REPEAT count) table as src/main.rs.
VARIANTS = [
    ("full", set(), 29),
    ("full-r1", set(), 1),
    ("no-noise", NOISE_OPS, 29),
    ("no-measure", MEASURE_OPS, 29),
    ("no-measure-r1", MEASURE_OPS, 1),
    ("gates-only", MEASURE_OPS | NOISE_OPS, 29),
]


def strip(src: str, ops: set[str], reps: int) -> str:
    out = []
    for line in src.splitlines():
        name = line.lstrip().replace("(", " ").split(" ")[0]
        if name in ops:
            continue
        out.append(f"REPEAT {reps} {{" if line == "REPEAT 29 {" else line)
    return "\n".join(out) + "\n"


def run_variant(name: str, src: str, shots: int) -> None:
    circuit = stim.Circuit(src)
    n_qubits = max(circuit.num_qubits, 1)
    build, execute, ones = [], [], []
    for seed in range(WARMUP + shots):
        t0 = time.perf_counter()
        sim = stim.TableauSimulator(seed=seed)
        sim.set_num_qubits(n_qubits)
        t1 = time.perf_counter()
        sim.do_circuit(circuit)
        t2 = time.perf_counter()
        if seed >= WARMUP:
            build.append(t1 - t0)
            execute.append(t2 - t1)
            ones.append(sum(sim.current_measurement_record()))
    ms = lambda s: f"{1e6 * s:10.1f}µs"
    print(f"\n[{name}] {n_qubits} qubits, {circuit.num_measurements} measurements")
    print(
        f" stim: construct {ms(statistics.median(build))}  execute median "
        f"{ms(statistics.median(execute))} (min {ms(min(execute))}, max {ms(max(execute))})"
        f"  | mean ones/shot {statistics.mean(ones):.1f}"
    )


def main() -> None:
    if len(sys.argv) > 2 and sys.argv[1] == "file":
        shots = int(sys.argv[3]) if len(sys.argv) > 3 else 20
        print(f"{sys.argv[2]}: stim {stim.__version__}, {shots} shots")
        run_variant("file", pathlib.Path(sys.argv[2]).read_text(), shots)
        return
    shots = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    only = sys.argv[2] if len(sys.argv) > 2 else None
    print(f"surface_d30: stim {stim.__version__}, {shots} shots")
    for name, ops, reps in VARIANTS:
        if only in (None, name):
            run_variant(name, strip(SRC, ops, reps), shots)


if __name__ == "__main__":
    main()
