"""Microbenchmark circuits for the inverse tableau, at the d=7 surface code's
size (118 qubits). Run with the `run_stim` example; per-operation cost is
(time per shot - base) / count.

- base: H on all qubits (per-shot overhead).
- cx: base plus 4096 random CX, in instructions of 16 disjoint pairs.
- mr_batched: 10 rounds of H on all, 128 random CX, one MR over all qubits
  (random collapses with fan-outs).
- mr_split: the same with one MR per instruction (no per-instruction batching).

Regenerate with `python3 gen.py` in this directory (seeded, deterministic).
"""
import random
n = 118
rng = random.Random(1)
h_all = "H " + " ".join(map(str, range(n)))
def cx_layer(pairs):
    lines = []
    for _ in range(pairs // 16):
        qs = rng.sample(range(n), 32)
        lines.append("CX " + " ".join(map(str, qs)))
    return lines
def write(name, lines):
    open(f"{name}.stim", "w").write("\n".join(lines) + "\n")
write("base", [h_all])
write("cx", [h_all] + cx_layer(4096))
rounds = []
for _ in range(10):
    rounds += [h_all] + cx_layer(128)
    rounds.append("MR " + " ".join(map(str, range(n))))
write("mr_batched", rounds)
write("mr_split", [l if not l.startswith("MR") else "\n".join(f"MR {q}" for q in range(n)) for l in rounds])
