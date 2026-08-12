# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""The monoprop side of the cross-library Pauli-propagation benchmark.

See `README.md` for the shared circuit definitions, the parameter contract, and
the CSV schema — every runner reads the same environment variables and prints the
same columns.

Three things about monoprop need saying, because each of them silently changes
the number this file prints:

**It is parallel by default.** With one MPI rank the engine takes one serial
partition per *physical core* (`resolve_partition_count_`), so an uncapped run on
this 14-core machine is a ~9.5x-CPU, ~3x-wall-faster run than the serial one, and
is not comparable with the other single-threaded engines. `monoprop_NUM_THREADS=1`
and `monoprop_PARTITIONS=off` pin it to one partition, and the runner **measures
its own CPU/wall ratio and fails** if the cap did not take. An environment
variable that is read once into a cached C++ static is exactly the kind of setting
that goes stale without anyone noticing.

**Its angles are not the spec's angles.** `ExpGate` applies `exp(+iθH)`, so the
spec's `exp(-iθ/2·G)` needs `θ_monoprop = -θ_spec/2`. All four other engines take
`θ_spec` directly.

**It reports more rows than it has terms.** The engine retains monomials whose
coefficient has cancelled to exactly zero, and `size()` counts those. At
`n=12` Heisenberg that is 3.06M tracked rows against 1.32M terms above the
truncation threshold. The `terms` column is the above-threshold support, which is
what the other four engines mean by it; the tracked size goes to stderr.

    MODEL=tfim QUBITS=8,16,24 STEPS=10 DT=0.1 JCOUP=1.0 HFIELD=1.0 ATOL=1e-6 \
      uv run --no-project --with monoprop python3 xbench_monoprop.py
"""

from __future__ import annotations

import os

# Must precede the import: the C++ side reads these once into a cached static
# (`monoprop::config::get()`), so a later assignment would be ignored. The driver
# sets them in the subprocess environment too; this makes a direct invocation of
# this file single-threaded as well, rather than quietly parallel.
os.environ.setdefault("monoprop_NUM_THREADS", "1")
os.environ.setdefault("monoprop_PARTITIONS", "off")

import resource
import sys
import time

from monoprop import (
    Circuit,
    ExpGate,
    Pauli,
    PauliOperator,
    PauliPropagator,
)

MODEL = os.environ.get("MODEL", "tfim")
STEPS = int(os.environ.get("STEPS", "10"))
DT = float(os.environ.get("DT", "0.1"))
JCOUP = float(os.environ.get("JCOUP", "1.0"))
HFIELD = float(os.environ.get("HFIELD", "1.0"))
ATOL = float(os.environ.get("ATOL", "1e-6"))
ITERS = int(os.environ.get("ITERS", "3"))
SEED = int(os.environ.get("SEED", "12345"))
# A serial run has CPU/wall ~1. Anything materially above it means the thread cap
# did not take and the timing is not comparable with the other engines.
CPU_WALL_MAX = float(os.environ.get("CPU_WALL_MAX", "1.5"))

THETA_BOND = 2 * JCOUP * DT
THETA_SITE = 2 * HFIELD * DT


def cpu_seconds() -> float:
    """Total CPU time (user+sys) charged to this process and its children."""
    me = resource.getrusage(resource.RUSAGE_SELF)
    kids = resource.getrusage(resource.RUSAGE_CHILDREN)
    return me.ru_utime + me.ru_stime + kids.ru_utime + kids.ru_stime


def seed_operator(n: int) -> PauliOperator:
    """`Σ_i Z_i` for TFIM, `Z_0` for Heisenberg and scramble."""
    if MODEL == "tfim":
        terms: dict[Pauli | str, float] = {Pauli("Z", (i,)): 1.0 for i in range(n)}
    else:
        terms = {Pauli("Z", (0,)): 1.0}
    return PauliOperator(terms, n)


class SplitMix64:
    """splitmix64, matching `SplitMix64` in `examples/xbench.rs` bit for bit.

    Both runners have to emit the *same* random circuit, and neither language's
    stdlib RNG is specified tightly enough to rely on. This is short enough to
    state in both and is checked by the term-for-term dump diff.
    """

    MASK = (1 << 64) - 1

    def __init__(self, seed: int) -> None:
        self.state = seed & self.MASK

    def next(self) -> int:
        self.state = (self.state + 0x9E3779B97F4A7C15) & self.MASK
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & self.MASK
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & self.MASK
        return z ^ (z >> 31)

    def unit(self) -> float:
        """Uniform in `[0, 1)`, from the top 53 bits."""
        return (self.next() >> 11) / float(1 << 53)


AXES = "XYZ"


def scramble_gates(n: int) -> list[tuple[str, tuple[int, ...], float]]:
    """`STEPS · n` random all-to-all two-qubit Pauli rotations.

    Draw order — pair, offset, axis, axis, angle — must match `scramble_gates`
    in `examples/xbench.rs` exactly, including the offset trick that avoids a
    rejection loop for `b != a` (a loop would consume a variable number of
    draws and desynchronise the two streams).
    """
    if n < 2:
        raise SystemExit("scramble needs at least 2 qubits")
    rng = SplitMix64(SEED)
    theta_max = 2 * JCOUP * DT
    gates = []
    for _ in range(STEPS * n):
        a = rng.next() % n
        b = (a + 1 + rng.next() % (n - 1)) % n
        pa = AXES[rng.next() % 3]
        pb = AXES[rng.next() % 3]
        theta = theta_max * rng.unit()
        # P_a ⊗ P_b is symmetric under swapping both, so normalise to ascending
        # qubit order -- monoprop's `Pauli` expects sorted supports.
        if a > b:
            a, b, pa, pb = b, a, pb, pa
        gates.append((pa + pb, (a, b), theta))
    return gates


def gate_list(n: int) -> list[tuple[str, tuple[int, ...], float]]:
    """The gate sequence in *application* order, per the spec."""
    if MODEL == "scramble":
        return scramble_gates(n)
    gates: list[tuple[str, tuple[int, ...], float]] = []
    for _ in range(STEPS):
        if MODEL == "tfim":
            gates += [("X", (i,), THETA_SITE) for i in range(n)]
            gates += [("ZZ", (i, i + 1), THETA_BOND) for i in range(n - 1)]
        else:
            for i in range(n - 1):
                gates += [
                    ("XX", (i, i + 1), THETA_BOND),
                    ("YY", (i, i + 1), THETA_BOND),
                    ("ZZ", (i, i + 1), THETA_BOND),
                ]
            gates += [("Z", (i,), THETA_SITE) for i in range(n)]
    return gates


def build_circuit(n: int) -> Circuit:
    """The spec's gate sequence as a monoprop `Circuit`.

    Two conversions happen here, and the term-for-term validation in
    `run_xbench.py` is what pinned both of them down:

    * **Angle.** `ExpGate` applies `exp(+iθH)` (its own docstring flags the
      positive sign as the difference from Qiskit's `r<P>`), while the shared
      spec's gate is `exp(-iθ_spec/2 · G)`. Hence `θ = -θ_spec/2`.
    * **Order.** Like Qiskit's `pauli-prop`, the propagator conjugates from the
      end of the gate list backwards, so the spec's application order is
      obtained by appending in reverse.

    Neither is cosmetic. All eight sign/order combinations were checked against
    the reference dump: forward order loses terms outright (109 against TFIM's
    124), and every wrong angle keeps the right support while moving
    coefficients by up to 1.3. Only this one lands within 5e-13.
    """
    gates = list(reversed(gate_list(n)))
    exp_gates = [
        ExpGate(PauliOperator({Pauli(string, qubits): 1.0}, n), index)
        for index, (string, qubits, _theta) in enumerate(gates)
    ]
    parameters = [-0.5 * theta for (_s, _q, theta) in gates]
    return Circuit(gates=exp_gates, parameters=parameters)


def propagate(n: int, circuit: Circuit) -> PauliPropagator:
    """Heisenberg-picture propagation with coefficient-only truncation.

    `cutoff` is monoprop's mandatory bound on retained Pauli *weight*, which the
    other four engines have no analogue of; `cutoff=n` is the whole register, so
    it never binds and `lower_atol` is left as the only truncation — the shared
    `|c| < atol` rule.
    """
    return PauliPropagator.from_circuit(
        circuit, seed_operator(n), cutoff=n, lower_atol=ATOL
    )


def support(mp: PauliPropagator, n: int) -> dict[str, float]:
    """`{word: coefficient}` with site 0 leftmost, exact-zero rows dropped.

    `atol=0.0` asks the engine for everything it holds, and the threshold is
    applied here so the filter is the same one the other engines apply during
    propagation.
    """
    operator = mp.evolved_operator(atol=0.0)
    terms: dict[str, float] = {}
    for pauli, coeff in operator.terms.items():
        value = float(coeff.real if isinstance(coeff, complex) else coeff)
        if abs(value) < ATOL:
            continue
        chars = ["I"] * n
        for qubit, letter in zip(pauli.qubits, pauli.string):
            chars[qubit] = letter
        terms["".join(chars)] = value
    return terms


def readout(terms: dict[str, float], n: int) -> float:
    """`⟨0…0|O|0…0⟩` for TFIM; the `Z_0` autocorrelator for Heisenberg."""
    if MODEL == "tfim":
        # ⟨0|Z|0⟩ = 1 and ⟨0|X|0⟩ = ⟨0|Y|0⟩ = 0, so only the X/Y-free terms survive.
        return sum(c for w, c in terms.items() if "X" not in w and "Y" not in w)
    return terms.get("Z" + "I" * (n - 1), 0.0)


def dump(n: int) -> None:
    """Print the propagated support as `word coefficient`, largest first."""
    terms = support(propagate(n, build_circuit(n)), n)
    print(f"# {len(terms)} terms")
    for word in sorted(terms, key=lambda w: (-abs(terms[w]), w)):
        print(f"{word} {terms[word]:+.12e}")


def main() -> None:
    qubits = [
        int(t) for t in os.environ.get("QUBITS", "8,12,16,20,24,28,32").split(",")
    ]
    if os.environ.get("DUMP"):
        dump(qubits[0])
        return
    print(
        f"monoprop {MODEL}: steps={STEPS} dt={DT} J={JCOUP} h={HFIELD} atol={ATOL} "
        f"iters={ITERS} threads={os.environ['monoprop_NUM_THREADS']} "
        f"partitions={os.environ['monoprop_PARTITIONS']} cores={os.cpu_count()}",
        file=sys.stderr,
    )
    print("model,library,qubits,steps,dt,atol,time_s,terms,observable")
    for n in qubits:
        # Circuit construction is not propagation; built once, outside the timed
        # region, as in every other runner. monoprop's `propagate` does expand
        # the gate list internally each call, but that is 0.0–2.2% of the total
        # at these widths (measured), not enough to move a ratio.
        circuit = build_circuit(n)

        best = float("inf")
        cpu_at_best = 0.0
        propagator = None
        for _ in range(ITERS):
            cpu0, wall0 = cpu_seconds(), time.perf_counter()
            propagator = propagate(n, circuit)
            wall = time.perf_counter() - wall0
            if wall < best:
                best, cpu_at_best = wall, cpu_seconds() - cpu0

        # A parallel run would report a wall time the other engines never had the
        # chance to compete with. Refuse to print it.
        ratio = cpu_at_best / best if best > 0 else 1.0
        if ratio > CPU_WALL_MAX:
            raise SystemExit(
                f"n={n}: monoprop used {ratio:.2f}x CPU per wall second "
                f"({cpu_at_best:.3f}s CPU in {best:.3f}s wall), so the thread cap "
                f"did not take and this timing is not comparable with the other "
                f"single-threaded engines. Check monoprop_NUM_THREADS / "
                f"monoprop_PARTITIONS reach the subprocess."
            )

        assert propagator is not None
        # Materializing the support is O(size) in Python and dwarfs the
        # propagation at the wide end, so it happens once, after timing.
        terms = support(propagator, n)
        obs = readout(terms, n)
        print(
            f"{MODEL},monoprop,{n},{STEPS},{DT},{ATOL},{best:.7g},{len(terms)},{obs!r}"
        )
        print(
            f"  n={n}  {best:.4f}s  {len(terms)} terms  obs={obs}  "
            f"[cpu/wall={ratio:.2f} tracked={propagator.size()}]",
            file=sys.stderr,
        )
        sys.stdout.flush()


if __name__ == "__main__":
    main()
