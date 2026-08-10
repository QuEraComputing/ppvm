# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""The Qiskit `pauli-prop` side of the cross-library Pauli-propagation benchmark.

See `README.md` for the shared circuit definitions, the parameter contract, and
the CSV schema — every runner reads the same environment variables and prints
the same columns.

`pauli-prop` truncates on `atol` *and* on a mandatory `max_terms` cap, which the
other three engines do not have. `MAX_TERMS` therefore defaults high enough to
be non-binding, and the runner asserts it never bound: if the propagated support
ever reaches the cap the row is not comparable and the run fails rather than
quietly reporting a differently-truncated number.

    MODEL=tfim QUBITS=8,16,24 STEPS=10 DT=0.1 JCOUP=1.0 HFIELD=1.0 ATOL=1e-6 \
      uv run --no-project --with pauli-prop python3 xbench_qiskit.py
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np
from pauli_prop import propagate_through_circuit
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp

MODEL = os.environ.get("MODEL", "tfim")
STEPS = int(os.environ.get("STEPS", "10"))
DT = float(os.environ.get("DT", "0.1"))
JCOUP = float(os.environ.get("JCOUP", "1.0"))
HFIELD = float(os.environ.get("HFIELD", "1.0"))
ATOL = float(os.environ.get("ATOL", "1e-6"))
ITERS = int(os.environ.get("ITERS", "3"))
MAX_TERMS = int(os.environ.get("MAX_TERMS", str(1 << 22)))

THETA_BOND = 2 * JCOUP * DT
THETA_SITE = 2 * HFIELD * DT


def seed_operator(n: int) -> SparsePauliOp:
    """`Σ_i Z_i` for TFIM, `Z_0` for Heisenberg.

    Qiskit Pauli labels are little-endian: the rightmost character is qubit 0.
    """
    if MODEL == "tfim":
        labels = ["I" * (n - 1 - i) + "Z" + "I" * i for i in range(n)]
    else:
        labels = ["I" * (n - 1) + "Z"]
    return SparsePauliOp(labels, coeffs=np.ones(len(labels)))


def gate_list(n: int) -> list[tuple[str, int]]:
    """`STEPS` first-order Trotter steps in *application* order, per the spec."""
    gates: list[tuple[str, int]] = []
    for _ in range(STEPS):
        if MODEL == "tfim":
            gates += [("rx", i) for i in range(n)]
            gates += [("rzz", i) for i in range(n - 1)]
        else:
            for i in range(n - 1):
                gates += [("rxx", i), ("ryy", i), ("rzz", i)]
            gates += [("rz", i) for i in range(n)]
    return gates


def build_circuit(n: int) -> QuantumCircuit:
    """The spec's gate sequence, as a circuit `pauli-prop` will apply in order.

    A Qiskit circuit whose instructions are appended `g_1 … g_k` denotes the
    unitary `U = g_k ⋯ g_1`, and in the Heisenberg frame `pauli-prop` computes
    `U† O U` by conjugating **from the end of the instruction list backwards**.
    So the gate `pauli-prop` applies first is the one appended last, and the
    spec's application order is obtained by appending in reverse.

    This is not cosmetic. Appending in forward order propagates a genuinely
    different operator: at `n=4, steps=3, atol=1e-14` it yields 108 terms
    against TFIM's reference 124 and 61 against Heisenberg's 64, with
    coefficients off by up to 0.1. Reversed, all four engines agree
    term-for-term to the last bit — which is what `run_xbench.py --validate`
    checks before it will report a timing.
    """
    qc = QuantumCircuit(n)
    for name, i in reversed(gate_list(n)):
        if name == "rx":
            qc.rx(THETA_SITE, i)
        elif name == "rz":
            qc.rz(THETA_SITE, i)
        elif name == "rxx":
            qc.rxx(THETA_BOND, i, i + 1)
        elif name == "ryy":
            qc.ryy(THETA_BOND, i, i + 1)
        elif name == "rzz":
            qc.rzz(THETA_BOND, i, i + 1)
        else:  # pragma: no cover - the list above is closed
            raise ValueError(f"unknown gate {name}")
    return qc


def support(op: SparsePauliOp) -> dict[str, float]:
    """`{word: coefficient}` with site 0 leftmost, duplicate Paulis summed.

    `propagate_through_circuit` returns an operator that may list the same Pauli
    more than once — `len(op)` counts rows, not distinct terms — so the support
    size comparable with the other three engines is the size of this mapping.
    Qiskit labels are little-endian, hence the reversal.
    """
    terms: dict[str, float] = {}
    for pauli, coeff in zip(op.paulis, np.asarray(op.coeffs)):
        word = str(pauli)[::-1]
        terms[word] = terms.get(word, 0.0) + float(np.real(coeff))
    return terms


def readout(op: SparsePauliOp, n: int) -> float:
    """`⟨0…0|O|0…0⟩` for TFIM; the `Z_0` autocorrelator for Heisenberg."""
    paulis = op.paulis
    coeffs = np.asarray(op.coeffs)
    if MODEL == "tfim":
        # ⟨0|Z|0⟩ = 1 and ⟨0|X|0⟩ = ⟨0|Y|0⟩ = 0, so only the X-free terms survive.
        diagonal = ~paulis.x.any(axis=1)
        return float(np.real(coeffs[diagonal].sum()))
    z0_x = np.zeros(n, dtype=bool)
    z0_z = np.zeros(n, dtype=bool)
    z0_z[0] = True
    hit = (paulis.x == z0_x).all(axis=1) & (paulis.z == z0_z).all(axis=1)
    return float(np.real(coeffs[hit].sum()))


def dump(n: int) -> None:
    """Print the propagated support as `word coefficient`, largest first.

    Site 0 leftmost, matching the Rust runner's `DUMP=1` output so the driver
    can diff the two term-for-term.
    """
    out, _bias = propagate_through_circuit(
        seed_operator(n), build_circuit(n), max_terms=MAX_TERMS, atol=ATOL, frame="h"
    )
    terms = support(out)
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
        f"pauli-prop {MODEL}: steps={STEPS} dt={DT} J={JCOUP} h={HFIELD} "
        f"atol={ATOL} iters={ITERS} max_terms={MAX_TERMS}",
        file=sys.stderr,
    )
    print("model,library,qubits,steps,dt,atol,time_s,terms,observable")
    for n in qubits:
        # Circuit construction is not propagation; built once, outside the
        # timed region, as in every other runner.
        circuit = build_circuit(n)
        seed = seed_operator(n)

        best = float("inf")
        terms = 0
        obs = float("nan")
        for _ in range(ITERS):
            t0 = time.perf_counter()
            out, _bias = propagate_through_circuit(
                seed, circuit, max_terms=MAX_TERMS, atol=ATOL, frame="h"
            )
            best = min(best, time.perf_counter() - t0)
            # Distinct Paulis, not the returned operator's row count.
            terms = len(support(out))
            obs = readout(out, n)
        if terms >= MAX_TERMS:
            raise SystemExit(
                f"n={n}: support hit the max_terms cap ({terms} >= {MAX_TERMS}); "
                "raise MAX_TERMS — this row would be truncated differently from "
                "the other engines and is not comparable"
            )
        print(f"{MODEL},pauli-prop,{n},{STEPS},{DT},{ATOL},{best:.7g},{terms},{obs!r}")
        print(f"  n={n}  {best:.4f}s  {terms} terms  obs={obs}", file=sys.stderr)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
