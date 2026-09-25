# SPDX-FileCopyrightText: 2026 The PPVM Authors
# SPDX-License-Identifier: Apache-2.0
"""Logical Pauli reference and fixtures, independent of physical word packing."""

from __future__ import annotations

import random

OPERATIONS = ("comm", "mul", "h", "s", "cnot", "circuit")
MASK64 = (1 << 64) - 1


def checksum(values):
    result = 14695981039346656037
    for value in values:
        result = ((result ^ value) * 1099511628211) & MASK64
    return f"{result:016x}"


def gate(columns, signs, name, q, t=0):
    """Conjugate all rows at once using Python's arbitrary precision integers."""
    x, z = columns[q]
    if name == "h":
        signs ^= x & z
        columns[q] = (z, x)
    elif name == "s":
        signs ^= x & z
        columns[q] = (x, z ^ x)
    elif name == "cnot":
        xt, zt = columns[t]
        signs ^= x & zt & ~(xt ^ z)
        columns[q] = (x, z ^ zt)
        columns[t] = (xt ^ x, zt)
    else:
        raise ValueError(name)
    return signs


def seed_gates(n):
    """A deterministic all-to-all Clifford circuit, shared with the public API."""
    rng = random.Random(0x204 + n)
    for _ in range(24):
        for q in range(n):
            if rng.randrange(2):
                yield ("h", q)
            if rng.randrange(2):
                yield ("s", q)
        order = list(range(n))
        rng.shuffle(order)
        for a, b in zip(order, order[1:] + order[:1]):
            yield ("cnot", a, b)


def fixture(n, directory):
    """Write a valid 2n-row frame and the exact gates that prepare it."""
    columns = [(1 << q, 1 << (n + q)) for q in range(n)]
    signs = 0
    gates = list(seed_gates(n))
    for instruction in gates:
        signs = gate(columns, signs, *instruction)
    rows = []
    path = directory / f"n{n:05}.txt"
    with path.open("w") as output:
        output.write(f"{n}\n")
        for r in range(2 * n):
            xs = "".join(str((x >> r) & 1) for x, _ in columns)
            zs = "".join(str((z >> r) & 1) for _, z in columns)
            phase = 2 * ((signs >> r) & 1)
            output.write(f"{phase} {xs} {zs}\n")
            rows.append((phase, int(xs[::-1], 2), int(zs[::-1], 2)))
    path.with_suffix(".gates").write_text(
        "".join(" ".join(map(str, instruction)) + "\n" for instruction in gates)
    )
    return rows


def arbitrary_fixture(n, directory):
    """Non-frame inputs cover odd phases and both symplectic parities."""
    rng = random.Random(n + 0xF00D)
    rows = [
        (rng.randrange(4), rng.getrandbits(n), rng.getrandbits(n)) for _ in range(2 * n)
    ]
    with (directory / f"n{n:05}.txt").open("w") as output:
        output.write(f"{n}\n")
        for p, x, z in rows:
            xs, zs = f"{x:0{n}b}"[::-1], f"{z:0{n}b}"[::-1]
            output.write(f"{p} {xs} {zs}\n")
    return rows


def multiply(left, right):
    """Use the X^x Z^z phase convention to derive a Hermitian-Pauli product."""
    lp, lx, lz = left
    rp, rx, rz = right
    x, z = lx ^ rx, lz ^ rz
    # tensor(Y) = i^(number of Y sites) X^x Z^z. Moving right X
    # past left Z adds a minus sign at every common site.
    phase = (
        lp
        + rp
        + (lx & lz).bit_count()
        + (rx & rz).bit_count()
        + 2 * (lz & rx).bit_count()
        - (x & z).bit_count()
    ) % 4
    return phase, x, z


def apply(rows, n, operation):
    rows = rows.copy()
    if operation == "comm":
        return checksum(
            (
                (rows[r][1] & rows[(r + 1) % len(rows)][2]).bit_count()
                + (rows[r][2] & rows[(r + 1) % len(rows)][1]).bit_count()
            )
            % 2
            for r in range(len(rows))
        )
    if operation == "mul":
        for r in range(len(rows)):
            rows[r] = multiply(rows[r], rows[(r + 1) % len(rows)])
    elif operation in ("h", "s", "cnot", "circuit"):
        # This scalar row representation differs from the fixture builder's
        # arbitrary-precision columns and the workers' packed word kernels.
        for q in range(n):
            instructions = (
                ("h", "s", "cnot") if operation == "circuit" else (operation,)
            )
            for name in instructions:
                for r, (p, x, z) in enumerate(rows):
                    a, b = (x >> q) & 1, (z >> q) & 1
                    if name == "h":
                        p = (p + 2 * a * b) % 4
                        x ^= (a ^ b) << q
                        z ^= (a ^ b) << q
                    elif name == "s":
                        p = (p + 2 * a * b) % 4
                        z ^= a << q
                    else:
                        t = (q + 1) % n
                        c, d = (x >> t) & 1, (z >> t) & 1
                        p = (p + 2 * a * d * (1 ^ c ^ b)) % 4
                        x ^= a << t
                        z ^= d << q
                    rows[r] = p, x, z
    elif operation != "transpose_roundtrip":
        raise ValueError(operation)
    return checksum(
        value
        for p, x, z in rows
        for value in (p, *(((x >> q) & 1) | (((z >> q) & 1) << 1) for q in range(n)))
    )


def validate_algebra():
    """The complex 2x2 Pauli matrices define all one-site signed products."""
    import itertools

    matrices = (
        ((1, 0), (0, 1)),
        ((0, 1), (1, 0)),
        ((1, 0), (0, -1)),
        ((0, -1j), (1j, 0)),
    )
    for a, b in itertools.product(range(4), repeat=2):
        p, x, z = multiply((0, a & 1, a >> 1), (0, b & 1, b >> 1))
        actual = tuple(
            tuple(
                sum(matrices[a][i][k] * matrices[b][k][j] for k in range(2))
                for j in range(2)
            )
            for i in range(2)
        )
        expected = tuple(
            tuple(1j**p * v for v in row) for row in matrices[x | (z << 1)]
        )
        assert actual == expected, (a, b, actual, expected)


if __name__ == "__main__":
    validate_algebra()
    print("All 16 one-site Pauli products match their complex matrices.")
