"""Validate Heisenberg-picture conventions against a numpy density-matrix sim.

`PauliSum` propagates an *observable* backwards through a circuit. Two
conventions follow from that and are easy to get wrong:

1. **Noise ordering.** A physical circuit applies a gate and *then* the
   noise that models that gate's imperfection (gate → noise). Because
   `PauliSum` runs backwards, the whole sequence must be reversed, so the
   noise call must come *before* the gate call in `PauliSum` code
   (noise → gate). Writing gate-then-noise in `PauliSum` code silently
   simulates noise-then-gate, which is a different circuit.

2. **Rotation angle sign.** `state.rx(theta)` already applies the adjoint
   conjugation ``U† O U`` with ``U = exp(-i theta/2 X)``. So backward
   propagation only requires *reversing the gate order* — the angles keep
   their original sign. Negating angles "because it runs backwards" is the
   classic mistake and gives the inverse circuit.

We check both by comparing `PauliSum` against a tiny first-principles
density-matrix simulator built from Kraus operators. The forward sim uses
standard textbook gate unitaries and noise Kraus maps; the `PauliSum` side
reverses the operation list. The identity being exercised is

    Tr[O · C(rho0)] == <0...0| C†(O) |0...0>

where ``C`` is the forward channel and ``C†(O)`` is what `PauliSum` stores.
"""

import functools

import numpy as np
import pytest

from ppvm import PauliSum

# --- single-qubit operators -------------------------------------------------

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.array([[1, 0], [0, 1j]], dtype=complex)
PAULIS = [I2, X, Y, Z]


def _embed1(op, q, n):
    """Place a 2x2 operator on qubit q (qubit 0 = leftmost kron factor)."""
    mats = [op if i == q else I2 for i in range(n)]
    return functools.reduce(np.kron, mats)  # ty: ignore[invalid-argument-type]


def _embed2(g4, q0, q1, n):
    """Embed a 4x4 two-qubit gate on (q0, q1) via its Pauli decomposition.

    Works for arbitrary (even non-adjacent) qubit positions.
    """
    out = np.zeros((2**n, 2**n), dtype=complex)
    for _, pa in enumerate(PAULIS):
        for _, pb in enumerate(PAULIS):
            coeff = np.trace(np.kron(pa, pb).conj().T @ g4) / 4
            if abs(coeff) > 1e-15:
                out += coeff * (_embed1(pa, q0, n) @ _embed1(pb, q1, n))
    return out


def _rot(pauli, theta):
    """exp(-i theta/2 P) for an involutory Pauli operator P (P**2 = I)."""
    return np.cos(theta / 2) * np.eye(pauli.shape[0]) - 1j * np.sin(theta / 2) * pauli


# --- forward (Schrodinger) circuit on a density matrix ----------------------


class DensityMatrixSim:
    """Minimal exact density-matrix simulator (forward / Schrodinger picture)."""

    def __init__(self, n):
        self.n = n
        dim = 2**n
        rho = np.zeros((dim, dim), dtype=complex)
        rho[0, 0] = 1.0  # |0...0><0...0|
        self.rho = rho

    def _unitary(self, u_full):
        self.rho = u_full @ self.rho @ u_full.conj().T

    def _kraus(self, ks_full):
        self.rho = sum(k @ self.rho @ k.conj().T for k in ks_full)

    # gates (standard textbook unitaries) -----------------------------------
    def h(self, q):
        self._unitary(_embed1(H, q, self.n))

    def s(self, q):
        self._unitary(_embed1(S, q, self.n))

    def x(self, q):
        self._unitary(_embed1(X, q, self.n))

    def rx(self, q, theta):
        self._unitary(_embed1(_rot(X, theta), q, self.n))

    def ry(self, q, theta):
        self._unitary(_embed1(_rot(Y, theta), q, self.n))

    def rz(self, q, theta):
        self._unitary(_embed1(_rot(Z, theta), q, self.n))

    def cnot(self, c, t):
        cnot4 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)
        self._unitary(_embed2(cnot4, c, t, self.n))

    def cz(self, q0, q1):
        cz4 = np.diag([1, 1, 1, -1]).astype(complex)
        self._unitary(_embed2(cz4, q0, q1, self.n))

    def rzz(self, q0, q1, theta):
        self._unitary(_embed2(_rot(np.kron(Z, Z), theta), q0, q1, self.n))

    def rxx(self, q0, q1, theta):
        self._unitary(_embed2(_rot(np.kron(X, X), theta), q0, q1, self.n))

    def ryy(self, q0, q1, theta):
        self._unitary(_embed2(_rot(np.kron(Y, Y), theta), q0, q1, self.n))

    # noise channels (first-principles Kraus operators) ---------------------
    def depolarize1(self, q, p):
        ks = [
            np.sqrt(1 - p) * I2,
            np.sqrt(p / 3) * X,
            np.sqrt(p / 3) * Y,
            np.sqrt(p / 3) * Z,
        ]
        self._kraus([_embed1(k, q, self.n) for k in ks])

    def pauli_error(self, q, p):
        px, py, pz = p
        ks = [
            np.sqrt(1 - px - py - pz) * I2,
            np.sqrt(px) * X,
            np.sqrt(py) * Y,
            np.sqrt(pz) * Z,
        ]
        self._kraus([_embed1(k, q, self.n) for k in ks])

    def amplitude_damping(self, q, gamma):
        k0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
        k1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
        self._kraus([_embed1(k0, q, self.n), _embed1(k1, q, self.n)])

    def expectation_sum_z(self):
        """<sum_i Z_i> = Tr[(sum_i Z_i) rho]."""
        obs = sum(_embed1(Z, i, self.n) for i in range(self.n))
        val = np.trace(obs @ self.rho)
        assert abs(val.imag) < 1e-12
        return val.real


# --- driving PauliSum from the same operation list --------------------------


def _apply_ppvm(ps, ops):
    """Apply an ordered op list to a PauliSum (`ops` already in call order)."""
    for name, args in ops:
        getattr(ps, name)(*args)


def _ppvm_sum_z(n, forward_ops):
    """Reverse the forward op list and propagate <sum_i Z_i> backwards."""
    ps = PauliSum.new(n, [f"Z{i}" for i in range(n)])
    _apply_ppvm(ps, list(reversed(forward_ops)))
    return ps.overlap_with_zero()


def _dm_sum_z(n, forward_ops):
    sim = DensityMatrixSim(n)
    for name, args in forward_ops:
        getattr(sim, name)(*args)
    return sim.expectation_sum_z()


# A forward circuit (physical order: each gate is followed by its noise).
# Mixes Clifford gates, rotations, an asymmetric Pauli error, depolarizing,
# and amplitude damping so that operation ordering genuinely matters.
FORWARD_CIRCUIT = [
    ("h", (0,)),
    ("cnot", (0, 1)),
    ("rx", (1, 0.7)),
    ("pauli_error", (1, [0.05, 0.0, 0.15])),  # asymmetric -> order matters
    ("rzz", (1, 2, 0.4)),
    ("depolarize1", (2, 0.1)),
    ("ry", (0, 1.1)),
    ("amplitude_damping", (0, 0.2)),  # non-unital -> order matters
]


def test_reversed_pauli_sum_matches_density_matrix():
    """Reversing the full op list reproduces the exact density-matrix result."""
    n = 3
    expected = _dm_sum_z(n, FORWARD_CIRCUIT)
    got = _ppvm_sum_z(n, FORWARD_CIRCUIT)
    assert got == pytest.approx(expected, abs=1e-9)


def test_gate_then_noise_in_pauli_sum_is_wrong():
    """Gate-then-noise in PauliSum code simulates the wrong circuit.

    Physical circuit: rx(theta) on |0>, then a pure-dephasing channel.
    Correct Heisenberg code calls noise *before* the gate; the buggy version
    calls gate *before* noise. Only the correct one matches the DM result, and
    the buggy one differs by a non-trivial amount (so this is a real trap, not
    a rounding artefact).
    """
    n = 1
    theta = 0.9
    forward = [("rx", (0, theta)), ("pauli_error", (0, [0.0, 0.0, 0.3]))]
    expected = _dm_sum_z(n, forward)

    correct = _ppvm_sum_z(n, forward)  # reversed -> noise then gate
    assert correct == pytest.approx(expected, abs=1e-9)

    buggy = PauliSum.new(n, ["Z0"])
    buggy.rx(0, theta)  # gate first ...
    buggy.pauli_error(0, [0.0, 0.0, 0.3])  # ... then noise (WRONG)
    assert abs(buggy.overlap_with_zero() - expected) > 1e-3


def test_rotation_angles_keep_their_sign_when_reversed():
    """Backward propagation reverses gate order but does NOT negate angles."""
    n = 1
    theta = 1.3
    forward = [("rx", (0, theta))]
    expected = _dm_sum_z(n, forward)

    # Same sign, reversed order (here a single gate) -> correct.
    same_sign = _ppvm_sum_z(n, forward)
    assert same_sign == pytest.approx(expected, abs=1e-9)

    # Negating the angle is the classic mistake -> wrong (gives the inverse).
    negated = PauliSum.new(n, ["Z0"])
    negated.rx(0, -theta)
    assert abs(negated.overlap_with_zero() - expected) > 1e-3
