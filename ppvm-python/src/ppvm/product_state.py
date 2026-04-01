from __future__ import annotations

import warnings
from typing import Sequence

import ppvm_python_native

from .paulisum import PauliSum


class ProductState:
    """A separable initial state ρ₀ = ⊗ᵢ ρ₀⁽ⁱ⁾, encoded as per-qubit Bloch vectors.

    Used to compute expectation values ⟨O(t)⟩ = Tr(ρ₀ O(t)) after Heisenberg-picture
    time evolution.  ρ₀ is never propagated — it is evaluated at output checkpoints only.

    The Bloch vector (bx, by, bz) for qubit i gives:
        Tr(ρ₀⁽ⁱ⁾ I) = 1,  Tr(ρ₀⁽ⁱ⁾ X) = bx,  Tr(ρ₀⁽ⁱ⁾ Y) = by,  Tr(ρ₀⁽ⁱ⁾ Z) = bz.
    Pure states satisfy |b|² = 1; mixed states |b|² < 1.
    """

    def __init__(self, bloch: list[float]) -> None:
        """Low-level constructor.  `bloch` is a flat array [bx₀,by₀,bz₀, bx₁,…]."""
        if len(bloch) % 3 != 0:
            raise ValueError(f"bloch must have length divisible by 3, got {len(bloch)}")
        self._bloch = bloch
        self._n_qubits = len(bloch) // 3

    # ------------------------------------------------------------------ constructors

    @staticmethod
    def all_zero(n_qubits: int) -> "ProductState":
        """All qubits in |0⟩: bz = +1."""
        return ProductState([v for _ in range(n_qubits) for v in (0.0, 0.0, 1.0)])

    @staticmethod
    def all_one(n_qubits: int) -> "ProductState":
        """All qubits in |1⟩: bz = -1."""
        return ProductState([v for _ in range(n_qubits) for v in (0.0, 0.0, -1.0)])

    @staticmethod
    def bitstring(bits: str | Sequence[int]) -> "ProductState":
        """Computational-basis state.

        Args:
            bits: String of '0'/'1' characters or a sequence of 0/1 integers.
                  bits[i] = 0 → |0⟩ (bz=+1);  bits[i] = 1 → |1⟩ (bz=-1).

        Example:
            ProductState.bitstring("0101")   # 4-qubit state |0101⟩
        """
        bloch = []
        for b in bits:
            bit = int(b)
            if bit not in (0, 1):
                raise ValueError(f"bitstring contains invalid character {b!r}")
            bz = 1.0 if bit == 0 else -1.0
            bloch.extend([0.0, 0.0, bz])
        return ProductState(bloch)

    @staticmethod
    def bloch_vectors(
        vectors: Sequence[tuple[float, float, float]],
    ) -> "ProductState":
        """Arbitrary product state via explicit per-qubit Bloch vectors.

        Args:
            vectors: List of (bx, by, bz) tuples, one per qubit.

        Example:
            ProductState.bloch_vectors([(0,0,1), (1,0,0)])  # |0⟩ ⊗ |+⟩
        """
        bloch = []
        for i, (bx, by, bz) in enumerate(vectors):
            norm_sq = bx**2 + by**2 + bz**2
            if norm_sq > 1.0 + 1e-9:
                warnings.warn(
                    f"Bloch vector for qubit {i} has |b|² = {norm_sq:.6g} > 1 "
                    f"(not a valid density matrix). Proceeding anyway.",
                    stacklevel=2,
                )
            bloch.extend([bx, by, bz])
        return ProductState(bloch)

    # ------------------------------------------------------------------ properties

    @property
    def n_qubits(self) -> int:
        return self._n_qubits

    # ------------------------------------------------------------------ expectation

    def expectation(self, observable: PauliSum) -> float:
        """Compute ⟨O⟩ = Tr(ρ₀ O) for a Heisenberg-picture observable O.

        Args:
            observable: A PauliSum representing the (possibly evolved) observable.

        Returns:
            Real-valued expectation value.
        """
        return ppvm_python_native.product_state_expectation(
            observable._interface, self._bloch
        )
