"""
Heisenberg-picture GHZ example.

Compute <ZZ> for the GHZ circuit ``H(0); CNOT(0, 1)`` by propagating ``ZZ``
backwards. Because ppvm runs Pauli propagation in the Heisenberg picture,
gates are applied in *reverse* order of the textbook circuit.
"""

from ppvm import PauliSum

state = PauliSum.new(n_qubits=2, terms=["ZZ"])

# Circuit is H(0); CNOT(0, 1) — propagate backwards.
state.cnot(0, 1)
state.h(0)

print(state)  # → 1.000 * IZ
print(state.overlap_with_zero())  # → 1.0
