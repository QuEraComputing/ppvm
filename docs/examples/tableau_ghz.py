"""
Generalized-tableau GHZ example.

Prepare a Bell/GHZ-on-2 state via ``H(0); CNOT(0, 1)`` and measure both
qubits. The outcomes are perfectly correlated, as expected.
"""

from ppvm import GeneralizedTableau

tab = GeneralizedTableau(n_qubits=2, seed=0)
tab.h(0)
tab.cnot(0, 1)

r0 = tab.measure(0)
r1 = tab.measure(1)

print(f"qubit 0: {r0}, qubit 1: {r1}")  # → qubit 0: 1, qubit 1: 1
print("correlated:", r0 == r1)  # → correlated: True
