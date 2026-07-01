"""Encode a message into a computational-basis state with GeneralizedTableau.

Each bit of the ASCII encoding of "Hello, world!" gets one qubit. Starting
from the all-zero state, we flip every qubit whose target bit is 1 with an X
gate, then measure all qubits back out and decode them into the string.
"""

from ppvm import GeneralizedTableau

MESSAGE = "Hello, world!"

# ASCII bytes -> bit list, MSB first within each byte.
data = MESSAGE.encode("ascii")
bits = [(byte >> (7 - i)) & 1 for byte in data for i in range(8)]

# One qubit per bit; flip every qubit whose target bit is 1.
tab = GeneralizedTableau(n_qubits=len(bits))
for q, bit in enumerate(bits):
    if bit:
        tab.x(q)

# Measure all qubits and decode the outcomes back into a string.
measured = [int(tab.measure(q)) for q in range(len(bits))]
decoded = bytes(
    int("".join(str(b) for b in measured[i : i + 8]), 2)
    for i in range(0, len(measured), 8)
).decode("ascii")

print(f"qubits:  {len(bits)}")
print(f"decoded: {decoded!r}")
assert decoded == MESSAGE, f"round-trip failed: {decoded!r}"
print("round-trip OK")
