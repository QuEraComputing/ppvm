# ppvm

**ppvm** is a fast Heisenberg-picture Pauli propagation simulator.

Instead of tracking a full quantum state vector, ppvm represents the quantum
state as a sum of Pauli operators with real coefficients — a Pauli Transfer
Matrix (PTM) representation. Gates are applied by conjugating each term in
the sum, and noise channels are applied as super-operators. This makes it
possible to simulate large and deep quantum circuits, including noisy ones,
at a fraction of the cost of full statevector simulation.

## Key features

- Clifford gates, arbitrary single- and two-qubit rotations, and common noise
  channels out of the box
- Automatic truncation: terms below a coefficient threshold or above a Pauli
  weight cutoff are dropped, controlling the approximation/performance trade-off
- Optional loss channel support via `LossyPauliSum`
- Python bindings with a clean, gate-by-gate API backed by a high-performance
  Rust core

## Quick example

```python
from ppvm import PauliSum

# Initialise a two-qubit state |00><00| = (I + Z) ⊗ (I + Z) / 4
# represented as the operator ZZ (up to normalisation)
state = PauliSum(initial_terms=["ZZ"])

# Apply a CNOT followed by a Hadamard — creates a Bell state
state.cnot(0, 1)
state.h(0)

# The state is now represented as IZ
print(state)  # 1.000 * IZ
```

## Getting started

- [Installation](install.md)
- [Examples](examples/index.md)
- [API Reference](api/index.md)
