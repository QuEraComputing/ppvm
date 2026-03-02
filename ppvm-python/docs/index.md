# ppvm

**ppvm** is a fast and easy-to-use Pauli propagation simulator.

Instead of tracking a full quantum state vector, ppvm represents the quantum
state as a sum of Pauli strings with real coefficients.
Gates are applied by conjugating each term in
the sum, and noise channels are applied as super-operators. This makes it
possible to simulate large and deep quantum circuits, including noisy ones,
at a fraction of the cost of full statevector simulation.

## Feature overview:

- Clifford gates, arbitrary single- and two-qubit rotations, and common noise
  channels.
- Automatic truncation: terms below a coefficient threshold or above a Pauli
  weight cutoff are dropped, controlling the approximation/performance trade-off
- Optional loss channel support via [`LossyPauliSum`](/api/ppvm/paulisum/#ppvm.paulisum.LossyPauliSum). This also includes a reset channel,
  allowing to reset a qubit from being lost to its zero-state.
- Python bindings backed by a high-performance Rust core.

## Quick installation

Just `pip install` the directly from git:

```bash
pip install git+https://github.com/QuEraComputing/ppvm.git#subdirectory=ppvm-python
```

## Short example

```python
from ppvm import PauliSum

# Start with the terms for which you want to compute the average
state = PauliSum.new(n_qubits = 2, terms = ["ZZ"])

# Create a GHZ state
# NOTE: since this is Pauli Propagation, we need to backwards propagate, hence
# the CNOT needs to precede the Hadamard
state.cnot(0, 1)
state.h(0)

# The state is now represented as IZ
print(state)  # 1.000 * IZ

# Compute the average value
print(state.overlap_with_zero())
```

## Getting started

- [Installation](install.md)
- [Examples](examples/index.md)
- [API Reference](api/index.md)
- [Rust API Reference](../rust-api/)


## Similar packages

- [PauliPropagation.jl](https://github.com/MSRudolph/PauliPropagation.jl) - Julia package for Pauli Propagation. They also have a nice paper that explains the underlying theory: [arXiv:2505.21606](https://arxiv.org/abs/2505.21606)
- [cuPauliProp](https://docs.nvidia.com/cuda/cuquantum/latest/cupauliprop/index.html) - CUDA package by NVIDIA
