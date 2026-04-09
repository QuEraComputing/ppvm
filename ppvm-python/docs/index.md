# ppvm

**ppvm** (Pauli Propagation and Virtual Machine) is a fast and easy-to-use quantum circuit simulator
backed by a high-performance Rust core with Python bindings.

It offers two complementary simulation approaches:

## Pauli Propagation

Instead of tracking a full quantum state vector, ppvm represents the quantum
state as a sum of Pauli strings with real coefficients.
Gates are applied by conjugating each term in
the sum, and noise channels are applied as super-operators. This makes it
possible to simulate large and deep quantum circuits, including noisy ones,
at a fraction of the cost of full statevector simulation.

- Clifford gates, arbitrary single- and two-qubit rotations, and common noise
  channels.
- Automatic truncation: terms below a coefficient threshold or above a Pauli
  weight cutoff are dropped, controlling the approximation/performance trade-off.
- Optional loss channel support via [`LossyPauliSum`][ppvm.paulisum.LossyPauliSum],
  including a reset channel to reset a qubit from being lost to its zero-state.

## Generalized Stabilizer Tableau

The [`GeneralizedTableau`][ppvm.generalized_tableau.GeneralizedTableau]
simulates quantum circuits in the Schrodinger picture using a generalized stabilizer decomposition.
Clifford operations are handled efficiently in the tableau, while non-Clifford gates (T gates, arbitrary
rotations) are tracked via a sparse coefficient vector that grows with the stabilizer rank.

- Clifford gates, T gates, and arbitrary single- and two-qubit rotations.
- Mid-circuit measurement returning [`MeasurementResult`][ppvm.generalized_tableau.MeasurementResult] (ZERO, ONE, or LOST).
- Noise channels: depolarizing, Pauli error, loss, and correlated loss.
- [STIM](https://github.com/quantumlib/Stim) circuit format support via `run_stim_string` and `run_stim_file`.
- `fork()` for branching into independent simulation trajectories.

## Quick installation

Just `pip install` directly from git:

```bash
pip install git+https://github.com/QuEraComputing/ppvm.git#subdirectory=ppvm-python
```

## Short examples

### Pauli Propagation

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

### Stabilizer Tableau

```python
from ppvm import GeneralizedTableau

tab = GeneralizedTableau(n_qubits=2)
tab.h(0)
tab.cnot(0, 1)

# Measure both qubits -- results are always correlated (Bell state)
r0 = tab.measure(0)
r1 = tab.measure(1)
print(f"Qubit 0: {r0}, Qubit 1: {r1}")
```

## Getting started

- [Installation](install.md)
- [Examples](examples/index.md)
- [API Reference](api/index.md)
- [Rust API Reference](../rust-api/)


## Similar packages

- [PauliPropagation.jl](https://github.com/MSRudolph/PauliPropagation.jl) - Julia package for Pauli Propagation. They also have a nice paper that explains the underlying theory: [arXiv:2505.21606](https://arxiv.org/abs/2505.21606)
- [cuPauliProp](https://docs.nvidia.com/cuda/cuquantum/latest/cupauliprop/index.html) - CUDA package by NVIDIA
- [SOFT Simulator](https://github.com/haoliri0/SOFT) - package for generalized tableau simulation with CUDA support. For background, see their paper, [arXiv:2512.23037](https://arxiv.org/pdf/2512.23037), and also the nice [report by Ted Yoder](https://www.scottaaronson.com/showcase2/report/ted-yoder.pdf)
