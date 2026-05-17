[![CI - python](https://github.com/QuEraComputing/ppvm/actions/workflows/python-ci.yml/badge.svg)](https://github.com/QuEraComputing/ppvm/actions/workflows/python-ci.yml)

# ppvm (Python)

**ppvm** (Pauli Propagation and Virtual Machine) is a fast quantum circuit simulator with Python bindings backed by a high-performance Rust core.

## Installation

We recommend using [uv](https://docs.astral.sh/uv/) to manage your Python
environment. Install directly from git:

```bash
uv add git+https://github.com/QuEraComputing/ppvm.git#subdirectory=ppvm-python
```

Or clone and install locally:

```bash
git clone https://github.com/QuEraComputing/ppvm
cd ppvm
uv sync --project ppvm-python
```

## Quick examples

### Pauli Propagation

```python
from ppvm import PauliSum

state = PauliSum.new(n_qubits=2, terms=["ZZ"])

# Backwards-propagate through a GHZ circuit
state.cnot(0, 1)
state.h(0)

print(state)  # 1.000 * IZ
print(state.overlap_with_zero())
```

### Generalized Stabilizer Tableau

```python
from ppvm import GeneralizedTableau

tab = GeneralizedTableau(n_qubits=2)
tab.h(0)
tab.cnot(0, 1)

r0 = tab.measure(0)
r1 = tab.measure(1)
print(f"Qubit 0: {r0}, Qubit 1: {r1}")  # always correlated
```

## Documentation

Full documentation is available at the [ppvm docs site](https://queracomputing.github.io/ppvm/).

## Examples

- [Trotter simulation](docs/examples/trotter.py) -- Trotterized time evolution of the XZZ Ising chain
- [Magic State Distillation](docs/examples/msd.py) -- 85-qubit MSD circuit with the generalized stabilizer tableau
- [demo/](demo/) -- Benchmarking scripts
