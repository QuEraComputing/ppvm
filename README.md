# Pauli Propagation Virtual Machine

A fast quantum circuit simulator written in Rust, with Python bindings.

[![Docs](https://img.shields.io/badge/docs-6437FF)](https://congenial-bassoon-l436wp3.pages.github.io/)
[![CI - rust](https://github.com/QuEraComputing/ppvm/actions/workflows/rust-ci.yml/badge.svg)](https://github.com/QuEraComputing/ppvm/actions/workflows/rust-ci.yml)
[![CI - python](https://github.com/QuEraComputing/ppvm/actions/workflows/python-ci.yml/badge.svg)](https://github.com/QuEraComputing/ppvm/actions/workflows/python-ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

## Install

**Python** (with [uv](https://docs.astral.sh/uv/)):

```bash
uv add git+https://github.com/QuEraComputing/ppvm.git#subdirectory=ppvm-python
```

**Rust**:

```toml
[dependencies]
ppvm-pauli-sum = { git = "https://github.com/QuEraComputing/ppvm" }
```

## Examples

Pauli propagation runs **backwards** (Heisenberg picture): write gates in reverse order.

```python
from ppvm import PauliSum

state = PauliSum.new(n_qubits=2, terms=["ZZ"])
state.cnot(0, 1)   # GHZ preparation, written in reverse
state.h(0)

print(state)                    # 1.000 * IZ
print(state.overlap_with_zero())
```

The generalized stabilizer tableau is itself a form of Pauli propagation — it
tracks stabilizer generators under Heisenberg evolution, extended to handle
non-Clifford gates and measurements:

```python
from ppvm import GeneralizedTableau

tab = GeneralizedTableau(n_qubits=2)
tab.h(0)
tab.cnot(0, 1)

r0, r1 = tab.measure(0), tab.measure(1)
print(f"Qubit 0: {r0}, Qubit 1: {r1}")  # always correlated
```

See the [documentation](https://congenial-bassoon-l436wp3.pages.github.io/) for the Rust API, Stim integration, and symbolic propagation.

## License & contributing

Licensed under [Apache 2.0](LICENSE); see [NOTICE](NOTICE) for attribution. Contributions are welcome — read [`CONTRIBUTING.md`](CONTRIBUTING.md) and the [CLA](CLA.md) before opening a pull request.
