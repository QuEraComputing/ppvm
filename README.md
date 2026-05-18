[![Documentation](https://img.shields.io/badge/Documentation-6437FF)](https://congenial-bassoon-l436wp3.pages.github.io/)

* Python build: [![CI - python](https://github.com/QuEraComputing/ppvm/actions/workflows/python-ci.yml/badge.svg)](https://github.com/QuEraComputing/ppvm/actions/workflows/python-ci.yml)
* Rust build: [![CI - rust](https://github.com/QuEraComputing/ppvm/actions/workflows/rust-ci.yml/badge.svg)](https://github.com/QuEraComputing/ppvm/actions/workflows/rust-ci.yml)

# ppvm - Pauli Propagation and Virtual Machine

**ppvm** is a fast quantum circuit simulator.
It is implemented in rust, but also offers [python bindings](#python).

# Short example

## Python

Install with [uv](https://docs.astral.sh/uv/):

```bash
uv add git+https://github.com/QuEraComputing/ppvm.git#subdirectory=ppvm-python
```

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

## Rust

Install with

```toml
[dependencies]
ppvm-runtime = { git = "https://github.com/QuEraComputing/ppvm" }
```

```rust
use ppvm_runtime::prelude::*;

fn main() {
    // Start with the terms for which you want to compute the average
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(2).build();
    state += ("ZZ", 1.0);

    // Create a GHZ state
    // NOTE: since this is Pauli Propagation, we need to backwards propagate, hence
    // the CNOT needs to precede the Hadamard
    state.cnot(0, 1);
    state.h(0);

    // The state is now represented as IZ
    println!("{}", state);  // 1.000 * IZ

    // Compute the average value
    let zero_state: PauliPattern = "Z?*".into();
    println!("{}", state.trace(&zero_state));  // 1
}
```

# License

ppvm is licensed under the [Apache License, Version 2.0](LICENSE).
See [NOTICE](NOTICE) for attribution requirements.

# Contributing

Contributions are welcome. By submitting a pull request to this repository,
you agree to license your contribution under the Apache License, Version 2.0
and to the terms of the [Contributor License Agreement](CLA.md).
