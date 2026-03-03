[![Documentation](https://img.shields.io/badge/Documentation-6437FF)](TODO)

* Python build: [![CI - python](https://github.com/QuEraComputing/bloqade/actions/workflows/ci.yml/badge.svg)](https://github.com/QuEraComputing/ppvm/actions/workflows/python-ci.yml)
* Rust build: [![CI - rust](https://github.com/QuEraComputing/bloqade/actions/workflows/ci.yml/badge.svg)](https://github.com/QuEraComputing/ppvm/actions/workflows/rust-ci.yml)

# ppvm - Pauli Propagation Virtual Machine

**ppvm** is a fast Pauli Propagation engine.
It is implemented in rust, but also offers [python bindings](#TODO).

# Short example

## Python

Install with

```bash
pip install git+https://github.com/QuEraComputing/ppvm.git#subdirectory=ppvm-python
```

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
