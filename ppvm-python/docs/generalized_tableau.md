# Generalized Stabilizer Tableau

## Overview

[`GeneralizedTableau`][ppvm.generalized_tableau.GeneralizedTableau] simulates quantum
circuits in the **Schrodinger picture** using a generalized stabilizer decomposition.

Clifford operations are tracked efficiently in the stabilizer tableau. Non-Clifford gates
(T gates, arbitrary rotations) expand the state into a superposition of stabilizer states,
tracked via a sparse coefficient vector. The cost scales exponentially only in the number
of non-Clifford gates, making it well-suited for circuits with few T gates
and many Clifford operations.

### Basic usage

```python
from ppvm import GeneralizedTableau

tab = GeneralizedTableau(n_qubits=2)
tab.h(0)
tab.cnot(0, 1)

r0 = tab.measure(0)
r1 = tab.measure(1)
print(f"Qubit 0: {r0}, Qubit 1: {r1}")  # always correlated (Bell state)
```

### Forking for multiple shots

Measurement mutates the tableau, so running multiple shots requires creating
independent copies. Use `fork()` to clone the quantum state with a fresh RNG:

```python
tab = GeneralizedTableau(n_qubits=2, seed=42)
tab.h(0)
tab.cnot(0, 1)

for shot in range(100):
    t = tab.fork(seed=shot)
    print(t.measure(0), t.measure(1))
```

To preserve the RNG state exactly (e.g. for checkpointing), use `copy.copy()` instead.

### STIM circuit support

Circuits in [STIM](https://github.com/quantumlib/Stim) format can be executed
directly via `run_stim_string` or `run_stim_file`. Only the gate/measurement
subset of STIM is supported; control flow instructions are not.

```python
from ppvm import GeneralizedTableau

tab = GeneralizedTableau(n_qubits=2)

stim_circuit = """
H 0
CX 0 1
M 0 1
"""

results = tab.run_stim_string(stim_circuit)
print(f"Bell state measurement: {results}")
```

For larger circuits stored as `.stim` files, use `run_stim_file`:

```python
tab = GeneralizedTableau(n_qubits=85)
results = tab.run_stim_file("path/to/circuit.stim")
```

## Noise and loss

`GeneralizedTableau` supports the same noise channels as `PauliSum`:

- **Depolarizing**: `depolarize(addr, p)` and `depolarize2(addr0, addr1, p)`.
- **Pauli error**: `pauli_error(addr, [p_x, p_y, p_z])`.
- **Loss**: `loss_channel(addr, p)` and `correlated_loss_channel(addr0, addr1, p)`.

When a qubit is lost, subsequent measurement returns
[`MeasurementResult.LOST`][ppvm.generalized_tableau.MeasurementResult]
instead of `ZERO` or `ONE`. You can check and reset loss state:

```python
tab = GeneralizedTableau(n_qubits=1, seed=0)
tab.loss_channel(0, 1.0)  # deterministic loss

print(tab.is_lost(0))   # True
print(tab.measure(0))    # MeasurementResult.LOST

tab.reset_loss_channel(0)
print(tab.is_lost(0))   # False
```

## Coefficient pruning

Each non-Clifford gate doubles the number of terms in the internal coefficient vector.
The `min_abs_coeff` parameter controls pruning of small coefficients:

```python
tab = GeneralizedTableau(n_qubits=4, min_abs_coeff=1e-8)
```
