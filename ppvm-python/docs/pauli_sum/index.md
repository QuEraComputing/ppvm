# Pauli Propagation

## Overview

[`PauliSum`][ppvm.paulisum.PauliSum] simulates quantum circuits in the **Heisenberg picture**.
Instead of evolving a state vector forward through the circuit, it propagates
*observables* backward. The observable is represented as a weighted sum of Pauli
strings, and each gate is applied analytically by conjugation.

This makes it possible to simulate large, deep circuits -- including noisy ones --
at a fraction of the cost of full statevector simulation.

### Basic usage

```python
from ppvm import PauliSum

# Observable: Z on each qubit
state = PauliSum.new(n_qubits=3, terms=[f"Z{i}" for i in range(3)])

# Apply gates in reverse circuit order
state.cnot(1, 2)
state.cnot(0, 1)
state.h(0)

# Expectation value with respect to |0...0>
print(state.overlap_with_zero())
```

!!! note
    Because this is Heisenberg-picture propagation, gates must be applied in
    **reverse circuit order**.

### Term notation

Terms can be specified as full Pauli strings (`"XZI"`) or in compact notation
(`"X0Z1"` -- Pauli + qubit index). Coefficients default to 1.0 but can be set
explicitly:

```python
ps = PauliSum.new(4, [("Z0Z1", 0.5), ("X2", 0.3)])
```

### Truncation

Two truncation strategies control the approximation/performance trade-off:

- **Coefficient truncation** (`min_abs_coeff`): drops terms with absolute
  coefficient below a threshold.
- **Weight truncation** (`max_pauli_weight`): drops terms with more non-identity
  Paulis than the cutoff.

```python
ps = PauliSum.new(10, "Z0", min_abs_coeff=1e-8, max_pauli_weight=5)
```

## Simulating loss

To simulate qubit loss, ppvm offers [`LossyPauliSum`][ppvm.paulisum.LossyPauliSum].
This is a dedicated class that behaves just like a [`PauliSum`][ppvm.paulisum.PauliSum],
but adds additional methods for the loss.

This separation exists because we need to extend the Pauli basis to
account for loss (see [Loss channel details](loss.md) for the full background). This comes at a storage overhead --
we now need 3 bits to represent a character in a Pauli string rather than 2.

Here is a small example:

```python
from ppvm import LossyPauliSum

ps = LossyPauliSum.new(n_qubits = 1, terms=["Z"])

# Reset at the end of the circuit
ps.reset_loss_channel(0)

# Loss after an X gate
ps.loss_channel(0, 0.1)

# Apply an X gate
ps.x(0)

z_exp = ps.overlap_with_zero()

# This will be -0.8: in 10% of cases we have <Z> = 1 instead of -1.
print(f"<Z>: {z_exp}")
```

A third truncation strategy is available for lossy simulations:

- **Loss weight truncation** (`max_loss_weight`): drops terms with more than a
  given number of `L` operators. Since the contribution of strings with `L` on
  many positions scales as $p_L^k$, this effectively controls the branching
  from loss and reset channels.

```python
ps = LossyPauliSum.new(3, "ZZZ", max_loss_weight=2)
```
