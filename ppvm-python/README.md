[![CI](https://github.com/QuEraComputing/ppvm/actions/workflows/ci.yml/badge.svg)](https://github.com/QuEraComputing/ppvm/actions/workflows/ci.yml)

# ppvm (Python)

**ppvm** (Programmable Pauli Virtual Machine) is a fast quantum circuit simulator with Python bindings backed by a high-performance Rust core.

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

### Stabilizer codes

Define a code by its signed commuting Pauli checks, including checks that mix
X, Y, and Z. The order of letters is physical-qubit order, starting with qubit
zero. Supply logical representatives to select a specific encoding convention.

```python
from ppvm import Code

five = Code.from_stabilizers(
    ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"],
    logical_xs=["XXXXX"], logical_zs=["ZZZZZ"],
)
assert (five.n, five.k) == (5, 1)
assert five.state_stabilizers([1])[-1] == "-ZZZZZ"

# One row of X bits, one row of Z bits, and a minus sign: the -ZZ code.
odd = Code.from_tableau([[0, 0]], [[1, 1]], signs=[1],
                        logical_xs=["XX"], logical_zs=["ZI"])
check, = odd.terminal_readout("ZZ").detectors
assert (check.indices, check.offset) == ((0, 1), True)
```

`Code.from_css(hx, hz)` accepts separate positive X and Z checks. All constructors
preserve redundant checks and reject inconsistent signs. `encoder_tableau()`
exports a forward Stim tableau with logical inputs last. `synthesize_state_prep()`,
`synthesize_encoder()`, and `synthesize_logical_clifford(logical_tableau)`
return verified Stim circuits for ideal Clifford operations. Terminal readout
covers one Pauli measurement per physical qubit; repeated syndrome rounds need
a separate measurement schedule.

## Documentation

Full documentation is available at the [ppvm docs site](https://queracomputing.github.io/ppvm/).

## Examples

- [Trotter simulation](docs/examples/trotter.py) -- Trotterized time evolution of the XZZ Ising chain
- [Magic State Distillation](docs/examples/msd.py) -- 85-qubit MSD circuit with the generalized stabilizer tableau
- [demo/](demo/) -- Benchmarking scripts
