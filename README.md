# Programmable Pauli Virtual Machine

A fast quantum circuit simulator written in Rust, with Python bindings.

[![Docs](https://img.shields.io/badge/docs-6437FF)](https://queracomputing.github.io/ppvm/)
[![CI](https://github.com/QuEraComputing/ppvm/actions/workflows/ci.yml/badge.svg)](https://github.com/QuEraComputing/ppvm/actions/workflows/ci.yml)
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

See the [documentation](https://queracomputing.github.io/ppvm/) for the Rust API, Stim integration, and symbolic propagation.

## Citation

If you use PPVM in your research, please cite [arXiv:2610.18225](https://arxiv.org/abs/2610.18225)
(BibTeX also in [`CITATION.bib`](CITATION.bib)):

```bibtex
@misc{luo2026ppvmcircuitsimulationframework,
      title={PPVM: A Circuit Simulation Framework for Digital Twins of Quantum Processors},
      author={Xiu-Zhe Luo and David Plankensteiner and Kai-Hsin Wu and Rob Patterson and Alexander Schuckert and Oriol Rubies-Bigorda and Refaat Ismail and Jonhas Colina and Neelay Fruitwala and Danylo Lykov and Tyson Jones and Stefan Krastanov and Rafael Haenel and Stefan Ostermann and Sheng-Tao Wang},
      year={2026},
      eprint={2610.18225},
      archivePrefix={arXiv},
      primaryClass={quant-ph},
      url={https://arxiv.org/abs/2610.18225},
}
```

## License & contributing

Licensed under [Apache 2.0](LICENSE); see [NOTICE](NOTICE) for attribution. Contributions are welcome — read [`CONTRIBUTING.md`](CONTRIBUTING.md) and the [CLA](CLA.md) before opening a pull request.
