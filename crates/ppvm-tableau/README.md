# ppvm-tableau

Generalized stabilizer-tableau simulator for ppvm.

`Tableau` simulates Clifford circuits. `GeneralizedTableau` adds non-Clifford gates by tracking a sparse vector of coefficients over bitstrings, so circuits with few non-Clifford gates stay cheap. Both evolve forward in circuit order (Schrödinger picture) and support measurement, expectation values, and noise channels including qubit loss. The optional `rayon` feature parallelises parts of the simulation.

Part of [ppvm](https://github.com/QuEraComputing/ppvm), the Programmable Pauli Virtual Machine.

## Usage

```sh
cargo add ppvm-tableau
```

- [Documentation](https://queracomputing.github.io/ppvm/)
- [API reference](https://docs.rs/ppvm-tableau)

## License

Apache-2.0. See [LICENSE](https://github.com/QuEraComputing/ppvm/blob/main/LICENSE) and [NOTICE](https://github.com/QuEraComputing/ppvm/blob/main/NOTICE).
