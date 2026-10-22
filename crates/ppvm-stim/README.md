# ppvm-stim

Parse and run [Stim](https://github.com/quantumlib/Stim) circuits on the ppvm `GeneralizedTableau`.

Parse a circuit once with `parse_extended`, then run a single shot with `execute` or many with `sample`. Parsing is done by `stim-parser`, so circuits can also use ppvm's extensions for non-Clifford gates, qubit loss and leakage. The optional `rayon` feature parallelises multi-shot sampling.

Part of [ppvm](https://github.com/QuEraComputing/ppvm), the Programmable Pauli Virtual Machine.

## Usage

```sh
cargo add ppvm-stim
```

- [Documentation](https://queracomputing.github.io/ppvm/)
- [API reference](https://docs.rs/ppvm-stim)

## License

Apache-2.0. See [LICENSE](https://github.com/QuEraComputing/ppvm/blob/main/LICENSE) and [NOTICE](https://github.com/QuEraComputing/ppvm/blob/main/NOTICE).
