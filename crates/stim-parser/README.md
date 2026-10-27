# stim-parser

Parser for the [Stim](https://github.com/quantumlib/Stim) circuit format.

It is written in pure Rust and does not depend on Stim, so you can parse, validate and print Stim circuits without the Stim C++ library or Python package.

`parse` returns the standard Stim syntax tree. `parse_extended` also understands ppvm's extensions, written as Stim instruction tags, for non-Clifford rotations, qubit loss, correlated loss and leakage. Errors are reported as diagnostics that point into the source. A printer turns either tree back into canonical Stim text, and parsing the printed output gives the same tree.

Part of [ppvm](https://github.com/QuEraComputing/ppvm), the Programmable Pauli Virtual Machine.

## Usage

```sh
cargo add stim-parser
```

- [Documentation](https://queracomputing.github.io/ppvm/)
- [API reference](https://docs.rs/stim-parser)

## License

Apache-2.0. See [LICENSE](https://github.com/QuEraComputing/ppvm/blob/main/LICENSE) and [NOTICE](https://github.com/QuEraComputing/ppvm/blob/main/NOTICE).
