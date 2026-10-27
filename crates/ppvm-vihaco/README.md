# ppvm-vihaco

Run [vihaco](https://crates.io/crates/vihaco) circuit programs on the ppvm simulators.

It parses and compiles `.sst` programs to `.ssb` bytecode, and runs either form. `run_shots` runs many independent shots of a compiled program; with the optional `rayon` feature the shots run in parallel. This is the library behind the `ppvm` command-line tool in `ppvm-cli`.

Part of [ppvm](https://github.com/QuEraComputing/ppvm), the Programmable Pauli Virtual Machine.

## Usage

```sh
cargo add ppvm-vihaco
```

- [Documentation](https://queracomputing.github.io/ppvm/)
- [API reference](https://docs.rs/ppvm-vihaco)

## License

Apache-2.0. See [LICENSE](https://github.com/QuEraComputing/ppvm/blob/main/LICENSE) and [NOTICE](https://github.com/QuEraComputing/ppvm/blob/main/NOTICE).
