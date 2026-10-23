# ppvm-pauli-word

Packed Pauli-string types for the ppvm simulators: `PauliWord`, `PhasedPauliWord`, `LossyPauliWord` and `PauliPattern`.

`PhasedPauliWord` carries a phase alongside the word. `LossyPauliWord` adds an `L` symbol for a lost qubit. `PauliPattern` matches sets of words; for example, `Z?*` matches any word made of `I` and `Z`.

Part of [ppvm](https://github.com/QuEraComputing/ppvm), the Programmable Pauli Virtual Machine.

## Usage

```sh
cargo add ppvm-pauli-word
```

- [Documentation](https://queracomputing.github.io/ppvm/)
- [API reference](https://docs.rs/ppvm-pauli-word)

## License

Apache-2.0. See [LICENSE](https://github.com/QuEraComputing/ppvm/blob/main/LICENSE) and [NOTICE](https://github.com/QuEraComputing/ppvm/blob/main/NOTICE).
