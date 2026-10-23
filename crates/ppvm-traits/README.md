# ppvm-traits

Core traits shared by the ppvm crates: the `Config` trait, the single-qubit `Pauli` alphabet, and the `ACMap`/`Trace` map implementations.

`Config` is a compile-time bundle that fixes the storage, coefficient type, truncation strategy, hasher, Pauli-word representation and map for one `PauliSum` type. The other ppvm crates are generic over it, so choosing a configuration costs nothing at runtime.

Most users will want a higher-level crate such as `ppvm-pauli-sum` or `ppvm-tableau` instead.

Part of [ppvm](https://github.com/QuEraComputing/ppvm), the Programmable Pauli Virtual Machine.

## Usage

```sh
cargo add ppvm-traits
```

- [Documentation](https://queracomputing.github.io/ppvm/)
- [API reference](https://docs.rs/ppvm-traits)

## License

Apache-2.0. See [LICENSE](https://github.com/QuEraComputing/ppvm/blob/main/LICENSE) and [NOTICE](https://github.com/QuEraComputing/ppvm/blob/main/NOTICE).
