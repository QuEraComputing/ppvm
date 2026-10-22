# ppvm-tableau-sum

Sums of generalized stabilizer tableaux for ppvm.

`GeneralizedTableauSum` represents a state as a weighted sum of `GeneralizedTableau`s. Measurement and noise split the state into branches; matching branches are merged and branches below a weight cutoff are dropped. Includes gates, measurement, noise channels and a `Sampler` for drawing shots.

Part of [ppvm](https://github.com/QuEraComputing/ppvm), the Programmable Pauli Virtual Machine.

## Usage

```sh
cargo add ppvm-tableau-sum
```

- [Documentation](https://queracomputing.github.io/ppvm/)
- [API reference](https://docs.rs/ppvm-tableau-sum)

## License

Apache-2.0. See [LICENSE](https://github.com/QuEraComputing/ppvm/blob/main/LICENSE) and [NOTICE](https://github.com/QuEraComputing/ppvm/blob/main/NOTICE).
