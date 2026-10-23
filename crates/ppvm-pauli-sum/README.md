# ppvm-pauli-sum

Pauli-propagation engine for ppvm: `PauliSum`, truncation strategies, and ready-made `Config` bundles. A `Config` built on `LossyPauliWord` makes `PauliSum` track qubit loss.

You start from an observable and apply the circuit's gates to it in reverse order (Heisenberg picture), then evaluate the result against an initial state. The number of terms can grow quickly, so truncation strategies such as `MaxPauliWeight`, `CoefficientThreshold` and `MaxLossWeight` drop terms as you go. `CombinedStrategy` applies two of them in sequence.

Part of [ppvm](https://github.com/QuEraComputing/ppvm), the Programmable Pauli Virtual Machine.

## Usage

```sh
cargo add ppvm-pauli-sum
```

- [Documentation](https://queracomputing.github.io/ppvm/)
- [API reference](https://docs.rs/ppvm-pauli-sum)

## License

Apache-2.0. See [LICENSE](https://github.com/QuEraComputing/ppvm/blob/main/LICENSE) and [NOTICE](https://github.com/QuEraComputing/ppvm/blob/main/NOTICE).
