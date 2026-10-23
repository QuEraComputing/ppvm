# vihaco-circuit-isa

Quantum-circuit instruction set for the [vihaco](https://crates.io/crates/vihaco) virtual machine, used by ppvm.

It defines the `Circuit` component that vihaco programs call into: Pauli and Clifford gates, `T`, single- and two-qubit rotations, `U3`, measurement and reset, Pauli and depolarizing noise, qubit loss and leakage, plus `Truncate` and `Trace`. This crate only describes the instructions; `ppvm-vihaco` executes them.

Part of [ppvm](https://github.com/QuEraComputing/ppvm), the Programmable Pauli Virtual Machine.

## Usage

```sh
cargo add vihaco-circuit-isa
```

- [Documentation](https://queracomputing.github.io/ppvm/)
- [API reference](https://docs.rs/vihaco-circuit-isa)

## License

Apache-2.0. See [LICENSE](https://github.com/QuEraComputing/ppvm/blob/main/LICENSE) and [NOTICE](https://github.com/QuEraComputing/ppvm/blob/main/NOTICE).
