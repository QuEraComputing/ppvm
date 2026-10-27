# ppvm-tui

Ratatui components and app state for the `ppvm` terminal UI.

It provides widgets for the program listing, simulator state, measurement record and command line, plus a line editor. The crate does not own a terminal or run an event loop, so the components can be embedded in another ratatui app. The `ppvm` binary itself is in `ppvm-cli`.

Part of [ppvm](https://github.com/QuEraComputing/ppvm), the Programmable Pauli Virtual Machine.

## Usage

```sh
cargo add ppvm-tui
```

- [Documentation](https://queracomputing.github.io/ppvm/)
- [API reference](https://docs.rs/ppvm-tui)

## License

Apache-2.0. See [LICENSE](https://github.com/QuEraComputing/ppvm/blob/main/LICENSE) and [NOTICE](https://github.com/QuEraComputing/ppvm/blob/main/NOTICE).
