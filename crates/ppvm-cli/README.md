# ppvm-cli

Command-line front-end for the Pauli-propagation virtual machine. Parses,
dumps, and runs `.sst` programs (and their compiled `.ssb` bytecode).

## Install

From this crate directory:

```sh
cargo install --path crates/ppvm-cli
```

This builds the release binary and copies it to `~/.cargo/bin/ppvm`. As long as
`~/.cargo/bin` is on your `PATH`, you can then invoke `ppvm` from anywhere.

During development you can skip the install and use `cargo run` instead — just
put CLI arguments after `--`:

```sh
cargo run -p ppvm-cli -- run examples/ghz.sst
```

## Run

`run` executes a program and prints its measurement record. The example
[`examples/ghz.sst`](examples/ghz.sst) prepares a 3-qubit GHZ state and measures
every qubit, so each shot reads `0 0 0` or `1 1 1`:

```sh
$ ppvm run examples/ghz.sst
Measurements: 1 1 1
```

Each measurement event is shown as a bit string (a lost qubit prints as `L`),
events separated by spaces. Use `-f debug` for the raw record, or `-q` to
suppress the output entirely:

```sh
$ ppvm run examples/ghz.sst -f debug
Measurement record:
[[One], [One], [One]]
```

## Dump

`dump` compiles a `.sst` program to `.ssb` bytecode. With no `-o`, it writes
next to the input (`ghz.sst` → `ghz.ssb`):

```sh
$ ppvm dump examples/ghz.sst
Bytecode written to examples/ghz.ssb
```

`run` auto-detects the format from the file's contents, so the bytecode runs the
same way as the source:

```sh
$ ppvm run examples/ghz.ssb
Measurements: 0 0 0
```

`dump` refuses to overwrite an existing file unless you pass `-f`/`--force`.
