# ppvm-cli

Command-line front-end for the Pauli-propagation virtual machine. Parses,
dumps, runs, and steps through `.sst` programs (and their compiled `.ssb`
bytecode).

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

## Debug

`debug` steps through a program interactively. At each pause it prints the
program counter, the next instruction, and the measurements so far, then waits
for a command (type the letter and press Enter; a bare Enter steps):

- `s` — step one instruction
- `c` — continue to the next breakpoint (or the end)
- `q` — quit

By default it pauses at `breakpoint` instructions in the program. Add one
wherever you want execution to stop:

```
fn @main() {
    const.u64 0
    gate h
    breakpoint        // execution pauses here
    const.u64 0
    gate measure
    ret
}
```

```sh
$ printf 's\nc\n' | ppvm debug program.sst
-- breakpoint hit --
pc=3  next: const.u64 0
measurements: (none)
> s step | c continue | q quit: pc=4  next: Measure
measurements: (none)
> s step | c continue | q quit: Program finished.
Measurements: 0
```

To step through a program that has no breakpoints, pass `-b`/`--break-at-start`
to pause before the very first instruction:

```sh
$ ppvm debug examples/ghz.sst -b
pc=0  next: const.u64 0
measurements: (none)
> s step | c continue | q quit:
```

Batch `run` ignores `breakpoint` instructions entirely, so the same file still
runs straight through with `ppvm run`.
