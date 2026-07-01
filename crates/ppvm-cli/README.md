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

`run` executes a program for one or more shots and prints the measurement
results. The example [`examples/ghz.sst`](examples/ghz.sst) prepares a 3-qubit
GHZ state and measures every qubit, so each shot reads `000` or `111`:

```sh
$ ppvm run examples/ghz.sst
000
```

Each shot is printed as a single flat bit string — `0`/`1`, with a lost qubit
shown as `2` — and shots are separated by newlines. Use `-s`/`--shots` to run
more than one:

```sh
$ ppvm run examples/ghz.sst --shots 5
000
000
111
000
111
```

Other options:

- `-t`/`--threads <N>` — run shots across `N` threads. More than one enables
  parallel execution (defaults to 1).
- `--seed <N>` — seed the RNG for reproducible results. The same seed yields the
  same shots regardless of the thread count.
- `-o`/`--output <FILE>` — write the results to a file (one shot per line)
  instead of stdout.
- `-f debug` — print the raw record for every shot instead of bit strings.
- `-q`/`--quiet` — run without printing anything.

```sh
$ ppvm run examples/ghz.sst --shots 2 -f debug
[[[Zero], [Zero], [Zero]], [[One], [One], [One]]]

$ ppvm run examples/ghz.sst --shots 1000 --threads 8 -o results.txt
Results written to results.txt
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
000
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
    circuit.h

    // execution pauses here
    breakpoint

    const.u64 0
    circuit.measure
    ret
}
```

```sh
$ printf 's\nc\n' | ppvm debug program.sst
-- breakpoint hit --
pc=3  next: const.u64 0
measurements:
> s step | c continue | q quit: pc=4  next: Measure
measurements:
> s step | c continue | q quit: Program finished.
Measurements: 0
```

To step through a program that has no breakpoints, pass `-b`/`--break-at-start`
to pause before the very first instruction:

```sh
$ ppvm debug examples/ghz.sst -b
pc=0  next: const.u64 0
measurements:
> s step | c continue | q quit:
```

Batch `run` ignores `breakpoint` instructions entirely, so the same file still
runs straight through with `ppvm run`.
