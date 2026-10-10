# Exact symbolic Pauli propagation

Compare PPVM's specialized trigonometric-polynomial coefficients with
PauliPropagation.jl's surrogate computation graph. Also compare each engine's
symbolic evaluation with its own repeated numerical propagation. No performance
tuning or production algorithm changes are included.

## Reproduce

```bash
cargo build --release --example symbolic_benchmark
JULIA_PKG_USE_CLI_GIT=true julia --project=benchmarks/symbolic -e 'using Pkg; Pkg.instantiate()'
uv run --python 3.12 python benchmarks/symbolic/run.py
uv run --with matplotlib python benchmarks/symbolic/plot.py
```

The Julia manifest pins PauliPropagation.jl 0.9.0 at
`b2125c7d2dc61626dbd50f44a47064f290dc5ce9`, including its dependency versions.
The default sweep uses 3 and 4 qubits, 1–6 Trotter steps, two parameter schemes,
seven trials, and batches of 2,000 evaluations. Override these with `--qubits`,
`--steps`, `--trials`, `--batch`, and `--out`; plotting selects four qubits by
default (`--qubits`). The Rust runner uses one byte of Pauli-word storage and
supports 1–8 qubits. Exact symbolic construction can grow exponentially with
depth: extend the sweep deliberately rather than using large circuits by default.

## Circuit and validation

An open chain has one forward Trotter step consisting of `RX(q, theta)` on
sites in ascending order, then `RZZ(q, q+1, theta)` on bonds in ascending order.
Rotations use `exp(-i theta P/2)`. Both engines propagate `Z_0` backwards through
the same forward circuit and compute its expectation in `|0...0>`.

- **Independent:** every gate has its own parameter.
- **Shared:** all RX gates use one parameter; all RZZ gates use another,
  including across steps. Julia receives the corresponding expanded gate-angle
  array outside timed regions; PPVM uses repeated variable IDs.

There is no Pauli-weight, sine-degree, or numerical-magnitude truncation.
PPVM's formal monomial arithmetic uses floating-point scalar coefficients;
its built-in roundoff handling is not a rigorous symbolic exact-arithmetic
backend. Here the formal scalar coefficients are integer-valued, and numerical
comparisons allow absolute roundoff error below `1e-10`.

Before timing every case, each runner compares *all* symbolic coefficients
against its numerical propagation at five deterministic parameter assignments.
The driver additionally compares full coefficient vectors across libraries,
including absent terms as zero, and checks the expectation against an independent
small dense-state simulation in forward gate order. The parameter generator is
`theta_i = 0.07 + ((17*i + 11*sample) % 53)/100`, indexed from zero. Timed batches
cycle through 53 different assignments; no evaluated-result cache is reused.

The Rust example includes analytic one-qubit and repeated-parameter substitution
checks, runnable with `cargo test --example symbolic_benchmark`.

## Timing contract

All runs use a single CPU thread, optimized Rust builds, and warmed Julia code.
Parsing, circuit generation, parameter assignment construction, correctness
checks, Julia startup and compilation are outside the reported timings.

- **Construction:** symbolic propagation plus readout preparation. PPVM forms
  the polynomial expectation; Julia's `zerofilter!` keeps the graph roots needed
  for the same expectation. Neither engine reevaluates the complete operator
  when the benchmark only asks for an expectation.
- **Evaluation:** evaluate that prepared expectation for each parameter set,
  including resetting Julia's graph evaluation flags. Report seconds per set.
- **Numeric:** propagate the original observable numerically and extract its
  expectation for each parameter set, using the same engine and no truncation.

Report medians of seven trials. Julia runs a full collection before each trial;
allocation and collection during measured execution are included. Construction
object destruction is excluded from construction timing. The numeric baseline
uses PPVM's FxHashMap configuration and Julia's dictionary-backed `PauliSum`,
not Julia's alternative vector storage. Thus these are comparisons of these
specific implementations and small circuits, not universal rankings.

The amortization count is `ceil(construction / (numeric - evaluation))`;
it is left blank when evaluation is not faster than numerical propagation.
Peak process RSS is not reported: runtime/JIT overhead would make it misleading
as a comparison of the expression representations themselves.

## Results and files

- `results/samples.csv`: every raw timing trial.
- `results/summary.csv`: medians, observable, support size, amortization count.
- `results/validation.csv`: maximum cross-library and dense-state errors.
- `results/metadata.json`: hardware, toolchain, revisions, runner checksums.
- `results/symbolic_propagation.pdf` and `.png`: construction/evaluation figure.
- `../../examples/symbolic_benchmark.rs`: Rust runner and analytic checks.
- `surrogate.jl`: Julia surrogate and numerical runner.
- `run.py`: validation and measurement orchestration (Python standard library).
- `plot.py`: paper-style plot in the QuEra palette (matplotlib).

PauliPropagation.jl's surrogate reference:
https://github.com/SparqleSim/PauliPropagation.jl/blob/b2125c7d2dc61626dbd50f44a47064f290dc5ce9/src/Surrogate/propagate.jl

### Recorded run (Apple M4 Pro)

Rust 1.96.0, Julia 1.12.6, PauliPropagation.jl 0.9.0; single-threaded.
All 24 cases passed full coefficient-vector validation (maximum cross-library
error `5.56e-16`) and independent dense-state validation (`1.23e-15`).
At four qubits and six steps:

| Parameters | Engine | Construction (µs) | Evaluation (µs/set) | Numeric (µs/set) | Amortization count |
|---|---|---:|---:|---:|---:|
| Shared | PPVM | 114.46 | 0.359 | 1.233 | 131 |
| Shared | PauliPropagation.jl | 55.75 | 26.884 | 13.408 | None |
| Independent | PPVM | 421.04 | 5.040 | 1.279 | None |
| Independent | PauliPropagation.jl | 53.92 | 26.968 | 13.225 | None |

Shared variables allow identical products from different propagation paths to
merge in PPVM's coefficient polynomials. At this depth, shared-angle evaluation
is 3.43 times faster than PPVM numerical propagation, but construction is
2.05 times slower than Julia's surrogate. Independent parameters expand many
more distinct monomials: PPVM construction is 7.81 times slower than the
surrogate, and evaluation is slower than PPVM's numerical propagation.
The surrogate evaluation entry point carries substantial fixed overhead on these
small problems; this experiment does not establish its behavior on larger
circuits. Both representations retain exact parameter dependence up to roundoff;
neither guarantees a benefit for every circuit or parameter sweep.
