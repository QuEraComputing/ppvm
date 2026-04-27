---
title: ppvm-stim Design
date: 2026-04-21
status: approved-in-chat
---

# ppvm-stim Design

## Goal

Extract the existing Stim circuit parser from `ppvm-tableau/src/stim.rs` into a new workspace crate, `ppvm-stim`, and redesign it so parsing, dialect normalization, and execution are clearly separated.

The extracted crate serves two audiences:

- ppvm's internal pipeline, which uses Stim as a circuit format plus a small set of ppvm-specific tags (`S[T]` → T gate, `I[R_X(theta=…)]` → Rx rotation, `I_ERROR[loss]` → loss channel, etc.) and executes against `GeneralizedTableau`.
- Pure-Stim consumers, who parse standard Stim circuits and walk the AST themselves without invoking internal-dialect resolution.

## Scope

Phase 1 (this spec):

- Extract parser + executor into `crates/ppvm-stim/`.
- Adopt a clean three-stage pipeline: `parse` → `normalize::to_tableau` → `execute`.
- Cover all Stim instructions that `GeneralizedTableau` already supports (see Coverage below).
- Large test suite, curated primarily from `quantumlib/Stim`, driving TDD.
- Wire `ppvm-python-native` directly to `ppvm-stim`; delete `ppvm-tableau/src/stim.rs`.
- Add a Rust-side `sample` function and Python-side `ppvm.sample_stim` for parse-once / normalize-once / execute-many-shots usage.

Phase 1 does not cover:

- Stim instructions not yet supported by `GeneralizedTableau` (`SWAP`, `MX`, `MRX`, `MPP`, `HERALDED_ERASE`, classical-feedback targets like `rec[-k]`, etc.). These fail at `normalize::to_tableau` with a specific `Unsupported` error; adding them is phase-2 work.
- Parser or executor benchmarks. Parse-over-text-file throughput is expected to be negligible relative to execution; we revisit only if measurement shows otherwise.
- Parallel shot sampling. A rayon-based `sample_parallel` is deferred.
- Round-tripping `Program` back to Stim source text.

## Guiding Principles

- **Readability over cleverness.** Prefer flat module layout, one big match, plain `Vec`, free functions. No abstraction without a measured reason.
- **Pure-Stim vs internal-dialect separated by design.** The parser has no internal-dialect knowledge; tag meaning is resolved by a separate normalizer.
- **Validation where Stim puts it.** Parse-time checks on instruction names, argument counts, and target arities — matching Stim's own `CircuitInstruction::validate()` flow. Normalizer handles semantics specific to our pipeline only.
- **Non-exhaustive public errors.** Error enums use `#[non_exhaustive]` so we can grow variants without SemVer breakage. This keeps open the path to richer typed errors when we eventually factor the parser into its own standalone crate.

## Architecture

### Crate layout

```
crates/ppvm-stim/
├── Cargo.toml              # deps: ppvm-tableau (rayon), chumsky = "0.12",
│                           #       thiserror, itertools
├── src/
│   ├── lib.rs              # public API re-exports; top-level Error;
│   │                       # run_string / run_file one-shot conveniences
│   ├── parser/
│   │   ├── mod.rs          # chumsky grammar + pub fn parse; one linear story
│   │   └── ast.rs          # Program, RawInstruction, GateName, NoiseName,
│   │                       # MeasureName, Tag, TagParam, AnnotationKind, ParseError;
│   │                       # GATE_TABLE / NOISE_TABLE / MEASURE_TABLE / ANNOTATION_TABLE
│   │                       # (arg/target arity metadata per instruction)
│   ├── normalize.rs        # Program → TableauProgram; dialect resolution;
│   │                       # phase-1 Unsupported rejection; NormalizeError
│   ├── tableau_program.rs  # TableauProgram, Instruction, GateKind, NoiseKind,
│   │                       # MeasureKind
│   └── executor.rs         # execute, sample, ExecError
├── tests/
│   ├── parser_*.rs         # one file per instruction family
│   ├── normalize.rs
│   ├── executor.rs
│   ├── stim_corpus.rs      # harness for tests/data/*.stim
│   └── data/               # curated .stim files from quantumlib/Stim
```

The parser lives under `src/parser/` with no dependency on `ppvm-tableau`. That boundary is what makes a future "lift parser into its own crate" refactor cheap. Everything above `parser` can use tableau types freely.

### Data flow

```
   source (&str / file)
          │ parse()
          ▼
     Program                 ← pure-Stim AST; tags preserved verbatim;
          │                    arg/target/name validation complete
          │ normalize::to_tableau()
          ▼
   TableauProgram            ← dialect resolved (S[T]→T, I[R_X]→RX,
          │                    I_ERROR[loss]→Loss, …); phase-1-unsupported
          │                    rejected; GateKind/NoiseKind/MeasureKind are
          │                    ready-to-execute enums
          │ execute(&prog, &mut tab)        ← single shot
          │ sample(&prog, num_shots, factory) ← N shots
          ▼
 Vec<Option<bool>>  /  Vec<Vec<Option<bool>>>
```

Each arrow is a free function. Functions take their inputs by immutable reference (except `execute`/`sample`, which mutate the tableau). The immutable borrow on `TableauProgram` is what enables shot reuse — parse and normalize happen once; execute runs thousands of times.

## Types

### Pure-Stim AST — `src/parser/ast.rs`

```rust
pub struct Program {
    pub instructions: Vec<RawInstruction>,
}

pub enum RawInstruction {
    Gate    { name: GateName,    tags: Vec<Tag>, args: Vec<f64>, targets: Vec<usize>, line: usize },
    Noise   { name: NoiseName,   tags: Vec<Tag>, args: Vec<f64>, targets: Vec<usize>, line: usize },
    Measure { name: MeasureName, tags: Vec<Tag>, args: Vec<f64>, targets: Vec<usize>, line: usize },
    Annotation { kind: AnnotationKind, args: Vec<f64>, targets: Vec<usize>, line: usize },
    Repeat  { count: u64, body: Vec<RawInstruction>, line: usize },
}

pub struct Tag {
    pub name: String,            // "T", "R_X", "loss", ...
    pub params: Vec<TagParam>,
}

pub enum TagParam {
    Positional(f64),
    Named { key: String, value: f64 },
}

pub enum GateName {    // phase-1 + Stim-valid-but-unsupported names known to the parser
    // Supported
    Reset, X, Y, Z, H, S, SDag, SqrtX, SqrtXDag, SqrtY, SqrtYDag,
    Identity, CX, CY, CZ,
    // Stim-valid, phase-1-unsupported (parser accepts, normalizer rejects)
    Swap, ISwap, ISwapDag, SqrtXX, SqrtYY, SqrtZZ,
    CXSwap, SwapCX, XCX, XCY, XCZ, YCX, YCY, YCZ,
    CXYZ, CZYX, HXY, HYZ,
    // add as we encounter them in the Stim corpus
}

pub enum NoiseName {
    // Supported
    Depolarize1, Depolarize2, PauliChannel1, PauliChannel2,
    XError, YError, ZError, IError,
    // Unsupported
    HeraldedErase, HeraldedPauliChannel1, CorrelatedError, ElseCorrelatedError,
}

pub enum MeasureName {
    // Supported
    M, MZ, MR,
    // Unsupported
    MX, MY, MRX, MRY, MXX, MYY, MZZ, MPP,
    // Note: MPAD is treated as an annotation, not a measurement (matches today's stim.rs).
}

pub enum AnnotationKind {
    Detector, MPad, ObservableInclude, QubitCoords, ShiftCoords, Tick,
}
```

A `GATE_TABLE` (and analogous `NOISE_TABLE`, `MEASURE_TABLE`) maps each name to its required arg count and target arity (`Any`, `Pairs`, `Quadruples`). The parser uses these to validate at parse time, matching Stim's own validation flow.

### Error types

```rust
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ParseError {
    #[error("syntax error at line {line}, col {col}: {message}")]
    Syntax { line: usize, col: usize, message: String },

    #[error("unknown instruction '{name}' at line {line}")]
    UnknownInstruction { name: String, line: usize },

    #[error("'{name}' at line {line} expected {expected} args, got {found}")]
    ArgCount { name: String, expected: usize, found: usize, line: usize },

    #[error("'{name}' at line {line} expected target count divisible by {divisor}, got {found}")]
    TargetCount { name: String, divisor: usize, found: usize, line: usize },
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum NormalizeError {
    #[error("unsupported instruction '{name}' at line {line} (phase 1)")]
    Unsupported { name: String, line: usize },

    #[error("invalid tag '{tag}' on '{instruction}' at line {line}: {message}")]
    InvalidTag { tag: String, instruction: String, line: usize, message: String },
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ExecError {
    // Empty in phase 1. Phase 2 will add cases like MeasurementRecordOutOfRange.
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    #[error(transparent)] Parse(#[from] ParseError),
    #[error(transparent)] Normalize(#[from] NormalizeError),
    #[error(transparent)] Exec(#[from] ExecError),
}
```

### Normalized AST — `src/tableau_program.rs`

```rust
pub struct TableauProgram {
    pub instructions: Vec<Instruction>,
    pub expected_measurement_count: usize,  // sum over M/MR, × REPEAT count; used to pre-size result buffers
}

pub enum Instruction {
    Gate    { kind: GateKind,   targets: Vec<usize>, line: usize },
    Noise   { kind: NoiseKind,  targets: Vec<usize>, args: Vec<f64>, line: usize },
    Measure { kind: MeasureKind, targets: Vec<usize>, line: usize },
    Annotation { /* phase-1 no-op */ },
    Repeat  { count: u64, body: Vec<Instruction>, line: usize },
}

pub enum GateKind {
    Reset,
    X, Y, Z, H, S, SDag, SqrtX, SqrtXDag, SqrtY, SqrtYDag,
    T, TDag,                                                    // from S[T] / S_DAG[T]
    RX { theta: f64 }, RY { theta: f64 }, RZ { theta: f64 },    // from I[R_X(theta=…)] etc.
    U3 { theta: f64, phi: f64, lambda: f64 },                   // from I[U3(…)]
    CX, CY, CZ,
}

pub enum NoiseKind {
    Depolarize1, Depolarize2,
    PauliChannel1, PauliChannel2,
    XError, YError, ZError,
    Loss,                 // from I_ERROR[loss]
    CorrelatedLoss,       // from I_ERROR[correlated_loss]
}

pub enum MeasureKind { M, MR }   // M and MZ both map to M at normalize time
```

## Parser — `src/parser/mod.rs`

Chumsky 0.12 grammar. Kept in a single file so a reader walks the grammar top-to-bottom:

```rust
fn grammar<'s>() -> impl Parser<'s, &'s str, Program, extra::Err<Rich<'s, char>>> {
    // atoms: ident, uint, pi_expr (replaces today's parse_pi_expr helper)
    // tags:  [ident], [ident(k=v,…)], [ident(f,…)]
    // args:  (f, f, …)
    // targets: uint (whitespace separated)
    // instruction: ident then tags? then args? then targets
    // repeat: "REPEAT" uint "{" instruction* "}"
    // program: (comment | blank | instruction | repeat)*
}
```

Parsing steps performed at this stage:

1. Grammar matches the line.
2. Instruction name is looked up in the corresponding family table (`GATE_TABLE` / `NOISE_TABLE` / `MEASURE_TABLE` / `ANNOTATION_TABLE`). Unknown name → `ParseError::UnknownInstruction`.
3. `args.len()` is validated against the table entry's expected arg count. Mismatch → `ParseError::ArgCount`.
4. `targets.len()` is validated against the table entry's target arity (`Any`, `Pairs`, `Quadruples`). Mismatch → `ParseError::TargetCount`.
5. Tags are parsed structurally but their meaning is not inspected. Tags are stored on the `RawInstruction`.
6. Line numbers are tracked as the grammar advances.

Chumsky's `Rich<char>` error is flattened into `ParseError::Syntax { line, col, message }` at the `pub fn parse` boundary. The rich structured spans are discarded; we can expose them in a future `ParseErrorDetailed` variant if needed when the parser moves to its own crate.

## Normalizer — `src/normalize.rs`

One function. All internal-dialect knowledge lives here:

```rust
pub fn to_tableau(program: &Program) -> Result<TableauProgram, NormalizeError>
```

Walks the `Program` once, building a `TableauProgram`. Key rules:

- `RawInstruction::Gate { name: S, tags: [Tag{T, []}], .. }` → `Instruction::Gate { kind: GateKind::T, .. }`.
- `RawInstruction::Gate { name: Identity, tags: [Tag{R_X, [Named{theta, θ}]}], .. }` → `Instruction::Gate { kind: GateKind::RX { theta: θ }, .. }`.
- `RawInstruction::Noise { name: IError, tags: [Tag{loss, []}], args: [p], .. }` → `Instruction::Noise { kind: NoiseKind::Loss, args: [p], .. }`.
- Unsupported names (`Swap`, `MX`, `MPP`, …) → `NormalizeError::Unsupported`.
- Malformed internal-dialect tags (`I[R_X]` without a theta, or `S[T, extra]`) → `NormalizeError::InvalidTag`.

`expected_measurement_count` is computed during this walk by adding `len(targets)` for each `M`/`MR`, multiplied by outer `REPEAT` counts.

## Executor — `src/executor.rs`

```rust
pub fn execute<T, I, C>(
    program: &TableauProgram,
    tab: &mut GeneralizedTableau<T, I, C>,
) -> Result<Vec<Option<bool>>, ExecError>
where /* same bounds as today's RunStim impl */
{
    let mut results = Vec::with_capacity(program.expected_measurement_count);
    execute_slice(&program.instructions, tab, &mut results)?;
    Ok(results)
}

fn execute_slice<T, I, C>(
    instructions: &[Instruction],
    tab: &mut GeneralizedTableau<T, I, C>,
    results: &mut Vec<Option<bool>>,
) -> Result<(), ExecError> {
    for instr in instructions {
        match instr {
            Instruction::Gate { kind, targets, .. } => match kind {
                GateKind::H  => targets.iter().for_each(|&q| tab.h(q)),
                GateKind::CX => targets.chunks_exact(2).for_each(|p| tab.cnot(p[0], p[1])),
                // ... one arm per GateKind variant
            },
            Instruction::Noise   { kind, targets, args, .. } => match kind { /* ... */ },
            Instruction::Measure { kind, targets, .. } => match kind {
                MeasureKind::M  => targets.iter().for_each(|&q| results.push(tab.measure(q))),
                MeasureKind::MR => targets.iter().for_each(|&q| {
                    let o = tab.measure(q);
                    if o == Some(true) { tab.x(q); }
                    results.push(o);
                }),
            },
            Instruction::Annotation { .. } => { /* phase-1 no-op */ }
            Instruction::Repeat { count, body, .. } => {
                for _ in 0..*count { execute_slice(body, tab, results)?; }
            }
        }
    }
    Ok(())
}
```

Key properties:

- `TableauProgram` is borrowed immutably, so shot reuse is zero-cost.
- No allocations in the hot loop beyond the pre-sized `results` vector.
- `chunks_exact(2)` and friends are safe because the parser already validated target parity.
- Trait bounds are verbatim from today's `RunStim for GeneralizedTableau` impl.

### Sample function

```rust
pub fn sample<T, I, C, F>(
    program: &TableauProgram,
    num_shots: usize,
    mut make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    F: FnMut() -> GeneralizedTableau<T, I, C>,
{
    (0..num_shots)
        .map(|_| {
            let mut tab = make_tableau();
            execute(program, &mut tab)
        })
        .collect()
}
```

Shot loop lives in Rust. Phase 2 may add a `sample_parallel` using rayon.

## Top-level API — `src/lib.rs`

```rust
// Primary three-stage pipeline:
pub fn parse(src: &str)                              -> Result<Program, ParseError>;
pub mod normalize { pub fn to_tableau(p: &Program) -> Result<TableauProgram, NormalizeError>; }
pub fn execute<T, I, C>(p: &TableauProgram, tab: &mut GeneralizedTableau<T, I, C>)
    -> Result<Vec<Option<bool>>, ExecError>;
pub fn sample<T, I, C, F>(p: &TableauProgram, num_shots: usize, make_tab: F)
    -> Result<Vec<Vec<Option<bool>>>, ExecError>
where F: FnMut() -> GeneralizedTableau<T, I, C>;

// One-shot convenience (re-parses every call; DO NOT use in shot loops):
pub fn run_string<T, I, C>(src: &str, tab: &mut GeneralizedTableau<T, I, C>)
    -> Result<Vec<Option<bool>>, Error>;
pub fn run_file<T, I, C>(path: &Path, tab: &mut GeneralizedTableau<T, I, C>)
    -> Result<Vec<Option<bool>>, Error>;
```

The crate-level rustdoc explicitly documents:

- The three-stage pipeline as the recommended path.
- `run_string` / `run_file` are for single-shot demos only, because they re-parse per call.
- A worked example of the N-shot pattern using `sample(...)`.

## Python bindings

`ppvm-python-native` gains a `ppvm-stim` dep and adds:

```python
from ppvm import StimProgram, GeneralizedTableau
import ppvm

# Parse + normalize once, return opaque handle
prog = StimProgram.from_file("surface_code.stim")
# or: StimProgram.parse(source_string)

# Single-shot on an existing tableau, no re-parse:
tab = GeneralizedTableau(n_qubits=50, tol=1e-10)
results = tab.run(prog)

# Multi-shot sampling — module-level function, shot loop in Rust:
shots = ppvm.sample_stim(prog, n_qubits=50, tol=1e-10, num_shots=10_000)
```

Under the hood:

- `PyStimProgram` wraps `ppvm_stim::TableauProgram`; one Python class, no generics to surface.
- `GeneralizedTableau_<N>.run(prog)` calls `ppvm_stim::execute(&prog.inner, &mut self.inner)`.
- `ppvm.sample_stim(prog, n_qubits, tol, num_shots)` dispatches to the appropriate generated tableau class (same dispatch already used by today's `GeneralizedTableau(n_qubits, tol)` constructor), then calls `ppvm_stim::sample(&prog.inner, num_shots, || GeneralizedTableau::new(n_qubits, tol))`.
- `parse` / `normalize` errors surface as a Python exception carrying the `Error::Display` string.

The existing Python API surface (`run_stim_string`, `run_stim_file` methods on the tableau) is replaced. Callers update to either `tab.run(StimProgram.parse(src))` or `ppvm.sample_stim(prog, ...)`. Existing Python tests in `ppvm-python/test/generalized_tableau/test_stim.py` are rewritten to the new API.

## Migration — existing `ppvm-tableau/src/stim.rs`

Deleted outright. Affected call sites:

- `crates/ppvm-tableau/src/lib.rs` — remove `pub mod stim;`.
- `crates/ppvm-tableau/tests/gates.rs` — Stim-driven tests migrate to `ppvm-stim/tests/`.
- `crates/ppvm-tableau/benches/tableau-msd-stim.rs` — update import to `ppvm_stim::run_string` (or the shot-based pattern if more appropriate).
- `crates/ppvm-python-native/src/interface_tableau.rs` — replace `run_stim_string` / `run_stim_file` methods with `run(prog)` and `sample_stim` dispatch.
- `crates/ppvm-python-native/ppvm_python_native.pyi` — update type stubs.
- `ppvm-python/src/ppvm/generalized_tableau.py` — update Python wrapper.
- `ppvm-python/test/generalized_tableau/test_stim.py` — update tests.
- `ppvm-python/docs/` — update `run_stim_*` references to `StimProgram` / `sample_stim` pattern.

Breaking these APIs is acceptable; the package is pre-1.0 and the whole pipeline is in flux.

## Coverage — phase 1

Inclusive (parser accepts, normalizer maps to a `GateKind`/`NoiseKind`/`MeasureKind`, executor dispatches):

- Reset: `R`, `RZ`.
- Single-qubit Clifford: `X`, `Y`, `Z`, `H`, `H_XZ`, `S`, `S_DAG`, `SQRT_Z`, `SQRT_Z_DAG`, `SQRT_X`, `SQRT_X_DAG`, `SQRT_Y`, `SQRT_Y_DAG`.
- Single-qubit tagged non-Clifford: `S[T]`, `S_DAG[T]`, `I[R_X(theta=…)]`, `I[R_Y(theta=…)]`, `I[R_Z(theta=…)]`, `I[U3(theta=…, phi=…, lambda=…)]`.
- Two-qubit Clifford: `CX`, `ZCX`, `CNOT`, `CY`, `ZCY`, `CZ`, `ZCZ`.
- Noise: `DEPOLARIZE1`, `DEPOLARIZE2`, `PAULI_CHANNEL_1`, `PAULI_CHANNEL_2`, `X_ERROR`, `Y_ERROR`, `Z_ERROR`.
- Loss: `I_ERROR[loss]`, `I_ERROR[correlated_loss]`.
- Measurement: `M`, `MZ`, `MR`.
- Annotation (no-op): `DETECTOR`, `MPAD`, `OBSERVABLE_INCLUDE`, `QUBIT_COORDS`, `SHIFT_COORDS`, `TICK`.
- Control flow: `REPEAT N { ... }` (parsed as nested AST; executor recurses; no parse-time inlining).

Excluded (parser accepts to the extent the Stim corpus exercises them, normalizer rejects with `NormalizeError::Unsupported`):

- Extra Cliffords: `SWAP`, `ISWAP`, `ISWAP_DAG`, `CXSWAP`, `SWAPCX`, `SQRT_XX`/`YY`/`ZZ`, `XCX`/`XCY`/`XCZ`, `YCX`/`YCY`/`YCZ`, `C_XYZ`, `C_ZYX`, `H_XY`, `H_YZ`.
- Measurements: `MX`, `MY`, `MRX`, `MRY`, `MPP`, `MXX`, `MYY`, `MZZ`.
- Heralded noise: `HERALDED_ERASE`, `HERALDED_PAULI_CHANNEL_1`.
- Classical control: `CX rec[-1] 0` style feedback, `CORRELATED_ERROR`, `ELSE_CORRELATED_ERROR`.
- Complex target types beyond plain qubit indices: `rec[-k]`, `sweep[k]`, Pauli-product targets like `X3*Y4`.

## Testing

TDD throughout. Each type/module has its tests written first.

```
tests/
├── parser_syntax.rs        # whitespace, comments, blank lines, REPEAT nesting,
│                           # tag bracket shape, args paren shape
├── parser_gates.rs         # every gate name parses
├── parser_noise.rs         # every noise name parses
├── parser_measure.rs       # every measurement name parses
├── parser_tags.rs          # tag shape: [ident], [ident(k=v)], [ident(f)], multi-tag
├── parser_errors.rs        # one test per ParseError variant
├── normalize.rs            # dialect resolution + Unsupported rejection
├── executor.rs             # end-to-end outcomes (GHZ parity, reset behavior,
│                           # tagged rotations, REPEAT correctness, sampling)
├── stim_corpus.rs          # walks tests/data/*.stim:
│                           #   - parse must succeed on all files
│                           #   - normalize succeeds OR returns a specific
│                           #     Unsupported error (asserted in a table)
│                           #   - files that normalize also execute
└── data/
    ├── README.md           # provenance: quantumlib/Stim commit SHA
    ├── ghz.stim            # hand-written
    ├── repetition_code_d3_r3.stim
    ├── surface_code_d3_r1.stim
    ├── color_code_d3.stim
    ├── tableau_basics.stim
    └── ... (~20–40 files)
```

### Corpus sourcing

Pull generously from `quantumlib/Stim`:

- `src/stim/io/stim_data_formats.test.cc` adjacent `.stim` fixtures.
- `src/stim/cmd/command_gen.test.cc` generated-circuit fixtures.
- `glue/stimcirq/test_circuits/` round-trip circuits.
- Small circuits from `doc/getting_started.ipynb`.

The `stim_corpus.rs` harness treats "file parses OK" and "file normalizes with expected result (OK or specific Unsupported)" as the two phase-1 assertions. Every phase-2 feature we add turns a skipped-corpus file into an executed one — free regression coverage over time.

### Implementation order (high-level)

1. Create crate, `Cargo.toml`, `ast.rs` skeleton with types only.
2. Parser grammar skeleton: single-gate parse (`H 0`), smallest possible `parser_syntax` test passes.
3. `GATE_TABLE` + arg/target validation, driven by `parser_errors.rs`.
4. Tag parsing, driven by `parser_tags.rs`.
5. REPEAT blocks, driven by nested `parser_syntax` tests.
6. Full Stim-gate/noise/measure/annotation coverage, test-family by test-family.
7. `normalize.rs`, driven by `normalize.rs` tests.
8. `executor.rs`, driven by `executor.rs` tests.
9. `sample`, driven by shot tests in `executor.rs`.
10. Wire `ppvm-python-native`; update Python tests and docs.
11. Delete `ppvm-tableau/src/stim.rs` and downstream references.
12. `stim_corpus.rs` as the final integration pass.

The implementation plan (produced by the `writing-plans` skill next) turns these into numbered tasks.

### Not tested in phase 1

- Parser benchmarks.
- `Program` → Stim text round-trip.
- Tableau gate correctness (that's `ppvm-tableau`'s surface).
- Parallel shot sampling.

## Open Questions for Phase 2

- Typed `ParseError` enum upgrade when the parser moves to a standalone crate.
- `sample_parallel` via rayon.
- Exposing chumsky's `Rich` error spans for richer diagnostics (would tie the public error type to chumsky version).
- Round-tripping `Program` → Stim source text.
- Unsupported-in-phase-1 gate coverage (`SWAP`, `ISWAP`, `MX`/`MY`, `MPP`, classical feedback, heralded noise).
