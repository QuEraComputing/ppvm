---
title: MPAD Implementation Design
date: 2026-04-30
branch: david/ppvm-stim-2
---

# MPAD Implementation Design

## Summary

`MPAD` is currently classified as an annotation in `ppvm-stim` and lowered to a no-op in the executor. This is incorrect: per the Stim spec, `MPAD` writes deterministic bits into the measurement record, and `rec[-k]` indices count those bits. Treating it as a no-op silently desynchronises the result vector relative to a reference Stim run and breaks any downstream consumer that depends on measurement-record alignment.

This spec promotes `MPAD` to a first-class instruction at every layer of the parse → normalize → execute pipeline:

- **Parser AST** gains `RawInstruction::MPad`; `AnnotationKind::MPad` is removed.
- **Parser table** reclassifies `"MPAD"` from `TableEntry::Annotation` to a new `TableEntry::MPad` with `ArgCount::Optional(1)` and `TargetArity::AtLeastOne`.
- **Normalized program** gains `Instruction::MPad { bits: Vec<bool>, noise: f64, line }`.
- **Normalize** validates each target is `0` or `1`, threads the optional probability arg into a `noise: f64` field (default `0.0`), and bumps `expected_measurement_count` by `bits.len() * enclosing_repeat_factor` (same formula as `Measure`).
- **Executor** appends one `Some(bit)` per target, optionally flipping each via `tab.bernoulli(noise)` to mirror the existing `M(p)` readout-noise plumbing.

## Context

### Current state

The legacy embedded parser at `crates/ppvm-tableau/src/stim.rs` (removed in `efbdf43`) treated `MPAD` as a no-op alongside `DETECTOR`, `OBSERVABLE_INCLUDE`, `QUBIT_COORDS`, `SHIFT_COORDS`, and `TICK`. The chumsky rewrite in `crates/ppvm-stim/` preserved that behaviour: `parser/ast.rs:158-164` lists `AnnotationKind::MPad`, `parser/table.rs:571-577` registers it as an `Annotation` with `ArgCount::None` / `TargetArity::Any`, and the executor's `Instruction::Annotation` arm at `executor.rs:242` is a no-op.

This is wrong on two axes:

1. **Argument shape.** Stim allows `MPAD(p) 0 1 0` where `p` is the probability that each padded bit is recorded flipped. The current table entry rejects any arg.
2. **Execution semantics.** `MPAD 0 1 0` should append three `Some(bool)` entries to the measurement record; the current implementation appends nothing.

Issue #44 lists `MPAD` among the instructions whose status to decide. The Stim spec (`Stim/doc/gates.md`, anchor `MPAD`) is unambiguous: MPAD is grouped textually with annotations but functionally writes to the measurement record.

### Authoritative Stim semantics

| Instruction | Affects measurement record? | Effect on simulator |
|---|---|---|
| `TICK` | No | None — visualizer hint |
| `QUBIT_COORDS` | No | None — drawing metadata |
| `SHIFT_COORDS` | No | None — accumulates a coord-offset for later annotations |
| `DETECTOR` | No (in measurement-sampling mode) | None — used in detector-sampling mode only |
| `OBSERVABLE_INCLUDE` | No (in measurement-sampling mode) | None — used in detector-sampling mode only |
| `MPAD` | **Yes** | Appends listed bits, optionally noisy |

Only `MPAD` requires a real implementation in measurement-sampling mode (which is the only mode `ppvm-stim` supports today). The other five remain legitimate no-ops in our pipeline.

## Decisions

| # | Decision |
|---|---|
| 1 | `MPAD` gets a dedicated AST variant at every layer rather than reusing `MeasureKind` or staying under `AnnotationKind`. |
| 2 | The optional probability argument is implemented in this round; deferral was rejected because reusing the existing `M(p)` plumbing costs essentially zero, and "spec-complete or rejected" is a cleaner contract than silent partial behaviour. |
| 3 | Non-`0`/`1` targets are caught in `normalize.rs`, not in the parser. The parser stays a syntactic-Stim recogniser; semantic validation lives where tag-validation already lives. |
| 4 | The parser carries targets as `Vec<usize>` and the normalized program carries `Vec<bool>`, with the typed conversion happening at the normalize step. This preserves the actual offending value for diagnostics. |
| 5 | `expected_measurement_count` is bumped for `MPAD` bits using the same `enclosing_repeat_factor` multiplier as `Measure`. This keeps the result-vector capacity hint aligned with what `execute` actually appends. |

## Architecture

### Pipeline

Unchanged shape:

```
src ──parse──▶ Program (RawInstruction)
            ──to_tableau──▶ TableauProgram (Instruction)
            ──execute──▶ Vec<Option<bool>>
```

`MPAD` joins as a parallel variant to `Measure` at each stage. The parser-side validation already in place (`ArgCount`, `TargetArity`) handles MPAD's shape requirements via a new `TableEntry::MPad` arm; no grammar change is needed because MPAD's lexical form (head, optional parens arg, integer targets) is already covered by the existing `instruction_line` combinator.

### Parser AST (`crates/ppvm-stim/src/parser/ast.rs`)

```rust
pub enum RawInstruction {
    // existing variants ...
    MPad {
        tags: Vec<Tag>,
        prob: Option<f64>,
        bits: Vec<usize>,
        line: usize,
    },
}
```

`AnnotationKind::MPad` is removed; the `canonical_name` arm for it is removed.

`bits` stays as `Vec<usize>` at the parser stage so that `NormalizeError::InvalidMPadTarget` can report the actual offending integer (e.g. `5`) rather than losing it to a parse-time bool conversion.

### Parser table (`crates/ppvm-stim/src/parser/table.rs`)

```rust
pub enum TableEntry {
    // existing variants ...
    MPad { args: ArgCount, targets: TargetArity },
}
```

The existing row:

```rust
("MPAD", TableEntry::Annotation { kind: AnnotationKind::MPad,
                                   args: ArgCount::None,
                                   targets: TargetArity::Any }),
```

is replaced with:

```rust
("MPAD", TableEntry::MPad { args: ArgCount::Optional(1),
                            targets: TargetArity::AtLeastOne }),
```

### Parser dispatch (`crates/ppvm-stim/src/parser/mod.rs`)

The `tolerate_non_numeric` flag at line 199 currently fires for any `TableEntry::Annotation`. Since `MPAD` no longer matches that variant, its targets must parse cleanly as `usize` — which is what we want (a stray `rec[-1]` on `MPAD` is a syntax error).

The `skip_arg_validation` flag at line 217 also currently triggers for annotations; `MPAD` falls through to the standard `ArgCount::Optional(1)` check, which is correct.

`build_instruction` (line 298) gains a new arm:

```rust
TableEntry::MPad { .. } => RawInstruction::MPad {
    tags,
    prob: args.into_iter().next(),
    bits: targets,
    line,
},
```

### Normalized program (`crates/ppvm-stim/src/tableau_program.rs`)

```rust
pub enum Instruction {
    // existing variants ...
    MPad {
        bits: Vec<bool>,
        noise: f64,
        line: usize,
    },
}
```

### Normalize (`crates/ppvm-stim/src/normalize.rs`)

`NormalizeError` gains:

```rust
#[error("'MPAD' at line {line} target #{index} = {value}, must be 0 or 1")]
InvalidMPadTarget { line: usize, index: usize, value: usize },
```

`normalize_slice` gains a new arm that:

1. Walks `bits` and converts each `usize` to `bool` (`0` → `false`, `1` → `true`), returning `InvalidMPadTarget` on the first non-bit value.
2. Sets `noise = prob.unwrap_or(0.0)`.
3. Bumps `*measure_count` by `bits.len() * (enclosing_repeat_factor as usize)` (saturating, mirroring the `Measure` arm).
4. Pushes `Instruction::MPad { bits, noise, line }`.

### Executor (`crates/ppvm-stim/src/executor.rs`)

`execute_slice` gains:

```rust
Instruction::MPad { bits, noise, .. } => {
    for &bit in bits {
        let recorded = if *noise > 0.0 && tab.bernoulli(*noise) {
            !bit
        } else {
            bit
        };
        results.push(Some(recorded));
    }
}
```

This matches the readout-noise pattern at `executor.rs:235` for `MR`, with the difference that the source bit is the deterministic pad value rather than a real measurement outcome. `tab.bernoulli` is the same RNG entry point already used by `MR`, so MPAD-recorded bits stay reproducible under the same shot seed.

### Compatibility

- `Instruction::Annotation` remains for the other five annotation kinds (`DETECTOR`, `OBSERVABLE_INCLUDE`, `QUBIT_COORDS`, `SHIFT_COORDS`, `TICK`). The executor's no-op arm for `Annotation` is unchanged.
- The stale comment at `parser/ast.rs:158` ("MPAD is treated as an annotation, not a measurement (matches today's `stim.rs`)") is deleted with `AnnotationKind::MPad`.
- The `test_stim_noop_instructions` comment at `tests/executor.rs:305` lists `MPAD` among no-op annotations even though the test source string does not include it; that stale mention is dropped from the comment.
- Public API (`pub use` in `lib.rs`) gains the new `Instruction::MPad` variant via the existing re-export of `Instruction`.

## Testing

New test cases in `crates/ppvm-stim/tests/`:

| Test | Purpose |
|---|---|
| `mpad_single_zero_appends_some_false` | `MPAD 0` → result vector ends with `Some(false)`. |
| `mpad_single_one_appends_some_true` | `MPAD 1` → result vector ends with `Some(true)`. |
| `mpad_multi_bit_in_order` | `MPAD 0 1 0 1` → four entries appended in order. |
| `mpad_inside_repeat_block` | `REPEAT 3 { MPAD 1 }` → three `Some(true)` entries; `expected_measurement_count` accounts for repeat factor. |
| `mpad_with_noise_distribution_within_3_sigma` | `MPAD(0.3) 0` over 4096 shots: flip rate within 3σ of `p`, mirroring `measure_noise_distribution_within_3_sigma` at `tests/executor.rs:383`. |
| `mpad_target_two_is_normalize_error` | `MPAD 2` → `Err(NormalizeError::InvalidMPadTarget { value: 2, .. })`. |
| `mpad_with_rec_target_is_syntax_error` | `MPAD rec[-1]` → `ParseError::Syntax(_)` (no longer tolerated). |
| `mpad_zero_targets_is_target_count_error` | `MPAD` (bare) → `ParseError::TargetCount` from the `TargetArity::AtLeastOne` rule. |
| `mpad_two_args_is_arg_count_error` | `MPAD(0.1, 0.2) 0` → `ParseError::ArgCount` from `ArgCount::Optional(1)`. |

The shot-loop noise test reuses `ppvm_stim::sample` with `GeneralizedTableau::new_with_seed`, matching the existing M-noise distribution test verbatim except for the input circuit.

## Out of scope

- Detector / observable sampling. `DETECTOR` and `OBSERVABLE_INCLUDE` remain no-ops; revisiting them is a separate spec.
- Confirming whether FLAIR currently emits `MPAD`. This implementation makes `MPAD` correct regardless; FLAIR-emission audit is a separate question tracked under #44's "Instruction coverage gaps" item.
- Any change to `tab.bernoulli` or the underlying RNG infrastructure.
- The other no-op annotations (`TICK`, `QUBIT_COORDS`, `SHIFT_COORDS`, `DETECTOR`, `OBSERVABLE_INCLUDE`) — their classification is correct and their `Instruction::Annotation` lowering stays unchanged.
