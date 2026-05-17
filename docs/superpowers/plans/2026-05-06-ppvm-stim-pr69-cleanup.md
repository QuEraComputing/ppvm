# PPVM-stim PR #69 Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Apply targeted simplifications to PR #69 (`david/ppvm-stim-2`): drop redundant Raw types in the parser, eliminate per-instruction clones, factor out duplicated executor logic, clean up the Python binding wrapper, and restore lost test coverage from the deleted `ppvm-tableau/tests/gates.rs`.

**Architecture:** Eleven independent commits, ordered bottom-up: parser/AST internals first (Steps 1–4), executor (Steps 5–7), bindings (Steps 8–9), restored tests (Step 10), then a final structural change to `ExtendedInstruction` with a documented bail criterion (Step 11). PR #57 territory (where-clauses, trait organization in `ppvm-runtime`/`ppvm-tableau`) is excluded.

**Tech Stack:** Rust 2024 edition, chumsky 0.12, PyO3 0.27, Python 3 (uv-managed), maturin.

---

## Pre-flight

Before starting, confirm the working tree is clean and on the right branch.

- [ ] **Step 0.1: Confirm branch and clean tree**

Run:
```bash
git status
git rev-parse --abbrev-ref HEAD
```
Expected: clean working tree, branch `david/ppvm-stim-2`.

- [ ] **Step 0.2: Baseline test run**

Run:
```bash
cargo test --workspace --exclude ppvm-python-native
```
Expected: all tests pass. Note any pre-existing failures so we don't blame them on our work later.

---

## Step 1: Drop `RawTag` / `RawTagParam`

**Goal:** Have the grammar produce `ast::Tag` / `ast::TagParam` directly. Delete the Raw types and the per-tag conversion loop in `validate_node`.

**Files:**
- Modify: `crates/stim-parser/src/parser.rs:41-57` (delete `RawTag`/`RawTagParam` types)
- Modify: `crates/stim-parser/src/parser.rs:135-148` (drop the conversion loop)
- Modify: `crates/stim-parser/src/parser.rs:26-39` (change `RawSyntaxNode::Instruction.tags` type)
- Modify: `crates/stim-parser/src/grammar.rs:97-122` (emit `Tag`/`TagParam`)
- Modify: `crates/stim-parser/src/parser.rs:340-351` (test helper `instr_with_tags`)
- Modify: `crates/stim-parser/src/parser.rs:480-506` (test `raw_tags_are_converted_to_ast_tags`)
- Modify: `crates/stim-parser/src/grammar.rs:322-347` (tests using `RawTag`/`RawTagParam`)

- [ ] **Step 1.1: Make grammar emit AST `Tag`/`TagParam` directly**

In `crates/stim-parser/src/grammar.rs`, replace the `use crate::parser::{RawTag, RawTagParam};` line and `tag_param`/`tag`/`tags_block` definitions:

```rust
use crate::ast::{Tag, TagParam};

/// `<ident>=<pi_expr>` (Named) or `<pi_expr>` (Positional).
pub(crate) fn tag_param<'src>() -> impl Parser<'src, &'src str, TagParam, Extra<'src>> + Clone {
    let named = ident()
        .then_ignore(inline_pad())
        .then_ignore(just('='))
        .then_ignore(inline_pad())
        .then(pi_expr())
        .map(|(key, value)| TagParam::Named { key, value });
    let positional = pi_expr().map(TagParam::Positional);
    choice((named, positional))
}

/// Tag: `<ident>` or `<ident>(<tag_param>, ...)`.
pub(crate) fn tag<'src>() -> impl Parser<'src, &'src str, Tag, Extra<'src>> + Clone {
    let params = tag_param()
        .separated_by(inline_pad().then(just(',')).then(inline_pad()))
        .allow_trailing()
        .collect::<Vec<_>>()
        .delimited_by(just('(').then(inline_pad()), inline_pad().then(just(')')));
    ident().then(params.or_not()).map(|(name, params)| Tag {
        name,
        params: params.unwrap_or_default(),
    })
}

/// `[tag, tag, ...]`.
pub(crate) fn tags_block<'src>() -> impl Parser<'src, &'src str, Vec<Tag>, Extra<'src>> + Clone {
    tag()
        .separated_by(inline_pad().then(just(',')).then(inline_pad()))
        .allow_trailing()
        .collect::<Vec<_>>()
        .delimited_by(just('[').then(inline_pad()), inline_pad().then(just(']')))
}
```

Then update `instruction_head`'s return type at `crates/stim-parser/src/grammar.rs:162-177` — change `Vec<RawTag>` to `Vec<Tag>`. Same for `instruction_line`.

- [ ] **Step 1.2: Update `RawSyntaxNode::Instruction` to carry `Tag`**

In `crates/stim-parser/src/parser.rs`, change line 29 from `tags: Vec<RawTag>` to `tags: Vec<Tag>` and add `Tag` to the imports at line 4 (it's already there). Then delete the `RawTag` and `RawTagParam` types at lines 41-51.

- [ ] **Step 1.3: Drop the conversion loop in `validate_node`**

In `crates/stim-parser/src/parser.rs`, delete the `let tags: Vec<Tag> = tags.into_iter().map(...)` block at lines 135-148. The `tags` variable destructured from `RawSyntaxNode::Instruction` now already has type `Vec<Tag>` (Step 1.2 changed the field type), so subsequent uses just work as-is. No replacement code, no comment needed.

- [ ] **Step 1.4: Update test helpers**

In `crates/stim-parser/src/parser.rs`, the `instr_with_tags` helper at line 340 takes `tags: Vec<RawTag>`. Update its signature:

```rust
    fn instr_with_tags(name: &str, tags: Vec<Tag>) -> RawSyntaxNode {
```

Update the test `raw_tags_are_converted_to_ast_tags` at line 480 — rename to `tags_pass_through_validator` and rewrite to construct `Tag` / `TagParam` directly:

```rust
    #[test]
    fn tags_pass_through_validator() {
        let nodes = vec![instr_with_tags(
            "H",
            vec![Tag {
                name: "R".to_string(),
                params: vec![
                    TagParam::Positional(0.5),
                    TagParam::Named {
                        key: "theta".to_string(),
                        value: 0.25,
                    },
                ],
            }],
        )];
        let result = validate_program(nodes, &lm()).unwrap();
        match &result[0] {
            RawInstruction::Gate { tags, .. } => {
                assert_eq!(tags[0].name, "R");
                assert!(matches!(tags[0].params[0], TagParam::Positional(0.5)));
                assert!(matches!(
                    &tags[0].params[1],
                    TagParam::Named { key, value } if key == "theta" && *value == 0.25
                ));
            }
            other => panic!("{other:?}"),
        }
    }
```

- [ ] **Step 1.5: Update grammar tests that match on `RawTagParam`**

In `crates/stim-parser/src/grammar.rs`, the test module imports `RawTagParam` indirectly through `tag()`. Two affected tests:
- `tag_with_positional_params` (line 329)
- `tag_with_named_param` (line 337)

Both pattern-match `RawTagParam::Positional(...)` / `RawTagParam::Named { ... }`. Change to `TagParam::Positional(...)` / `TagParam::Named { ... }`. Add `use crate::ast::TagParam;` near the top of the test module (line ~257).

- [ ] **Step 1.6: Run parser tests**

Run:
```bash
cargo test -p stim-parser
```
Expected: all tests pass. If a test fails to compile because it still references `RawTag`/`RawTagParam`, update it to use `Tag`/`TagParam`.

- [ ] **Step 1.7: Run full workspace tests**

Run:
```bash
cargo test --workspace --exclude ppvm-python-native
```
Expected: all tests pass.

- [ ] **Step 1.8: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/stim-parser/src/parser.rs crates/stim-parser/src/grammar.rs
git commit -m "$(cat <<'EOF'
refactor(stim-parser): drop RawTag/RawTagParam, emit ast::Tag directly

Grammar combinators now produce `ast::Tag` / `ast::TagParam` directly.
Removes the crate-private mirror types and the per-tag conversion loop
in `validate_node`.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Step 2: Stop cloning in `interpret_one`

**Goal:** Change `interpret(prog: &Program)` to take `Program` by value, then move fields out of each `RawInstruction` instead of cloning.

**Files:**
- Modify: `crates/stim-parser/src/extended/interpret.rs:7-21` (interpret takes ownership)
- Modify: `crates/stim-parser/src/extended/interpret.rs:23-84` (interpret_one takes owned RawInstruction)
- Modify: `crates/stim-parser/src/extended/parser.rs:26-29` (parse_extended passes owned)
- Modify: `crates/stim-parser/src/extended/interpret.rs:86-271` (helpers take owned slices/vecs as appropriate)

- [ ] **Step 2.1: Make `interpret_one` take owned `RawInstruction`**

Replace the entire `interpret_one` function in `crates/stim-parser/src/extended/interpret.rs:23-84` with:

```rust
fn interpret_one(raw: RawInstruction) -> Result<ExtendedInstruction, ExtendedParseError> {
    match raw {
        RawInstruction::Gate {
            name,
            tags,
            args,
            targets,
            line,
        } => interpret_gate(name, tags, args, targets, line),
        RawInstruction::Noise {
            name,
            tags,
            args,
            targets,
            line,
        } => interpret_noise(name, tags, args, targets, line),
        RawInstruction::Measure {
            name,
            tags,
            args,
            targets,
            line,
        } => Ok(ExtendedInstruction::Measure {
            name,
            tags,
            args,
            targets,
            line,
        }),
        RawInstruction::Annotation {
            kind,
            args,
            targets,
            line,
        } => Ok(ExtendedInstruction::Annotation {
            kind,
            args,
            targets,
            line,
        }),
        RawInstruction::MPad {
            tags,
            prob,
            bits,
            line,
        } => Ok(ExtendedInstruction::MPad {
            tags,
            prob,
            bits: convert_mpad_bits(&bits, line)?,
            line,
        }),
        RawInstruction::Repeat { count, body, line } => {
            let mut inner = Vec::with_capacity(body.len());
            interpret_slice(body, &mut inner)?;
            Ok(ExtendedInstruction::Repeat {
                count,
                body: inner,
                line,
            })
        }
    }
}
```

- [ ] **Step 2.2: Update `interpret_slice` and `interpret` to consume**

Replace `interpret` and `interpret_slice` at `crates/stim-parser/src/extended/interpret.rs:7-21` with:

```rust
pub(crate) fn interpret(prog: Program) -> Result<ExtendedProgram, ExtendedParseError> {
    let mut out = Vec::with_capacity(prog.instructions.len());
    interpret_slice(prog.instructions, &mut out)?;
    Ok(ExtendedProgram { instructions: out })
}

fn interpret_slice(
    src: Vec<RawInstruction>,
    out: &mut Vec<ExtendedInstruction>,
) -> Result<(), ExtendedParseError> {
    for raw in src {
        out.push(interpret_one(raw)?);
    }
    Ok(())
}
```

- [ ] **Step 2.3: Update helper signatures to take owned where it eliminates a clone**

The owned-friendly version of `interpret_gate` takes `tags: Vec<Tag>`, `args: Vec<f64>`, `targets: Vec<usize>` instead of slices, so it can move them into `ExtendedInstruction::Gate` without `to_vec()`.

Replace `interpret_gate` at `crates/stim-parser/src/extended/interpret.rs:86-133` with:

```rust
fn interpret_gate(
    name: GateName,
    tags: Vec<Tag>,
    args: Vec<f64>,
    targets: Vec<usize>,
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    use GateName::*;

    match (name, tags.as_slice()) {
        (S, [t]) if t.name == "T" => {
            require_no_params(t, "S", line)?;
            Ok(ExtendedInstruction::T { targets, line })
        }
        (SDag, [t]) if t.name == "T" => {
            require_no_params(t, "S_DAG", line)?;
            Ok(ExtendedInstruction::TDag { targets, line })
        }
        (S, [t]) => Err(invalid_tag(t.name.clone(), "S", line, "expected [T]")),
        (SDag, [t]) => Err(invalid_tag(t.name.clone(), "S_DAG", line, "expected [T]")),
        (S, _) | (SDag, _) if !tags.is_empty() => Err(invalid_tag(
            tags[0].name.clone(),
            if matches!(name, S) { "S" } else { "S_DAG" },
            line,
            "expected exactly one tag",
        )),
        (Identity, [t]) => interpret_identity_tag(t, targets, line),
        (Identity, _) if !tags.is_empty() => Err(invalid_tag(
            tags[0].name.clone(),
            "I",
            line,
            "expected exactly one tag",
        )),
        _ => Ok(ExtendedInstruction::Gate {
            name,
            tags,
            args,
            targets,
            line,
        }),
    }
}
```

Note: `tags.as_slice()` is used for matching since match patterns can't bind on `Vec<T>`. The `Identity` arm passes `t` (a `&Tag` borrowed from the slice) and `targets` (owned) into `interpret_identity_tag`.

- [ ] **Step 2.4: Update `interpret_identity_tag` to take owned `targets`**

Replace at `crates/stim-parser/src/extended/interpret.rs:147-186`:

```rust
fn interpret_identity_tag(
    tag: &Tag,
    targets: Vec<usize>,
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    let axis = match tag.name.as_str() {
        "R_X" => Some(Axis::X),
        "R_Y" => Some(Axis::Y),
        "R_Z" => Some(Axis::Z),
        _ => None,
    };

    if let Some(axis) = axis {
        let [theta] = exact_named_params(tag, ["theta"], "I", line)?;
        return Ok(ExtendedInstruction::Rotation {
            axis,
            theta,
            targets,
            line,
        });
    }

    if tag.name == "U3" {
        let [theta, phi, lambda] = exact_named_params(tag, ["theta", "phi", "lambda"], "I", line)?;
        return Ok(ExtendedInstruction::U3 {
            theta,
            phi,
            lambda,
            targets,
            line,
        });
    }

    Err(invalid_tag(
        tag.name.clone(),
        "I",
        line,
        "unrecognized tag (expected R_X / R_Y / R_Z / U3)",
    ))
}
```

- [ ] **Step 2.5: Update `interpret_noise` to take owned**

Replace at `crates/stim-parser/src/extended/interpret.rs:188-271`:

```rust
fn interpret_noise(
    name: NoiseName,
    tags: Vec<Tag>,
    args: Vec<f64>,
    targets: Vec<usize>,
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    use NoiseName::*;

    match (name, tags.as_slice()) {
        (IError, [t]) if t.name == "loss" => {
            require_no_params(t, "I_ERROR", line)?;
            if args.len() != 1 {
                return Err(invalid_tag(
                    "loss",
                    "I_ERROR",
                    line,
                    format!("[loss] expects 1 arg, got {}", args.len()),
                ));
            }
            Ok(ExtendedInstruction::Loss {
                p: args[0],
                targets,
                line,
            })
        }
        (IError, [t]) if t.name == "correlated_loss" => {
            require_no_params(t, "I_ERROR", line)?;
            if targets.is_empty() || !targets.len().is_multiple_of(2) {
                return Err(invalid_tag(
                    "correlated_loss",
                    "I_ERROR",
                    line,
                    format!(
                        "[correlated_loss] expects a nonzero even target count, got {}",
                        targets.len()
                    ),
                ));
            }
            let ps = match args.len() {
                1 => [args[0], 0.0, 0.0],
                3 => [args[0], args[1], args[2]],
                n => {
                    return Err(invalid_tag(
                        "correlated_loss",
                        "I_ERROR",
                        line,
                        format!("[correlated_loss] expects 1 or 3 args, got {n}"),
                    ));
                }
            };
            Ok(ExtendedInstruction::CorrelatedLoss {
                ps,
                targets: pair_targets(&targets),
                line,
            })
        }
        (IError, []) => Err(invalid_tag(
            "",
            "I_ERROR",
            line,
            "I_ERROR requires a [loss] or [correlated_loss] tag",
        )),
        (IError, [t]) => Err(invalid_tag(
            t.name.clone(),
            "I_ERROR",
            line,
            "expected [loss] or [correlated_loss]",
        )),
        (IError, _) => Err(invalid_tag(
            tags[0].name.clone(),
            "I_ERROR",
            line,
            "expected exactly one tag",
        )),
        _ => Ok(ExtendedInstruction::Noise {
            name,
            tags,
            args,
            targets,
            line,
        }),
    }
}
```

- [ ] **Step 2.6: Update `parse_extended` caller**

In `crates/stim-parser/src/extended/parser.rs:26-29`, change:

```rust
pub fn parse_extended(src: &str) -> Result<ExtendedProgram, ExtendedParseError> {
    let prog = crate::parser::parse(src)?;
    interpret(prog)
}
```

(Drop the `&` before `prog`.)

- [ ] **Step 2.7: Run parser tests**

Run:
```bash
cargo test -p stim-parser
```
Expected: all tests pass.

- [ ] **Step 2.8: Run full workspace tests**

Run:
```bash
cargo test --workspace --exclude ppvm-python-native
```
Expected: all tests pass.

- [ ] **Step 2.9: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/stim-parser/src/extended/interpret.rs crates/stim-parser/src/extended/parser.rs
git commit -m "$(cat <<'EOF'
refactor(stim-parser): consume Program in interpret pass

`interpret` now takes `Program` by value and moves `tags`/`args`/`targets`
into `ExtendedInstruction` instead of cloning per pass-through instruction.
Only caller is `parse_extended`, which constructs and discards `prog`.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Step 3: Restructure `interpret_gate` for S/SDag/Identity

**Goal:** Replace the `if !tags.is_empty()` guard pattern with a tags-shape match. The guard is load-bearing-but-non-obvious; a structural match makes the dispatch read clearly.

**Files:**
- Modify: `crates/stim-parser/src/extended/interpret.rs` (the `interpret_gate` function rewritten in Step 2)

- [ ] **Step 3.1: Replace S/SDag arm shape**

In `interpret_gate` (which Step 2 rewrote), replace the four S/SDag arms and the two Identity arms. The new structure dispatches on `tags.as_slice()` shape instead of using guards.

Replace the body of `interpret_gate` (in `crates/stim-parser/src/extended/interpret.rs`) with:

```rust
fn interpret_gate(
    name: GateName,
    tags: Vec<Tag>,
    args: Vec<f64>,
    targets: Vec<usize>,
    line: usize,
) -> Result<ExtendedInstruction, ExtendedParseError> {
    use GateName::*;

    match name {
        S | SDag => match tags.as_slice() {
            [] => Ok(ExtendedInstruction::Gate {
                name,
                tags,
                args,
                targets,
                line,
            }),
            [t] if t.name == "T" => {
                require_no_params(t, name.canonical_name(), line)?;
                Ok(if matches!(name, S) {
                    ExtendedInstruction::T { targets, line }
                } else {
                    ExtendedInstruction::TDag { targets, line }
                })
            }
            [t] => Err(invalid_tag(
                t.name.clone(),
                name.canonical_name(),
                line,
                "expected [T]",
            )),
            _ => Err(invalid_tag(
                tags[0].name.clone(),
                name.canonical_name(),
                line,
                "expected exactly one tag",
            )),
        },
        Identity => match tags.as_slice() {
            [] => Ok(ExtendedInstruction::Gate {
                name,
                tags,
                args,
                targets,
                line,
            }),
            [t] => interpret_identity_tag(t, targets, line),
            _ => Err(invalid_tag(
                tags[0].name.clone(),
                "I",
                line,
                "expected exactly one tag",
            )),
        },
        _ => Ok(ExtendedInstruction::Gate {
            name,
            tags,
            args,
            targets,
            line,
        }),
    }
}
```

This drops the `(S, [t]) if t.name == "T"`-style mixed guards, the `if !tags.is_empty()` guard, and the stringly-typed `if matches!(name, S) { "S" } else { "S_DAG" }` (replaced by `name.canonical_name()`).

- [ ] **Step 3.2: Run parser tests**

Run:
```bash
cargo test -p stim-parser
```
Expected: all tests pass — the existing extended-dialect tests cover `S[T]`, `S_DAG[T]`, `S[U]` (rejected), `I[R_X(theta=...)]`, and the multi-tag rejection cases.

- [ ] **Step 3.3: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/stim-parser/src/extended/interpret.rs
git commit -m "$(cat <<'EOF'
refactor(stim-parser): restructure interpret_gate dispatch

Match on `tags.as_slice()` shape ([], [t], _) instead of guarded arms.
Drops the `if !tags.is_empty()` guard (load-bearing but non-obvious)
and the stringly-typed `if matches!(name, S) { "S" } else { "S_DAG" }`
in favor of `name.canonical_name()`.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Step 4: `TableEntry` accessors

**Goal:** Add `arg_count`/`target_arity`/`canonical` inherent methods on `TableEntry`. Collapse `arity_of` into call-site uses; keep `build_instruction` (it constructs the typed enum) but simplify around it.

**Files:**
- Modify: `crates/stim-parser/src/table.rs` (add accessors)
- Modify: `crates/stim-parser/src/parser.rs:78-102` (delete `arity_of`)
- Modify: `crates/stim-parser/src/parser.rs:133` (use accessors directly)

- [ ] **Step 4.1: Add accessors to `TableEntry`**

In `crates/stim-parser/src/table.rs`, after the `TableEntry` enum definition (around line 38), add an `impl` block:

```rust
impl TableEntry {
    pub fn arg_count(&self) -> ArgCount {
        match self {
            TableEntry::Gate { args, .. }
            | TableEntry::Noise { args, .. }
            | TableEntry::Measure { args, .. }
            | TableEntry::Annotation { args, .. }
            | TableEntry::MPad { args, .. } => *args,
        }
    }

    pub fn target_arity(&self) -> TargetArity {
        match self {
            TableEntry::Gate { targets, .. }
            | TableEntry::Noise { targets, .. }
            | TableEntry::Measure { targets, .. }
            | TableEntry::Annotation { targets, .. }
            | TableEntry::MPad { targets, .. } => *targets,
        }
    }

    pub fn canonical(&self) -> &'static str {
        match self {
            TableEntry::Gate { name, .. } => name.canonical_name(),
            TableEntry::Noise { name, .. } => name.canonical_name(),
            TableEntry::Measure { name, .. } => name.canonical_name(),
            TableEntry::Annotation { kind, .. } => kind.canonical_name(),
            TableEntry::MPad { .. } => "MPAD",
        }
    }
}
```

- [ ] **Step 4.2: Delete `arity_of` and use accessors at call site**

In `crates/stim-parser/src/parser.rs`, delete the `arity_of` function entirely (lines 78-102).

At line 133, replace:
```rust
            let (arg_rule, target_rule, canonical) = arity_of(entry);
```
with:
```rust
            let arg_rule = entry.arg_count();
            let target_rule = entry.target_arity();
            let canonical = entry.canonical();
```

Everything else in `validate_node` stays the same — it already uses `arg_rule`, `target_rule`, `canonical`.

- [ ] **Step 4.3: Run parser tests**

Run:
```bash
cargo test -p stim-parser
```
Expected: all tests pass. The accessors are also covered transitively by every parser test.

- [ ] **Step 4.4: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/stim-parser/src/table.rs crates/stim-parser/src/parser.rs
git commit -m "$(cat <<'EOF'
refactor(stim-parser): add TableEntry accessors

`arg_count()`, `target_arity()`, `canonical()` replace the standalone
`arity_of` helper. `validate_node` calls each accessor directly.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Step 5: Consolidate `XError`/`YError`/`ZError` arms

**Goal:** Three near-identical match arms in the executor's noise dispatch become one.

**Files:**
- Modify: `crates/ppvm-stim/src/executor.rs:231-251` (the three error arms)

- [ ] **Step 5.1: Replace the three arms with one**

In `crates/ppvm-stim/src/executor.rs`, find the noise match arms for `XError`, `YError`, `ZError` (currently lines 231-251). Replace them with a single combined arm:

```rust
                NoiseName::XError | NoiseName::YError | NoiseName::ZError => {
                    debug_assert_eq!(args.len(), 1);
                    let p: T::Coeff = args[0].into();
                    let zero = T::Coeff::zero();
                    let ps: [T::Coeff; 3] = match name {
                        NoiseName::XError => [p, zero.clone(), zero],
                        NoiseName::YError => [zero.clone(), p, zero],
                        NoiseName::ZError => [zero.clone(), zero, p],
                        _ => unreachable!(),
                    };
                    for &q in targets {
                        tab.pauli_error(q, ps.clone());
                    }
                }
```

- [ ] **Step 5.2: Run ppvm-stim tests**

Run:
```bash
cargo test -p ppvm-stim
```
Expected: all tests pass. The X_ERROR / Y_ERROR / Z_ERROR opcodes are exercised by the executor tests and the corpus harness.

- [ ] **Step 5.3: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/ppvm-stim/src/executor.rs
git commit -m "$(cat <<'EOF'
refactor(ppvm-stim): consolidate XError/YError/ZError arms

Single match arm with name-driven slot selection replaces three
near-identical arms. ~12 lines removed; behavior unchanged.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Step 6: `flip_with_prob` helper for measurement readout

**Goal:** Extract the `if noise > 0.0 && bernoulli(noise) { !b } else { b }` pattern shared by MR and MPad into a single helper.

**Files:**
- Modify: `crates/ppvm-tableau/src/measure.rs` (add `flip_with_prob` method)
- Modify: `crates/ppvm-stim/src/executor.rs:289-326` (MR and MPad arms use it)

- [ ] **Step 6.1: Add `flip_with_prob` to `GeneralizedTableau`**

In `crates/ppvm-tableau/src/measure.rs`, in the `impl<T, I, C> GeneralizedTableau<T, I, C>` block that contains `measure_noisy` and `bernoulli` (around lines 203-249), add:

```rust
    /// Flip `bit` with probability `p`. Used by Stim MR/MPad readout-noise
    /// dispatch in `ppvm-stim`. Returns `bit` unchanged when `p <= 0.0`.
    pub fn flip_with_prob(&mut self, bit: bool, p: f64) -> bool {
        if p > 0.0 && self.bernoulli(p) { !bit } else { bit }
    }
```

Place it adjacent to `bernoulli` so the relationship is obvious.

- [ ] **Step 6.2: Use it from MR**

In `crates/ppvm-stim/src/executor.rs`, replace the MR arm body at lines 289-301 with:

```rust
                    MeasureName::MR => {
                        for &q in targets {
                            let true_outcome = tab.measure(q);
                            if true_outcome == Some(true) {
                                tab.x(q);
                            }
                            let recorded = true_outcome.map(|b| tab.flip_with_prob(b, noise));
                            results.push(recorded);
                        }
                    }
```

- [ ] **Step 6.3: Use it from MPad**

In `crates/ppvm-stim/src/executor.rs`, replace the MPad arm body at lines 316-326 with:

```rust
            ExtendedInstruction::MPad { bits, prob, .. } => {
                let noise = prob.unwrap_or(0.0);
                for &bit in bits {
                    results.push(Some(tab.flip_with_prob(bit, noise)));
                }
            }
```

- [ ] **Step 6.4: Run tests**

Run:
```bash
cargo test -p ppvm-tableau -p ppvm-stim
```
Expected: all tests pass. MR and MPad readout-noise behavior is covered by `crates/ppvm-stim/tests/executor.rs`.

- [ ] **Step 6.5: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/ppvm-tableau/src/measure.rs crates/ppvm-stim/src/executor.rs
git commit -m "$(cat <<'EOF'
refactor: extract flip_with_prob for MR/MPad readout noise

Shared helper on `GeneralizedTableau` replaces the inline
`if noise > 0.0 && bernoulli(noise) { !b } else { b }` pattern in
both MR and MPad executor arms.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Step 7: Investigate `measure_noisy` split

**Goal:** Decide whether `measure_noisy` can be split into `measure` (already exists) + `flip_with_prob` so MZ becomes the same `measure → flip_with_prob` chain that MR uses, removing the divergence noted in commit 212dc41. This is a **decision step** — it may end with no code change.

**Files (potentially):**
- Read: `crates/ppvm-tableau/src/measure.rs:235-242` (`measure_noisy`)
- Read: PR commit 212dc41 message
- Possibly modify: `crates/ppvm-stim/src/executor.rs:279-301` (MZ and MR arms)

- [ ] **Step 7.1: Read `measure_noisy` body**

Run:
```bash
git log --oneline -1 212dc41
git show 212dc41
```
Re-read `crates/ppvm-tableau/src/measure.rs:235-242`. The current body is:

```rust
pub fn measure_noisy(&mut self, addr0: usize, flip_prob: f64) -> Option<bool> {
    let outcome = self.measure(addr0)?;
    if flip_prob > 0.0 && self.tableau.rng.random::<f64>() < flip_prob {
        Some(!outcome)
    } else {
        Some(outcome)
    }
}
```

Note this is `measure` (already returns `Option<bool>`) followed by an inline-bernoulli flip on the recorded bit. The flip is purely readout-side: it does **not** mutate quantum state.

- [ ] **Step 7.2: Decide**

If the flip is purely readout-side (which the source confirms — only the returned value is altered, not `self.coefficients` or `self.tableau`), then `measure_noisy` is equivalent to `measure().map(|b| flip_with_prob(b, p))`. Proceed to step 7.3. If you discover any side effect (e.g., a stabilizer-tracking subtlety), abandon — close out this step with no commit and skip to Step 8.

- [ ] **Step 7.3: Reduce `measure_noisy` to delegate**

In `crates/ppvm-tableau/src/measure.rs`, replace `measure_noisy` body with:

```rust
pub fn measure_noisy(&mut self, addr0: usize, flip_prob: f64) -> Option<bool> {
    let outcome = self.measure(addr0)?;
    Some(self.flip_with_prob(outcome, flip_prob))
}
```

This makes the relationship between MZ, MR, and MPad explicit: all three go through `flip_with_prob`. The MR commit-212dc41 divergence is now narrated structurally — MR cannot delegate to `measure_noisy` because it interleaves a `tab.x(q)` between `measure` and `flip_with_prob`, but it uses the same readout helper.

- [ ] **Step 7.4: Update the MR comment**

In `crates/ppvm-stim/src/executor.rs:284-288`, the comment currently says:

```
// Note: MR cannot delegate to `tab.measure_noisy` like MZ
// because the reset must be conditioned on the *true*
// outcome. The recorded bit gets a separate Bernoulli
// flip — same observable distribution as MZ, but a
// distinct RNG draw.
```

Update to:
```
// MR cannot delegate to `measure_noisy` because the reset must use
// the *true* outcome — but it shares the readout flip via
// `flip_with_prob`, so the RNG-draw shape matches MZ exactly.
```

- [ ] **Step 7.5: Run tests**

Run:
```bash
cargo test -p ppvm-tableau -p ppvm-stim
```
Expected: all tests pass. MZ behavior with non-zero readout-noise probability is covered by the executor and corpus tests.

- [ ] **Step 7.6: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/ppvm-tableau/src/measure.rs crates/ppvm-stim/src/executor.rs
git commit -m "$(cat <<'EOF'
refactor(tableau): measure_noisy delegates to flip_with_prob

MZ and MR now share the same readout-flip helper, narrating the
divergence (MR's reset interleaves `tab.x(q)`) more clearly.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

If Step 7.2 ended with "abandon", skip Steps 7.3–7.6 and move to Step 8.

---

## Step 8: `PyStimProgram` newtype

**Goal:** Replace the named `program: ExtendedProgram` field with a tuple newtype + `Deref<Target=ExtendedProgram>`. Call sites lose one `.program` indirection.

**Files:**
- Modify: `crates/ppvm-python-native/src/stim_program.rs` (newtype + Deref)
- Modify: `crates/ppvm-python-native/src/interface_tableau.rs:170-212` (lose `.program` indirection)

- [ ] **Step 8.1: Convert `PyStimProgram` to a tuple struct with `Deref`**

Replace `crates/ppvm-python-native/src/stim_program.rs` entirely with:

```rust
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::ops::Deref;

use ppvm_stim::{ExtendedProgram, parse_extended, prepare};

/// Python-facing wrapper around a prepared extended Stim program.
#[pyclass(name = "StimProgram", module = "ppvm_python_native")]
pub struct PyStimProgram(pub ExtendedProgram);

#[pymethods]
impl PyStimProgram {
    /// Parse and prepare a Stim circuit string.
    #[staticmethod]
    pub fn parse(src: &str) -> PyResult<Self> {
        let program = parse_extended(src).map_err(stim_to_pyerr)?;
        prepare(&program).map_err(stim_to_pyerr_exec)?;
        Ok(Self(program))
    }

    /// Read a `.stim` file and parse it.
    #[staticmethod]
    pub fn from_file(path: &str) -> PyResult<Self> {
        let src = std::fs::read_to_string(path)
            .map_err(|e| PyIOError::new_err(format!("failed to read {path}: {e}")))?;
        Self::parse(&src)
    }

    fn __repr__(&self) -> String {
        format!(
            "<StimProgram instructions={} measurements={}>",
            self.0.instructions.len(),
            self.0.measurement_count()
        )
    }
}

impl Deref for PyStimProgram {
    type Target = ExtendedProgram;
    fn deref(&self) -> &ExtendedProgram {
        &self.0
    }
}

fn stim_to_pyerr(e: ppvm_stim::ExtendedParseError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}

fn stim_to_pyerr_exec(e: ppvm_stim::ExecError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}
```

- [ ] **Step 8.2: Update interface_tableau.rs `run` and `sample` callers**

In `crates/ppvm-python-native/src/interface_tableau.rs`, change `run`:

```rust
            // STIM integration — runs a parsed and prepared StimProgram.
            pub fn run(
                &mut self,
                prog: &crate::stim_program::PyStimProgram,
            ) -> pyo3::PyResult<Vec<Option<bool>>> {
                let mut results = Vec::with_capacity(prog.measurement_count());
                ppvm_stim::execute_prepared(
                    &prog.instructions,
                    &mut self.inner,
                    &mut results,
                );
                Ok(results)
            }
```

(Replaced `prog.program.measurement_count()` → `prog.measurement_count()` and `prog.program.instructions` → `prog.instructions`. `Deref` provides both.)

In `sample` at line 187, change `&prog.program` to `&prog`:

```rust
                ppvm_stim::sample(&prog, num_shots, || {
```

(Note: `&PyStimProgram` derefs to `&ExtendedProgram`, which is what `ppvm_stim::sample` expects.)

Also delete the change-narrating `// STIM integration —` comment on the `run` method, and the `// some python niceties` comment higher up — they're noise. Remove the dead 14-line commented-out `__richcmp__` block at lines 241-254 (it's not part of this PR's value-add and just got pulled along).

- [ ] **Step 8.3: Build the Python native crate**

Run:
```bash
cargo build -p ppvm-python-native
```
Expected: clean build.

- [ ] **Step 8.4: Run workspace tests**

Run:
```bash
cargo test --workspace --exclude ppvm-python-native
```
Expected: all tests pass.

- [ ] **Step 8.5: Run the Python test suite**

This requires the maturin-built native module. Per `memory/reference_python_native_build.md`:

```bash
cd ppvm-python
uv sync --reinstall-package ppvm-python-native
uv run pytest
```
Expected: all Python tests pass.

- [ ] **Step 8.6: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/ppvm-python-native/src/stim_program.rs crates/ppvm-python-native/src/interface_tableau.rs
git commit -m "$(cat <<'EOF'
refactor(python-native): make PyStimProgram a Deref newtype

Tuple struct + `Deref<Target=ExtendedProgram>` replaces the named
`program:` field. Call sites lose one indirection. Also drops dead
`__richcmp__` comment block and change-narrating comments.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Step 9: `MeasurementResult` mapping at the PyO3 boundary

**Goal:** PyO3 binding returns `Vec<u8>` (0=ZERO, 1=ONE, 2=LOST). Python wrapper uses `[MeasurementResult(x) for x in raw]` instead of `_from_raw` per element.

**Files:**
- Modify: `crates/ppvm-python-native/src/interface_tableau.rs:170-212` (run/sample return Vec<u8>)
- Modify: `crates/ppvm-python-native/ppvm_python_native.pyi:143-151` (typing)
- Modify: `ppvm-python/src/ppvm/generalized_tableau.py:25-31, 215, 234` (drop `_from_raw`)

- [ ] **Step 9.1: Add a small mapping helper in PyO3**

In `crates/ppvm-python-native/src/interface_tableau.rs`, near the top of the file (after the imports, before the macro), add:

```rust
fn measurement_to_u8(m: Option<bool>) -> u8 {
    match m {
        Some(false) => 0,
        Some(true) => 1,
        None => 2,
    }
}
```

(The values match `ppvm.MeasurementResult.ZERO/ONE/LOST` — keep them in sync with the Python enum at `ppvm-python/src/ppvm/generalized_tableau.py:18-23`.)

- [ ] **Step 9.2: Change `run` return type to `Vec<u8>`**

In the same file, change `run`. Since `run` is generated inside `macro_rules! create_interface!` (expanded at top level of `interface_tableau.rs`), use the explicit `crate::interface_tableau::measurement_to_u8` path — that's the canonical reference and avoids any `macro_rules!` hygiene surprises:

```rust
            pub fn run(
                &mut self,
                prog: &crate::stim_program::PyStimProgram,
            ) -> pyo3::PyResult<Vec<u8>> {
                let mut results = Vec::with_capacity(prog.measurement_count());
                ppvm_stim::execute_prepared(
                    &prog.instructions,
                    &mut self.inner,
                    &mut results,
                );
                Ok(results.into_iter().map(crate::interface_tableau::measurement_to_u8).collect())
            }
```

- [ ] **Step 9.3: Change `sample` return type to `Vec<Vec<u8>>`**

In the same file, update `sample`:

```rust
            #[staticmethod]
            #[pyo3(signature = (prog, n_qubits, min_abs_coeff = 1e-10, num_shots = 1, seed = None))]
            pub fn sample(
                prog: &crate::stim_program::PyStimProgram,
                n_qubits: usize,
                min_abs_coeff: f64,
                num_shots: usize,
                seed: Option<u64>,
            ) -> pyo3::PyResult<Vec<Vec<u8>>> {
                let mut next_seed = seed;
                let raw = ppvm_stim::sample(&prog, num_shots, || {
                    let s = next_seed;
                    if let Some(ref mut v) = next_seed {
                        *v = v.wrapping_add(1);
                    }
                    match s {
                        Some(s) => GeneralizedTableau::<$type, $indexType>::new_with_seed(
                            n_qubits,
                            min_abs_coeff,
                            s,
                        ),
                        None => {
                            GeneralizedTableau::<$type, $indexType>::new(n_qubits, min_abs_coeff)
                        }
                    }
                })
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
                Ok(raw
                    .into_iter()
                    .map(|shot| shot.into_iter().map(crate::interface_tableau::measurement_to_u8).collect())
                    .collect())
            }
```

(The `?` on `.map_err(...)` propagates the error before the `Ok(...)` mapping.)

- [ ] **Step 9.4: Update the .pyi stub**

In `crates/ppvm-python-native/ppvm_python_native.pyi`, change lines 143-151:

```python
    def run(self, prog: StimProgram) -> list[int]: ...
    @staticmethod
    def sample(
        prog: StimProgram,
        n_qubits: int,
        min_abs_coeff: float = 1e-10,
        num_shots: int = 1,
        seed: int | None = None,
    ) -> list[list[int]]: ...
```

(Quotes around `StimProgram` are no longer necessary in this position — Python 3 PEP 563 / `from __future__ import annotations` would resolve them lazily, but since `StimProgram` is declared later in the same file, leave them as `"StimProgram"` if the existing style does; otherwise drop them. Match the rest of the file.)

- [ ] **Step 9.5: Update Python wrapper to construct `MeasurementResult` directly**

In `ppvm-python/src/ppvm/generalized_tableau.py`:

Delete the `_from_raw` staticmethod at lines 25-31. The IntEnum constructor `MeasurementResult(x)` for `x ∈ {0, 1, 2}` returns the same value.

Update `measure` at lines 139-150 — `measure` still returns `bool | None` from the underlying interface (it's a single-qubit measurement, not a Stim run), so this method is unaffected. Leave its body, but inline the conversion:

```python
    def measure(self, addr0: int) -> MeasurementResult:
        """Measure the specified qubit in the Z basis. ..."""
        m = self._interface.measure(addr0)
        if m is None:
            return MeasurementResult.LOST
        return MeasurementResult.ONE if m else MeasurementResult.ZERO
```

Update `run` at lines 206-215:

```python
    def run(self, prog: StimProgram) -> list[MeasurementResult]:
        """Execute a parsed Stim program against this tableau (single shot).

        .. note::
            This **mutates** the tableau in place. For independent shots use
            :meth:`fork` or the :func:`ppvm.sample_stim` / :meth:`sample`
            helpers (which build a fresh tableau per shot).
        """
        raw = self._interface.run(prog)
        return [MeasurementResult(x) for x in raw]
```

Update `sample` classmethod at lines 217-234:

```python
    @classmethod
    def sample(
        cls,
        prog: StimProgram,
        n_qubits: int,
        min_abs_coeff: float = 1e-10,
        num_shots: int = 1,
        seed: int | None = None,
    ) -> list[list[MeasurementResult]]:
        """Run ``num_shots`` shots of ``prog`` and return all measurement results.

        Each shot starts from a fresh tableau, so this is the right entry
        point for multi-shot sampling.
        """
        N_interface = math.ceil(n_qubits / 64.0)
        native_cls = getattr(ppvm_python_native, f"GeneralizedTableau{N_interface}")
        raw = native_cls.sample(prog, n_qubits, min_abs_coeff, num_shots, seed)
        return [[MeasurementResult(x) for x in shot] for shot in raw]
```

- [ ] **Step 9.6: Build the Python native crate**

Run:
```bash
cargo build -p ppvm-python-native
```
Expected: clean build.

- [ ] **Step 9.7: Rebuild and run Python tests**

Run:
```bash
cd ppvm-python
uv sync --reinstall-package ppvm-python-native
uv run pytest
```
Expected: all Python tests pass. The tests in `test_stim.py` (run/sample) and `test_loss.py` (LOST measurements) cover the conversion.

- [ ] **Step 9.8: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/ppvm-python-native/src/interface_tableau.rs crates/ppvm-python-native/ppvm_python_native.pyi ppvm-python/src/ppvm/generalized_tableau.py
git commit -m "$(cat <<'EOF'
refactor: map MeasurementResult at PyO3 boundary

`tab.run` / `cls.sample` return `Vec<u8>` enum-coded values (0/1/2)
instead of `Vec<Option<bool>>`. Python wrapper batch-converts via
`MeasurementResult(x)` (IntEnum constructor, C-speed) instead of
the `_from_raw` per-element 3-way branch.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Step 10: Restore lost gate test coverage

**Goal:** Recreate `crates/ppvm-tableau/tests/gates.rs` with only the missing-coverage subset: cross-axis 2q rotations and direct stabilizer-row assertions for basic Cliffords. Skip composition identities (redundant) and the lost-qubit/Stim-parser sections (covered elsewhere now).

**Files:**
- Create: `crates/ppvm-tableau/tests/gates.rs`

The original file is at `git show main:crates/ppvm-tableau/tests/gates.rs` — copy specific tests from there verbatim.

- [ ] **Step 10.1: Recreate the test file with the targeted subset**

Create `crates/ppvm-tableau/tests/gates.rs` with the following content (this is a strict subset of the deleted file — the targeted tests from sections 1, 2 (basic), and 4 (cross-axis rotations) only):

```rust
//! Targeted gate tests restored from the pre-PR-69 ppvm-tableau test suite.
//!
//! Covers two areas not exercised by the surviving src-file tests:
//! - Tableau direct Clifford row assertions (X, Y, Z, H, S, S†, CNOT, CZ).
//! - Cross-axis 2q rotations (rxy, rxz, ryx, ryz, rzx, rzy) and rzz on
//!   computational basis states.

use ppvm_runtime::config::fxhash::ByteF64;
use ppvm_tableau::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI};

type TC = ByteF64<1>;
type Tab = Tableau<TC>;
type GTab = GeneralizedTableau<TC>;

fn stab1(t: &Tab) -> String {
    t.stabilizers()[0].to_string()
}

fn destab1(t: &Tab) -> String {
    t.destabilizers()[0].to_string()
}

fn stab2(t: &Tab) -> (String, String) {
    (
        t.stabilizers()[0].to_string(),
        t.stabilizers()[1].to_string(),
    )
}

fn destab2(t: &Tab) -> (String, String) {
    (
        t.destabilizers()[0].to_string(),
        t.destabilizers()[1].to_string(),
    )
}

// ============================================================
// Tableau direct Clifford gate tests
// ============================================================

#[test]
fn test_tableau_x_gate() {
    let mut t: Tab = Tableau::new(1);
    t.x(0);
    assert_eq!(stab1(&t), "-Z");
    assert_eq!(destab1(&t), "+X");
}

#[test]
fn test_tableau_y_gate() {
    let mut t: Tab = Tableau::new(1);
    t.y(0);
    assert_eq!(stab1(&t), "-Z");
    assert_eq!(destab1(&t), "-X");
}

#[test]
fn test_tableau_z_gate() {
    let mut t: Tab = Tableau::new(1);
    t.z(0);
    assert_eq!(stab1(&t), "+Z");
    assert_eq!(destab1(&t), "-X");
}

#[test]
fn test_tableau_h_gate() {
    let mut t: Tab = Tableau::new(1);
    t.h(0);
    assert_eq!(stab1(&t), "+X");
    assert_eq!(destab1(&t), "+Z");
}

#[test]
fn test_tableau_s_gate() {
    let mut t: Tab = Tableau::new(1);
    t.s(0);
    assert_eq!(stab1(&t), "+Z");
    assert_eq!(destab1(&t), "+Y");
}

#[test]
fn test_tableau_s_adj_gate() {
    let mut t: Tab = Tableau::new(1);
    t.s_adj(0);
    assert_eq!(stab1(&t), "+Z");
    assert_eq!(destab1(&t), "-Y");
}

#[test]
fn test_tableau_cnot_on_00() {
    let mut t: Tab = Tableau::new(2);
    t.cnot(0, 1);
    assert_eq!(stab2(&t), ("+ZI".to_string(), "+ZZ".to_string()));
    assert_eq!(destab2(&t), ("+XX".to_string(), "+IX".to_string()));
}

#[test]
fn test_tableau_cz_on_00() {
    let mut t: Tab = Tableau::new(2);
    t.cz(0, 1);
    assert_eq!(stab2(&t), ("+ZI".to_string(), "+IZ".to_string()));
    assert_eq!(destab2(&t), ("+XZ".to_string(), "+ZX".to_string()));
}

#[test]
fn test_tableau_bell_state_via_h_cnot() {
    let mut t: Tab = Tableau::new(2);
    t.h(0);
    t.cnot(0, 1);
    assert_eq!(stab2(&t), ("+XX".to_string(), "+ZZ".to_string()));
}

// ============================================================
// Cross-axis 2q rotation tests
// ============================================================

/// rxy(π) on |00⟩: XY|00⟩ = X|0⟩ ⊗ Y|0⟩ = |1⟩ ⊗ (i|1⟩) = i|11⟩.
/// So rxy(π)|00⟩ = |11⟩ (deterministic, no branching).
#[test]
fn test_rxy_pi_flips_both() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.rxy(0, 1, PI);
    assert_eq!(g.coefficients.len(), 1, "rxy(π) should not branch");
    assert!(g.measure(0).unwrap());
    assert!(g.measure(1).unwrap());
}

#[test]
fn test_rxy_half_pi_branches() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.rxy(0, 1, FRAC_PI_2);
    assert_eq!(g.coefficients.len(), 2, "rxy(π/2) should create 2 branches");
}

/// rxz(π) on |00⟩: XZ|00⟩ = X|0⟩ ⊗ Z|0⟩ = |1⟩ ⊗ |0⟩ = |10⟩.
/// So rxz(π)|00⟩ = -i|10⟩.
#[test]
fn test_rxz_pi_flips_first() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.rxz(0, 1, PI);
    assert_eq!(g.coefficients.len(), 1, "rxz(π) should not branch");
    assert!(g.measure(0).unwrap());
    assert!(!g.measure(1).unwrap());
}

#[test]
fn test_rxz_half_pi_branches() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.rxz(0, 1, FRAC_PI_2);
    assert_eq!(g.coefficients.len(), 2, "rxz(π/2) should create 2 branches");
}

/// ryx(π) on |00⟩: YX|00⟩ = Y|0⟩ ⊗ X|0⟩ = (i|1⟩) ⊗ |1⟩ = i|11⟩.
/// So ryx(π)|00⟩ = |11⟩.
#[test]
fn test_ryx_pi_flips_both() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.ryx(0, 1, PI);
    assert_eq!(g.coefficients.len(), 1, "ryx(π) should not branch");
    assert!(g.measure(0).unwrap());
    assert!(g.measure(1).unwrap());
}

#[test]
fn test_ryx_half_pi_branches() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.ryx(0, 1, FRAC_PI_2);
    assert_eq!(g.coefficients.len(), 2, "ryx(π/2) should create 2 branches");
}

/// ryz(π) on |00⟩: YZ|00⟩ = Y|0⟩ ⊗ Z|0⟩ = (i|1⟩) ⊗ |0⟩ = i|10⟩.
/// So ryz(π)|00⟩ = |10⟩.
#[test]
fn test_ryz_pi_flips_first() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.ryz(0, 1, PI);
    assert_eq!(g.coefficients.len(), 1, "ryz(π) should not branch");
    assert!(g.measure(0).unwrap());
    assert!(!g.measure(1).unwrap());
}

#[test]
fn test_ryz_half_pi_branches() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.ryz(0, 1, FRAC_PI_2);
    assert_eq!(g.coefficients.len(), 2, "ryz(π/2) should create 2 branches");
}

/// rzx(π) on |00⟩: ZX|00⟩ = Z|0⟩ ⊗ X|0⟩ = |0⟩ ⊗ |1⟩ = |01⟩.
/// So rzx(π)|00⟩ = -i|01⟩.
#[test]
fn test_rzx_pi_flips_second() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.rzx(0, 1, PI);
    assert_eq!(g.coefficients.len(), 1, "rzx(π) should not branch");
    assert!(!g.measure(0).unwrap());
    assert!(g.measure(1).unwrap());
}

#[test]
fn test_rzx_half_pi_branches() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.rzx(0, 1, FRAC_PI_2);
    assert_eq!(g.coefficients.len(), 2, "rzx(π/2) should create 2 branches");
}

/// rzy(π) on |00⟩: ZY|00⟩ = Z|0⟩ ⊗ Y|0⟩ = |0⟩ ⊗ (i|1⟩) = i|01⟩.
/// So rzy(π)|00⟩ = |01⟩.
#[test]
fn test_rzy_pi_flips_second() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.rzy(0, 1, PI);
    assert_eq!(g.coefficients.len(), 1, "rzy(π) should not branch");
    assert!(!g.measure(0).unwrap());
    assert!(g.measure(1).unwrap());
}

#[test]
fn test_rzy_half_pi_branches() {
    let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
    g.rzy(0, 1, FRAC_PI_2);
    assert_eq!(g.coefficients.len(), 2, "rzy(π/2) should create 2 branches");
}

/// rzz on computational basis never branches (ZZ is diagonal in Z basis).
#[test]
fn test_rzz_never_branches_on_comp_basis() {
    for state in [(false, false), (true, false), (false, true), (true, true)] {
        let mut g: GTab = GeneralizedTableau::new(2, 1e-12);
        if state.0 {
            g.x(0);
        }
        if state.1 {
            g.x(1);
        }
        g.rzz(0, 1, 0.7);
        assert_eq!(
            g.coefficients.len(),
            1,
            "rzz should not branch on |{}{}⟩",
            state.0 as u8,
            state.1 as u8
        );
    }
}
```

- [ ] **Step 10.2: Run the new tests**

Run:
```bash
cargo test -p ppvm-tableau --test gates
```
Expected: all 22 tests pass. If any fail, double-check that `Tableau::new(1)` and the `stabilizers()`/`destabilizers()` accessors are still public and behave the same — they're load-bearing for Step 10.

- [ ] **Step 10.3: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/ppvm-tableau/tests/gates.rs
git commit -m "$(cat <<'EOF'
test(ppvm-tableau): restore Clifford row + cross-axis rot2 coverage

Targeted subset of the pre-PR-69 gates.rs: direct stabilizer/destabilizer
row assertions for basic Cliffords (X, Y, Z, H, S, S†, CNOT, CZ) and
cross-axis 2q rotations (rxy, rxz, ryx, ryz, rzx, rzy, rzz). Skips the
composition identities (redundant with row assertions) and the
lost-qubit/Stim-parser sections (covered by gates/clifford.rs and
ppvm-stim's tests respectively).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Step 11: `ExtendedInstruction::Raw(RawInstruction)` — last, with bail

**Goal:** Replace the six pass-through variants of `ExtendedInstruction` with a single `Raw(RawInstruction)` variant. This consolidates the vanilla-Stim AST shape into one source of truth.

**Bail criterion:** After applying the change locally, scan `crates/ppvm-stim/src/executor.rs` and `crates/ppvm-stim/src/prepare.rs`. If the dispatch reads meaningfully harder to follow — not "one wrapper noisier" but "I have to think to find the gate arm" — abort: `git reset --hard HEAD~1` (only the pre-Step-11 commit) and skip to "Final verification."

**Files:**
- Modify: `crates/stim-parser/src/extended/ast.rs` (collapse variants)
- Modify: `crates/stim-parser/src/extended/interpret.rs` (wrap pass-through arms)
- Modify: `crates/ppvm-stim/src/executor.rs` (match on `Raw(...)`)
- Modify: `crates/ppvm-stim/src/prepare.rs` (match on `Raw(...)`)

**Design note:** The collapse covers `Gate`, `Noise`, `Measure`, `Annotation` — these pass through unchanged from `RawInstruction`. `MPad` keeps its own variant because the extended dialect validates its bits (`Vec<usize>` → `Vec<bool>`). `Repeat` keeps its own variant because its body type changes (`Vec<RawInstruction>` → `Vec<ExtendedInstruction>`) — folding it into `Raw` would force the executor to dispatch over both raw and extended bodies, which kills the abstraction. So `Raw` covers exactly the four pass-through-without-modification cases.

- [ ] **Step 11.1: Replace `ExtendedInstruction` definition**

In `crates/stim-parser/src/extended/ast.rs`, replace the entire enum (lines 11-82) with:

```rust
#[derive(Debug, Clone, PartialEq)]
pub enum ExtendedInstruction {
    /// Pass-through from vanilla Stim — covers `Gate`, `Noise`, `Measure`,
    /// `Annotation`. `MPad` and `Repeat` are NOT in `Raw` because their
    /// bits/body shapes diverge between dialects.
    Raw(RawInstruction),

    // --- Extended-dialect sugar variants ---
    T {
        targets: Vec<usize>,
        line: usize,
    },
    TDag {
        targets: Vec<usize>,
        line: usize,
    },
    Rotation {
        axis: Axis,
        theta: f64,
        targets: Vec<usize>,
        line: usize,
    },
    U3 {
        theta: f64,
        phi: f64,
        lambda: f64,
        targets: Vec<usize>,
        line: usize,
    },
    Loss {
        p: f64,
        targets: Vec<usize>,
        line: usize,
    },
    CorrelatedLoss {
        ps: [f64; 3],
        targets: Vec<(usize, usize)>,
        line: usize,
    },
    /// Extended-dialect MPad with validated bits (each ∈ {0, 1}).
    MPad {
        tags: Vec<Tag>,
        prob: Option<f64>,
        bits: Vec<bool>,
        line: usize,
    },
    Repeat {
        count: u64,
        body: Vec<ExtendedInstruction>,
        line: usize,
    },
}
```

Replace the import line at the top of the file with:

```rust
use crate::ast::{RawInstruction, Tag};
```

(`AnnotationKind`, `GateName`, `MeasureName`, `NoiseName` are no longer pattern-matched directly in this file — they reach `RawInstruction` arms via `Raw(_)`.)

- [ ] **Step 11.2: Update `count_in_slice` for the new shape**

In `crates/stim-parser/src/extended/ast.rs`, replace `count_in_slice` (lines 99-124) with:

```rust
fn count_in_slice(instructions: &[ExtendedInstruction], factor: u64) -> usize {
    let mut total = 0usize;
    for instr in instructions {
        match instr {
            ExtendedInstruction::Raw(RawInstruction::Measure { targets, .. }) => {
                total = total.saturating_add(targets.len().saturating_mul(factor as usize));
            }
            ExtendedInstruction::MPad { bits, .. } => {
                total = total.saturating_add(bits.len().saturating_mul(factor as usize));
            }
            ExtendedInstruction::Repeat { count, body, .. } => {
                total = total.saturating_add(count_in_slice(body, factor.saturating_mul(*count)));
            }
            ExtendedInstruction::Raw(_)
            | ExtendedInstruction::T { .. }
            | ExtendedInstruction::TDag { .. }
            | ExtendedInstruction::Rotation { .. }
            | ExtendedInstruction::U3 { .. }
            | ExtendedInstruction::Loss { .. }
            | ExtendedInstruction::CorrelatedLoss { .. } => {}
        }
    }
    total
}
```

Note: the `Raw(_)` wildcard arm sweeps up `Raw(Gate)`, `Raw(Noise)`, `Raw(Annotation)`, and (defensively) `Raw(MPad)` / `Raw(Repeat)` — none of which produce measurement bits. The interpret pass guarantees the latter two never appear, but matching them as `Raw(_) => {}` keeps the match exhaustive without an extra arm.

- [ ] **Step 11.3: Update `interpret_one` to wrap pass-through arms**

In `crates/stim-parser/src/extended/interpret.rs`, replace `interpret_one` to wrap `Measure` and `Annotation` (and inner `Repeat` body recursion stays the same), but `Gate` and `Noise` go through `interpret_gate`/`interpret_noise` which return either a sugar variant or `Raw`.

```rust
fn interpret_one(raw: RawInstruction) -> Result<ExtendedInstruction, ExtendedParseError> {
    match raw {
        RawInstruction::Gate {
            name,
            tags,
            args,
            targets,
            line,
        } => interpret_gate(name, tags, args, targets, line),
        RawInstruction::Noise {
            name,
            tags,
            args,
            targets,
            line,
        } => interpret_noise(name, tags, args, targets, line),
        m @ RawInstruction::Measure { .. } => Ok(ExtendedInstruction::Raw(m)),
        a @ RawInstruction::Annotation { .. } => Ok(ExtendedInstruction::Raw(a)),
        RawInstruction::MPad {
            tags,
            prob,
            bits,
            line,
        } => Ok(ExtendedInstruction::MPad {
            tags,
            prob,
            bits: convert_mpad_bits(&bits, line)?,
            line,
        }),
        RawInstruction::Repeat { count, body, line } => {
            let mut inner = Vec::with_capacity(body.len());
            interpret_slice(body, &mut inner)?;
            Ok(ExtendedInstruction::Repeat {
                count,
                body: inner,
                line,
            })
        }
    }
}
```

In `interpret_gate`, the default (no-tag, non-S/SDag/Identity) arm changes from `ExtendedInstruction::Gate { ... }` to:

```rust
        _ => Ok(ExtendedInstruction::Raw(RawInstruction::Gate {
            name,
            tags,
            args,
            targets,
            line,
        })),
```

Same for the `S | SDag` empty-tags arm and `Identity` empty-tags arm — they all now wrap.

In `interpret_noise`, the default arm becomes:

```rust
        _ => Ok(ExtendedInstruction::Raw(RawInstruction::Noise {
            name,
            tags,
            args,
            targets,
            line,
        })),
```

- [ ] **Step 11.4: Update `executor.rs` dispatch**

In `crates/ppvm-stim/src/executor.rs`, the top-level match now handles `Raw(RawInstruction::Gate { .. })`, `Raw(RawInstruction::Noise { .. })`, `Raw(RawInstruction::Measure { .. })`, and `Raw(RawInstruction::Annotation { .. })`. The sugar variants stay flat.

Replace the top-level match (around line 127) with:

```rust
    for instr in instructions {
        match instr {
            ExtendedInstruction::Raw(RawInstruction::Gate { name, targets, .. }) => match name {
                // ... gate dispatch arms unchanged ...
            },
            ExtendedInstruction::T { targets, .. } => targets.iter().for_each(|&q| tab.t(q)),
            ExtendedInstruction::TDag { targets, .. } => {
                targets.iter().for_each(|&q| tab.t_adj(q));
            }
            ExtendedInstruction::Rotation { axis, theta, targets, .. } => match axis {
                Axis::X => targets.iter().for_each(|&q| tab.rx(q, *theta)),
                Axis::Y => targets.iter().for_each(|&q| tab.ry(q, *theta)),
                Axis::Z => targets.iter().for_each(|&q| tab.rz(q, *theta)),
            },
            ExtendedInstruction::U3 { theta, phi, lambda, targets, .. } => targets
                .iter()
                .for_each(|&q| tab.u3(q, (*theta).into(), (*phi).into(), (*lambda).into())),
            ExtendedInstruction::Raw(RawInstruction::Noise { name, targets, args, .. }) => match name {
                // ... noise dispatch arms unchanged ...
            },
            ExtendedInstruction::Loss { p, targets, .. } => {
                for &q in targets {
                    tab.loss_channel(q, (*p).into());
                }
            }
            ExtendedInstruction::CorrelatedLoss { ps, targets, .. } => {
                let ps: [T::Coeff; 3] = [ps[0].into(), ps[1].into(), ps[2].into()];
                for &(a, b) in targets {
                    tab.correlated_loss_channel(a, b, ps.clone());
                }
            }
            ExtendedInstruction::Raw(RawInstruction::Measure { name, args, targets, .. }) => {
                let noise = args.first().copied().unwrap_or(0.0);
                match name {
                    // ... measure arms unchanged ...
                }
            }
            ExtendedInstruction::MPad { bits, prob, .. } => {
                let noise = prob.unwrap_or(0.0);
                for &bit in bits {
                    results.push(Some(tab.flip_with_prob(bit, noise)));
                }
            }
            ExtendedInstruction::Raw(RawInstruction::Annotation { .. }) => { /* no-op */ }
            ExtendedInstruction::Raw(RawInstruction::MPad { .. }) => {
                unreachable!("MPad is consumed into the extended dialect; never reaches Raw");
            }
            ExtendedInstruction::Raw(RawInstruction::Repeat { .. }) => {
                unreachable!("Repeat is consumed into the extended dialect; never reaches Raw");
            }
            ExtendedInstruction::Repeat { count, body, .. } => {
                for _ in 0..*count {
                    execute_prepared(body, tab, results);
                }
            }
        }
    }
```

The two `unreachable!` arms guard the invariant that `MPad` and `Repeat` always go through their dedicated extended variants — the interpret pass never produces `Raw(RawInstruction::MPad)` or `Raw(RawInstruction::Repeat)`. Keep those guards for type-system completeness.

Add `use stim_parser::ast::RawInstruction;` to the imports near the top of `executor.rs`.

- [ ] **Step 11.5: Update `prepare.rs` dispatch**

In `crates/ppvm-stim/src/prepare.rs`, replace `validate_slice` (around lines 14-34) with:

```rust
fn validate_slice(instructions: &[ExtendedInstruction]) -> Result<(), ExecError> {
    for instr in instructions {
        match instr {
            ExtendedInstruction::Raw(RawInstruction::Gate { name, line, .. }) => {
                check_gate_supported(*name, *line)?;
            }
            ExtendedInstruction::Raw(RawInstruction::Noise { name, line, .. }) => {
                check_noise_supported(*name, *line)?;
            }
            ExtendedInstruction::Raw(RawInstruction::Measure { name, line, .. }) => {
                check_measure_supported(*name, *line)?;
            }
            ExtendedInstruction::Repeat { body, .. } => validate_slice(body)?,
            ExtendedInstruction::Raw(RawInstruction::Annotation { .. })
            | ExtendedInstruction::MPad { .. }
            | ExtendedInstruction::T { .. }
            | ExtendedInstruction::TDag { .. }
            | ExtendedInstruction::Rotation { .. }
            | ExtendedInstruction::U3 { .. }
            | ExtendedInstruction::Loss { .. }
            | ExtendedInstruction::CorrelatedLoss { .. } => {}
            ExtendedInstruction::Raw(RawInstruction::MPad { .. }) => {
                unreachable!("MPad never appears as Raw — interpret pass consumes it");
            }
            ExtendedInstruction::Raw(RawInstruction::Repeat { .. }) => {
                unreachable!("Repeat never appears as Raw — interpret pass consumes it");
            }
        }
    }
    Ok(())
}
```

Add `use stim_parser::ast::RawInstruction;` to the imports near the top of `prepare.rs`.

- [ ] **Step 11.6: Run all tests**

Run:
```bash
cargo test --workspace --exclude ppvm-python-native
```
Expected: all tests pass. The parser, executor, prepare, and corpus tests all exercise the new shape.

- [ ] **Step 11.7: Bail check — read the executor + prepare diffs**

Run:
```bash
git diff --stat HEAD
```
Then read `crates/ppvm-stim/src/executor.rs` start-to-finish (the dispatch top-level match) and `crates/ppvm-stim/src/prepare.rs::validate_slice`. Ask:

- Is the gate arm easy to find from a cold read?
- Does the `Raw(RawInstruction::Foo { .. })` wrapper feel cosmetic or does it obscure intent?
- Are the two `unreachable!` arms (Raw MPad / Raw Repeat) carrying weight, or do they read as junk?

**If the dispatch reads meaningfully harder to follow:** abandon. Run:
```bash
git restore --staged .
git restore .
```
Then move on to the final verification step. The cost (per-site noise, two unreachable! arms) was real and didn't pay off. This is an acceptable outcome — the bail criterion existed for exactly this reason.

**If the dispatch reads about the same or improved:** proceed to Step 11.8.

- [ ] **Step 11.8: Run Python tests**

Run:
```bash
cd ppvm-python
uv sync --reinstall-package ppvm-python-native
uv run pytest
```
Expected: all Python tests pass.

- [ ] **Step 11.9: Format and commit**

Run:
```bash
cargo fmt --all
git add crates/stim-parser/src/extended/ast.rs crates/stim-parser/src/extended/interpret.rs crates/ppvm-stim/src/executor.rs crates/ppvm-stim/src/prepare.rs
git commit -m "$(cat <<'EOF'
refactor(stim-parser): collapse pass-through ExtendedInstruction variants

ExtendedInstruction::{Gate, Noise, Measure, Annotation} fold into a
single `Raw(RawInstruction)` variant — single source of truth for
vanilla Stim shapes. MPad and Repeat keep dedicated variants because
their bodies diverge between dialects (validated bits, nested
ExtendedInstruction body).

Two `unreachable!` arms guard the invariant that interpret never
produces `Raw(MPad)` or `Raw(Repeat)`.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Final verification

- [ ] **Step 12.1: Full workspace test run**

Run:
```bash
cargo test --workspace --exclude ppvm-python-native
```
Expected: all tests pass.

- [ ] **Step 12.2: Build Python native crate**

Run:
```bash
cargo build -p ppvm-python-native
```
Expected: clean build.

- [ ] **Step 12.3: Run Python tests**

Run:
```bash
cd ppvm-python
uv sync --reinstall-package ppvm-python-native
uv run pytest
```
Expected: all Python tests pass.

- [ ] **Step 12.4: Clippy clean**

Run:
```bash
cargo clippy --workspace --exclude ppvm-python-native -- -D warnings
```
Expected: no warnings.

- [ ] **Step 12.5: Confirm commit log**

Run:
```bash
git log --oneline main..HEAD | head -20
```
Expected: see the 9–11 cleanup commits added on top of the existing PR #69 history (depending on whether Step 7 and Step 11 were committed or abandoned).

---

## Summary of changes

| Step | Files touched | Net LOC | Status |
|------|---------------|---------|--------|
| 1 | `stim-parser/src/parser.rs`, `grammar.rs` | ~−30 | Required |
| 2 | `stim-parser/src/extended/interpret.rs`, `parser.rs` | ~−20 | Required |
| 3 | `stim-parser/src/extended/interpret.rs` | ~0 | Required |
| 4 | `stim-parser/src/table.rs`, `parser.rs` | ~−20 | Required |
| 5 | `ppvm-stim/src/executor.rs` | ~−12 | Required |
| 6 | `ppvm-tableau/src/measure.rs`, `ppvm-stim/src/executor.rs` | ~−8 | Required |
| 7 | `ppvm-tableau/src/measure.rs`, `ppvm-stim/src/executor.rs` | ~−5 | Conditional |
| 8 | `ppvm-python-native/src/stim_program.rs`, `interface_tableau.rs` | ~−15 | Required |
| 9 | `ppvm-python-native/src/interface_tableau.rs`, `.pyi`, `generalized_tableau.py` | ~−5 | Required |
| 10 | `ppvm-tableau/tests/gates.rs` | ~+155 (tests) | Required |
| 11 | `stim-parser/src/extended/ast.rs`, `interpret.rs`, `ppvm-stim/src/executor.rs`, `prepare.rs` | ~−30 | Conditional (bail) |
