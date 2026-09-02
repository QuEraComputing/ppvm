# PPVM vihaco 0.4.0 changes

This changelog describes the completed migration of `ppvm-vihaco` and its
supporting ISA crate from vihaco 0.1.1 to vihaco 0.4.0. It is an informative
reference for users and contributors who need to understand what changed,
especially in `.sst` instruction syntax.

PPVM now uses the vihaco 0.4 model:

```text
source instruction -> surface instruction -> runtime instruction
```

The lowering step resolves names, types, literals, and addresses while
preserving one source instruction per runtime instruction. Syntax sugar that
expands one source instruction into several runtime instructions is not part
of the language.

## Migration summary

The migration is implemented across the PPVM runtime, circuit ISA, parser,
bytecode codec, CLI/TUI call sites, and `.sst` fixtures. The implementation
uses v0.4 surface and runtime instruction types, generated component
instruction sets, explicit one-to-one lowering, and bytecode v2 metadata for
functions, labels, signatures, and the entry point.

The full workspace test suite is the semantic verification gate. It covers
source execution, bytecode round trips, function calls, branching,
measurements, all PPVM backends, traces, loss, truncation, and source/bytecode
parity.

The implementation was guided by Stellarscope PR #75, commit `b9e24258`.

## Semantic preservation

This was an API and syntax migration, not a simulator-behavior migration.
The implementation preserves PPVM semantics. In particular, it preserves:

- Heisenberg/Pauli propagation direction or gate ordering;
- operand stack ordering, including the order in which qubit addresses and
  floating-point parameters are popped;
- measurement and reset behavior;
- tableau, PauliSum, and LossyPauliSum backend results;
- truncation thresholds, maximum Pauli weights, and loss-channel behavior;
- branch and conditional-branch behavior after label resolution;
- function-call behavior, return-value counts, and function boundaries;
- string-table values used by `trace`, `print`, and other string-consuming
  instructions;
- seeded randomness and shot-level reproducibility;
- bytecode round-trip behavior and equivalence between source and bytecode
  execution.

Source/runtime instruction conversion was verified in two ways:

1. **Implementation review:** the surface-to-runtime lowering preserves
   operand values, operand order, instruction order, control-flow targets,
   and backend dispatch. Any source-only construct, such as a label, is
   explicitly identified as metadata rather than silently dropped.
2. **Existing behavioral tests:** the tests already in this repository remain
   the semantic oracle. Their source spellings and expected enum names were
   migrated without weakening or replacing their behavioral assertions.

The completed implementation passes the existing tests for
GHZ/Bell circuits, measurements, branching, function calls, traces,
PauliSum/LossyPauliSum behavior, truncation, reset, bytecode serialization,
and source/bytecode parity. Parser tests additionally assert the one-to-one
source-to-runtime mapping, while the existing execution tests assert that the
mapping preserves behavior.

## Dependency changes

### Workspace and PPVM dependencies

PPVM now uses the published 0.4.0 crates consistently:

- `vihaco = "0.4.0"`
- `vihaco-cpu = "0.4.0"`
- `vihaco-parser = "0.4.0"`
- `vihaco-parser-derive = "0.4.0"`

The migration removed direct dependencies on:

- `vihaco-parser-core`
- `vihaco-derive`
- the old 0.1.1 vihaco crates

`vihaco-circuit-isa` was migrated as well. The workspace now has one vihaco
version, and its instruction type is usable as a v0.4 component instruction
set.

### Parser imports and derives

Old parser-core imports were replaced:

```rust
use vihaco_parser_core::Parse;
```

with:

```rust
use vihaco_parser::Parse;
```

`vihaco_parser_derive::Parse` is used for derives required by the 0.4 API.
The 0.4 derive uses `#[pattern = ...]` and syntax-class metadata; the old
`#[head]`, `#[token]`, and `#[delimiters]` form is no longer the general
instruction-definition mechanism.

## Component and instruction-set changes

### Circuit component

The circuit instruction set is defined with `vihaco::component!`. Its generated
module exposes separate surface and runtime types:

```rust
vihaco::component! {
    pub component Circuit {
        // state fields
    }

    instruction {
        X,
        H,
        // ...
    }
}

pub use circuit::{runtime, syntax};
```

The generated runtime instruction type is used by `#[dispatch]` and by the
composite. The generated syntax instruction type is used by the parser. This
implements `HasInstructionSet`, which the v0.4 `#[composite]` macro requires.

The hand-written `vihaco_circuit_isa::CircuitInstruction` was replaced by the
generated v0.4 instruction type. `CircuitMessage` and `CircuitEffect` remain
application-level payload types.

### CPU instruction types

vihaco 0.4 distinguishes:

- `vihaco_cpu::SurfaceInstruction`: parsed source instructions;
- `vihaco_cpu::RuntimeInstruction`: executable instructions;
- `vihaco_cpu::SurfaceType` and `SurfaceValue`: source-level typed operands.

PPVM parses into surface instructions and lowers them to
`RuntimeInstruction`, resolving labels, function names, types, and
string-table indices along the way.

### Composite instruction type

The composite-generated instruction enum is the runtime instruction type
generated by the composite, for example:

```rust
pub type PPVMInstruction = ppvm::runtime::Instruction;
```

Conversions now exist from the runtime CPU and runtime circuit instruction
types into that enum. The loader uses the v0.4 `ProgramImage` and
`LocalModule`, replacing the removed `Module`/old loader API combination.

## Parser and resolver changes

The old `ParsedModule`/`BodyItem`/`RawForm` pipeline was replaced by the v0.4
pipeline:

```rust
ParsedModule<PPVMSurfaceInstruction, SurfaceType, Vec<PPVMHeader>>
```

and:

```rust
impl Resolve<PPVMSurfaceInstruction, SurfaceType, PPVMHeader> for PPVMResolver {
    type Module = LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo>;
}
```

The resolver now performs these operations:

1. Applies device headers to `PPVMDeviceInfo`.
2. Assigns function IDs and addresses.
3. Records source labels as metadata and resolves branch/call targets.
4. Converts CPU surface instructions to runtime instructions.
5. Converts circuit surface instructions to circuit runtime instructions.
6. Interns string constants into the module string table.
7. Populates functions, labels, `main_function`, code, strings, and device info.

There is no generic fallback for arbitrary `RawForm` values. Unknown
instructions and malformed operands fail during surface parsing.

## Instruction syntax changes

The v0.4 composite syntax qualifies an instruction with both the composite
field and the component dialect:

```text
<composite field>::<component dialect>.<mnemonic>
```

For PPVM this means CPU instructions use `cpu::cpu.*`, and circuit
instructions use `circuit::circuit.*`.

### CPU instructions

| Old PPVM spelling | v0.4 spelling | Runtime result |
| --- | --- | --- |
| `const.u64 0` | `cpu::cpu.const u64, 0` | `RuntimeInstruction::Const(Type::U64, Value::U64(0))` |
| `const.u32 1` | `cpu::cpu.const u32, 1` | `RuntimeInstruction::Const(Type::U32, Value::U32(1))` |
| `const.i64 -1` | `cpu::cpu.const i64, -1` | typed `Const` |
| `const.f64 0.5` | `cpu::cpu.const f64, 0.5` | typed `Const` |
| `const.bool true` | `cpu::cpu.const bool, true` | typed `Const` |
| `const.str "Z?*"` | `cpu::cpu.const str, "Z?*"` | one typed `Const`, with string interning |
| `ret` | `cpu::cpu.ret 0` | `RuntimeInstruction::Return(0)` |
| `ret 1` | `cpu::cpu.ret 1` | `RuntimeInstruction::Return(1)` |
| `br @done` | `cpu::cpu.br @done` | one patched `Branch` |
| `cond_br @yes, @no` | `cpu::cpu.cond_br @yes, @no` | one patched `ConditionalBranch` |
| `call 0, @measure` | `cpu::cpu.call 0, measure` | one patched `Call` |
| `@label:` | `cpu::cpu.label @label` | source label metadata; no runtime instruction |
| `breakpoint` | `cpu::cpu.breakpoint` | one `Breakpoint` |
| `halt` | `cpu::cpu.halt` | one `Halt` |
| `print` | `cpu::cpu.print` | one `Print` |
| `heap_alloc 4` | `cpu::cpu.heap_alloc 4` | one `HeapAlloc` |
| `load ...` | `cpu::cpu.load ...` | one typed `Load` |
| `store ...` | `cpu::cpu.store ...` | one typed `Store` |
| arithmetic/logical CPU op | `cpu::cpu.<mnemonic> ...` | one runtime CPU instruction |

The exact operand spelling for typed CPU instructions follows the v0.4 CPU
surface grammar. In particular, types and values are separate source
operands: `const f64, 0.5`, not `const.f64 0.5`.

`const.str` is not a multi-instruction expansion. It is a source-level typed
constant whose string value is interned while producing one runtime `Const`.
The old spelling was removed because it encoded the type in the mnemonic.

### Circuit instructions

Circuit instructions do not take qubit or numeric operands in their mnemonic.
PPVM pushes operands with CPU constants, then executes one circuit instruction
that consumes the required values from the CPU stack.

| Old PPVM spelling | v0.4 spelling | Values consumed from stack |
| --- | --- | --- |
| `circuit.x` | `circuit::circuit.x` | qubit |
| `circuit.y` | `circuit::circuit.y` | qubit |
| `circuit.z` | `circuit::circuit.z` | qubit |
| `circuit.h` | `circuit::circuit.h` | qubit |
| `circuit.s` / `circuit.s_adj` | `circuit::circuit.s` / `circuit::circuit.s_adj` | qubit |
| `circuit.sqrt_x` / `circuit.sqrt_x_adj` | `circuit::circuit.sqrt_x` / `circuit::circuit.sqrt_x_adj` | qubit |
| `circuit.sqrt_y` / `circuit.sqrt_y_adj` | `circuit::circuit.sqrt_y` / `circuit::circuit.sqrt_y_adj` | qubit |
| `circuit.t` / `circuit.t_adj` | `circuit::circuit.t` / `circuit::circuit.t_adj` | qubit |
| `circuit.cnot` | `circuit::circuit.cnot` | two qubits |
| `circuit.cz` | `circuit::circuit.cz` | two qubits |
| `circuit.rx`, `ry`, `rz` | `circuit::circuit.rx`, `ry`, `rz` | qubit, float |
| `circuit.rxx`, `ryy`, `rzz` | `circuit::circuit.rxx`, `ryy`, `rzz` | two qubits, float |
| `circuit.r` | `circuit::circuit.r` | qubit, two floats |
| `circuit.u3` | `circuit::circuit.u3` | qubit, three floats |
| `circuit.measure` | `circuit::circuit.measure` | qubit |
| `circuit.reset` | `circuit::circuit.reset` | qubit |
| `circuit.depolarize` | `circuit::circuit.depolarize` | qubit, float |
| `circuit.depolarize2` | `circuit::circuit.depolarize2` | two qubits, float |
| `circuit.paulierror` | `circuit::circuit.pauli_error` | qubit, three probabilities |
| `circuit.two_qubit_pauli_error` | `circuit::circuit.two_qubit_pauli_error` | two qubits, fifteen probabilities |
| `circuit.loss` | `circuit::circuit.loss` | qubit, float |
| `circuit.correlated_loss` | `circuit::circuit.correlated_loss` | two qubits, three probabilities |
| `circuit.trace` | `circuit::circuit.trace` | pattern string |
| `circuit.truncate` | `circuit::circuit.truncate` | none |

The generated component defines the names explicitly where needed. Each
listed source instruction maps to exactly one generated circuit runtime
instruction.

### Device headers

The PPVM device headers remain conceptually the same and are now parsed with
the v0.4 header syntax:

```text
device circuit.n_qubits 2;
device circuit.backend paulisum;
device circuit.observable ZZ;
device circuit.coefficient_threshold 1e-10;
device circuit.max_pauli_weight 8;
```

They are metadata, not runtime instructions, and their meaning is unchanged.

## One-to-one lowering behavior

The implementation follows these rules:

- A typed `const` produces one runtime `Const`.
- `ret`, branch, conditional branch, and call each produce one runtime CPU
  instruction.
- A circuit mnemonic produces one runtime circuit instruction.
- A source label records a label address but produces no executable runtime
  instruction. Labels are metadata, not instructions.
- Symbolic branch and call targets may be resolved in a separate pass, but
  resolution replaces operands on the same instruction rather than emitting
  additional instructions.
- No `play`, `poly`, `digi`, combined gate, or measurement shorthand expands
  into a constant plus another instruction.
- Stack setup is written explicitly in the source. For example, a
  measurement requires an explicit `cpu::cpu.const u64, 0` followed by
  `circuit::circuit.measure`.

## Bytecode and runtime changes

Bytecode serialization and tests were updated for v0.4 runtime instructions:

- the codec uses `RuntimeInstruction`, not the removed
  `vihaco_cpu::Instruction`;
- typed `Const` instructions include their `Type` argument;
- modules use `LocalModule` and the v0.4 module fields;
- the existing PPVM string table and device-info serialization is preserved;
- labels and function metadata are serialized so source and bytecode execution have
  the same resolved control-flow and call behavior. PPVM bytecode v2 stores
  function signatures, function ranges, labels, and `main_function`;
- round-trip tests cover every runtime instruction family used by PPVM.

## Call-site and fixture changes

All embedded programs and `.sst` fixtures were updated, including:

- `crates/ppvm-vihaco/tests/*.sst`;
- `crates/ppvm-vihaco/src` test programs;
- `crates/ppvm-cli/examples/*.sst`;
- CLI and TUI breakpoint/trace programs;
- README and documentation examples.

For every fixture, verification covers both:

1. the new source parses into the intended surface instruction; and
2. resolution emits exactly one corresponding runtime instruction per source
   instruction, excluding source-only labels.

The fixtures run through the existing execution assertions, which compare
measurements, traces, states, branch outcomes, and bytecode results with the
pre-migration behavior. A fixture that parses but changes a result fails the
migration's semantic-preservation requirement.

## Verification

The completed migration was verified with:

```bash
cargo fmt --all -- --check
cargo check -p ppvm-vihaco
cargo test -p vihaco-circuit-isa
cargo test -p ppvm-vihaco
cargo test -p ppvm-cli
cargo test --workspace
```

The workspace test suite passes, including PPVM source/bytecode fixtures,
CLI/TUI tests, circuit ISA tests, and the existing simulator behavior tests.

The dependency tree was also checked to contain only the v0.4 vihaco family:

```bash
cargo tree -i vihaco@0.1.1
```

That command reports no PPVM-related dependency path for the migrated
`vihaco-circuit-isa`.
