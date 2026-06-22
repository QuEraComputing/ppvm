# Mirroring the `stim.TableauSimulator` API in ppvm

**Date:** 2026-06-17
**Status:** Design — pending review
**Topic:** Reshape ppvm's gate/measurement/noise API surface to mirror `stim.TableauSimulator`

---

## 1. Goal & motivation

Make ppvm's simulator API **muscle-memory compatible** with `stim.TableauSimulator`, so
that anyone familiar with stim can use ppvm without relearning method names or calling
conventions.

The goal is **drop-in familiarity**, not a literal runnable-stim-code clone. Python adopts
stim's *naming* and *calling conventions* (broadcasting over targets, keyword-only noise
probabilities, stim method names), while Rust adopts the names but keeps scalar method
calls scalar and exposes broadcasting through explicit batch helpers. We deliberately keep
ppvm's semantics where they genuinely differ (loss-aware measurement, fixed qubit count,
no pure-stabilizer state inspection).

## 2. Decisions (locked during brainstorming)

| Decision | Choice |
|---|---|
| Motivation | Drop-in familiarity |
| Placement | **Replace** the existing API (migrate methods + all callers), not an additive facade |
| Rust vs Python | **Python mirrors stim calling conventions**; Rust uses stim names plus explicit `*_batch` methods |
| Trait scope | **Change the shared traits in `ppvm-traits`** — both the tableau backend and the Pauli-propagation (`PauliSum`) backend move together |
| Method families in scope | Gates + broadcasting; Measurement + reset; Noise |
| State inspection | **Out of scope** (`peek_*`, `current_inverse_tableau`, `canonical_stabilizers`, `state_vector`, `measure_kickback`, `set_state_from_*`, `postselect_*`) |
| Rust target shape | **Scalar methods stay scalar** — batches use explicit `*_batch` helpers |
| `current_measurement_record` | **Include it** — add measurement-record state |

## 3. Architecture context (why the blast radius is what it is)

The gate traits are defined once in `crates/ppvm-traits/src/traits/` and implemented by
**both** simulator backends:

- **Pauli propagation:** `PauliSum<T>` — `crates/ppvm-pauli-sum/src/sum/` (Heisenberg picture).
- **Tableau:** `Tableau<T>` and `GeneralizedTableau<T, ...>` — `crates/ppvm-tableau/`
  (Schrödinger picture). `ppvm-tableau` depends on `ppvm-traits`, so the traits sit below
  both.

Because we are changing the **shared traits**, the changes ripple into both backends. The
families split by where they live:

| Trait surface | Implementers | In this refactor |
|---|---|---|
| `Clifford`, `CliffordExtensions` | `PauliSum` + `Tableau` + `GeneralizedTableau` | renamed + explicit batches |
| `RotationOne` / `RotationTwo` | `PauliSum` + `GeneralizedTableau` | scalar names + explicit batches |
| Noise (`PauliError`, `Depolarizing`, `Depolarizing2`, `TwoQubitPauliError`, loss channels) | `PauliSum` + tableaux | renamed + explicit batches |
| `TGate` (`t`, `t_dag`) | `GeneralizedTableau` only | renamed (`t_adj`→`t_dag`) |
| `Measure` / `LossyMeasure` / `Reset` | tableaux only | reshaped + record |

## 4. Rust: scalar names plus explicit batches

stim broadcasts over `*targets`, and the Python API mirrors that convention. Rust keeps
the scalar method names scalar because `gate(q)` and `gate_many(qs)` are semantically
different operations in Rust API design. Batched application uses explicit `*_batch`
helpers:

```rust
tab.h(0);                         // single
tab.h_batch(&[0, 1, 2]);          // explicit batch
tab.cnot(0, 1);                   // single pair
tab.cnot_batch(&[(0, 1), (2, 3)]); // explicit pair batch
```

### Broadcasting semantics

- **Python:** follows stim target-list broadcasting.
- **Rust scalar methods:** accept one target or one pair.
- **Rust `*_batch` methods:** apply targets or pairs in order and are the fast path for
  Python's broadcast wrappers.

## 5. Naming map

| Current ppvm | New (stim) name | Aliases added | Backends touched |
|---|---|---|---|
| `s_adj` | `s_dag` | — | all three |
| `sqrt_x_adj` | `sqrt_x_dag` | — | all three |
| `sqrt_y_adj` | `sqrt_y_dag` | — | all three |
| `t_adj` | `t_dag` | — | GeneralizedTableau |
| `cnot` | `cnot` | `cx`, `zcx` | all three |
| `cz` | `cz` | `zcz` | all three |
| `cy` | `cy` | `zcy` | all three |
| `depolarize` | `depolarize1` | — | all three |
| `depolarize2` | `depolarize2` | — | all three |
| `pauli_error([px,py,pz])` | retained as ppvm extension | adds `x_error`/`y_error`/`z_error` | all three |
| `two_qubit_pauli_error` | retained as ppvm extension | — | all three |
| `reset` | `reset` | `reset_z` (alias); add `reset_x`, `reset_y` | tableaux |

Aliases are thin forwarders (Rust: `#[inline]` methods or default-trait methods; Python:
assignment / thin wrappers).

`x_error`/`y_error`/`z_error` are convenience wrappers over `pauli_error` with the relevant
single-axis probability. `reset_x`/`reset_y` are basis-change + `reset_z` (`reset`).

## 6. Measurement & reset

- `measure(target) -> MeasurementResult` — **single target, no broadcast** (mirrors stim's
  deliberate footgun-avoidance). Return type stays `MeasurementResult` (`ZERO`/`ONE`/`LOST`)
  in Python and `Option<bool>` in Rust — `LOST` is a real outcome ppvm must express and
  cannot collapse to `bool`.
- `measure_many(targets) -> Vec<MeasurementResult>` / `list[MeasurementResult]` — broadcast.
- Python `reset` / `reset_x` / `reset_y` / `reset_z` broadcast over targets; Rust keeps
  scalar reset methods plus explicit `*_batch` helpers.

### Measurement record (new state)

Add a measurement log to the tableau state:

- Rust: a `Vec<Option<bool>>` (or `Vec<MeasurementResult>`) field on the tableau, appended on
  every `measure` / `measure_many`, and on measurements performed by `run`/`do`.
- Exposed as `current_measurement_record() -> Vec<...>` / `list[MeasurementResult]`.
- **Copied** by `fork`, `__copy__`, `__deepcopy__` (record travels with the state).
- Cleared by an explicit reset of the simulator only (not by qubit `reset`, matching stim).

## 7. Noise

stim-shaped, **keyword-only `p`** in Python:

- `x_error(*targets, p=...)`, `y_error`, `z_error`
- `depolarize1(*targets, p=...)`, `depolarize2(*targets, p=...)` (pairs)

Rust scalar equivalents take one target (or one target pair) plus `p`; broadcast behavior
uses explicit `*_batch` methods. These wrap the existing `pauli_error`/`depolarize`
machinery (which is implemented via the `impl_tableau_noise!` macro delegating to
`TableauLike` methods — so the per-type impl blocks need few edits; the trait definition
and `TableauLike` impl carry the changes).

ppvm-specific channels with no stim equivalent are **retained** under their current names:
`pauli_error`, `two_qubit_pauli_error`, `loss_channel`, `correlated_loss_channel`,
`reset_loss_channel`, and the `two_qubit_pauli_error_probabilities` helper.

## 8. Python layer

- `*targets` varargs collected into a sequence and passed to native wrappers that dispatch
  to Rust `*_batch` methods. Single-qubit, two-qubit-pair, and odd-count-error semantics
  are enforced consistently.
- Method renames mirror the Rust renames; aliases (`cx`, `zcx`, `zcz`, `zcy`, `reset_z`)
  added.
- Noise methods expose keyword-only `p`.
- `MeasurementResult` enum unchanged.
- Existing `run(StimProgram)` / `sample(...)` retained; `run` appends to the measurement
  record. Optionally add `do` as an alias of `run` for stim familiarity (low cost — included).
- Update `.pyi` stubs and mixins (`CliffordMixin`, `CliffordExtensionMixin`, `RotationsMixin`,
  `NoiseMixin`) to the new surface.

## 9. Explicitly NOT mirrored (and why)

| stim API | Reason omitted |
|---|---|
| `peek_bloch`, `peek_x/y/z`, `peek_observable_expectation` | Assume a pure stabilizer state; ambiguous over a coefficient-weighted superposition |
| `current_inverse_tableau`, `canonical_stabilizers`, `set_state_from_*` | Pure-stabilizer concepts; not well-defined for `GeneralizedTableau` |
| `state_vector` | Well-defined but expensive; out of selected scope |
| `measure_kickback`, `postselect_*` | Out of scope; stabilizer-state assumptions |
| auto-growing `num_qubits` | ppvm fixes `n_qubits` at construction; tableau + sparse coefficient vector are not cheap to resize |
| `measure` returning plain `bool` | ppvm has a third outcome (`LOST`) |

## 10. Migration plan & blast radius

Order of operations:

1. **`ppvm-traits/src/traits/`** — rename trait methods; keep scalar signatures; add
   explicit batch defaults + aliases; reshape `Measure`/`Reset`/noise traits.
2. **`ppvm-pauli-sum/src/sum/`, `ppvm-pauli-word/src/{phase,word,loss}/`** — update `PauliSum` and blanket
   `PauliWordTrait` impls to the new trait surface.
3. **`ppvm-tableau/`** — update `Tableau` and `GeneralizedTableau` impls
   (`gates/clifford.rs`, `tgate.rs`, `rot1.rs`, `rot2.rs`, `measure.rs`, `gates/reset.rs`,
   `noise.rs`); add the measurement record to the tableau state; wire `*_batch` as the
   broadcast fast path.
4. **`ppvm-stim/src/executor.rs`** — update call sites; append measurements to the record.
5. **`ppvm-python-native/`** — update PyO3 bindings (`interface.rs`, `interface_tableau.rs`)
   for the new names, sequence-accepting targets, keyword `p`, and `current_measurement_record`.
6. **`ppvm-python/src/ppvm/`** — update mixins, `generalized_tableau.py`, `.pyi` stubs;
   add aliases; add `do`.
7. **Tests / benches / examples** — ~40 Rust files reference the renamed methods
   (`ppvm-traits`, `ppvm-pauli-sum`, `ppvm-tableau`, `ppvm-stim` tests/benches/examples)
   plus Python tests.
   Update all call sites.

**No backward-compatibility shims** — this is a clean replace per the placement decision. A
single rename pass + compiler-driven fixup is the mechanism; the type system surfaces every
stale call site.

## 11. Testing

- Existing gate/measure/noise/loss tests updated to new names — they remain the correctness
  oracle (behavior must be unchanged; only the surface moves).
- New tests for **batching / Python broadcasting**: Rust `h_batch(&[0,1,2])` ≡ three
  `h` calls, pair-batch for two-qubit gates, Python odd-count error.
- New tests for the **measurement record**: ordering, `measure_many`, propagation through
  `fork`/copy, population by `run`.
- Python parity tests asserting the wrapper broadcasts identically to the native layer.
- Both backends (`PauliSum` and tableaux) exercised for the shared renamed Clifford/rotation/
  noise surface.

## 12. Open considerations (non-blocking)

- Whether the measurement record stores `MeasurementResult` (keeps `LOST`) or `bool` (stim
  parity, but lossy). **Proposed:** keep `Option<bool>`/`MeasurementResult` to preserve
  loss information.
- Whether Rust should grow higher-level builders for circuit-style bulk application in
  addition to `*_batch`. **Proposed:** keep the PR scoped to explicit batch helpers.
