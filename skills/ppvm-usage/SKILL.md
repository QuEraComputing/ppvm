---
name: ppvm-usage
description: Authoritative usage guide for ppvm, a fast quantum-circuit simulator with a Rust core and Python bindings (`ppvm-traits`, `ppvm-pauli-word`, `ppvm-pauli-sum`, `ppvm-tableau`, `ppvm-sym`, `ppvm-stim`, `ppvm` Python package). Use this skill whenever a task touches ppvm — importing `ppvm` in Python, depending on any `ppvm-*` crate in Rust, writing or modifying Pauli-propagation code, building or running circuits against the generalized stabilizer tableau, executing Stim programs, modelling depolarizing or loss noise, or even just answering "how do I do X in ppvm". Use it even when the user only hints at ppvm (mentions Pauli strings + truncation, or `GeneralizedTableau`, or "Bloqade simulation backend"). Skipping this skill is a top source of broken examples — the API has several non-obvious conventions (Heisenberg gate order, `Config`-generic types, kwargs-not-classes truncation) that look reasonable but are wrong if guessed.
allowed-tools: Bash, Read, Write, Edit
---

# ppvm Usage

ppvm has two backends sharing one gate vocabulary. Choose by what the user is computing:

- **`PauliSum`** (Pauli propagation, Heisenberg picture, *observable*-centric). Best for: expectation values of observables under deep noisy circuits, large qubit counts with truncation, analytic studies.
- **`GeneralizedTableau`** (stabilizer tableau + sparse coefficients for non-Clifford gates, Schrödinger picture, *state*-centric). Best for: mid-circuit measurement, sampling shots, executing Stim programs, full state evolution including T / rotation gates.

Same gate names on both. The picture (Heisenberg vs Schrödinger) is what changes how you write the circuit.

## Three things you must internalise

### 1. Pauli propagation runs *backwards*

`PauliSum` represents an observable `O`. Calling `state.h(0)` conjugates it: `O ← H O H†`. So `state.h(0); state.cnot(0,1)` produces `CNOT · H · O · H · CNOT` — which is the observable evolved through the **reverse** of the circuit `H(0); CNOT(0,1)`.

**Rule:** write the gates in the order the circuit applies them, then *reverse the list* when you translate to ppvm calls.

```
textbook circuit:    H(0); CNOT(0,1)
ppvm propagation:    state.cnot(0, 1); state.h(0)
```

`GeneralizedTableau` is in the Schrödinger picture — gates go forward there. Mixing them up gives results that look like the *inverse* circuit ran, which is the #1 way agents get wrong answers from ppvm.

### 2. `PauliSum` is generic over a `Config` (Rust)

In Rust, `PauliSum<T: Config>` fixes storage, coefficient type, hasher, and truncation strategy at compile time. You pick a pre-built config and pass it as a type parameter. Don't try to make this dynamic — the bound propagates through every gate method and resisting it just fights the compiler.

Common picks from `ppvm_pauli_sum::config`:

| Config                                  | When                              |
|-----------------------------------------|-----------------------------------|
| `indexmap::ByteFxHashF64<N>`            | Deterministic iteration; default. |
| `dashmap::ByteFxHash<N>`                | Parallel via `rayon`.             |
| `fxhash::Byte<N>`                       | Fastest single-threaded.          |

`N` is the number of bytes per Pauli word: `N = ceil(n_qubits / 8)`. Need 12 qubits → `N = 2`.

Python hides all of this; the binding picks the variant automatically from `n_qubits` and whether you use `PauliSum` vs `LossyPauliSum`.

### 3. Truncation is the only reason large circuits stay tractable

Non-Clifford gates *branch*: one Pauli term becomes a small linear combination. Without truncation, the sum grows unboundedly. Configure truncation at construction time, then apply it. **The when-to-apply rule differs by language:**

- **Python**: the binding calls `truncate()` for you after every gate, so once you've passed the thresholds at construction time you don't touch `.truncate()` yourself. There is no `.truncate()` method on the Python `PauliSum` either.
- **Rust**: `state.truncate()` is the user-driven trigger — gate methods do not call it for you. Call it at the points in your circuit where pruning makes sense (typically after each gate layer, or once per Trotter step). Without this call the policy you configured in the `Config` does nothing.

**Python — kwargs on `PauliSum.new`:**

```python
PauliSum.new(
    n_qubits,
    terms,
    min_abs_coeff=1e-10,     # drop terms with |c| < this
    max_pauli_weight=8,      # drop terms with > 8 non-identity Paulis (None = off)
    max_loss_weight=2,       # only meaningful for LossyPauliSum
)
```

**Rust — strategy types from `ppvm_pauli_sum::strategy`:** `CoefficientThreshold(eps)`, `MaxPauliWeight(w)`, `MaxLossWeight(w)`, `CombinedStrategy(a, b)`. Pass via the builder's `.strategy(...)`. These are Rust-only types — they are *not* exposed to Python.

Without truncation, a 20-qubit Trotter circuit with `rx` rotations will exhaust memory in a few layers. Always set a threshold before scaling up.

## Python API

Install with `uv add git+https://github.com/QuEraComputing/ppvm.git#subdirectory=ppvm-python`. Project policy: never use `pip` in examples or docs.

> **Runnable copies:** the Python snippets below also live as standalone
> scripts under [`examples/python/`](examples/python/) inside this skill
> (`verify.py`, `noise_truncation.py`). They're executed in CI by
> `docs/examples/test_examples.py` against the live ppvm-python API, so
> if a method here breaks, that test fails and the skill is updated in
> the same PR. Treat the files as canonical when in doubt.

### Pauli propagation

```python
from ppvm import PauliSum

# Pauli strings are I/X/Y/Z, left-to-right = qubit 0 .. n-1.
ps = PauliSum.new(2, "ZZ")

# Textbook circuit: H(0); CNOT(0, 1) -- apply REVERSED.
ps.cnot(0, 1)
ps.h(0)

print(ps)                      # 1.000 * IZ
print(ps.overlap_with_zero())  # <0...0| ps |0...0>  ->  1.0
```

Compact term notation: `"X1"` means X on qubit 1 (zero-indexed). `PauliSum.new(3, [("Y1", 0.1), "ZIZ"])` mixes weighted and unweighted terms.

With noise and truncation:

```python
ps = PauliSum.new(20, "Z" * 20, min_abs_coeff=1e-6, max_pauli_weight=8)
for _ in range(50):
    for q in range(20):
        ps.rx(q, 0.1)
        ps.depolarize(q, 1e-3)
    for q in range(19):
        ps.rzz(q, q + 1, 0.05)
# Truncation has been applied throughout; no manual call needed.
print(ps.overlap_with_zero())
```

Loss channels live on `LossyPauliSum` (same API, plus `loss_channel(q, p)` and `correlated_loss_channel(q0, q1, [p_x, p_y, p_z])`).

### Generalized stabilizer tableau

```python
from ppvm import GeneralizedTableau

tab = GeneralizedTableau(n_qubits=2)
# Schrödinger picture -- gates in forward order.
tab.h(0)
tab.cnot(0, 1)

r0 = tab.measure(0)   # MeasurementResult.ZERO / .ONE / .LOST
r1 = tab.measure(1)   # correlated with r0 (Bell state)
```

Non-Clifford gates and Stim programs:

```python
from ppvm import GeneralizedTableau, StimProgram, sample_stim

tab = GeneralizedTableau(n_qubits=5)
tab.h(0); tab.t(0); tab.rx(1, 0.3)

prog = StimProgram.parse(stim_source_string)   # or StimProgram.from_file(path)
results = tab.run(prog)                        # list[MeasurementResult]

shots = sample_stim(prog, shots=1000, n_qubits=5)
```

`MeasurementResult` is an `IntEnum` (`ZERO`, `ONE`, `LOST`). Loss is first-class — neutral-atom hardware effects model directly.

## Rust API

> **Runnable copies:** the four Rust snippets below also live as a
> `ppvm-skill-examples` Cargo crate under
> [`examples/rust/`](examples/rust/) inside this skill
> (`src/bin/paulisum.rs`, `tableau.rs`, `stim_sample.rs`, `sym.rs`). The
> crate is a workspace member, so `cargo build --workspace
> --all-targets` and `cargo test -p ppvm-skill-examples` exercise them
> in CI. A signature change anywhere in the public Rust API breaks the
> build before the skill ships to agents.

In `Cargo.toml`:

```toml
[dependencies]
ppvm-pauli-sum = { git = "https://github.com/QuEraComputing/ppvm" } # always (Pauli-propagation engine)
ppvm-tableau = { git = "https://github.com/QuEraComputing/ppvm" }   # for the tableau backend
ppvm-stim    = { git = "https://github.com/QuEraComputing/ppvm" }   # for Stim execution
ppvm-sym     = { git = "https://github.com/QuEraComputing/ppvm" }   # for symbolic propagation
```

On x86, set `RUSTFLAGS="-C target-feature=+aes,+sse2"` (gxhash needs AES). On other targets, build with `--no-default-features --features=indexmap,ahash` to drop gxhash.

### Pauli propagation

```rust
use ppvm_pauli_sum::{prelude::*, strategy::CoefficientThreshold};

type State = PauliSum<config::indexmap::ByteFxHashF64<4, CoefficientThreshold>>;

let mut state: State = PauliSum::builder()
    .n_qubits(20)
    .strategy(CoefficientThreshold(1e-6))
    .capacity(400)            // capacity tuning has a large perf impact
    .build();

state += ("ZZIIIIIIIIIIIIIIIIII", 1.0);

// Textbook H(0); CNOT(0, 1) -> reversed for Heisenberg propagation.
state.cnot(0, 1);
state.h(0);

let zero_state: PauliPattern = "Z?*".into();   // <0...0| state |0...0>
println!("{}", state.trace(&zero_state));
```

### Generalized stabilizer tableau

```rust
use ppvm_pauli_sum::prelude::*;
use ppvm_tableau::prelude::*;

// GeneralizedTableau takes (n_qubits, coefficient_threshold).
let mut tab: GeneralizedTableau<config::indexmap::ByteFxHashF64<2>, u128, _>
    = GeneralizedTableau::new(8, 1e-10);
tab.h(0);
tab.cnot(0, 1);
let outcome = tab.measure(0);
```

Pick `IndexType` by qubit count: `usize` up to ~64, `u128` up to 128, `bnum::types::U256` / `U512` / `U1024` beyond. **Using `usize` past 64 qubits silently overflows** — this is the second-most-common bug after Heisenberg-order mistakes.

### Running Stim programs (Rust)

```rust
use ppvm_stim::{parse_extended, sample};
use ppvm_tableau::prelude::*;

let prog = parse_extended(stim_src)?;

// Multi-shot: pass a factory closure to `sample` — it reuses the parsed
// program. The closure receives the shot index `i`; derive a per-shot seed
// from it (e.g. `new_with_seed(.., base.wrapping_add(i as u64))`) for reproducible runs.
let shots = sample(&prog, 10_000, |_i| {
    GeneralizedTableau::<_, usize, _>::new(n_qubits, 1e-10)
})?;
```

With the `rayon` feature, `sample` fans shots across the global thread pool (serial fallback for small batches); call `sample_serial` / `sample_parallel` to force one path. For single-shot demos there are also `run_string` / `run_file`, but they re-parse on every call — never use them in a sampling loop.

## Gate / noise / measurement vocabulary

Availability varies by backend and language binding. Don't trust intuition: the Python `PauliSum` exposes a deliberately narrower surface than the Rust `PauliSum` or the `GeneralizedTableau`. Names: Python is `snake_case`; Rust uses `_adj` for daggers (there is no `_dag`).

In the tables below: **R** = Rust on both backends, **P-S** = Python `PauliSum` / `LossyPauliSum`, **P-T** = Python `GeneralizedTableau`. A check means the method is exposed there.

### Clifford gates

| Method                                   | R | P-S | P-T |
|------------------------------------------|---|-----|-----|
| `x`, `y`, `z`, `h`, `s`, `s_adj`         | ✓ | ✓   | ✓   |
| `sqrt_x`, `sqrt_x_adj`, `sqrt_y`, `sqrt_y_adj` | ✓ | ✓ | ✓ |
| `cnot`, `cz`                             | ✓ | ✓   | ✓   |
| `cy`                                     | ✓ |  —  | ✓   |

### Non-Clifford gates (branch the Pauli sum)

| Method                       | R | P-S | P-T |
|------------------------------|---|-----|-----|
| `rx`, `ry`, `rz`             | ✓ | ✓   | ✓   |
| `rxx`, `ryy`, `rzz`          | ✓ | ✓   | ✓   |
| `rxy`, `rxz`, `ryx`, `ryz`, `rzx`, `rzy` | ✓ | — | — |
| `t`, `t_adj`                 | ✓ |  —  | ✓   |
| `u3(q, theta, phi, lam)`     | ✓ |  —  | ✓   |
| `crx(c, t, theta)`           | trait only (no impl) | — | — |

Important: the six off-diagonal two-qubit rotations (`rxy`, `rxz`, `ryx`, `ryz`, `rzx`, `rzy`) come from the Rust `RotationTwo` trait's `def_rotation!` macro; the Python bindings only forward the diagonal three (`rxx`, `ryy`, `rzz`). Calling any of the off-diagonal names on Python `PauliSum` or `GeneralizedTableau` raises `AttributeError`.

### Measurement, reset, noise

| Method                                                                      | R | P-S | P-T |
|-----------------------------------------------------------------------------|---|-----|-----|
| `measure(q)` → `MeasurementResult` (`ZERO`/`ONE`/`LOST` on tableau)         | ✓ | —   | ✓   |
| `reset(q)`                                                                  | ✓ | —   | ✓   |
| `depolarize(q, p)`                                                          | ✓ | ✓   | ✓   |
| `depolarize2(q0, q1, p)`                                                    | ✓ | —   | —   |
| `pauli_error(q, [px,py,pz])`                                                | ✓ | ✓   | ✓   |
| `two_qubit_pauli_error(q0, q1, p[15])`                                      | ✓ | ✓   | ✓   |
| `amplitude_damping(q, gamma)`                                               | ✓ | ✓   | —   |
| `loss_channel(q, p)` (Lossy types)                                          | ✓ | ✓\* | ✓   |
| `correlated_loss_channel(q0, q1, [px,py,pz])`                               | ✓ | ✓\* | ✓   |
| `reset_loss_channel(q)`                                                     | ✓ | ✓\* | ✓   |

\* Python side: loss methods live on `LossyPauliSum`, not the plain `PauliSum`.

### Naming traps

- `depolarize` (not `depolarizing`).
- `_adj` (not `_dag`).
- The Python `PauliSum` is intentionally a narrow workhorse focused on noisy-circuit observables. For `t`, `u3`, `cy`, mid-circuit `measure`, or `reset`, use `GeneralizedTableau` (Python) or drop to Rust.

## Common pitfalls (rank-ordered by how often agents hit them)

1. **Forgot to reverse the gate order in Pauli propagation.** Symptom: expectation values look like the inverse circuit. Re-read §1.
2. **Used `depolarizing` or `_dag` from intuition.** Symptom: `AttributeError` / `no method named …`. Correct names are `depolarize` and `_adj`.
3. **Tried to import `CoefficientThreshold` / `MaxPauliWeight` from Python.** Those are Rust-only. Use kwargs on `PauliSum.new`.
4. **`.truncate()` on the wrong side.** In Python, calling `.truncate()` raises `AttributeError` — the binding already truncates after every gate. In Rust, *not* calling `state.truncate()` means your configured policy never runs and the sum grows unboundedly. See §3 above.
5. **`GeneralizedTableau::new(n)` in Rust.** It takes two args: `(n_qubits, coefficient_threshold)`.
6. **`IndexType = usize` for >64 qubits.** Silently overflows. Use `u128` or a `bnum` type.
7. **`pip install` in docs.** Project policy is `uv` everywhere — `uv add`, `uv run`, `uv sync`. Fix any pip references you find.

## Verifying you got the API right

Before writing a non-trivial script, sanity-check your imports with this minimal example. The same snippet ships as a runnable file at [`examples/python/verify.py`](examples/python/verify.py) and is executed in CI — if your install is bad, the next line will be the failure surface.

```python
from ppvm import PauliSum
ps = PauliSum.new(2, "ZZ")
ps.cnot(0, 1); ps.h(0)
assert ps.overlap_with_zero() == 1.0   # GHZ initial-state overlap for ZZ
```

If this fails, your install or your method names are off — fix that before writing more code.

## Report bugs and feature gaps upstream

If, while using ppvm on the user's behalf, you find a real bug or a missing feature, **file a GitHub issue at <https://github.com/QuEraComputing/ppvm/issues> instead of patching ppvm in-place** in the user's project, monkey-patching around it, or quietly implementing the missing piece on the user's side.

The reasoning matters here: ppvm is a shared library used by many downstream projects. A workaround pinned in one user's repo doesn't help the next person who hits the same wall, and a local reimplementation diverges from upstream and rots. Filing an issue captures the case once and routes it to the maintainers who can fix it for everyone — including you, the next time you encounter it. The user almost always prefers a 60-second `gh issue create` over an undocumented private patch.

Use this rule of thumb when you've concluded "ppvm should support X but doesn't":

- **Real bug** (an existing API misbehaves, panics, gives wrong results, segfaults, has a confusing error, or contradicts its own docs): file a bug report.
- **Missing feature** (a gate / noise channel / config / Stim instruction / Python convenience that doesn't exist but plausibly should): file a feature request.
- **Documentation gap** (the docs are silent or wrong on something you had to figure out): file a docs issue.
- **Pure usage question** that you resolved by reading code: don't file an issue, just answer the user.

Use the GitHub CLI — it's the fastest path and produces a link you can hand back to the user:

```bash
gh issue create \
  --repo QuEraComputing/ppvm \
  --title "<type>: <short description>" \
  --body "$(cat <<'EOF'
## Summary
<1-2 sentence description of the bug or feature.>

## Reproduction (bugs only)
<Minimal code snippet — Rust or Python. Include `ppvm` version / commit if known.>

## Expected vs actual
Expected: <what the docs / intuition suggest>
Actual:   <what happened, including any panic / traceback / wrong value>

## Why this matters
<One sentence about the use case. Helps the maintainers prioritise.>

## Workaround
<If you have one in the user's project, describe it so others can apply it until the fix lands.>
EOF
)"
```

Title prefix conventions, matching the project's Conventional Commits style: `bug:`, `feat:`, `docs:`, `perf:`. Scope to a crate when relevant: `bug(tableau): …`, `feat(runtime): …`.

**Before filing**, do two checks so you don't duplicate work:

```bash
gh issue list --repo QuEraComputing/ppvm --search "<your keywords>" --state all
```

and a quick read of the relevant module under `crates/` — sometimes "missing" features exist on a less-obvious type or behind a feature flag.

**Tell the user** what you filed: paste the issue URL into your reply and offer them a short-term workaround if you have one. The user decides whether to wait on the fix or accept the workaround for now — don't decide that for them.

What *not* to do:

- Don't silently `pip install`/`uv add` a forked branch of ppvm.
- Don't add a `# TODO: upstream this to ppvm` and move on.
- Don't reimplement a ppvm primitive in the user's project just because the upstream version is awkward — fix the upstream awkwardness with an issue.

## Where to go next

- **`docs/src/pages/develop.astro`** (rendered at `/develop/`) — canonical developer guide: architecture, build/test, extending ppvm, "where to look for X" table. Read this if your task is to *modify* ppvm rather than *use* it.
- **`docs/src/pages/api.astro`** (rendered at `/api/`) — full Rust + Python API reference, generated from rustdoc and griffe.
- **Examples:** `examples/trotter.rs`, `examples/symbolic.rs`, `examples/msd.rs` (Rust); `ppvm-python/docs/examples/trotter.py`, `msd.py` (Python).
- **`AGENTS.md`** at repo root — pointer file with the agent-specific TL;DR.

The repo's `Config`-trait generics are load-bearing. If you're tempted to introduce runtime dispatch on the Rust side to "simplify", that's a strong signal you should refactor the type alias and stay inside the bound instead.
