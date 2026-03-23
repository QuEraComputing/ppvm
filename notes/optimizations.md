# Optimization Ideas

> **Profiling basis:** All findings and estimates below are derived from flamegraph profiling of the **85-qubit MSD circuit** (`target/release/examples/flamegraph.svg`, 10k iterations). Speedup projections are specific to that workload and may not generalize to other circuit shapes.

## 1. Masked single-qubit Clifford gates

### Background

The MSD flamegraph shows:

| Hotspot | % of runtime |
|---|---|
| H gate (`PhasedPauliWord::h` + `PauliWord::h`) | ~37% |
| Measurement (`LossyMeasure::measure`) | ~27% |
| CZ gate | ~14% |
| S gate | ~13% |
| `MulAssign` | ~7% |
| `compute_decomposition` | ~6% |

H and S together dominate at ~50%. The root cause: `sqrt_y`, `sqrt_x`, and their adjoints decompose into sequences of 5 single-qubit H/S calls. Each call iterates over all `2n` tableau rows. Applying `sqrt_y` to 17 qubits in a block means 85 separate row-scan passes.

### The optimization

Add `h_mask` and `s_mask` methods that apply the gate to a set of qubits (given as a `BitArray<A>` mask) in a single pass over tableau rows, using word-level operations.

**Within each row**, the word-level math is:
- H: `phase += count_ones(xw & zw & mask) * 2`, then swap `xw ↔ zw` under the mask
- S: `phase += count_ones(xw & zw & mask) * 2` (Tableau convention), then `zw ^= xw & mask`

Phase accumulation uses POPCNT — the same trick as the `MulAssign` speedup.

**Across rows**, gate fusion collapses k separate row-scan loops into one.

For a layer of k simultaneous gates these two effects multiply. For k=17 (one code block):

| k | Loop reduction | Word-level inner | Combined (rough) |
|---|---|---|---|
| 4 | ~4× | ~3× | ~10× |
| 8 | ~8× | ~4× | ~25× |
| 17 | ~17× | ~5× | ~50× (theoretical) |

Realistically **10–20× on the H/S portion** → **~5–10× end-to-end** on the MSD workload.

### Implementation scope

~100 lines across 3 files, ~1 day of work:

- `word/data.rs`: `PauliWord::h_mask`, `PauliWord::s_mask` (~16 lines). Same raw-slice pattern as `MulAssign`, needs `<A as BitView>::Store: PrimInt`.
- `phase/clifford.rs`: `PhasedPauliWord::h_mask`, `PhasedPauliWord::s_mask` (~16 lines). POPCNT phase + delegate to word.
- `tableau/clifford.rs`: `Tableau::h_mask`, `Tableau::s_mask`, `GeneralizedTableau` wrappers (~30 lines).

**Key gotcha**: `Tableau::s_mask` must **not** delegate to `PhasedPauliWord::s_mask`. `Tableau::s` uses the forward-propagation phase convention (`x & z`), while `PhasedPauliWord::s` uses backward-propagation (`x & !z`). Inline the phase as `count_ones(xw & zw & mask) * 2`, then call `word.s_mask(mask)` for the bit op only.

No trait changes required to start — inherent methods on `Tableau`/`GeneralizedTableau` are sufficient since call sites use concrete types.

### CZ / CNOT

Not easily vectorizable word-level (cross-position bit interactions between control and target). Gate fusion (batching k pairs in one row pass) still applies but requires a different API (`cz_pairs(&[(usize, usize)])`). Lower priority.
