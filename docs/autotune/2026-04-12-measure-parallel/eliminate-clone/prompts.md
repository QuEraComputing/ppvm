# Eliminate coefficient clone in measurement

## Hypothesis
The `measure()` method in `crates/ppvm-tableau/src/measure.rs` (the `LossyMeasure` impl for `GeneralizedTableau`) clones `self.coefficients` into a `FxHashMap` at line 92-97. This clone accounts for ~25-36% of measurement time. However, `self.coefficients` is never read again after the clone — it's overwritten at line 192 (case a) or via `trim_coefficients_for_measurement` (case b). We can use `std::mem::replace` to take ownership instead of cloning.

## What to change
File: `crates/ppvm-tableau/src/measure.rs`

1. Replace the clone:
```rust
// BEFORE:
let mut coeff_map: HashMap<I, Complex<T::Coeff>> = self
    .coefficients
    .clone()
    .into_iter()
    .map(|(v, i)| (i, v))
    .collect();

// AFTER:
let mut coeff_map: HashMap<I, Complex<T::Coeff>> = std::mem::replace(&mut self.coefficients, C::new())
    .into_iter()
    .map(|(v, i)| (i, v))
    .collect();
```

2. For case_b (line 205+), the current code calls `self.trim_coefficients_for_measurement(destab_anticomm_bits, outcome, z_sign)` which operates on `self.coefficients`. Since we took ownership of the coefficients, we need to implement the trim logic directly on `coeff_map`. Look at the `trim_coefficients_for_measurement` method in `data.rs` to understand the logic:
- It iterates over coefficients, keeps only entries where `(symplectic_inner(alpha, destab_anticomm_bits) % 2 != 0) ^ z_sign) == outcome`
- Then normalizes if anything was removed

For case_b, convert coeff_map back to self.coefficients with filtering:
```rust
// In case_b:
let z_sign = phase_decomp == 2;
let old_len = coeff_map.len();
self.coefficients = C::new();
for (idx, coeff) in coeff_map {
    let parity = symplectic_inner(idx, destab_anticomm_bits) % 2 != 0;
    if (parity ^ z_sign) == outcome {
        self.coefficients.unsafe_insert(idx, coeff);
    }
}
if self.coefficients.len() < old_len {
    self.coefficients.normalize();
}
```

## Files in scope
- `crates/ppvm-tableau/src/measure.rs` — the only file to modify

## Files NOT to touch
- All benchmarks, tests, examples, docs

## Expected impact
Eliminate ~25-36% of first-measurement time by avoiding a full coefficient vector clone+allocation.

## Commit instruction
Commit all changes with message: "perf(tableau): eliminate coefficient clone in measurement"
