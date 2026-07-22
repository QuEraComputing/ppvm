# PR 180 Translation-Symmetry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split PR 180's translation-symmetry implementation into focused modules and correct group validation, traversal, momentum projection, and sector checking without changing existing public import paths.

**Architecture:** Keep `ppvm_pauli_sum::symmetry` as the public facade. Put group construction and a shared odometer traversal in `group.rs`, unnormalized real merging in `merge.rs`, and exact-character momentum operations in `momentum.rs`; `mod.rs` re-exports the existing API. Treat each declared generator order as exact while allowing relations between different generators, and use exact modular character numerators to handle stabilizers.

**Tech Stack:** Rust 2024, `PauliWord`, `PauliSum<T: Config>`, `fxhash`, `num::Complex`, Cargo workspace tests.

---

## File map

- Delete: `crates/ppvm-pauli-sum/src/symmetry.rs`
- Create: `crates/ppvm-pauli-sum/src/symmetry/mod.rs`
- Create: `crates/ppvm-pauli-sum/src/symmetry/group.rs`
- Create: `crates/ppvm-pauli-sum/src/symmetry/merge.rs`
- Create: `crates/ppvm-pauli-sum/src/symmetry/momentum.rs`
- Create: `crates/ppvm-pauli-sum/src/symmetry/tests.rs`
- Create: `crates/ppvm-pauli-sum/tests/symmetry_api.rs`
- Preserve: `crates/ppvm-pauli-sum/src/lib.rs` (`pub mod symmetry;` remains unchanged)

Every new Rust file starts with the repository's SPDX header.

## Task 1: Split the module without changing behavior

**Files:**

- Delete: `crates/ppvm-pauli-sum/src/symmetry.rs`
- Create: `crates/ppvm-pauli-sum/src/symmetry/mod.rs`
- Create: `crates/ppvm-pauli-sum/src/symmetry/group.rs`
- Create: `crates/ppvm-pauli-sum/src/symmetry/merge.rs`
- Create: `crates/ppvm-pauli-sum/src/symmetry/momentum.rs`
- Create: `crates/ppvm-pauli-sum/src/symmetry/tests.rs`
- Create: `crates/ppvm-pauli-sum/tests/symmetry_api.rs`

- [ ] **Step 1: Establish the behavioral baseline**

Run:

```bash
cargo test -p ppvm-pauli-sum symmetry::tests
```

Expected: all 15 existing symmetry unit tests pass.

- [ ] **Step 2: Add an external public-surface smoke test**

Create `crates/ppvm-pauli-sum/tests/symmetry_api.rs`:

```rust
// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::Complex;
use ppvm_pauli_sum::symmetry::{
    TranslationGroup, canonicalize_pauli_sum, canonicalize_pauli_sum_complex,
    check_momentum_sector,
};
use ppvm_pauli_word::word::PauliWord;

type W = PauliWord<[u8; 1], fxhash::FxBuildHasher, true>;

#[test]
fn public_symmetry_imports_remain_available() {
    let group = TranslationGroup::chain_1d(2);

    let mut real_basis: Vec<W> = vec![W::from("XI"), W::from("IX")];
    let mut real_coeffs = vec![1.0, 1.0];
    canonicalize_pauli_sum(&mut real_basis, &mut real_coeffs, &group);
    assert_eq!(real_basis.len(), 1);

    let mut complex_basis: Vec<W> = vec![W::from("ZI"), W::from("IZ")];
    let mut complex_coeffs = vec![Complex::new(1.0, 0.0); 2];
    assert!(check_momentum_sector(
        &complex_basis,
        &complex_coeffs,
        &group,
        &[0],
        1e-12,
    )
    .is_ok());
    canonicalize_pauli_sum_complex(
        &mut complex_basis,
        &mut complex_coeffs,
        &group,
        &[0],
    );
    assert_eq!(complex_basis.len(), 1);
}
```

- [ ] **Step 3: Move symbols to their semantic owners**

Use this exact routing, making no body or documentation changes yet:

| Destination | Symbols/content |
| --- | --- |
| `mod.rs` | Existing module-level docs, module declarations, public re-exports |
| `group.rs` | `TranslationGroup` and its complete current inherent implementation |
| `merge.rs` | `canonicalize_pauli_sum`, `symmetry_merge_pauli_sum` |
| `momentum.rs` | `character` inherent impl, `canonicalize_pauli_sum_complex`, `check_momentum_sector`, `SectorCheckError` and formatting impls |
| `tests.rs` | Existing `tests` module contents without the outer `mod tests { ... }` wrapper |

Start `tests.rs` with the imports that were previously inherited from the
monolithic parent module:

```rust
use super::*;
use fxhash::FxHashMap;
use num::Complex;
use ppvm_pauli_word::word::PauliWord;
use std::f64::consts::PI;
```

`mod.rs` must end with this facade:

```rust
mod group;
mod merge;
mod momentum;

pub use group::TranslationGroup;
pub use merge::{canonicalize_pauli_sum, symmetry_merge_pauli_sum};
pub use momentum::{
    SectorCheckError, canonicalize_pauli_sum_complex, check_momentum_sector,
};

#[cfg(test)]
mod tests;
```

Move `character` out of the main `TranslationGroup` impl into this additional
impl in `momentum.rs`:

```rust
impl TranslationGroup {
    pub fn character(&self, k_modes: &[i32], counter: &[u32]) -> Complex<f64> {
        // Move the existing body verbatim in this task.
    }
}
```

- [ ] **Step 4: Verify the behavior-neutral split**

Run:

```bash
cargo fmt --all
cargo test -p ppvm-pauli-sum symmetry
cargo test -p ppvm-pauli-sum --test symmetry_api
```

Expected: the 15 unit tests and the public-surface integration test pass.

- [ ] **Step 5: Commit the structural split**

```bash
git add crates/ppvm-pauli-sum/src/symmetry.rs \
  crates/ppvm-pauli-sum/src/symmetry \
  crates/ppvm-pauli-sum/tests/symmetry_api.rs
git commit -m "refactor(pauli-sum): split translation symmetry module"
```

## Task 2: Enforce the translation-group construction contract

**Files:**

- Modify: `crates/ppvm-pauli-sum/src/symmetry/group.rs`
- Modify: `crates/ppvm-pauli-sum/src/symmetry/tests.rs`

- [ ] **Step 1: Add failing constructor-validation tests**

Append these tests to `tests.rs`:

```rust
#[test]
#[should_panic(expected = "generator 0 order must be nonzero")]
fn rejects_zero_generator_order() {
    TranslationGroup::from_generators(2, vec![vec![1, 0]], vec![0]);
}

#[test]
#[should_panic(expected = "declared order 4 != exact permutation order 2")]
fn rejects_inflated_generator_order() {
    TranslationGroup::from_generators(2, vec![vec![1, 0]], vec![4]);
}

#[test]
#[should_panic(expected = "generators 0 and 1 do not commute")]
fn rejects_noncommuting_generators() {
    let swap_01 = vec![1, 0, 2];
    let swap_12 = vec![0, 2, 1];
    TranslationGroup::from_generators(3, vec![swap_01, swap_12], vec![2, 2]);
}

#[test]
fn rejects_zero_lattice_dimensions() {
    assert!(std::panic::catch_unwind(|| TranslationGroup::chain_1d(0)).is_err());
    assert!(std::panic::catch_unwind(|| TranslationGroup::torus_2d(0, 2)).is_err());
    assert!(std::panic::catch_unwind(|| TranslationGroup::torus_3d(2, 0, 2)).is_err());
    assert!(std::panic::catch_unwind(|| TranslationGroup::ladder(2, 0)).is_err());
}

#[test]
fn rejects_dimension_product_overflow_before_allocation() {
    assert!(
        std::panic::catch_unwind(|| TranslationGroup::torus_2d(usize::MAX, 2)).is_err()
    );
    assert!(
        std::panic::catch_unwind(|| TranslationGroup::ladder(usize::MAX, 2)).is_err()
    );
}

#[test]
fn rejects_group_order_overflow() {
    let orders = if usize::BITS == 64 {
        vec![u32::MAX, u32::MAX, u32::MAX]
    } else {
        vec![u32::MAX, u32::MAX]
    };
    assert!(std::panic::catch_unwind(|| {
        super::group::checked_group_order(&orders)
    })
    .is_err());
}
```

- [ ] **Step 2: Run the new tests and confirm the old behavior fails**

Run:

```bash
cargo test -p ppvm-pauli-sum symmetry::tests::rejects_
```

Expected: at least the zero-order, inflated-order, and noncommuting tests fail.

- [ ] **Step 3: Add checked arithmetic and exact permutation-order helpers**

Add these private helpers to `group.rs`:

```rust
fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn checked_lcm(a: usize, b: usize, context: &str) -> usize {
    a.checked_div(gcd(a, b))
        .and_then(|q| q.checked_mul(b))
        .unwrap_or_else(|| panic!("{context} overflow"))
}

fn permutation_order(perm: &[u32], generator: usize) -> u32 {
    let mut seen = vec![false; perm.len()];
    let mut order = 1usize;
    for start in 0..perm.len() {
        if seen[start] {
            continue;
        }
        let mut length = 0usize;
        let mut q = start;
        loop {
            assert!(!seen[q], "generator {generator} contains a malformed cycle");
            seen[q] = true;
            length += 1;
            q = perm[q] as usize;
            if q == start {
                break;
            }
        }
        order = checked_lcm(order, length, "permutation order");
    }
    u32::try_from(order).unwrap_or_else(|_| {
        panic!("generator {generator} exact permutation order does not fit in u32")
    })
}

fn permutations_commute(left: &[u32], right: &[u32]) -> bool {
    (0..left.len()).all(|q| {
        left[right[q] as usize] == right[left[q] as usize]
    })
}

pub(super) fn checked_group_order(orders: &[u32]) -> usize {
    orders.iter().enumerate().fold(1usize, |acc, (g, &value)| {
        acc.checked_mul(value as usize)
            .unwrap_or_else(|| panic!("group order overflows usize at generator {g}"))
    })
}
```

- [ ] **Step 4: Cache validated group metadata**

Extend the private fields:

```rust
pub struct TranslationGroup {
    n_qubits: usize,
    perms: Vec<Vec<u32>>,
    orders: Vec<u32>,
    order: usize,
    phase_modulus: usize,
}
```

After validating each permutation's length/range/uniqueness, validate its exact
order and pairwise commutativity. Compute metadata with checked folds:

```rust
for (g, &declared) in orders.iter().enumerate() {
    assert!(declared != 0, "generator {g} order must be nonzero");
    let exact = permutation_order(&perms[g], g);
    assert_eq!(
        declared, exact,
        "generator {g} declared order {declared} != exact permutation order {exact}",
    );
}
for left in 0..perms.len() {
    for right in left + 1..perms.len() {
        assert!(
            permutations_commute(&perms[left], &perms[right]),
            "generators {left} and {right} do not commute",
        );
    }
}
let order = checked_group_order(&orders);
let phase_modulus = orders.iter().fold(1usize, |acc, &value| {
    checked_lcm(acc, value as usize, "character phase modulus")
});
```

Return cached `self.order` from `order()`.

- [ ] **Step 5: Make lattice constructors fail before arithmetic/allocation**

At the start of each constructor, assert every dimension is positive. Convert
each generator order with `u32::try_from`, form site counts with `checked_mul`,
and convert every generated target index with `u32::try_from`. Use messages that
name the constructor and offending dimension/product.

- [ ] **Step 6: Run validation tests and the crate suite**

```bash
cargo fmt --all
cargo test -p ppvm-pauli-sum symmetry::tests::rejects_
cargo test -p ppvm-pauli-sum symmetry
```

Expected: all tests pass.

- [ ] **Step 7: Commit group validation**

```bash
git add crates/ppvm-pauli-sum/src/symmetry/group.rs \
  crates/ppvm-pauli-sum/src/symmetry/tests.rs
git commit -m "fix(pauli-sum): validate translation groups"
```

## Task 3: Replace duplicated reconstruction with one odometer traversal

**Files:**

- Modify: `crates/ppvm-pauli-sum/src/symmetry/group.rs`
- Modify: `crates/ppvm-pauli-sum/src/symmetry/tests.rs`

- [ ] **Step 1: Add failing traversal and release-check tests**

```rust
#[test]
fn odometer_yields_expected_counter_order() {
    let group = TranslationGroup::torus_2d(2, 3);
    let counters: Vec<Vec<u32>> = group
        .orbit_with_counters(&word("XIIIII"))
        .map(|(_, counter)| counter)
        .collect();
    assert_eq!(
        counters,
        vec![
            vec![0, 0], vec![1, 0], vec![0, 1],
            vec![1, 1], vec![0, 2], vec![1, 2],
        ],
    );
}

#[test]
fn traversal_matches_brute_force_composition() {
    let group = TranslationGroup::torus_2d(2, 3);
    let source = word("XYZIII");
    for (candidate, counter) in group.orbit_with_counters(&source) {
        let mut brute = source;
        for (g, &count) in counter.iter().enumerate() {
            for _ in 0..count {
                brute = group.apply_generator(&brute, g);
            }
        }
        assert_eq!(candidate, brute);
    }
}

#[test]
fn public_word_width_checks_are_not_debug_only() {
    let group = TranslationGroup::chain_1d(4);
    let short = word("XI");
    assert!(std::panic::catch_unwind(|| group.canonicalize(&short)).is_err());
    assert!(std::panic::catch_unwind(|| group.canonicalize_with_shift(&short)).is_err());
    assert!(std::panic::catch_unwind(|| group.orbit(&short)).is_err());
}
```

- [ ] **Step 2: Run the tests and confirm the traversal API is absent**

```bash
cargo test -p ppvm-pauli-sum symmetry::tests::odometer_yields_expected_counter_order
```

Expected: compilation fails because `orbit_with_counters` does not exist.

- [ ] **Step 3: Implement the shared iterator**

Add this iterator to `group.rs` and make `apply_generator` `pub(super)`:

```rust
pub(super) struct GroupOrbit<'a, A, S, const R: bool>
where
    A: PauliStorage,
{
    group: &'a TranslationGroup,
    current: PauliWord<A, S, R>,
    counter: Vec<u32>,
    remaining: usize,
}

impl<A, S, const R: bool> Iterator for GroupOrbit<'_, A, S, R>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    type Item = (PauliWord<A, S, R>, Vec<u32>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }
        let item = (self.current, self.counter.clone());
        self.remaining -= 1;
        if self.remaining == 0 {
            return Some(item);
        }
        for g in 0..self.group.orders.len() {
            if self.group.orders[g] == 1 {
                continue;
            }
            self.current = self.group.apply_generator(&self.current, g);
            self.counter[g] += 1;
            if self.counter[g] < self.group.orders[g] {
                break;
            }
            self.counter[g] = 0;
        }
        Some(item)
    }
}
```

Add the constructor, asserting width before returning the lazy iterator:

```rust
pub(super) fn orbit_with_counters<'a, A, S, const R: bool>(
    &'a self,
    word: &'a PauliWord<A, S, R>,
) -> GroupOrbit<'a, A, S, R>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    assert_eq!(word.n_qubits(), self.n_qubits, "word and group must agree on n_qubits");
    GroupOrbit {
        group: self,
        current: *word,
        counter: vec![0; self.orders.len()],
        remaining: self.order,
    }
}
```

- [ ] **Step 4: Rewrite all three consumers**

Use these bodies so all consumers share traversal and preserve the first
counter on equal representatives:

```rust
// canonicalize
let mut traversal = self.orbit_with_counters(word);
let (mut best, _) = traversal.next().expect("a finite group contains the identity");
for (candidate, _) in traversal {
    if candidate < best {
        best = candidate;
    }
}
best

// canonicalize_with_shift
let mut traversal = self.orbit_with_counters(word);
let (mut best, mut counter_from_word) =
    traversal.next().expect("a finite group contains the identity");
for (candidate, counter) in traversal {
    if candidate < best {
        best = candidate;
        counter_from_word = counter;
    }
}
let counter_to_word = counter_from_word
    .iter()
    .zip(self.orders.iter())
    .map(|(&counter, &order)| (order - counter) % order)
    .collect();
(best, counter_to_word)

// orbit
self.orbit_with_counters(word)
    .map(|(candidate, _)| candidate)
```

- [ ] **Step 5: Run traversal, symmetry, and public API tests**

```bash
cargo fmt --all
cargo test -p ppvm-pauli-sum symmetry::tests::odometer_
cargo test -p ppvm-pauli-sum symmetry::tests::traversal_matches_
cargo test -p ppvm-pauli-sum symmetry
cargo test -p ppvm-pauli-sum --test symmetry_api
```

Expected: all tests pass.

- [ ] **Step 6: Commit traversal centralization**

```bash
git add crates/ppvm-pauli-sum/src/symmetry/group.rs \
  crates/ppvm-pauli-sum/src/symmetry/tests.rs
git commit -m "refactor(pauli-sum): centralize symmetry traversal"
```

## Task 4: Add exact character arithmetic

**Files:**

- Modify: `crates/ppvm-pauli-sum/src/symmetry/group.rs`
- Modify: `crates/ppvm-pauli-sum/src/symmetry/momentum.rs`
- Modify: `crates/ppvm-pauli-sum/src/symmetry/tests.rs`

- [ ] **Step 1: Add exact-character tests**

```rust
#[test]
fn character_numerator_normalizes_negative_modes() {
    let group = TranslationGroup::chain_1d(4);
    assert_eq!(group.character_numerator(&[-1], &[1]), 3);
    assert_eq!(group.character_numerator(&[3], &[1]), 3);
    assert!((group.character(&[-1], &[1]) - Complex::new(0.0, -1.0)).norm() < 1e-12);
}

#[test]
fn exact_character_detects_cross_generator_kernel() {
    let swap = vec![1, 0];
    let group = TranslationGroup::from_generators(
        2,
        vec![swap.clone(), swap],
        vec![2, 2],
    );
    assert_ne!(group.character_numerator(&[1, 0], &[1, 1]), 0);
    assert_eq!(group.character_numerator(&[1, 1], &[1, 1]), 0);
}

#[test]
fn character_checks_slice_lengths_in_release_builds() {
    let group = TranslationGroup::chain_1d(4);
    assert!(std::panic::catch_unwind(|| group.character(&[], &[0])).is_err());
    assert!(std::panic::catch_unwind(|| group.character(&[0], &[])).is_err());
}
```

- [ ] **Step 2: Run the tests and confirm the exact helper is absent**

```bash
cargo test -p ppvm-pauli-sum symmetry::tests::character_numerator_
```

Expected: compilation fails because `character_numerator` does not exist.

- [ ] **Step 3: Implement exact numerator calculation**

Add a `pub(super)` getter for `phase_modulus` if momentum code needs it, then
implement this method in the `TranslationGroup` impl in `momentum.rs`:

```rust
pub(super) fn character_numerator(&self, k_modes: &[i32], counter: &[u32]) -> usize {
    assert_eq!(k_modes.len(), self.n_generators(), "k_modes length mismatch");
    assert_eq!(counter.len(), self.n_generators(), "counter length mismatch");
    let modulus = self.phase_modulus() as u128;
    let mut numerator = 0u128;
    for g in 0..self.n_generators() {
        let order = self.generator_order(g);
        let k = (k_modes[g] as i64).rem_euclid(order as i64) as u128;
        let count = (counter[g] % order) as u128;
        let reduced = (k * count) % order as u128;
        let factor = self.phase_modulus() as u128 / order as u128;
        numerator = (numerator + reduced * factor) % modulus;
    }
    numerator as usize
}
```

Rewrite `character` to convert only this exact numerator:

```rust
let numerator = self.character_numerator(k_modes, counter);
let phase = 2.0 * PI * numerator as f64 / self.phase_modulus() as f64;
Complex::from_polar(1.0, phase)
```

- [ ] **Step 4: Run exact-character and full symmetry tests**

```bash
cargo fmt --all
cargo test -p ppvm-pauli-sum symmetry::tests::character_
cargo test -p ppvm-pauli-sum symmetry
```

Expected: all tests pass.

- [ ] **Step 5: Commit exact character arithmetic**

```bash
git add crates/ppvm-pauli-sum/src/symmetry/group.rs \
  crates/ppvm-pauli-sum/src/symmetry/momentum.rs \
  crates/ppvm-pauli-sum/src/symmetry/tests.rs
git commit -m "fix(pauli-sum): compute symmetry characters exactly"
```

## Task 5: Correct momentum projection for stabilizers

**Files:**

- Modify: `crates/ppvm-pauli-sum/src/symmetry/momentum.rs`
- Modify: `crates/ppvm-pauli-sum/src/symmetry/tests.rs`

- [ ] **Step 1: Add projection regressions**

```rust
#[test]
fn period_two_k_zero_round_trip_preserves_rep_coefficient() {
    let group = TranslationGroup::chain_1d(4);
    let mut basis = vec![word("XIXI"), word("IXIX")];
    let mut coeffs = vec![Complex::new(1.0, 0.0); 2];
    canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &group, &[0]);
    assert_eq!(basis.len(), 1);
    assert!((coeffs[0] - Complex::new(1.0, 0.0)).norm() < 1e-12);
}

#[test]
fn period_two_compatible_k_two_round_trip_preserves_rep_coefficient() {
    let group = TranslationGroup::chain_1d(4);
    let rep = group.canonicalize(&word("XIXI"));
    let mut members: FxHashMap<W, Complex<f64>> = FxHashMap::default();
    for (member, counter) in group.orbit_with_counters(&rep) {
        members.entry(member).or_insert_with(|| group.character(&[2], &counter).conj());
    }
    let (mut basis, mut coeffs): (Vec<W>, Vec<Complex<f64>>) = members.into_iter().unzip();
    canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &group, &[2]);
    assert_eq!(basis, vec![rep]);
    assert!((coeffs[0] - Complex::new(1.0, 0.0)).norm() < 1e-12);
}

#[test]
fn incompatible_stabilizer_projects_orbit_to_zero() {
    let group = TranslationGroup::chain_1d(4);
    let mut basis = vec![word("XXXX")];
    let mut coeffs = vec![Complex::new(1.0, 0.0)];
    canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &group, &[1]);
    assert!(basis.is_empty());
    assert!(coeffs.is_empty());
}

#[test]
fn partial_period_two_orbit_is_averaged_with_missing_member_zero() {
    let group = TranslationGroup::chain_1d(4);
    let mut basis = vec![word("XIXI")];
    let mut coeffs = vec![Complex::new(1.0, 0.0)];
    canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &group, &[0]);
    assert_eq!(basis.len(), 1);
    assert!((coeffs[0] - Complex::new(0.5, 0.0)).norm() < 1e-12);
}
```

- [ ] **Step 2: Run the regressions and verify current normalization fails**

```bash
cargo test -p ppvm-pauli-sum symmetry::tests::period_two_
cargo test -p ppvm-pauli-sum symmetry::tests::incompatible_stabilizer_projects_
```

Expected: the period-two and incompatible-stabilizer tests fail.

- [ ] **Step 3: Implement one stabilizer-aware projection per represented orbit**

Replace the fixed `inv_g` loop with this data flow:

```rust
let mut input: FxHashMap<PauliWord<A, S, R>, Complex<f64>> = FxHashMap::default();
for (word, &coeff) in basis.iter().zip(coeffs.iter()) {
    *input.entry(*word).or_insert(Complex::new(0.0, 0.0)) += coeff;
}
let reps: FxHashSet<_> = input.keys().map(|word| group.canonicalize(word)).collect();
let mut projected = FxHashMap::default();

for rep in reps {
    let mut members: FxHashMap<_, (Vec<u32>, usize)> = FxHashMap::default();
    let mut compatible = true;
    for (member, counter) in group.orbit_with_counters(&rep) {
        let numerator = group.character_numerator(k_modes, &counter);
        match members.entry(member) {
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert((counter, numerator));
            }
            std::collections::hash_map::Entry::Occupied(entry) => {
                if entry.get().1 != numerator {
                    compatible = false;
                    break;
                }
            }
        }
    }
    if !compatible {
        continue;
    }
    let orbit_size = members.len() as f64;
    let mut rep_coeff = Complex::new(0.0, 0.0);
    for (member, (counter, _)) in members {
        let coeff = input
            .get(&member)
            .copied()
            .unwrap_or(Complex::new(0.0, 0.0));
        rep_coeff += group.character(k_modes, &counter) * coeff / orbit_size;
    }
    projected.insert(rep, rep_coeff);
}
```

Use `fxhash::FxHashSet`. Replace `basis` and `coeffs` from `projected` exactly
as the current function does from `merged`. Update the rustdoc to call `k=0`
an orbit average rather than a plain merge. Rename
`momentum_zero_complex_merge_matches_real_merge` to
`momentum_zero_complex_projection_is_orbit_average` and retain its assertions
that the real result is `10.0` while the complex result is `2.5`.

- [ ] **Step 4: Run projection and existing momentum tests**

```bash
cargo fmt --all
cargo test -p ppvm-pauli-sum symmetry::tests::period_two_
cargo test -p ppvm-pauli-sum symmetry::tests::incompatible_stabilizer_projects_
cargo test -p ppvm-pauli-sum symmetry::tests::momentum_
```

Expected: all tests pass; the renamed `k=0` test still demonstrates real sum
`10.0` versus complex average `2.5`.

- [ ] **Step 5: Commit projection correctness**

```bash
git add crates/ppvm-pauli-sum/src/symmetry/momentum.rs \
  crates/ppvm-pauli-sum/src/symmetry/tests.rs
git commit -m "fix(pauli-sum): normalize momentum projection by orbit"
```

## Task 6: Make sector validation complete and diagnostic

**Files:**

- Modify: `crates/ppvm-pauli-sum/src/symmetry/momentum.rs`
- Modify: `crates/ppvm-pauli-sum/src/symmetry/tests.rs`

- [ ] **Step 1: Add validator regressions**

```rust
#[test]
fn sector_check_rejects_missing_orbit_members() {
    let group = TranslationGroup::chain_1d(4);
    let basis = vec![word("ZIII")];
    let coeffs = vec![Complex::new(1.0, 0.0)];
    assert!(matches!(
        check_momentum_sector(&basis, &coeffs, &group, &[0], 1e-12),
        Err(SectorCheckError::CoefficientMismatch { .. })
    ));
}

#[test]
fn sector_check_rejects_incompatible_stabilizer() {
    let group = TranslationGroup::chain_1d(4);
    let basis = vec![word("XXXX")];
    let coeffs = vec![Complex::new(1.0, 0.0)];
    assert!(matches!(
        check_momentum_sector(&basis, &coeffs, &group, &[1], 1e-12),
        Err(SectorCheckError::IncompatibleStabilizer { .. })
    ));
}

#[test]
fn sector_check_rejects_invalid_numeric_inputs() {
    let group = TranslationGroup::chain_1d(2);
    let basis = vec![word("ZI")];
    assert!(matches!(
        check_momentum_sector(&basis, &[Complex::new(1.0, 0.0)], &group, &[0], f64::NAN),
        Err(SectorCheckError::InvalidTolerance { .. })
    ));
    assert!(matches!(
        check_momentum_sector(&basis, &[Complex::new(f64::NAN, 0.0)], &group, &[0], 1e-12),
        Err(SectorCheckError::NonFiniteCoefficient { .. })
    ));
}

#[test]
fn sector_error_display_names_the_words() {
    let group = TranslationGroup::chain_1d(2);
    let basis = vec![word("ZI")];
    let coeffs = vec![Complex::new(1.0, 0.0)];
    let message = check_momentum_sector(&basis, &coeffs, &group, &[0], 1e-12)
        .unwrap_err()
        .to_string();
    assert!(message.contains("ZI") || message.contains("IZ"));
}
```

- [ ] **Step 2: Run the regressions and confirm current validation accepts bad input**

```bash
cargo test -p ppvm-pauli-sum symmetry::tests::sector_check_rejects_
```

Expected: the missing-member, incompatible-stabilizer, and NaN tests fail.

- [ ] **Step 3: Replace the error struct with explicit variants**

```rust
pub enum SectorCheckError<A: PauliStorage, S, const R: bool> {
    InvalidTolerance { tol: f64 },
    NonFiniteCoefficient {
        pauli: PauliWord<A, S, R>,
        coeff: Complex<f64>,
    },
    CoefficientMismatch {
        rep: PauliWord<A, S, R>,
        offending_pauli: PauliWord<A, S, R>,
        expected: Complex<f64>,
        actual: Complex<f64>,
        shift: Vec<u32>,
    },
    IncompatibleStabilizer {
        rep: PauliWord<A, S, R>,
        shift: Vec<u32>,
    },
}
```

Implement `Debug` and `Display` with
`S: BuildHasher + Clone + Default + HashFinalize`; format `rep`, `pauli`, and
`offending_pauli` with `{}` rather than placeholders.

- [ ] **Step 4: Validate the complete represented orbit**

At function entry, return `InvalidTolerance` unless `tol.is_finite() && tol >=
0.0`, and return `NonFiniteCoefficient` for any non-finite component. Coalesce
duplicate words and remove exactly-zero totals. Build the set of canonical
representatives. For each representative:

1. enumerate `orbit_with_counters` into a map from distinct word to its first
   `(counter, exact_numerator)`;
2. return `IncompatibleStabilizer` if the same word reappears with a different
   numerator;
3. select the first present nonzero member and infer
   `rep_coeff = character(counter) * actual`;
4. for every distinct member compute
   `expected = character(counter).conj() * rep_coeff`;
5. load absent members as `Complex::new(0.0, 0.0)` and return
   `CoefficientMismatch` when
   `(actual - expected).norm() > tol * rep_coeff.norm().max(1.0)`.

Keep the existing assertions for basis/coeff and momentum-mode lengths.
Implement the core loop as follows after those assertions:

```rust
if !tol.is_finite() || tol < 0.0 {
    return Err(SectorCheckError::InvalidTolerance { tol });
}
let mut input: FxHashMap<PauliWord<A, S, R>, Complex<f64>> = FxHashMap::default();
for (pauli, &coeff) in basis.iter().zip(coeffs.iter()) {
    if !coeff.re.is_finite() || !coeff.im.is_finite() {
        return Err(SectorCheckError::NonFiniteCoefficient {
            pauli: *pauli,
            coeff,
        });
    }
    *input.entry(*pauli).or_insert(Complex::new(0.0, 0.0)) += coeff;
}
input.retain(|_, coeff| *coeff != Complex::new(0.0, 0.0));
let reps: FxHashSet<_> = input.keys().map(|pauli| group.canonicalize(pauli)).collect();

for rep in reps {
    let mut members: FxHashMap<_, (Vec<u32>, usize)> = FxHashMap::default();
    for (member, counter) in group.orbit_with_counters(&rep) {
        let numerator = group.character_numerator(k_modes, &counter);
        match members.entry(member) {
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert((counter, numerator));
            }
            std::collections::hash_map::Entry::Occupied(entry) => {
                if entry.get().1 != numerator {
                    return Err(SectorCheckError::IncompatibleStabilizer {
                        rep,
                        shift: counter,
                    });
                }
            }
        }
    }

    let (reference_word, (reference_counter, _)) = members
        .iter()
        .find(|(member, _)| input.contains_key(*member))
        .expect("represented orbit has a nonzero member");
    let rep_coeff = group.character(k_modes, reference_counter) * input[reference_word];

    for (member, (counter, _)) in members {
        let expected = group.character(k_modes, &counter).conj() * rep_coeff;
        let actual = input
            .get(&member)
            .copied()
            .unwrap_or(Complex::new(0.0, 0.0));
        if (actual - expected).norm() > tol * rep_coeff.norm().max(1.0) {
            return Err(SectorCheckError::CoefficientMismatch {
                rep,
                offending_pauli: member,
                expected,
                actual,
                shift: counter,
            });
        }
    }
}
Ok(())
```

- [ ] **Step 5: Run validator, projection, and public API tests**

```bash
cargo fmt --all
cargo test -p ppvm-pauli-sum symmetry::tests::sector_
cargo test -p ppvm-pauli-sum symmetry::tests::momentum_
cargo test -p ppvm-pauli-sum --test symmetry_api
```

Expected: all tests pass, including the pre-existing valid `k=1` eigenstate.

- [ ] **Step 6: Commit complete validation**

```bash
git add crates/ppvm-pauli-sum/src/symmetry/momentum.rs \
  crates/ppvm-pauli-sum/src/symmetry/tests.rs
git commit -m "fix(pauli-sum): validate complete momentum sectors"
```

## Task 7: Finish documentation and verify the branch

**Files:**

- Modify: `crates/ppvm-pauli-sum/src/symmetry/mod.rs`
- Modify: `crates/ppvm-pauli-sum/src/symmetry/group.rs`
- Modify: `crates/ppvm-pauli-sum/src/symmetry/momentum.rs`
- Verify: `crates/ppvm-pauli-sum/tests/symmetry_api.rs`

- [ ] **Step 1: Align rustdoc with the implemented contracts**

Update module and item docs to state all of the following explicitly:

- generator orders are exact, but the combined action may have a kernel;
- `order()` is the abstract product order and `orbit()` may repeat words;
- `canonicalize_with_shift` returns the first valid counter and counters are not
  unique in the presence of stabilizers;
- complex `k=0` projection averages a distinct orbit while real merging sums;
- incompatible stabilizer sectors project to zero; and
- `check_momentum_sector` treats absent orbit members as zero.

- [ ] **Step 2: Run formatting, lint, tests, docs, and wasm verification**

Run each command and stop on the first failure:

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo doc --workspace --no-deps
cargo build --target wasm32-unknown-unknown --workspace \
  --exclude ppvm-python-native --exclude ppvm-cli --exclude ppvm-tui
```

Expected: every command exits successfully with no warnings or failed tests.

- [ ] **Step 3: Inspect the final public diff**

```bash
git diff origin/split/1-translation-symmetry...HEAD --stat
git diff origin/split/1-translation-symmetry...HEAD -- crates/ppvm-pauli-sum/src/lib.rs
git status --short
```

Expected: `lib.rs` is unchanged, the old `symmetry.rs` is replaced by the five
semantic module files, the integration test is present, and the worktree is
clean except for intentional documentation changes.

- [ ] **Step 4: Commit final rustdoc changes if Step 1 produced a diff**

```bash
git add crates/ppvm-pauli-sum/src/symmetry
git commit -m "docs(pauli-sum): clarify symmetry projection semantics"
```

If Step 1 was already completed in earlier implementation commits and `git
status --short` is empty, do not create an empty commit.

## Task 8: Prepare the stacked-PR handoff

**Files:** None on this branch.

- [ ] **Step 1: Record the descendant rebase order in the handoff**

Use this exact order after the implementation branch is merged into
`split/1-translation-symmetry`:

```text
split/2-ctpp-core
split/3-symmetric-evolution
split/4-autotune-ledgers
kossakowski-dissipator (also based on split/3-symmetric-evolution)
```

- [ ] **Step 2: Identify PR 182 follow-up edits without applying them here**

After rebasing PR 182, update its fixed-`1/|G|` comments in
`crates/ppvm-python-native/src/interface.rs`, preserve the Rust re-export paths,
and add a Python period-two regression to
`ppvm-python/test/test_momentum_merge.py`. Run:

```bash
uv run --project ppvm-python --group dev pytest ppvm-python/test/test_momentum_merge.py
```

Expected: the rebased Python momentum tests pass. Do not copy those binding or
Python changes into this PR 180 implementation branch.

- [ ] **Step 3: Prepare review-thread responses**

Group duplicate Copilot discussions under the commits that fix generator
validation, zero dimensions, shared traversal, runtime shape checks, `k=0`
documentation, stabilizer normalization, and error diagnostics. Resolve a
thread only after its regression test and implementation are both present.
