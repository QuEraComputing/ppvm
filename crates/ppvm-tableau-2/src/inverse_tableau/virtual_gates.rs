// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Virtual gates, i.e. changes of basis `C ← C·V`: they rewrite columns.

use ppvm_pauli_kernels::heisenberg;

use super::InverseTableau;
use super::bitword::{BitWord, ZERO, bit_position, mask_if};

impl InverseTableau {
    /// Change the basis of virtual qubit `v` by the one-qubit Clifford `V`
    /// whose Heisenberg rule is `rule`: `C ← C·V`.
    ///
    /// Updates the tableau only. The state `C (|φ⟩ ⊗ |0…0⟩)` stays the same if
    /// the caller also applies `V†` to `φ`, or if `V` acts trivially on it (e.g.
    /// `v` dormant and `V` fixes `|0⟩`).
    ///
    /// Every row `R` becomes `V† R V`: the rule runs on each row's bits in
    /// column `v`. Strided: one word per row.
    ///
    /// Used by: dormant-pivot collapse (`H` on the pivot, `X` to fix the
    /// outcome), active collapse (the final rotation onto `Z_h`), reset.
    pub fn change_basis1(&mut self, v: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        debug_assert!(v < self.n);
        let (w, m) = bit_position(v);
        for (r, (x, z)) in self.row_offsets().enumerate() {
            let (x, z) = (x + w, z + w);
            let old = [self.words[x] & m != ZERO, self.words[z] & m != ZERO];
            let (new, flip) = rule(old[0], old[1]);
            self.words[x] ^= m & mask_if(new[0] != old[0]);
            self.words[z] ^= m & mask_if(new[1] != old[1]);
            let (sw, sm) = bit_position(r);
            self.signs[sw] ^= sm & mask_if(flip);
        }
    }

    /// Change the basis of virtual qubits `a`, `b` by the two-qubit Clifford
    /// `V` whose Heisenberg rule is `rule` (`a` is the control for controlled
    /// gates): `C ← C·V`. Same contract as [`change_basis1`](Self::change_basis1).
    ///
    /// Used by: T-gate localization (`CX(p, v)` collects the X bits of `O` onto
    /// the pivot `p`), dormant-pivot collapse (`CX(p, v)`, `CZ(p, c)`), active
    /// collapse (the `CX` / `CZ` that reduce `P` to one qubit).
    pub fn change_basis2(
        &mut self,
        a: usize,
        b: usize,
        rule: impl Fn([bool; 4]) -> ([bool; 4], bool),
    ) {
        debug_assert!(a < self.n && b < self.n && a != b);
        let [(wa, ma), (wb, mb)] = [bit_position(a), bit_position(b)];
        for (r, (x, z)) in self.row_offsets().enumerate() {
            let words = [x + wa, z + wa, x + wb, z + wb];
            let masks = [ma, ma, mb, mb];
            let old: [bool; 4] = std::array::from_fn(|j| self.words[words[j]] & masks[j] != ZERO);
            let (new, flip) = rule(old);
            for j in 0..4 {
                self.words[words[j]] ^= masks[j] & mask_if(new[j] != old[j]);
            }
            let (sw, sm) = bit_position(r);
            self.signs[sw] ^= sm & mask_if(flip);
        }
    }

    /// Change basis by a pivot fan-out in one pass over the rows:
    /// `C ← C·V` with `V = CX(p, t₁)⋯CX(p, t_k) · CZ(p, c₁)⋯CZ(p, c_m) · L(p)`,
    /// where `L` is the one-qubit Clifford whose Heisenberg rule is `last`.
    ///
    /// Same result and contract as the corresponding [`change_basis2`] /
    /// [`change_basis1`] calls in that order, but each row loads and stores the
    /// pivot's bits and its sign once instead of once per gate. Targets must be
    /// distinct and differ from the pivot.
    ///
    /// Used by: T-gate localization (CX only, `last` the identity), both
    /// random measurement collapses (CX, CZ, then `H`/`S·H` and an `X`).
    ///
    /// [`change_basis1`]: Self::change_basis1
    /// [`change_basis2`]: Self::change_basis2
    pub fn fan_out(
        &mut self,
        pivot: usize,
        cx_targets: &[usize],
        cz_targets: &[usize],
        last: impl Fn(bool, bool) -> ([bool; 2], bool),
    ) {
        debug_assert!(
            cx_targets
                .iter()
                .chain(cz_targets)
                .all(|&t| t < self.n && t != pivot)
        );
        let (pw, pm) = bit_position(pivot);
        let cx = cx_targets
            .iter()
            .map(|&t| bit_position(t))
            .collect::<Vec<_>>();
        let cz = cz_targets
            .iter()
            .map(|&t| bit_position(t))
            .collect::<Vec<_>>();
        for (r, (x, z)) in self.row_offsets().enumerate() {
            let old_pivot = [
                self.words[x + pw] & pm != ZERO,
                self.words[z + pw] & pm != ZERO,
            ];
            let (mut pivot_bits, mut flip) = (old_pivot, false);
            for &(tw, tm) in &cx {
                flip ^= target_step(
                    &mut self.words,
                    [x + tw, z + tw],
                    tm,
                    &mut pivot_bits,
                    heisenberg::cnot,
                );
            }
            for &(tw, tm) in &cz {
                flip ^= target_step(
                    &mut self.words,
                    [x + tw, z + tw],
                    tm,
                    &mut pivot_bits,
                    heisenberg::cz,
                );
            }
            let (new_pivot, f) = last(pivot_bits[0], pivot_bits[1]);
            self.words[x + pw] ^= pm & mask_if(new_pivot[0] != old_pivot[0]);
            self.words[z + pw] ^= pm & mask_if(new_pivot[1] != old_pivot[1]);
            let (sw, sm) = bit_position(r);
            self.signs[sw] ^= sm & mask_if(flip ^ f);
        }
    }
}

/// Run a two-qubit `rule` on the pivot's bits (held in `pivot_bits`) and one
/// target's bits (at `words[at]`, bit `mask`); store the target's changed bits.
/// Returns the sign flip.
#[inline(always)]
fn target_step(
    words: &mut [BitWord],
    at: [usize; 2],
    mask: BitWord,
    pivot_bits: &mut [bool; 2],
    rule: impl Fn([bool; 4]) -> ([bool; 4], bool),
) -> bool {
    let old = at.map(|i| words[i] & mask != ZERO);
    let ([px, pz, tx, tz], flip) = rule([pivot_bits[0], pivot_bits[1], old[0], old[1]]);
    words[at[0]] ^= mask & mask_if(tx != old[0]);
    words[at[1]] ^= mask & mask_if(tz != old[1]);
    *pivot_bits = [px, pz];
    flip
}

#[cfg(test)]
mod tests {
    use super::super::test_support::*;
    use super::*;

    #[test]
    fn change_basis_preserves_commutation() {
        for (i, n) in [1, 2, 7, 65].into_iter().enumerate() {
            let mut t = InverseTableau::new(n);
            for gate in random_circuit(n, 100, 0xba5e + i as u64) {
                change_basis(&mut t, gate);
                assert_canonical_commutation(&t);
            }
        }
    }

    /// `(G·C)·V = G·(C·V)`: the order of physical gates relative to changes of
    /// basis does not matter.
    #[test]
    fn physical_gates_commute_with_changes_of_basis() {
        for (i, n) in SIZES.into_iter().filter(|&n| n > 1).enumerate() {
            let physical = random_circuit(n, 100, 0xc0de + i as u64);
            let basis = random_circuit(n, 100, 0xfeed + i as u64);
            let (mut first, mut second) = (InverseTableau::new(n), InverseTableau::new(n));
            physical.iter().for_each(|&g| apply(&mut first, g));
            basis.iter().for_each(|&g| change_basis(&mut first, g));
            basis.iter().for_each(|&g| change_basis(&mut second, g));
            physical.iter().for_each(|&g| apply(&mut second, g));
            assert_eq!(first, second, "n = {n}");
        }
    }

    type Rule1 = fn(bool, bool) -> ([bool; 2], bool);

    /// The fused fan-out equals the separate changes of basis it replaces.
    #[test]
    fn fan_out_matches_separate_changes_of_basis() {
        use ppvm_pauli_kernels::heisenberg as rules;
        let lasts: [Rule1; 6] = [
            |x, z| ([x, z], false),
            rules::h,
            rules::s,
            rules::x,
            rules::sqrt_x,
            rules::sqrt_y_dag,
        ];
        let mut seed = 0xfa40_u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as usize
        };
        for n in SIZES.into_iter().filter(|&n| n > 1) {
            for round in 0..20 {
                let mut start = InverseTableau::new(n);
                random_circuit(n, 60, next() as u64)
                    .into_iter()
                    .for_each(|g| apply(&mut start, g));
                random_circuit(n, 20, next() as u64)
                    .into_iter()
                    .for_each(|g| change_basis(&mut start, g));
                let pivot = next() % n;
                // Each other qubit is a CX target, a CZ target, or neither.
                let (mut cx, mut cz) = (Vec::new(), Vec::new());
                for t in (0..n).filter(|&t| t != pivot) {
                    match next() % 4 {
                        0 => cx.push(t),
                        1 => cz.push(t),
                        _ => {}
                    }
                }
                let last = lasts[round % lasts.len()];
                let mut fused = start.clone();
                fused.fan_out(pivot, &cx, &cz, last);
                let mut separate = start;
                cx.iter()
                    .for_each(|&t| separate.change_basis2(pivot, t, rules::cnot));
                cz.iter()
                    .for_each(|&t| separate.change_basis2(pivot, t, rules::cz));
                separate.change_basis1(pivot, last);
                assert_eq!(fused, separate, "n = {n}, round {round}");
            }
        }
    }
}
