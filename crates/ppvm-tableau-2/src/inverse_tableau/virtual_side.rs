// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Virtual gates, i.e. changes of basis `C ← C·V`: they rewrite columns.

use super::InverseTableau;
use super::bitword::{ZERO, bit_position, mask_if};

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
}
