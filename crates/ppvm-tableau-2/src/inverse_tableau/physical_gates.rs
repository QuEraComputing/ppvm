// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Physical gates, `C ← G·C`: they mix whole rows.

use ppvm_pauli_kernels::{row_product, word_product};
use ppvm_traits_2::{ApplyGateRule, CliffordBatch};

use super::InverseTableau;
use super::bitword::{BitWord, ZERO};

/// Where a gate sends one of its generators: `G† P G = (−1)^flip · P(factors)`.
///
/// The generators are X then Z for each gate qubit; bit `j` of `factors` says
/// whether the Hermitian Pauli `P(factors)` contains generator `j`. Plain
/// integer masks, so that for a given gate everything below folds to constants.
#[derive(Clone, Copy)]
struct Image<const M: usize> {
    factors: u8,
    flip: bool,
}

impl<const M: usize> Image<M> {
    fn new((bits, flip): ([bool; M], bool)) -> Self {
        let mut factors = 0;
        for (j, &b) in bits.iter().enumerate() {
            factors |= u8::from(b) << j;
        }
        Self { factors, flip }
    }

    fn has(&self, j: usize) -> bool {
        self.factors >> j & 1 == 1
    }

    /// Whether generator `g` maps to itself up to sign, so its row keeps its bits.
    fn keeps_bits(&self, g: usize) -> bool {
        self.factors == 1 << g
    }

    /// Whether row `g` becomes itself times rows that keep their bits (`changed`
    /// marks the rows that do not), so it can be multiplied in place (CX, CZ, S,
    /// √X and their inverses).
    fn multiplies_in_place(&self, g: usize, changed: u8) -> bool {
        self.has(g) && self.factors & !(1 << g) & changed == 0
    }

    /// Number of Y factors: an X and a Z on the same qubit.
    fn y_count(&self) -> u32 {
        (self.factors & self.factors >> 1 & 0b0101_0101).count_ones()
    }

    /// One word of the product of the old rows in `factors`, and its phase.
    #[inline(always)]
    fn multiply(&self, old_words: &[[BitWord; 2]; M]) -> ([BitWord; 2], u32) {
        let mut word = [ZERO; 2];
        let mut phase = 0u32;
        for (j, old) in old_words.iter().enumerate() {
            if self.has(j) {
                let (product, k) = word_product(word[0], word[1], old[0], old[1]);
                word = product;
                phase = phase.wrapping_add(k);
            }
        }
        (word, phase)
    }
}

/// Bit `g` set for each row the gate changes.
fn changed_rows<const M: usize>(images: &[Image<M>; M]) -> u8 {
    let mut changed = 0;
    for (g, image) in images.iter().enumerate() {
        changed |= u8::from(!image.keeps_bits(g)) << g;
    }
    changed
}

/// Rewrite row `g` (X words `x[g]`, Z words `z[g]`) as `C⁻¹ (G† P_g G) C`.
///
/// Write the image as `(−1)^flip · i^#Y` times its generators in order (a Y is
/// `i·X·Z`). Conjugating by `C` turns each generator into its old row, so the
/// new row is the product of those old rows times `(−1)^flip · i^#Y`. For
/// example `S` sends `X ↦ −Y = −i·X·Z`, so `ix ← −i · ix · iz`, and keeps `Z`,
/// so `iz` stays.
///
/// The product is taken site by site, so the new words at index `w` depend
/// only on the old words at `w`: each index is loaded, multiplied in registers
/// and stored once. `phases[g]` (a power of `i`) comes in with the image's own
/// factors and signs and collects each word's product phase.
#[inline(always)]
fn recombine<const M: usize>(
    x: [&mut [BitWord]; M],
    z: [&mut [BitWord]; M],
    images: &[Image<M>; M],
    phases: &mut [u32; M],
) {
    let len = x[0].len();
    assert!(x.iter().chain(&z).all(|s| s.len() == len));
    for w in 0..len {
        let old_words: [[BitWord; 2]; M] = std::array::from_fn(|j| [x[j][w], z[j][w]]);
        for (g, image) in images.iter().enumerate() {
            if image.keeps_bits(g) {
                continue;
            }
            let ([new_x, new_z], phase) = image.multiply(&old_words);
            x[g][w] = new_x;
            z[g][w] = new_z;
            phases[g] = phases[g].wrapping_add(phase);
        }
    }
}

/// `(&mut v[g], &v[j])` for `g ≠ j`.
fn pair<T>(v: &mut [T], g: usize, j: usize) -> (&mut T, &T) {
    if g < j {
        let (lo, hi) = v.split_at_mut(j);
        (&mut lo[g], &hi[0])
    } else {
        let (lo, hi) = v.split_at_mut(g);
        (&mut hi[0], &lo[j])
    }
}

/// The in-place path: each changed row `g` is multiplied by its other factors
/// with `row_product`, like Stim's `zs[t] *= zs[c]`.
///
/// That computes `R_g · R_j …` with `R_g` first, while the image is a product
/// in generator order; moving `R_g` past an earlier factor costs `−1` when the
/// two anticommute, i.e. when they are the X and Z rows of one qubit.
#[inline(always)]
fn multiply_in_place<const M: usize>(
    mut x: [&mut [BitWord]; M],
    mut z: [&mut [BitWord]; M],
    images: &[Image<M>; M],
    phases: &mut [u32; M],
) {
    for (g, image) in images.iter().enumerate() {
        if image.keeps_bits(g) {
            continue;
        }
        for j in 0..M {
            if j == g || !image.has(j) {
                continue;
            }
            let (dst_x, src_x) = pair(&mut x, g, j);
            let (dst_z, src_z) = pair(&mut z, g, j);
            let k = row_product(dst_x, dst_z, src_x, src_z);
            let reorder = if j < g && j / 2 == g / 2 { 2 } else { 0 };
            phases[g] = phases[g].wrapping_add(u32::from(k) + reorder);
        }
    }
}

impl InverseTableau {
    /// Starting phase of each new row, as a power of `i`: one `i` per Y in the
    /// image, `−1` for its flip and `−1` per negative old row it multiplies.
    /// A row whose image contains the row itself leaves its own sign out: its
    /// phase is then the *change* of its sign (see `set_signs`).
    #[inline(always)]
    fn initial_phases<const M: usize>(&self, rows: [usize; M], images: &[Image<M>; M]) -> [u32; M] {
        let mut phases = [0; M];
        for (g, image) in images.iter().enumerate() {
            phases[g] = image.y_count() + 2 * u32::from(image.flip);
            for (j, &r) in rows.iter().enumerate() {
                if j != g && image.has(j) && self.negative(r) {
                    phases[g] += 2;
                }
            }
        }
        phases
    }

    /// Set each row's sign from its final phase (`i⁰ = +`, `i² = −`): relative
    /// to its old sign when its image contains the row itself, else absolute.
    #[inline(always)]
    fn set_signs<const M: usize>(
        &mut self,
        rows: [usize; M],
        images: &[Image<M>; M],
        phases: [u32; M],
    ) {
        for (g, &r) in rows.iter().enumerate() {
            debug_assert!(phases[g].is_multiple_of(2), "row is not Hermitian");
            let negative = phases[g] % 4 == 2;
            if images[g].has(g) {
                if negative {
                    self.flip_sign(r);
                }
            } else if self.negative(r) != negative {
                self.flip_sign(r);
            }
        }
    }

    /// Replace each row of `qubits` by its image's product of old rows: in place
    /// when the pattern allows it, else with the fused pass. The choice depends
    /// only on the images, which are constants for a given gate.
    #[inline(always)]
    fn rewrite_rows<const M: usize>(&mut self, qubits: &[usize], images: &[Image<M>; M]) {
        let rows = self.qubit_rows(qubits);
        let mut phases = self.initial_phases(rows, images);
        let in_place = {
            let changed = changed_rows(images);
            let mut all = true;
            for (g, image) in images.iter().enumerate() {
                all &= changed >> g & 1 == 0 || image.multiplies_in_place(g, changed);
            }
            all
        };
        let (x, z) = self.qubit_rows_mut(qubits);
        if in_place {
            multiply_in_place(x, z, images, &mut phases);
        } else {
            recombine(x, z, images, &mut phases);
        }
        self.set_signs(rows, images, phases);
    }
}

/// Physical gates, `C ← G·C`: the new row for generator `P` is `C⁻¹ (G† P G) C`,
/// so the rule evaluated on `P` names the old rows to multiply.
impl ApplyGateRule for InverseTableau {
    type B = bool;

    fn apply1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        let images = [rule(true, false), rule(false, true)].map(Image::new);
        self.rewrite_rows(&[q], &images);
    }

    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        let images = [
            rule([true, false, false, false]),
            rule([false, true, false, false]),
            rule([false, false, true, false]),
            rule([false, false, false, true]),
        ]
        .map(Image::new);
        self.rewrite_rows(&[a, b], &images);
    }

    fn apply_sign1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        let ((x_bits, x_flip), (z_bits, z_flip)) = (rule(true, false), rule(false, true));
        debug_assert!(x_bits == [true, false] && z_bits == [false, true]);
        if x_flip {
            self.flip_sign(q);
        }
        if z_flip {
            self.flip_sign(self.n + q);
        }
    }
}

/// The default loops are right here: gates on different qubits touch disjoint
/// rows, so a batch has no shared work to fuse.
impl CliffordBatch for InverseTableau {}

#[cfg(test)]
mod tests {
    use ppvm_pauli_word_2::PauliWord;
    use ppvm_phased_pauli_word_2::Phased;
    use ppvm_traits_2::{PauliBits, Phase};

    use super::super::RowRef;
    use super::super::test_support::*;
    use super::*;

    type Reference = Phased<PauliWord<[u64; 4]>>;

    /// `C⁻¹ P C` for `C = G_m ⋯ G_1`: propagate `P` through the gates in reverse.
    fn heisenberg(n: usize, q: usize, z: bool, circuit: &[Gate]) -> Reference {
        let mut word = PauliWord::new(n);
        if z {
            word.set_z_bit(q, true);
        } else {
            word.set_x_bit(q, true);
        }
        let mut p = Phased::new(word);
        for &gate in circuit.iter().rev() {
            apply(&mut p, gate);
        }
        p
    }

    fn assert_row_eq(row: RowRef, expected: &Reference, n: usize) {
        for j in 0..n {
            assert_eq!(row.x_bit(j), expected.word().x_bit(j), "x bit {j}");
            assert_eq!(row.z_bit(j), expected.word().z_bit(j), "z bit {j}");
        }
        let negative = match expected.phase() {
            Phase::Pos1 => false,
            Phase::Neg1 => true,
            phase => panic!("non-Hermitian reference phase {phase:?}"),
        };
        assert_eq!(row.negative, negative, "sign");
    }

    /// Rows must equal `C⁻¹ X_q C` / `C⁻¹ Z_q C` computed by Heisenberg
    /// propagation of single Pauli words.
    #[test]
    fn physical_gates_match_heisenberg_propagation() {
        for (i, n) in SIZES.into_iter().enumerate() {
            let circuit = random_circuit(n, 200, 0x9e37_79b9 + i as u64);
            let mut t = InverseTableau::new(n);
            for &gate in &circuit {
                apply(&mut t, gate);
            }
            for q in 0..n {
                assert_row_eq(t.ix(q), &heisenberg(n, q, false, &circuit), n);
                assert_row_eq(t.iz(q), &heisenberg(n, q, true, &circuit), n);
            }
        }
    }

    #[test]
    fn physical_gates_preserve_commutation() {
        for (i, n) in [1, 2, 7, 65].into_iter().enumerate() {
            let mut t = InverseTableau::new(n);
            for gate in random_circuit(n, 100, 0x5eed + i as u64) {
                apply(&mut t, gate);
                assert_canonical_commutation(&t);
            }
        }
    }
}
