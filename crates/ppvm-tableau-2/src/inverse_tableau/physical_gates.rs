// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Physical gates, `C ← G·C`: they mix whole rows.

use ppvm_pauli_kernels::word_product;
use ppvm_traits_2::ApplyGateRule;

use super::InverseTableau;
use super::bitword::{BitWord, ZERO};

/// Where a gate sends one of its generators: `G† P G = (−1)^flip · P(factors)`.
///
/// The generators are X then Z for each gate qubit; `factors[j]` says whether
/// the Hermitian Pauli `P(factors)` contains generator `j`.
#[derive(Clone, Copy)]
struct Image<const M: usize> {
    factors: [bool; M],
    flip: bool,
}

impl<const M: usize> Image<M> {
    fn new((factors, flip): ([bool; M], bool)) -> Self {
        Self { factors, flip }
    }

    /// Whether generator `g` maps to itself up to sign, so its row keeps its bits.
    fn keeps_bits(&self, g: usize) -> bool {
        (0..M).all(|j| self.factors[j] == (j == g))
    }

    /// Number of Y factors: an X and a Z on the same qubit.
    fn y_count(&self) -> u32 {
        self.factors.chunks(2).filter(|xz| xz[0] && xz[1]).count() as u32
    }

    /// One word of the product of the old rows in `factors`, and its phase.
    #[inline(always)]
    fn multiply(&self, old_words: &[[BitWord; 2]; M]) -> ([BitWord; 2], u32) {
        let mut word = [ZERO; 2];
        let mut phase = 0u32;
        for j in (0..M).filter(|&j| self.factors[j]) {
            let (product, k) = word_product(word[0], word[1], old_words[j][0], old_words[j][1]);
            word = product;
            phase = phase.wrapping_add(k);
        }
        (word, phase)
    }
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
        for (g, image) in images
            .iter()
            .enumerate()
            .filter(|(g, image)| !image.keeps_bits(*g))
        {
            let ([new_x, new_z], phase) = image.multiply(&old_words);
            x[g][w] = new_x;
            z[g][w] = new_z;
            phases[g] = phases[g].wrapping_add(phase);
        }
    }
}

impl InverseTableau {
    /// Starting phase of each new row, as a power of `i`: one `i` per Y in the
    /// image, `−1` for its flip and `−1` per negative old row it multiplies.
    fn initial_phases<const M: usize>(&self, rows: [usize; M], images: &[Image<M>; M]) -> [u32; M] {
        images.map(|image| {
            let negatives = (0..M).filter(|&j| image.factors[j] && self.negative(rows[j]));
            image.y_count() + 2 * image.flip as u32 + 2 * negatives.count() as u32
        })
    }

    /// Set each row's sign from its final phase (`i⁰ = +`, `i² = −`).
    fn set_signs<const M: usize>(&mut self, rows: [usize; M], phases: [u32; M]) {
        for (g, &r) in rows.iter().enumerate() {
            debug_assert!(phases[g].is_multiple_of(2), "row is not Hermitian");
            if self.negative(r) != (phases[g] % 4 == 2) {
                self.flip_sign(r);
            }
        }
    }
}

/// Physical gates, `C ← G·C`: the new row for generator `P` is `C⁻¹ (G† P G) C`,
/// so the rule evaluated on `P` names the old rows to multiply.
impl ApplyGateRule for InverseTableau {
    type B = bool;

    fn apply1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        let images = [rule(true, false), rule(false, true)].map(Image::new);
        let rows = [q, self.n + q];
        let mut phases = self.initial_phases(rows, &images);
        let (x, z) = self.rows_mut(rows);
        recombine(x, z, &images, &mut phases);
        self.set_signs(rows, phases);
    }

    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        let generator = |j: usize| std::array::from_fn(|i| i == j);
        let images = std::array::from_fn(|j| Image::new(rule(generator(j))));
        let rows = [a, self.n + a, b, self.n + b];
        let mut phases = self.initial_phases(rows, &images);
        let (x, z) = self.rows_mut(rows);
        recombine(x, z, &images, &mut phases);
        self.set_signs(rows, phases);
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
