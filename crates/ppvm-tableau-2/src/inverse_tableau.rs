// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The inverse stabilizer tableau `C⁻¹` of a Clifford `C`.
//!
//! Row `ix_q` is `C⁻¹ X_q C` and row `iz_q` is `C⁻¹ Z_q C`, each stored as
//! packed X and Z bits plus one sign bit (rows are Hermitian, so their phase
//! is ±1). Rows are indexed `ix_q = q`, `iz_q = n + q`.
//!
//! Words are touched only through [`BitWord`], [`BITWORD_BITS`], [`ZERO`],
//! [`get_bit`], [`flip_bit`] and `Bits`, so a wider storage word can be tried
//! by changing the alias. Where rows sit in the arena is decided only by
//! `x_at` / `z_at`.

use std::ops::Range;

use ppvm_pauli_kernels::word_product;
use ppvm_traits_2::ApplyGateRule;

/// Storage word. Swap to benchmark other widths.
type BitWord = u64;
const BITWORD_BITS: usize = BitWord::BITS as usize;
const ZERO: BitWord = 0;

/// Bit `i` of a packed bit string.
#[inline(always)]
fn get_bit(words: &[BitWord], i: usize) -> bool {
    (words[i / BITWORD_BITS] >> (i % BITWORD_BITS)) & 1 == 1
}

#[inline(always)]
fn flip_bit(words: &mut [BitWord], i: usize) {
    words[i / BITWORD_BITS] ^= 1 << (i % BITWORD_BITS);
}

/// The inverse tableau of an `n`-qubit Clifford.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InverseTableau {
    n: usize,
    /// Words per row and plane: `n.div_ceil(BITWORD_BITS)`.
    stride: usize,
    /// X and Z words of all `2n` rows, arranged by `x_at` / `z_at`.
    words: Vec<BitWord>,
    /// Bit `r` is the sign of row `r` (set = −).
    signs: Vec<BitWord>,
}

// NOTE: we could impl Word for RowRef to make it a true PauliWord if we ever need it
/// A borrowed row: its X words, Z words and sign.
#[derive(Clone, Copy, Debug)]
pub struct RowRef<'a> {
    pub x: &'a [BitWord],
    pub z: &'a [BitWord],
    pub negative: bool,
}

impl RowRef<'_> {
    /// X bit of the row at qubit `q`.
    pub fn x_bit(&self, q: usize) -> bool {
        get_bit(self.x, q)
    }

    /// Z bit of the row at qubit `q`.
    pub fn z_bit(&self, q: usize) -> bool {
        get_bit(self.z, q)
    }
}

impl InverseTableau {
    /// The identity Clifford, i.e. the state `|0…0⟩`: `ix_q = X_q`, `iz_q = Z_q`.
    pub fn new(n: usize) -> Self {
        let stride = n.div_ceil(BITWORD_BITS);
        let mut t = Self {
            n,
            stride,
            words: vec![ZERO; 4 * n * stride],
            signs: vec![ZERO; (2 * n).div_ceil(BITWORD_BITS)],
        };
        for q in 0..n {
            let (ix, iz) = (t.x_at(q), t.z_at(n + q));
            flip_bit(&mut t.words[ix..], q);
            flip_bit(&mut t.words[iz..], q);
        }
        t
    }

    pub fn num_qubits(&self) -> usize {
        self.n
    }

    /// Row `ix_q = C⁻¹ X_q C`.
    pub fn ix(&self, q: usize) -> RowRef<'_> {
        debug_assert!(q < self.n);
        self.row(q)
    }

    /// Row `iz_q = C⁻¹ Z_q C`.
    pub fn iz(&self, q: usize) -> RowRef<'_> {
        debug_assert!(q < self.n);
        self.row(self.n + q)
    }

    fn row(&self, r: usize) -> RowRef<'_> {
        let (x, z) = (self.x_at(r), self.z_at(r));
        RowRef {
            x: &self.words[x..x + self.stride],
            z: &self.words[z..z + self.stride],
            negative: get_bit(&self.signs, r),
        }
    }

    fn flip_sign(&mut self, r: usize) {
        flip_bit(&mut self.signs, r);
    }

    /// Starting phase of each new row, as a power of `i`: one `i` per Y in the
    /// image, `−1` for its flip and `−1` per negative old row it multiplies.
    fn initial_phases<const M: usize>(&self, rows: [usize; M], images: &[Image<M>; M]) -> [u32; M] {
        images.map(|image| {
            let negatives = (0..M).filter(|&j| image.factors[j] && self.row(rows[j]).negative);
            image.y_count() + 2 * image.flip as u32 + 2 * negatives.count() as u32
        })
    }

    /// Set each row's sign from its final phase (`i⁰ = +`, `i² = −`).
    fn set_signs<const M: usize>(&mut self, rows: [usize; M], phases: [u32; M]) {
        for (g, &r) in rows.iter().enumerate() {
            debug_assert!(phases[g].is_multiple_of(2), "row is not Hermitian");
            if self.row(r).negative != (phases[g] % 4 == 2) {
                self.flip_sign(r);
            }
        }
    }

    /// Positions of row `r`'s X words and Z words in the arena.
    fn xz(&self, r: usize) -> [Range<usize>; 2] {
        let (x, z) = (self.x_at(r), self.z_at(r));
        [x..x + self.stride, z..z + self.stride]
    }

    // Arrangement: blocks by row kind, ix-X, ix-Z, iz-X, iz-Z (#204's order).
    #[inline(always)]
    fn x_at(&self, r: usize) -> usize {
        let (kind, i) = (r / self.n, r % self.n);
        ((2 * kind) * self.n + i) * self.stride
    }

    #[inline(always)]
    fn z_at(&self, r: usize) -> usize {
        let (kind, i) = (r / self.n, r % self.n);
        ((2 * kind + 1) * self.n + i) * self.stride
    }
}

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

/// Physical gates, `C ← G·C`: the new row for generator `P` is `C⁻¹ (G† P G) C`,
/// so the rule evaluated on `P` names the old rows to multiply.
impl ApplyGateRule for InverseTableau {
    type B = bool;

    fn apply1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        let images = [rule(true, false), rule(false, true)].map(Image::new);
        let rows = [q, self.n + q];
        let mut phases = self.initial_phases(rows, &images);
        let [[ix_x, ix_z], [iz_x, iz_z]] = rows.map(|r| self.xz(r));
        let [ix_x, ix_z, iz_x, iz_z] = self
            .words
            .get_disjoint_mut([ix_x, ix_z, iz_x, iz_z])
            .unwrap();
        recombine([ix_x, iz_x], [ix_z, iz_z], &images, &mut phases);
        self.set_signs(rows, phases);
    }

    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        let generator = |j: usize| std::array::from_fn(|i| i == j);
        let images = std::array::from_fn(|j| Image::new(rule(generator(j))));
        let rows = [a, self.n + a, b, self.n + b];
        let mut phases = self.initial_phases(rows, &images);
        let [[r0x, r0z], [r1x, r1z], [r2x, r2z], [r3x, r3z]] = rows.map(|r| self.xz(r));
        let [r0x, r0z, r1x, r1z, r2x, r2z, r3x, r3z] = self
            .words
            .get_disjoint_mut([r0x, r0z, r1x, r1z, r2x, r2z, r3x, r3z])
            .unwrap();
        recombine(
            [r0x, r1x, r2x, r3x],
            [r0z, r1z, r2z, r3z],
            &images,
            &mut phases,
        );
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
    use super::*;
    use ppvm_pauli_kernels::Bits;
    use ppvm_pauli_word_2::PauliWord;
    use ppvm_phased_pauli_word_2::Phased;
    use ppvm_traits_2::{Clifford, PauliBits, Phase};

    const SIZES: [usize; 7] = [1, 2, 7, 63, 64, 65, 130];

    /// Whether two rows anticommute: odd symplectic product.
    fn anticommute(a: RowRef, b: RowRef) -> bool {
        let mut parity = 0;
        for i in 0..a.x.len() {
            parity ^= ((a.x[i] & b.z[i]) ^ (a.z[i] & b.x[i])).popcount() & 1;
        }
        parity == 1
    }

    #[test]
    fn identity_rows_are_single_paulis() {
        for n in SIZES {
            let t = InverseTableau::new(n);
            for q in 0..n {
                for j in 0..n {
                    assert_eq!(t.ix(q).x_bit(j), j == q);
                    assert!(!t.ix(q).z_bit(j));
                    assert!(!t.iz(q).x_bit(j));
                    assert_eq!(t.iz(q).z_bit(j), j == q);
                }
                assert!(!t.ix(q).negative && !t.iz(q).negative);
            }
        }
    }

    /// Conjugation preserves commutation: `ix_a`, `iz_b` anticommute iff `a == b`.
    fn assert_canonical_commutation(t: &InverseTableau) {
        let n = t.num_qubits();
        for a in 0..n {
            for b in 0..n {
                assert_eq!(anticommute(t.ix(a), t.iz(b)), a == b);
                assert!(!anticommute(t.ix(a), t.ix(b)));
                assert!(!anticommute(t.iz(a), t.iz(b)));
            }
        }
    }

    #[test]
    fn identity_has_canonical_commutation() {
        for n in SIZES {
            assert_canonical_commutation(&InverseTableau::new(n));
        }
    }

    #[derive(Clone, Copy, Debug)]
    enum Gate {
        One(usize, usize),
        Two(usize, usize, usize),
    }

    fn apply<T: Clifford>(t: &mut T, gate: Gate) {
        match gate {
            Gate::One(g, q) => match g {
                0 => t.x(q),
                1 => t.y(q),
                2 => t.z(q),
                3 => t.h(q),
                4 => t.s(q),
                5 => t.s_dag(q),
                6 => t.sqrt_x(q),
                7 => t.sqrt_x_dag(q),
                8 => t.sqrt_y(q),
                _ => t.sqrt_y_dag(q),
            },
            Gate::Two(g, a, b) => match g {
                0 => t.cnot(a, b),
                1 => t.cz(a, b),
                _ => t.cy(a, b),
            },
        }
    }

    /// A random circuit from a xorshift seed.
    fn random_circuit(n: usize, len: usize, mut seed: u64) -> Vec<Gate> {
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as usize
        };
        (0..len)
            .map(|_| {
                if n < 2 || next() % 3 == 0 {
                    Gate::One(next() % 10, next() % n)
                } else {
                    let a = next() % n;
                    let b = (a + 1 + next() % (n - 1)) % n;
                    Gate::Two(next() % 3, a, b)
                }
            })
            .collect()
    }

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
