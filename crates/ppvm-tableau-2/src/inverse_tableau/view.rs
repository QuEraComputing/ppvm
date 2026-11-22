// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`VirtualView`]: the part of the tableau a virtual state may use.

use ppvm_pauli_kernels::row_product;

use super::bitword::{BITWORD_BITS, BitWord, ZERO, bit_position, get_bit, mask_if};
use super::{InverseTableau, RowRef};

/// A set of qubits as packed bits, e.g. the active virtual qubits, so rows
/// can be classified against it word by word.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QubitMask {
    pub(super) words: Vec<BitWord>,
}

impl QubitMask {
    /// The empty set over `n` qubits.
    pub fn new(n: usize) -> Self {
        Self {
            words: vec![ZERO; n.div_ceil(BITWORD_BITS)],
        }
    }

    pub fn contains(&self, q: usize) -> bool {
        get_bit(&self.words, q)
    }

    pub fn set(&mut self, q: usize, member: bool) {
        let (w, m) = bit_position(q);
        self.words[w] = self.words[w] & !m | m & mask_if(member);
    }
}

/// An owned row, e.g. a product of rows: its X words, Z words and sign.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OwnedRow {
    x: Vec<BitWord>,
    z: Vec<BitWord>,
    negative: bool,
}

impl OwnedRow {
    pub fn row(&self) -> RowRef<'_> {
        RowRef {
            x: &self.x,
            z: &self.z,
            negative: self.negative,
        }
    }
}

/// The virtual side of an [`InverseTableau`]: row reads and changes of basis
/// `C ← C·V`, but no physical gates. A virtual state gets this view, so it can
/// change basis together with its own state while `|ψ⟩` stays the same.
pub struct VirtualView<'a> {
    t: &'a mut InverseTableau,
}

impl InverseTableau {
    pub fn virtual_view(&mut self) -> VirtualView<'_> {
        VirtualView { t: self }
    }
}

impl VirtualView<'_> {
    pub fn num_qubits(&self) -> usize {
        self.t.num_qubits()
    }

    /// Row `ix_q = C⁻¹ X_q C`.
    pub fn ix(&self, q: usize) -> RowRef<'_> {
        self.t.ix(q)
    }

    /// Row `iz_q = C⁻¹ Z_q C`.
    pub fn iz(&self, q: usize) -> RowRef<'_> {
        self.t.iz(q)
    }

    /// `O = C⁻¹ (∏ Z_q) C = ∏ iz_q` over `qubits`, which must be distinct. The
    /// `Z_q` commute, so `O` is a Hermitian Pauli string with sign ±1.
    pub fn iz_product(&self, qubits: &[usize]) -> OwnedRow {
        let stride = self.t.stride;
        let (mut x, mut z) = (vec![ZERO; stride], vec![ZERO; stride]);
        // Power of `i` in front of the product of the bare rows.
        let mut phase = 0u8;
        for &q in qubits {
            let r = self.t.iz(q);
            let k = row_product(&mut x, &mut z, r.x, r.z) + 2 * u8::from(r.negative);
            phase = phase.wrapping_add(k);
        }
        debug_assert_eq!(phase % 2, 0, "product of commuting rows is Hermitian");
        OwnedRow {
            x,
            z,
            negative: phase % 4 == 2,
        }
    }

    /// See [`InverseTableau::change_basis1`].
    pub fn change_basis1(&mut self, v: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        self.t.change_basis1(v, rule);
    }

    /// See [`InverseTableau::change_basis2`].
    pub fn change_basis2(
        &mut self,
        a: usize,
        b: usize,
        rule: impl Fn([bool; 4]) -> ([bool; 4], bool),
    ) {
        self.t.change_basis2(a, b, rule);
    }

    /// See [`InverseTableau::fan_out`].
    pub fn fan_out(
        &mut self,
        pivot: usize,
        cx_targets: &[usize],
        cz_targets: &[usize],
        last: impl Fn(bool, bool) -> ([bool; 2], bool),
    ) {
        self.t.fan_out(pivot, cx_targets, cz_targets, last);
    }

    /// Measure `Z_q` with every virtual qubit dormant, i.e. on `C |0…0⟩`. See
    /// [`InverseTableau::measure`]; the collapse is a change of basis only.
    pub fn measure(&mut self, q: usize, coin: impl FnOnce() -> bool) -> bool {
        self.t.measure(q, coin)
    }

    /// See [`InverseTableau::collapse_dormant`].
    pub fn collapse_dormant(&mut self, q: usize, active: &QubitMask, m: bool) {
        self.t.collapse_dormant(q, active, m);
    }
}
