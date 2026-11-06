// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Phase-tracking Clifford conjugation. Gates touching lost sites are no-ops.
//! The base word must implement both `PauliBits` and `LossInfo`.

use ppvm_traits_2::{Clifford, PauliBits, Phase, loss::LossInfo};

use crate::data::Phased;

impl<W: PauliBits + LossInfo> Phased<W> {
    #[inline]
    fn flip_sign_if(&mut self, cond: bool) {
        if cond {
            self.phase *= Phase::Neg1;
        }
    }
}

impl<W: PauliBits + LossInfo> Clifford for Phased<W> {
    #[inline]
    fn x(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        self.flip_sign_if(self.word.z_bit(q));
    }

    #[inline]
    fn y(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        self.flip_sign_if(self.word.x_bit(q) ^ self.word.z_bit(q));
    }

    #[inline]
    fn z(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        self.flip_sign_if(self.word.x_bit(q));
    }

    #[inline]
    fn h(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        let x = self.word.x_bit(q);
        let z = self.word.z_bit(q);
        self.word.set_xz_bits(q, z, x);
        self.flip_sign_if(x && z);
    }

    // Heisenberg convention: S† X S = -Y.
    #[inline]
    fn s(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        let x = self.word.x_bit(q);
        let z = self.word.z_bit(q);
        self.word.set_z_bit(q, z ^ x);
        self.flip_sign_if(x && !z);
    }

    #[inline]
    fn cnot(&mut self, ctrl: usize, tgt: usize) {
        debug_assert_ne!(ctrl, tgt, "CNOT requires distinct qubits");
        debug_assert!(
            ctrl < self.word.n_sites() && tgt < self.word.n_sites(),
            "qubit out of bounds"
        );
        if self.word.is_lost(ctrl) || self.word.is_lost(tgt) {
            return;
        }
        let xc = self.word.x_bit(ctrl);
        let zc = self.word.z_bit(ctrl);
        let xt = self.word.x_bit(tgt);
        let zt = self.word.z_bit(tgt);
        self.word.set_x_bit_and_z_bit(tgt, xt ^ xc, ctrl, zc ^ zt);
        self.flip_sign_if(xc && zt && (xt == zc));
    }

    #[inline]
    fn cz(&mut self, a: usize, b: usize) {
        debug_assert_ne!(a, b, "CZ requires distinct qubits");
        debug_assert!(
            a < self.word.n_sites() && b < self.word.n_sites(),
            "qubit out of bounds"
        );
        if self.word.is_lost(a) || self.word.is_lost(b) {
            return;
        }
        let xa = self.word.x_bit(a);
        let za = self.word.z_bit(a);
        let xb = self.word.x_bit(b);
        let zb = self.word.z_bit(b);
        self.word.set_z_bit_pair(a, za ^ xb, b, zb ^ xa);
        self.flip_sign_if(xa && xb && (za ^ zb));
    }

    #[inline]
    fn s_dag(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        let x = self.word.x_bit(q);
        let z = self.word.z_bit(q);
        self.word.set_z_bit(q, z ^ x);
        self.flip_sign_if(x && z);
    }

    #[inline]
    fn sqrt_x(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        let x = self.word.x_bit(q);
        let z = self.word.z_bit(q);
        self.word.set_x_bit(q, x ^ z);
        self.flip_sign_if(x && z);
    }

    #[inline]
    fn sqrt_x_dag(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        let x = self.word.x_bit(q);
        let z = self.word.z_bit(q);
        self.word.set_x_bit(q, x ^ z);
        self.flip_sign_if(!x && z);
    }

    #[inline]
    fn sqrt_y(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        let x = self.word.x_bit(q);
        let z = self.word.z_bit(q);
        self.word.set_xz_bits(q, z, x);
        self.flip_sign_if(!x && z);
    }

    #[inline]
    fn sqrt_y_dag(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        let x = self.word.x_bit(q);
        let z = self.word.z_bit(q);
        self.word.set_xz_bits(q, z, x);
        self.flip_sign_if(x && !z);
    }

    #[inline]
    fn cy(&mut self, control: usize, target: usize) {
        debug_assert_ne!(control, target, "CY requires distinct qubits");
        debug_assert!(
            control < self.word.n_sites() && target < self.word.n_sites(),
            "qubit out of bounds"
        );
        if self.word.is_lost(control) || self.word.is_lost(target) {
            return;
        }
        let xc = self.word.x_bit(control);
        let zc = self.word.z_bit(control);
        let xt = self.word.x_bit(target);
        let zt = self.word.z_bit(target);
        self.word
            .set_xz_bits2(control, xc, zc ^ xt ^ zt, target, xt ^ xc, zt ^ xc);
        self.flip_sign_if(xc && (xt ^ zt) && !(zc ^ zt));
    }
}

#[cfg(test)]
mod tests {
    use crate::{Phased, PhasedPauliWord};
    use ppvm_lossy_pauli_word_2::LossyPauliWord;
    use ppvm_traits_2::Clifford;

    type TestWord = Phased<LossyPauliWord>;

    #[test]
    fn gates_touching_lost_sites_preserve_word_and_phase() {
        use ppvm_traits_2::Phase;

        let single: [fn(&mut TestWord, usize); 10] = [
            TestWord::x,
            TestWord::y,
            TestWord::z,
            TestWord::h,
            TestWord::s,
            TestWord::s_dag,
            TestWord::sqrt_x,
            TestWord::sqrt_x_dag,
            TestWord::sqrt_y,
            TestWord::sqrt_y_dag,
        ];
        let pairs: [fn(&mut TestWord, usize, usize); 3] =
            [TestWord::cnot, TestWord::cz, TestWord::cy];
        for phase in [Phase::Pos1, Phase::PosI, Phase::Neg1, Phase::NegI] {
            for input in ["XL", "LY", "LL"] {
                let original = Phased::with_phase(LossyPauliWord::from(input), phase);
                let lost = if input.starts_with('L') { 0 } else { 1 };
                for gate in single {
                    let mut word = original.clone();
                    gate(&mut word, lost);
                    assert_eq!(word, original);
                }
                for gate in pairs {
                    for (a, b) in [(0, 1), (1, 0)] {
                        let mut word = original.clone();
                        gate(&mut word, a, b);
                        assert_eq!(word, original);
                    }
                }
            }
        }
    }

    fn conj(input: &str, gate: impl Fn(&mut PhasedPauliWord)) -> String {
        let mut w = PhasedPauliWord::from(input);
        gate(&mut w);
        w.to_string()
    }

    #[test]
    fn ordinary_word_reports_no_loss() {
        use crate::PauliWord;
        use ppvm_traits_2::loss::LossInfo;

        let word: PauliWord = PauliWord::from("XYZ");
        assert_eq!(word.loss_weight(), 0);
        for i in 0..3 {
            assert!(!word.is_lost(i));
        }
    }

    #[test]
    fn single_qubit_gates_track_sign() {
        // (gate, input → target) ported from the old phase/clifford.rs tests.
        let x = [("+I", "+I"), ("+X", "+X"), ("+Y", "-Y"), ("+Z", "-Z")];
        let y = [("+I", "+I"), ("+X", "-X"), ("+Y", "+Y"), ("+Z", "-Z")];
        let z = [("+I", "+I"), ("+X", "-X"), ("+Y", "-Y"), ("+Z", "+Z")];
        let h = [("+I", "+I"), ("+X", "+Z"), ("+Y", "-Y"), ("+Z", "+X")];
        let s = [("+I", "+I"), ("+X", "-Y"), ("+Y", "+X"), ("+Z", "+Z")];
        for (i, t) in x {
            assert_eq!(conj(i, |w| w.x(0)), t, "X {i}");
        }
        for (i, t) in y {
            assert_eq!(conj(i, |w| w.y(0)), t, "Y {i}");
        }
        for (i, t) in z {
            assert_eq!(conj(i, |w| w.z(0)), t, "Z {i}");
        }
        for (i, t) in h {
            assert_eq!(conj(i, |w| w.h(0)), t, "H {i}");
        }
        for (i, t) in s {
            assert_eq!(conj(i, |w| w.s(0)), t, "S {i}");
        }
    }

    #[test]
    fn cnot_tracks_sign() {
        for (i, t) in [
            ("+II", "+II"),
            ("+IX", "+IX"),
            ("+IZ", "+ZZ"),
            ("+IY", "+ZY"),
            ("+XI", "+XX"),
            ("+XX", "+XI"),
            ("+XY", "+YZ"),
            ("+XZ", "-YY"),
            ("+ZI", "+ZI"),
            ("+ZX", "+ZX"),
            ("+ZY", "+IY"),
            ("+ZZ", "+IZ"),
            ("+YI", "+YX"),
            ("+YX", "+YI"),
            ("+YY", "-XZ"),
            ("+YZ", "+XY"),
        ] {
            assert_eq!(conj(i, |w| w.cnot(0, 1)), t, "CNOT {i}");
        }
    }

    #[test]
    fn cz_tracks_sign() {
        for (i, t) in [
            ("+II", "+II"),
            ("+IX", "+ZX"),
            ("+IY", "+ZY"),
            ("+IZ", "+IZ"),
            ("+XI", "+XZ"),
            ("+XX", "+YY"),
            ("+XY", "-YX"),
            ("+XZ", "+XI"),
            ("+ZI", "+ZI"),
            ("+ZX", "+IX"),
            ("+ZY", "+IY"),
            ("+ZZ", "+ZZ"),
            ("+YI", "+YZ"),
            ("+YX", "-XY"),
            ("+YY", "+XX"),
            ("+YZ", "+YI"),
        ] {
            assert_eq!(conj(i, |w| w.cz(0, 1)), t, "CZ {i}");
        }
    }

    #[test]
    fn phase_carries_through_gate() {
        // A pre-existing +i factor is preserved under a signed conjugation:
        // H (−iY) = −i·(−Y) = +iY.
        assert_eq!(conj("-iY", |w| w.h(0)), "+iY");
    }

    #[test]
    fn extension_gates_track_sign() {
        // The single-qubit table of `ppvm-traits/src/traits/clifford.rs`.
        let s_dag = [("+I", "+I"), ("+X", "+Y"), ("+Y", "-X"), ("+Z", "+Z")];
        let sqrt_x = [("+I", "+I"), ("+X", "+X"), ("+Y", "-Z"), ("+Z", "+Y")];
        let sqrt_x_dag = [("+I", "+I"), ("+X", "+X"), ("+Y", "+Z"), ("+Z", "-Y")];
        let sqrt_y = [("+I", "+I"), ("+X", "+Z"), ("+Y", "+Y"), ("+Z", "-X")];
        let sqrt_y_dag = [("+I", "+I"), ("+X", "-Z"), ("+Y", "+Y"), ("+Z", "+X")];
        for (i, t) in s_dag {
            assert_eq!(conj(i, |w| w.s_dag(0)), t, "S† {i}");
        }
        for (i, t) in sqrt_x {
            assert_eq!(conj(i, |w| w.sqrt_x(0)), t, "√X {i}");
        }
        for (i, t) in sqrt_x_dag {
            assert_eq!(conj(i, |w| w.sqrt_x_dag(0)), t, "√X† {i}");
        }
        for (i, t) in sqrt_y {
            assert_eq!(conj(i, |w| w.sqrt_y(0)), t, "√Y {i}");
        }
        for (i, t) in sqrt_y_dag {
            assert_eq!(conj(i, |w| w.sqrt_y_dag(0)), t, "√Y† {i}");
        }
    }

    #[test]
    fn cy_tracks_sign() {
        // The 16-entry CY table of `ppvm-traits/src/traits/clifford.rs`.
        for (i, t) in [
            ("+II", "+II"),
            ("+IX", "+ZX"),
            ("+IY", "+IY"),
            ("+IZ", "+ZZ"),
            ("+XI", "+XY"),
            ("+XX", "-YZ"),
            ("+XY", "+XI"),
            ("+XZ", "+YX"),
            ("+YI", "+YY"),
            ("+YX", "+XZ"),
            ("+YY", "+YI"),
            ("+YZ", "-XX"),
            ("+ZI", "+ZI"),
            ("+ZX", "+IX"),
            ("+ZY", "+ZY"),
            ("+ZZ", "+IZ"),
        ] {
            assert_eq!(conj(i, |w| w.cy(0, 1)), t, "CY {i}");
            assert_eq!(conj(i, |w| w.zcy(0, 1)), t, "ZCY {i}");
        }
    }
}
