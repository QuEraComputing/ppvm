// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Phase-tracking Clifford conjugation. Gates touching lost sites are no-ops.
//!
//! `Phased<W>` implements `ApplyGateRule`, so `Clifford` comes from the blanket
//! impl. The bit update is delegated to the word's own `ApplyGateRule`, so each
//! word type keeps its fastest write path; `Phased` only multiplies its `ℤ₄`
//! phase by −1 when the rule (from `ppvm_pauli_kernels::heisenberg`,
//! `S† X S = −Y`) reports a sign flip.

use std::cell::Cell;

use ppvm_traits_2::{ApplyGateRule, PauliBits, Phase, loss::LossInfo};

use crate::data::Phased;

impl<W> Phased<W> {
    /// Branch-free: the flip is a coin toss on random Paulis, so a branch here
    /// mispredicts often (measured ~1 ns per gate on arm64).
    #[inline]
    fn flip_sign_if(&mut self, cond: bool) {
        self.phase = Phase::from_exponent(self.phase.exponent() + 2 * cond as u8);
    }
}

/// A word that skips a gate (e.g. on a lost site) never calls the rule, so the
/// recorded flip stays `false` and the phase is untouched. This relies on `W`
/// calling the rule once per applied gate, as `ApplyGateRule` requires.
impl<W> ApplyGateRule for Phased<W>
where
    W: ApplyGateRule<B = bool> + PauliBits + LossInfo,
{
    type B = bool;

    #[inline]
    fn apply1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        let flip = Cell::new(false);
        self.word.apply1(q, |x, z| {
            let out = rule(x, z);
            flip.set(out.1);
            out
        });
        self.flip_sign_if(flip.get());
    }

    #[inline]
    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        let flip = Cell::new(false);
        self.word.apply2(a, b, |bits| {
            let out = rule(bits);
            flip.set(out.1);
            out
        });
        self.flip_sign_if(flip.get());
    }

    /// Words may skip sign-only rules entirely, so the flip is computed here.
    #[inline]
    fn apply_sign1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        if self.word.is_lost(q) {
            return;
        }
        let (_, flip) = rule(self.word.x_bit(q), self.word.z_bit(q));
        self.flip_sign_if(flip);
    }
}

#[cfg(test)]
mod tests {
    use crate::{Phased, PhasedPauliWord};
    use ppvm_lossy_pauli_word_2::LossyPauliWord;
    use ppvm_traits_2::Clifford;

    type TestWord = Phased<LossyPauliWord>;

    /// Both word types under `Phased` agree on random circuits, phase included.
    #[test]
    fn word_types_agree_on_bits_and_phase() {
        let mut seed = 0xD1B5_4A32_D192_ED03u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for _ in 0..64 {
            let start: String = (0..64)
                .map(|_| ['I', 'X', 'Y', 'Z'][(next() % 4) as usize])
                .collect();
            let mut plain = PhasedPauliWord::new(start.as_str().into());
            let mut lossy = TestWord::new(start.as_str().into());
            for _ in 0..256 {
                let (g, a) = (next() % 13, (next() % 64) as usize);
                let b = (a + 1 + (next() % 63) as usize) % 64;
                for w in [&mut plain as &mut dyn Clifford, &mut lossy] {
                    match g {
                        0 => w.x(a),
                        1 => w.y(a),
                        2 => w.z(a),
                        3 => w.h(a),
                        4 => w.s(a),
                        5 => w.s_dag(a),
                        6 => w.sqrt_x(a),
                        7 => w.sqrt_x_dag(a),
                        8 => w.sqrt_y(a),
                        9 => w.sqrt_y_dag(a),
                        10 => w.cnot(a, b),
                        11 => w.cz(a, b),
                        _ => w.cy(a, b),
                    }
                }
            }
            assert_eq!(plain.to_string(), lossy.to_string());
        }
    }

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
