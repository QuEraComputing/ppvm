// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::sum::PauliSum;
use ppvm_traits::char::Pauli;
use ppvm_traits::config::Config;
use ppvm_traits::traits::*;

/// Reset a live qubit to `|0⟩` in the Heisenberg picture.
///
/// The adjoint of the reset channel maps `A_q ⊗ B` to `⟨0|A|0⟩ · I_q ⊗ B`:
/// `I` and `Z` on the reset qubit become `I`, while `X` and `Y` drop out.
/// Terms that collide after relabelling have their coefficients added.
/// Lost qubits are left untouched, matching the tableau, which skips
/// inactive qubits on state reset.
impl<T: Config> ResetState for PauliSum<T> {
    fn reset_state(&mut self, addr0: usize) {
        self.data_mut()
            .retain(|k, _| k.get_lbit(addr0) || !k.get_xbit(addr0));
        self.map_add(|k, v| {
            if !k.get_lbit(addr0) && k.get_zbit(addr0) {
                (k.set_new(addr0, Pauli::I), v.clone())
            } else {
                (k.clone(), v.clone())
            }
        });
    }
}

/// Full reset to `|0⟩` in the Heisenberg picture, term by term via
/// [`WordReset`]. For lossy words this also reloads lost qubits.
impl<T: Config> Reset for PauliSum<T>
where
    T::PauliWordType: WordReset,
{
    fn reset(&mut self, addr0: usize) {
        let (data, aux) = self.data_aux_mut();
        aux.clear();
        data.map_flat_add_assign(aux, |k, v| {
            let v = v.clone();
            k.word_reset(addr0).map(move |w| (w, v.clone()))
        });
        self.swap();
    }

    // The trait defaults run reset first (Schrödinger order); the Heisenberg
    // picture needs the basis change applied before the reset instead.
    fn reset_x(&mut self, addr0: usize) {
        self.h(addr0);
        self.reset(addr0);
    }

    fn reset_y(&mut self, addr0: usize) {
        self.s(addr0);
        self.h(addr0);
        self.reset(addr0);
    }
}

#[cfg(test)]
mod tests {
    use crate::config::fxhash::ByteF64;

    use super::*;
    use ppvm_pauli_word::loss::LossyPauliWord;

    type Sum = PauliSum<ByteF64<2>>;

    fn sum(n_qubits: usize, terms: &[(&str, f64)]) -> Sum {
        let mut s: Sum = PauliSum::builder().n_qubits(n_qubits).build();
        for &(word, coef) in terms {
            s += (word, coef);
        }
        s
    }

    #[test]
    fn test_reset_single_qubit_paulis() {
        let mut s = sum(1, &[("I", 0.5)]);
        s.reset(0);
        assert_eq!(s, sum(1, &[("I", 0.5)]));

        let mut s = sum(1, &[("Z", 0.5)]);
        s.reset(0);
        assert_eq!(s, sum(1, &[("I", 0.5)]));

        let mut s = sum(1, &[("X", 0.5), ("Y", 0.25)]);
        s.reset(0);
        assert!(s.is_empty());
    }

    #[test]
    fn test_reset_merges_colliding_terms() {
        // Qubit 0 is the leftmost character.
        let mut s = sum(2, &[("ZX", 0.5), ("IX", 0.25), ("YZ", 1.0)]);
        s.reset(0);
        assert_eq!(s, sum(2, &[("IX", 0.75)]));
    }

    #[test]
    fn test_reset_only_touches_target_qubit() {
        let mut s = sum(2, &[("XZ", 1.0), ("ZY", 1.0)]);
        s.reset(1);
        assert_eq!(s, sum(2, &[("XI", 1.0)]));
    }

    type Lossy = PauliSum<
        crate::config::fxhash::Byte<
            1,
            f64,
            ppvm_traits::traits::NoStrategy,
            LossyPauliWord<[u8; 1], fxhash::FxBuildHasher>,
        >,
    >;

    fn lossy(n_qubits: usize, terms: &[(&str, f64)]) -> Lossy {
        let mut s: Lossy = PauliSum::builder().n_qubits(n_qubits).build();
        for &(word, coef) in terms {
            s += (word, coef);
        }
        s
    }

    #[test]
    fn test_reset_gives_zero_state_even_from_loss() {
        // ⟨Z⟩ = 1 after reset on any state: Z maps to the full identity I + L.
        let mut o = lossy(1, &[("Z", 1.0)]);
        o.reset(0);
        assert_eq!(o, lossy(1, &[("I", 1.0), ("L", 1.0)]));
    }

    #[test]
    fn test_reset_leaves_no_loss_probability() {
        // The lost-projector expectation vanishes after a full reset.
        let mut o = lossy(2, &[("LZ", 1.0)]);
        o.reset(0);
        assert!(o.data().iter().all(|(_, v)| *v == 0.0));
    }

    #[test]
    fn test_reset_state_leaves_lost_qubit_untouched() {
        let mut s: Lossy = PauliSum::builder().n_qubits(2).build();
        s += ("LZ", 1.0);
        s += ("ZZ", 0.5);
        s.reset_state(0);
        let mut expect: Lossy = PauliSum::builder().n_qubits(2).build();
        expect += ("LZ", 1.0);
        expect += ("IZ", 0.5);
        assert_eq!(s, expect);
    }

    #[test]
    fn test_reset_after_x_gives_zero_state() {
        // X then reset: ⟨Z⟩ = +1. Heisenberg order: reset, then X.
        let mut o = sum(1, &[("Z", 1.0)]);
        o.reset(0);
        o.x(0);
        assert_eq!(o, sum(1, &[("I", 1.0)]));
    }

    #[test]
    fn test_reset_x_prepares_plus_state() {
        // reset_x then measure X: ⟨X⟩ = +1, ⟨Z⟩ = 0.
        let mut o = sum(1, &[("X", 1.0)]);
        o.reset_x(0);
        assert_eq!(o, sum(1, &[("I", 1.0)]));

        let mut o = sum(1, &[("Z", 1.0)]);
        o.reset_x(0);
        assert!(o.is_empty());
    }

    #[test]
    fn test_reset_y_prepares_plus_i_state() {
        // reset_y then measure Y: ⟨Y⟩ = +1 on |+i⟩.
        let mut o = sum(1, &[("Y", 1.0)]);
        o.reset_y(0);
        assert_eq!(o, sum(1, &[("I", 1.0)]));

        let mut o = sum(1, &[("X", 1.0)]);
        o.reset_y(0);
        assert!(o.is_empty());
    }
}
