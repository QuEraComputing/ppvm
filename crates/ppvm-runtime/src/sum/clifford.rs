// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::{
    config::Config,
    phase::PhasedPauliWord,
    sum::PauliSum,
    traits::{Clifford, CliffordExtensions, PauliWordTrait, Targets},
};

// Two-qubit Clifford action on a `PauliSum`, broadcast over consecutive
// pairs. The per-pair sign rules are reused from the `PhasedPauliWord` impl.
macro_rules! map_word {
    ($name:ident) => {
        fn $name(&mut self, targets: impl Targets) {
            for (a, b) in targets.pairs() {
                self.map_add(|k, v| {
                    let mut p: PhasedPauliWord<T::Storage, T::BuildHasher, <T as Config>::PauliWordType> = k.clone().into();
                    p.$name([a, b]);
                    if p.is_positive() {
                        (p.word, v.clone())
                    } else {
                        (p.word, -v.clone())
                    }
                })
            }
        }
    };
}

// Single- and two-qubit Clifford action on a `PauliSum`. Single-qubit gates
// are specialised to bit-level updates so we avoid round-tripping through
// `PhasedPauliWord`; the two-qubit gates still go through the macro because
// their sign rules are easier to reuse from the `PhasedPauliWord` impl.
impl<T: Config> Clifford for PauliSum<T> {
    #[inline]
    fn x(&mut self, targets: impl Targets) {
        // X conjugation: only flips sign for Z and Y (zbit set); word unchanged
        for index in targets.each() {
            self.scale(|k, v| {
                if !k.get_lbit(index) && k.get_zbit(index) {
                    *v *= -1.0;
                }
            })
        }
    }

    #[inline]
    fn y(&mut self, targets: impl Targets) {
        // Y conjugation: flips sign for X and Z (xbit XOR zbit); word unchanged
        for index in targets.each() {
            self.scale(|k, v| {
                if !k.get_lbit(index) && (k.get_xbit(index) ^ k.get_zbit(index)) {
                    *v *= -1.0;
                }
            })
        }
    }

    #[inline]
    fn z(&mut self, targets: impl Targets) {
        // Z conjugation: only flips sign for X and Y (xbit set); word unchanged
        for index in targets.each() {
            self.scale(|k, v| {
                if !k.get_lbit(index) && k.get_xbit(index) {
                    *v *= -1.0;
                }
            })
        }
    }

    #[inline]
    fn h(&mut self, targets: impl Targets) {
        // H swaps x and z bits; phase flip only for Y (both bits set)
        for index in targets.each() {
            self.map_add(|k, v| {
                if k.get_lbit(index) {
                    return (k.clone(), v.clone());
                }
                let xbit = k.get_xbit(index);
                let zbit = k.get_zbit(index);
                let mut new_word = k.clone();
                new_word.set_xbit(index, zbit);
                new_word.set_zbit(index, xbit);
                new_word.rehash();
                if xbit & zbit {
                    (new_word, -v.clone())
                } else {
                    (new_word, v.clone())
                }
            })
        }
    }

    #[inline]
    fn s(&mut self, targets: impl Targets) {
        // S: zbit = xbit XOR zbit; phase flip only for X (xbit set, zbit clear)
        for index in targets.each() {
            self.map_add(|k, v| {
                if k.get_lbit(index) {
                    return (k.clone(), v.clone());
                }
                let xbit = k.get_xbit(index);
                let zbit = k.get_zbit(index);
                let mut new_word = k.clone();
                new_word.set_zbit(index, xbit ^ zbit);
                new_word.rehash();
                if xbit & !zbit {
                    (new_word, -v.clone())
                } else {
                    (new_word, v.clone())
                }
            })
        }
    }

    map_word!(cnot);
    map_word!(cz);
}

impl<T: Config> CliffordExtensions for PauliSum<T> {
    #[inline]
    fn s_dag(&mut self, targets: impl Targets) {
        // S†: same bit map as S; phase flip for Y (both bits set)
        for addr0 in targets.each() {
            self.map_add(|k, v| {
                if k.get_lbit(addr0) {
                    return (k.clone(), v.clone());
                }
                let xbit = k.get_xbit(addr0);
                let zbit = k.get_zbit(addr0);
                let mut new_word = k.clone();
                new_word.set_zbit(addr0, xbit ^ zbit);
                new_word.rehash();
                if xbit & zbit {
                    (new_word, -v.clone())
                } else {
                    (new_word, v.clone())
                }
            })
        }
    }

    #[inline]
    fn sqrt_x(&mut self, targets: impl Targets) {
        // √X: xbit = xbit XOR zbit; phase flip for Y (both bits set)
        for addr0 in targets.each() {
            self.map_add(|k, v| {
                if k.get_lbit(addr0) {
                    return (k.clone(), v.clone());
                }
                let xbit = k.get_xbit(addr0);
                let zbit = k.get_zbit(addr0);
                let mut new_word = k.clone();
                new_word.set_xbit(addr0, xbit ^ zbit);
                new_word.rehash();
                if xbit & zbit {
                    (new_word, -v.clone())
                } else {
                    (new_word, v.clone())
                }
            })
        }
    }

    #[inline]
    fn sqrt_y(&mut self, targets: impl Targets) {
        // √Y: swap x and z bits; phase flip for Z (zbit set, xbit clear)
        for addr0 in targets.each() {
            self.map_add(|k, v| {
                if k.get_lbit(addr0) {
                    return (k.clone(), v.clone());
                }
                let xbit = k.get_xbit(addr0);
                let zbit = k.get_zbit(addr0);
                let mut new_word = k.clone();
                new_word.set_xbit(addr0, zbit);
                new_word.set_zbit(addr0, xbit);
                new_word.rehash();
                if !xbit & zbit {
                    (new_word, -v.clone())
                } else {
                    (new_word, v.clone())
                }
            })
        }
    }

    #[inline]
    fn sqrt_x_dag(&mut self, targets: impl Targets) {
        // √X†: xbit = xbit XOR zbit; phase flip for Z (zbit set, xbit clear)
        for addr0 in targets.each() {
            self.map_add(|k, v| {
                if k.get_lbit(addr0) {
                    return (k.clone(), v.clone());
                }
                let xbit = k.get_xbit(addr0);
                let zbit = k.get_zbit(addr0);
                let mut new_word = k.clone();
                new_word.set_xbit(addr0, xbit ^ zbit);
                new_word.rehash();
                if !xbit & zbit {
                    (new_word, -v.clone())
                } else {
                    (new_word, v.clone())
                }
            })
        }
    }

    #[inline]
    fn sqrt_y_dag(&mut self, targets: impl Targets) {
        // √Y†: swap x and z bits; phase flip for X (xbit set, zbit clear)
        for addr0 in targets.each() {
            self.map_add(|k, v| {
                if k.get_lbit(addr0) {
                    return (k.clone(), v.clone());
                }
                let xbit = k.get_xbit(addr0);
                let zbit = k.get_zbit(addr0);
                let mut new_word = k.clone();
                new_word.set_xbit(addr0, zbit);
                new_word.set_zbit(addr0, xbit);
                new_word.rehash();
                if xbit & !zbit {
                    (new_word, -v.clone())
                } else {
                    (new_word, v.clone())
                }
            })
        }
    }

    map_word!(cy);
}
