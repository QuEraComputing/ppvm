use crate::{
    config::Config,
    phase::PhasedPauliWord,
    sum::PauliSum,
    traits::{Clifford, CliffordExtensions, PauliWordTrait},
};

macro_rules! map_word {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            self.map_add(|k, v| {
                let mut p: PhasedPauliWord<T::Storage, T::BuildHasher, <T as Config>::PauliWordType> = k.clone().into();
                p.$name($($index),*);
                if p.is_positive() {
                    (p.word, v.clone())
                } else {
                    (p.word, -v.clone())
                }
            })
        }
    };
}

// Single- and two-qubit Clifford action on a `PauliSum`. Single-qubit gates
// are specialised to bit-level updates so we avoid round-tripping through
// `PhasedPauliWord`; the two-qubit gates still go through the macro because
// their sign rules are easier to reuse from the `PhasedPauliWord` impl.
impl<T: Config> Clifford for PauliSum<T> {
    #[inline]
    fn x(&mut self, index: usize) {
        // X conjugation: only flips sign for Z and Y (zbit set); word unchanged
        self.scale(|k, v| {
            if !k.get_lbit(index) && k.get_zbit(index) {
                *v *= -1.0;
            }
        })
    }

    #[inline]
    fn y(&mut self, index: usize) {
        // Y conjugation: flips sign for X and Z (xbit XOR zbit); word unchanged
        self.scale(|k, v| {
            if !k.get_lbit(index) && (k.get_xbit(index) ^ k.get_zbit(index)) {
                *v *= -1.0;
            }
        })
    }

    #[inline]
    fn z(&mut self, index: usize) {
        // Z conjugation: only flips sign for X and Y (xbit set); word unchanged
        self.scale(|k, v| {
            if !k.get_lbit(index) && k.get_xbit(index) {
                *v *= -1.0;
            }
        })
    }

    #[inline]
    fn h(&mut self, index: usize) {
        // H swaps x and z bits; phase flip only for Y (both bits set)
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

    #[inline]
    fn s(&mut self, index: usize) {
        // S: zbit = xbit XOR zbit; phase flip only for X (xbit set, zbit clear)
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

    map_word!(cnot, a, b);
    map_word!(cz, a, b);
}

impl<T: Config> CliffordExtensions for PauliSum<T> {
    #[inline]
    fn s_adj(&mut self, addr0: usize) {
        // S†: same bit map as S; phase flip for Y (both bits set)
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

    #[inline]
    fn sqrt_x(&mut self, addr0: usize) {
        // √X: xbit = xbit XOR zbit; phase flip for Y (both bits set)
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

    #[inline]
    fn sqrt_y(&mut self, addr0: usize) {
        // √Y: swap x and z bits; phase flip for Z (zbit set, xbit clear)
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

    #[inline]
    fn sqrt_x_adj(&mut self, addr0: usize) {
        // √X†: xbit = xbit XOR zbit; phase flip for Z (zbit set, xbit clear)
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

    #[inline]
    fn sqrt_y_adj(&mut self, addr0: usize) {
        // √Y†: swap x and z bits; phase flip for X (xbit set, zbit clear)
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

    map_word!(cy, a, b);
}
