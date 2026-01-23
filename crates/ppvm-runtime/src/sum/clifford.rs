use crate::{config::Config, phase::PhasedPauliWord, sum::PauliSum, traits::Clifford};

macro_rules! map_scale {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            self.scale(|k, v| {
                let mut p: PhasedPauliWord<T::Storage, T::BuildHasher> = k.clone().into();
                let l = p.word.lbits;
                if $( l[$index] )||* {
                    *v *= 0.0;
                    return;
                }
                p.$name($($index),*);
                if !p.is_positive() {
                    *v *= -1.0;
                }
            })
        }
    };
}

macro_rules! map_word {
    ($name:ident, $index:ident) => {
        fn $name(&mut self, $index: usize) {
            self.map_add(|k, v| {
                let mut p: PhasedPauliWord<T::Storage, T::BuildHasher> = k.clone().into();

                let l = p.word.lbits;
                if l[$index] {
                    return (p.word, v.clone() * 0.0);
                }

                p.$name($index);
                if p.is_positive() {
                    (p.word, v.clone())
                } else {
                    (p.word, -v.clone())
                }
            })
        }
    };
    ($name:ident, $first:ident, $second:ident) => {
        fn $name(&mut self, $first: usize, $second: usize) {
            self.map_add(|k, v| {
                let mut p: PhasedPauliWord<T::Storage, T::BuildHasher> = k.clone().into();

                let l = p.word.lbits;
                if l[$second] {
                    return (p.word, v.clone() * 0.0);
                }

                if l[$first] {
                    return (p.word, v.clone());
                }

                p.$name($first, $second);
                if p.is_positive() {
                    (p.word, v.clone())
                } else {
                    (p.word, -v.clone())
                }
            })
        }
    };
}

impl<T: Config> Clifford for PauliSum<T> {
    map_scale!(x, index);
    map_scale!(y, index);
    map_scale!(z, index);
    map_word!(h, index);
    map_word!(s, index);
    map_word!(cnot, a, b);
    map_word!(cz, a, b);
}
