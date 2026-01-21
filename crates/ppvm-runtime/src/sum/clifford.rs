use crate::{config::Config, phase::PhasedPauliWord, sum::PauliSum, traits::Clifford};

macro_rules! map_scale {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            self.scale(|k, v| {
                let mut p: PhasedPauliWord<T::Storage, T::BuildHasher> = k.clone().into();
                p.$name($($index),*);
                if !p.is_positive() {
                    *v *= -1.0;
                }
            })
        }
    };
}

macro_rules! map_word {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            self.map_add(|k, v| {
                let mut p: PhasedPauliWord<T::Storage, T::BuildHasher> = k.clone().into();
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

impl<T: Config> Clifford for PauliSum<T> {
    map_scale!(x, index);
    map_scale!(y, index);
    map_scale!(z, index);
    map_word!(h, index);
    map_word!(s, index);
    map_word!(s_dagger, index);
    map_word!(cnot, a, b);
    map_word!(cz, a, b);
}
