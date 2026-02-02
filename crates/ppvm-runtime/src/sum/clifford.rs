use crate::{
    config::Config, phase::PhasedPauliWord, sum::PauliSum, traits::Clifford, traits::PauliStorage,
    word::PauliWord,
};
use std::hash::BuildHasher;

macro_rules! map_scale {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            self.scale(|k, v| {
                let mut p: PhasedPauliWord<T::Storage, T::BuildHasher, <T as Config>::PauliWordType> = k.clone().into();
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

// NOTE: impl for PauliWord only; not a blanket, since PhasedPauliWord Clifford also isn't
impl<S, H, T> Clifford for PauliSum<T>
where
    S: PauliStorage,
    H: BuildHasher + Clone + Default,
    T: Config<Storage = S, BuildHasher = H, PauliWordType = PauliWord<S, H>>,
{
    map_scale!(x, index);
    map_scale!(y, index);
    map_scale!(z, index);
    map_word!(h, index);
    map_word!(s, index);
    map_word!(cnot, a, b);
    map_word!(cz, a, b);
}
