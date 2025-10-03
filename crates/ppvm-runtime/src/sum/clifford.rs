use crate::{
    config::Config,
    phase::PhasedPauliWord,
    sum::PauliSum,
    traits::{ACMapAddAssign, ACMapBase, Clifford},
};

macro_rules! map_word {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            let (data, aux) = self.data_aux_mut();
            aux.clear();

            data.map_add_assign(aux, |k, v| {
                let mut p: PhasedPauliWord<T::Storage, T::BuildHasher> = k.clone().into();
                p.$name($($index),*);
                if p.is_positive() {
                    (p.word, v.clone())
                } else {
                    (p.word, -v.clone())
                }
            });
            self.swap();
        }
    };
}

impl<T: Config> Clifford for PauliSum<T> {
    map_word!(x, index);
    map_word!(y, index);
    map_word!(z, index);
    map_word!(h, index);
    map_word!(s, index);
    map_word!(cnot, a, b);
    map_word!(cz, a, b);
}
