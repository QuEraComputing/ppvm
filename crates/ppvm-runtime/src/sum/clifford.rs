use crate::{
    config::Config,
    phase::PhasedPauliWord,
    sum::PauliSum,
    traits::{ACMap, ACMapAddAssign, Clifford},
};

macro_rules! map_word {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            let (data, aux) = self.data_aux_mut();
            aux.clear();

            data.map_add_assign(aux, |k, v| {
                let mut p: PhasedPauliWord<T::Storage> = (*k).into();
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

impl<T: Config> Clifford for PauliSum<T>
where
    T::Coeff: std::ops::AddAssign + Clone + std::ops::Neg<Output = T::Coeff>,
    T::Map: ACMap<T::Storage, T::Coeff> + ACMapAddAssign<T::Storage, T::Coeff>,
{
    map_word!(x, index);
    map_word!(y, index);
    map_word!(z, index);
    map_word!(h, index);
    map_word!(s, index);
    map_word!(cnot, a, b);
    map_word!(cz, a, b);
}
