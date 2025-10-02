use crate::traits::*;
use crate::{char::Pauli, config::Config, sum::PauliSum};

impl<T: Config> Projection for PauliSum<T>
where
    T::Coeff: std::ops::MulAssign + std::ops::Neg<Output = T::Coeff> + Clone,
    T::Map: ACMapInsert<T::Storage, T::Coeff> + ACMapCombineUnique,
{
    fn p0(&mut self, pos: usize) {
        self.map_insert(|k, v| {
            let half = v.half();
            match k.get(pos) {
                Pauli::I => {
                    *v *= half;
                    let nk = k.set_new(pos, Pauli::Z);
                    return Some((nk, v.clone()));
                }
                Pauli::Z => {
                    *v *= half;
                    let nk = k.set_new(pos, Pauli::I);
                    return Some((nk, v.clone()));
                }
                _ => {
                    return None;
                }
            }
        });
    }

    fn p1(&mut self, pos: usize) {
        self.map_insert(|k, v| {
            let half = v.half();
            match k.get(pos) {
                Pauli::I => {
                    *v *= half;
                    let nk = k.set_new(pos, Pauli::Z);
                    return Some((nk, -v.clone()));
                }
                Pauli::Z => {
                    *v *= half;
                    let nk = k.set_new(pos, Pauli::I);
                    return Some((nk, -v.clone()));
                }
                _ => {
                    return None;
                }
            }
        });
    }
}
