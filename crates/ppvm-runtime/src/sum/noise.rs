use crate::char::Pauli;
use crate::traits::*;
use crate::{config::Config, sum::PauliSum};

impl<T: Config> PauliError<T> for PauliSum<T>
where
    f64: std::ops::Mul<T::Coeff, Output = T::Coeff>
        + std::ops::Add<T::Coeff, Output = T::Coeff>
        + std::ops::Sub<T::Coeff, Output = T::Coeff>,
{
    fn pauli_error(&mut self, addr0: usize, p: [<T as Config>::Coeff; 3]) {
        self.scale(|k, v| {
            match k.get(addr0) {
                Pauli::I => {}
                Pauli::X => {
                    *v *= 1.0f64 - 2.0f64 * p[1].clone() - 2.0f64 * p[2].clone();
                }
                Pauli::Y => {
                    *v *= 1.0f64 - 2.0f64 * p[0].clone() - 2.0f64 * p[2].clone();
                }
                Pauli::Z => {
                    *v *= 1.0f64 - 2.0f64 * p[0].clone() - 2.0f64 * p[1].clone();
                }
            };
        });
    }
}


impl<T: Config> LossChannel<T> for PauliSum<T>
where
    f64: std::ops::Mul<T::Coeff, Output = T::Coeff>
        + std::ops::Add<T::Coeff, Output = T::Coeff>
        + std::ops::Sub<T::Coeff, Output = T::Coeff>
        + std::ops::MulAssign<<T as Config>::Coeff>,
{
    fn loss_channel(&mut self, _addr0: usize, p: T::Coeff) {
        self.scale(|_, v| {
            *v *= 1.0f64 - p.clone();
        });
    }
}
