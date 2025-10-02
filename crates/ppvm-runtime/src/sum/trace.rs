use crate::{config::Config, sum::PauliSum, traits::Trace, word::PauliWord};
use num::Zero;

impl<'a, T: Config, Rhs> Trace<'a, Rhs> for PauliSum<T>
where
    <T as Config>::Coeff: Zero + Clone + std::ops::AddAssign + 'a,
    <T as Config>::Storage: 'a,
    <T as Config>::Map: Trace<'a, Rhs, Output = <T as Config>::Coeff>,
    <T as Config>::BuildHasher: 'a,
    Rhs: Trace<'a, PauliWord<T::Storage, T::BuildHasher>, Output = bool> + 'a,
{
    type Output = T::Coeff;
    fn trace(&'a self, value: &'a Rhs) -> Self::Output {
        self.data().trace(value)
    }
}
