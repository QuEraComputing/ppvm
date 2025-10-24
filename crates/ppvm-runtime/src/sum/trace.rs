use crate::{config::Config, sum::PauliSum, traits::Trace, word::PauliWord};
use num::Zero;

impl<'a, T: Config, Rhs> Trace<'a, Rhs, T::Coeff> for PauliSum<T>
where
    <T as Config>::Coeff: Zero + Clone + std::ops::AddAssign + 'a,
    <T as Config>::Storage: 'a,
    <T as Config>::Map: Trace<'a, Rhs, T::Coeff>,
    <T as Config>::BuildHasher: 'a,
    Rhs: Trace<'a, PauliWord<T::Storage, T::BuildHasher>, T::Coeff> + 'a,
{
    fn trace(&'a self, value: &'a Rhs) -> T::Coeff {
        self.data().trace(value)
    }
}
