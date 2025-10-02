use super::contains::Contains;
use super::data::PauliPattern;

use crate::traits::{PauliStorage, Trace};
use crate::word::PauliWord;
use std::hash::BuildHasher;

impl<'a, A, H> Trace<'a, PauliPattern> for PauliWord<A, H>
where
    A: PauliStorage + 'a,
    H: Default + BuildHasher + Clone + 'a,
{
    type Output = bool;
    fn trace(&'a self, value: &'a PauliPattern) -> Self::Output {
        value.contains(&self)
    }
}

impl<'a, A, H> Trace<'a, PauliWord<A, H>> for PauliPattern
where
    A: PauliStorage + 'a,
    H: Default + BuildHasher + Clone + 'a,
{
    type Output = bool;
    fn trace(&'a self, value: &'a PauliWord<A, H>) -> Self::Output {
        self.contains(value)
    }
}
