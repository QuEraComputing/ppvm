use super::contains::Contains;
use super::data::PauliPattern;

use crate::traits::{PauliStorage, Trace};
use crate::word::PauliWord;
use std::hash::BuildHasher;

impl<'a, A, H, T> Trace<'a, PauliPattern, T> for PauliWord<A, H>
where
    A: PauliStorage + 'a,
    H: Default + BuildHasher + Clone + 'a,
    T: From<bool>,
{
    fn trace(&'a self, value: &'a PauliPattern) -> T {
        T::from(value.contains(&self))
    }
}

impl<'a, A, H, T> Trace<'a, PauliWord<A, H>, T> for PauliPattern
where
    A: PauliStorage + 'a,
    H: Default + BuildHasher + Clone + 'a,
    T: From<bool>,
{
    fn trace(&'a self, value: &'a PauliWord<A, H>) -> T {
        T::from(self.contains(value))
    }
}
