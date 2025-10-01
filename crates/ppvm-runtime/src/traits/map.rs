use crate::{
    traits::{Coefficient, PauliStorage},
    word::PauliWord,
};

pub trait ACMap<S: PauliStorage, V: Coefficient>
{
    fn with_capacity(capacity: usize) -> Self;
    fn len(&self) -> usize;
}

pub trait ACMapIter<'a, S: PauliStorage, V: Coefficient> {
    type Item;
    type Iter: Iterator<Item = Self::Item>;
    fn iter(&'a self) -> Self::Iter;
}

pub trait ACMapIterMut<'a, S: PauliStorage, V: Coefficient> {
    type Item;
    type IterMut: Iterator<Item = Self::Item>;
    fn iter_mut(&'a mut self) -> Self::IterMut;
}

pub trait ACMapAddAssign<S: PauliStorage, V: Coefficient> {
    fn add_assign(&mut self, key: PauliWord<S>, value: V);
}

pub trait ACMapMulAssign<V: Coefficient> {
    fn mul_assign(&mut self, value: V);
}
