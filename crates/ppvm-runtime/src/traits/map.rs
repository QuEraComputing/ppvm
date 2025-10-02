use crate::{
    traits::{Coefficient, PauliStorage},
    word::PauliWord,
};

pub trait ACMapBase {
    fn with_capacity(capacity: usize) -> Self;
    fn len(&self) -> usize;
    fn clear(&mut self);
}

pub trait ACMapIter<'a> {
    type Item;
    type Iter: Iterator<Item = Self::Item>;
    fn iter(&'a self) -> Self::Iter;
}

pub trait ACMapAddAssign<S: PauliStorage, V: Coefficient> {
    fn add_assign(&mut self, key: PauliWord<S>, value: V);
    fn map_add_assign<F>(&self, dest: &mut Self, f: F)
    where
        F: Fn(&PauliWord<S>, &V) -> (PauliWord<S>, V) + Sync + Send;
}

pub trait ACMapMulAssign<V: Coefficient> {
    fn mul_assign(&mut self, value: V);
}

pub trait ACMapInsert<S: PauliStorage, V: Coefficient> {
    /// modify in place and insert some new entry into dest based on
    /// existing entries in self.
    fn map_insert<F>(&mut self, dest: &mut Self, f: F)
    where
        F: Fn(&PauliWord<S>, &mut V) -> Option<(PauliWord<S>, V)> + Sync + Send;
}

pub trait ACMapCombineUnique {
    /// combine two maps, assuming they have unique keys
    fn combine_unique(&mut self, dest: &mut Self);
}

pub trait ACMap<S: PauliStorage, V: Coefficient>:
    ACMapBase + ACMapAddAssign<S, V> + ACMapMulAssign<V> + ACMapInsert<S, V>
{
}

impl<T, Storage, Coeff> ACMap<Storage, Coeff> for T
where
    Storage: PauliStorage,
    Coeff: Coefficient,
    T: ACMapBase
        + ACMapAddAssign<Storage, Coeff>
        + ACMapMulAssign<Coeff>
        + ACMapInsert<Storage, Coeff>,
{
}
