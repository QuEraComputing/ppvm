use std::hash::BuildHasher;

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

pub trait ACMapAddAssign<S: PauliStorage, V: Coefficient, H: BuildHasher + Clone + Default> {
    fn add_assign(&mut self, key: PauliWord<S, H>, value: V);
    fn map_add_assign<F>(&self, dest: &mut Self, f: F)
    where
        F: Fn(&PauliWord<S, H>, &V) -> (PauliWord<S, H>, V) + Sync + Send;
}

pub trait ACMapMulAssign<V: Coefficient, H: BuildHasher + Clone + Default> {
    fn mul_assign(&mut self, value: V);
}

pub trait ACMapInsert<S: PauliStorage, V: Coefficient, H: BuildHasher + Clone + Default> {
    /// modify in place and insert some new entry into dest based on
    /// existing entries in self.
    fn map_insert<F>(&mut self, dest: &mut Self, f: F)
    where
        F: Fn(&PauliWord<S, H>, &mut V) -> Option<(PauliWord<S, H>, V)> + Sync + Send;
}

pub trait ACMapContains<S: PauliStorage, V: Coefficient, H: BuildHasher + Clone + Default> {
    fn contains(&self, key: &PauliWord<S, H>, value: &V) -> bool {
        self.contains_with(key, |v| v == value)
    }
    fn contains_with<F>(&self, key: &PauliWord<S, H>, f: F) -> bool
    where
        F: Fn(&V) -> bool;
}

pub trait ACMapConsume {
    /// consume dest into self, guranteeing accumulation of values with the same key.
    fn consume(&mut self, dest: &mut Self);
}

pub trait ACMapScale<S: PauliStorage, V: Coefficient, H: BuildHasher + Clone + Default> {
    fn scale<F>(&mut self, f: F)
    where
        F: Fn(&PauliWord<S, H>, &mut V) + Sync + Send;
}

pub trait ACMap<S: PauliStorage, V: Coefficient, H: BuildHasher + Clone + Default>:
    Clone
    + ACMapBase
    + ACMapAddAssign<S, V, H>
    + ACMapMulAssign<V, H>
    + ACMapInsert<S, V, H>
    + ACMapContains<S, V, H>
    + ACMapScale<S, V, H>
    + ACMapConsume
{
}

impl<T, Storage, Coeff, Hasher> ACMap<Storage, Coeff, Hasher> for T
where
    Storage: PauliStorage,
    Coeff: Coefficient,
    Hasher: BuildHasher + Clone + Default,
    T: Clone
        + ACMapBase
        + ACMapAddAssign<Storage, Coeff, Hasher>
        + ACMapMulAssign<Coeff, Hasher>
        + ACMapInsert<Storage, Coeff, Hasher>
        + ACMapScale<Storage, Coeff, Hasher>
        + ACMapContains<Storage, Coeff, Hasher>
        + ACMapConsume,
{
}
