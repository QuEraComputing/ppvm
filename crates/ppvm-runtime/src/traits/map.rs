use std::hash::BuildHasher;

use crate::traits::{Coefficient, PauliStorage, PauliWordTrait};

pub trait ACMapBase {
    fn with_capacity(capacity: usize) -> Self;
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    fn clear(&mut self);
}

pub trait ACMapIter<'a> {
    type Item;
    type Iter: Iterator<Item = Self::Item>;
    fn iter(&'a self) -> Self::Iter;
}

pub trait ACMapAddAssign<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>
{
    fn add_assign(&mut self, key: W, value: V);
    fn map_add_assign<F>(&self, dest: &mut Self, f: F)
    where
        F: Fn(&W, &V) -> (W, V) + Sync + Send;
}

pub trait ACMapMulAssign<V: Coefficient, H: BuildHasher + Clone + Default> {
    fn mul_assign(&mut self, value: V);
}

pub trait ACMapInsert<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>
{
    /// modify in place and insert some new entry into dest based on
    /// existing entries in self.
    fn map_insert<F>(&mut self, dest: &mut Self, f: F)
    where
        F: Fn(&W, &mut V) -> Option<(W, V)> + Sync + Send;

    /// modify in place and insert multiple new entries into dest based on
    /// existing entries in self.
    fn map_insert_multiple<F>(&mut self, dest: &mut Self, f: F)
    where
        F: Fn(&W, &mut V) -> Option<Vec<(W, V)>> + Sync + Send;
}

pub trait ACMapContains<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>
{
    fn contains(&self, key: &W, value: &V) -> bool {
        self.contains_with(key, |v| v == value)
    }
    fn contains_with<F>(&self, key: &W, f: F) -> bool
    where
        F: Fn(&V) -> bool;
}

pub trait ACMapConsume {
    /// consume dest into self, guaranteeing accumulation of values with the same key.
    fn consume(&mut self, dest: &mut Self);
}

pub trait ACMapScale<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>
{
    fn scale<F>(&mut self, f: F)
    where
        F: Fn(&W, &mut V) + Sync + Send;
}

pub trait ACMapRetain<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>
{
    fn retain<F>(&mut self, f: F)
    where
        F: Fn(&W, &V) -> bool + Sync + Send;
}

pub trait ACMap<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>:
    Clone
    + ACMapBase
    + ACMapAddAssign<S, V, H, W>
    + ACMapMulAssign<V, H>
    + ACMapInsert<S, V, H, W>
    + ACMapContains<S, V, H, W>
    + ACMapScale<S, V, H, W>
    + ACMapRetain<S, V, H, W>
    + ACMapConsume
{
}

impl<T, Storage, Coeff, Hasher, Word> ACMap<Storage, Coeff, Hasher, Word> for T
where
    Storage: PauliStorage,
    Coeff: Coefficient,
    Hasher: BuildHasher + Clone + Default,
    Word: PauliWordTrait,
    T: Clone
        + ACMapBase
        + ACMapAddAssign<Storage, Coeff, Hasher, Word>
        + ACMapMulAssign<Coeff, Hasher>
        + ACMapInsert<Storage, Coeff, Hasher, Word>
        + ACMapScale<Storage, Coeff, Hasher, Word>
        + ACMapContains<Storage, Coeff, Hasher, Word>
        + ACMapRetain<Storage, Coeff, Hasher, Word>
        + ACMapConsume,
{
}
