use std::collections::HashMap;
use std::hash::BuildHasher;

use crate::Tableau;

pub trait TableauMapBase {
    fn with_capacity(capacity: usize) -> Self;
    fn len(&self) -> usize;
    fn clear(&mut self);
}

pub trait TableauMapIter<'a, C: 'a> {
    type Iter: Iterator<Item = (&'a Tableau, &'a C)>;
    fn iter(&'a self) -> Self::Iter;
}

pub trait TableauMapDrain<C> {
    type Drain<'a>: Iterator<Item = (Tableau, C)>
    where
        Self: 'a;

    fn drain(&mut self) -> Self::Drain<'_>;
}

pub trait TableauMapAddAssign<C, H> {
    fn add_assign(&mut self, key: Tableau, value: C);
}

pub trait TableauMapGet<C> {
    fn get(&self, key: &Tableau) -> Option<&C>;
}

pub trait TableauMap<C, H>:
    Clone
    + TableauMapBase
    + for<'a> TableauMapIter<'a, C>
    + TableauMapDrain<C>
    + TableauMapAddAssign<C, H>
    + TableauMapGet<C>
{
}

impl<T, C, H> TableauMap<C, H> for T where
    T: Clone
        + TableauMapBase
        + for<'a> TableauMapIter<'a, C>
        + TableauMapDrain<C>
        + TableauMapAddAssign<C, H>
        + TableauMapGet<C>
{
}

impl<C, H> TableauMapBase for HashMap<Tableau, C, H>
where
    H: BuildHasher + Clone + Default,
{
    fn with_capacity(capacity: usize) -> Self {
        Self::with_capacity_and_hasher(capacity, H::default())
    }

    fn len(&self) -> usize {
        HashMap::len(self)
    }

    fn clear(&mut self) {
        HashMap::clear(self)
    }
}

impl<'a, C: 'a, H> TableauMapIter<'a, C> for HashMap<Tableau, C, H>
where
    H: BuildHasher + Clone + Default,
{
    type Iter = std::collections::hash_map::Iter<'a, Tableau, C>;

    fn iter(&'a self) -> Self::Iter {
        HashMap::iter(self)
    }
}

impl<C, H> TableauMapDrain<C> for HashMap<Tableau, C, H>
where
    H: BuildHasher + Clone + Default,
{
    type Drain<'a> = std::collections::hash_map::Drain<'a, Tableau, C> where Self: 'a;

    fn drain(&mut self) -> Self::Drain<'_> {
        HashMap::drain(self)
    }
}

impl<C, H> TableauMapAddAssign<C, H> for HashMap<Tableau, C, H>
where
    C: Clone + std::ops::AddAssign,
    H: BuildHasher + Clone + Default,
{
    fn add_assign(&mut self, key: Tableau, value: C) {
        self.entry(key)
            .and_modify(|v| *v += value.clone())
            .or_insert(value);
    }
}

impl<C, H> TableauMapGet<C> for HashMap<Tableau, C, H>
where
    H: BuildHasher + Clone + Default,
{
    fn get(&self, key: &Tableau) -> Option<&C> {
        HashMap::get(self, key)
    }
}

#[cfg(feature = "indexmap")]
impl<C, H> TableauMapBase for indexmap::IndexMap<Tableau, C, H>
where
    H: BuildHasher + Clone + Default,
{
    fn with_capacity(capacity: usize) -> Self {
        indexmap::IndexMap::with_capacity_and_hasher(capacity, H::default())
    }

    fn len(&self) -> usize {
        indexmap::IndexMap::len(self)
    }

    fn clear(&mut self) {
        indexmap::IndexMap::clear(self)
    }
}

#[cfg(feature = "indexmap")]
impl<'a, C: 'a, H> TableauMapIter<'a, C> for indexmap::IndexMap<Tableau, C, H>
where
    H: BuildHasher + Clone + Default,
{
    type Iter = indexmap::map::Iter<'a, Tableau, C>;

    fn iter(&'a self) -> Self::Iter {
        indexmap::IndexMap::iter(self)
    }
}

#[cfg(feature = "indexmap")]
impl<C, H> TableauMapDrain<C> for indexmap::IndexMap<Tableau, C, H>
where
    H: BuildHasher + Clone + Default,
{
    type Drain<'a> = indexmap::map::Drain<'a, Tableau, C> where Self: 'a;

    fn drain(&mut self) -> Self::Drain<'_> {
        self.drain(..)
    }
}

#[cfg(feature = "indexmap")]
impl<C, H> TableauMapAddAssign<C, H> for indexmap::IndexMap<Tableau, C, H>
where
    C: Clone + std::ops::AddAssign,
    H: BuildHasher + Clone + Default,
{
    fn add_assign(&mut self, key: Tableau, value: C) {
        self.entry(key)
            .and_modify(|v| *v += value.clone())
            .or_insert(value);
    }
}

#[cfg(feature = "indexmap")]
impl<C, H> TableauMapGet<C> for indexmap::IndexMap<Tableau, C, H>
where
    H: BuildHasher + Clone + Default,
{
    fn get(&self, key: &Tableau) -> Option<&C> {
        indexmap::IndexMap::get(self, key)
    }
}
