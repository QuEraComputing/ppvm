use crate::config::Config;
use crate::traits::{self, ACMapBase, ACMapCombineUnique, ACMapInsert, ACMapIter};
use crate::word::PauliWord;

#[derive(Clone, Debug)]
pub struct PauliSum<T: Config> {
    map: (T::Map, T::Map),
    aux: bool,
    n_qubits: usize,
    capacity: usize,
}

impl<T: Config> PauliSum<T> {
    pub fn with_capacity(n_qubits: usize, capacity: usize) -> Self {
        Self {
            map: (
                T::Map::with_capacity(capacity),
                T::Map::with_capacity(capacity),
            ),
            aux: false,
            n_qubits,
            capacity,
        }
    }

    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    pub fn capacity(&self) -> usize {
        self.capacity
    }

    #[inline(always)]
    pub fn data(&self) -> &T::Map {
        if self.aux { &self.map.1 } else { &self.map.0 }
    }

    #[inline(always)]
    pub fn data_mut(&mut self) -> &mut T::Map {
        if self.aux {
            &mut self.map.1
        } else {
            &mut self.map.0
        }
    }

    #[inline(always)]
    pub fn aux(&self) -> &T::Map {
        if self.aux { &self.map.0 } else { &self.map.1 }
    }

    #[inline(always)]
    pub fn aux_mut(&mut self) -> &mut T::Map {
        if self.aux {
            &mut self.map.0
        } else {
            &mut self.map.1
        }
    }

    #[inline(always)]
    pub fn data_aux(&self) -> (&T::Map, &T::Map) {
        if self.aux {
            (&self.map.1, &self.map.0)
        } else {
            (&self.map.0, &self.map.1)
        }
    }

    #[inline(always)]
    pub fn data_aux_mut(&mut self) -> (&mut T::Map, &mut T::Map) {
        if self.aux {
            (&mut self.map.1, &mut self.map.0)
        } else {
            (&mut self.map.0, &mut self.map.1)
        }
    }

    #[inline(always)]
    pub fn swap(&mut self) {
        self.aux = !self.aux;
    }

    pub fn len(&self) -> usize {
        self.data().len()
    }
}

impl<'a, T: Config> PauliSum<T>
where
    T::Map: ACMapIter<'a>,
{
    pub fn iter(
        &'a self,
    ) -> <<T as Config>::Map as traits::ACMapIter<'a>>::Iter {
        self.data().iter()
    }
}

impl<T: Config> PauliSum<T>
where
    T::Map: ACMapCombineUnique,
{
    /// combine entries with the same key assuming unique keys
    /// in either data or aux. The combined entries are stored in `.data()`.
    pub fn combine_unique(&mut self) {
        let (data, aux) = self.data_aux_mut();
        if aux.len() > data.len() {
            aux.combine_unique(data);
            self.swap();
        } else {
            data.combine_unique(aux);
        }
    }
}

impl<T: Config> PauliSum<T>
where
    T::Map: ACMapInsert<T::Storage, T::Coeff> + ACMapCombineUnique,
{
    /// modify in place existing entries and insert some new entries
    /// if `f` return Some((k,v)) for an existing entry (k0,v0), then
    /// the existing entry is modified by `f` and a new entry (k,v) is added.
    /// if `f` return None, then the existing entry is only modified.
    /// finally, all entries are combined assuming unique keys.
    pub fn map_insert<F>(&mut self, f: F)
    where
        F: Fn(&PauliWord<T::Storage>, &mut T::Coeff) -> Option<(PauliWord<T::Storage>, T::Coeff)>
            + Sync
            + Send,
    {
        let (data, aux) = self.data_aux_mut();
        aux.clear();
        data.map_insert(aux, f);
        self.combine_unique();
    }
}

impl<T: Config> IntoIterator for PauliSum<T>
where
    T::Map: IntoIterator,
{
    type Item = <T::Map as IntoIterator>::Item;
    type IntoIter = <T::Map as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        if self.aux {
            self.map.1.into_iter()
        } else {
            self.map.0.into_iter()
        }
    }
}

impl<T: Config> PartialEq for PauliSum<T>
where
    T::Map: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.n_qubits == other.n_qubits && self.data() == other.data()
    }
}

impl<T: Config> Extend<(PauliWord<T::Storage>, T::Coeff)> for PauliSum<T>
where
    T::Map: Extend<(PauliWord<T::Storage>, T::Coeff)>,
{
    fn extend<I: IntoIterator<Item = (PauliWord<T::Storage>, T::Coeff)>>(&mut self, iter: I) {
        self.data_mut().extend(iter);
    }
}
