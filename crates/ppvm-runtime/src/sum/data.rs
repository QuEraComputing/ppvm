use crate::config::Config;
use crate::traits::{self, ACMap, ACMapIter};
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
    T::Map: ACMapIter<'a, T::Storage, T::Coeff>,
{
    pub fn iter(
        &'a self,
    ) -> <<T as Config>::Map as traits::ACMapIter<'a, T::Storage, T::Coeff>>::Iter {
        self.data().iter()
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
