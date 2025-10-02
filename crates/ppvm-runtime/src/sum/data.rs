use crate::config::Config;
use crate::traits::*;
use crate::word::PauliWord;

#[derive(Clone)]
pub struct PauliSum<T: Config> {
    map: (T::Map, T::Map),
    aux: bool,
    n_qubits: usize,
    capacity: usize,
}

impl<T: Config> PauliSum<T> {
    pub fn new(n_qubits: usize) -> Self {
        Self::with_capacity(n_qubits, 1 << n_qubits)
    }

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

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn contains(&self, key: &PauliWord<T::Storage, T::BuildHasher>, value: &T::Coeff) -> bool {
        self.data().contains(key, value)
    }

    /// combine entries with the same key
    /// in either data or aux. The combined entries are stored in `.data()`.
    pub fn consume(&mut self) {
        let (data, aux) = self.data_aux_mut();
        if aux.len() > data.len() {
            aux.consume(data);
            self.swap();
        } else {
            data.consume(aux);
        }
    }

    /// modify in place existing entries and insert some new entries
    /// if `f` return Some((k,v)) for an existing entry (k0,v0), then
    /// the existing entry is modified by `f` and a new entry (k,v) is added.
    /// if `f` return None, then the existing entry is only modified.
    /// finally, all entries are combined assuming unique keys.
    pub fn map_insert<F>(&mut self, f: F)
    where
        F: Fn(
                &PauliWord<T::Storage, T::BuildHasher>,
                &mut T::Coeff,
            ) -> Option<(PauliWord<T::Storage, T::BuildHasher>, T::Coeff)>
            + Sync
            + Send,
    {
        let (data, aux) = self.data_aux_mut();
        aux.clear();
        data.map_insert(aux, f);
        self.consume();
    }
}

impl<'a, T: Config> PauliSum<T>
where
    T::Map: ACMapIter<'a>,
{
    pub fn iter(&'a self) -> <<T as Config>::Map as ACMapIter<'a>>::Iter {
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

impl<T: Config> Extend<(PauliWord<T::Storage, T::BuildHasher>, T::Coeff)> for PauliSum<T>
where
    T::Map: Extend<(PauliWord<T::Storage, T::BuildHasher>, T::Coeff)>,
{
    fn extend<I: IntoIterator<Item = (PauliWord<T::Storage, T::BuildHasher>, T::Coeff)>>(
        &mut self,
        iter: I,
    ) {
        self.data_mut().extend(iter);
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_yaml_snapshot;

    use super::*;
    use crate::config::fxhash::ByteF64;

    #[test]
    fn test_pauli_sum_creation() {
        let word = PauliWord::<[u8; 2]>::new(4);
        let mut sum: PauliSum<ByteF64<2>> = PauliSum::new(word.n_qubits());
        assert!(sum.data().is_empty());
        sum += "IIII";
        assert!(!sum.data().is_empty());
        assert_yaml_snapshot!(sum.to_string());
        sum += ("IIII", 2.0);
        assert_yaml_snapshot!(sum.to_string());
        sum += ("XIII", 1.0);
        assert_yaml_snapshot!(sum.to_string());
        sum += ("XIII", 2.0);
        assert_yaml_snapshot!(sum.to_string());
        sum += "IYII";
        assert_yaml_snapshot!(sum.to_string());
        assert!(sum.contains(&PauliWord::from("IIII"), &3.0));
        assert!(sum.contains(&PauliWord::from("XIII"), &3.0));
        assert!(sum.contains(&PauliWord::from("IYII"), &1.0));
    }

    #[test]
    fn test_pauli_sum_top_bottom() {
        let mut sum: PauliSum<ByteF64<2>> = PauliSum::new(4);
        assert!(sum.is_empty());
        sum += ("IIII", 1.0);
        assert!(!sum.is_empty());
        sum += ("IIII", 1.0);
    }
}
