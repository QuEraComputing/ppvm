use crate::traits::{Coefficient, PauliStorage};

pub trait ACMap<S: PauliStorage, V: Coefficient> {
    type Iter: Iterator<Item = (S, V)>;
    fn with_capacity(capacity: usize) -> Self;
    fn len(&self) -> usize;
    fn iter(&self) -> Self::Iter;
}
