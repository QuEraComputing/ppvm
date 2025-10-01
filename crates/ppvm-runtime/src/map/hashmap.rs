use std::collections::HashMap;

use crate::{
    traits::{ACMapIterMut, Coefficient, PauliStorage},
    word::PauliWord,
};

// impl<'a, S, V> ACMapIterMut<'a, S, V> for HashMap<PauliWord<S>, V>
// where
//     S: PauliStorage + 'a,
//     V: Coefficient + 'a,
// {
//     type IterMut = std::collections::hash_map::IterMut<'a, PauliWord<S>, V>;
//     fn iter_mut(&'a mut self) -> Self::IterMut {
//         self.iter_mut()
//     }
// }
