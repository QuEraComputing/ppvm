// use crate::{
//     config::Config, phase::PhasedPauliWord, sum::PauliSum, traits::ACMapIter, word::PauliWord,
// };

// impl<T: Config> std::ops::MulAssign<PauliWord<T::Storage, T::BuildHasher>> for PauliSum<T>
// where
//     T::Map: for<'a> ACMapIter<'a, Item = (&'a PauliWord<T::Storage, T::BuildHasher>, &'a T::Coeff)>,
//     T::BuildHasher: Sync + Send,
// {
//     fn mul_assign(&mut self, rhs: PauliWord<T::Storage, T::BuildHasher>) {
//         let phased_rhs = PhasedPauliWord {
//             word: rhs,
//             phase: 0,
//         };
//         self.map_add(|word, coeff| {
//             let phased_word = PhasedPauliWord {
//                 word: word.clone(),
//                 phase: 0,
//             };
//             let new_phased_word = phased_word * phased_rhs.clone();
//             let new_coeff = if new_phased_word.is_positive() {
//                 coeff.clone()
//             } else {
//                 -coeff.clone()
//             };
//             (new_phased_word.word, new_coeff)
//         });
//     }
// }

// impl<T: Config> std::ops::Mul<PauliWord<T::Storage, T::BuildHasher>> for PauliSum<T>
// where
//     T::Map: for<'a> ACMapIter<'a, Item = (&'a PauliWord<T::Storage, T::BuildHasher>, &'a T::Coeff)>,
//     T::BuildHasher: Sync + Send,
// {
//     type Output = PauliSum<T>;

//     fn mul(self, rhs: PauliWord<T::Storage, T::BuildHasher>) -> Self::Output {
//         let phased_rhs = PhasedPauliWord {
//             word: rhs,
//             phase: 0,
//         };
//         let mut output = self.clone();
//         output.map_insert(|word, coeff| {
//             let phased_word = PhasedPauliWord {
//                 word: word.clone(),
//                 phase: 0,
//             };
//             let new_phased_word = phased_word * phased_rhs.clone();
//             Some((new_phased_word.word, coeff.clone()))
//         });
//         output
//     }
// }
