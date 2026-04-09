use std::hash::BuildHasher;

use num::Integer;

use crate::char::Pauli;
use crate::traits::{PauliStorage, PauliWordTrait};
use crate::word::PauliWord;

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct PhasedPauliWord<
    A: PauliStorage,
    H = fxhash::FxBuildHasher,
    W: PauliWordTrait = PauliWord<A, H>,
> {
    pub word: W,
    /// 0: +1, 1: +i, 2: -1, 3: -i
    ///   sign imag
    /// +1: 0   0
    /// +i: 0   1
    /// -1: 1   0
    /// -i: 1   1
    pub phase: u8,
    _phantom: std::marker::PhantomData<(A, H)>,
}

impl<A: PauliStorage, H: BuildHasher, W: PauliWordTrait> PhasedPauliWord<A, H, W> {
    pub fn new(n_qubits: usize) -> Self {
        Self {
            word: W::new(n_qubits),
            phase: 0, // Default phase is +1
            _phantom: std::marker::PhantomData,
        }
    }

    pub fn build_from_word(word: W, phase: u8) -> Self {
        Self {
            word,
            phase,
            _phantom: std::marker::PhantomData,
        }
    }

    pub fn n_qubits(&self) -> usize {
        self.word.n_qubits()
    }

    pub fn is_positive(&self) -> bool {
        // is second bit 0
        (self.phase & 0b10) == 0
    }

    #[inline(always)]
    pub fn add_phase(&mut self, phase: u8) {
        self.phase = (self.phase + phase) % 4;
    }

    #[inline(always)]
    pub fn get(&self, index: usize) -> Pauli {
        self.word.get(index)
    }

    #[inline(always)]
    pub fn set(&mut self, index: usize, pauli: Pauli) -> &mut Self {
        self.word.set(index, pauli);
        self
    }

    #[inline(always)]
    pub fn set_new(&self, index: usize, pauli: Pauli) -> Self {
        let new_words = self.word.set_new(index, pauli);
        Self {
            word: new_words,
            phase: self.phase,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<A: PauliStorage, H: BuildHasher + Default + Clone, W: PauliWordTrait> From<W>
    for PhasedPauliWord<A, H, W>
{
    fn from(words: W) -> Self {
        Self {
            word: words,
            phase: 0,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<H: BuildHasher + Default + Clone, W: PauliWordTrait> From<String>
    for PhasedPauliWord<u64, H, W>
{
    fn from(s: String) -> Self {
        let mut chars = s.chars();
        let phase: u8 = match (chars.next(), chars.next()) {
            (Some('+'), Some('i')) => 1, // +i
            (Some('-'), Some('i')) => 3, // -i
            (Some('+'), _) => 0,         // +1
            (Some('-'), _) => 2,         // -1
            _ => panic!("Invalid phase format"),
        };
        // Remaining characters are the Pauli string
        let s: String = s.chars().skip(if phase.is_odd() { 2 } else { 1 }).collect();
        let words: W = s.into();
        Self {
            word: words,
            phase,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<H: BuildHasher + Default + Clone, W: PauliWordTrait> From<&str>
    for PhasedPauliWord<u64, H, W>
{
    fn from(s: &str) -> Self {
        s.to_string().into()
    }
}

impl<A: PauliStorage, H: BuildHasher + Default + Clone, const REHASH: bool> std::fmt::Display
    for PhasedPauliWord<A, H, PauliWord<A, H, REHASH>>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.phase {
            0 => write!(f, "+")?,
            1 => write!(f, "+i")?,
            2 => write!(f, "-")?,
            3 => write!(f, "-i")?,
            _ => unreachable!("Invalid phase value: {}", self.phase),
        };
        write!(f, "{}", self.word)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_pauli_string_with_phase_creation() {
        let ps = PhasedPauliWord::<u64>::new(4);
        assert_eq!(ps.word.n_qubits(), 4);
    }

    #[test]
    fn test_pauli_string_with_phase_set_get() {
        let ps = PhasedPauliWord::<u64>::new(4);
        let ps = ps.set_new(0, Pauli::X);
        assert_eq!(ps.get(0), Pauli::X);
        let ps = ps.set_new(1, Pauli::Y);
        assert_eq!(ps.get(1), Pauli::Y);
        let ps = ps.set_new(2, Pauli::Z);
        assert_eq!(ps.get(2), Pauli::Z);
        let ps = ps.set_new(3, Pauli::I);
        assert_eq!(ps.get(3), Pauli::I);
    }

    #[test]
    fn test_pauli_string_with_phase_display() {
        let ps = PhasedPauliWord::<u64>::new(4)
            .set_new(0, Pauli::X)
            .set_new(1, Pauli::Y)
            .set_new(2, Pauli::Z)
            .set_new(3, Pauli::I);
        assert_eq!(ps.to_string(), "+XYZI");
        let ps = PhasedPauliWord::<u64>::new(4)
            .set_new(0, Pauli::X)
            .set_new(1, Pauli::Y)
            .set_new(2, Pauli::Z)
            .set_new(3, Pauli::I);
        assert_eq!(ps.to_string(), "+XYZI");
        let ps = PhasedPauliWord::<u64>::new(4)
            .set_new(0, Pauli::X)
            .set_new(1, Pauli::Y)
            .set_new(2, Pauli::Z)
            .set_new(3, Pauli::I);
        assert_eq!(ps.to_string(), "+XYZI");
        let ps = PhasedPauliWord::<u64>::new(4)
            .set_new(0, Pauli::X)
            .set_new(1, Pauli::Y)
            .set_new(2, Pauli::Z)
            .set_new(3, Pauli::I);
        assert_eq!(ps.to_string(), "+XYZI");
    }
    #[test]
    fn test_pauli_string_with_phase_from_string() {
        let ps: PhasedPauliWord<u64> = "+XYZI".to_string().into();
        assert_eq!(ps.word.n_qubits(), 4);
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        assert_eq!(ps.phase, 0);
        let ps: PhasedPauliWord<u64> = "-XYZI".into();
        assert_eq!(ps.word.n_qubits(), 4);
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        assert_eq!(ps.phase, 2);
        let ps: PhasedPauliWord<u64> = "+iXYZI".to_string().into();
        assert_eq!(ps.word.n_qubits(), 4);
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        assert_eq!(ps.phase, 1);
        let ps: PhasedPauliWord<u64> = "-iXYZI".to_string().into();
        assert_eq!(ps.word.n_qubits(), 4);
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        assert_eq!(ps.phase, 3);
    }

    #[test]
    fn test_is_positive() {
        let ps: PhasedPauliWord<u64> = "+XYZI".into();
        assert!(ps.is_positive());
        let ps: PhasedPauliWord<u64> = "-XYZI".into();
        assert!(!ps.is_positive());
        let ps: PhasedPauliWord<u64> = "+iXYZI".into();
        assert!(ps.is_positive());
        let ps: PhasedPauliWord<u64> = "-iXYZI".into();
        assert!(!ps.is_positive());
    }
}
