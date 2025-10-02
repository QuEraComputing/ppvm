use std::hash::BuildHasher;

use num::Integer;

use crate::char::Pauli;
use crate::traits::PauliStorage;
use crate::word::PauliWord;

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct PhasedPauliWord<A: PauliStorage, H = fxhash::FxBuildHasher> {
    pub word: PauliWord<A, H>,
    pub phase: u8, // 0: +1, 1: -1, 2: +i, 3: -i
}

impl<A: PauliStorage, H: BuildHasher + Default + Clone> PhasedPauliWord<A, H> {
    pub fn new(n_qubits: usize) -> Self {
        Self {
            word: PauliWord::new(n_qubits),
            phase: 0, // Default phase is +1
        }
    }

    pub fn is_positive(&self) -> bool {
        self.phase.is_even()
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
        }
    }
}

impl<A: PauliStorage, H: BuildHasher + Default + Clone> From<PauliWord<A, H>> for PhasedPauliWord<A, H> {
    fn from(words: PauliWord<A, H>) -> Self {
        Self {
            word: words,
            phase: 0,
        }
    }
}

impl<S: AsRef<str>, H: BuildHasher + Default + Clone> From<S> for PhasedPauliWord<u64, H> {
    fn from(s: S) -> Self {
        let mut chars = s.as_ref().chars();
        let phase: u8 = match (chars.next(), chars.next()) {
            (Some('+'), Some('i')) => 2, // +i
            (Some('-'), Some('i')) => 3, // -i
            (Some('+'), _) => 0,         // +1
            (Some('-'), _) => 1,         // -1
            _ => panic!("Invalid phase format"),
        };
        // Remaining characters are the Pauli string
        let s: String = s.as_ref().chars().skip((phase / 2 + 1) as usize).collect();
        let words = PauliWord::from(s);
        Self { word: words, phase }
    }
}

impl<A: PauliStorage, H: BuildHasher + Default + Clone> std::fmt::Display for PhasedPauliWord<A, H> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.phase {
            0 => write!(f, "+")?,
            1 => write!(f, "-")?,
            2 => write!(f, "+i")?,
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
        let ps: PhasedPauliWord<u64> = "-XYZI".to_string().into();
        assert_eq!(ps.word.n_qubits(), 4);
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        assert_eq!(ps.phase, 1);
        let ps: PhasedPauliWord<u64> = "+iXYZI".to_string().into();
        assert_eq!(ps.word.n_qubits(), 4);
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        assert_eq!(ps.phase, 2);
        let ps: PhasedPauliWord<u64> = "-iXYZI".to_string().into();
        assert_eq!(ps.word.n_qubits(), 4);
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        assert_eq!(ps.phase, 3);
    }
}
