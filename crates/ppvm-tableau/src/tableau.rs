use ppvm_runtime::traits::Clifford;

use crate::simd::{
    xor_phase_with, xor_phase_with_scalar, xor_phase_with_xor, xor_phase_with_xor_scalar,
};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Tableau {
    pub(crate) n_qubits: usize,
    pub(crate) n_words: usize,
    pub(crate) x: Vec<u64>,
    pub(crate) z: Vec<u64>,
    pub(crate) phase: Vec<u64>,
}

impl Tableau {
    pub fn new(n_qubits: usize) -> Self {
        let n_words = words_for(n_qubits);
        let x = vec![0u64; n_qubits * n_words];
        let mut z = vec![0u64; n_qubits * n_words];
        for i in 0..n_qubits {
            set_column_bit(&mut z, n_words, i, i, true);
        }
        let phase = vec![0u64; n_words];
        Self {
            n_qubits,
            n_words,
            x,
            z,
            phase,
        }
    }

    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    pub fn x(&mut self, qubit: usize) {
        let base = qubit * self.n_words;
        let z = &self.z[base..base + self.n_words];
        xor_phase_with(&mut self.phase, z);
    }

    pub fn x_scalar(&mut self, qubit: usize) {
        let base = qubit * self.n_words;
        let z = &self.z[base..base + self.n_words];
        xor_phase_with_scalar(&mut self.phase, z);
    }

    pub fn y(&mut self, qubit: usize) {
        let base = qubit * self.n_words;
        let x = &self.x[base..base + self.n_words];
        let z = &self.z[base..base + self.n_words];
        xor_phase_with_xor(&mut self.phase, x, z);
    }

    pub fn y_scalar(&mut self, qubit: usize) {
        let base = qubit * self.n_words;
        let x = &self.x[base..base + self.n_words];
        let z = &self.z[base..base + self.n_words];
        xor_phase_with_xor_scalar(&mut self.phase, x, z);
    }

    pub fn z(&mut self, qubit: usize) {
        let base = qubit * self.n_words;
        let x = &self.x[base..base + self.n_words];
        xor_phase_with(&mut self.phase, x);
    }

    pub fn z_scalar(&mut self, qubit: usize) {
        let base = qubit * self.n_words;
        let x = &self.x[base..base + self.n_words];
        xor_phase_with_scalar(&mut self.phase, x);
    }

    pub fn h(&mut self, qubit: usize) {
        let base = qubit * self.n_words;
        for word in 0..self.n_words {
            let idx = base + word;
            let x = self.x[idx];
            let z = self.z[idx];
            self.phase[word] ^= x & z;
            self.x[idx] = z;
            self.z[idx] = x;
        }
    }

    pub fn h_scalar(&mut self, qubit: usize) {
        let base = qubit * self.n_words;
        for word in 0..self.n_words {
            let idx = base + word;
            let x = self.x[idx];
            let z = self.z[idx];
            self.phase[word] ^= x & z;
            self.x[idx] = z;
            self.z[idx] = x;
        }
    }

    pub fn s(&mut self, qubit: usize) {
        let base = qubit * self.n_words;
        for word in 0..self.n_words {
            let idx = base + word;
            let x = self.x[idx];
            let z = self.z[idx];
            self.phase[word] ^= x & z;
            self.z[idx] = z ^ x;
        }
    }

    pub fn s_scalar(&mut self, qubit: usize) {
        let base = qubit * self.n_words;
        for word in 0..self.n_words {
            let idx = base + word;
            let x = self.x[idx];
            let z = self.z[idx];
            self.phase[word] ^= x & z;
            self.z[idx] = z ^ x;
        }
    }
}

impl Clifford for Tableau {
    fn x(&mut self, index: usize) {
        Tableau::x(self, index);
    }

    fn y(&mut self, index: usize) {
        Tableau::y(self, index);
    }

    fn z(&mut self, index: usize) {
        Tableau::z(self, index);
    }

    fn h(&mut self, index: usize) {
        Tableau::h(self, index);
    }

    fn s(&mut self, index: usize) {
        Tableau::s(self, index);
    }

    fn cnot(&mut self, _control: usize, _target: usize) {
        unimplemented!("CNOT not yet implemented for Tableau");
    }

    fn cz(&mut self, _control: usize, _target: usize) {
        unimplemented!("CZ not yet implemented for Tableau");
    }
}

fn words_for(n_qubits: usize) -> usize {
    (n_qubits + 63) / 64
}

fn set_column_bit(words: &mut [u64], n_words: usize, col: usize, row: usize, value: bool) {
    let word = row / 64;
    let mask = 1u64 << (row % 64);
    let idx = col * n_words + word;
    if value {
        words[idx] |= mask;
    } else {
        words[idx] &= !mask;
    }
}
