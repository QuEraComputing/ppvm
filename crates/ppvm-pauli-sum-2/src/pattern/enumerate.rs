// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::marker::PhantomData;

use ppvm_pauli_word_2::{PauliStorage, PauliWord};
use ppvm_traits_2::{Pauli, PauliBits};

use super::data::{Decorated, OpPattern, PauliPattern};

impl PauliPattern {
    /// Enumerate matching words at a fixed width.
    ///
    /// As in the original implementation, any star is rejected with a panic:
    /// a finite word width does not implicitly bind an unbounded pattern.
    pub fn enumerate_matches<A: PauliStorage>(
        &self,
        n_qubits: usize,
    ) -> EnumMatchesPauliPattern<A> {
        let mut choices = Vec::new();
        let mut start = 0;
        for decorated in &self.0 {
            match decorated {
                Decorated::Position(op, position) => {
                    choices.extend((start..*position).map(|_| vec![Pauli::I]));
                    choices.push(op_choices(*op));
                    start = *position + 1;
                }
                Decorated::Repeat(op, count) => {
                    choices.extend((0..*count).map(|_| op_choices(*op)));
                }
                Decorated::Star(_) => panic!("Star patterns are not supported"),
            }
        }
        choices.extend((choices.len()..n_qubits).map(|_| vec![Pauli::I]));
        EnumMatchesPauliPattern {
            n_qubits,
            indices: vec![0; choices.len()],
            choices,
            first: true,
            done: false,
            storage: PhantomData,
        }
    }
}

fn op_choices(pattern: OpPattern) -> Vec<Pauli> {
    match pattern {
        OpPattern::Identity => vec![Pauli::I],
        OpPattern::Single(pauli) => vec![pauli],
        OpPattern::Double(a, b) => vec![a, b],
        OpPattern::AnyNonIdentity => vec![Pauli::X, Pauli::Z, Pauli::Y],
        OpPattern::SingleOrIdentity(pauli) => vec![pauli, Pauli::I],
        OpPattern::DoubleOrIdentity(a, b) => vec![a, b, Pauli::I],
        // Old starts this iterator in encoded order. Unlike old's missing stop
        // condition, this finite implementation ends after the four valid values.
        OpPattern::AnyPauliOrIdentity => vec![Pauli::I, Pauli::X, Pauli::Z, Pauli::Y],
    }
}

/// Iterator over the concrete words accepted by a bounded pattern.
#[derive(Debug, Clone)]
pub struct EnumMatchesPauliPattern<A: PauliStorage> {
    n_qubits: usize,
    choices: Vec<Vec<Pauli>>,
    indices: Vec<usize>,
    first: bool,
    done: bool,
    storage: PhantomData<A>,
}

impl<A: PauliStorage> Iterator for EnumMatchesPauliPattern<A> {
    type Item = PauliWord<A>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        if self.first {
            self.first = false;
        } else if !self.advance() {
            self.done = true;
            return None;
        }

        let mut word = PauliWord::new(self.n_qubits);
        for (position, (&choice, options)) in self.indices.iter().zip(&self.choices).enumerate() {
            let (x, z) = match options[choice] {
                Pauli::I => (false, false),
                Pauli::X => (true, false),
                Pauli::Z => (false, true),
                Pauli::Y => (true, true),
            };
            word.set_x_bit(position, x);
            word.set_z_bit(position, z);
        }
        Some(word)
    }
}

impl<A: PauliStorage> EnumMatchesPauliPattern<A> {
    fn advance(&mut self) -> bool {
        for position in (0..self.indices.len()).rev() {
            self.indices[position] += 1;
            if self.indices[position] < self.choices[position].len() {
                return true;
            }
            self.indices[position] = 0;
        }
        false
    }
}
