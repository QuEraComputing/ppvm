// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod clifford;
mod data;
mod product;

// ppvm-phased-pauli-word-2: The phased Pauli word wrapper.
//
// `Phased<W>` wraps a base word `W` with an explicit `ℤ₄` phase — the standalone
// phased Pauli operator (`i^φ · g(w)`). Unlike the bare PauliWord / LossyPauliWord,
// whose clifford drops every conjugation sign, Phased<W> carries a fused clifford impl
// that tracks the ℤ₄ sign (`HYH = −Y`, `SXS† = Y`, the CNOT/CZ deltas).

pub use data::Phased;
pub type PhasedPauliWord = Phased<ppvm_pauli_word_2::PauliWord>;
pub use ppvm_pauli_word_2::PauliWord;
