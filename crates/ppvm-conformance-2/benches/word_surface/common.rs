// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

pub type Storage = [u8; 32];

pub type NewWord = ppvm_pauli_word_2::PauliWord<Storage>;
pub type OldWord = ppvm_pauli_word::word::PauliWord<Storage>;
pub type NewLossy = ppvm_lossy_pauli_word_2::LossyPauliWord<Storage>;
pub type OldLossy = ppvm_pauli_word::loss::LossyPauliWord<Storage>;
pub type NewPhased = ppvm_phased_pauli_word_2::Phased<NewWord>;
pub type OldPhased = ppvm_pauli_word::phase::PhasedPauliWord<Storage>;

pub const WIDTH: usize = 256;
pub const SITE: usize = 127;
pub const SITE2: usize = 191;

pub fn ordinary_string(width: usize) -> String {
    "IXYZ".chars().cycle().take(width).collect()
}

pub fn lossy_string(width: usize) -> String {
    "IXYZL".chars().cycle().take(width).collect()
}

pub fn phased_string(width: usize) -> String {
    format!("+i{}", ordinary_string(width))
}

pub fn old_pauli(p: ppvm_traits_2::Pauli) -> ppvm_traits::char::Pauli {
    match p {
        ppvm_traits_2::Pauli::I => ppvm_traits::char::Pauli::I,
        ppvm_traits_2::Pauli::X => ppvm_traits::char::Pauli::X,
        ppvm_traits_2::Pauli::Y => ppvm_traits::char::Pauli::Y,
        ppvm_traits_2::Pauli::Z => ppvm_traits::char::Pauli::Z,
    }
}
