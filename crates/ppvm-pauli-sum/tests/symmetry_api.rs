// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::Complex;
use ppvm_pauli_sum::symmetry::{
    TranslationGroup, canonicalize_pauli_sum, canonicalize_pauli_sum_complex, check_momentum_sector,
};
use ppvm_pauli_word::word::PauliWord;

type W = PauliWord<[u8; 1], fxhash::FxBuildHasher, true>;

#[test]
fn public_symmetry_imports_remain_available() {
    let group = TranslationGroup::chain_1d(2);

    let mut real_basis: Vec<W> = vec![W::from("XI"), W::from("IX")];
    let mut real_coeffs = vec![1.0, 1.0];
    canonicalize_pauli_sum(&mut real_basis, &mut real_coeffs, &group);
    assert_eq!(real_basis.len(), 1);

    let mut complex_basis: Vec<W> = vec![W::from("ZI"), W::from("IZ")];
    let mut complex_coeffs = vec![Complex::new(1.0, 0.0); 2];
    assert!(check_momentum_sector(&complex_basis, &complex_coeffs, &group, &[0], 1e-12,).is_ok());
    canonicalize_pauli_sum_complex(&mut complex_basis, &mut complex_coeffs, &group, &[0]);
    assert_eq!(complex_basis.len(), 1);
}
