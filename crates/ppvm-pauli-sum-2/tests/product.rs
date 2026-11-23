// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod common;

use common::{C, CSum, assert_matrix_eq, assert_same_operator, c, dense, matmul, random_sum, word};
use ppvm_pauli_sum_2::{CoefficientThreshold, PauliSum, PauliWord};

#[test]
fn single_qubit_products_follow_the_pauli_table() {
    let (one, i) = (c(1.0, 0.0), c(0.0, 1.0));
    let table = [
        ("X", "Y", "Z", i),
        ("Y", "Z", "X", i),
        ("Z", "X", "Y", i),
        ("Y", "X", "Z", -i),
        ("Z", "Y", "X", -i),
        ("X", "Z", "Y", -i),
        ("X", "X", "I", one),
        ("Y", "Y", "I", one),
        ("Z", "Z", "I", one),
        ("I", "Y", "Y", one),
    ];
    for (a, b, product, phase) in table {
        let result = &word(a) * &word(b);
        assert_eq!(result.len(), 1, "{a}{b}");
        assert_eq!(result.get(&PauliWord::from(product)), Some(phase), "{a}{b}");
    }
}

#[test]
fn products_match_dense_matrix_multiplication() {
    for seed in 1..=8 {
        let a = random_sum(3, 7, seed);
        let b = random_sum(3, 7, seed + 100);
        let product = &a * &b;
        assert_matrix_eq(&dense(&product), &matmul(&dense(&a), &dense(&b)));
        let mut assigned = a.clone();
        assigned *= &b;
        assert_eq!(assigned, product);
    }
}

#[test]
fn products_satisfy_ring_identities() {
    let (a, b, d) = (
        random_sum(3, 6, 1),
        random_sum(3, 6, 2),
        random_sum(3, 6, 3),
    );
    assert_same_operator(&(&(&a * &b) * &d), &(&a * &(&b * &d)));
    assert_same_operator(&(&a * &(&b + &d)), &(&(&a * &b) + &(&a * &d)));
    assert_same_operator(&(&(&a + &b) * &d), &(&(&a * &d) + &(&b * &d)));
    let s = c(0.3, -1.2);
    assert_same_operator(&(&(&a * s) * &b), &(&(&a * &b) * s));
    assert_same_operator(&(&a * &(&b * s)), &(&(&a * &b) * s));
    let identity = word("III");
    assert_eq!(&identity * &a, a);
    assert_eq!(&a * &identity, a);
    let zero = CSum::new(3);
    assert!((&a * &zero).is_empty() && (&zero * &a).is_empty());
}

#[test]
fn adjoint_reverses_products() {
    // Pauli words are Hermitian, so the adjoint conjugates coefficients only.
    let adjoint =
        |sum: &CSum| CSum::from_terms(sum.n_sites(), sum.iter().map(|(k, v)| (k, v.conj())));
    let (a, b) = (random_sum(3, 6, 4), random_sum(3, 6, 5));
    assert_same_operator(&adjoint(&(&a * &b)), &(&adjoint(&b) * &adjoint(&a)));
}

#[test]
fn products_keep_zeros_and_the_left_policy_without_truncating() {
    let mut a = PauliSum::<C, _>::with_policy(1, CoefficientThreshold { threshold: 10.0 });
    a += ("X", c(1.0, 0.0));
    a += ("Z", c(1.0, 0.0));
    a.preserve_key("X".into());
    let mut b = CSum::new(1);
    b += ("X", c(1.0, 0.0));
    b += ("Z", c(-1.0, 0.0));
    // (X + Z)(X - Z) = (I - I) + (ZX - XZ) = 0·I + 2i·Y.
    let mut product = &a * &b;
    assert_eq!(product.len(), 2);
    assert_eq!(product.get(&"I".into()), Some(c(0.0, 0.0)));
    assert_eq!(product.get(&"Y".into()), Some(c(0.0, 2.0)));
    assert_eq!(product.policy().threshold, 10.0);
    assert!(product.preserved_keys().contains(&"X".into()));
    product.truncate();
    assert!(product.is_empty());
}

#[test]
#[should_panic(expected = "product requires equal widths")]
fn products_reject_mismatched_widths() {
    let _ = &word("XI") * &word("X");
}
