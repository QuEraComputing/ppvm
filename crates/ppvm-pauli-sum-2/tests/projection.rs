// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod common;

use common::{CSum, assert_same_operator, c, random_sum};
use ppvm_pauli_sum_2::{CoefficientThreshold, LossyPauliSum, LossyPauliWord, PauliSum};
use ppvm_traits_2::{PauliBits, Projection};

type Terms = &'static [(&'static str, f64)];

#[test]
fn projectors_map_each_site_value() {
    let cases: [(&str, Terms, Terms); 5] = [
        ("I", &[("I", 0.5), ("Z", 0.5)], &[("I", 0.5), ("Z", -0.5)]),
        ("Z", &[("I", 0.5), ("Z", 0.5)], &[("I", -0.5), ("Z", 0.5)]),
        ("X", &[("X", 0.0)], &[("X", 0.0)]),
        ("Y", &[("Y", 0.0)], &[("Y", 0.0)]),
        ("L", &[("L", 0.0)], &[("L", 0.0)]),
    ];
    for (word, p0, p1) in cases {
        for (project, expected) in [(0, p0), (1, p1)] {
            let mut sum = LossyPauliSum::<f64>::new(1);
            sum += (word, 1.0);
            if project == 0 {
                sum.p0(0)
            } else {
                sum.p1(0)
            }
            let expected: LossyPauliSum = LossyPauliSum::from_terms(
                1,
                expected.iter().map(|&(w, v)| (LossyPauliWord::from(w), v)),
            );
            assert_eq!(sum, expected, "p{project}({word})");
        }
    }
}

/// Every two-site lossy word, with distinct coefficients.
fn all_lossy_words() -> LossyPauliSum {
    let mut sum = LossyPauliSum::<f64>::new(2);
    for (i, a) in "IXYZL".chars().enumerate() {
        for (j, b) in "IXYZL".chars().enumerate() {
            sum += (format!("{a}{b}").as_str(), (5 * i + j + 1) as f64);
        }
    }
    sum
}

#[test]
fn projectors_are_idempotent_orthogonal_and_complete() {
    let o = all_lossy_words();
    for q in 0..2 {
        let project = |sum: &LossyPauliSum, one: bool| {
            let mut out = sum.clone();
            if one {
                out.p1(q)
            } else {
                out.p0(q)
            }
            out
        };
        let (p0, p1) = (project(&o, false), project(&o, true));
        assert_same_operator(&project(&p0, false), &p0);
        assert_same_operator(&project(&p1, true), &p1);
        assert!(project(&p0, true).iter().all(|(_, v)| v == 0.0));
        assert!(project(&p1, false).iter().all(|(_, v)| v == 0.0));
        // P0 O P0 + P1 O P1 dephases the site: X, Y, and lost terms vanish.
        let dephased: LossyPauliSum =
            LossyPauliSum::from_terms(2, o.iter().filter(|(k, _)| !k.is_lost(q) && !k.x_bit(q)));
        assert_same_operator(&(&p0 + &p1), &dephased);
    }
}

#[test]
fn projection_matches_conjugation_by_the_projector() {
    let o = random_sum(3, 12, 9);
    for q in 0..3 {
        for one in [false, true] {
            let z: String = (0..3).map(|i| if i == q { 'Z' } else { 'I' }).collect();
            let sign = if one { -0.5 } else { 0.5 };
            let mut projector = CSum::new(3);
            projector += ("III", c(0.5, 0.0));
            projector += (z.as_str(), c(sign, 0.0));
            let mut actual = o.clone();
            if one {
                actual.p1(q)
            } else {
                actual.p0(q)
            }
            assert_same_operator(&actual, &(&(&projector * &o) * &projector));
        }
    }
}

#[test]
fn projection_keeps_zero_terms_and_ignores_policy() {
    let mut sum = PauliSum::<f64, _>::with_policy(2, CoefficientThreshold { threshold: 10.0 });
    sum += ("XI", 1.0);
    sum += ("ZI", 1.0);
    sum.p0(0);
    assert_eq!(sum.len(), 3);
    assert_eq!(sum.get(&"XI".into()), Some(0.0));
    sum.reduce();
    assert_eq!(sum.len(), 2);
    sum.truncate();
    assert!(sum.is_empty());
}

#[test]
fn invalid_projection_target_panics_before_mutation() {
    use std::panic::{AssertUnwindSafe, catch_unwind};
    let mut sum = PauliSum::<f64>::new(1);
    sum += ("Z", 1.0);
    let original = sum.clone();
    assert!(catch_unwind(AssertUnwindSafe(|| sum.p0(1))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| sum.p1(1))).is_err());
    assert_eq!(sum, original);
}
