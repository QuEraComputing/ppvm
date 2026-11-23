// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Shared helpers: deterministic random sums, a dense-matrix oracle, and
//! operator comparison that treats missing keys as zero.
#![allow(dead_code)]

use approx::AbsDiffEq;
use num::{Complex, Zero};
use ppvm_pauli_sum_2::{PauliSum, PauliWord, Sum};
use ppvm_traits_2::{PauliBits, Support};

pub type C = Complex<f64>;
pub type CSum = PauliSum<C>;

pub fn c(re: f64, im: f64) -> C {
    Complex::new(re, im)
}

/// A one-term sum with coefficient one.
pub fn word(w: &str) -> CSum {
    CSum::from_terms(w.len(), [(PauliWord::from(w), c(1.0, 0.0))])
}

/// A deterministic random complex sum; repeated words accumulate.
pub fn random_sum(n_sites: usize, n_terms: usize, seed: u64) -> CSum {
    let mut state = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1;
    let mut below = move |bound: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % bound
    };
    let mut sum = CSum::new(n_sites);
    for _ in 0..n_terms {
        let w: String = (0..n_sites)
            .map(|_| ['I', 'X', 'Y', 'Z'][below(4) as usize])
            .collect();
        let mut unit = || below(2001) as f64 / 1000.0 - 1.0;
        sum += (w.as_str(), c(unit(), unit()));
    }
    sum
}

/// Dense `2^n × 2^n` matrix, row-major; site 0 is the most significant factor.
pub fn dense(sum: &CSum) -> Vec<C> {
    let n = sum.n_sites();
    let dim = 1 << n;
    let one = c(1.0, 0.0);
    let i = c(0.0, 1.0);
    let o = C::zero();
    // Indexed by packed code: I, X, Z, Y.
    let paulis = [
        [[one, o], [o, one]],
        [[o, one], [one, o]],
        [[one, o], [o, -one]],
        [[o, -i], [i, o]],
    ];
    let mut out = vec![C::zero(); dim * dim];
    for (key, coeff) in sum.iter() {
        for row in 0..dim {
            for col in 0..dim {
                let entry = (0..n).fold(coeff, |acc, site| {
                    let bit = |x: usize| (x >> (n - 1 - site)) & 1;
                    acc * paulis[key.pauli_code(site) as usize][bit(row)][bit(col)]
                });
                out[row * dim + col] += entry;
            }
        }
    }
    out
}

pub fn matmul(a: &[C], b: &[C]) -> Vec<C> {
    let dim = (a.len() as f64).sqrt() as usize;
    let mut out = vec![C::zero(); dim * dim];
    for row in 0..dim {
        for k in 0..dim {
            for col in 0..dim {
                out[row * dim + col] += a[row * dim + k] * b[k * dim + col];
            }
        }
    }
    out
}

pub fn assert_matrix_eq(actual: &[C], expected: &[C]) {
    assert_eq!(actual.len(), expected.len());
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!((a - e).norm() < 1e-12, "entry {i}: {a} != {e}");
    }
}

/// Equal as operators: every key's coefficient agrees, absent keys counting as zero.
pub fn assert_same_operator<S, P1, P2>(actual: &Sum<S, P1>, expected: &Sum<S, P2>)
where
    S: Support,
    S::Key: std::fmt::Debug,
    S::Coeff: AbsDiffEq<Epsilon = f64> + std::fmt::Debug,
{
    assert_eq!(actual.n_sites(), expected.n_sites());
    let check = |key: &S::Key| {
        let a = actual.get(key).unwrap_or_else(S::Coeff::zero);
        let e = expected.get(key).unwrap_or_else(S::Coeff::zero);
        assert!(a.abs_diff_eq(&e, 1e-12), "{key:?}: {a:?} != {e:?}");
    };
    actual.for_each_ref(|key, _| check(key));
    expected.for_each_ref(|key, _| check(key));
}
