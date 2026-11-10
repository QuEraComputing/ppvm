// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The row product against the matrix-model site table, and group invariants.

mod common;
#[allow(dead_code)]
#[path = "generated/tables.rs"]
#[rustfmt::skip]
mod tables;

use common::Rng;
use ppvm_pauli_kernels::{row_product, row_product_phase, word_product};

/// A Pauli word on `64 · len` qubits as X and Z planes.
#[derive(Clone, PartialEq, Debug)]
struct Row {
    x: Vec<u64>,
    z: Vec<u64>,
}

impl Row {
    fn random(rng: &mut Rng, len: usize) -> Self {
        Self {
            x: rng.words(len),
            z: rng.words(len),
        }
    }

    /// `self ← self · other`, returning the phase exponent.
    fn mul_assign(&mut self, other: &Row) -> u8 {
        row_product(&mut self.x, &mut self.z, &other.x, &other.z)
    }
}

/// Row lengths in words, on both sides of the kernel's short/long-row switch (8 words).
const LENS: [usize; 7] = [1, 2, 5, 8, 9, 16, 64];

fn site_exponent(a: bool, b: bool, c: bool, d: bool) -> u8 {
    tables::SITE_PRODUCT
        .iter()
        .find(|(k, _)| *k == [a, b, c, d])
        .unwrap()
        .1
}

#[test]
fn single_site_matches_matrix_model() {
    for ([a, b, c, d], e) in tables::SITE_PRODUCT {
        let (mut x, mut z) = ([a as u8], [b as u8]);
        let k = row_product(&mut x, &mut z, &[c as u8], &[d as u8]);
        assert_eq!(k, e, "site {:?}", [a, b, c, d]);
        assert_eq!((x[0], z[0]), ((a ^ c) as u8, (b ^ d) as u8));
    }
}

/// One word's contribution, mod 4, is the sum of its 64 site exponents.
#[test]
fn word_product_matches_site_phases() {
    let mut rng = Rng::new(19);
    for _ in 0..256 {
        let [a, b, c, d] = [
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
        ];
        let expected: u32 = (0..64)
            .map(|i| {
                let bit = |w: u64| (w >> i) & 1 == 1;
                site_exponent(bit(a), bit(b), bit(c), bit(d)) as u32
            })
            .sum();
        let ([x, z], e) = word_product(a, b, c, d);
        assert_eq!(e % 4, expected % 4);
        assert_eq!([x, z], [a ^ c, b ^ d]);
    }
}

#[test]
fn phase_is_sum_of_site_phases() {
    let mut rng = Rng::new(3);
    for len in LENS {
        let (p, q) = (Row::random(&mut rng, len), Row::random(&mut rng, len));
        let mut expected = 0u32;
        for w in 0..len {
            for i in 0..64 {
                let bit = |v: &[u64]| (v[w] >> i) & 1 == 1;
                expected += site_exponent(bit(&p.x), bit(&p.z), bit(&q.x), bit(&q.z)) as u32;
            }
        }
        let mut prod = p.clone();
        let k = prod.mul_assign(&q);
        assert_eq!(k as u32, expected % 4);
        let xor = |a: &[u64], b: &[u64]| a.iter().zip(b).map(|(a, b)| a ^ b).collect::<Vec<_>>();
        assert_eq!(
            prod,
            Row {
                x: xor(&p.x, &q.x),
                z: xor(&p.z, &q.z)
            }
        );
        assert_eq!(row_product_phase(&p.x, &p.z, &q.x, &q.z), k);
    }
}

#[test]
fn product_is_associative_with_phase() {
    let mut rng = Rng::new(5);
    for len in LENS {
        let (a, b, c) = (
            Row::random(&mut rng, len),
            Row::random(&mut rng, len),
            Row::random(&mut rng, len),
        );
        let mut left = a.clone();
        let k1 = left.mul_assign(&b);
        let k2 = left.mul_assign(&c);
        let mut bc = b.clone();
        let k3 = bc.mul_assign(&c);
        let mut right = a.clone();
        let k4 = right.mul_assign(&bc);
        assert_eq!(left, right);
        assert_eq!((k1 + k2) % 4, (k3 + k4) % 4);
    }
}

#[test]
fn hermitian_paulis_square_to_identity() {
    let mut rng = Rng::new(9);
    for len in LENS {
        let p = Row::random(&mut rng, len);
        let mut sq = p.clone();
        assert_eq!(sq.mul_assign(&p), 0);
        assert!(sq.x.iter().chain(&sq.z).all(|&w| w == 0));
    }
}

#[test]
fn commutator_phase_is_symplectic_form() {
    let mut rng = Rng::new(13);
    for len in LENS {
        let (a, b) = (Row::random(&mut rng, len), Row::random(&mut rng, len));
        let omega: u32 = (0..len)
            .map(|w| ((a.x[w] & b.z[w]) ^ (a.z[w] & b.x[w])).count_ones())
            .sum::<u32>()
            % 2;
        let kab = row_product_phase(&a.x, &a.z, &b.x, &b.z) as u32;
        let kba = row_product_phase(&b.x, &b.z, &a.x, &a.z) as u32;
        assert_eq!((kab + 4 - kba) % 4, 2 * omega);
    }
}

/// Zero words change nothing. Padding a row from 8 to 9 or 16 words also moves it
/// from the short-row loop to the long-row loop, so this compares the two.
#[test]
fn padding_bits_contribute_nothing() {
    let mut rng = Rng::new(17);
    for len in LENS {
        let (p, q) = (Row::random(&mut rng, len), Row::random(&mut rng, len));
        let k = row_product_phase(&p.x, &p.z, &q.x, &q.z);
        for padded in [len + 1, 2 * len, len + 9] {
            let pad = |v: &[u64]| {
                let mut w = v.to_vec();
                w.resize(padded, 0);
                w
            };
            let (mut px, mut pz) = (pad(&p.x), pad(&p.z));
            assert_eq!(row_product_phase(&px, &pz, &pad(&q.x), &pad(&q.z)), k);
            assert_eq!(row_product(&mut px, &mut pz, &pad(&q.x), &pad(&q.z)), k);
            assert!(px[len..].iter().chain(&pz[len..]).all(|&w| w == 0));
        }
    }
}
