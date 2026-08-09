// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Gate-surface parity between the hash and columnar storage backends.

use num::Complex;
use ppvm_pauli_sum_2::{
    ColumnPauliSum, ColumnStore, HashMapStore, LossyPauliSum, LossyPauliWord, NoPolicy, PauliSum,
    PauliWord, Sum,
};
use ppvm_traits_2::{
    AmplitudeDamping, CliffordBatch, CliffordExtensions, CliffordExtensionsBatch,
    CorrelatedLossChannel, Depolarizing, Depolarizing2, LossChannel, Projection, ResetLossChannel,
    RotXY, RotationTwo, TwoQubitPauliError,
};
use rand::SeedableRng;
use rand::rngs::SmallRng;

/// The `-2` sum backends are density-matrix-like: every channel scales
/// coefficients analytically and never draws. The injected RNG is threaded only
/// to satisfy the trait surface, so a fixed seed suffices everywhere here.
fn rng() -> SmallRng {
    SmallRng::seed_from_u64(0)
}

fn sorted<S>(sum: &S) -> Vec<(String, f64)>
where
    S: Terms,
{
    let mut out = sum.terms();
    out.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    out
}

trait Terms {
    fn terms(&self) -> Vec<(String, f64)>;
}

impl<Store> Terms for Sum<Store, NoPolicy>
where
    Store: ppvm_traits_2::Accumulate<Coeff = f64>,
    Store::Key: ppvm_traits_2::Word + ppvm_traits_2::Indexable + std::fmt::Display,
{
    fn terms(&self) -> Vec<(String, f64)> {
        self.iter()
            .map(|(key, coeff)| (key.to_string(), coeff))
            .collect()
    }
}

fn assert_close(a: &impl Terms, b: &impl Terms) {
    let (a, b) = (sorted(a), sorted(b));
    assert_eq!(a.len(), b.len(), "{a:?}\n{b:?}");
    for ((ak, av), (bk, bv)) in a.iter().zip(&b) {
        assert_eq!(ak, bk);
        assert!(
            (av - bv).abs() <= 1e-12 * av.abs().max(1.0),
            "{ak}: {av} vs {bv}"
        );
    }
}

fn ordinary_pair() -> (PauliSum, ColumnPauliSum) {
    let terms = [
        ("IXYZ", 1.0),
        ("XYZI", -0.25),
        ("ZZII", 0.75),
        ("IIII", 2.0),
    ];
    (
        PauliSum::from_terms(4, terms.map(|(w, c)| (PauliWord::from(w), c))),
        ColumnPauliSum::from_terms(4, terms.map(|(w, c)| (PauliWord::from(w), c))),
    )
}

#[test]
fn complete_gate_surface_matches_between_backends() {
    let (mut hash, mut column) = ordinary_pair();
    let mut rng = rng();
    macro_rules! both {
        ($method:ident($($arg:expr),* $(,)?)) => {{
            hash.$method($($arg),*);
            column.$method($($arg),*);
            assert_close(&hash, &column);
        }};
    }

    both!(s_dag(0));
    both!(sqrt_x(1));
    both!(sqrt_y_dag(2));
    both!(cy(0, 3));
    both!(h_many(&[0, 2]));
    both!(sqrt_x_many(&[1, 3]));
    both!(cy_many(&[(0, 1), (2, 3)]));
    both!(rotate_2([1, 0], [0, 1], 0, 2, 0.31));
    both!(rxx(1, 3, -0.27));
    both!(ryy(0, 2, 0.19));
    both!(rzz(2, 3, 0.41));
    both!(r(1, 0.7, 0.23));
    both!(two_qubit_pauli_error(
        0,
        3,
        std::array::from_fn(|i| 0.001 * (i + 1) as f64),
        &mut rng
    ));
    both!(depolarize1(2, 0.13, &mut rng));
    both!(depolarize2(0, 1, 0.09, &mut rng));
    both!(amplitude_damping(3, 0.17));
    both!(p0(0));
    both!(p1(2));
}

#[test]
fn hermitian_overlap_matches_between_backends() {
    type H = Sum<HashMapStore<PauliWord, Complex<f64>>, NoPolicy>;
    type C = Sum<ColumnStore<PauliWord, Complex<f64>>, NoPolicy>;
    let terms = [
        (PauliWord::from("XI"), Complex::new(1.0, 2.0)),
        (PauliWord::from("IZ"), Complex::new(-0.5, 0.25)),
    ];
    let hash = H::from_terms(2, terms);
    let column = C::from_terms(2, terms);
    assert_eq!(
        hash.hermitian_overlap(&hash),
        column.hermitian_overlap(&column)
    );
}

#[test]
fn lossy_channels_match_on_columnar_lossy_words() {
    type ColumnLossy = Sum<ColumnStore<LossyPauliWord, f64>, NoPolicy>;
    let terms = [("LI", 1.0), ("IL", 2.0), ("ZZ", 3.0), ("II", -0.5)];
    let mut hash: LossyPauliSum =
        Sum::from_terms(2, terms.map(|(w, c)| (LossyPauliWord::from(w), c)));
    let mut column: ColumnLossy =
        Sum::from_terms(2, terms.map(|(w, c)| (LossyPauliWord::from(w), c)));

    let mut rng = rng();
    hash.reset_loss_channel(0);
    column.reset_loss_channel(0);
    hash.loss_channel(1, 0.13, &mut rng);
    column.loss_channel(1, 0.13, &mut rng);
    hash.correlated_loss_channel(0, 1, [0.03, 0.07, 0.11], &mut rng);
    column.correlated_loss_channel(0, 1, [0.03, 0.07, 0.11], &mut rng);
    assert_close(&hash, &column);
}
