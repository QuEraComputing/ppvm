// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::complex::Complex64;
use ppvm_pauli_sum::config::fxhash::ByteF64;
use ppvm_tableau_2::GeneralizedTableauMixture;
use ppvm_tableau_sum::data::GeneralizedTableauSum as OldMixture;
use ppvm_tableau_sum::storage::EntryStore;
use ppvm_tableau_sum::storage::vec::VecStorage;

pub type Old = OldMixture<
    ByteF64<2>,
    u128,
    Vec<(Complex64, u128)>,
    VecStorage<ByteF64<2>, u128, Vec<(Complex64, u128)>>,
>;
pub type New = GeneralizedTableauMixture<[u8; 2], u128>;

#[derive(Clone, Debug)]
pub struct StateSnapshot {
    pub rows: Vec<([u8; 2], [u8; 2], u8)>,
    pub amplitudes: Vec<(u128, Complex64)>,
    pub loss: Vec<bool>,
    pub probability: f64,
}

pub fn old(seed: u64, cutoff: f64) -> Old {
    Old::new_with_seed(12, 1e-12, cutoff, seed)
}

pub fn new(seed: u64, cutoff: f64) -> New {
    New::new_with_seed(12, 1e-12, cutoff, seed)
}

pub fn old_snapshot(mixture: &Old) -> Vec<StateSnapshot> {
    mixture
        .entries
        .iter()
        .map(|(tab, probability)| {
            let rows = tab
                .tableau
                .data
                .iter()
                .map(|row| (row.word.xbits.data, row.word.zbits.data, row.phase))
                .collect();
            let mut amplitudes: Vec<_> = tab
                .coefficients
                .iter()
                .map(|(value, index)| (*index, *value))
                .collect();
            amplitudes.sort_by_key(|entry| entry.0);
            StateSnapshot {
                rows,
                amplitudes,
                loss: tab.is_lost.clone(),
                probability: *probability,
            }
        })
        .collect()
}

pub fn new_snapshot(mixture: &New) -> Vec<StateSnapshot> {
    mixture
        .iter()
        .map(|(tab, probability)| {
            let rows = tab.tableau.rows().collect();
            let mut amplitudes: Vec<_> = tab
                .coefficients
                .iter()
                .map(|(value, index)| (*index, *value))
                .collect();
            amplitudes.sort_by_key(|entry| entry.0);
            StateSnapshot {
                rows,
                amplitudes,
                loss: tab.is_lost.clone(),
                probability: *probability,
            }
        })
        .collect()
}

pub fn assert_snapshots_close(mut old: Vec<StateSnapshot>, mut new: Vec<StateSnapshot>) {
    old.sort_by_key(|state| format!("{:?}{:?}", state.rows, state.loss));
    new.sort_by_key(|state| format!("{:?}{:?}", state.rows, state.loss));
    assert_eq!(old.len(), new.len());
    for (old, new) in old.iter().zip(new.iter()) {
        assert_eq!(old.rows, new.rows);
        assert_eq!(old.loss, new.loss);
        assert_eq!(old.amplitudes.len(), new.amplitudes.len());
        assert!((old.probability - new.probability).abs() < 1e-11);
        for ((oi, ov), (ni, nv)) in old.amplitudes.iter().zip(&new.amplitudes) {
            assert_eq!(oi, ni);
            assert!((*ov - *nv).norm() < 1e-11);
        }
    }
}
