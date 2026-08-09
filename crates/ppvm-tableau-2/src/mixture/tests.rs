// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::complex::Complex64;
use ppvm_traits_2::Clifford;

use super::GeneralizedTableauMixture;
use super::fingerprint::fingerprint;

type Mixture = GeneralizedTableauMixture<u64, usize>;

#[test]
fn constructor_obeys_strict_sum_cutoff() {
    assert_eq!(Mixture::new_with_seed(1, 1e-12, 0.999, 1).len(), 1);
    assert!(Mixture::new_with_seed(1, 1e-12, 1.0, 1).is_empty());
    assert!(Mixture::new_with_seed(1, 1e-12, 2.0, 1).is_empty());
}

#[test]
fn fingerprint_collision_still_checks_full_structure() {
    let mut mixture = Mixture::new_with_seed(2, 1e-12, 0.0, 3);
    let colliding_fp = mixture.fingerprints[0];
    let mut different = mixture.entries[0].0.clone();
    different.h(0);
    assert!(!mixture.insert_branches(vec![(different, 0.5, colliding_fp)]));
    assert_eq!(mixture.len(), 2);
}

#[test]
fn approximate_amplitudes_merge_without_eq_or_hash() {
    let mut mixture = Mixture::new_with_seed(1, 1e-6, 0.0, 5);
    let mut close = mixture.entries[0].0.clone();
    close
        .coefficients
        .mul_element_by(0, Complex64::new(1.0 + 1e-7, 0.0));
    let fp = fingerprint(&close);
    assert!(!mixture.insert_branches(vec![(close, 0.25, fp)]));
    assert_eq!(mixture.len(), 1);
    assert_eq!(mixture.entries[0].1, 1.25);
}

#[test]
fn record_scratch_and_probability_are_not_identity() {
    let mut mixture = Mixture::new_with_seed(1, 1e-12, 0.0, 7);
    let mut same_state = mixture.entries[0].0.fork();
    same_state.append_measurement_record(Some(true));
    let fp = fingerprint(&same_state);
    assert!(!mixture.insert_branches(vec![(same_state, 0.5, fp)]));
    assert_eq!(mixture.len(), 1);
    assert_eq!(mixture.entries[0].1, 1.5);
}

#[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
#[test]
fn parallel_sampler_matches_precomputed_serial_stream() {
    let mut serial_mixture = Mixture::new_with_seed(3, 1e-12, 0.0, 11);
    serial_mixture.h(0);
    serial_mixture.cnot(0, 1);
    let mut parallel_mixture = serial_mixture.clone();
    let serial = serial_mixture.sampler().sample_shots_serial(128);
    let parallel = parallel_mixture.sampler().sample_shots_parallel(128);
    assert_eq!(serial, parallel);
}
