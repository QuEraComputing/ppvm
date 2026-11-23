// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod common;

use approx::assert_abs_diff_eq;
use common::{dense, random_sum};
use ppvm_pauli_sum_2::{LossyPauliSum, LossyPauliWord, Pattern, PauliSum, SiteSet};
use ppvm_traits_2::{Clifford, Depolarizing, LossChannel, Projection, RotationOne};

#[test]
fn patterns_parse_and_match() {
    let pattern: Pattern = "[IZ]_XL".parse().unwrap();
    assert_eq!(pattern.n_sites(), 4);
    let sites = vec![
        SiteSet::I | SiteSet::Z,
        SiteSet::ANY,
        SiteSet::X,
        SiteSet::L,
    ];
    assert_eq!(pattern, Pattern::from(sites));
    for (word, matches) in [
        ("ZYXL", true),
        ("ILXL", true),
        ("XYXL", false),
        ("IIXI", false),
    ] {
        let key: LossyPauliWord = word.into();
        assert_eq!(pattern.matches(&key), matches, "{word}");
    }
    assert_eq!(Pattern::zero_state(2), "[IZ][IZ]".parse().unwrap());
    for bad in ["Q", "[XQ]", "[X", "[]", "X]"] {
        assert!(bad.parse::<Pattern>().is_err(), "{bad}");
    }
}

#[test]
fn trace_sums_matching_coefficients() {
    let mut sum = LossyPauliSum::<f64>::new(2);
    for (word, coeff) in [
        ("II", 1.0),
        ("ZI", 2.0),
        ("ZZ", 4.0),
        ("XI", 8.0),
        ("LZ", 16.0),
    ] {
        sum += (word, coeff);
    }
    assert_eq!(sum.trace(&Pattern::zero_state(2)), 7.0);
    assert_eq!(sum.trace(&"[IX]I".parse().unwrap()), 9.0);
    assert_eq!(sum.trace(&"L_".parse().unwrap()), 16.0);
    assert_eq!(sum.trace(&"YY".parse().unwrap()), 0.0);
    assert_eq!(
        LossyPauliSum::<f64>::new(2).trace(&Pattern::zero_state(2)),
        0.0
    );
}

#[test]
#[should_panic(expected = "pattern width must match sum width")]
fn trace_rejects_mismatched_widths() {
    PauliSum::<f64>::new(2).trace(&Pattern::zero_state(3));
}

#[test]
fn basis_state_readout_matches_the_dense_diagonal() {
    let o = random_sum(3, 16, 11);
    let matrix = dense(&o);
    for b in 0..8 {
        // <b|O|b> = <0|X^b O X^b|0>.
        let mut flipped = o.clone();
        for q in (0..3).filter(|q| b >> (2 - q) & 1 == 1) {
            flipped.x(q);
        }
        let readout = flipped.trace(&Pattern::zero_state(3));
        assert!((readout - matrix[b * 8 + b]).norm() < 1e-12, "b = {b:03b}");
    }
}

/// Observables propagate backwards: apply the circuit's gates in reverse order.
#[test]
fn propagation_workloads_read_out_known_values() {
    let zero = |n| Pattern::zero_state(n);
    let observable = |w: &str| PauliSum::<f64>::from_terms(w.len(), [(w.into(), 1.0)]);

    let theta = 0.7;
    let mut z = observable("Z");
    z.rx(0, theta);
    assert_abs_diff_eq!(z.trace(&zero(1)), theta.cos(), epsilon = 1e-12);

    // GHZ: H(0), CNOT(0, 1), CNOT(1, 2).
    for (w, expected) in [
        ("ZZI", 1.0),
        ("IZZ", 1.0),
        ("XXX", 1.0),
        ("YYX", -1.0),
        ("ZII", 0.0),
    ] {
        let mut o = observable(w);
        o.cnot(1, 2);
        o.cnot(0, 1);
        o.h(0);
        assert_eq!(o.trace(&zero(3)), expected, "{w}");
    }

    let p = 0.3;
    let mut z = observable("Z");
    z.depolarize1(0, p);
    assert_abs_diff_eq!(z.trace(&zero(1)), 1.0 - 4.0 * p / 3.0, epsilon = 1e-12);

    // Probability of outcome 0 after H on |0>: the observable is P0 itself.
    let mut born = observable("I");
    born.p0(0);
    born.h(0);
    assert_abs_diff_eq!(born.trace(&zero(1)), 0.5, epsilon = 1e-12);

    // Loss from a live |0>: L reads out p, Z reads out 1 - p; a lost start reads L as 1.
    let lossy = |w: &str| {
        let mut sum = LossyPauliSum::<f64>::from_terms(1, [(w.into(), 1.0)]);
        sum.loss_channel(0, p);
        sum
    };
    assert_abs_diff_eq!(lossy("L").trace(&zero(1)), p, epsilon = 1e-12);
    assert_abs_diff_eq!(lossy("Z").trace(&zero(1)), 1.0 - p, epsilon = 1e-12);
    assert_eq!(lossy("L").trace(&Pattern::uniform(1, SiteSet::L)), 1.0);
}
