// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The Phase-6 gate for the columnar backend: **a backend swap must be
//! observationally identical** (implementation-plan §"Phase 6 — `ColumnStore`
//! (SoA backend)").
//!
//! Every test here runs the *same* driver body twice — once on [`PauliSum`] (the
//! `HashMapStore` hash-join backend) and once on [`ColumnPauliSum`] (the
//! `ColumnStore` structure-of-arrays backend) — and compares the full support
//! term for term. The body is duplicated by a macro rather than made generic on
//! purpose: a generic driver would need a `where` clause naming every capability
//! trait, and the point is that *user* code written against the concrete alias
//! behaves the same either way.
//!
//! The comparison is on the **sorted** `(word, coefficient)` list, since neither
//! backend specifies its iteration order (the hash map's is bucket order, the
//! column store's is insertion order). Coefficients are compared **bit-exactly**:
//! every path exercised here applies the same float operations to the same term
//! in the same order on both backends. The two places where that is genuinely
//! not true — the pairings (`overlap`) and the L4 product, whose *summation*
//! order follows the support order — are compared at a relative tolerance and
//! are called out where they appear.

use ppvm_pauli_sum_2::{
    CoefficientThreshold, ColumnPauliSum, CombinedPolicy, MaxPauliWeight, NoPolicy, PauliPattern,
    PauliSum, PauliWord, Policy, Sum,
};
use ppvm_traits_2::{
    Accumulate, Clifford, Indexable, PauliError, Projection, RotationOne, Trace, Word,
};

fn pw(s: &str) -> PauliWord {
    PauliWord::from(s)
}

/// The support as a sorted `(word, coefficient)` list — backend-order-free.
fn sorted<S, P>(sum: &Sum<S, P>) -> Vec<(String, f64)>
where
    S: Accumulate<Key = PauliWord, Coeff = f64>,
    P: Policy<PauliWord, f64>,
    PauliWord: Word + Indexable,
{
    let mut v: Vec<(String, f64)> = sum.iter().map(|(k, c)| (format!("{k}"), c)).collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

/// Bit-exact support equality.
fn assert_identical(a: &[(String, f64)], b: &[(String, f64)]) {
    assert_eq!(
        a.len(),
        b.len(),
        "support size diverged:\n hash  ={a:?}\n column={b:?}"
    );
    for ((ka, ca), (kb, cb)) in a.iter().zip(b.iter()) {
        assert_eq!(ka, kb, "key sets diverged:\n hash  ={a:?}\n column={b:?}");
        assert_eq!(
            ca.to_bits(),
            cb.to_bits(),
            "coefficient diverged at {ka}: {ca} vs {cb}"
        );
    }
}

/// Expand one driver body against both backends and return the two sorted
/// supports.
macro_rules! both_backends {
    (n = $n:expr, policy = $policy:expr, seed = $seed:expr, |$sum:ident| $body:block) => {{
        let hash = {
            let mut $sum: PauliSum<f64, _> = PauliSum::with_policy($n, $policy);
            for (w, c) in $seed {
                $sum += (pw(w), c);
            }
            $body
            sorted(&$sum)
        };
        let column = {
            let mut $sum: ColumnPauliSum<f64, _> = ColumnPauliSum::with_policy($n, $policy);
            for (w, c) in $seed {
                $sum += (pw(w), c);
            }
            $body
            sorted(&$sum)
        };
        (hash, column)
    }};
}

#[test]
fn clifford_and_rotation_propagation_is_identical() {
    let (a, b) = both_backends!(
        n = 4,
        policy = CoefficientThreshold { threshold: 1e-9 },
        seed = [("ZIII", 1.0), ("IZII", 1.0), ("IIZI", 1.0), ("IIIZ", 1.0)],
        |sum| {
            for q in 0..4 {
                sum.h(q);
                sum.rz(q, 1.1);
                sum.ry(q, 2.1);
                sum.s(q);
                sum.x(q);
                sum.y(q);
                sum.z(q);
            }
            for q in 0..3 {
                sum.cnot(q, q + 1);
                sum.rx(q, 0.37);
                sum.cz(q, q + 1);
            }
            sum.truncate();
            for q in 0..4 {
                sum.pauli_error(q, [1e-3, 2e-3, 3e-3]);
                sum.ry(q, 0.9);
            }
            sum.truncate();
        }
    );
    assert!(
        a.len() > 8,
        "the workload should build a real support, got {}",
        a.len()
    );
    assert_identical(&a, &b);
}

#[test]
fn truncation_and_zero_preservation_are_identical() {
    let (a, b) = both_backends!(
        n = 4,
        policy = CombinedPolicy(CoefficientThreshold { threshold: 1e-6 }, MaxPauliWeight(3)),
        seed = [
            ("XIII", 1.0),
            ("ZIII", 1.0),
            ("IZII", 1.0),
            ("IIZI", 1.0),
            ("IIIZ", 1.0)
        ],
        |sum| {
            // λ_X = 0 here, so the X-carrying term is scaled to exactly 0.0 and
            // MUST stay in the support on both backends (old has no `reduce`).
            sum.pauli_error(0, [0.0, 0.25, 0.25]);
            // Sub-threshold branches that must be allowed to merge before any
            // truncation sees them.
            sum.rx(1, 0.03);
            sum.rx(1, 0.03);
            sum.truncate();
            // An identity rotation still emits a 0.0-coefficient branch, which
            // both backends insert.
            sum.rz(2, 0.0);
            sum.truncate();
        }
    );
    assert_identical(&a, &b);
}

#[test]
fn zero_coefficients_survive_every_gate_on_both_backends() {
    // The sharpest no-implicit-reduce probe: no `truncate` anywhere, so the only
    // thing that could drop these terms is a backend deciding to. `pauli_error`
    // with `p_X = 0` gives λ_X = 0 (the X term is scaled to exactly 0.0) and the
    // identity rotation `rz(_, 0.0)` emits a 0.0-coefficient branch that old's
    // `map_insert` inserts.
    let (a, b) = both_backends!(
        n = 4,
        policy = NoPolicy,
        seed = [("XIII", 1.0), ("ZIII", 1.0), ("IIXI", 1.0)],
        |sum| {
            sum.pauli_error(0, [0.0, 0.25, 0.25]);
            sum.rz(2, 0.0);
        }
    );
    assert_identical(&a, &b);
    assert_eq!(
        b.iter().filter(|(_, c)| *c == 0.0).count(),
        2,
        "both the zeroed X term and the 0.0 identity-rotation branch must \
         survive on the columnar backend: {b:?}"
    );
    assert!(
        b.iter().any(|(k, c)| k == "XIII" && *c == 0.0),
        "the zeroed term keeps its key: {b:?}"
    );
}

#[test]
fn lean_correct_projection_is_identical_on_both_backends() {
    let (a, b) = both_backends!(
        n = 1,
        policy = NoPolicy,
        seed = [("I", 2.0), ("X", 3.0)],
        |sum| {
            sum.p0(0);
            sum.p0(0);
        }
    );
    assert_identical(&a, &b);
    assert_eq!(
        b,
        vec![("I".into(), 1.0), ("X".into(), 0.0), ("Z".into(), 1.0)]
    );
}

#[test]
fn deep_untruncated_fanout_is_identical() {
    // The random-circuit shape: no truncation, so the support grows by pure
    // fan-out — this is what stresses column growth and index rebuilds.
    let (a, b) = both_backends!(n = 6, policy = NoPolicy, seed = [("ZZIIII", 1.0)], |sum| {
        for _layer in 0..4 {
            for q in 0..6 {
                sum.rz(q, 1.1);
                sum.ry(q, 2.1);
                sum.rz(q, 1.1);
            }
            for q in 0..6 {
                sum.cnot(q, (q + 1) % 6);
            }
        }
    });
    assert!(
        a.len() > 200,
        "the fan-out workload should grow a large support, got {}",
        a.len()
    );
    assert_identical(&a, &b);
}

#[test]
fn preserve_set_restore_is_identical() {
    let seed = [
        ("ZII", 1e-6),
        ("IZI", 1e-6),
        ("IIZ", 1e-6),
        ("XYZ", 1e-6),
        ("XXX", 0.7),
    ];
    let policy = CoefficientThreshold { threshold: 0.5 };
    let keep = ["ZII", "IZI", "IIZ"].map(pw);

    let mut hash: PauliSum<f64, _> = PauliSum::with_policy(3, policy);
    let mut column: ColumnPauliSum<f64, _> = ColumnPauliSum::with_policy(3, policy);
    for (w, c) in seed {
        hash += (pw(w), c);
        column += (pw(w), c);
    }
    hash = hash.preserving(keep);
    column = column.preserving(keep);

    hash.truncate();
    column.truncate();
    assert_identical(&sorted(&hash), &sorted(&column));

    // Old's semantics: a preserved key is restored at its PRE-truncate value,
    // never doubled, and a non-preserved sub-threshold key is gone.
    assert_eq!(column.get(&pw("ZII")), Some(1e-6));
    assert_eq!(column.get(&pw("XYZ")), None);
    assert_eq!(column.get(&pw("XXX")), Some(0.7));
}

#[test]
fn reduce_is_caller_driven_on_the_columnar_backend() {
    let mut sum: ColumnPauliSum =
        ColumnPauliSum::from_terms(2, [(pw("XI"), 1.0), (pw("IZ"), 2.0), (pw("XI"), -1.0)]);
    // The exact cancellation stays until the caller asks.
    assert_eq!(sum.len(), 2);
    assert_eq!(sum.get(&pw("XI")), Some(0.0));
    sum.reduce();
    assert_eq!(sum.len(), 1);
    assert_eq!(sum.get(&pw("XI")), None);
    assert_eq!(sum.get(&pw("IZ")), Some(2.0));
}

#[test]
fn pairings_and_trace_agree_within_tolerance() {
    // The pairings sum in support order, which differs between the backends, so
    // the bar here is relative rather than bit-exact.
    let build_hash = || {
        let mut s: PauliSum = PauliSum::new(4);
        s += (pw("ZIII"), 1.0);
        s += (pw("IZII"), 1.0);
        s.cnot(0, 1);
        s.h(0);
        s.ry(2, 0.7);
        s
    };
    let build_column = || {
        let mut s: ColumnPauliSum = ColumnPauliSum::new(4);
        s += (pw("ZIII"), 1.0);
        s += (pw("IZII"), 1.0);
        s.cnot(0, 1);
        s.h(0);
        s.ry(2, 0.7);
        s
    };
    let (h, c) = (build_hash(), build_column());
    assert_identical(&sorted(&h), &sorted(&c));

    let hh = h.overlap(&build_hash());
    let cc = c.overlap(&build_column());
    assert!((hh - cc).abs() <= 1e-12 * hh.abs().max(1.0), "{hh} vs {cc}");

    let pattern = PauliPattern::zero_state();
    let th: f64 = h.trace(&pattern);
    let tc: f64 = c.trace(&pattern);
    assert!((th - tc).abs() <= 1e-12 * th.abs().max(1.0), "{th} vs {tc}");
}

#[test]
fn ghz_backward_trace_is_exactly_one_on_the_columnar_backend() {
    // The frozen scalar from old `tests/ghz.rs`, through the SoA backend.
    let mut sum: ColumnPauliSum = ColumnPauliSum::from_terms(2, [(pw("ZZ"), 1.0)]);
    sum.cnot(0, 1);
    sum.h(0);
    let t: f64 = sum.trace(&PauliPattern::zero_state());
    assert_eq!(t, 1.0);
}

/// The L4 operator product on the columnar backend: the fresh-accumulator form,
/// the aux-backed in-place form, and agreement with the hash-join backend.
///
/// Coefficients here are compared at a relative tolerance, not bit-exactly: the
/// convolution accumulates `|A|·|B|` contributions **in support order**, which
/// differs between the backends, so the float summation reassociates.
#[test]
fn operator_product_agrees_with_the_hash_backend() {
    use num::Complex;

    fn c(re: f64, im: f64) -> Complex<f64> {
        Complex::new(re, im)
    }

    let terms = [
        ("XIZ", c(0.5, 0.25)),
        ("IYZ", c(-1.5, 0.0)),
        ("ZZI", c(0.0, 2.0)),
        ("YXY", c(1.0, -1.0)),
    ];
    let rhs = [
        ("ZIX", c(2.0, 0.0)),
        ("IYI", c(0.0, -0.5)),
        ("XXX", c(1.0, 1.0)),
    ];

    let hash_a: PauliSum<Complex<f64>> = PauliSum::from_terms(3, terms.map(|(w, v)| (pw(w), v)));
    let hash_b: PauliSum<Complex<f64>> = PauliSum::from_terms(3, rhs.map(|(w, v)| (pw(w), v)));
    let col_a: ColumnPauliSum<Complex<f64>> =
        ColumnPauliSum::from_terms(3, terms.map(|(w, v)| (pw(w), v)));
    let col_b: ColumnPauliSum<Complex<f64>> =
        ColumnPauliSum::from_terms(3, rhs.map(|(w, v)| (pw(w), v)));

    let want = hash_a.multiply(&hash_b);
    let got = col_a.multiply(&col_b);
    assert_eq!(want.len(), got.len(), "product support size diverged");
    for (k, v) in want.iter() {
        let g = got.get(&k).unwrap_or_else(|| panic!("missing key {k}"));
        assert!(
            (v - g).norm() <= 1e-12 * v.norm().max(1.0),
            "{k}: {v} vs {g}"
        );
    }

    // The in-place form goes through the store's persistent aux double-buffer;
    // it must land on the same map, and leave the store reusable afterwards.
    let mut in_place = col_a.clone();
    in_place.multiply_in_place(&col_b);
    assert_eq!(in_place.len(), got.len());
    for (k, v) in got.iter() {
        let g = in_place
            .get(&k)
            .unwrap_or_else(|| panic!("missing key {k}"));
        assert!(
            (v - g).norm() <= 1e-12 * v.norm().max(1.0),
            "{k}: {v} vs {g}"
        );
    }
    // A second product on the SAME store exercises the swapped-in aux. Cloning
    // here would replace the workspace and let a broken reuse path pass.
    let twice_expected = got.multiply(&col_b);
    in_place.multiply_in_place(&col_b);
    assert_eq!(in_place.len(), twice_expected.len());
    for (k, want) in twice_expected.iter() {
        let actual = in_place
            .get(&k)
            .unwrap_or_else(|| panic!("missing key {k} after reused-aux product"));
        assert!(
            (actual - want).norm() <= 1e-12 * want.norm().max(1.0),
            "{k}: {actual} vs {want}"
        );
    }
}

/// The headline integration workload's shape (noisy TFIM Trotter propagation,
/// Heisenberg picture, `truncate()` after every operation) run on both backends
/// and compared term for term. This is the circuit the Phase-6 perf gate times,
/// so it is also the one whose numerics must be pinned: a backend that were
/// quietly dropping terms would look fast *and* wrong.
#[test]
fn tfim_trotter_propagation_is_identical() {
    const NQ: usize = 6;
    const STEPS: usize = 4;
    const THETA_X: f64 = 0.1;
    const THETA_ZZ: f64 = 0.0125;
    const NOISE: [f64; 3] = [2.5e-5, 2.5e-5, 2.5e-5];

    let seed: Vec<(&str, f64)> = vec![
        ("ZIIIII", 1.0),
        ("IZIIII", 1.0),
        ("IIZIII", 1.0),
        ("IIIZII", 1.0),
        ("IIIIZI", 1.0),
        ("IIIIIZ", 1.0),
    ];

    let (a, b) = both_backends!(
        n = NQ,
        policy = CombinedPolicy(
            CoefficientThreshold { threshold: 1e-6 },
            MaxPauliWeight(usize::MAX)
        ),
        seed = seed.clone(),
        |sum| {
            for _ in 0..STEPS {
                for i in 0..NQ {
                    sum.pauli_error(i, NOISE);
                    sum.truncate();
                    sum.rx(i, THETA_X);
                    sum.truncate();
                }
                for i in 0..NQ - 1 {
                    sum.pauli_error(i + 1, NOISE);
                    sum.truncate();
                    sum.pauli_error(i, NOISE);
                    sum.truncate();
                    sum.cnot(i, i + 1);
                    sum.rz(i + 1, THETA_ZZ);
                    sum.cnot(i, i + 1);
                    sum.truncate();
                }
            }
        }
    );
    assert!(
        a.len() > 50,
        "the Trotter run should hold a real support, got {}",
        a.len()
    );
    assert_identical(&a, &b);
}
