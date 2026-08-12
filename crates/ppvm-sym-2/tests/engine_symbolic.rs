// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The point of the crate, end to end: instantiate `ppvm-pauli-sum-2::Sum` over
//! the symbolic coefficient ring and propagate a parametric circuit through it,
//! exactly as old `PauliSum<config::fxhash::Byte<2, Term>>` could.
//!
//! `parametric_trace_golden_master` is a verbatim port of `examples/symbolic.rs`
//! (the `sym.trace.parametric` integration workload) against the **hard golden
//! master captured from the old crate on this tree**:
//!
//! ```text
//! $ cargo run --release --example symbolic
//! Trace expression: [1.000 * cos^3(%1) + -1.000 * sin^1(%0) sin^2(%1) cos^1(%1) + …]
//! Trace: 0.18803675917759355
//! ```
//!
//! The `Display` string is asserted **byte for byte**: symbolic output is
//! user-facing text, so a changed rendering is a behaviour change. It is stable
//! across runs because the monomial table's hasher is seed-free (integration
//! baseline, perf feature 4) and, on this workload, every monomial has a distinct
//! `(sin_pow, cos_pow)` sort key.

use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, PauliPattern, PauliWord, Sum};
use ppvm_sym_2::{Inner, Term};
use ppvm_traits_2::{Clifford, RotationOne, Trace};

type Key = PauliWord<[u8; 2]>;
type SymSum = Sum<HashMapStore<Key, Term>, NoPolicy>;

/// The old golden master, captured by running `examples/symbolic.rs` on this
/// tree (deterministic across runs).
const GOLDEN_DISPLAY: &str = "[1.000 * cos^3(%1) + -1.000 * sin^1(%0) sin^2(%1) cos^1(%1) + 1.000 * sin^3(%1) cos^1(%0) cos^1(%1) + 1.000 * sin^1(%0) sin^3(%1) cos^1(%0) + 1.000 * sin^1(%0) sin^3(%1) cos^1(%0) cos^1(%1)]";
const GOLDEN_VALUE: f64 = 0.188_036_759_177_593_55;

/// The `examples/symbolic.rs` circuit, run on the new engine over `Term`.
fn parametric_trace() -> Term {
    let mut sum: SymSum = SymSum::new(2);
    sum += (Key::from("ZZ"), Term::from(1.0));

    sum.rz(0, Term::var(0));
    sum.ry(0, Term::var(1));
    sum.rz(0, Term::var(0));

    sum.rz(1, Term::var(0));
    sum.ry(1, Term::var(1));
    sum.rz(1, Term::var(0));

    sum.cnot(0, 1);

    sum.rx(0, Term::var(1));
    sum.rx(1, Term::var(1));

    // Old's `PauliPattern::from("Z?*")`: every site is `I` or `Z`.
    sum.trace(&PauliPattern::zero_state())
}

#[test]
fn parametric_trace_golden_master() {
    let trace = parametric_trace();

    // (a) The exact old `Display` string (snapshot contract).
    assert_eq!(trace.to_string(), GOLDEN_DISPLAY);

    // (b) The numeric value. Tolerance 1e-12: the trace fold and the monomial
    // table are walked in a hash order that legitimately differs from old's.
    let value = trace.eval(&[1.1, 2.1]).unwrap();
    assert!(
        (value - GOLDEN_VALUE).abs() < 1e-12,
        "eval = {value}, golden = {GOLDEN_VALUE}"
    );

    // (c) The symbolic structure: the same five monomials with the same
    // coefficients, and `c0 == 0.0`.
    let Inner::Sum(s) = trace.inner() else {
        panic!("trace should be a map-backed sum, got {:?}", trace.inner());
    };
    assert_eq!(s.c0(), 0.0);
    assert_eq!(s.len(), 5);
    let mut got: Vec<(String, f64)> = s.iter().map(|(p, c)| (p.to_string(), c)).collect();
    got.sort_by(|a, b| a.0.cmp(&b.0));
    let mut want: Vec<(String, f64)> = vec![
        ("cos^3(%1)".into(), 1.0),
        ("sin^1(%0) sin^2(%1) cos^1(%1)".into(), -1.0),
        ("sin^3(%1) cos^1(%0) cos^1(%1)".into(), 1.0),
        ("sin^1(%0) sin^3(%1) cos^1(%0)".into(), 1.0),
        ("sin^1(%0) sin^3(%1) cos^1(%0) cos^1(%1)".into(), 1.0),
    ];
    want.sort_by(|a, b| a.0.cmp(&b.0));
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want.iter()) {
        assert_eq!(g.0, w.0);
        assert!((g.1 - w.1).abs() < 1e-12, "{}: {} vs {}", g.0, g.1, w.1);
    }
}

#[test]
fn trace_readout_resets_the_truncation_parameters() {
    // Behavioural contract 10: `std::iter::Sum for Term` folds from
    // `Term::from_f64(0.0)`, whose parameters are the defaults — so a user who
    // seeded `max_sin` on the propagated coefficients gets a trace accumulated
    // with NO `max_sin` bound.
    let mut c = Term::from(1.0);
    c.set_max_sin(1);
    c.set_min_eps(1e-12);

    let mut sum: SymSum = SymSum::new(2);
    sum += (Key::from("ZZ"), c);
    sum.rx(0, Term::var(0));
    sum.rx(1, Term::var(1));

    let trace = sum.trace(&PauliPattern::zero_state());
    assert_eq!(trace.max_sin(), usize::MAX);
    assert_eq!(trace.min_eps(), f64::EPSILON);
}

#[test]
fn seeded_max_sin_propagates_through_a_whole_circuit() {
    // Behavioural contract 1 + 2: the parameters are seeded on the INITIAL
    // observable coefficient and travel with it through every gate, because the
    // engine always writes `v.clone() * sin` / `*v *= cos` with `v` on the left.
    const K: usize = 2;
    let mut c = Term::from(1.0);
    c.set_max_sin(K);
    c.set_min_eps(1e-12);

    let mut sum: SymSum = SymSum::new(4);
    sum += (Key::from("ZIII"), c);
    for layer in 0..4u32 {
        for q in 0..4 {
            sum.rx(q, Term::var(2 * layer));
        }
        for q in 0..3 {
            // rzz(q, q+1, θ) = cnot; rz; cnot
            sum.cnot(q, q + 1);
            sum.rz(q + 1, Term::var(2 * layer + 1));
            sum.cnot(q, q + 1);
        }
    }

    // Golden master mined from OLD on this tree (`examples/symbolic.rs`'s crate,
    // same circuit, `PauliSum<config::fxhash::Byte<2, Term>>`):
    //
    //   support = 8, one_form = 1, sum_form = 7, monomials = 22,
    //   worst_monomial_sin_pow = 7
    //
    // Note `worst = 7 > K`: old's `max_sin` bound is consulted only inside
    // `Sum::add_term`/`Sum::mul_term`, i.e. only once a coefficient has been
    // promoted to the map-backed form. The `One × One → One` and
    // `Const × One → One` fast arms never truncate, so a coefficient that stays
    // a single monomial for the whole circuit escapes the bound entirely. That
    // is old behaviour and is preserved deliberately; the assertions below pin
    // both halves of it.
    let (mut one_form, mut sum_form, mut monomials, mut worst) = (0usize, 0usize, 0usize, 0usize);
    for (_, coeff) in sum.iter() {
        assert_eq!(coeff.max_sin(), K, "max_sin did not travel with the clone");
        match coeff.inner() {
            Inner::One(p, _) => {
                one_form += 1;
                monomials += 1;
                worst = worst.max(p.sin_pow());
            }
            Inner::Sum(s) => {
                sum_form += 1;
                monomials += s.len();
                for (p, _) in s.iter() {
                    assert!(
                        p.sin_pow() <= K,
                        "a map-backed monomial survived with sin_pow {} > {K}",
                        p.sin_pow()
                    );
                    worst = worst.max(p.sin_pow());
                }
            }
            other => panic!("unexpected coefficient form {other:?}"),
        }
    }
    assert_eq!(sum.len(), 8, "support size");
    assert_eq!(one_form, 1, "single-monomial coefficients");
    assert_eq!(sum_form, 7, "map-backed coefficients");
    assert_eq!(monomials, 22, "total monomials");
    assert_eq!(worst, 7, "the un-truncated `One`-form escapee");
}

#[test]
fn sum_level_truncation_is_inert_on_symbolic_coefficients() {
    // Behavioural contract 3: however tiny a symbolic coefficient's monomial
    // coefficients are, `CoefficientThreshold` never drops it — but a literal
    // constant coefficient is dropped.
    use ppvm_pauli_sum_2::CoefficientThreshold;

    type ThreshSum = Sum<HashMapStore<Key, Term>, CoefficientThreshold>;
    let policy = CoefficientThreshold { threshold: 1e-6 };
    let mut sum: ThreshSum = ThreshSum::with_policy(2, policy);
    sum += (Key::from("ZI"), Term::var(0).sin() * 1e-30);
    sum += (Key::from("IZ"), Term::from(1e-30));
    sum.truncate();

    assert!(
        sum.contains_key(&Key::from("ZI")),
        "symbolic term was dropped"
    );
    assert!(
        !sum.contains_key(&Key::from("IZ")),
        "constant term below threshold survived"
    );
}

// ---------------------------------------------------------------------------
// Per-key golden master, captured from OLD `ppvm-sym` on this tree.
// ---------------------------------------------------------------------------

/// `key => Display` for every surviving term of the circuit in
/// [`seeded_max_sin_propagates_through_a_whole_circuit`], captured by replaying
/// it on old `PauliSum<config::fxhash::Byte<2, Term>>`.
///
/// Compared **monomial-set-wise**, not byte-wise: old's `Sum` `Display` sorts by
/// `(sin_pow, cos_pow)`, a *non-total* order whose ties fall back to `FxHashMap`
/// iteration order. The new monomial layout (a packed sorted factor vector
/// instead of two `BTreeMap`s) and the `mul_term` aux double-buffer both change
/// bucket occupancy, so the tie order legitimately differs — see the crate-level
/// notes. Everything that is *not* tie order — the surviving key set, the
/// monomial set per key, each coefficient, the empty-`Sum` `[]` results of the
/// `mul_term` zero shortcut, and the single un-promoted `One` form — is asserted
/// exactly.
const OLD_GOLDEN: &[(&str, &str)] = &[
    (
        "XXXY",
        "1.000 * sin^1(%0) sin^1(%1) sin^1(%2) sin^1(%3) sin^1(%4) sin^1(%5) sin^1(%6)",
    ),
    ("XXXZ", "[]"),
    ("XXYI", "[]"),
    ("XXZI", "[]"),
    ("XYII", "[]"),
    (
        "XZII",
        "[1.000 * sin^1(%6) sin^1(%7) cos^1(%0) cos^1(%2) cos^1(%4) + 1.000 * sin^1(%4) sin^1(%7) cos^1(%0) cos^1(%2) cos^1(%5) cos^1(%6) + 1.000 * sin^1(%4) sin^1(%5) cos^1(%0) cos^1(%2) cos^1(%6) cos^1(%7) + 1.000 * sin^1(%2) sin^1(%3) cos^1(%0) cos^1(%4) cos^1(%5) cos^1(%6) cos^1(%7) + 1.000 * sin^1(%2) sin^1(%5) cos^1(%0) cos^1(%3) cos^1(%4) cos^1(%6) cos^1(%7) + 1.000 * sin^1(%2) sin^1(%7) cos^1(%0) cos^1(%3) cos^1(%4) cos^1(%5) cos^1(%6) + 1.000 * sin^1(%0) sin^1(%5) cos^1(%1) cos^1(%2) cos^1(%3) cos^1(%4) cos^1(%6) cos^1(%7) + 1.000 * sin^1(%0) sin^1(%1) cos^1(%2) cos^1(%3) cos^1(%4) cos^1(%5) cos^1(%6) cos^1(%7) + 1.000 * sin^1(%0) sin^1(%3) cos^1(%1) cos^1(%2) cos^1(%4) cos^1(%5) cos^1(%6) cos^1(%7) + 1.000 * sin^1(%0) sin^1(%7) cos^1(%1) cos^1(%2) cos^1(%3) cos^1(%4) cos^1(%5) cos^1(%6)]",
    ),
    (
        "YIII",
        "[1.000 * sin^1(%6) cos^1(%0) cos^1(%2) cos^1(%4) cos^1(%7) + 1.000 * sin^1(%4) cos^1(%0) cos^1(%2) cos^1(%5) cos^1(%6) cos^1(%7) + 1.000 * sin^1(%2) cos^1(%0) cos^1(%3) cos^1(%4) cos^1(%5) cos^1(%6) cos^1(%7) + 1.000 * sin^1(%0) cos^1(%1) cos^1(%2) cos^1(%3) cos^1(%4) cos^1(%5) cos^1(%6) cos^1(%7)]",
    ),
    (
        "ZIII",
        "[1.000 * cos^1(%0) cos^1(%2) cos^1(%4) cos^1(%6) + -1.000 * sin^1(%0) sin^1(%2) cos^1(%1) cos^1(%4) cos^1(%6) + -1.000 * sin^1(%4) sin^1(%6) cos^1(%0) cos^1(%2) cos^1(%5) + -1.000 * sin^1(%2) sin^1(%4) cos^1(%0) cos^1(%3) cos^1(%6) + -1.000 * sin^1(%0) sin^1(%4) cos^1(%1) cos^1(%2) cos^1(%3) cos^1(%6) + -1.000 * sin^1(%2) sin^1(%6) cos^1(%0) cos^1(%3) cos^1(%4) cos^1(%5) + -1.000 * sin^1(%0) sin^1(%6) cos^1(%1) cos^1(%2) cos^1(%3) cos^1(%4) cos^1(%5)]",
    ),
];

/// `key => eval(&[0.3 + 0.17·i])`, same replay, captured from old.
const OLD_GOLDEN_EVAL: &[(&str, f64)] = &[
    ("XXXY", 4.251_052_176_955_197e-2),
    ("XXXZ", -0.0),
    ("XXYI", 0.0),
    ("XXZI", -0.0),
    ("XYII", 0.0),
    ("XZII", 5.268_178_107_896_917e-1),
    ("YIII", 4.103_488_084_661_929e-2),
    ("ZIII", -3.975_817_799_003_077_4e-1),
];

/// Split a `Term` `Display` string into its sorted monomial list, so two
/// renderings that differ only in the (non-total) tie order compare equal.
fn monomial_set(s: &str) -> Vec<String> {
    let body = s
        .strip_prefix('[')
        .and_then(|b| b.strip_suffix(']'))
        .unwrap_or(s);
    if body.is_empty() {
        return Vec::new();
    }
    let mut v: Vec<String> = body.split(" + ").map(|m| m.trim().to_string()).collect();
    v.sort();
    v
}

/// Rebuild the `seeded_max_sin_propagates_through_a_whole_circuit` state.
fn seeded_trotter_state() -> SymSum {
    let mut c = Term::from(1.0);
    c.set_max_sin(2);
    c.set_min_eps(1e-12);

    let mut sum: SymSum = SymSum::new(4);
    sum += (Key::from("ZIII"), c);
    for layer in 0..4u32 {
        for q in 0..4 {
            sum.rx(q, Term::var(2 * layer));
        }
        for q in 0..3 {
            sum.cnot(q, q + 1);
            sum.rz(q + 1, Term::var(2 * layer + 1));
            sum.cnot(q, q + 1);
        }
    }
    sum
}

#[test]
fn trotter_replay_matches_the_old_golden_master() {
    let sum = seeded_trotter_state();

    let mut got: Vec<(String, String)> = sum
        .iter()
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect();
    got.sort();
    assert_eq!(got.len(), OLD_GOLDEN.len(), "support size: {got:?}");

    for ((gk, gv), (wk, wv)) in got.iter().zip(OLD_GOLDEN.iter()) {
        assert_eq!(gk, wk, "support key mismatch");
        assert_eq!(
            monomial_set(gv),
            monomial_set(wv),
            "monomial set differs for {gk}\n new: {gv}\n old: {wv}"
        );
        // The empty-`Sum` results of the `mul_term` zero shortcut must print as
        // `[]` on both sides (behavioural contract 9), and the un-promoted
        // `One` form must stay bracket-free.
        assert_eq!(
            gv.starts_with('['),
            wv.starts_with('['),
            "representation form differs for {gk}"
        );
    }
}

#[test]
fn trotter_replay_evaluates_to_the_old_values() {
    let sum = seeded_trotter_state();
    let vals: Vec<f64> = (0..8).map(|i| 0.3 + 0.17 * i as f64).collect();

    let mut got: Vec<(String, f64)> = sum
        .iter()
        .map(|(k, v)| (k.to_string(), v.eval(&vals).unwrap()))
        .collect();
    got.sort_by(|a, b| a.0.cmp(&b.0));
    assert_eq!(got.len(), OLD_GOLDEN_EVAL.len());
    for ((gk, gv), (wk, wv)) in got.iter().zip(OLD_GOLDEN_EVAL.iter()) {
        assert_eq!(gk, wk);
        assert!(
            (gv - wv).abs() < 1e-12,
            "{gk}: new {gv:e} vs old {wv:e} (Σ-fold order differs, but not by this much)"
        );
    }
}

#[test]
fn eval_error_semantics_match_old() {
    // Behavioural contract 7 / the `sym.expectation.grid` numeric bar: a `vals`
    // slice shorter than the highest variable id is an `Err` with old's message
    // shape — not a panic and not a silent zero.
    let sum = seeded_trotter_state();
    let short = [0.1, 0.2];
    let mut saw_err = false;
    for (_, v) in sum.iter() {
        if let Err(e) = v.eval(&short) {
            saw_err = true;
            let msg = e.to_string();
            assert!(
                msg.starts_with("variable %") && msg.ends_with(" not found"),
                "unexpected message: {msg}"
            );
        }
    }
    assert!(
        saw_err,
        "no coefficient referenced an out-of-range variable"
    );
}

#[test]
fn l4_operator_product_runs_over_the_symbolic_ring() {
    // Capability parity with old's `impl Mul<PauliSum> for PauliSum` (which
    // required `ComplexCoefficient`, and which `ppvm-sym::Term` implemented).
    // The `-2` spelling is `Sum::multiply`, bounded on `ImaginaryUnit`.
    //
    // `X · Y = i·Z`, so the product coefficient must DENOTE `i`. Old could not
    // express that: its `Prod` multiply dropped the phase and its `eval` returned
    // `f64`, so the same product evaluated as `+1` (`oldSuspectedBugs` #2/#4).
    // The oracle here is Lean — `lean/PPVM/Pauli/Phase.lean` `phaseExp_eq_ref`
    // and `lean/PPVM/Algebra/Twisted.lean` `iPow_add`.
    let mut x: SymSum = SymSum::new(2);
    x += (Key::from("XI"), Term::from(1.0));
    let mut y: SymSum = SymSum::new(2);
    y += (Key::from("YI"), Term::from(1.0));

    let p = x.multiply(&y);
    let c = p.get(&Key::from("ZI")).expect("X·Y should land on Z");
    let v = c.eval_complex(&[]).unwrap();
    assert!(v.re.abs() < 1e-15, "X·Y must be purely imaginary, got {v}");
    assert!((v.im - 1.0).abs() < 1e-15, "X·Y must be +i·Z, got {v}");

    // …and it stays symbolic: the same product with a symbolic coefficient
    // carries both the monomial and the phase.
    let mut xs: SymSum = SymSum::new(2);
    xs += (Key::from("XI"), Term::var(0).sin());
    let ps = xs.multiply(&y);
    let cs = ps.get(&Key::from("ZI")).unwrap();
    let vs = cs.eval_complex(&[0.7]).unwrap();
    assert!(vs.re.abs() < 1e-15);
    assert!((vs.im - 0.7f64.sin()).abs() < 1e-12, "{vs}");
}
