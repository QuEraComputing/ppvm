// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Phase-5 **differential** suite: new `ppvm-sym-2` vs old `ppvm-sym`.
//!
//! Structure:
//!
//! 1. the coefficient ring itself — construction, `add`/`sub`/`mul`/`neg`/
//!    `zero`, evaluation/substitution, `Display`, and the `Coefficient` trait
//!    obligations (`magnitude` vs old `cutoff`, `mul_sign`, `is_zero`);
//! 2. the twelve **behavioural contracts** from the integration baseline — the
//!    prime-directive parity checks, including *when* truncation fires;
//! 3. the deliberate **divergences** (old bugs the Lean oracle adjudicates),
//!    pinned explicitly so they can never drift silently; and
//! 4. the six `sym.*` **integration workloads**, propagated end to end through
//!    BOTH `ppvm-pauli-sum-2::Sum<…, Term2>` and old
//!    `PauliSum<config::fxhash::Byte<8, Term>>` under an identical gate
//!    sequence, identical storage width and identical (coefficient-intrinsic)
//!    truncation policy.
//!
//! Old's `Sum`/`Prod` fields are `pub(crate)`, so the structural comparison goes
//! through `Display` (parsed by [`TermView`]) plus `eval` — see the
//! `ppvm_conformance_2::sym` module docs for why that is a *stronger* contract
//! than poking at private fields, not a weaker one.

use std::collections::{BTreeMap, BTreeSet};

use ppvm_conformance_2::sym::*;
use ppvm_conformance_2::{GateOp, random_circuit, seeded_rng};
use rand::RngExt;
use rand::rngs::StdRng;

// Old trait surface.
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_pauli_word::pattern::PauliPattern as OldPauliPattern;
use ppvm_sym::{Prod as OldProd, Sum as OldSumTy};
use ppvm_traits::traits::{
    Clifford as OldClifford, Coefficient as OldCoefficient,
    ComplexCoefficient as OldComplexCoefficient, RotationOne as OldRotationOne,
};

// New trait surface.
use ppvm_sym_2::{Inner, Prod as NewProd, Sum as NewSumTy};
use ppvm_traits_2::{
    Angle as NewAngle, Clifford as NewClifford, Coefficient as NewCoefficient,
    RotationOne as NewRotationOne,
};

use num::Zero;

/// Structural + numeric equality of an old and a new `Term`.
///
/// * the parsed `Display` view (representation form, `c0`, and the
///   monomial → printed-coefficient map) — **exact**;
/// * `eval` at every supplied angle vector — within `1e-12`.
#[track_caller]
fn assert_terms_match(old: &OldTerm, new: &NewTerm, angles: &[Vec<f64>], what: &str) {
    let ov = old_view(old);
    let nv = new_view(new);
    assert_eq!(
        ov.form, nv.form,
        "{what}: representation form differs\n old: {old}\n new: {new}"
    );
    assert_eq!(
        ov.c0, nv.c0,
        "{what}: constant part differs\n old: {old}\n new: {new}"
    );
    assert_eq!(
        ov.monomials, nv.monomials,
        "{what}: monomial set differs\n old: {old}\n new: {new}"
    );
    for vals in angles {
        match (old.eval(vals), new.eval(vals)) {
            (Ok(o), Ok(n)) => assert!(
                (o - n).abs() < 1e-12,
                "{what}: eval differs at {vals:?}: old {o} vs new {n}"
            ),
            (Err(o), Err(n)) => assert_eq!(
                o.to_string(),
                n.to_string(),
                "{what}: eval error message differs"
            ),
            (o, n) => panic!("{what}: eval Ok/Err disagree: old {o:?} vs new {n:?}"),
        }
    }
}

// ===========================================================================
// 1. The coefficient ring: construction, the ring surface, eval, Display.
// ===========================================================================

/// A replayable symbolic expression, interpreted once per crate so both see
/// literally the same construction sequence.
///
/// `+ f64` / `- f64` are deliberately **absent**: they route through
/// `AddAssign<f64> for Term`, whose `One` receiver arm is old bug #1 (a silent
/// no-op) and is a documented divergence, pinned separately by
/// [`divergence_add_f64_to_a_single_monomial`].
#[derive(Debug, Clone)]
enum Expr {
    Sin(u32),
    Cos(u32),
    Const(f64),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Neg(Box<Expr>),
    ScaleF(Box<Expr>, f64),
}

fn random_expr(rng: &mut StdRng, n_vars: u32, depth: usize) -> Expr {
    if depth == 0 {
        return match rng.random_range(0..3usize) {
            0 => Expr::Sin(rng.random_range(0..n_vars)),
            1 => Expr::Cos(rng.random_range(0..n_vars)),
            _ => Expr::Const(rng.random_range(-2.0..2.0)),
        };
    }
    match rng.random_range(0..7usize) {
        0 => Expr::Sin(rng.random_range(0..n_vars)),
        1 => Expr::Cos(rng.random_range(0..n_vars)),
        2 => Expr::Add(
            Box::new(random_expr(rng, n_vars, depth - 1)),
            Box::new(random_expr(rng, n_vars, depth - 1)),
        ),
        3 => Expr::Sub(
            Box::new(random_expr(rng, n_vars, depth - 1)),
            Box::new(random_expr(rng, n_vars, depth - 1)),
        ),
        4 => Expr::Mul(
            Box::new(random_expr(rng, n_vars, depth - 1)),
            Box::new(random_expr(rng, n_vars, depth - 1)),
        ),
        5 => Expr::Neg(Box::new(random_expr(rng, n_vars, depth - 1))),
        _ => Expr::ScaleF(
            Box::new(random_expr(rng, n_vars, depth - 1)),
            rng.random_range(-2.0..2.0),
        ),
    }
}

fn eval_old(e: &Expr, max_sin: usize, min_eps: f64) -> OldTerm {
    match e {
        Expr::Sin(u) => OldTerm::var(*u).sin(),
        Expr::Cos(u) => OldTerm::var(*u).cos(),
        Expr::Const(c) => {
            let mut t = OldTerm::from(*c);
            t.set_max_sin(max_sin);
            t.set_min_eps(min_eps);
            t
        }
        Expr::Add(a, b) => eval_old(a, max_sin, min_eps) + eval_old(b, max_sin, min_eps),
        Expr::Sub(a, b) => eval_old(a, max_sin, min_eps) - eval_old(b, max_sin, min_eps),
        Expr::Mul(a, b) => eval_old(a, max_sin, min_eps) * eval_old(b, max_sin, min_eps),
        Expr::Neg(a) => -eval_old(a, max_sin, min_eps),
        Expr::ScaleF(a, c) => eval_old(a, max_sin, min_eps) * *c,
    }
}

fn eval_new(e: &Expr, max_sin: usize, min_eps: f64) -> NewTerm {
    match e {
        Expr::Sin(u) => NewTerm::var(*u).sin(),
        Expr::Cos(u) => NewTerm::var(*u).cos(),
        Expr::Const(c) => {
            let mut t = NewTerm::from(*c);
            t.set_max_sin(max_sin);
            t.set_min_eps(min_eps);
            t
        }
        Expr::Add(a, b) => eval_new(a, max_sin, min_eps) + eval_new(b, max_sin, min_eps),
        Expr::Sub(a, b) => eval_new(a, max_sin, min_eps) - eval_new(b, max_sin, min_eps),
        Expr::Mul(a, b) => eval_new(a, max_sin, min_eps) * eval_new(b, max_sin, min_eps),
        Expr::Neg(a) => -eval_new(a, max_sin, min_eps),
        Expr::ScaleF(a, c) => eval_new(a, max_sin, min_eps) * *c,
    }
}

#[test]
fn construction_matches_old() {
    let angles = vec![vec![0.3, 1.1, -2.2, 0.0], vec![1.7, -0.4, 2.9, 3.0]];
    // Bare variables, the sin/cos wrappers, and the four `From` conversions.
    assert_eq!(OldTerm::var(3).to_string(), NewTerm::var(3).to_string());
    assert_terms_match(
        &OldTerm::var(1).sin(),
        &NewTerm::var(1).sin(),
        &angles,
        "sin(x1)",
    );
    assert_terms_match(
        &OldTerm::var(2).cos(),
        &NewTerm::var(2).cos(),
        &angles,
        "cos(x2)",
    );
    // Constant folding on a constant argument.
    assert_terms_match(
        &OldTerm::from(0.5).sin(),
        &NewTerm::from(0.5).sin(),
        &angles,
        "sin(0.5)",
    );
    assert_terms_match(
        &OldTerm::from(0.5).cos(),
        &NewTerm::from(0.5).cos(),
        &angles,
        "cos(0.5)",
    );
    for (o, n) in [
        (OldTerm::from(2.0f64), NewTerm::from(2.0f64)),
        (OldTerm::from(2.0f32), NewTerm::from(2.0f32)),
        (OldTerm::from(2i32), NewTerm::from(2i32)),
        (OldTerm::from(2i64), NewTerm::from(2i64)),
        (OldTerm::from_f64(-3.25), NewTerm::from_f64(-3.25)),
    ] {
        assert_terms_match(&o, &n, &angles, "From conversions");
    }
    // Every conversion lands on the SAME value, on both crates.
    assert_eq!(OldTerm::from(2i32), OldTerm::from(2.0f64));
    assert_eq!(NewTerm::from(2i32), NewTerm::from(2.0f64));
    assert_eq!(NewTerm::from(2i64), NewTerm::from(2.0f64));
    assert_eq!(NewTerm::from(2.0f32), NewTerm::from(2.0f64));
}

#[test]
fn ring_surface_matches_old_on_seeded_random_expressions() {
    // add / sub / mul / neg / scalar-scale over sin/cos/const atoms, replayed on
    // both crates, at the default truncation parameters and at a seeded
    // `max_sin`/`min_eps` (so the drop-at-accumulate path is exercised too).
    let mut rng = seeded_rng(0x5111_0DE7);
    let n_vars = 4u32;
    let angle_rng = &mut seeded_rng(0xA9_1E5);
    let angles = angle_grid(angle_rng, n_vars as usize, 8);

    for case in 0..256 {
        let e = random_expr(&mut rng, n_vars, 3);
        for (max_sin, min_eps) in [(usize::MAX, f64::EPSILON), (2, 1e-12), (1, 1e-6)] {
            let o = eval_old(&e, max_sin, min_eps);
            let n = eval_new(&e, max_sin, min_eps);
            assert_terms_match(
                &o,
                &n,
                &angles,
                &format!("case {case} (max_sin={max_sin}, min_eps={min_eps:e}): {e:?}"),
            );
        }
    }
}

#[test]
fn zero_neg_and_sub_match_old() {
    let angles = vec![vec![0.4, -1.2]];
    let o: OldTerm = num::Zero::zero();
    let n: NewTerm = num::Zero::zero();
    assert_terms_match(&o, &n, &angles, "zero()");
    assert!(o.is_zero() && n.is_zero());

    assert_terms_match(
        &(-OldTerm::from(2.0)),
        &(-NewTerm::from(2.0)),
        &angles,
        "neg const",
    );
    assert_eq!(-OldTerm::from(2.0), OldTerm::from(-2.0));
    assert_eq!(-NewTerm::from(2.0), NewTerm::from(-2.0));

    let o = OldTerm::var(0).sin() - OldTerm::var(1).cos();
    let n = NewTerm::var(0).sin() - NewTerm::var(1).cos();
    assert_terms_match(&o, &n, &angles, "sub");

    // `num::One` is a NEW addition (old had no impl); it must agree with
    // `From<f64>(1.0)`, the only spelling consistent with `zero()`.
    assert_eq!(<NewTerm as num::One>::one(), NewTerm::from(1.0));
}

#[test]
fn display_snapshots_match_old() {
    // Old's own two `display::tests` cases, replayed on both crates byte for
    // byte (behavioural contract 8).
    let mut op = OldProd::new();
    op.mul_sin(1);
    op.mul_sin(1);
    op.mul_cos(2);
    let mut np = NewProd::new();
    np.mul_sin(1);
    np.mul_sin(1);
    np.mul_cos(2);
    assert_eq!(op.to_string(), "sin^2(%1) cos^1(%2)");
    assert_eq!(np.to_string(), op.to_string());

    let mut os = OldTerm::from_f64(3.0);
    os += OldTerm::var(1).sin();
    os += OldTerm::var(2).cos();
    let mut ns = NewTerm::from_f64(3.0);
    ns += NewTerm::var(1).sin();
    ns += NewTerm::var(2).cos();
    assert_eq!(
        os.to_string(),
        "[3.000 + 1.000 * cos^1(%2) + 1.000 * sin^1(%1)]"
    );
    assert_eq!(ns.to_string(), os.to_string());

    os.set_max_sin(2);
    ns.set_max_sin(2);
    os *= OldTerm::var(2).sin();
    ns *= NewTerm::var(2).sin();
    os *= OldTerm::var(1).sin();
    ns *= NewTerm::var(1).sin();
    assert_eq!(
        os.to_string(),
        "[3.000 * sin^1(%1) sin^1(%2) + 1.000 * sin^1(%1) sin^1(%2) cos^1(%2)]"
    );
    assert_eq!(ns.to_string(), os.to_string());

    // The four small forms.
    assert_eq!(OldTerm::var(7).to_string(), "%7");
    assert_eq!(NewTerm::var(7).to_string(), "%7");
    // `Const` prints with NO precision spec.
    assert_eq!(OldTerm::from_f64(1.5).to_string(), "1.5");
    assert_eq!(NewTerm::from_f64(1.5).to_string(), "1.5");
    assert_eq!(OldTerm::var(0).sin().to_string(), "1.000 * sin^1(%0)");
    assert_eq!(NewTerm::var(0).sin().to_string(), "1.000 * sin^1(%0)");

    // Old's unbalanced rendering of a `Sum` with a non-zero `c0` and an empty
    // table: it returns *before* the closing bracket. Reproduced exactly.
    let o = OldTerm::from(1.0) + OldTerm::var(0).sin() * 1e-300;
    let n = NewTerm::from(1.0) + NewTerm::var(0).sin() * 1e-300;
    assert_eq!(o.to_string(), "[1.000 ");
    assert_eq!(n.to_string(), o.to_string());
}

#[test]
fn display_tie_break_is_reported_not_silently_changed() {
    // Behavioural contract 8: the `Sum` sort key `(sin_pow, cos_pow)` is a
    // NON-total order; ties fall back to monomial-table iteration order, which
    // is a function of the (seed-free) hash values and the table's capacity
    // history. Both demanded perf features change that — the packed-vector
    // monomial layout changes every digest and the `mul_term` aux double-buffer
    // changes the capacity history — so the tie order legitimately differs.
    //
    // What must NOT differ is anything the ordering is *about*. This test builds
    // two monomials that TIE on `(sin_pow, cos_pow)` and asserts the monomial
    // SETS agree exactly while recording whether the byte strings do.
    let mut o = OldTerm::from(0.0);
    o += OldTerm::var(0).sin();
    o += OldTerm::var(1).sin();
    let mut n = NewTerm::from(0.0);
    n += NewTerm::var(0).sin();
    n += NewTerm::var(1).sin();

    let ov = old_view(&o);
    let nv = new_view(&n);
    assert_eq!(
        ov.monomials, nv.monomials,
        "tied monomials must still be the same SET"
    );
    assert_eq!(ov.form, nv.form);
    // On this input the tie order happens to agree as well; assert it so a
    // future layout change that *does* reorder shows up as a reported
    // divergence rather than passing unnoticed.
    assert_eq!(
        o.to_string(),
        n.to_string(),
        "tie-break order diverged — a Display divergence, report it"
    );
}

#[test]
fn eval_fold_order_is_bit_identical_to_old() {
    // Behavioural contract 7: within a monomial the sin factors are folded in
    // ascending variable order, then the cos factors — any reassociation moves
    // the last ulp, so this is a BIT comparison.
    let mut op = OldProd::sin(0);
    op.mul_sin(0);
    op.mul_sin(0);
    op.mul_sin(1);
    op.mul_sin(1);
    op.mul_cos(0);
    let mut np = NewProd::sin(0);
    np.mul_sin(0);
    np.mul_sin(0);
    np.mul_sin(1);
    np.mul_sin(1);
    np.mul_cos(0);

    let vals = [1.1, 2.1];
    let o = op.eval(&vals).unwrap();
    let n = np.eval(&vals).unwrap();
    assert_eq!(o.to_bits(), n.to_bits(), "fold order changed: {o} vs {n}");

    // The empty product is `1.0` on both.
    assert_eq!(OldProd::new().eval(&[]).unwrap(), 1.0);
    assert_eq!(NewProd::new().eval(&[]).unwrap(), 1.0);
}

#[test]
fn eval_error_shape_matches_old() {
    let o = OldTerm::var(5).sin().eval(&[0.1]).unwrap_err();
    let n = NewTerm::var(5).sin().eval(&[0.1]).unwrap_err();
    assert_eq!(o.to_string(), "variable %5 not found");
    assert_eq!(n.to_string(), o.to_string());

    // A bare variable evaluates to its binding on both.
    assert_eq!(OldTerm::var(3).eval(&[0., 1., 2., 9.]).unwrap(), 9.0);
    assert_eq!(NewTerm::var(3).eval(&[0., 1., 2., 9.]).unwrap(), 9.0);
    // …and errors identically when unbound.
    assert_eq!(
        OldTerm::var(3).eval(&[0.0]).unwrap_err().to_string(),
        NewTerm::var(3).eval(&[0.0]).unwrap_err().to_string()
    );

    // Constant folding needs no variables at all.
    assert_eq!(
        OldTerm::from_f64(0.5).sin().eval(&[]).unwrap(),
        NewTerm::from_f64(0.5).sin().eval(&[]).unwrap()
    );
}

// ===========================================================================
// 2. `Coefficient` trait obligations: magnitude/cutoff, mul_sign, is_zero.
// ===========================================================================

#[test]
fn magnitude_reproduces_old_cutoff_exactly() {
    // Behavioural contract 3. Old: `cutoff(t)` is `|c| < t` for `Const` and
    // `false` for every symbolic form. New: the keep-rule is
    // `magnitude() >= threshold`, so `magnitude()` must be `|c|` for `Const`
    // and something that always passes for the rest.
    let cases: Vec<(OldTerm, NewTerm, &str)> = vec![
        (OldTerm::from(1e-30), NewTerm::from(1e-30), "tiny const"),
        (OldTerm::from(0.0), NewTerm::from(0.0), "zero const"),
        (OldTerm::from(3.5), NewTerm::from(3.5), "big const"),
        (
            OldTerm::var(0).sin() * 1e-30,
            NewTerm::var(0).sin() * 1e-30,
            "tiny One",
        ),
        (
            OldTerm::from(1e-30) + OldTerm::var(0).sin() * 1e-30,
            NewTerm::from(1e-30) + NewTerm::var(0).sin() * 1e-30,
            "tiny Sum",
        ),
    ];
    for threshold in [1e-6, 1e-12, 1.0, 10.0] {
        for (o, n, what) in &cases {
            let old_drops = OldCoefficient::cutoff(o, threshold);
            let new_drops = NewCoefficient::magnitude(n) < threshold;
            assert_eq!(
                old_drops, new_drops,
                "{what} @ {threshold:e}: old drops {old_drops}, new drops {new_drops}"
            );
        }
    }
    // …spelled out: symbolic forms are NEVER droppable.
    assert_eq!(
        NewCoefficient::magnitude(&(NewTerm::var(0).sin())),
        f64::INFINITY
    );
}

#[test]
fn mul_sign_matches_old() {
    let angles = vec![vec![0.4, 1.2]];
    for sign in [1i8, -1] {
        for (o, n, what) in [
            (OldTerm::from(2.5), NewTerm::from(2.5), "const"),
            (
                OldTerm::var(0).sin(),
                NewTerm::var(0).sin(),
                "single monomial",
            ),
            (
                OldTerm::var(0).sin() + OldTerm::var(1).cos(),
                NewTerm::var(0).sin() + NewTerm::var(1).cos(),
                "map-backed sum",
            ),
        ] {
            assert_terms_match(
                &OldCoefficient::mul_sign(&o, sign),
                &NewCoefficient::mul_sign(&n, sign),
                &angles,
                &format!("mul_sign({sign}) on {what}"),
            );
        }
    }
}

#[test]
fn is_zero_matches_old_on_every_form() {
    // Behavioural contract 4: `is_zero` is true ONLY for a constant under
    // `min_eps`; every non-`Const` form is false, including an empty `Sum`
    // (which denotes 0) and a `Sum` whose coefficients all cancelled to 0.
    assert!(OldTerm::from(0.0).is_zero() && NewTerm::from(0.0).is_zero());
    assert!(!OldTerm::from(1.0).is_zero() && !NewTerm::from(1.0).is_zero());

    // Exact cancellation leaves a zero-coefficient monomial in the table.
    let mut o = OldTerm::var(0).sin();
    o += OldTerm::var(0).sin() * -1.0;
    let mut n = NewTerm::var(0).sin();
    n += NewTerm::var(0).sin() * -1.0;
    assert!(!o.is_zero(), "old: cancelled sum must not be is_zero");
    assert!(!n.is_zero(), "new: cancelled sum must not be is_zero");
    assert_ne!(o, OldTerm::from(0.0));
    assert_ne!(n, NewTerm::from(0.0));
    assert_eq!(old_view(&o).n_monomials(), 1);
    assert_eq!(new_view(&n).n_monomials(), 1);
    assert!(o.to_string().contains("0.000"));
    assert_eq!(n.to_string(), o.to_string());

    // The empty `Sum` left by `mul_term`'s zero shortcut (contract 9).
    let (o, n) = zero_shortcut_pair();
    assert!(!o.is_zero() && !n.is_zero());
}

// ===========================================================================
// 3. Behavioural contracts 1..12.
// ===========================================================================

/// Contract 9's fixture: a map-backed `Term` with `max_sin = 1` multiplied by
/// `sin(x0)·sin(x1)`, which trips `mul_term`'s whole-sum-to-zero shortcut.
fn zero_shortcut_pair() -> (OldTerm, NewTerm) {
    let mut o = OldTerm::from(1.0) + OldTerm::var(1).cos();
    o.set_max_sin(1);
    o *= OldTerm::var(0).sin() * OldTerm::var(1).sin();
    let mut n = NewTerm::from(1.0) + NewTerm::var(1).cos();
    n.set_max_sin(1);
    n *= NewTerm::var(0).sin() * NewTerm::var(1).sin();
    (o, n)
}

#[test]
fn contract_1_truncation_parameters_are_inherited_from_the_lhs_only() {
    // A map-backed receiver is needed to observe it: the `Const × One` and
    // `One × One` fast arms never consult `max_sin` at all, on either crate.
    let mut oa = OldTerm::from(1.0) + OldTerm::var(1).cos();
    oa.set_max_sin(0);
    let ob = OldTerm::var(0).sin(); // max_sin = usize::MAX
    let mut na = NewTerm::from(1.0) + NewTerm::var(1).cos();
    na.set_max_sin(0);
    let nb = NewTerm::var(0).sin();

    let o_ab = oa.clone() * ob.clone();
    let n_ab = na.clone() * nb.clone();
    let o_ba = ob * oa;
    let n_ba = nb * na;

    assert_eq!(old_view(&o_ab).n_monomials(), 0, "old: LHS max_sin=0 wins");
    assert_eq!(new_view(&n_ab).n_monomials(), 0, "new: LHS max_sin=0 wins");
    assert!(
        old_view(&o_ba).n_monomials() > 0,
        "old: LHS max_sin=MAX keeps"
    );
    assert!(
        new_view(&n_ba).n_monomials() > 0,
        "new: LHS max_sin=MAX keeps"
    );
    assert_eq!(o_ab.to_string(), n_ab.to_string());
    assert_eq!(o_ba.to_string(), n_ba.to_string());
}

#[test]
fn contract_2_truncation_is_always_on_inside_the_coefficient() {
    // No `truncate()`/`reduce()` call anywhere — the monomials are already gone.
    let (o, n) = zero_shortcut_pair();
    assert_eq!(old_view(&o).n_monomials(), 0);
    assert_eq!(new_view(&n).n_monomials(), 0);
    assert_eq!(o.to_string(), "[]");
    assert_eq!(n.to_string(), "[]");

    // Conversely, the SUM-level `truncate()` must not touch a symbolic
    // coefficient's monomials.
    let mut os = new_old_sum(2);
    os += ("ZI", OldTerm::var(0).sin() + OldTerm::var(1).cos());
    let mut ns = new_new_sum(2);
    ns += (
        NewSymKey::from("ZI"),
        NewTerm::var(0).sin() + NewTerm::var(1).cos(),
    );
    let before_o = old_sym_support(&os);
    let before_n = new_sym_support(&ns);
    os.truncate();
    ns.truncate();
    assert_eq!(old_sym_support(&os).len(), before_o.len());
    assert_eq!(new_sym_support(&ns).len(), before_n.len());
    assert_eq!(
        old_view(&old_sym_support(&os)[0].1),
        old_view(&before_o[0].1)
    );
    assert_eq!(
        new_view(&new_sym_support(&ns)[0].1),
        new_view(&before_n[0].1)
    );
}

#[test]
fn contract_2_min_eps_drops_at_insert_and_is_not_a_post_pass() {
    // The `min_eps` half of contract 2, which the integration `min_eps` sweep
    // CANNOT reach: every coefficient the Trotter workload produces has
    // magnitude >= 1 (each rotation branch multiplies by a `sin`/`cos` whose
    // `f64` weight is exactly 1), so no `min_eps` below 1.0 ever fires there and
    // the sweep is vacuous by construction. This test drives the branch head-on.
    //
    // `max_sin` is a monomial ideal, so dropping at insert equals truncating the
    // finished product (`mulMono_drop_at_insert_eq_drop_at_end`,
    // `lean/PPVM/Instantiations/Symbolic.lean`). `min_eps` reads the COEFFICIENT
    // and is NOT interchangeable with a post-pass — `eps_drop_at_insert_ne_drop_at_end`:
    // two sub-threshold contributions to the same monomial are each dropped on
    // arrival, so the monomial vanishes even though their sum is above the
    // threshold and a post-pass `retain` would have kept it. That distinction is
    // exactly the branch a redesign is most likely to relocate, so it is pinned
    // on both crates.
    let build_old = |eps: f64, contrib: f64, times: usize| {
        let mut t = OldTerm::from(1.0) + OldTerm::var(1).cos();
        t.set_min_eps(eps);
        for _ in 0..times {
            t += OldTerm::var(0).sin() * contrib;
        }
        t
    };
    let build_new = |eps: f64, contrib: f64, times: usize| {
        let mut t = NewTerm::from(1.0) + NewTerm::var(1).cos();
        t.set_min_eps(eps);
        for _ in 0..times {
            t += NewTerm::var(0).sin() * contrib;
        }
        t
    };

    // Two contributions of 6e-4 sum to 1.2e-3 > 1e-3: both are dropped ON
    // ARRIVAL, so `sin^1(%0)` never enters the table at all.
    let o = build_old(1e-3, 6e-4, 2);
    let n = build_new(1e-3, 6e-4, 2);
    let ov = old_view(&o);
    assert_eq!(ov, new_view(&n), "min_eps drop-at-insert diverged");
    assert!(
        !ov.monomials.contains_key("sin^1(%0)"),
        "a post-pass would have kept the accumulated 1.2e-3; old drops at insert: {ov:?}"
    );
    assert_eq!(ov.n_monomials(), 1, "only cos^1(%1) survives");

    // Positive control: one contribution ABOVE the threshold does survive, so
    // the assertion above is about the threshold and not about the fixture.
    let o = build_old(1e-3, 2e-3, 1);
    let n = build_new(1e-3, 2e-3, 1);
    let ov = old_view(&o);
    assert_eq!(ov, new_view(&n));
    assert!(ov.monomials.contains_key("sin^1(%0)"));

    // …and the accumulate-then-fall-below case is NOT re-checked: `add_term`
    // does `*entry += coeff` and never re-tests the running total, so a monomial
    // that cancels to 0.0 stays with coefficient 0.0 (contract 4).
    let o = build_old(1e-3, 2e-3, 1) + OldTerm::var(0).sin() * -2e-3;
    let n = build_new(1e-3, 2e-3, 1) + NewTerm::var(0).sin() * -2e-3;
    let ov = old_view(&o);
    assert_eq!(ov, new_view(&n));
    assert_eq!(
        ov.monomials.get("sin^1(%0)").map(String::as_str),
        Some("0.000"),
        "a cancelled monomial stays in the table with coefficient 0.0"
    );

    // The `min_eps` arm of `mul_term`'s whole-sum-to-zero shortcut (perf
    // feature 8 / contract 9): a sub-threshold MULTIPLIER clears the whole
    // table in one `clear()` and leaves an EMPTY `Sum`, not a `Const(0.0)`.
    let mut o = OldTerm::from(1.0) + OldTerm::var(1).cos();
    o.set_min_eps(1e-3);
    o *= OldTerm::var(0).sin() * 1e-9;
    let mut n = NewTerm::from(1.0) + NewTerm::var(1).cos();
    n.set_min_eps(1e-3);
    n *= NewTerm::var(0).sin() * 1e-9;
    assert_eq!(o.to_string(), "[]");
    assert_eq!(n.to_string(), "[]");
    assert!(!o.is_zero() && !n.is_zero());
    assert_ne!(o, OldTerm::from(0.0));
    assert_ne!(n, NewTerm::from(0.0));
}

#[test]
fn contract_3_sum_level_threshold_is_inert_on_symbolic_coefficients() {
    use ppvm_pauli_sum::config::fxhash::Byte as OldByte;
    use ppvm_pauli_sum::strategy::CoefficientThreshold as OldThreshold;
    use ppvm_pauli_sum_2::{CoefficientThreshold as NewThreshold, HashMapStore, Sum as NewSum};

    type OldThreshSum = OldPauliSum<OldByte<8, OldTerm, OldThreshold>>;
    type NewThreshSum = NewSum<HashMapStore<NewSymKey, NewTerm>, NewThreshold>;

    let mut o: OldThreshSum = OldPauliSum::builder()
        .n_qubits(2)
        .strategy(OldThreshold(1e-6))
        .build();
    o += ("ZI", OldTerm::var(0).sin() * 1e-30);
    o += ("IZ", OldTerm::from(1e-30));
    o.truncate();

    let mut n: NewThreshSum = NewThreshSum::with_policy(2, NewThreshold { threshold: 1e-6 });
    n += (NewSymKey::from("ZI"), NewTerm::var(0).sin() * 1e-30);
    n += (NewSymKey::from("IZ"), NewTerm::from(1e-30));
    n.truncate();

    let old_keys: BTreeSet<String> = o.data().keys().map(|k| k.to_string()).collect();
    let new_keys: BTreeSet<String> = n.iter().map(|(k, _)| k.to_string()).collect();
    assert_eq!(old_keys, new_keys, "old {old_keys:?} vs new {new_keys:?}");
    assert!(old_keys.contains("ZI"), "the symbolic term must survive");
    assert!(
        !old_keys.contains("IZ"),
        "the constant term must be dropped"
    );
}

#[test]
fn contract_5_partial_eq_is_representational() {
    // All three of these DENOTE 1; none of them are equal, on either crate.
    let o_const = OldTerm::from(1.0);
    let o_one = OldComplexCoefficient::mul_phase(&OldTerm::from(1.0), 0);
    let o_sum = OldTerm::from(1.0) + OldTerm::var(0).sin() * 1e-300;
    let n_const = NewTerm::from(1.0);
    let n_one = NewTerm::from(1.0).mul_phase(0);
    let n_sum = NewTerm::from(1.0) + NewTerm::var(0).sin() * 1e-300;

    // The three forms render distinguishably, identically on both crates.
    assert_eq!(o_const.to_string(), n_const.to_string());
    assert_eq!(o_one.to_string(), n_one.to_string());
    assert_eq!(o_sum.to_string(), n_sum.to_string());
    assert!(matches!(n_const.inner(), Inner::Const(_)));
    assert!(matches!(n_one.inner(), Inner::One(..)));
    assert!(matches!(n_sum.inner(), Inner::Sum(_)));

    assert!(o_const != o_one && o_one != o_sum && o_const != o_sum);
    assert!(n_const != n_one && n_one != n_sum && n_const != n_sum);

    // `Prod`'s `Eq`/`Hash` include the phase byte, so `sin(x)` and `i·sin(x)`
    // are distinct monomials that never coalesce — on both crates.
    let op = OldProd::sin(0);
    let mut oq = OldProd::sin(0);
    oq.add_phase(1);
    assert_ne!(op, oq);
    assert_ne!(fx_hash(&op), fx_hash(&oq));

    let np = NewProd::sin(0);
    let mut nq = NewProd::sin(0);
    nq.add_phase(1);
    assert_ne!(np, nq);
    assert_ne!(fx_hash(&np), fx_hash(&nq));
}

fn fx_hash<T: std::hash::Hash>(v: &T) -> u64 {
    use std::hash::Hasher;
    let mut h = fxhash::FxHasher64::default();
    v.hash(&mut h);
    h.finish()
}

#[test]
fn contract_6_sin_cos_are_partial_and_constant_folding() {
    // Positive halves already covered by `construction_matches_old`; the panic
    // halves are asserted per-crate below (a panicking closure per crate, so a
    // single test can cover both sides).
    for (what, old_call, new_call) in [
        (
            "sin of a compound",
            (|| {
                let _ = OldTerm::var(0).sin().sin();
            }) as fn(),
            (|| {
                let _ = NewTerm::var(0).sin().sin();
            }) as fn(),
        ),
        (
            "cos of a compound",
            (|| {
                let _ = OldTerm::var(0).cos().cos();
            }) as fn(),
            (|| {
                let _ = NewTerm::var(0).cos().cos();
            }) as fn(),
        ),
        (
            "arithmetic on a bare variable",
            (|| {
                let _ = OldTerm::var(0) * OldTerm::from(2.0);
            }) as fn(),
            (|| {
                let _ = NewTerm::var(0) * NewTerm::from(2.0);
            }) as fn(),
        ),
        (
            "mul_sign on a bare variable",
            (|| {
                let _ = OldCoefficient::mul_sign(&OldTerm::var(0), 1);
            }) as fn(),
            (|| {
                let _ = NewCoefficient::mul_sign(&NewTerm::var(0), 1);
            }) as fn(),
        ),
        (
            "mul_phase on a bare variable",
            (|| {
                let _ = OldComplexCoefficient::mul_phase(&OldTerm::var(0), 1);
            }) as fn(),
            (|| {
                let _ = NewTerm::var(0).mul_phase(1);
            }) as fn(),
        ),
        (
            "sin_cos of a compound angle",
            (|| {
                let _ = OldCoefficient::sin_cos(&OldTerm::var(0).sin());
            }) as fn(),
            (|| {
                let _ = NewAngle::<NewTerm>::sin_cos(&NewTerm::var(0).sin());
            }) as fn(),
        ),
    ] {
        let om = panic_message(old_call);
        let nm = panic_message(new_call);
        assert_eq!(om, nm, "{what}: panic message differs");
        assert!(om.is_some(), "{what}: old did not panic");
    }
}

/// Run `f`, returning its panic message (if it panicked) with the panic hook
/// silenced so the test output stays readable.
fn panic_message(f: fn()) -> Option<String> {
    let prev = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(f));
    std::panic::set_hook(prev);
    match r {
        Ok(()) => None,
        Err(e) => Some(
            e.downcast_ref::<String>()
                .cloned()
                .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_default(),
        ),
    }
}

#[test]
fn contract_9_zero_shortcut_leaves_an_empty_sum_not_a_constant_zero() {
    let (o, n) = zero_shortcut_pair();
    assert_eq!(o.to_string(), "[]");
    assert_eq!(n.to_string(), "[]");
    assert!(!o.is_zero() && !n.is_zero());
    assert_ne!(o, OldTerm::from(0.0));
    assert_ne!(n, NewTerm::from(0.0));
    assert_eq!(old_view(&o).form, Form::MapBacked);
    assert_eq!(new_view(&n).form, Form::MapBacked);
    assert!(matches!(n.inner(), Inner::Sum(_)));
    // …and it still evaluates to 0.
    assert_eq!(o.eval(&[0.3, 0.9]).unwrap(), 0.0);
    assert_eq!(n.eval(&[0.3, 0.9]).unwrap(), 0.0);
}

/// The exported-macro API surface: `pauli_sum *= Term` and `pauli_sum * Term`.
///
/// Old `ppvm-sym` took a real (non-dev) dependency on `ppvm-pauli-sum` for the
/// single line `impl_op_mul_assign_coefficient!(Term)`, so `sum *= Term::from(2)`
/// compiled on a symbolic Pauli sum. `ppvm-sym-2` reproduces it with
/// `ppvm_pauli_sum_2::impl_scalar_mul!(Term)`. This test is as much a
/// *compile*-level assertion as a value one: if the instantiation is dropped,
/// this file stops building.
///
/// Semantics on both sides are "scale every coefficient, remove nothing" —
/// old's `ACMapMulAssign` (`*v *= value.clone()` per entry), so a symbolic
/// coefficient is multiplied by the `Term` in the coefficient ring and a
/// zeroed-out term stays in the support.
#[test]
fn scalar_mul_by_a_term_matches_old() {
    let angles = vec![vec![0.3, 1.1, -0.7], vec![2.2, -1.4, 0.9]];

    let mut os = new_old_sum(4);
    os += ("ZIII", OldTerm::var(0).sin());
    os += ("IZII", OldTerm::from(1.0) + OldTerm::var(1).cos());
    os *= OldTerm::from(2.0);
    os *= OldTerm::var(2).cos();

    let mut ns = new_new_sum(4);
    ns += (NewSymKey::from("ZIII"), NewTerm::var(0).sin());
    ns += (
        NewSymKey::from("IZII"),
        NewTerm::from(1.0) + NewTerm::var(1).cos(),
    );
    ns *= NewTerm::from(2.0);
    ns *= NewTerm::var(2).cos();

    let (ok, nk) = (old_sym_support(&os), new_sym_support(&ns));
    assert_eq!(
        ok.iter().map(|(k, _)| k).collect::<Vec<_>>(),
        nk.iter().map(|(k, _)| k).collect::<Vec<_>>(),
        "support differs after `*= Term`"
    );
    for ((k, oc), (_, nc)) in ok.iter().zip(nk.iter()) {
        assert_terms_match(oc, nc, &angles, &format!("*= Term at {k}"));
    }

    // Scaling by a zero constant keeps the whole key set (nothing is removed),
    // exactly as `*= 0.0` does on an `f64` sum.
    let mut oz = new_old_sum(4);
    oz += ("ZIII", OldTerm::var(0).sin());
    oz *= OldTerm::from(0.0);
    let mut nz = new_new_sum(4);
    nz += (NewSymKey::from("ZIII"), NewTerm::var(0).sin());
    nz *= NewTerm::from(0.0);
    assert_eq!(oz.len(), 1);
    assert_eq!(nz.len(), 1);
    assert_terms_match(
        &old_sym_support(&oz)[0].1,
        &new_sym_support(&nz)[0].1,
        &angles,
        "*= Term::from(0.0)",
    );

    // The by-value `sum * Term` form the macro also emits.
    let mut nm = new_new_sum(4);
    nm += (NewSymKey::from("ZIII"), NewTerm::var(0).sin());
    let nm = nm * NewTerm::from(2.0);
    assert_eq!(nm.len(), 1);
}

#[test]
fn contract_10_trace_readout_resets_truncation_to_the_defaults() {
    // `std::iter::Sum for Term` (what `trace`/`overlap` fold with) starts from
    // `Term::from_f64(0.0)`, whose parameters are the DEFAULTS — so a user who
    // set `max_sin = 1` on the propagated coefficients gets a read-out
    // accumulated with NO `max_sin` bound and a `1e-16` epsilon.
    //
    // The observable consequence: a `One`-form coefficient that escaped the
    // bound (the `One × One → One` fast arm never consults `max_sin`, on either
    // crate) is admitted into the read-out, where an accumulator carrying the
    // user's `max_sin = 1` would have dropped it.
    let build = |k: usize| {
        let mut a = OldTerm::var(0).sin() * OldTerm::var(1).sin() * OldTerm::var(2).sin();
        a.set_max_sin(k);
        let mut b = OldTerm::from(1.0) + OldTerm::var(3).cos();
        b.set_max_sin(k);
        (a, b)
    };
    let build_new = |k: usize| {
        let mut a = NewTerm::var(0).sin() * NewTerm::var(1).sin() * NewTerm::var(2).sin();
        a.set_max_sin(k);
        let mut b = NewTerm::from(1.0) + NewTerm::var(3).cos();
        b.set_max_sin(k);
        (a, b)
    };

    let (oa, ob) = build(1);
    let (na, nb) = build_new(1);
    assert_eq!(old_view(&oa).max_sin_pow(), 3, "the One-form escapee");
    assert_eq!(new_view(&na).max_sin_pow(), 3);

    // The read-out fold: `std::iter::Sum`, i.e. from `Const(0.0)` with MAX.
    let ot: OldTerm = vec![ob.clone(), oa.clone()].into_iter().sum();
    let nt: NewTerm = vec![nb.clone(), na.clone()].into_iter().sum();
    assert_eq!(old_view(&ot).monomials, new_view(&nt).monomials);
    assert!(
        old_view(&ot).max_sin_pow() > 1,
        "the read-out must NOT have applied max_sin = 1: {ot}"
    );
    assert_eq!(nt.max_sin(), usize::MAX, "the read-out resets max_sin");
    assert_eq!(nt.min_eps(), f64::EPSILON, "…and min_eps");

    // Contrast: folding into an accumulator that DOES carry `max_sin = 1` drops
    // the same monomial — on both crates, so the reset is what makes the
    // difference and not a representation quirk.
    let mut oacc = OldTerm::from(0.0);
    oacc.set_max_sin(1);
    oacc += ob;
    oacc += oa;
    let mut nacc = NewTerm::from(0.0);
    nacc.set_max_sin(1);
    nacc += nb;
    nacc += na;
    assert_eq!(
        old_view(&oacc).max_sin_pow(),
        0,
        "a max_sin = 1 accumulator drops the sin³ escapee: {oacc}"
    );
    assert_eq!(new_view(&nacc).max_sin_pow(), 0);
    assert_eq!(old_view(&oacc).monomials, new_view(&nacc).monomials);

    // …and end to end through the engines: old and new agree on the traced
    // `Term`, and the new one reports the reset parameters.
    let spec = TrotterSpec {
        n: 4,
        layers: 3,
        max_sin: 1,
        min_eps: 1e-12,
        observable: "ZIII",
    };
    let os = trotter_old(&spec);
    let ns = trotter_new(&spec);
    let ot = ppvm_traits::traits::Trace::trace(&os, &OldPauliPattern::from("Z?*"));
    let nt = ppvm_traits_2::Trace::trace(&ns, &ppvm_pauli_sum_2::PauliPattern::zero_state());
    assert_eq!(
        old_view(&ot).monomials,
        new_view(&nt).monomials,
        "old {ot}\nnew {nt}"
    );
    assert_eq!(nt.max_sin(), usize::MAX);
    assert_eq!(nt.min_eps(), f64::EPSILON);
}

#[test]
fn contract_11_operator_form_additions_bypass_truncation() {
    // `Sum += f64` does a bare `c0 += rhs` with NO `min_eps` check, while the
    // method `Sum::add_const` DOES drop below `min_eps`.
    let mut o = OldSumTy::new();
    o += 1e-9;
    let mut n = NewSumTy::new();
    n += 1e-9;
    assert_eq!(o.eval(&[]).unwrap(), 1e-9);
    assert_eq!(n.eval(&[]).unwrap(), 1e-9);

    let mut o = OldSumTy::new();
    o.add_const(1e-9, 1e-3);
    let mut n = NewSumTy::new();
    n.add_const(1e-9, 1e-3);
    assert_eq!(o.eval(&[]).unwrap(), 0.0);
    assert_eq!(n.eval(&[]).unwrap(), 0.0);

    // `Sum += Prod` inserts with NO check at all — even a `pow() == 0` monomial,
    // violating the invariant that `terms` only holds `pow > 0`.
    let mut o = OldSumTy::new();
    o += OldProd::new();
    let mut n = NewSumTy::new();
    n += NewProd::new();
    assert_eq!(o.to_string(), "[1.000 * ]");
    assert_eq!(n.to_string(), o.to_string());
    assert_eq!(n.len(), 1);

    // The `Term`-level `+= Term::from(c)` path routes through `add_const` and
    // therefore DOES drop.
    let mut o = OldTerm::from(1.0) + OldTerm::var(0).sin();
    o.set_min_eps(1e-3);
    o += OldTerm::from(1e-9);
    let mut n = NewTerm::from(1.0) + NewTerm::var(0).sin();
    n.set_min_eps(1e-3);
    n += NewTerm::from(1e-9);
    assert_eq!(old_view(&o).c0, "1.000");
    assert_eq!(new_view(&n).c0, old_view(&o).c0);
}

#[test]
fn contract_12_trait_defaults_and_conversions() {
    // Old's `half` is the CONSTANT 0.5 regardless of `self` — a `Halvable` law
    // violation (`x.half() + x.half() == x` fails for every `x != 1`). The new
    // crate deliberately does not implement `Halvable` at all (implementation
    // plan §Phase 5), which is a compile-time absence and therefore cannot be
    // asserted here; what IS asserted is that old's value is the law-violating
    // one, so the divergence is documented rather than inferred.
    assert_eq!(
        OldCoefficient::half(&OldTerm::from(3.0)),
        OldTerm::from(0.5)
    );
    let h = OldCoefficient::half(&OldTerm::from(3.0));
    assert_ne!(
        h.clone() + h,
        OldTerm::from(3.0),
        "old `half` violates the Halvable law, which is why -2 omits it"
    );
    // Nothing on the live `-2` path can reach it: `ppvm-pauli-sum-2` has no
    // projection kernel, so `Halvable` is never required of a symbolic sum.

    // `sin_cos` is two clones of the angle, and is partial in exactly old's way.
    let (os, oc) = OldCoefficient::sin_cos(&OldTerm::var(0));
    let (ns, nc) = NewAngle::<NewTerm>::sin_cos(&NewTerm::var(0));
    assert_eq!(os.to_string(), ns.to_string());
    assert_eq!(oc.to_string(), nc.to_string());
    let (os, oc) = OldCoefficient::sin_cos(&OldTerm::from(0.3));
    let (ns, nc) = NewAngle::<NewTerm>::sin_cos(&NewTerm::from(0.3));
    assert_eq!(os.eval(&[]).unwrap(), ns.eval(&[]).unwrap());
    assert_eq!(oc.eval(&[]).unwrap(), nc.eval(&[]).unwrap());
}

// ===========================================================================
// 4. The deliberate divergences (old is wrong; Lean adjudicates).
// ===========================================================================

#[test]
fn divergence_add_f64_to_a_single_monomial() {
    // `oldSuspectedBugs` #1: old's `AddAssign<f64> for Term` `One` arm builds the
    // promoted `Sum` into a local and never assigns it, so `t += 2.0` is a
    // SILENT NO-OP on the most common non-constant form.
    //
    // Lean: `lean/PPVM/Algebra/GradedMap.lean` `accumulate_comm`/
    // `accumulate_assoc` require the coefficient ring to be an additive monoid,
    // so `x + c == x` for `c != 0` is not a permissible reading. The new crate
    // is correct; the divergence is pinned here so it can never drift silently.
    let mut o = OldTerm::var(0).sin();
    o += 2.0;
    let mut n = NewTerm::var(0).sin();
    n += 2.0;

    assert_eq!(
        o.eval(&[0.5]).unwrap(),
        0.5f64.sin(),
        "old is expected to DROP the addend (bug #1)"
    );
    assert_eq!(
        n.eval(&[0.5]).unwrap(),
        2.0 + 0.5f64.sin(),
        "new must add it"
    );
    // The `Const` and `Sum` receivers agree (old's bug is only in the `One` arm).
    let mut o = OldTerm::from(1.0);
    o += 2.0;
    let mut n = NewTerm::from(1.0);
    n += 2.0;
    assert_eq!(o.eval(&[]).unwrap(), n.eval(&[]).unwrap());
    let mut o = OldTerm::from(1.0) + OldTerm::var(0).sin();
    o += 2.0;
    let mut n = NewTerm::from(1.0) + NewTerm::var(0).sin();
    n += 2.0;
    assert_eq!(o.eval(&[0.5]).unwrap(), n.eval(&[0.5]).unwrap());
}

#[test]
fn divergence_monomial_multiplication_composes_the_phase() {
    // `oldSuspectedBugs` #2: old's `MulAssign<Prod> for Prod` merges the
    // sin/cos maps but never combines `phase`. Lean:
    // `lean/PPVM/Algebra/Twisted.lean` `iPow_add` (`i^a · i^b = i^{a+b}`), which
    // `tmul_assoc` needs.
    let mut oa = OldProd::sin(0);
    oa.add_phase(1);
    let mut ob = OldProd::cos(0);
    ob.add_phase(1);
    let oc = oa * ob;

    let mut na = NewProd::sin(0);
    na.add_phase(1);
    let mut nb = NewProd::cos(0);
    nb.add_phase(1);
    let nc = na * nb;

    // Old cannot report a phase (no accessor), but its *keying* is observable:
    // `i·P · i·Q` and `i·P · Q` must be DIFFERENT monomials, and on old they are
    // not.
    let mut ob2 = OldProd::cos(0);
    ob2.add_phase(0);
    let oc2 = {
        let mut a = OldProd::sin(0);
        a.add_phase(1);
        a * ob2
    };
    assert_eq!(
        oc, oc2,
        "old is expected to mis-key: dropping rhs.phase makes i·P·(i·Q) == i·P·Q"
    );

    let mut nb2 = NewProd::cos(0);
    nb2.add_phase(0);
    let nc2 = {
        let mut a = NewProd::sin(0);
        a.add_phase(1);
        a * nb2
    };
    assert_ne!(nc, nc2, "new must compose the phase (i² ≠ i)");
    assert_eq!(nc.phase(), 2, "i·i = i²");
    assert_eq!(nc2.phase(), 1);
}

#[test]
fn divergence_mul_phase_keeps_the_constant_part_phased() {
    // `oldSuspectedBugs` #3: old's `add_term` short-circuits on `pow() == 0`
    // alone and folds the value into `c0`, throwing the phase away, so
    // multiplying a symbolic sum by `i` leaves its constant term unphased.
    // Lean: `twistedConv_add_left`/`twistedConv_add_right` — the product is
    // additive in each argument, so the constant summand must be phased too.
    let o = OldTerm::from(2.0) + OldTerm::var(0).sin();
    let o = OldComplexCoefficient::mul_phase(&o, 1);
    let n = NewTerm::from(2.0) + NewTerm::var(0).sin();
    let n = n.mul_phase(1);

    // Old's real-valued `eval` sees the unphased constant and the (also
    // unphased, since `eval` ignores the phase) monomial: it returns the whole
    // real value.
    assert_eq!(o.eval(&[0.5]).unwrap(), 2.0 + 0.5f64.sin());
    // The new crate's real `eval` still ignores the phase (parity — every
    // real-valued golden master is untouched)…
    assert_eq!(n.eval(&[0.5]).unwrap(), 2.0 + 0.5f64.sin());
    // …but the phase-aware readout shows the whole value is imaginary, which is
    // what old could not represent.
    let v = n.eval_complex(&[0.5]).unwrap();
    assert!(v.re.abs() < 1e-15, "{v}");
    assert!((v.im - (2.0 + 0.5f64.sin())).abs() < 1e-12, "{v}");

    // Old's structural symptom: the constant summand was absorbed into `c0`,
    // so its rendering still shows a constant part; the new one does not.
    assert_eq!(old_view(&o).c0, "2.000");
    assert_eq!(new_view(&n).c0, "0");
}

#[test]
fn divergence_real_eval_ignores_the_phase_on_both_crates() {
    // `oldSuspectedBugs` #4: neither `Prod::eval` nor `Display for Prod` applies
    // the phase. The new crate PRESERVES that for `eval`/`Display` (so real
    // golden masters are byte-identical) and ADDS `eval_complex` rather than
    // changing `eval`'s type.
    let mut op = OldProd::sin(0);
    op.add_phase(1);
    let mut np = NewProd::sin(0);
    np.add_phase(1);
    assert_eq!(op.eval(&[0.7]).unwrap(), 0.7f64.sin());
    assert_eq!(np.eval(&[0.7]).unwrap(), 0.7f64.sin());
    assert_eq!(op.to_string(), np.to_string());
    assert_eq!(op.to_string(), "sin^1(%0)", "the phase is not printed");
    // The added readout.
    assert_eq!(
        np.eval_complex(&[0.7]).unwrap(),
        num::Complex::new(0.0, 0.7f64.sin())
    );
}

// ===========================================================================
// 5. The `sym.*` integration workloads, end to end.
// ===========================================================================

/// `sym.trace.parametric` — `examples/symbolic.rs`, verbatim, on both engines.
#[test]
fn integration_sym_trace_parametric() {
    const GOLDEN_DISPLAY: &str = "[1.000 * cos^3(%1) + -1.000 * sin^1(%0) sin^2(%1) cos^1(%1) + 1.000 * sin^3(%1) cos^1(%0) cos^1(%1) + 1.000 * sin^1(%0) sin^3(%1) cos^1(%0) + 1.000 * sin^1(%0) sin^3(%1) cos^1(%0) cos^1(%1)]";
    const GOLDEN_VALUE: f64 = 0.188_036_759_177_593_55;

    let o = parametric_trace_old();
    let n = parametric_trace_new();

    // (a) The captured golden master still describes OLD (so the bar is real).
    assert_eq!(o.to_string(), GOLDEN_DISPLAY);
    assert!((o.eval(&[1.1, 2.1]).unwrap() - GOLDEN_VALUE).abs() < 1e-12);

    // (b) New matches old, byte for byte on `Display` and to 1e-12 on `eval`.
    assert_eq!(n.to_string(), GOLDEN_DISPLAY);
    let mut rng = seeded_rng(0x5E1_1E5);
    let mut angles = vec![vec![1.1, 2.1]];
    angles.extend(angle_grid(&mut rng, 2, 32));
    assert_terms_match(&o, &n, &angles, "sym.trace.parametric");

    // (c) The five-monomial structure with `c0 == 0.0`.
    let v = new_view(&n);
    assert_eq!(v.n_monomials(), 5);
    assert_eq!(v.c0, "0");
    let want: BTreeMap<String, String> = [
        ("cos^3(%1)", "1.000"),
        ("sin^1(%0) sin^2(%1) cos^1(%1)", "-1.000"),
        ("sin^3(%1) cos^1(%0) cos^1(%1)", "1.000"),
        ("sin^1(%0) sin^3(%1) cos^1(%0)", "1.000"),
        ("sin^1(%0) sin^3(%1) cos^1(%0) cos^1(%1)", "1.000"),
    ]
    .into_iter()
    .map(|(a, b)| (a.to_string(), b.to_string()))
    .collect();
    assert_eq!(v.monomials, want);
    assert!((n.eval(&[1.1, 2.1]).unwrap() - GOLDEN_VALUE).abs() < 1e-12);
}

/// Compare two propagated symbolic sums key for key.
///
/// Asserts (1) identical support, exact; (2) identical monomial set and printed
/// coefficient per key, exact; (3) identical representation form per key; (4)
/// identical TOTAL monomial count, exact (a differing count means a truncation
/// or coalescing rule moved even if the numbers happen to agree); (5) `eval`
/// agreement within `1e-12` at every supplied angle vector.
#[track_caller]
fn assert_sums_match(
    old: &OldSymSum,
    new: &NewSymSum,
    angles: &[Vec<f64>],
    what: &str,
) -> (usize, usize) {
    let os = old_sym_support(old);
    let ns = new_sym_support(new);
    let ok: Vec<&String> = os.iter().map(|(k, _)| k).collect();
    let nk: Vec<&String> = ns.iter().map(|(k, _)| k).collect();
    assert_eq!(ok, nk, "{what}: support differs");

    let (mut total, mut peak) = (0usize, 0usize);
    for ((k, oc), (_, nc)) in os.iter().zip(ns.iter()) {
        assert_terms_match(oc, nc, angles, &format!("{what} @ {k}"));
        let m = new_view(nc).n_monomials();
        total += m;
        peak = peak.max(m);
    }
    (total, peak)
}

/// `sym.tfim.trotter` — the headline deep symbolic Trotter workload.
#[test]
fn integration_sym_tfim_trotter() {
    let mut rng = seeded_rng(0x77_07_7E);
    for k in [3usize, 4] {
        for observable in ["ZIIIII", "ZZIIII"] {
            let spec = TrotterSpec {
                observable,
                ..TrotterSpec::headline(k)
            };
            let mut angles = vec![fixed_angles(spec.n_vars())];
            angles.extend(angle_grid(&mut rng, spec.n_vars(), 32));

            let os = trotter_old(&spec);
            let ns = trotter_new(&spec);
            let (total, peak) = assert_sums_match(
                &os,
                &ns,
                &angles,
                &format!("sym.tfim.trotter k={k} obs={observable}"),
            );
            assert!(total > 0);
            println!(
                "[sym.tfim.trotter] k={k} obs={observable}: support={} monomials={total} peak={peak}",
                ns.len()
            );

            // (4) The seeded `max_sin` really propagated: every MAP-BACKED
            // monomial obeys it. The `One`-form fast arms never consult
            // `max_sin` — old behaviour, preserved deliberately — so a
            // single-monomial coefficient may exceed it. Both crates are checked
            // through the same view, so a divergence in *which* form a
            // coefficient ends up in would already have failed above.
            for (key, c) in new_sym_support(&ns) {
                let v = new_view(&c);
                if v.form == Form::MapBacked {
                    assert!(
                        v.max_sin_pow() <= k,
                        "{key}: a map-backed monomial survived with sin_pow {} > {k}",
                        v.max_sin_pow()
                    );
                }
            }
        }
    }
}

/// `sym.truncation.sweep` — the `max_sin` / `min_eps` cost-curve gate.
#[test]
fn integration_sym_truncation_sweep() {
    let mut rng = seeded_rng(0x5_9EE7);
    let base = TrotterSpec {
        n: 5,
        layers: 4,
        max_sin: 0,
        min_eps: 1e-12,
        observable: "ZIIII",
    };
    let mut angles = vec![fixed_angles(base.n_vars())];
    angles.extend(angle_grid(&mut rng, base.n_vars(), 16));

    // --- the `max_sin` sweep --------------------------------------------------
    let mut per_k: Vec<BTreeMap<String, BTreeSet<String>>> = Vec::new();
    for k in 1..=5usize {
        let spec = TrotterSpec { max_sin: k, ..base };
        let os = trotter_old(&spec);
        let ns = trotter_new(&spec);
        let (total, _) = assert_sums_match(&os, &ns, &angles, &format!("sweep k={k}"));
        println!(
            "[sym.truncation.sweep] k={k}: support={} monomials={total}",
            ns.len()
        );

        per_k.push(
            new_sym_support(&ns)
                .into_iter()
                .map(|(key, c)| (key, new_view(&c).monomials.keys().cloned().collect()))
                .collect(),
        );
    }
    // Monotonicity: the k-truncated monomials must be a SUBSET of the
    // (k+1)-truncated ones, per key.
    for k in 0..per_k.len() - 1 {
        for (key, small) in &per_k[k] {
            let big = per_k[k + 1].get(key).unwrap_or_else(|| {
                panic!("key {key} vanished going from k={} to {}", k + 1, k + 2)
            });
            assert!(
                small.is_subset(big),
                "k={} monomials for {key} are not a subset of k={}",
                k + 1,
                k + 2
            );
        }
    }

    // --- the `min_eps` sweep --------------------------------------------------
    // Old drops a monomial the instant its accumulated coefficient falls under
    // `min_eps` *at insert time* (not as a post-pass), so the surviving set must
    // match EXACTLY, not merely within tolerance.
    for min_eps in [f64::EPSILON, 1e-12, 1e-6] {
        let spec = TrotterSpec {
            max_sin: 3,
            min_eps,
            ..base
        };
        let os = trotter_old(&spec);
        let ns = trotter_new(&spec);
        let (total, _) =
            assert_sums_match(&os, &ns, &angles, &format!("sweep min_eps={min_eps:e}"));
        println!("[sym.truncation.sweep] min_eps={min_eps:e}: monomials={total}");
    }
}

/// `sym.truncation.sweep`, the `min_eps` half — driven END TO END through both
/// engines on a circuit where the threshold genuinely fires.
///
/// The pure-symbolic Trotter sweep above cannot reach this branch: every
/// coefficient it produces has magnitude ≥ 1, so no `min_eps` below `1.0` ever
/// drops anything and the sweep is vacuous. Mixing SMALL `f64` angles into the
/// symbolic circuit is what makes the coefficients decay — `rz(q, 0.05)`
/// constant-folds to `Const(sin 0.05)` ≈ `5e-2`, and after a few layers the
/// branch weights fall through `1e-6`. Old drops such a monomial the instant its
/// accumulated coefficient falls under `min_eps` *at insert time*, so the
/// surviving set must match EXACTLY — not merely within tolerance.
#[test]
fn integration_sym_min_eps_sweep_fires_end_to_end() {
    let n = 4usize;
    let layers = 5u32;
    // `max_sin` is bounded so the monomial count stays polynomial in the depth;
    // this test is about the OTHER truncation axis, and an unbounded sine degree
    // would make it run for minutes without exercising `min_eps` any harder.
    let max_sin = 3usize;
    let mut rng = seeded_rng(0xE9_5EE9);
    let angles = angle_grid(&mut rng, 2 * layers as usize, 8);

    let run = |min_eps: f64| -> (OldSymSum, NewSymSum) {
        let mut os = new_old_sum(n);
        os += ("ZIII", old_seed_coeff(max_sin, min_eps));
        let mut ns = new_new_sum(n);
        ns += (NewSymKey::from("ZIII"), new_seed_coeff(max_sin, min_eps));
        for l in 0..layers {
            for q in 0..n {
                // A symbolic rotation: builds genuine monomials.
                OldRotationOne::rx(&mut os, q, OldTerm::var(l));
                NewRotationOne::rx(&mut ns, q, NewTerm::var(l));
                // A SMALL concrete rotation: constant-folds to a weight ≈ 5e-2,
                // so the branch coefficients decay geometrically with depth.
                OldRotationOne::rz(&mut os, q, OldTerm::from(0.05));
                NewRotationOne::rz(&mut ns, q, 0.05);
            }
            for q in 0..n - 1 {
                OldClifford::cnot(&mut os, q, q + 1);
                NewClifford::cnot(&mut ns, q, q + 1);
            }
        }
        (os, ns)
    };

    let mut monomials_by_eps = Vec::new();
    for min_eps in [f64::EPSILON, 1e-12, 1e-6, 1e-3] {
        let (os, ns) = run(min_eps);
        let (total, _) = assert_sums_match(
            &os,
            &ns,
            &angles,
            &format!("sym.truncation.sweep min_eps={min_eps:e}"),
        );
        println!(
            "[sym.truncation.sweep/min_eps] {min_eps:e}: support={} monomials={total}",
            ns.len()
        );
        monomials_by_eps.push(total);
    }

    // The branch really fired: a coarser threshold retains strictly fewer
    // monomials. Without this the parity assertions above would be vacuous.
    assert!(
        monomials_by_eps[3] < monomials_by_eps[0],
        "min_eps never dropped anything — the sweep is vacuous: {monomials_by_eps:?}"
    );
    // …and it is monotone in the threshold.
    for w in monomials_by_eps.windows(2) {
        assert!(
            w[1] <= w[0],
            "min_eps retention is not monotone: {monomials_by_eps:?}"
        );
    }
}

/// `sym.expectation.grid` — pay symbolic propagation once, sweep angles cheaply.
#[test]
fn integration_sym_expectation_grid() {
    let spec = TrotterSpec::headline(3);
    let ot = ppvm_traits::traits::Trace::trace(&trotter_old(&spec), &OldPauliPattern::from("Z?*"));
    let nt = ppvm_traits_2::Trace::trace(
        &trotter_new(&spec),
        &ppvm_pauli_sum_2::PauliPattern::zero_state(),
    );

    let mut rng = seeded_rng(0x9_81D);
    let grid = angle_grid(&mut rng, spec.n_vars(), 1000);
    let mut worst = 0.0f64;
    for vals in &grid {
        let o = ot.eval(vals).unwrap();
        let n = nt.eval(vals).unwrap();
        worst = worst.max((o - n).abs());
    }
    println!("[sym.expectation.grid] max |old-new| over 1000 points = {worst:e}");
    assert!(worst < 1e-12, "grid divergence {worst:e}");

    // Error/index semantics: a `vals` slice shorter than the highest variable id
    // is an `Err` with old's message shape — not a panic, not a silent zero.
    // Both must fail, with old's message SHAPE. (Which variable id gets named
    // is a function of monomial-table iteration order, unspecified on both
    // crates; the exact message is pinned deterministically for a
    // single-monomial term in `eval_error_shape_matches_old`.)
    let short = vec![0.1, 0.2];
    let oe = ot.eval(&short).unwrap_err().to_string();
    let ne = nt.eval(&short).unwrap_err().to_string();
    for m in [&oe, &ne] {
        assert!(
            m.starts_with("variable %") && m.ends_with(" not found"),
            "unexpected message: {m}"
        );
    }
}

/// `sym.random.circuit` — deep heterogeneous replay, with the per-phase split
/// (Clifford-only prefix vs rotation-heavy suffix) the baseline asks for.
#[test]
fn integration_sym_random_circuit() {
    let n = 6usize;
    let depth = 150usize;
    let n_vars = 8u32;
    let mut rng = seeded_rng(0xC0FFEE);
    let circuit = random_sym_circuit(&mut rng, n, depth, n_vars);

    let mut os = new_old_sum(n);
    os += ("ZIIIII", old_seed_coeff(3, 1e-12));
    let mut ns = new_new_sum(n);
    ns += (NewSymKey::from("ZIIIII"), new_seed_coeff(3, 1e-12));

    let mut angles = vec![fixed_angles(n_vars as usize)];
    angles.extend(angle_grid(&mut seeded_rng(0xBEEF), n_vars as usize, 16));

    // Phase 1: the Clifford-only prefix (cheap re-key path, coefficient
    // untouched). Phase 2: the full circuit. Both are diffed, so an attribution
    // split is available without a separate microbench.
    let prefix: Vec<_> = circuit
        .iter()
        .copied()
        .filter(|g| {
            matches!(
                g,
                SymGate::H(_) | SymGate::S(_) | SymGate::Cnot(..) | SymGate::Cz(..)
            )
        })
        .collect();
    replay_old(&mut os, &prefix);
    replay_new(&mut ns, &prefix);
    let (total, _) = assert_sums_match(&os, &ns, &angles, "sym.random.circuit/clifford-prefix");
    println!(
        "[sym.random.circuit] prefix: support={} monomials={total}",
        ns.len()
    );

    replay_old(&mut os, &circuit);
    replay_new(&mut ns, &circuit);
    let (total, peak) = assert_sums_match(&os, &ns, &angles, "sym.random.circuit/full");
    println!(
        "[sym.random.circuit] full: support={} monomials={total} peak={peak}",
        ns.len()
    );
    assert!(ns.len() > 1, "the replay should have grown the support");
}

/// `sym.exact.multiply`, the half where old's semantics are sound: the
/// **single-word** operator product `A ← A·P`.
///
/// Old's *sum × sum* `MulAssign` is not a parity target — it calls `map_add`
/// once per rhs term and therefore computes the product **chain** rather than
/// the bilinear sum (already recorded by `pauli_sum_multiply_diff.rs`), on top
/// of the three independent phase bugs. Old's `MulAssign<PauliWord>` path IS a
/// single `map_add` over a bijection and is sound, so it is diffed here on a
/// **real, phase-free** product (`Z`-type operands, so every emitted `iᵏ` is
/// `+1` and old's dropped phase cannot bite).
#[test]
fn integration_sym_exact_multiply_single_word_matches_old() {
    use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, Sum as NewSum};
    use ppvm_pauli_word::word::PauliWord as OldWordT;

    let mut oa = new_old_sum(4);
    oa += ("ZIII", OldTerm::var(0).sin());
    oa += ("IZII", OldTerm::var(1).cos());
    oa += ("ZZII", OldTerm::var(0).sin() + OldTerm::var(2).cos());
    oa *= OldWordT::<[u8; 8]>::from("ZIII");

    type NewSumT = NewSum<HashMapStore<NewSymKey, NewTerm>, NoPolicy>;
    let mut na: NewSumT = NewSumT::new(4);
    na += (NewSymKey::from("ZIII"), NewTerm::var(0).sin());
    na += (NewSymKey::from("IZII"), NewTerm::var(1).cos());
    na += (
        NewSymKey::from("ZZII"),
        NewTerm::var(0).sin() + NewTerm::var(2).cos(),
    );
    na.mul_word_assign(&NewSymKey::from("ZIII"));

    let ok: Vec<String> = old_sym_support(&oa).into_iter().map(|(k, _)| k).collect();
    let nk: Vec<String> = new_sym_support(&na).into_iter().map(|(k, _)| k).collect();
    assert_eq!(ok, nk, "real single-word product support differs");

    let angles = vec![vec![0.3, 1.1, -0.7], vec![2.2, -1.4, 0.9]];
    for ((k, oc), (_, nc)) in old_sym_support(&oa).iter().zip(new_sym_support(&na).iter()) {
        for vals in &angles {
            let o = oc.eval(vals).unwrap();
            let n = nc.eval(vals).unwrap();
            assert!((o - n).abs() < 1e-12, "{k}: old {o} vs new {n}");
        }
        // The monomial set AND the representation form must agree: the engine
        // folds the emitted `iᵏ` through the coefficient's own `mul_phase` (via
        // `ImaginaryUnit::mul_i_pow`), so old's unconditional `Const → One(i⁰,c)`
        // promotion is preserved rather than short-circuited by a `Pos1 =>
        // c.clone()` arm. See
        // `single_word_product_phase_fold_keeps_olds_representation`.
        let (ov, nv) = (old_view(oc), new_view(nc));
        assert_eq!(ov.form, nv.form, "{k}: representation form differs");
        assert_eq!(
            ov.monomials.keys().collect::<Vec<_>>(),
            nv.monomials.keys().collect::<Vec<_>>(),
            "{k}: monomial set differs"
        );
    }
}

/// **Closed divergence** (was a gap, now a parity pin): on the single-word
/// operator product the two engines land in the *same representation*, not
/// merely the same value.
///
/// Old folds the emitted phase with `Coefficient::mul_phase(k)`, whose `Const`
/// arm builds `One(Prod::new() · iᵏ, c)` **unconditionally** — even for `k = 0`.
/// The `-2` engine folds it with `Phase::apply`, which now delegates to
/// `ImaginaryUnit::mul_i_pow`; `Term` overrides that to its own `mul_phase`, so
/// the `i⁰` promotion is preserved. Before that override the `Pos1` arm was
/// `c.clone()` and the coefficient stayed a `Const`, which `Term`'s
/// representational `PartialEq` (behavioural contract 5) and `Display`
/// (`"2.000 * "` vs `"2"`, contract 8) both exposed.
///
/// This pins every arm of the fold, not just `k = 0`: `k = 2` used to take
/// `Phase::apply`'s `-c` arm, which folds the sign onto the `f64` instead of into
/// the monomial's phase byte, and is the same class of divergence.
#[test]
fn single_word_product_phase_fold_keeps_olds_representation() {
    use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, Sum as NewSum};
    use ppvm_pauli_word::word::PauliWord as OldWordT;

    type NewSumT = NewSum<HashMapStore<NewSymKey, NewTerm>, NoPolicy>;

    // `(obs, rhs)`: `Z·Z = +I` emits `i⁰`, `Z·X = iY` emits `i¹`, `Y·X = -iZ`
    // emits `i³` — three of the four phase arms, on a `Const` coefficient.
    for (obs, rhs) in [("ZIII", "ZIII"), ("ZIII", "XIII"), ("YIII", "XIII")] {
        let mut oa = new_old_sum(4);
        oa += (obs, OldTerm::from(2.0));
        oa *= OldWordT::<[u8; 8]>::from(rhs);

        let mut na: NewSumT = NewSumT::new(4);
        na += (NewSymKey::from(obs), NewTerm::from(2.0));
        na.mul_word_assign(&NewSymKey::from(rhs));

        let oc = old_sym_support(&oa)[0].1.clone();
        let nc = new_sym_support(&na)[0].1.clone();

        // The value agrees (old's real-only `eval` ignores the phase byte on both
        // sides, so this is the weak half of the assertion) …
        assert_eq!(oc.eval(&[]).unwrap(), nc.eval(&[]).unwrap());
        assert_eq!(oc.eval(&[]).unwrap(), 2.0);
        // … and so does the representation, which is the point.
        let (ov, nv) = (old_view(&oc), new_view(&nc));
        assert_eq!(ov.form, Form::SingleMonomial, "{obs}·{rhs}: old promotes");
        assert_eq!(nv.form, ov.form, "{obs}·{rhs}: form differs");
        assert_eq!(nv.monomials, ov.monomials, "{obs}·{rhs}: monomials differ");
        assert_eq!(nc.to_string(), oc.to_string(), "{obs}·{rhs}: Display");
        assert_eq!(oc.to_string(), "2.000 * ");
    }
}

/// **Reported divergence** (already documented in `ppvm-sym-2`'s crate docs):
/// `Display`'s tie order among monomials sharing `(sin_pow, cos_pow)` differs
/// from old, because that sort key is *non-total* and ties fall back to the
/// monomial table's iteration order — a function of the digests (changed by the
/// packed-vector layout) and the capacity history (changed by the `mul_term` aux
/// double-buffer).
///
/// This test measures how often it actually bites on the real workloads and
/// asserts that everything the ordering is *about* is identical: the monomial
/// set, the printed coefficients, the representation form.
#[test]
fn display_tie_order_divergence_rate_on_the_real_workloads() {
    let mut differ = 0usize;
    let mut total = 0usize;
    let mut check = |oc: &OldTerm, nc: &NewTerm| {
        total += 1;
        let (ov, nv) = (old_view(oc), new_view(nc));
        assert_eq!(ov.form, nv.form);
        assert_eq!(ov.c0, nv.c0);
        assert_eq!(ov.monomials, nv.monomials);
        if oc.to_string() != nc.to_string() {
            differ += 1;
        }
    };

    for k in [3usize, 4] {
        let spec = TrotterSpec::headline(k);
        let os = trotter_old(&spec);
        let ns = trotter_new(&spec);
        for ((_, oc), (_, nc)) in old_sym_support(&os).iter().zip(new_sym_support(&ns).iter()) {
            check(oc, nc);
        }
    }

    let n = 6usize;
    let circuit = random_sym_circuit(&mut seeded_rng(0xC0FFEE), n, 150, 8);
    let mut os = new_old_sum(n);
    os += ("ZIIIII", old_seed_coeff(3, 1e-12));
    let mut ns = new_new_sum(n);
    ns += (NewSymKey::from("ZIIIII"), new_seed_coeff(3, 1e-12));
    replay_old(&mut os, &circuit);
    replay_new(&mut ns, &circuit);
    for ((_, oc), (_, nc)) in old_sym_support(&os).iter().zip(new_sym_support(&ns).iter()) {
        check(oc, nc);
    }

    println!(
        "[sym Display tie order] byte-differs on {differ} of {total} coefficients \
         ({:.1}%); monomial sets, coefficients and forms identical throughout",
        100.0 * differ as f64 / total as f64
    );
    // The `examples/symbolic.rs` snapshot (whose monomials all have distinct
    // sort keys) must still be byte-identical — that is the contract users
    // actually see.
    assert_eq!(
        parametric_trace_old().to_string(),
        parametric_trace_new().to_string()
    );
}

// ===========================================================================
// 6. Behaviour parity: WHEN the side effects fire.
// ===========================================================================

#[test]
fn parity_no_gate_truncates_on_its_own() {
    // The prime directive's canonical check: apply gates and, WITHOUT calling
    // `truncate()`/`reduce()` anywhere, assert new support == old support after
    // EVERY gate. A new engine that auto-truncated (or auto-reduced) where old
    // does not would show up as a support-size divergence at the first gate that
    // cancels a key.
    let n = 5usize;
    let mut rng = seeded_rng(0x9A7E5);
    let circuit = random_circuit(&mut rng, n, 40);

    let mut os = new_old_sum(n);
    os += ("ZIIII", old_seed_coeff(2, 1e-12));
    let mut ns = new_new_sum(n);
    ns += (NewSymKey::from("ZIIII"), new_seed_coeff(2, 1e-12));

    // The shared `GateOp` generator carries `f64` angles; on a symbolic sum the
    // `f64` widens to a constant-folded `Term` on both crates (old via
    // `impl Into<Coeff>`, new via `Angle<Term> for f64`) — which is itself the
    // parity being checked.
    for (i, g) in circuit.iter().enumerate() {
        match *g {
            GateOp::H(q) => {
                OldClifford::h(&mut os, q);
                NewClifford::h(&mut ns, q);
            }
            GateOp::S(q) => {
                OldClifford::s(&mut os, q);
                NewClifford::s(&mut ns, q);
            }
            GateOp::Cnot(a, b) => {
                OldClifford::cnot(&mut os, a, b);
                NewClifford::cnot(&mut ns, a, b);
            }
            GateOp::Rx(q, theta) => {
                OldRotationOne::rx(&mut os, q, OldTerm::from(theta));
                NewRotationOne::rx(&mut ns, q, theta);
            }
            GateOp::Rz(q, theta) => {
                OldRotationOne::rz(&mut os, q, OldTerm::from(theta));
                NewRotationOne::rz(&mut ns, q, theta);
            }
        }
        let ok: Vec<String> = old_sym_support(&os).into_iter().map(|(k, _)| k).collect();
        let nk: Vec<String> = new_sym_support(&ns).into_iter().map(|(k, _)| k).collect();
        assert_eq!(
            ok, nk,
            "support diverged after gate {i} ({g:?}) — did a gate truncate on its own?"
        );
    }
    assert!(os.data().len() > 1, "the circuit should have branched");
}

#[test]
fn parity_exact_cancellation_keeps_the_key_on_both_crates() {
    // `is_zero()` is false for every non-`Const` form, and an exactly-cancelling
    // monomial stays in the table with coefficient `0.0`, so neither crate may
    // drop the key.
    let mut os = new_old_sum(2);
    os += ("ZI", OldTerm::var(0).sin());
    os += ("ZI", OldTerm::var(0).sin() * -1.0);
    let mut ns = new_new_sum(2);
    ns += (NewSymKey::from("ZI"), NewTerm::var(0).sin());
    ns += (NewSymKey::from("ZI"), NewTerm::var(0).sin() * -1.0);

    assert_eq!(os.data().len(), 1);
    assert_eq!(ns.len(), 1);
    let oc = old_sym_support(&os)[0].1.clone();
    let nc = new_sym_support(&ns)[0].1.clone();
    assert_terms_match(&oc, &nc, &[vec![0.7]], "exact cancellation");
    assert_eq!(old_view(&oc).n_monomials(), 1, "the zero monomial survives");

    // …and a sum-level `truncate()` still does not remove it, because a
    // symbolic coefficient is never below any threshold (contract 3).
    os.truncate();
    ns.truncate();
    assert_eq!(os.data().len(), 1);
    assert_eq!(ns.len(), 1);
}
