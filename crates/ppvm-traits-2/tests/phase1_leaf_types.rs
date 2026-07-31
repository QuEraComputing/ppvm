// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lean-oracle-style unit tests pinning the concrete leaf types and provided
//! impls of `ppvm-traits-2` (implementation-plan Phase 1).
//!
//! There is no old behavioral twin for pure trait defs, so these tests encode
//! the algebraic laws the design/Lean spec state directly on the shipped
//! `f64` / `Complex<f64>` witnesses and the leaf types (`Phase`, hashers), and
//! type-check the `Clifford` blanket impl against a stub.

use num::{Complex, One, Zero};
use std::hash::{BuildHasher, Hasher};

use ppvm_traits_2::algebra::{Conjugate, ImaginaryUnit, Phase};
use ppvm_traits_2::coefficient::{Angle, Coefficient, Halvable};
use ppvm_traits_2::gates::Clifford;
use ppvm_traits_2::hash::{IdentityBuildHasher, IdentityHasher};
use ppvm_traits_2::pauli::{PhaseTrack, SymplecticColumns};

type C = Complex<f64>;

// ---------------------------------------------------------------------------
// ImaginaryUnit: i * i == -one()  (lean Matrix.lean `iU_sq`)
// ---------------------------------------------------------------------------

#[test]
fn imaginary_unit_squares_to_neg_one() {
    let i = <C as ImaginaryUnit>::imaginary_unit();
    assert_eq!(i * i, -<C as One>::one());
    // i is a *primitive* 4th root of unity: i^4 == 1, i^2 != 1.
    assert_eq!(i * i * i * i, <C as One>::one());
    assert_ne!(i * i, <C as One>::one());
    // The distinguished element really is (0, 1).
    assert_eq!(i, C::new(0.0, 1.0));
}

// ---------------------------------------------------------------------------
// Conjugate: involution laws over Complex<f64>, identity over f64,
// conj(i) == -i  (lean Matrix.lean `star_iU`).
// ---------------------------------------------------------------------------

#[test]
fn conjugate_is_an_involution_on_complex() {
    for &(re, im) in &[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (3.0, -2.5), (-7.0, 4.0)] {
        let x = C::new(re, im);
        assert_eq!(Conjugate::conj(&Conjugate::conj(&x)), x);
        // conj negates the imaginary part.
        assert_eq!(Conjugate::conj(&x), C::new(re, -im));
    }
}

#[test]
fn conjugate_is_a_ring_homomorphism_on_complex() {
    let a = C::new(2.0, -3.0);
    let b = C::new(-1.0, 5.0);
    assert_eq!(
        Conjugate::conj(&(a + b)),
        Conjugate::conj(&a) + Conjugate::conj(&b)
    );
    assert_eq!(
        Conjugate::conj(&(a * b)),
        Conjugate::conj(&a) * Conjugate::conj(&b)
    );
}

#[test]
fn conjugate_is_identity_on_reals() {
    for &x in &[0.0_f64, 1.0, -3.5, 42.0] {
        assert_eq!(Conjugate::conj(&x), x);
    }
}

#[test]
fn conjugate_of_i_is_neg_i() {
    let i = <C as ImaginaryUnit>::imaginary_unit();
    assert_eq!(Conjugate::conj(&i), -i);
}

// ---------------------------------------------------------------------------
// Phase behaves as Z/4 and matches the {1, i, -1, -i} interpretation.
//
// The group structure is exercised through the first-class `Phase::one` /
// `Phase::compose` / `Phase::inverse` API (and its `Mul`/`MulAssign` spelling),
// cross-checked against the packed Z/4 exponent and against multiplication of
// the concrete {1, i, -1, -i} values produced by `Phase::apply`.
// ---------------------------------------------------------------------------

const ALL_PHASES: [Phase; 4] = [Phase::Pos1, Phase::PosI, Phase::Neg1, Phase::NegI];

/// The {1, i, -1, -i} value a phase denotes.
fn value(p: Phase) -> C {
    p.apply(&<C as One>::one())
}

#[test]
fn phase_exponent_roundtrips() {
    for p in ALL_PHASES {
        assert_eq!(Phase::from_exponent(p.exponent()), p);
    }
    for k in 0u8..8 {
        assert_eq!(Phase::from_exponent(k).exponent(), k & 3);
    }
}

#[test]
fn phase_value_interpretation() {
    assert_eq!(value(Phase::Pos1), C::new(1.0, 0.0));
    assert_eq!(value(Phase::PosI), C::new(0.0, 1.0));
    assert_eq!(value(Phase::Neg1), C::new(-1.0, 0.0));
    assert_eq!(value(Phase::NegI), C::new(0.0, -1.0));
}

#[test]
fn phase_is_z_mod_4_group() {
    let one = Phase::one();
    assert_eq!(one, Phase::Pos1);
    for a in ALL_PHASES {
        // Identity element.
        assert_eq!(a.compose(one), a);
        assert_eq!(one.compose(a), a);

        // Inverse: exponent 4 - k (mod 4).
        let inv = a.inverse();
        assert_eq!(inv, Phase::from_exponent(4 - a.exponent()));
        assert_eq!(a.compose(inv), one);
        assert_eq!(inv.compose(a), one);

        for b in ALL_PHASES {
            // `compose`, `Mul`, and mod-4 exponent addition all agree.
            let composed = a.compose(b);
            assert_eq!(composed, a * b);
            assert_eq!(composed, Phase::from_exponent(a.exponent() + b.exponent()));

            // `MulAssign` mirrors `compose`.
            let mut acc = a;
            acc *= b;
            assert_eq!(acc, composed);

            // Composition matches multiplication of the denoted values.
            assert_eq!(value(composed), value(a) * value(b));
            // Commutative (Z/4 is abelian).
            assert_eq!(composed, b.compose(a));

            // Associativity.
            for c in ALL_PHASES {
                assert_eq!(a.compose(b).compose(c), a.compose(b.compose(c)));
            }
        }
    }
}

#[test]
fn phase_apply_folds_onto_coefficient() {
    // apply(c) == value(phase) * c for an arbitrary coefficient.
    let c = C::new(2.0, -3.0);
    for p in ALL_PHASES {
        assert_eq!(p.apply(&c), value(p) * c);
    }
}

// ---------------------------------------------------------------------------
// Angle<f64> for f64: sin_cos matches f64::sin_cos.
// ---------------------------------------------------------------------------

#[test]
fn angle_sin_cos_matches_std() {
    use std::f64::consts::PI;
    for &theta in &[
        0.0_f64,
        0.5,
        1.0,
        -1.0,
        PI / 6.0,
        PI / 4.0,
        PI,
        2.0 * PI,
        123.456,
    ] {
        let (s, c) = Angle::<f64>::sin_cos(&theta);
        let (rs, rc) = theta.sin_cos();
        assert_eq!(s, rs);
        assert_eq!(c, rc);
    }
}

// ---------------------------------------------------------------------------
// IdentityHasher / IdentityBuildHasher: write_u64(n) then finish() == n.
// ---------------------------------------------------------------------------

#[test]
fn identity_hasher_is_pass_through() {
    for &n in &[0u64, 1, 42, u64::MAX, 0x0123_4567_89ab_cdef, 1u64 << 63] {
        let mut h = IdentityHasher::default();
        h.write_u64(n);
        assert_eq!(h.finish(), n);
    }
}

#[test]
#[should_panic(expected = "Indexable keys write exactly one u64")]
fn identity_hasher_rejects_byte_writes() {
    // The pass-through contract only admits `write_u64` (an `Indexable` key
    // writes its finalized digest as a single u64); any other write path is a
    // contract violation and must not silently corrupt the digest.
    let mut h = IdentityHasher::default();
    h.write(&[1, 2, 3]);
}

#[test]
fn identity_build_hasher_produces_pass_through_hashers() {
    let bh = IdentityBuildHasher;
    let n = 0xdead_beef_cafe_f00d;
    let mut h = bh.build_hasher();
    h.write_u64(n);
    assert_eq!(h.finish(), n);

    // Independent hashers from the same builder do not share state.
    let mut h2 = bh.build_hasher();
    assert_eq!(h2.finish(), 0);
    h2.write_u64(7);
    assert_eq!(h2.finish(), 7);
}

// ---------------------------------------------------------------------------
// Coefficient ring smoke test over f64 (a trivial reference impl).
// ---------------------------------------------------------------------------

#[test]
fn coefficient_ring_ops_on_f64() {
    let a = 3.0_f64;
    let b = -2.0_f64;

    // Ring identities.
    assert_eq!(a + <f64 as Zero>::zero(), a);
    assert_eq!(a + (-a), <f64 as Zero>::zero());
    assert_eq!((a + b) - b, a);
    assert_eq!(a * (b + 1.0), a * b + a); // distributivity

    // Provided methods.
    assert_eq!(Coefficient::mul_sign(&a, -1), -a);
    assert_eq!(Coefficient::mul_sign(&a, 1), a);
    assert_eq!(Halvable::half(&a), a / 2.0);
    assert_eq!(Halvable::half(&a) + Halvable::half(&a), a); // exact
    assert_eq!(Coefficient::magnitude(&b), 2.0);
    assert_eq!(Coefficient::magnitude(&a), 3.0);

    // AddAssign / MulAssign / Sum.
    let mut acc = 0.0_f64;
    acc += a;
    acc *= 2.0;
    assert_eq!(acc, 6.0);
    let s: f64 = [1.0_f64, 2.0, 3.0].into_iter().sum();
    assert_eq!(s, 6.0);
}

#[test]
fn coefficient_ring_ops_on_complex() {
    let a = C::new(1.0, 2.0);
    assert_eq!(Coefficient::mul_sign(&a, -1), -a);
    assert_eq!(Halvable::half(&a), a / 2.0);
    assert_eq!(Halvable::half(&a) + Halvable::half(&a), a); // exact
    assert_eq!(Coefficient::magnitude(&C::new(3.0, 4.0)), 5.0);
}

// ---------------------------------------------------------------------------
// Compile-only: a stub SymplecticColumns + PhaseTrack automatically gets
// Clifford via the blanket impl in pauli.rs. Recording each primitive call
// lets us also assert the blanket impl fires the documented sequence.
// ---------------------------------------------------------------------------

#[derive(Default)]
struct StubPauli {
    log: Vec<String>,
}

impl SymplecticColumns for StubPauli {
    fn n_qubits(&self) -> usize {
        1
    }
    fn swap_xz(&mut self, q: usize) {
        self.log.push(format!("swap_xz({q})"));
    }
    fn xor_z_from_x(&mut self, q: usize) {
        self.log.push(format!("xor_z_from_x({q})"));
    }
    fn xor_x_col(&mut self, ctrl: usize, tgt: usize) {
        self.log.push(format!("xor_x_col({ctrl},{tgt})"));
    }
    fn xor_z_col(&mut self, tgt: usize, ctrl: usize) {
        self.log.push(format!("xor_z_col({tgt},{ctrl})"));
    }
    fn cz_bits(&mut self, a: usize, b: usize) {
        self.log.push(format!("cz_bits({a},{b})"));
    }
}

impl PhaseTrack for StubPauli {
    fn flip_phase_where_xz(&mut self, q: usize) {
        self.log.push(format!("flip_phase_where_xz({q})"));
    }
    fn s_phase(&mut self, q: usize) {
        self.log.push(format!("s_phase({q})"));
    }
    fn cnot_phase(&mut self, ctrl: usize, tgt: usize) {
        self.log.push(format!("cnot_phase({ctrl},{tgt})"));
    }
    fn cz_phase(&mut self, a: usize, b: usize) {
        self.log.push(format!("cz_phase({a},{b})"));
    }
    fn x_phase(&mut self, q: usize) {
        self.log.push(format!("x_phase({q})"));
    }
    fn y_phase(&mut self, q: usize) {
        self.log.push(format!("y_phase({q})"));
    }
    fn z_phase(&mut self, q: usize) {
        self.log.push(format!("z_phase({q})"));
    }
}

/// Free function usable only if `StubPauli: Clifford` — the compile-only guard.
fn requires_clifford<T: Clifford>(_t: &T) {}

#[test]
fn stub_gets_clifford_blanket_impl() {
    let mut s = StubPauli::default();
    requires_clifford(&s);

    // H = flip_phase_where_xz then swap_xz.
    s.h(0);
    assert_eq!(s.log, vec!["flip_phase_where_xz(0)", "swap_xz(0)"]);

    // S = s_phase then xor_z_from_x.
    s.log.clear();
    s.s(0);
    assert_eq!(s.log, vec!["s_phase(0)", "xor_z_from_x(0)"]);

    // CNOT = cnot_phase, xor_x_col, xor_z_col.
    s.log.clear();
    s.cnot(0, 1);
    assert_eq!(
        s.log,
        vec!["cnot_phase(0,1)", "xor_x_col(0,1)", "xor_z_col(1,0)"]
    );

    // CZ = cz_phase then cz_bits.
    s.log.clear();
    s.cz(0, 1);
    assert_eq!(s.log, vec!["cz_phase(0,1)", "cz_bits(0,1)"]);

    // Pure-sign Paulis touch only the phase.
    s.log.clear();
    s.x(0);
    s.y(0);
    s.z(0);
    assert_eq!(s.log, vec!["x_phase(0)", "y_phase(0)", "z_phase(0)"]);

    // stim aliases route to cnot/cz.
    s.log.clear();
    s.cx(0, 1);
    s.zcx(0, 1);
    s.zcz(0, 1);
    assert_eq!(
        s.log,
        vec![
            "cnot_phase(0,1)",
            "xor_x_col(0,1)",
            "xor_z_col(1,0)",
            "cnot_phase(0,1)",
            "xor_x_col(0,1)",
            "xor_z_col(1,0)",
            "cz_phase(0,1)",
            "cz_bits(0,1)",
        ]
    );
}

// ---------------------------------------------------------------------------
// Exact-ring witness: ℤ[i] (Gaussian integers) as a `Coefficient` +
// `ImaginaryUnit` + `Conjugate`, with NO floats.
//
// This mirrors the Lean `GaussianInt` instance (`lean/PPVM/Pauli/Matrix.lean`:
// `iU := ⟨0, 1⟩`, `StarRing GaussianInt`) and pins the design's headline claim
// that dropping `Mul<f64>` lets an *exact* ring be a `Coefficient` (L4 admits
// exact rings). Crucially ℤ[i] implements `Coefficient` WITHOUT implementing
// `Halvable`: halving (`0.5·x`) is not total on ℤ[i] — `half(1+i)` would leave
// the ring — so it was moved to the `Halvable` capability (like `sin_cos`→`Angle`
// and `i`→`ImaginaryUnit`), and propagation never needs it. A `Coefficient`-only
// generic below accepts `Zi`, pinning that halving is off the propagation path.
// ---------------------------------------------------------------------------

/// A Gaussian integer `re + im·i` ∈ ℤ[i]. Exact; no floating point in its ring
/// arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Zi {
    re: i64,
    im: i64,
}

impl Zi {
    const fn new(re: i64, im: i64) -> Self {
        Zi { re, im }
    }
}

impl num::Zero for Zi {
    fn zero() -> Self {
        Zi::new(0, 0)
    }
    fn is_zero(&self) -> bool {
        self.re == 0 && self.im == 0
    }
}

impl num::One for Zi {
    fn one() -> Self {
        Zi::new(1, 0)
    }
}

impl std::ops::Neg for Zi {
    type Output = Zi;
    fn neg(self) -> Zi {
        Zi::new(-self.re, -self.im)
    }
}

impl std::ops::Add for Zi {
    type Output = Zi;
    fn add(self, rhs: Zi) -> Zi {
        Zi::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl std::ops::Sub for Zi {
    type Output = Zi;
    fn sub(self, rhs: Zi) -> Zi {
        Zi::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl std::ops::Mul for Zi {
    type Output = Zi;
    fn mul(self, rhs: Zi) -> Zi {
        // (a + bi)(c + di) = (ac − bd) + (ad + bc)i
        Zi::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl std::ops::AddAssign for Zi {
    fn add_assign(&mut self, rhs: Zi) {
        *self = *self + rhs;
    }
}

impl std::ops::MulAssign for Zi {
    fn mul_assign(&mut self, rhs: Zi) {
        *self = *self * rhs;
    }
}

impl std::iter::Sum for Zi {
    fn sum<I: Iterator<Item = Zi>>(iter: I) -> Zi {
        iter.fold(<Zi as num::Zero>::zero(), |a, b| a + b)
    }
}

impl Coefficient for Zi {
    fn mul_sign(&self, sign: i8) -> Self {
        Zi::new((sign as i64) * self.re, (sign as i64) * self.im)
    }
    // NOTE: no `half`. `Halvable` (`0.5·x`) is a separate capability precisely
    // because ℤ[i] cannot implement it exactly — `half(1+i)` would leave the
    // ring. An exact `Coefficient` need not be `Halvable`.
    fn magnitude(&self) -> f64 {
        ((self.re * self.re + self.im * self.im) as f64).sqrt()
    }
}

impl ImaginaryUnit for Zi {
    fn imaginary_unit() -> Self {
        Zi::new(0, 1)
    }
}

impl Conjugate for Zi {
    fn conj(&self) -> Self {
        Zi::new(self.re, -self.im)
    }
}

#[test]
fn gaussian_integer_is_an_exact_coefficient_ring() {
    use num::{One, Zero};
    let one = <Zi as One>::one();
    let i = <Zi as ImaginaryUnit>::imaginary_unit();

    // Headline law, exact and float-free: i·i == −1 and i⁴ == 1, i² != 1.
    assert_eq!(i * i, -one);
    assert_eq!(i * i * i * i, one);
    assert_ne!(i * i, one);

    // *-ring involution laws, exact.
    let a = Zi::new(3, -2);
    let b = Zi::new(-1, 5);
    assert_eq!(Conjugate::conj(&Conjugate::conj(&a)), a);
    assert_eq!(Conjugate::conj(&(a + b)), a.conj() + b.conj());
    assert_eq!(Conjugate::conj(&(a * b)), a.conj() * b.conj());
    // conj(i) == −i  (Lean Matrix.lean `star_iU`).
    assert_eq!(Conjugate::conj(&i), -i);

    // Ring identities.
    assert_eq!(a + <Zi as Zero>::zero(), a);
    assert_eq!(a + (-a), <Zi as Zero>::zero());
    assert_eq!((a + b) - b, a);
    assert_eq!(a * (b + one), a * b + a); // distributivity
    assert_eq!(a * b, b * a); // commutative ring (required by L4/tmul_assoc)

    // mul_sign / magnitude behave.
    assert_eq!(Coefficient::mul_sign(&a, -1), -a);
    assert_eq!(Coefficient::mul_sign(&a, 1), a);
    assert_eq!(Coefficient::magnitude(&Zi::new(3, 4)), 5.0);

    // Phase::apply folds an exact iᵏ onto an exact coefficient — no floats.
    assert_eq!(Phase::NegI.apply(&one), -i);
    assert_eq!(Phase::PosI.apply(&a), a * i);
    assert_eq!(Phase::Neg1.apply(&a), -a);
}

/// Accepts any value-domain `Coefficient`, with NO `Halvable` bound — the
/// propagation path (Clifford / rotation / twisted product) never halves. That
/// this instantiates at `Zi` is the compile-time witness that an exact ring with
/// no total `0.5·x` is a full `Coefficient`.
fn accepts_any_coefficient<T: Coefficient>(a: &T, b: &T) -> T {
    a.clone() * b.clone() + a.clone()
}

#[test]
fn gaussian_integer_is_a_coefficient_without_being_halvable() {
    // `Zi` has no exact `half`, yet it is a `Coefficient`: the value domain and
    // the `Halvable` capability are cleanly separated (the fix for the retained
    // `half()` that would otherwise force a lossy `/2` on ℤ[i]).
    let a = Zi::new(1, 1);
    let b = Zi::new(2, -3);
    assert_eq!(accepts_any_coefficient(&a, &b), a * b + a);

    // The float/complex measurement domains keep an *exact* `Halvable`, so the
    // `(I ± Z)/2` projector kernel still type-checks there.
    assert_eq!(Halvable::half(&4.0_f64) + Halvable::half(&4.0_f64), 4.0);
    let c = C::new(1.0, 2.0);
    assert_eq!(Halvable::half(&c) + Halvable::half(&c), c);
}
