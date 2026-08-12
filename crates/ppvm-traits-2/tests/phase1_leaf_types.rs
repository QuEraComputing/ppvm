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

use ppvm_traits_2::algebra::{Conjugate, ImaginaryUnit, KeyProduct, Phase};
use ppvm_traits_2::coefficient::{Angle, Coefficient, Halvable};
use ppvm_traits_2::gates::Clifford;
use ppvm_traits_2::hash::{IdentityBuildHasher, IdentityHasher, Indexable};
use ppvm_traits_2::pauli::{BlanketClifford, PhaseTrack, SymplecticColumns};
use ppvm_traits_2::word::{Pauli, PauliBits, Word};

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
// Angle<Complex<f64>> for Complex<f64>: behaviour parity with the old
// `Coefficient::sin_cos` impl for `Complex<f64>` — the angle is the real part
// and the amplitudes come back purely real.
// ---------------------------------------------------------------------------

#[test]
fn complex_angle_sin_cos_matches_old_coefficient_impl() {
    use std::f64::consts::PI;
    for &(re, im) in &[
        (0.0_f64, 0.0_f64),
        (0.5, 0.0),
        (PI / 4.0, 0.0),
        (PI, 0.0),
        // The old impl discarded the imaginary part; so must this one.
        (1.0, 2.5),
        (-1.0, -0.25),
    ] {
        let theta = C::new(re, im);
        let (s, c) = Angle::<C>::sin_cos(&theta);
        let (rs, rc) = re.sin_cos();
        assert_eq!(s, C::new(rs, 0.0));
        assert_eq!(c, C::new(rc, 0.0));
    }
}

// ---------------------------------------------------------------------------
// Angle<Complex<f64>> for f64: the *real-angle-over-complex-coefficient*
// instantiation. This impl exists purely for behaviour parity — the old
// rotation methods took `theta: impl Into<T::Coeff>` with `Coefficient:
// From<f64>`, so `sum.rx(0, 0.1)` compiled on a `Complex<f64>` sum by widening
// the `f64` to `Complex::new(0.1, 0.0)` and then dropping the (zero) imaginary
// part again in `sin_cos`. The redesigned `RotationOne` cannot take
// `impl Into<A>`, so the caller spelling is preserved by this impl instead; if
// its body ever drifted from "widen, then old sin_cos", every `f64`-angle call
// on a complex-coefficient sum would silently change amplitude.
// ---------------------------------------------------------------------------

/// The old `ppvm_traits::Coefficient::sin_cos` for `Complex<f64>`, verbatim
/// (`crates/ppvm-traits/src/traits/coefficient.rs`).
fn old_complex_sin_cos(theta: C) -> (C, C) {
    let (s, c) = num::traits::Float::sin_cos(theta.re);
    (C::new(s, 0.0), C::new(c, 0.0))
}

#[test]
fn real_angle_over_complex_coefficients_matches_the_old_widening() {
    use std::f64::consts::PI;
    for &theta in &[0.0_f64, 0.1, -0.1, 0.5, PI / 4.0, PI, -2.5, 123.456] {
        // The old call path: `theta.into()` (`Complex::from`) then the old
        // `Coefficient::sin_cos`.
        let want = old_complex_sin_cos(C::from(theta));
        let got = Angle::<C>::sin_cos(&theta);
        assert_eq!(got, want, "f64 angle over Complex broke at {theta}");

        // Equivalently: routing the *widened* value through the `A = C`
        // instantiation gives the same amplitudes, so a caller may pass either
        // spelling and get identical numbers.
        assert_eq!(got, Angle::<C>::sin_cos(&C::new(theta, 0.0)));

        // The amplitudes are purely real and on the unit circle (`rot_norm_sq`).
        assert_eq!(got.0.im, 0.0);
        assert_eq!(got.1.im, 0.0);
        assert!((got.0.re * got.0.re + got.1.re * got.1.re - 1.0).abs() < 1e-12);
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
// Compile-only: a stub SymplecticColumns + PhaseTrack that opts into
// BlanketClifford automatically gets Clifford via the blanket impl in pauli.rs.
// Recording each primitive call lets us also assert the blanket impl fires the
// documented sequence.
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

// Opt into the shared blanket Clifford.
impl BlanketClifford for StubPauli {}

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

// ---------------------------------------------------------------------------
// Behaviour parity: `Phase::apply` reproduces the old
// `ppvm_traits::ComplexCoefficient::mul_phase` packing verbatim.
//
// The old crate packed a fourth root of unity as a `u8` with the table
// `0 → +1, 1 → +i, 2 → −1, 3 → −i` (`crates/ppvm-traits/src/traits/coefficient.rs`)
// and multiplied by hand on the real/imaginary parts. `Phase` replaces that
// `u8` with an enum whose `exponent()` MUST use the same numbering (the Lean
// `Phase` numbering of `lean/PPVM/Pauli/Phase.lean` is likewise
// `Pos1 = 0, PosI = 1, Neg1 = 2, NegI = 3`), otherwise every phased word in the
// port silently changes sign.
// ---------------------------------------------------------------------------

/// The old `ComplexCoefficient::mul_phase`, copied verbatim from
/// `crates/ppvm-traits/src/traits/coefficient.rs` as the parity oracle.
fn old_mul_phase(c: C, phase: u8) -> C {
    match phase % 4 {
        0 => c,
        1 => C::new(-c.im, c.re),
        2 => -c,
        3 => C::new(c.im, -c.re),
        _ => unreachable!(),
    }
}

/// Bit-exact comparison, so that `+0.0` vs `-0.0` (equal under `==`, but
/// distinguishable through `Display`/serialization) and the *placement* of a
/// `NaN` component are both caught. `NaN` is compared by payload-agnostic
/// "is it NaN", per component.
fn same_complex(a: C, b: C) -> bool {
    fn same_f64(x: f64, y: f64) -> bool {
        if x.is_nan() || y.is_nan() {
            x.is_nan() && y.is_nan()
        } else {
            x.to_bits() == y.to_bits()
        }
    }
    same_f64(a.re, b.re) && same_f64(a.im, b.im)
}

#[test]
fn phase_apply_matches_old_mul_phase_encoding() {
    for &c in &[
        C::new(0.0, 0.0),
        C::new(1.0, 0.0),
        C::new(0.0, 1.0),
        C::new(2.0, -3.0),
        C::new(-1.5, 0.25),
        // Signed zeros: `c * i` computes `re·0 − im·1`, which turns `+0` into
        // `-0` (and vice versa) where the old component swap did not.
        C::new(-0.0, 0.0),
        C::new(0.0, -0.0),
        C::new(1.0, -0.0),
        // Non-finite components: `inf·0` and `NaN·0` are `NaN`, so a generic
        // complex multiply would spread the contamination into the *other*
        // component. The old `mul_phase` swap kept it where it was — and since
        // `magnitude() < threshold` is false for `NaN`, truncation never removes
        // such a term, so the contamination would be permanent.
        C::new(f64::INFINITY, 0.0),
        C::new(f64::NEG_INFINITY, 0.0),
        C::new(0.0, f64::INFINITY),
        C::new(f64::NAN, 1.0),
        C::new(1.0, f64::NAN),
        C::new(f64::INFINITY, f64::NAN),
    ] {
        for p in ALL_PHASES {
            // Same exponent numbering AND the same resulting value.
            let (got, want) = (p.apply(&c), old_mul_phase(c, p.exponent()));
            assert!(same_complex(got, want), "{p:?} on {c}: {got} != {want}");
        }
        // The u8 table itself round-trips through the enum.
        for k in 0u8..4 {
            let (got, want) = (Phase::from_exponent(k).apply(&c), old_mul_phase(c, k));
            assert!(same_complex(got, want), "i^{k} on {c}: {got} != {want}");
        }
    }
}

/// `ImaginaryUnit::mul_i` is `·i` — the law the default body states — checked on
/// the finite values where the generic product is itself well behaved. This pins
/// that the `Complex<f64>` override (a component swap) is not a *different*
/// function, only a total one.
#[test]
fn mul_i_agrees_with_multiplying_by_the_imaginary_unit_on_finite_values() {
    for &c in &[
        C::new(0.0, 0.0),
        C::new(1.0, 0.0),
        C::new(0.0, 1.0),
        C::new(2.0, -3.0),
        C::new(-1.5, 0.25),
    ] {
        assert_eq!(c.mul_i(), c * <C as ImaginaryUnit>::imaginary_unit());
        // Four applications are the identity: `i⁴ = 1`.
        assert_eq!(c.mul_i().mul_i().mul_i().mul_i(), c);
    }
}

// ---------------------------------------------------------------------------
// Behaviour parity: `magnitude()` thresholded by a policy reproduces the old
// `Coefficient::cutoff(threshold)` decision exactly.
//
// The old trait owned the cutoff *decision* (`self.abs() < threshold` for f64,
// `self.norm() < threshold` for Complex); the redesign moves the decision to
// `Policy` and leaves only the value property here. The composed predicate must
// still agree bit for bit, including at the boundary and on non-finite inputs
// (old `<` is false for NaN, so a NaN coefficient was never cut — preserved).
// ---------------------------------------------------------------------------

/// The old `Coefficient::cutoff` for `f64`, copied verbatim as the oracle.
fn old_cutoff_f64(x: f64, threshold: f64) -> bool {
    x.abs() < threshold
}

/// The old `Coefficient::cutoff` for `Complex<f64>`, copied verbatim.
fn old_cutoff_complex(x: C, threshold: f64) -> bool {
    x.norm() < threshold
}

#[test]
fn magnitude_threshold_matches_old_cutoff_f64() {
    let thresholds = [0.0_f64, 1e-12, 0.5, 1.0, 1e9];
    let values = [
        0.0_f64,
        -0.0,
        1e-13,
        -0.5,
        0.5,
        -1.0,
        3.25,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ];
    for &t in &thresholds {
        for &v in &values {
            let new_cut = Coefficient::magnitude(&v) < t;
            assert_eq!(
                new_cut,
                old_cutoff_f64(v, t),
                "magnitude/cutoff parity broke at v={v}, t={t}"
            );
        }
    }
    // The boundary is exclusive on both sides: |x| == t is *kept*.
    assert_eq!(Coefficient::magnitude(&0.5_f64), 0.5);
    assert!(!old_cutoff_f64(0.5, 0.5));
}

#[test]
fn magnitude_threshold_matches_old_cutoff_complex() {
    let thresholds = [0.0_f64, 1e-12, 1.0, 5.0, 6.0];
    let values = [
        C::new(0.0, 0.0),
        C::new(3.0, 4.0),
        C::new(-3.0, -4.0),
        C::new(1e-13, 0.0),
        C::new(0.0, -1.0),
        C::new(f64::NAN, 0.0),
    ];
    for &t in &thresholds {
        for &v in &values {
            assert_eq!(
                Coefficient::magnitude(&v) < t,
                old_cutoff_complex(v, t),
                "magnitude/cutoff parity broke at v={v}, t={t}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// KeyProduct: the twisted (2-cocycle) product, pinned against the Lean oracle.
//
// `ppvm-traits-2` ships no concrete Pauli key, so the trait's *contract* is
// pinned on a single-qubit stub key whose phase is the packed boolean formula
// of `crates/ppvm-pauli-word/src/phase/mul.rs` — the same formula Lean names
// `PPVM.PauliPhase.phaseExp` (`2·sign + imag`). The tests check the three
// machine-checked claims of `lean/PPVM/Pauli/Phase.lean` in *trait terms*:
//
//   * `phaseExp_eq_ref`   — the booleans equal the analytic matrix exponent
//                           `ab + cd + 2bc − (a⊕c)(b⊕d)` (and, below, a genuine
//                           ℤ[i] 2×2 matrix product, `PauliMatrix.pauliMat_mul`);
//   * `phaseExp_cocycle`  — the residual phases `key_mul` emits compose (via the
//                           `Phase` group) into an associative twisted product;
//   * `phaseExp_sub_comm` — the phase asymmetry is `i^{2ω(p,q)}`, i.e.
//                           `P·Q = (−1)^{ω} Q·P`.
// ---------------------------------------------------------------------------

/// A single-qubit Pauli key `g(x, z) = i^{xz} Xˣ Zᶻ` — a point of GF(2)².
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PauliKey {
    x: bool,
    z: bool,
}

impl PauliKey {
    const fn new(x: bool, z: bool) -> Self {
        PauliKey { x, z }
    }
}

/// `sign` bit of the per-qubit product phase (Lean `PPVM.PauliPhase.signBit`).
///
/// Written verbatim as in `phase/mul.rs` / the Lean `signBit`, so the oracle is
/// syntactically comparable to the spec; clippy's "minimal" rewrite would break
/// that correspondence.
#[allow(clippy::nonminimal_bool)]
fn sign_bit(a: bool, b: bool, c: bool, d: bool) -> bool {
    (a && b && c && !d) || (a && !b && !c && d) || (!a && b && c && d)
}

/// `imag` bit of the per-qubit product phase (Lean `PPVM.PauliPhase.imagBit`).
///
/// Verbatim as in the Lean `imagBit`; see the note on [`sign_bit`].
#[allow(clippy::nonminimal_bool)]
fn imag_bit(a: bool, b: bool, c: bool, d: bool) -> bool {
    (a && !b && d) || (a && !c && d) || (!a && b && c) || (b && c && !d)
}

impl KeyProduct for PauliKey {
    fn key_mul(&self, other: &Self) -> (Self, Phase) {
        let (a, b, c, d) = (self.x, self.z, other.x, other.z);
        // Sp part: vector addition on GF(2)² (Lean `mulBits`).
        let key = PauliKey::new(a ^ c, b ^ d);
        // Extension part: 2·sign + imag (Lean `phaseExp`).
        let k = 2 * u8::from(sign_bit(a, b, c, d)) + u8::from(imag_bit(a, b, c, d));
        (key, Phase::from_exponent(k))
    }
}

const ALL_KEYS: [PauliKey; 4] = [
    PauliKey::new(false, false), // I
    PauliKey::new(true, false),  // X
    PauliKey::new(true, true),   // Y
    PauliKey::new(false, true),  // Z
];

/// Lean `phaseRef`: `ab + cd + 2bc − (a⊕c)(b⊕d)` in ℤ/4.
fn phase_ref(a: bool, b: bool, c: bool, d: bool) -> u8 {
    let n = |v: bool| i32::from(v);
    let e = n(a) * n(b) + n(c) * n(d) + 2 * (n(b) * n(c)) - n(a ^ c) * n(b ^ d);
    e.rem_euclid(4) as u8
}

#[test]
fn key_product_bits_and_phase_match_lean_reference() {
    for p in ALL_KEYS {
        for q in ALL_KEYS {
            let (k, phase) = p.key_mul(&q);
            // Sp part: ⊕ on the symplectic bits.
            assert_eq!(k, PauliKey::new(p.x ^ q.x, p.z ^ q.z));
            // Extension part: packed booleans == the analytic reference.
            assert_eq!(phase.exponent(), phase_ref(p.x, p.z, q.x, q.z));
        }
    }
    // Spot-checks against the textbook table: X·Z = −i·Y, Z·X = +i·Y, P·P = I.
    let (x, y, z) = (ALL_KEYS[1], ALL_KEYS[2], ALL_KEYS[3]);
    assert_eq!(x.key_mul(&z), (y, Phase::NegI));
    assert_eq!(z.key_mul(&x), (y, Phase::PosI));
    for p in ALL_KEYS {
        assert_eq!(p.key_mul(&p), (ALL_KEYS[0], Phase::Pos1));
    }
}

#[test]
fn key_product_phase_is_a_2_cocycle() {
    // Lean `phaseExp_cocycle`: the residual phases compose associatively, which
    // is exactly what makes the twisted product `C[K]` well defined.
    for p in ALL_KEYS {
        for q in ALL_KEYS {
            for r in ALL_KEYS {
                let (pq, a1) = p.key_mul(&q);
                let (pq_r, a2) = pq.key_mul(&r);

                let (qr, b1) = q.key_mul(&r);
                let (p_qr, b2) = p.key_mul(&qr);

                assert_eq!(pq_r, p_qr);
                // Accumulated with the first-class `Phase` group op.
                assert_eq!(a1.compose(a2), b1.compose(b2));
            }
        }
    }
}

#[test]
fn key_product_phase_asymmetry_is_the_symplectic_form() {
    // Lean `phaseExp_sub_comm`: phaseExp p q − phaseExp q p = 2·ω(p, q).
    for p in ALL_KEYS {
        for q in ALL_KEYS {
            let (_, pq) = p.key_mul(&q);
            let (_, qp) = q.key_mul(&p);
            let omega = (p.x && q.z) ^ (p.z && q.x);
            assert_eq!(
                pq.compose(qp.inverse()),
                Phase::from_exponent(2 * u8::from(omega))
            );
            // Equivalently: P·Q = (−1)^ω Q·P.
            if omega {
                assert_eq!(pq, qp.compose(Phase::Neg1));
            } else {
                assert_eq!(pq, qp);
            }
        }
    }
}

// --- The ℤ[i] matrix model (Lean `PPVM.PauliMatrix.pauliMat_mul`) -----------

/// A 2×2 matrix over ℤ[i] — exact, no floats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mat2([[Zi; 2]; 2]);

impl Mat2 {
    fn mul(self, rhs: Mat2) -> Mat2 {
        let mut out = [[<Zi as num::Zero>::zero(); 2]; 2];
        for (i, row) in out.iter_mut().enumerate() {
            for (j, slot) in row.iter_mut().enumerate() {
                *slot = self.0[i][0] * rhs.0[0][j] + self.0[i][1] * rhs.0[1][j];
            }
        }
        Mat2(out)
    }

    fn scale(self, s: Zi) -> Mat2 {
        Mat2([
            [self.0[0][0] * s, self.0[0][1] * s],
            [self.0[1][0] * s, self.0[1][1] * s],
        ])
    }
}

/// `g(x, z) = i^{xz} Xˣ Zᶻ` as an exact ℤ[i] matrix (the Lean normalization).
fn pauli_matrix(k: PauliKey) -> Mat2 {
    let zero = Zi::new(0, 0);
    let one = Zi::new(1, 0);
    let i = Zi::new(0, 1);
    match (k.x, k.z) {
        (false, false) => Mat2([[one, zero], [zero, one]]), // I
        (true, false) => Mat2([[zero, one], [one, zero]]),  // X
        (false, true) => Mat2([[one, zero], [zero, -one]]), // Z
        (true, true) => Mat2([[zero, -i], [i, zero]]),      // Y = i·X·Z
    }
}

#[test]
fn key_product_matches_the_gaussian_integer_matrix_model() {
    // `g(p)·g(q) == iᵏ · g(p·q)` with `k` the phase `key_mul` emits — the
    // grounding of the boolean formula in a genuine matrix product.
    for p in ALL_KEYS {
        for q in ALL_KEYS {
            let (k, phase) = p.key_mul(&q);
            let lhs = pauli_matrix(p).mul(pauli_matrix(q));
            let rhs = pauli_matrix(k).scale(phase.apply(&<Zi as num::One>::one()));
            assert_eq!(lhs, rhs, "matrix model disagrees for {p:?}·{q:?}");
        }
    }
}

// ---------------------------------------------------------------------------
// L4 `Multiply` over the exact ring: the twisted product `C[K]` is expressible
// with exactly the bounds the trait declares (`Key: KeyProduct`,
// `Coeff: ImaginaryUnit`) and is associative — `lean/PPVM/Algebra/Twisted.lean`
// (`tmul_assoc`) at the map level, and `multiply_single` in
// `lean/PPVM/Algebra/GradedMap.lean` for the basis monomials.
//
// The crate ships no L4 impl (`containers.rs` stops at L3 + `Retain`), so this
// instantiates the whole L0/L1/L4 stack on a local newtype over the shipped
// `Vec` layout, keyed by the stub Pauli key with exact ℤ[i] coefficients: a
// compile-level check that the L4 bounds are *usable* as declared (no extra
// bound is silently needed, and `ImaginaryUnit` really is enough — no
// `Complex<f64>`) plus a value-level check of the algebra laws, float-free.
//
// The newtype is needed only because a test is a separate crate and the orphan
// rule forbids `impl Multiply for Vec<_>` outside `ppvm-traits-2` — the same
// constraint `containers.rs` documents for L0–L3.
// ---------------------------------------------------------------------------

use ppvm_traits_2::batch::TermBatch;
use ppvm_traits_2::graded::{Accumulate, Multiply, Support};

/// `ℤ[i][PauliKey]` on the coordinate-list layout.
#[derive(Debug, Default, Clone, PartialEq)]
struct ZiSum(Vec<(PauliKey, Zi)>);

impl Support for ZiSum {
    type Key = PauliKey;
    type Coeff = Zi;

    fn len(&self) -> usize {
        self.0.len()
    }

    fn get(&self, key: &PauliKey) -> Option<Zi> {
        // `as_slice()` on purpose: `Support::iter` is in scope for `Vec`, and a
        // bare `self.0.iter()` would resolve to *it* rather than the slice
        // iterator (the same reason `containers.rs` spells `as_slice()`).
        self.0
            .as_slice()
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, c)| *c)
    }

    fn iter(&self) -> impl Iterator<Item = (PauliKey, Zi)> {
        self.0.as_slice().iter().copied()
    }
}

impl Accumulate for ZiSum {
    fn accumulate_batch(&mut self, terms: &TermBatch<PauliKey, Zi>) {
        for (k, c) in terms.iter() {
            if let Some(slot) = self.0.as_mut_slice().iter_mut().find(|(ek, _)| ek == k) {
                slot.1 += *c;
            } else {
                self.0.push((*k, *c));
            }
        }
    }

    fn reduce(&mut self) {
        self.0.retain(|(_, c)| !num::Zero::is_zero(c));
    }
}

impl Multiply for ZiSum {
    fn multiply_into(&self, other: &Self, acc: &mut Self) {
        for (p, a) in self.0.as_slice().iter() {
            for (q, b) in other.0.as_slice().iter() {
                let (key, phase) = p.key_mul(q);
                // The residual phase is absorbed by the coefficient (`iPow`).
                acc.accumulate(key, phase.apply(&(*a * *b)));
            }
        }
    }
}

fn zi_sum(terms: &[(PauliKey, Zi)]) -> ZiSum {
    let mut v = ZiSum::default();
    for (k, c) in terms {
        v.accumulate(*k, *c);
    }
    v
}

fn zi_mul(a: &ZiSum, b: &ZiSum) -> ZiSum {
    let mut acc = ZiSum::default();
    a.multiply_into(b, &mut acc);
    acc.reduce();
    acc
}

/// Canonical (sorted) view, so two supports compare independent of term order.
fn sorted(v: &ZiSum) -> Vec<((bool, bool), (i64, i64))> {
    let mut out: Vec<((bool, bool), (i64, i64))> = Support::iter(v)
        .map(|(k, c)| ((k.x, k.z), (c.re, c.im)))
        .collect();
    out.sort();
    out
}

#[test]
fn multiply_realizes_the_twisted_product_on_basis_monomials() {
    let one = Zi::new(1, 0);
    let (id, x, y, z) = (ALL_KEYS[0], ALL_KEYS[1], ALL_KEYS[2], ALL_KEYS[3]);

    // X·Z = −i·Y (`multiply_single`: the product of two basis monomials is the
    // single term `phase · (k₁·k₂)`).
    let prod = zi_mul(&zi_sum(&[(x, one)]), &zi_sum(&[(z, one)]));
    assert_eq!(sorted(&prod), sorted(&zi_sum(&[(y, Zi::new(0, -1))])));

    // (X + Z)·(X + Z) = 2·I — the ∓i·Y cross terms cancel and `reduce` drops
    // them (never inline: before `reduce` the zero key is still supported).
    let s = zi_sum(&[(x, one), (z, one)]);
    let mut raw = ZiSum::default();
    s.multiply_into(&s, &mut raw);
    assert_eq!(Support::len(&raw), 2); // I and the cancelled Y
    assert_eq!(Support::get(&raw, &y), Some(Zi::new(0, 0)));
    let sq = zi_mul(&s, &s);
    assert_eq!(sorted(&sq), sorted(&zi_sum(&[(id, Zi::new(2, 0))])));
}

// ---------------------------------------------------------------------------
// `Coefficient::magnitude` is an **absolute value**, not merely a nonnegative
// number. This is the standing impl obligation the trait doc states and the
// hypothesis the whole truncation guarantee consumes:
// `lean/PPVM/Algebra/Truncation.lean` (`l1_bound_abv` over any coefficient ring
// with such an `N`, specialized by `l1_bound_norm` / `l1_bound_complex`).
// `l1_bound_needs_subadditive` there shows nonnegativity alone does NOT suffice
// (`N(x) = x²` is nonnegative, vanishes only at 0, and is multiplicative, yet
// breaks the ℓ¹ bound) — so subadditivity is checked explicitly below, on every
// shipped witness plus the exact ℤ[i] one.
// ---------------------------------------------------------------------------

/// Check the four absolute-value laws on a sample of ring elements.
fn assert_absolute_value<T: Coefficient>(samples: &[T]) {
    let eps = 1e-9;
    for a in samples {
        let n = a.magnitude();
        // N(x) >= 0.
        assert!(n >= 0.0, "magnitude went negative");
        // N(x) == 0 iff x == 0.
        assert_eq!(n == 0.0, num::Zero::is_zero(a), "vanishing set is wrong");
        for b in samples {
            // Subadditive: N(x + y) <= N(x) + N(y).
            let sum = (a.clone() + b.clone()).magnitude();
            assert!(
                sum <= n + b.magnitude() + eps,
                "subadditivity broke: {sum} > {n} + {}",
                b.magnitude()
            );
            // Multiplicative: N(x * y) == N(x) * N(y).
            let prod = (a.clone() * b.clone()).magnitude();
            assert!(
                (prod - n * b.magnitude()).abs() <= eps * (1.0 + prod),
                "multiplicativity broke: {prod} != {n} * {}",
                b.magnitude()
            );
        }
    }
}

#[test]
fn magnitude_is_an_absolute_value_on_every_shipped_coefficient() {
    assert_absolute_value(&[0.0_f64, 1.0, -1.0, 0.5, -3.25, 7.0]);
    assert_absolute_value(&[
        C::new(0.0, 0.0),
        C::new(1.0, 0.0),
        C::new(0.0, -1.0),
        C::new(3.0, 4.0),
        C::new(-0.5, 0.25),
    ]);
    // The exact ring too: `magnitude` is the ℤ[i] norm √(a²+b²), still an
    // absolute value, so ℤ[i] coefficients get the same ℓ¹ truncation bound.
    assert_absolute_value(&[
        Zi::new(0, 0),
        Zi::new(1, 0),
        Zi::new(0, -1),
        Zi::new(3, 4),
        Zi::new(-2, 1),
    ]);
}

#[test]
fn magnitude_rejects_the_squared_pseudo_norm() {
    // Guard on the guard: the Lean counterexample `l1_bound_needs_subadditive`
    // is that `N(x) = x²` satisfies every *other* clause. Pin that the shipped
    // `magnitude` is not that function — i.e. the test above has teeth.
    let squared = |x: f64| x * x;
    let (a, b) = (3.0_f64, 4.0_f64);
    assert!(squared(a + b) > squared(a) + squared(b)); // not subadditive
    assert!(
        Coefficient::magnitude(&(a + b)) <= Coefficient::magnitude(&a) + Coefficient::magnitude(&b)
    );
}

// ---------------------------------------------------------------------------
// The angle domain drives a norm-preserving, angle-additive 2-D rotation —
// `lean/PPVM/Instantiations/Rotation.lean` (`rot_norm_sq`, `rot_rot`). The
// crate ships no rotation kernel, so the laws are checked on the `(sin, cos)`
// the shipped [`Angle`] impls hand it: the pair is a point of the unit circle
// and composing two rotations adds the angles.
// ---------------------------------------------------------------------------

/// Apply the 2-D rotation `[[cos, −sin], [sin, cos]]` built from an angle.
fn rot(theta: f64, v: (f64, f64)) -> (f64, f64) {
    let (s, c) = Angle::<f64>::sin_cos(&theta);
    (c * v.0 - s * v.1, s * v.0 + c * v.1)
}

#[test]
fn angle_sin_cos_is_norm_preserving_and_angle_additive() {
    use std::f64::consts::PI;
    let eps = 1e-12;
    let angles = [0.0_f64, 0.25, -0.75, PI / 6.0, PI / 3.0, PI, 2.5];

    for &t in &angles {
        // `rot_norm_sq`: sin²θ + cos²θ == 1, so the rotation is an isometry.
        let (s, c) = Angle::<f64>::sin_cos(&t);
        assert!((s * s + c * c - 1.0).abs() < eps, "not on the unit circle");
        let v = (2.0_f64, -3.0);
        let w = rot(t, v);
        assert!(
            ((w.0 * w.0 + w.1 * w.1) - (v.0 * v.0 + v.1 * v.1)).abs() < 1e-10,
            "rotation changed the norm"
        );

        // `rot_rot`: rot(t2)∘rot(t1) == rot(t1 + t2).
        for &u in &angles {
            let composed = rot(u, rot(t, v));
            let direct = rot(t + u, v);
            assert!(
                (composed.0 - direct.0).abs() < 1e-10 && (composed.1 - direct.1).abs() < 1e-10,
                "angle additivity broke at ({t}, {u})"
            );
        }
    }
}

#[test]
fn complex_angle_sin_cos_is_also_on_the_unit_circle() {
    // The `A = C` instantiation returns purely real amplitudes taken from the
    // real part, so the same `rot_norm_sq` law holds there (and the imaginary
    // part of theta is inert, as the old `Coefficient::sin_cos` had it).
    for &(re, im) in &[(0.0_f64, 0.0_f64), (0.6, 3.0), (-2.2, -1.0)] {
        let (s, c) = Angle::<C>::sin_cos(&C::new(re, im));
        assert_eq!(s.im, 0.0);
        assert_eq!(c.im, 0.0);
        assert!((s.re * s.re + c.re * c.re - 1.0).abs() < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// `Phase::apply` is a ℤ/4 *group action* on the coefficient ring, not just a
// per-element table: this is the `iPow` fold of
// `lean/PPVM/Algebra/Twisted.lean`, and it is what lets a `KeyProduct` chain
// accumulate residual phases with `compose` and fold onto the coefficient only
// once at the end (the deferral the `Phase::compose` doc promises).
// ---------------------------------------------------------------------------

#[test]
fn phase_apply_is_a_group_action() {
    // Over Complex<f64>...
    for &c in &[C::new(0.0, 0.0), C::new(2.0, -3.0), C::new(-1.5, 0.25)] {
        assert_eq!(Phase::one().apply(&c), c);
        for a in ALL_PHASES {
            for b in ALL_PHASES {
                // Deferred fold (compose then apply) == eager fold (apply twice).
                assert_eq!(a.compose(b).apply(&c), a.apply(&b.apply(&c)));
            }
            // The inverse phase undoes the fold.
            assert_eq!(a.inverse().apply(&a.apply(&c)), c);
        }
    }
    // ...and exactly, over ℤ[i], with no floating point at all.
    for &z in &[Zi::new(0, 0), Zi::new(3, -2), Zi::new(-1, 5)] {
        assert_eq!(Phase::one().apply(&z), z);
        for a in ALL_PHASES {
            for b in ALL_PHASES {
                assert_eq!(a.compose(b).apply(&z), a.apply(&b.apply(&z)));
            }
            assert_eq!(a.inverse().apply(&a.apply(&z)), z);
        }
    }
}

// ---------------------------------------------------------------------------
// The pass-through storage contract end to end: an `Indexable` key's `Hash`
// impl writes exactly `write_u64(key_hash())`, so a `BuildHasher`-driven hash
// of the key — the path hashbrown actually takes — equals `key_hash()` verbatim.
// The unit test above only pins `write_u64` → `finish`; this pins that no byte
// is added by the `Hash`/`BuildHasher` plumbing in between.
//
// Design: §"The pass-through storage contract".
// ---------------------------------------------------------------------------

/// An `Indexable` key with an avalanche-quality digest distinct from its body,
/// so a pass-through failure cannot be masked by the two happening to agree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DigestKey(u64);

impl std::hash::Hash for DigestKey {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // The contract: exactly one `write_u64`, of the finalized digest.
        state.write_u64(Indexable::key_hash(self));
    }
}

impl Indexable for DigestKey {
    fn key_hash(&self) -> u64 {
        // A cheap avalanche so low bits and the top-7 control tag both move.
        let mut h = self.0 ^ 0x9e37_79b9_7f4a_7c15;
        h ^= h >> 33;
        h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h ^= h >> 29;
        h
    }
}

#[test]
fn build_hasher_hash_of_an_indexable_key_is_its_key_hash() {
    let bh = IdentityBuildHasher;
    for raw in [0u64, 1, 7, 0xdead_beef, u64::MAX] {
        let k = DigestKey(raw);
        // `hash_one` is the exact path hashbrown drives.
        assert_eq!(bh.hash_one(k), Indexable::key_hash(&k));
        // And the manual spelling agrees.
        let mut h = bh.build_hasher();
        std::hash::Hash::hash(&k, &mut h);
        assert_eq!(h.finish(), Indexable::key_hash(&k));
        // Equal keys ⇒ equal digests (the second `Indexable` contract clause).
        assert_eq!(
            Indexable::key_hash(&k),
            Indexable::key_hash(&DigestKey(raw))
        );
    }
}

// ---------------------------------------------------------------------------
// `PauliBits::is_lost` is the crate's one remaining provided body outside the
// gate traits, and its default is load-bearing: an *ordinary* (non-lossy) word
// must report every site present, so the loss guard the lossy word adds inside
// its column primitives is a pure extension and never fires on `PauliWord`.
// A stub over the shipped `Pauli` alphabet pins the default (and `Word`'s own
// `n_sites`/`get`/`weight`/`iter` reads it is defined against).
// ---------------------------------------------------------------------------

struct BitsStub {
    x: Vec<bool>,
    z: Vec<bool>,
}

impl Word for BitsStub {
    type Site = Pauli;

    fn n_sites(&self) -> usize {
        self.x.len()
    }

    fn get(&self, index: usize) -> Pauli {
        match (self.x[index], self.z[index]) {
            (false, false) => Pauli::I,
            (true, false) => Pauli::X,
            (true, true) => Pauli::Y,
            (false, true) => Pauli::Z,
        }
    }

    fn weight(&self) -> usize {
        (0..self.n_sites())
            .filter(|&i| self.get(i) != Pauli::I)
            .count()
    }

    fn iter(&self) -> impl Iterator<Item = Pauli> {
        (0..self.n_sites()).map(|i| self.get(i))
    }
}

impl PauliBits for BitsStub {
    fn x_bit(&self, i: usize) -> bool {
        self.x[i]
    }
    fn z_bit(&self, i: usize) -> bool {
        self.z[i]
    }
    fn set_x_bit(&mut self, i: usize, v: bool) {
        self.x[i] = v;
    }
    fn set_z_bit(&mut self, i: usize, v: bool) {
        self.z[i] = v;
    }
    // NOTE: `is_lost` deliberately not overridden — the default is under test.
}

#[test]
fn pauli_bits_is_lost_defaults_to_present_everywhere() {
    let mut w = BitsStub {
        x: vec![false, true, true],
        z: vec![true, false, true],
    };
    // The `Word` reads agree with the X/Z bits.
    assert_eq!(w.n_sites(), 3);
    assert_eq!(
        w.iter().collect::<Vec<_>>(),
        vec![Pauli::Z, Pauli::X, Pauli::Y]
    );
    assert_eq!(w.weight(), 3);

    // The provided `is_lost` reports every site present, including out of the
    // written range — a non-lossy word has no loss plane to consult.
    for i in 0..5 {
        assert!(!PauliBits::is_lost(&w, i));
    }

    // Bit mutation does not make a site lost either.
    w.set_x_bit(0, true);
    w.set_z_bit(0, false);
    assert_eq!(w.get(0), Pauli::X);
    assert_eq!(w.weight(), 3);
    assert!(!PauliBits::is_lost(&w, 0));

    // Clearing both bits is the identity site, still present, weight drops.
    w.set_x_bit(0, false);
    assert_eq!(w.get(0), Pauli::I);
    assert_eq!(w.weight(), 2);
    assert!(!PauliBits::is_lost(&w, 0));
}

// ---------------------------------------------------------------------------
// Coherence guard for the `BlanketClifford` opt-in marker.
//
// `pauli.rs` gates the blanket `impl Clifford` on the empty `BlanketClifford`
// marker precisely so a type that wants a *fused* single-pass `Clifford` can
// write its own impl without tripping E0119. That guarantee is a compile-time
// property, so this stub is the regression test: it implements
// `SymplecticColumns + PhaseTrack` (which would satisfy an *unconditional*
// blanket), does NOT opt in, and supplies its own `Clifford`. If the marker
// gate were ever dropped from the blanket, this file would stop compiling —
// which is exactly the failure `Phased<W>` in `ppvm-phased-pauli-word-2` would
// hit.
// ---------------------------------------------------------------------------

#[derive(Default)]
struct FusedStub {
    log: Vec<String>,
}

impl SymplecticColumns for FusedStub {
    fn n_qubits(&self) -> usize {
        1
    }
    fn swap_xz(&mut self, _q: usize) {}
    fn xor_z_from_x(&mut self, _q: usize) {}
    fn xor_x_col(&mut self, _ctrl: usize, _tgt: usize) {}
    fn xor_z_col(&mut self, _tgt: usize, _ctrl: usize) {}
    fn cz_bits(&mut self, _a: usize, _b: usize) {}
}

impl PhaseTrack for FusedStub {
    fn flip_phase_where_xz(&mut self, _q: usize) {}
    fn s_phase(&mut self, _q: usize) {}
    fn cnot_phase(&mut self, _ctrl: usize, _tgt: usize) {}
    fn cz_phase(&mut self, _a: usize, _b: usize) {}
    fn x_phase(&mut self, _q: usize) {}
    fn y_phase(&mut self, _q: usize) {}
    fn z_phase(&mut self, _q: usize) {}
}

// Deliberately NOT `impl BlanketClifford for FusedStub {}` — this type supplies
// its own fused `Clifford` instead.
impl Clifford for FusedStub {
    fn x(&mut self, q: usize) {
        self.log.push(format!("fused_x({q})"));
    }
    fn y(&mut self, q: usize) {
        self.log.push(format!("fused_y({q})"));
    }
    fn z(&mut self, q: usize) {
        self.log.push(format!("fused_z({q})"));
    }
    fn h(&mut self, q: usize) {
        self.log.push(format!("fused_h({q})"));
    }
    fn s(&mut self, q: usize) {
        self.log.push(format!("fused_s({q})"));
    }
    fn cnot(&mut self, c: usize, t: usize) {
        self.log.push(format!("fused_cnot({c},{t})"));
    }
    fn cz(&mut self, a: usize, b: usize) {
        self.log.push(format!("fused_cz({a},{b})"));
    }
}

#[test]
fn a_non_opting_in_type_may_supply_its_own_fused_clifford() {
    let mut f = FusedStub::default();
    requires_clifford(&f);

    // The fused impl runs — the blanket does not shadow or duplicate it, and
    // `h` fires exactly one fused step rather than the two blanket primitives.
    f.h(0);
    f.s(1);
    f.cnot(0, 1);
    f.cz(1, 0);
    assert_eq!(
        f.log,
        vec![
            "fused_h(0)",
            "fused_s(1)",
            "fused_cnot(0,1)",
            "fused_cz(1,0)"
        ]
    );

    // The defaulted stim aliases still route through the *fused* required
    // methods (they are `Clifford` defaults, not blanket-specific).
    f.log.clear();
    f.cx(2, 3);
    f.zcz(3, 2);
    assert_eq!(f.log, vec!["fused_cnot(2,3)", "fused_cz(3,2)"]);
}

#[test]
fn multiply_is_associative_over_the_exact_ring() {
    // `tmul_assoc`: associativity of the twisted product, i.e. the 2-cocycle
    // property of the phase lifted to the whole algebra.
    let a = zi_sum(&[(ALL_KEYS[1], Zi::new(2, 0)), (ALL_KEYS[3], Zi::new(0, 1))]);
    let b = zi_sum(&[(ALL_KEYS[2], Zi::new(1, -1)), (ALL_KEYS[0], Zi::new(3, 0))]);
    let c = zi_sum(&[(ALL_KEYS[3], Zi::new(-1, 2)), (ALL_KEYS[1], Zi::new(1, 1))]);

    let left = zi_mul(&zi_mul(&a, &b), &c);
    let right = zi_mul(&a, &zi_mul(&b, &c));
    assert_eq!(sorted(&left), sorted(&right));
    assert!(!Support::is_empty(&left));
}
