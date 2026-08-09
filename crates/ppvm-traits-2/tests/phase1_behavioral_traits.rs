// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Tests pinning the *contracts* of the behavioral gate/noise traits
//! ([`RotationOne`], [`PauliError`], [`Measure`]) on stubs, since the concrete
//! implementors land downstream (`ppvm-pauli-sum-2`, `ppvm-tableau-2`).
//!
//! What is genuinely executable in this crate — the trait *defaults* and the
//! angle-domain plumbing — is pinned here:
//!
//!   * the **defaulted angle domain** `RotationOne<C, A = C>` — a stub must be
//!     instantiable both at the old `rx(theta: C)` shape (`A = C`, the
//!     behaviour-parity case) *and* at a symbolic angle `A ≠ C`, the capability
//!     the `Coefficient::sin_cos` → [`Angle`] split buys (design §"Coefficient,
//!     angle, and truncation"), *and* at the old caller spelling
//!     `sum.rx(0, 0.1)` on a complex-coefficient sum (`A = f64`, `C =
//!     Complex<f64>`), which the old `impl Into<T::Coeff>` used to provide;
//!   * the nine named [`RotationTwo`] gates, whose defaults each carry a literal
//!     `[x, z]` axis pair transcribed from the old `def_rotation!`;
//!   * the [`TGate`] batch loops; and
//!   * the **provided** [`Measure::measure_many`] loop, copied verbatim from
//!     the old `LossyMeasure::measure_many`: one result per target, in order,
//!     `None` for a lost qubit. Sharing the `Option<bool>` result type replaced
//!     the old `Measure -> bool` / `LossyMeasure -> Option<bool>` split.

use std::f64::consts::TAU;

use ppvm_traits_2::coefficient::Angle;
use ppvm_traits_2::gates::{Measure, PauliError, RotationOne, RotationTwo, TGate};
use ppvm_traits_2::word::Pauli;
use rand::SeedableRng;
use rand::rngs::SmallRng;

// ---------------------------------------------------------------------------
// RotationOne: defaulted angle (A = C) and a symbolic angle domain (A != C).
// ---------------------------------------------------------------------------

/// A symbolic angle measured in *turns*, over `f64` coefficients — an `A` that
/// is deliberately not the coefficient type.
#[derive(Debug, Clone, Copy)]
struct Turns(f64);

impl Angle<f64> for Turns {
    fn sin_cos(&self) -> (f64, f64) {
        (self.0 * TAU).sin_cos()
    }
}

#[derive(Default)]
struct StubSum {
    /// `(axis, qubit, sin, cos)` per rotation, in call order.
    log: Vec<(char, usize, f64, f64)>,
    /// `(qubit, probabilities)` per noise call.
    noise: Vec<(usize, [f64; 3])>,
}

/// The defaulted case `A = C`: exactly today's `rx(theta: C)`.
///
/// Only the *required* [`RotationOne::rotate_1`] is written, so `rx`/`ry`/`rz`
/// and the `*_many` batch forms below exercise the old crate's provided bodies
/// (`ppvm-traits/src/traits/branch/rot1.rs`).
impl RotationOne<f64> for StubSum {
    fn rotate_1(&mut self, axis: Pauli, qubit: usize, theta: f64) {
        let (s, c) = Angle::<f64>::sin_cos(&theta);
        let axis = match axis {
            Pauli::I => 'i',
            Pauli::X => 'x',
            Pauli::Y => 'y',
            Pauli::Z => 'z',
        };
        self.log.push((axis, qubit, s, c));
    }
}

/// The same coefficient domain driven by a *different* angle domain — and the
/// override point: `rx`/`ry`/`rz` are supplied directly, as `PauliSum` does with
/// its per-axis fast paths, so `rotate_1` is never consulted for those axes.
impl RotationOne<f64, Turns> for StubSum {
    fn rotate_1(&mut self, _axis: Pauli, _qubit: usize, _theta: Turns) {
        unreachable!("the per-axis overrides below must be preferred");
    }
    fn rx(&mut self, qubit: usize, theta: Turns) {
        let (s, c) = theta.sin_cos();
        self.log.push(('x', qubit, s, c));
    }
    fn ry(&mut self, qubit: usize, theta: Turns) {
        let (s, c) = theta.sin_cos();
        self.log.push(('y', qubit, s, c));
    }
    fn rz(&mut self, qubit: usize, theta: Turns) {
        let (s, c) = theta.sin_cos();
        self.log.push(('z', qubit, s, c));
    }
}

impl PauliError<f64> for StubSum {
    fn pauli_error<R: rand::Rng + ?Sized>(
        &mut self,
        qubit: usize,
        probabilities: [f64; 3],
        _rng: &mut R,
    ) {
        self.noise.push((qubit, probabilities));
    }
}

#[test]
fn rotation_one_defaults_its_angle_to_the_coefficient() {
    let mut s = StubSum::default();
    // No turbofish needed: `A` defaults to `C = f64`, i.e. the old `rx(theta: C)`.
    s.rx(0, 0.25);
    s.ry(1, -1.0);
    s.rz(2, 0.0);

    let (s0, c0) = 0.25_f64.sin_cos();
    let (s1, c1) = (-1.0_f64).sin_cos();
    assert_eq!(
        s.log,
        vec![('x', 0, s0, c0), ('y', 1, s1, c1), ('z', 2, 0.0, 1.0)]
    );
}

#[test]
fn rotation_one_admits_a_symbolic_angle_domain() {
    let mut s = StubSum::default();
    // Half a turn = π: the angle domain converts, the coefficients stay f64.
    RotationOne::<f64, Turns>::rx(&mut s, 3, Turns(0.5));

    let (sin, cos) = std::f64::consts::PI.sin_cos();
    assert_eq!(s.log, vec![('x', 3, sin, cos)]);

    // A quarter turn is a genuine π/2 rotation (sin ≈ 1, cos ≈ 0).
    s.log.clear();
    RotationOne::<f64, Turns>::rz(&mut s, 0, Turns(0.25));
    let (sin, cos) = (s.log[0].2, s.log[0].3);
    assert!((sin - 1.0).abs() < 1e-12);
    assert!(cos.abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// Behaviour parity: the old `RotationOne` required `rotate_1(axis, …)` and
// defaulted `rx`/`ry`/`rz` onto it, then defaulted `rx_many`/`ry_many`/`rz_many`
// onto those (`ppvm-traits/src/traits/branch/rot1.rs`). Those defaults are the
// entry points `rotate_2` and the Python bindings call, so they are pinned here.
// ---------------------------------------------------------------------------

#[test]
fn rx_ry_rz_default_onto_rotate_1() {
    let mut s = StubSum::default();
    s.rx(0, 0.25);
    s.ry(1, 0.25);
    s.rz(2, 0.25);

    // The axis reaching `rotate_1` is exactly the old default's.
    let mut direct = StubSum::default();
    direct.rotate_1(Pauli::X, 0, 0.25);
    direct.rotate_1(Pauli::Y, 1, 0.25);
    direct.rotate_1(Pauli::Z, 2, 0.25);

    assert_eq!(s.log, direct.log);
}

#[test]
fn rotate_1_accepts_the_identity_axis() {
    // The old `rotate_1` panicked only on `Pauli::L` (no longer a `Pauli`); an
    // `I` axis was accepted and commuted with everything.
    let mut s = StubSum::default();
    s.rotate_1(Pauli::I, 7, 0.5);
    let (sin, cos) = 0.5_f64.sin_cos();
    assert_eq!(s.log, vec![('i', 7, sin, cos)]);
}

#[test]
fn rotation_batch_forms_loop_in_target_order() {
    let mut s = StubSum::default();
    s.rx_many(&[0, 2], 0.25);
    s.ry_many(&[1], -0.5);
    s.rz_many(&[3, 3], 1.0);

    let mut expect = StubSum::default();
    expect.rx(0, 0.25);
    expect.rx(2, 0.25);
    expect.ry(1, -0.5);
    expect.rz(3, 1.0);
    expect.rz(3, 1.0);

    assert_eq!(s.log, expect.log);
}

#[test]
fn rotation_batch_forms_reach_a_per_axis_override() {
    // `*_many` loops over `rx`/`ry`/`rz`, so an override wins there too — the
    // `Turns` impl would `unreachable!()` if the batch bypassed it.
    let mut s = StubSum::default();
    RotationOne::<f64, Turns>::rx_many(&mut s, &[4, 5], Turns(0.5));

    let (sin, cos) = std::f64::consts::PI.sin_cos();
    assert_eq!(s.log, vec![('x', 4, sin, cos), ('x', 5, sin, cos)]);
}

// ---------------------------------------------------------------------------
// Behaviour parity: the old caller spelling `sum.rx(0, 0.1)` on a
// **complex-coefficient** sum.
//
// The old methods took `theta: impl Into<T::Coeff>` and `Coefficient: From<f64>`,
// so a bare `f64` literal was accepted on a `Complex<f64>` sum and widened. `A`
// is a free trait parameter here, so `impl Into<A>` is not available; the
// spelling is preserved by `Angle<Complex<f64>> for f64` instead. This stub is
// the *call-site* witness — it implements `RotationOne<Complex<f64>, f64>` only,
// so `s.rx(0, 0.1)` must resolve with no turbofish, no `.into()`, and no type
// annotation, and must produce the amplitudes the old widening produced.
// ---------------------------------------------------------------------------

type C = num::Complex<f64>;

#[derive(Default)]
struct ComplexStubSum {
    /// `(axis, qubit, sin, cos)` per rotation, in call order.
    log: Vec<(char, usize, C, C)>,
}

impl RotationOne<C, f64> for ComplexStubSum {
    fn rotate_1(&mut self, axis: Pauli, qubit: usize, theta: f64) {
        let (s, c) = Angle::<C>::sin_cos(&theta);
        let axis = match axis {
            Pauli::I => 'i',
            Pauli::X => 'x',
            Pauli::Y => 'y',
            Pauli::Z => 'z',
        };
        self.log.push((axis, qubit, s, c));
    }
}

#[test]
fn a_real_f64_angle_still_drives_a_complex_coefficient_rotation() {
    let mut s = ComplexStubSum::default();
    // Exactly the old spelling: an `f64` literal on a `Complex<f64>` sum.
    s.rx(0, 0.1);
    s.ry(1, -0.25);
    s.rz(2, 0.0);
    // The batch defaults clone the angle (`A: Clone`), as the old `*_many` did
    // after converting once.
    s.rx_many(&[3, 4], 0.1);

    /// The old path: `theta.into()` (`Complex::from`) then `Coefficient::sin_cos`.
    fn old(theta: f64) -> (C, C) {
        let widened = C::from(theta);
        let (sin, cos) = widened.re.sin_cos();
        (C::new(sin, 0.0), C::new(cos, 0.0))
    }

    let (s0, c0) = old(0.1);
    let (s1, c1) = old(-0.25);
    let (s2, c2) = old(0.0);
    assert_eq!(
        s.log,
        vec![
            ('x', 0, s0, c0),
            ('y', 1, s1, c1),
            ('z', 2, s2, c2),
            ('x', 3, s0, c0),
            ('x', 4, s0, c0),
        ]
    );
}

#[test]
fn pauli_error_takes_three_probabilities_in_the_coefficient_domain() {
    let mut s = StubSum::default();
    let mut rng = SmallRng::seed_from_u64(0);
    s.pauli_error(2, [0.01, 0.02, 0.03], &mut rng);
    assert_eq!(s.noise, vec![(2, [0.01, 0.02, 0.03])]);
}

// ---------------------------------------------------------------------------
// RotationTwo: the nine named gates are *defaults* whose bodies carry a literal
// `[x, z]` axis pair each, transcribed from the old `def_rotation!` invocations
// (`ppvm-traits/src/traits/branch/rot2.rs`). A transposed pair (`rxy` vs `ryx`)
// would silently rotate about the wrong generator with no type error anywhere,
// so the whole table is pinned against the old constants, along with the batch
// loops (which the old crate spelled `impl Into<T::Coeff>` + one conversion, and
// which now clone an `A: Clone`).
//
// Encoding (unchanged): `[0,0]` = I, `[1,0]` = X, `[0,1]` = Z, `[1,1]` = Y.
// ---------------------------------------------------------------------------

/// One recorded `rotate_2` call: `(axis_a, axis_b, a, b, theta)`.
type Rot2Call = ([u8; 2], [u8; 2], usize, usize, f64);

#[derive(Default)]
struct Rot2Stub {
    /// Every `rotate_2` the named defaults emitted, in call order.
    log: Vec<Rot2Call>,
}

impl RotationTwo<f64> for Rot2Stub {
    fn rotate_2(&mut self, axis_a: [u8; 2], axis_b: [u8; 2], a: usize, b: usize, theta: f64) {
        self.log.push((axis_a, axis_b, a, b, theta));
    }
}

const X: [u8; 2] = [1, 0];
const Y: [u8; 2] = [1, 1];
const Z: [u8; 2] = [0, 1];

#[test]
fn named_two_qubit_rotations_carry_the_old_axis_encodings() {
    /// One named two-qubit rotation: its name, the call, and the axis pair the
    /// default must forward to `rotate_2`.
    type NamedRotation = (&'static str, fn(&mut Rot2Stub), [u8; 2], [u8; 2]);

    let cases: [NamedRotation; 9] = [
        ("rxx", |s| s.rxx(0, 1, 0.5), X, X),
        ("rxy", |s| s.rxy(0, 1, 0.5), X, Y),
        ("rxz", |s| s.rxz(0, 1, 0.5), X, Z),
        ("ryx", |s| s.ryx(0, 1, 0.5), Y, X),
        ("ryy", |s| s.ryy(0, 1, 0.5), Y, Y),
        ("ryz", |s| s.ryz(0, 1, 0.5), Y, Z),
        ("rzx", |s| s.rzx(0, 1, 0.5), Z, X),
        ("rzy", |s| s.rzy(0, 1, 0.5), Z, Y),
        ("rzz", |s| s.rzz(0, 1, 0.5), Z, Z),
    ];
    for (name, call, axis_a, axis_b) in cases {
        let mut s = Rot2Stub::default();
        call(&mut s);
        assert_eq!(s.log, vec![(axis_a, axis_b, 0, 1, 0.5)], "{name}");
    }
}

#[test]
fn two_qubit_rotation_batch_forms_loop_in_pair_order() {
    let mut s = Rot2Stub::default();
    s.rzz_many(&[(0, 1), (2, 3), (0, 1)], 0.25);

    let mut expect = Rot2Stub::default();
    expect.rzz(0, 1, 0.25);
    expect.rzz(2, 3, 0.25);
    expect.rzz(0, 1, 0.25);
    assert_eq!(s.log, expect.log);

    // An empty pair list is a no-op, as the old loop was.
    let mut s = Rot2Stub::default();
    s.rxy_many(&[], 0.25);
    assert!(s.log.is_empty());
}

// ---------------------------------------------------------------------------
// TGate: `t_many` / `t_dag_many` are the old crate's loop defaults
// (`ppvm-traits/src/traits/branch/tgate.rs`), and `T` takes no numeric argument,
// so the new trait carries no coefficient parameter at all.
// ---------------------------------------------------------------------------

#[derive(Default)]
struct TStub {
    log: Vec<String>,
}

impl TGate for TStub {
    fn t(&mut self, qubit: usize) {
        self.log.push(format!("t({qubit})"));
    }
    fn t_dag(&mut self, qubit: usize) {
        self.log.push(format!("t_dag({qubit})"));
    }
}

#[test]
fn t_gate_batch_defaults_loop_over_the_required_methods() {
    let mut s = TStub::default();
    s.t_many(&[0, 2, 0]);
    s.t_dag_many(&[1]);
    s.t_many(&[]);
    assert_eq!(s.log, vec!["t(0)", "t(2)", "t(0)", "t_dag(1)"]);
}

// ---------------------------------------------------------------------------
// Measure::measure_many — the provided loop, ported verbatim from the old
// `LossyMeasure::measure_many`.
// ---------------------------------------------------------------------------

struct StubMeasure {
    /// Qubits that have been lost; measuring one yields `None`.
    lost: Vec<usize>,
    /// Every qubit measured, in order (measurement is destructive, so the
    /// *order* and *count* of side effects is part of the contract).
    seen: Vec<usize>,
}

impl Measure for StubMeasure {
    fn measure<R: rand::Rng + ?Sized>(&mut self, qubit: usize, _rng: &mut R) -> Option<bool> {
        self.seen.push(qubit);
        if self.lost.contains(&qubit) {
            None
        } else {
            Some(qubit % 2 == 1)
        }
    }
}

#[test]
fn measure_many_returns_one_result_per_target_in_order() {
    let mut m = StubMeasure {
        lost: vec![2],
        seen: Vec::new(),
    };
    let mut rng = SmallRng::seed_from_u64(0);
    let out = m.measure_many(&[0, 1, 2, 3, 1], &mut rng);
    assert_eq!(
        out,
        vec![Some(false), Some(true), None, Some(true), Some(true)]
    );
    // One `measure` call per target, in the given order — repeats included.
    assert_eq!(m.seen, vec![0, 1, 2, 3, 1]);

    // Empty target list measures nothing.
    m.seen.clear();
    assert!(m.measure_many(&[], &mut rng).is_empty());
    assert!(m.seen.is_empty());
}
