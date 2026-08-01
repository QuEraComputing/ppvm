// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The shared behavioral gate/noise traits: [`Clifford`] (+ [`CliffordExtensions`]
//! and the batched forms), the rotation family ([`RotationOne`], [`RotationTwo`],
//! [`RotXY`], [`CRx`], [`U3Gate`]), [`TGate`], [`Projection`], [`Reset`],
//! [`Measure`], and the noise channel family headed by [`PauliError`].
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Behavioral traits". These
//! describe operations, not representation layout; Clifford gates need no
//! coefficient parameter, numeric operations take the coefficient type directly.
//!
//! # Behaviour parity with `ppvm-traits`
//!
//! Every trait, method name, defaulted stim alias, and default body below is the
//! old crate's (`ppvm-traits/src/traits/{clifford,reset,noise}.rs` and
//! `ppvm-traits/src/traits/branch/*.rs`) — the *only* differences are the ones
//! the design's compatibility table already records for the whole crate:
//!
//! * the global `Config` parameter is replaced by the coefficient type itself
//!   (`PauliError<T: Config>` → `PauliError<C: Coefficient>`, and likewise for
//!   every channel); [`ResetLossChannel`] consumes no coefficient at all, so it
//!   loses the parameter entirely, exactly as `Clifford` did;
//! * `addr0`/`addr1` argument names read `qubit`/`control`/`target`;
//! * a rotation angle is the angle domain `A` (defaulting to `C`) rather than
//!   `impl Into<T::Coeff>`, and the batch defaults that clone it say so with an
//!   explicit `where A: Clone` (the old crate got `Clone` free from
//!   `Coefficient: Clone`). See [`RotationOne`] for why `impl Into<A>` is not
//!   available here.
//!
//! Nothing else moves: the defaults (`reset_x` = `reset` then `h`, `x_error(p)` =
//! `pauli_error([p, 0, 0])`, `rx` = `rotate_1(Pauli::X, …)`, every `*_many` loop)
//! are byte-for-byte the old bodies, so a caller sees the same behaviour and the
//! same override points.

use crate::coefficient::{Angle, Coefficient};
use crate::word::Pauli;

/// The Clifford gate set, applied in the Heisenberg picture.
///
/// `Clifford` is a *derived* behavioral trait: it is **not** implemented by
/// hand on each type but blanket-implemented once over the Pauli algebra
/// primitives ([`crate::pauli::SymplecticColumns`] + [`crate::pauli::PhaseTrack`],
/// see `pauli.rs`), and separately on `Sum` in terms of its key's `Clifford`.
///
/// Design: §"Behavioral traits" and §"Pauli algebra traits". Each generator is
/// an `Sp(2n, 2)` isometry on the symplectic bits — machine-checked per
/// generator in `lean/PPVM/Pauli/Symplectic.lean`
/// (`hAct_isometry`/`sAct_isometry`/`cnotAct_isometry`/`czAct_isometry`) — with
/// the sign action of `lean/PPVM/Pauli/Conjugation.lean` (`conjH_Y`: `HYH = −Y`,
/// etc.).
pub trait Clifford {
    /// Apply Pauli `X` to one qubit.
    fn x(&mut self, qubit: usize);
    /// Apply Pauli `Y` to one qubit.
    fn y(&mut self, qubit: usize);
    /// Apply Pauli `Z` to one qubit.
    fn z(&mut self, qubit: usize);
    /// Apply Hadamard `H` to one qubit.
    fn h(&mut self, qubit: usize);
    /// Apply the phase gate `S` to one qubit.
    fn s(&mut self, qubit: usize);
    /// Apply `CNOT` to one `(control, target)` pair.
    fn cnot(&mut self, control: usize, target: usize);
    /// Apply `CZ` to one qubit pair.
    fn cz(&mut self, qubit0: usize, qubit1: usize);

    /// stim alias for [`cnot`](Clifford::cnot).
    fn cx(&mut self, control: usize, target: usize) {
        self.cnot(control, target)
    }
    /// stim alias for [`cnot`](Clifford::cnot).
    fn zcx(&mut self, control: usize, target: usize) {
        self.cnot(control, target)
    }
    /// stim alias for [`cz`](Clifford::cz).
    fn zcz(&mut self, qubit0: usize, qubit1: usize) {
        self.cz(qubit0, qubit1)
    }
}

/// Additional Clifford gates beyond the minimal set: `S†`, `√X`, `√X†`, `√Y`,
/// `√Y†`, and `CY`.
///
/// Same shape as the old `ppvm_traits::traits::CliffordExtensions`: the six
/// generators are required methods (so a tableau can fuse them into one pass over
/// its bit planes) and `zcy` is the defaulted stim alias.
///
/// Like [`Clifford`], this is a *derived* trait — it is blanket-implemented once
/// over the Pauli algebra primitives for the types that opt in via
/// [`BlanketClifford`](crate::pauli::BlanketClifford) (see `pauli.rs`), which is
/// where the old crate's `impl<T: PauliWordTrait> CliffordExtensions for T`
/// blanket lands.
///
/// # Conjugation table (backward Heisenberg convention)
///
/// Each gate conjugates as `P ↦ U†PU` — the convention the phased word and the
/// Lean `conjSdag` pin (`lean/PPVM/Pauli/Conjugation.lean`, the `S`/`S†` note).
/// It is the table of the old crate's blanket impl, unchanged:
///
/// | Gate | `X` | `Y` | `Z` |
/// |:---:|:---:|:---:|:---:|
/// | `s` | `-Y` | `X` | `Z` |
/// | `s_dag` | `Y` | `-X` | `Z` |
/// | `sqrt_x` | `X` | `-Z` | `Y` |
/// | `sqrt_x_dag` | `X` | `Z` | `-Y` |
/// | `sqrt_y` | `Z` | `Y` | `-X` |
/// | `sqrt_y_dag` | `-Z` | `Y` | `X` |
pub trait CliffordExtensions: Clifford {
    /// Apply `S†` to one qubit.
    fn s_dag(&mut self, qubit: usize);
    /// Apply `√X` to one qubit.
    fn sqrt_x(&mut self, qubit: usize);
    /// Apply `(√X)†` to one qubit.
    fn sqrt_x_dag(&mut self, qubit: usize);
    /// Apply `√Y` to one qubit.
    fn sqrt_y(&mut self, qubit: usize);
    /// Apply `(√Y)†` to one qubit.
    fn sqrt_y_dag(&mut self, qubit: usize);
    /// Apply `CY` to one `(control, target)` pair.
    fn cy(&mut self, control: usize, target: usize);
    /// stim alias for [`cy`](CliffordExtensions::cy).
    fn zcy(&mut self, control: usize, target: usize) {
        self.cy(control, target)
    }
}

/// Batched Clifford gates: apply the same gate to many qubits in one call.
///
/// The default implementations loop over the corresponding single-qubit (or
/// single-pair) [`Clifford`] method — identical to the old crate, including the
/// convention that a type with no specialization opts in with an empty `impl`.
/// Types like the stabilizer `Tableau` override with a fused inner loop or a
/// bitmask sweep.
pub trait CliffordBatch: Clifford {
    /// Apply Pauli `X` to every qubit in `indices`.
    fn x_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.x(q);
        }
    }
    /// Apply Pauli `Y` to every qubit in `indices`.
    fn y_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.y(q);
        }
    }
    /// Apply Pauli `Z` to every qubit in `indices`.
    fn z_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.z(q);
        }
    }
    /// Apply Hadamard `H` to every qubit in `indices`.
    fn h_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.h(q);
        }
    }
    /// Apply the phase gate `S` to every qubit in `indices`.
    fn s_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.s(q);
        }
    }
    /// Apply `CNOT` to every `(control, target)` pair.
    fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cnot(c, t);
        }
    }
    /// Apply `CZ` to every `(control, target)` pair.
    fn cz_many(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cz(c, t);
        }
    }
}

/// Batched form of [`CliffordExtensions`], with the same loop defaults as the old
/// crate.
pub trait CliffordExtensionsBatch: CliffordExtensions + CliffordBatch {
    /// Apply `S†` to every qubit in `indices`.
    fn s_dag_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.s_dag(q);
        }
    }
    /// Apply `√X` to every qubit in `indices`.
    fn sqrt_x_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_x(q);
        }
    }
    /// Apply `(√X)†` to every qubit in `indices`.
    fn sqrt_x_dag_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_x_dag(q);
        }
    }
    /// Apply `√Y` to every qubit in `indices`.
    fn sqrt_y_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_y(q);
        }
    }
    /// Apply `(√Y)†` to every qubit in `indices`.
    fn sqrt_y_dag_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_y_dag(q);
        }
    }
    /// Apply `CY` to every `(control, target)` pair.
    fn cy_many(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cy(c, t);
        }
    }
}

/// Single-qubit rotations parameterized by an angle domain `A` that yields
/// amplitudes in coefficient domain `C`.
///
/// The angle defaults to the coefficient (`A = C`), recovering today's
/// `rx(theta: C)` while permitting a symbolic/parametric angle over an
/// `f64`-coefficient sum.
///
/// Design: §"Behavioral traits" (`RotationOne`). The branch each rotation stages
/// — `c·P → cos·c·P + sin·c·(iGP)` — produces exactly one genuinely-new term and
/// is a norm-preserving, angle-additive 2-D rotation on the coefficient pair,
/// machine-checked in `lean/PPVM/Instantiations/Rotation.lean`
/// (`anticommute_new_key`, `rot_norm_sq`, `rot_rot`).
///
/// # Behaviour parity with `ppvm-traits`
///
/// The required/defaulted split is the old
/// `ppvm_traits::RotationOne`'s (`ppvm-traits/src/traits/branch/rot1.rs`):
/// [`rotate_1`](RotationOne::rotate_1) is the axis-generic entry point that a
/// two-qubit rotation composes with (old `ppvm-pauli-sum/src/sum/rot2.rs`,
/// `ppvm-tableau-sum/src/gates/rot1.rs`), `rx`/`ry`/`rz` default onto it (a
/// backend that has a per-axis fast path overrides them, as `PauliSum` does),
/// and the `*_many` forms are the batch loops the Python bindings call
/// (`ppvm-python-native/src/interface.rs`).
///
/// Two edits are forced by the `<T: Config>` → `<C, A>` shape change:
///
/// * the old `rotate_1` panicked on the `Pauli::L` axis. `L` is not a variant of
///   the redesigned [`Pauli`] (loss is a `LossySite`, not a Pauli letter), so
///   that panic is unrepresentable rather than removed. `Pauli::I` behaves as it
///   always did: it commutes with everything, so the pass is a no-op.
/// * the old methods took `theta: impl Into<T::Coeff>`, which the `Coefficient:
///   From<f64>` bound made useful (`sum.rx(0, 0.1)` on a `Complex<f64>` sum).
///   Here the angle `A` is a *free trait parameter*, so `impl Into<A>` would
///   leave `A` unconstrained at the call site and inference would fail; the
///   parameter therefore stays `A`. The one instantiation callers actually used
///   is preserved by [`Angle<Complex<f64>> for f64`](crate::coefficient), which
///   lets `sum.rx(0, 0.1)` keep working on a complex-coefficient sum.
pub trait RotationOne<C: Coefficient, A: Angle<C> = C> {
    /// Rotate about `axis` (one of `X`, `Y`, `Z`) on `qubit` by `theta`.
    ///
    /// `Pauli::I` commutes with every term, so an `I` axis is a no-op.
    fn rotate_1(&mut self, axis: Pauli, qubit: usize, theta: A);

    /// Rotate about `X` on `qubit` by `theta`.
    fn rx(&mut self, qubit: usize, theta: A) {
        self.rotate_1(Pauli::X, qubit, theta)
    }
    /// Rotate about `Y` on `qubit` by `theta`.
    fn ry(&mut self, qubit: usize, theta: A) {
        self.rotate_1(Pauli::Y, qubit, theta)
    }
    /// Rotate about `Z` on `qubit` by `theta`.
    fn rz(&mut self, qubit: usize, theta: A) {
        self.rotate_1(Pauli::Z, qubit, theta)
    }

    /// Explicit batched `RX(θ)`.
    fn rx_many(&mut self, targets: &[usize], theta: A)
    where
        A: Clone,
    {
        for &q in targets {
            self.rx(q, theta.clone())
        }
    }
    /// Explicit batched `RY(θ)`.
    fn ry_many(&mut self, targets: &[usize], theta: A)
    where
        A: Clone,
    {
        for &q in targets {
            self.ry(q, theta.clone())
        }
    }
    /// Explicit batched `RZ(θ)`.
    fn rz_many(&mut self, targets: &[usize], theta: A)
    where
        A: Clone,
    {
        for &q in targets {
            self.rz(q, theta.clone())
        }
    }
}

/// One named two-qubit rotation plus its batch form, defaulted onto
/// [`RotationTwo::rotate_2`] — the old crate's `def_rotation!`
/// (`ppvm-traits/src/traits/branch/rot2.rs`) with `impl Into<T::Coeff>` replaced
/// by the angle domain `A`.
macro_rules! def_rotation_two {
    ($name:ident, $batch:ident, $x_a:expr, $z_a:expr, $x_b:expr, $z_b:expr, $doc:literal) => {
        #[doc = $doc]
        fn $name(&mut self, a: usize, b: usize, theta: A) {
            self.rotate_2([$x_a, $z_a], [$x_b, $z_b], a, b, theta)
        }

        #[doc = concat!("Explicit batched form of [`RotationTwo::", stringify!($name), "`].")]
        fn $batch(&mut self, pairs: &[(usize, usize)], theta: A)
        where
            A: Clone,
        {
            for &(a, b) in pairs {
                self.$name(a, b, theta.clone())
            }
        }
    };
}

/// Two-qubit Pauli rotations, generated by `P_a ⊗ P_b` for any pair of
/// non-identity Paulis: `exp(−i·θ/2·P_a ⊗ P_b)`. Provides the named convenience
/// methods `rxx`, `rxy`, …, `rzz` on top of the generic
/// [`rotate_2`](RotationTwo::rotate_2).
///
/// Ported from `ppvm_traits::RotationTwo` (`ppvm-traits/src/traits/branch/rot2.rs`)
/// with the crate-wide shape edit only: the `[x, z]` axis encoding, the nine
/// named gates, their argument order, and the `*_many` loops are the old ones.
/// The batch defaults need `A: Clone`, which the old crate got for free from
/// `Coefficient: Clone`.
pub trait RotationTwo<C: Coefficient, A: Angle<C> = C> {
    /// Two-qubit Pauli rotation `exp(−i·θ/2·P_a ⊗ P_b)`.
    ///
    /// Each axis is encoded as `[x, z]` bits:
    /// `[0,0]` = I, `[1,0]` = X, `[0,1]` = Z, `[1,1]` = Y.
    fn rotate_2(&mut self, axis_a: [u8; 2], axis_b: [u8; 2], a: usize, b: usize, theta: A);

    //                 x, z, x, z
    def_rotation_two!(rxx, rxx_many, 1, 0, 1, 0, "`exp(-i θ/2 · X_a X_b)`.");
    def_rotation_two!(rxy, rxy_many, 1, 0, 1, 1, "`exp(-i θ/2 · X_a Y_b)`.");
    def_rotation_two!(rxz, rxz_many, 1, 0, 0, 1, "`exp(-i θ/2 · X_a Z_b)`.");

    def_rotation_two!(ryx, ryx_many, 1, 1, 1, 0, "`exp(-i θ/2 · Y_a X_b)`.");
    def_rotation_two!(ryy, ryy_many, 1, 1, 1, 1, "`exp(-i θ/2 · Y_a Y_b)`.");
    def_rotation_two!(ryz, ryz_many, 1, 1, 0, 1, "`exp(-i θ/2 · Y_a Z_b)`.");

    def_rotation_two!(rzx, rzx_many, 0, 1, 1, 0, "`exp(-i θ/2 · Z_a X_b)`.");
    def_rotation_two!(rzy, rzy_many, 0, 1, 1, 1, "`exp(-i θ/2 · Z_a Y_b)`.");
    def_rotation_two!(rzz, rzz_many, 0, 1, 0, 1, "`exp(-i θ/2 · Z_a Z_b)`.");
}

/// Rotation about an axis in the x/y plane:
/// `R(axis_angle, θ) = exp(−i·θ/2·(cos(axis_angle)·X + sin(axis_angle)·Y))`.
///
/// The in-plane axis is `X` rotated about `Z` by `axis_angle`, so
/// `R(axis_angle, θ) = RZ(axis_angle)·RX(θ)·RZ(−axis_angle)` — the old
/// `ppvm_traits::RotXY`, unchanged apart from the angle domain.
pub trait RotXY<C: Coefficient, A: Angle<C> = C> {
    /// `R(axis_angle, θ)` on `qubit`.
    fn r(&mut self, qubit: usize, axis_angle: A, theta: A);
}

/// Controlled `RX` rotation — the old `ppvm_traits::CRx`.
pub trait CRx<C: Coefficient, A: Angle<C> = C> {
    /// Apply `CRX(θ)` with the given control and target.
    fn crx(&mut self, control: usize, target: usize, theta: A);
}

/// The general single-qubit `U3(θ, φ, λ)` gate — the old
/// `ppvm_traits::U3Gate`.
pub trait U3Gate<C: Coefficient, A: Angle<C> = C> {
    /// Apply `U3(θ, φ, λ)` to `qubit`.
    fn u3(&mut self, qubit: usize, theta: A, phi: A, lambda: A);
}

/// The non-Clifford `T` gate and its adjoint, `T = diag(1, e^{iπ/4})`.
///
/// Ported from `ppvm_traits::TGate` with the same required/defaulted split and
/// the same `*_many` loop bodies. The old trait carried a `<T: Config>`
/// parameter it never used; `T` takes no numeric argument, so — by the rule that
/// already unparameterized `Clifford` and `ResetLossChannel` — it carries none.
pub trait TGate {
    /// Apply `T` (`diag(1, e^{iπ/4})`) to one qubit.
    fn t(&mut self, qubit: usize);
    /// Apply `T†` to one qubit.
    fn t_dag(&mut self, qubit: usize);

    /// Explicit batched `T`.
    fn t_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.t(q);
        }
    }

    /// Explicit batched `T†`.
    fn t_dag_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.t_dag(q);
        }
    }
}

/// Projective Z-basis projectors `|0⟩⟨0|` and `|1⟩⟨1|` — the old
/// `ppvm_traits::Projection`, verbatim (it never carried a `Config`).
pub trait Projection {
    /// Project `qubit` onto `|0⟩`.
    fn p0(&mut self, qubit: usize);
    /// Project `qubit` onto `|1⟩`.
    fn p1(&mut self, qubit: usize);
}

/// A unital single-qubit Pauli error channel `P ↦ λ_P·P`.
///
/// Design: §"Behavioral traits" (`PauliError`). Acting diagonally in the Pauli
/// basis, its transfer eigenvalue collapses (using `Σ_Q p_Q = 1`) to
/// `λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q`, machine-checked in
/// `lean/PPVM/Algebra/Noise.lean` (`pauli_channel_eigenvalue`, and
/// `pauli_channel_eigenvalue_omega` tying anticommutation to
/// `PPVM.Symplectic.omega`).
pub trait PauliError<C: Coefficient> {
    /// Apply a single-qubit Pauli channel with `X`, `Y`, `Z` probabilities.
    fn pauli_error(&mut self, qubit: usize, probabilities: [C; 3]);

    /// stim `X_ERROR(p)` — apply `X` with probability `p` to one qubit.
    fn x_error(&mut self, qubit: usize, p: C) {
        let zero = C::zero();
        self.pauli_error(qubit, [p, zero.clone(), zero])
    }

    /// stim `Y_ERROR(p)` — apply `Y` with probability `p` to one qubit.
    fn y_error(&mut self, qubit: usize, p: C) {
        let zero = C::zero();
        self.pauli_error(qubit, [zero.clone(), p, zero])
    }

    /// stim `Z_ERROR(p)` — apply `Z` with probability `p` to one qubit.
    fn z_error(&mut self, qubit: usize, p: C) {
        let zero = C::zero();
        self.pauli_error(qubit, [zero.clone(), zero, p])
    }

    /// Explicit batched Pauli-error channel.
    fn pauli_error_many(&mut self, targets: &[usize], p: [C; 3]) {
        for &q in targets {
            self.pauli_error(q, p.clone());
        }
    }

    /// Explicit batched `X_ERROR(p)`.
    fn x_error_many(&mut self, targets: &[usize], p: C) {
        for &q in targets {
            self.x_error(q, p.clone());
        }
    }

    /// Explicit batched `Y_ERROR(p)`.
    fn y_error_many(&mut self, targets: &[usize], p: C) {
        for &q in targets {
            self.y_error(q, p.clone());
        }
    }

    /// Explicit batched `Z_ERROR(p)`.
    fn z_error_many(&mut self, targets: &[usize], p: C) {
        for &q in targets {
            self.z_error(q, p.clone());
        }
    }
}

/// Apply the same single-qubit Pauli error channel uniformly to every qubit in
/// the system.
pub trait PauliErrorAll<C: Coefficient> {
    /// Apply the Pauli channel `p = [p_x, p_y, p_z]` to every qubit.
    fn pauli_error_all(&mut self, p: [C; 3]);
}

/// Two-qubit Pauli error channel.
pub trait TwoQubitPauliError<C: Coefficient> {
    /// Apply a two-qubit Pauli-error channel to one pair. Probabilities are given
    /// in the order
    /// `{IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ}`.
    fn two_qubit_pauli_error(&mut self, qubit0: usize, qubit1: usize, p: [C; 15]);

    /// Explicit batched two-qubit Pauli-error channel.
    fn two_qubit_pauli_error_many(&mut self, pairs: &[(usize, usize)], p: [C; 15]) {
        for &(a, b) in pairs {
            self.two_qubit_pauli_error(a, b, p.clone());
        }
    }
}

/// Single-qubit depolarizing channel.
pub trait Depolarizing<C: Coefficient> {
    /// Depolarize one qubit with probability `p`.
    fn depolarize1(&mut self, qubit: usize, p: C);

    /// Explicit batched single-qubit depolarizing channel.
    fn depolarize1_many(&mut self, targets: &[usize], p: C) {
        for &q in targets {
            self.depolarize1(q, p.clone());
        }
    }
}

/// Two-qubit depolarizing channel.
pub trait Depolarizing2<C: Coefficient> {
    /// Depolarize one qubit pair with probability `p`.
    fn depolarize2(&mut self, qubit0: usize, qubit1: usize, p: C);

    /// Explicit batched two-qubit depolarizing channel.
    fn depolarize2_many(&mut self, pairs: &[(usize, usize)], p: C) {
        for &(a, b) in pairs {
            self.depolarize2(a, b, p.clone());
        }
    }
}

/// Amplitude-damping channel (single qubit).
pub trait AmplitudeDamping<C: Coefficient> {
    /// Apply amplitude damping with damping parameter `gamma`.
    fn amplitude_damping(&mut self, qubit: usize, gamma: C);
}

/// Single-qubit loss channel — with probability `p`, mark the qubit as lost
/// ([`LossySite::Lost`](crate::word::LossySite::Lost)).
pub trait LossChannel<C: Coefficient> {
    /// Apply a loss channel to `qubit` with loss probability `p`.
    fn loss_channel(&mut self, qubit: usize, p: C);
}

/// Correlated two-qubit loss channel.
pub trait CorrelatedLossChannel<C: Coefficient> {
    /// Apply a correlated loss channel to `qubit0` and `qubit1`.
    ///
    /// The three probabilities are:
    /// * `p[0]`: losing both qubits simultaneously when both are in the qubit
    ///   subspace.
    /// * `p[1]`: losing either one qubit when both are in the qubit subspace.
    /// * `p[2]`: losing one qubit when the other has already been lost prior to
    ///   the channel.
    fn correlated_loss_channel(&mut self, qubit0: usize, qubit1: usize, p: [C; 3]);
}

/// Reset the loss bit on a qubit — models a re-cooling / re-loading event that
/// brings a previously-lost atom back.
///
/// Unlike the old `ResetLossChannel<T: Config>`, this takes no coefficient
/// parameter: clearing a loss bit consumes no probability, and the design's rule
/// is that an operation with no numeric parameter carries none (as `Clifford`
/// does).
pub trait ResetLossChannel {
    /// Clear the loss bit at `qubit`.
    fn reset_loss_channel(&mut self, qubit: usize);
}

/// State-dependent ("asymmetric") single-qubit loss channel: a qubit is lost from
/// `|0⟩` with probability `p0` and from `|1⟩` with probability `p1`. Unlike
/// [`LossChannel`], the total loss probability depends on the qubit's
/// populations, so the channel reads the current `⟨Z⟩`.
pub trait AsymmetricLossChannel<C: Coefficient> {
    /// Apply asymmetric loss to `qubit`, with `p0` / `p1` the loss probabilities
    /// from `|0⟩` / `|1⟩`. See the backend impl for the trajectory approximation
    /// used (the survival back-action is omitted).
    fn asymmetric_loss_channel(&mut self, qubit: usize, p0: C, p1: C);
}

/// Loss-aware projective computational-basis measurement.
///
/// `Some(false)` and `Some(true)` denote the `|0⟩`/`|1⟩` outcomes and `None`
/// denotes a lost qubit — the former `Measure -> bool` / `LossyMeasure ->
/// Option<bool>` split is removed. Sharing the result type does not share the
/// algorithm: `Tableau` uses the pure Clifford procedure, `GeneralizedTableau`
/// the coefficient-aware `O(n²)` decomposition.
///
/// Design: §"Behavioral traits" (`Measure`). The deterministic-vs-random
/// dichotomy the pivot search rests on is machine-checked in
/// `lean/PPVM/Tableau/Frame.lean` (`measurement_dichotomy`,
/// `measure_deterministic_iff_xfree`).
pub trait Measure {
    /// Measure `qubit`; `None` if the qubit has been lost.
    fn measure(&mut self, qubit: usize) -> Option<bool>;

    /// Measure each target in order, one result per target.
    fn measure_many(&mut self, targets: &[usize]) -> Vec<Option<bool>> {
        targets.iter().map(|&q| self.measure(q)).collect()
    }
}

/// Reset one qubit to a computational/Pauli basis state.
///
/// The supertrait bound and every default body are the old crate's
/// (`ppvm-traits/src/traits/reset.rs`): only `reset` (stim `R`/`RZ`) is required,
/// and the basis variants are *defined* as `reset` followed by the basis-change
/// Cliffords — `reset_x` = `reset` then `h`, `reset_y` = `reset` then `h` then
/// `s`. Those compositions are behaviour, not an implementation detail, so they
/// are reproduced call-for-call rather than re-derived.
pub trait Reset: Clifford + CliffordExtensions {
    /// Reset one qubit to `|0⟩` (stim `R`/`RZ`).
    fn reset(&mut self, qubit: usize);

    /// stim `RZ` alias — reset to `|0⟩`.
    fn reset_z(&mut self, qubit: usize) {
        self.reset(qubit)
    }

    /// stim `RX` — reset to `|+⟩`.
    fn reset_x(&mut self, qubit: usize) {
        self.reset(qubit);
        self.h(qubit);
    }

    /// stim `RY` — reset to `|i⟩`.
    fn reset_y(&mut self, qubit: usize) {
        self.reset(qubit);
        self.h(qubit);
        self.s(qubit);
    }

    /// Explicit batched reset to `|0⟩`.
    fn reset_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset(q);
        }
    }

    /// Explicit batched `RZ` alias.
    fn reset_z_many(&mut self, targets: &[usize]) {
        self.reset_many(targets)
    }

    /// Explicit batched `RX`.
    fn reset_x_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset_x(q);
        }
    }

    /// Explicit batched `RY`.
    fn reset_y_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset_y(q);
        }
    }
}
