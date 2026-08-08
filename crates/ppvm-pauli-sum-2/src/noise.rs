// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`PauliError`] propagation for a Pauli-keyed [`Sum`] (the `PauliSum` alias):
//! the **diagonal**, in-place unital Pauli channel.
//!
//! `pauli_error(qubit, [pX, pY, pZ])` is the unital single-qubit Pauli channel
//! `ρ ↦ Σ_Q p_Q · Q ρ Q` (with `p_I = 1 − pX − pY − pZ` implicit). In the
//! Heisenberg picture it acts **diagonally** in the Pauli basis: each term `P` is
//! scaled by its real transfer eigenvalue
//!
//! ```text
//! λ_P = Σ_Q p_Q (−1)^{ω(P, Q)} = 1 − 2·Σ_{Q anticommutes with P at the qubit} p_Q,
//! ```
//!
//! a factor depending only on `P`'s Pauli at `qubit` (the collapse uses
//! `Σ_Q p_Q = 1`). So — no branching, no rebuild — this is a pure per-key
//! coefficient scale through the in-place [`Sum::scale_by_key`] fast path,
//! restoring the old crate's in-place `scale` (`ppvm-pauli-sum::sum::noise`).
//! Writing `X`/`Y`/`Z` for `P`'s Pauli at the qubit:
//!
//! ```text
//! λ_X = 1 − 2(pY + pZ),   λ_Y = 1 − 2(pX + pZ),   λ_Z = 1 − 2(pX + pY),   λ_I = 1.
//! ```
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Behavioral traits"
//! (`PauliError`). The eigenvalue — and its tie to anticommutation via
//! `PPVM.Symplectic.omega` — is machine-checked in `lean/PPVM/Algebra/Noise.lean`
//! (`pauli_channel_eigenvalue`, `pauli_channel_eigenvalue_omega`).
//!
//! # Friction: the eigenvalue needs a ring `1`, which `Coefficient` omits
//!
//! `λ_P = 1 − 2·(…)` needs the coefficient ring's multiplicative identity to
//! build, but [`Coefficient`](ppvm_traits_2::Coefficient) bounds only `num::Zero`
//! (the additive identity), having dropped every `1`-flavored bound — `Mul<f64>`,
//! `From<f64>` — that excluded exact rings. The minimal fix is a `C: num::One`
//! bound on **this impl** (not the trait), which every intended noise coefficient
//! (`f64`, `Complex<f64>`) satisfies; the doubling `2·x` is formed as `x + x`, so
//! no bare-`f64` scaling sneaks back in.

use ppvm_traits_2::{
    Accumulate, AmplitudeDamping, Coefficient, Depolarizing, Depolarizing2, Indexable, PauliBits,
    PauliError, Retain, TwoQubitPauliError, Word,
};

use crate::store::{RotateInPlace, ScaleByKey};
use crate::sum::Sum;

/// Diagonal unital Pauli-channel propagation on a Pauli-keyed `Sum`: an in-place
/// per-term scale by the channel's real transfer eigenvalue.
impl<S, P, W, C> PauliError<C> for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + ScaleByKey<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient + num::One,
    P: crate::policy::Policy<W, C>,
{
    /// Scale each term in place by `λ_P` for its Pauli at `qubit`. The three
    /// non-trivial eigenvalues are computed once; the identity term (and any lost
    /// qubit) is left exactly untouched.
    ///
    /// Mirrors the old `ppvm-pauli-sum::sum::noise`'s `self.scale(|k, v| ...)`
    /// verbatim in shape as well as arithmetic: an in-place mutation that cannot
    /// remove, so a zero eigenvalue leaves a zero-coefficient term in the support
    /// exactly as old does.
    #[inline(always)]
    fn pauli_error(&mut self, qubit: usize, probabilities: [C; 3]) {
        let [x_factor, z_factor, y_factor] = C::pauli_error_factors(probabilities);
        self.scale_pauli_error(qubit, x_factor, z_factor, y_factor);
    }

    #[inline(always)]
    fn pauli_error_many(&mut self, targets: &[usize], probabilities: [C; 3]) {
        let [x_factor, z_factor, y_factor] = C::pauli_error_factors(probabilities);
        self.scale_by_key(move |key, coeff| {
            for &qubit in targets {
                if key.is_lost(qubit) {
                    continue;
                }
                match key.pauli_code(qubit) {
                    0 => {}
                    1 => *coeff *= x_factor.clone(),
                    2 => *coeff *= z_factor.clone(),
                    3 => *coeff *= y_factor.clone(),
                    _ => unreachable!(),
                }
            }
        });
    }

    #[inline(always)]
    fn x_error_many(&mut self, targets: &[usize], p: C) {
        let zero = C::zero();
        self.pauli_error_many(targets, [p, zero.clone(), zero]);
    }

    #[inline(always)]
    fn y_error_many(&mut self, targets: &[usize], p: C) {
        let zero = C::zero();
        self.pauli_error_many(targets, [zero.clone(), p, zero]);
    }

    #[inline(always)]
    fn z_error_many(&mut self, targets: &[usize], p: C) {
        let zero = C::zero();
        self.pauli_error_many(targets, [zero.clone(), zero, p]);
    }
}

/// `1 − 2·(p[i₀] + p[i₁] + …)` over a hand-written index list — old's
/// `one_minus_two_sum` (`ppvm-pauli-sum/src/sum/noise.rs`), with the doubling
/// formed as the ring addition `x + x` instead of a bare-`f64` `* 2.0` (see the
/// module's `num::One` friction note). The accumulation order is old's, so the
/// float rounding is bit-identical.
#[inline]
fn one_minus_two_sum<C: Coefficient + num::One, const N: usize>(
    p: &[C; 15],
    indices: [usize; N],
) -> C {
    let mut acc = C::one();
    for i in indices {
        acc = acc - (p[i].clone() + p[i].clone());
    }
    acc
}

/// The two-qubit Pauli channel: a **diagonal** in-place scale by the pair's
/// transfer eigenvalue `λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q`.
///
/// The 16 index lists are ported **verbatim** from
/// `ppvm-pauli-sum/src/sum/noise.rs`; each names the 8 two-qubit Paulis that
/// anticommute with the observed pair, in the probability order
/// `{IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ}`. They were
/// written by hand in old with no derivation in the source, and only one-hot
/// probability vectors are covered by old's tests, so a transposed index would be
/// invisible on a mixed vector — the enumeration is flagged for the Lean oracle
/// (suspected old bug 6) rather than re-derived here, since re-deriving it *is* a
/// behaviour change if old is right.
///
/// Rides [`Sum::scale_by_key`] (behavioural contract 4): a pure in-place walk that
/// can neither insert nor remove, whatever the factors — including a factor of
/// zero. A lost site on either qubit silently applies **no** noise (old's `_ =>`
/// arm: "if just one atom is lost, then there is no well-defined noise channel on
/// the other atom").
impl<S, P, W, C> TwoQubitPauliError<C> for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + ScaleByKey<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient + num::One,
    P: crate::policy::Policy<W, C>,
{
    fn two_qubit_pauli_error(&mut self, qubit0: usize, qubit1: usize, p: [C; 15]) {
        self.scale_by_key(move |k: &W, c: &mut C| {
            if k.is_lost(qubit0) || k.is_lost(qubit1) {
                return;
            }
            // 2-bit Pauli code (x, z): 0 I, 1 X, 2 Z, 3 Y.
            let a = (k.x_bit(qubit0) as u8) | ((k.z_bit(qubit0) as u8) << 1);
            let b = (k.x_bit(qubit1) as u8) | ((k.z_bit(qubit1) as u8) << 1);
            match (a, b) {
                (0, 0) => {}
                (0, 1) => *c *= one_minus_two_sum(&p, [1, 10, 13, 14, 2, 5, 6, 9]),
                (0, 3) => *c *= one_minus_two_sum(&p, [0, 10, 12, 14, 2, 4, 6, 8]),
                (0, 2) => *c *= one_minus_two_sum(&p, [0, 1, 12, 13, 4, 5, 8, 9]),
                (1, 0) => *c *= one_minus_two_sum(&p, [10, 11, 12, 13, 14, 7, 8, 9]),
                (1, 1) => *c *= one_minus_two_sum(&p, [1, 11, 12, 2, 5, 6, 7, 8]),
                (1, 3) => *c *= one_minus_two_sum(&p, [0, 11, 13, 2, 4, 6, 7, 9]),
                (1, 2) => *c *= one_minus_two_sum(&p, [0, 1, 10, 11, 14, 4, 5, 7]),
                (3, 0) => *c *= one_minus_two_sum(&p, [11, 12, 13, 14, 3, 4, 5, 6]),
                (3, 1) => *c *= one_minus_two_sum(&p, [1, 10, 11, 12, 2, 3, 4, 9]),
                (3, 3) => *c *= one_minus_two_sum(&p, [0, 10, 11, 13, 2, 3, 5, 8]),
                (3, 2) => *c *= one_minus_two_sum(&p, [0, 1, 11, 14, 3, 6, 8, 9]),
                (2, 0) => *c *= one_minus_two_sum(&p, [10, 3, 4, 5, 6, 7, 8, 9]),
                (2, 1) => *c *= one_minus_two_sum(&p, [1, 13, 14, 2, 3, 4, 7, 8]),
                (2, 3) => *c *= one_minus_two_sum(&p, [0, 12, 14, 2, 3, 5, 7, 9]),
                (2, 2) => *c *= one_minus_two_sum(&p, [0, 1, 10, 12, 13, 3, 6, 7]),
                _ => unreachable!("2-bit Pauli code is 0..=3"),
            }
        });
    }
}

/// Single-qubit depolarizing: scale every **non-identity** Pauli at the qubit by
/// `1 − 4p/3`, in place ([`Sum::scale_by_key`], contract 4).
///
/// # Friction: `4p/3` needs a rational scale the `Coefficient` ring does not carry
///
/// Old computes `T::Coeff::from(1.0) - p * (4.0 / 3.0)`, leaning on the old
/// `Coefficient: From<f64> + Mul<f64>` bounds. `-2`'s [`Coefficient`] drops both —
/// deliberately, so exact rings qualify — and `4/3` is not expressible from
/// `One`/`Zero`/ring ops at all. The bound is therefore re-added **on this impl
/// only**, as `Mul<f64, Output = C>`: `f64` and `Complex<f64>` satisfy it, an
/// exact ring rightly does not (a depolarizing channel at `4p/3` leaves `ℤ[i]`
/// regardless).
impl<S, P, W, C> Depolarizing<C> for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + ScaleByKey<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient + num::One + std::ops::Mul<f64, Output = C>,
    P: crate::policy::Policy<W, C>,
{
    #[inline(always)]
    fn depolarize1(&mut self, qubit: usize, p: C) {
        let factor = C::one() - p * (4.0 / 3.0);
        self.scale_by_key(move |k: &W, c: &mut C| {
            if !k.is_lost(qubit) && (k.x_bit(qubit) || k.z_bit(qubit)) {
                *c *= factor.clone();
            }
        });
    }

    #[inline(always)]
    fn depolarize1_many(&mut self, targets: &[usize], p: C) {
        let factor = C::one() - p * (4.0 / 3.0);
        self.scale_by_key(move |k: &W, c: &mut C| {
            for &qubit in targets {
                if !k.is_lost(qubit) && (k.x_bit(qubit) || k.z_bit(qubit)) {
                    *c *= factor.clone();
                }
            }
        });
    }
}

/// Two-qubit depolarizing: scale by `1 − 16p/15` every term that is non-identity
/// on **either** qubit, in place. See [`Depolarizing`] for the `Mul<f64>` bound.
impl<S, P, W, C> Depolarizing2<C> for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + ScaleByKey<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient + num::One + std::ops::Mul<f64, Output = C>,
    P: crate::policy::Policy<W, C>,
{
    #[inline(always)]
    fn depolarize2(&mut self, qubit0: usize, qubit1: usize, p: C) {
        let factor = C::one() - p * (16.0 / 15.0);
        self.scale_by_key(move |k: &W, c: &mut C| {
            if k.is_lost(qubit0) || k.is_lost(qubit1) {
                return;
            }
            let a = k.x_bit(qubit0) || k.z_bit(qubit0);
            let b = k.x_bit(qubit1) || k.z_bit(qubit1);
            if a || b {
                *c *= factor.clone();
            }
        });
    }

    #[inline(always)]
    fn depolarize2_many(&mut self, pairs: &[(usize, usize)], p: C) {
        let factor = C::one() - p * (16.0 / 15.0);
        self.scale_by_key(move |k: &W, c: &mut C| {
            for &(qubit0, qubit1) in pairs {
                if k.is_lost(qubit0) || k.is_lost(qubit1) {
                    continue;
                }
                let a = k.x_bit(qubit0) || k.z_bit(qubit0);
                let b = k.x_bit(qubit1) || k.z_bit(qubit1);
                if a || b {
                    *c *= factor.clone();
                }
            }
        });
    }
}

/// Amplitude damping — the one **branching** channel on the non-lossy word:
/// `X`/`Y ↦ √(1−γ)`, `I ↦ I`, and `Z ↦ (1−γ)·Z + γ·I`.
///
/// The `Z → I` branch rides the fused [`Sum::rotate_in_place`] merge, whose
/// [`AddTerm`](crate::AddTerm) is old's accumulating `add_assign` — so on a sum
/// that already carries the `I` term the branch **accumulates** onto it
/// (`I_new = I_old + γ·Z_old`) and never overwrites (behavioural contract 3(c)).
///
/// Old bounds this impl on `T::Coeff: Float` (it needs `sqrt`), i.e. it exists for
/// `f64` and not for `Complex<f64>`; the same bound is kept, so the port neither
/// widens nor narrows the surface.
impl<S, P, W, C> AmplitudeDamping<C> for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + RotateInPlace<W, C> + Retain<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient + num::Float,
    P: crate::policy::Policy<W, C>,
{
    fn amplitude_damping(&mut self, qubit: usize, gamma: C) {
        let survive = (C::one() - gamma).sqrt();
        let keep = C::one() - gamma;
        self.rotate_in_place(move |k: &W, c: &mut C| {
            if k.is_lost(qubit) {
                return None;
            }
            match (k.x_bit(qubit), k.z_bit(qubit)) {
                // I — untouched.
                (false, false) => None,
                // X / Y — a pure in-place damping, no branch.
                (true, _) => {
                    *c *= survive;
                    None
                }
                // Z — branches γ·I off, keeping (1−γ)·Z.
                (false, true) => {
                    let branch = *c * gamma;
                    *c *= keep;
                    Some((k.toggled_bits(qubit, false, true), branch))
                }
            }
        });
    }
}
