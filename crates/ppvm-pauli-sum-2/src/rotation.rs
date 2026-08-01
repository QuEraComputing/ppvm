// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`RotationOne`] propagation for a Pauli-keyed [`Sum`] (the `PauliSum` alias):
//! the non-Clifford **branch**.
//!
//! A single-qubit rotation `exp(−i·θ/2·G)` (with axis Pauli `G` = `X`/`Y`/`Z` for
//! `rx`/`ry`/`rz`) conjugates each stored term `(P, c)` to
//! `c·cosθ·P + c·sinθ·(iGP)` when `G` anticommutes with `P` at the qubit, and
//! leaves it unchanged when they commute. This is a genuine **fan-out** (1 or 2
//! output terms per input), but the fan-out is *lopsided*: the `cos·P` diagonal
//! keeps its key, so only the `sin·(iGP)` branch is a genuinely-new term. Rather
//! than pay [`Sum::apply`]'s batch round-trip — which clones every key, resets the
//! map, and re-hashes the *whole* `2N` fan-out including the untouched diagonals —
//! this drives the fused in-place [`Sum::rotate_in_place`] fast path: it scales
//! each diagonal coefficient by `cosθ` **where it sits** (cached hash intact, no
//! bucket move) and hashes/merges only the ≤`N` branch terms, aggregating any
//! collision (a branch's `iGP` may land on another term) and dropping exact
//! cancellations before the policy truncates. This mirrors the old crate's
//! single-pass `map_insert` (`ppvm-pauli-sum::sum::rot1`), the pure-sign
//! [`Sum::flip_sign_by_key`], and the diagonal [`Sum::scale_by_key`] paths.
//!
//! The `iGP` branch key is a **real** Pauli — the single-qubit anticommuting
//! product `GP = ±iP'` carries one factor of `i`, which the leading `i` of `iGP`
//! cancels — so the branch coefficient stays real and its `±1` sign is drained
//! through [`Coefficient::mul_sign`] exactly as the old crate's
//! `sin.mul_sign(eps)` (`ppvm-pauli-sum::sum::rot1`).
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Behavioral traits"
//! (`RotationOne`) and §"Every gate is a producer feeding `accumulate`". The
//! branch is a norm-preserving, angle-additive 2-D rotation on the coefficient
//! pair whose new key is genuinely new, machine-checked in
//! `lean/PPVM/Instantiations/Rotation.lean` (`anticommute_new_key`, `rot_norm_sq`,
//! `rot_rot`).

use std::hash::BuildHasher;

use ppvm_pauli_word_2::{HashFinalize, PauliStorage, PauliWord};
use ppvm_traits_2::{Accumulate, Angle, Coefficient, PauliBits, Retain, RotationOne};

use crate::store::RotateInPlace;
use crate::sum::Sum;

/// Single-qubit rotation propagation on a Pauli-keyed `Sum`. Each rotation drives
/// the fused in-place [`Sum::rotate_in_place`] fast path: it scales every
/// diagonal coefficient by `cosθ` where it sits and merges only the anticommuting
/// `iGP` branch terms, then the policy truncates.
///
/// The per-axis commute test, flipped bits, and `±1` sign `ε` are ported
/// bit-for-bit from `ppvm-pauli-sum::sum::rot1` (`rx`/`ry`/`rz`). The `ε` column
/// is derived from the Pauli phase of `iGP` in
/// `lean/PPVM/Instantiations/Rotation.lean`
/// (`rx_eps_from_product`/`ry_eps_from_product`/`rz_eps_from_product`, real by
/// `branchExp_isRealPhase`).
impl<S, P, A, H, C, Ang> RotationOne<C, Ang> for Sum<S, P>
where
    S: Accumulate<Key = PauliWord<A, H>, Coeff = C>
        + RotateInPlace<PauliWord<A, H>, C>
        + Retain<PauliWord<A, H>, C>,
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
    C: Coefficient,
    Ang: Angle<C>,
    P: crate::policy::Policy<PauliWord<A, H>, C>,
{
    /// Rotate about `X` on `qubit` by `theta`. A term commutes iff its `z` bit at
    /// `qubit` is clear (`I`/`X`); an anticommuting `Z`/`Y` (flip `x`) branches to
    /// `Y`/`Z` with `ε = −1` if `x` else `+1`.
    #[inline]
    fn rx(&mut self, qubit: usize, theta: Ang) {
        let (sin, cos) = theta.sin_cos();
        self.rotate_in_place(move |k: &PauliWord<A, H>, c: &mut C| {
            // A lost qubit has no rotation action (matches `get_lbit` early-out).
            if k.is_lost(qubit) || !k.z_bit(qubit) {
                return None;
            }
            let x = k.x_bit(qubit);
            let branch = c.clone() * sin.mul_sign(if x { -1 } else { 1 });
            *c *= cos.clone();
            let mut new_key = k.clone();
            new_key.set_x_bit(qubit, !x);
            Some((new_key, branch))
        });
    }

    /// Rotate about `Y` on `qubit` by `theta`. A term commutes iff its `x` and `z`
    /// bits at `qubit` agree (`I`/`Y`); an anticommuting `X`/`Z` (flip `x` and `z`)
    /// branches to `Z`/`X` with `ε = −1` if `z` else `+1`.
    #[inline]
    fn ry(&mut self, qubit: usize, theta: Ang) {
        let (sin, cos) = theta.sin_cos();
        self.rotate_in_place(move |k: &PauliWord<A, H>, c: &mut C| {
            if k.is_lost(qubit) {
                return None;
            }
            let x = k.x_bit(qubit);
            let z = k.z_bit(qubit);
            if x == z {
                return None;
            }
            let branch = c.clone() * sin.mul_sign(if z { -1 } else { 1 });
            *c *= cos.clone();
            let mut new_key = k.clone();
            new_key.set_x_bit(qubit, !x);
            new_key.set_z_bit(qubit, !z);
            Some((new_key, branch))
        });
    }

    /// Rotate about `Z` on `qubit` by `theta`. A term commutes iff its `x` bit at
    /// `qubit` is clear (`I`/`Z`); an anticommuting `X`/`Y` (flip `z`) branches to
    /// `Y`/`X` with `ε = +1` if `z` else `−1`.
    #[inline]
    fn rz(&mut self, qubit: usize, theta: Ang) {
        let (sin, cos) = theta.sin_cos();
        self.rotate_in_place(move |k: &PauliWord<A, H>, c: &mut C| {
            if k.is_lost(qubit) || !k.x_bit(qubit) {
                return None;
            }
            let z = k.z_bit(qubit);
            let branch = c.clone() * sin.mul_sign(if z { 1 } else { -1 });
            *c *= cos.clone();
            let mut new_key = k.clone();
            new_key.set_z_bit(qubit, !z);
            Some((new_key, branch))
        });
    }
}
