// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`Projection`] — the computational-basis projectors `p0`/`p1` on a Pauli-keyed
//! [`Sum`], ported from `ppvm-pauli-sum/src/sum/proj.rs`.
//!
//! In the Heisenberg picture `|0⟩⟨0|` acts on a Pauli term as
//! `P ↦ ½(P + Z_q P)`, which is a branch exactly when `P` is `I` or `Z` at the
//! qubit (`X`/`Y` anticommute with `Z_q` and the two halves cancel — old returns
//! `None`, leaving the coefficient untouched). `p1` is the same branch with the
//! branch term negated. Both ride the fused [`Sum::rotate_in_place`](crate::Sum)
//! path, so — as everywhere in this crate — nothing truncates, nothing reduces,
//! and a branch colliding with an existing key **accumulates** (old's
//! `add_assign`).
//!
//! # Suspected old bug: the halving is quadratic in the coefficient
//!
//! Old's kernel reads
//!
//! ```text
//! let half = v.half();          // half == v/2, a *value*, not the constant ½
//! match k.get(pos) {
//!     Pauli::I => { *v *= half; let nk = k.set_new(pos, Pauli::Z); Some((nk, v.clone())) }
//!     ...
//! ```
//!
//! so it computes `c ↦ c²/2` on both the survivor and the branch instead of
//! `c ↦ c/2`. The map is therefore **not linear** (and not idempotent): it agrees
//! with the projector only where `c == 1`, which is why old's own usage — GHZ-style
//! unit-coefficient stabilizer sums — never exposed it, and no old test or
//! benchmark covers `p0`/`p1` at all.
//!
//! Per the prime directive a behaviour change is allowed only when the Lean
//! oracle adjudicates old wrong. **The oracle has ruled against old.**
//! `lean/PPVM/Instantiations/Projector.lean` proves the intended map is linear
//! (`projLin_add`, `projLin_smul`) and idempotent (`projLin_idem`), while old's
//! `c ↦ c²/2` is neither (`oldStep_not_additive`, `oldStep_not_homogeneous`,
//! `oldProj_not_idem`) and coincides with the correct map exactly on `c ∈ {0, 1}`
//! (`oldStep_eq_half_iff`) — i.e. only on the unit-coefficient stabilizer sums
//! that were old's sole usage. The Lean-correct value is `c/2`; the fix is to
//! build the ring's `½` (`C::one().half()`) **once outside the walk** and multiply
//! by *that*.
//!
//! # Second, independent divergence: the `X`/`Y` arm
//!
//! The same Lean file grounds the projector in genuine `ℤ[i]` matrices with
//! `2Π = I + Z`: `twoProj_conj_I`/`twoProj_conj_Z` agree with the intended map on
//! the `I`/`Z` block, but `twoProj_conj_X`/`twoProj_conj_Y` give
//! `Π X Π = Π Y Π = 0` — the projector *annihilates* `X` and `Y`, where old's
//! `_ => None` leaves them untouched. `projLin_p0_add_p1`
//! exhibits the observable consequence: `p0 + p1` is the identity on `I`/`Z` but
//! **doubles** `X`/`Y`, whereas completeness `Π₀ + Π₁ = 1` forces the dephasing
//! channel. The implementation therefore zeros `X`/`Y` in place. It deliberately
//! keeps the zero-coefficient key: explicit [`Sum::reduce`](crate::Sum::reduce)
//! remains the only structural zero-removal operation.

use ppvm_traits_2::{
    Accumulate, Coefficient, Halvable, Indexable, PauliBits, Projection, Retain, Word,
};

use crate::store::RotateInPlace;
use crate::sum::Sum;

impl<S, P, W, C> Projection for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + RotateInPlace<W, C> + Retain<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient + Halvable + num::One,
    P: crate::policy::Policy<W, C>,
{
    /// Project `qubit` onto `|0⟩`: `I ↦ ½(I + Z)`, `Z ↦ ½(Z + I)`, and
    /// `X, Y ↦ 0`.
    fn p0(&mut self, qubit: usize) {
        let half = C::one().half();
        self.rotate_in_place(move |k: &W, c: &mut C| {
            if k.is_lost(qubit) {
                return None;
            }
            if k.x_bit(qubit) {
                *c *= C::zero();
                return None;
            }
            // I ↦ branch on Z, Z ↦ branch on I: either way the z bit toggles.
            *c *= half.clone();
            Some((k.toggled_bits(qubit, false, true), c.clone()))
        });
    }

    /// Project `qubit` onto `|1⟩`: as [`p0`](Projection::p0) with the branch term
    /// negated.
    fn p1(&mut self, qubit: usize) {
        let half = C::one().half();
        self.rotate_in_place(move |k: &W, c: &mut C| {
            if k.is_lost(qubit) {
                return None;
            }
            if k.x_bit(qubit) {
                *c *= C::zero();
                return None;
            }
            *c *= half.clone();
            Some((k.toggled_bits(qubit, false, true), -c.clone()))
        });
    }
}
