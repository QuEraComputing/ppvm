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

use ppvm_pauli_word_2::{HashFinalize, PauliStorage, PauliWord};
use ppvm_traits_2::{Accumulate, Coefficient, PauliBits, PauliError};
use std::hash::BuildHasher;

use crate::store::ScaleByKey;
use crate::sum::Sum;

/// `1 − 2(a + b)`, formed without any bare-`f64` scaling: the doubling is the
/// ring addition `x + x`.
#[inline]
fn eigenvalue<C: Coefficient + num::One>(a: &C, b: &C) -> C {
    C::one() - (a.clone() + a.clone()) - (b.clone() + b.clone())
}

/// Diagonal unital Pauli-channel propagation on a Pauli-keyed `Sum`: an in-place
/// per-term scale by the channel's real transfer eigenvalue.
impl<S, P, A, H, C> PauliError<C> for Sum<S, P>
where
    S: Accumulate<Key = PauliWord<A, H>, Coeff = C> + ScaleByKey<PauliWord<A, H>, C>,
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
    C: Coefficient + num::One,
    P: crate::policy::Policy<PauliWord<A, H>, C>,
{
    /// Scale each term in place by `λ_P` for its Pauli at `qubit`. The three
    /// non-trivial eigenvalues are computed once; the identity term (and any lost
    /// qubit) is left exactly untouched.
    ///
    /// Mirrors the old `ppvm-pauli-sum::sum::noise`'s `self.scale(|k, v| ...)`
    /// verbatim in shape as well as arithmetic: an in-place mutation that cannot
    /// remove, so a zero eigenvalue leaves a zero-coefficient term in the support
    /// exactly as old does.
    #[inline]
    fn pauli_error(&mut self, qubit: usize, probabilities: [C; 3]) {
        let [px, py, pz] = probabilities;
        // λ_X = 1 − 2(pY + pZ), λ_Z = 1 − 2(pX + pY), λ_Y = 1 − 2(pX + pZ).
        let x_factor = eigenvalue(&py, &pz);
        let z_factor = eigenvalue(&px, &py);
        let y_factor = eigenvalue(&px, &pz);
        self.scale_by_key(move |k: &PauliWord<A, H>, c: &mut C| {
            if k.is_lost(qubit) {
                return;
            }
            // 2-bit Pauli code (x, z): (0,0) I, (1,0) X, (0,1) Z, (1,1) Y.
            match (k.x_bit(qubit), k.z_bit(qubit)) {
                (false, false) => {} // I: λ_I = 1, an exact no-op.
                (true, false) => *c *= x_factor.clone(),
                (false, true) => *c *= z_factor.clone(),
                (true, true) => *c *= y_factor.clone(),
            }
        });
    }
}
