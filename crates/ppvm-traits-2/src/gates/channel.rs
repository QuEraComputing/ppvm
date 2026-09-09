// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::arithmetic::Coefficient;

/// Coefficients that can calculate single-qubit Pauli-channel factors.
///
/// This optional channel capability keeps noise formulas out of scalar
/// arithmetic. Implementers may use the generic default or specialize it.
pub trait PauliErrorFactors: Coefficient + num::One {
    /// Transfer eigenvalues `(λ_X, λ_Z, λ_Y)` for Pauli probabilities
    /// `(p_X, p_Y, p_Z)`.
    #[inline(always)]
    fn pauli_error_factors(probabilities: [Self; 3]) -> [Self; 3] {
        let [px, py, pz] = probabilities;
        let one = Self::one();
        [
            one.clone() - py.doubled() - pz.doubled(),
            one.clone() - px.doubled() - py.doubled(),
            one - px.doubled() - pz.doubled(),
        ]
    }
}

impl PauliErrorFactors for f64 {
    #[inline(always)]
    fn pauli_error_factors([px, py, pz]: [Self; 3]) -> [Self; 3] {
        [
            1.0 - py * 2.0 - pz * 2.0,
            1.0 - px * 2.0 - py * 2.0,
            1.0 - px * 2.0 - pz * 2.0,
        ]
    }
}

impl PauliErrorFactors for num::Complex<f64> {
    #[inline(always)]
    fn pauli_error_factors([px, py, pz]: [Self; 3]) -> [Self; 3] {
        let one = Self::new(1.0, 0.0);
        [
            one - py * 2.0 - pz * 2.0,
            one - px * 2.0 - py * 2.0,
            one - px * 2.0 - pz * 2.0,
        ]
    }
}

/// A unital single-qubit Pauli error channel `P ↦ λ_P·P`.
pub trait PauliError<C: Coefficient> {
    /// Apply a single-qubit Pauli channel with `X`, `Y`, `Z` probabilities.
    fn pauli_error<R: rand::Rng + ?Sized>(
        &mut self,
        qubit: usize,
        probabilities: [C; 3],
        rng: &mut R,
    );

    /// stim `X_ERROR(p)` — apply `X` with probability `p` to one qubit.
    fn x_error<R: rand::Rng + ?Sized>(&mut self, qubit: usize, p: C, rng: &mut R) {
        let zero = C::zero();
        self.pauli_error(qubit, [p, zero.clone(), zero], rng)
    }

    /// stim `Y_ERROR(p)` — apply `Y` with probability `p` to one qubit.
    fn y_error<R: rand::Rng + ?Sized>(&mut self, qubit: usize, p: C, rng: &mut R) {
        let zero = C::zero();
        self.pauli_error(qubit, [zero.clone(), p, zero], rng)
    }

    /// stim `Z_ERROR(p)` — apply `Z` with probability `p` to one qubit.
    fn z_error<R: rand::Rng + ?Sized>(&mut self, qubit: usize, p: C, rng: &mut R) {
        let zero = C::zero();
        self.pauli_error(qubit, [zero.clone(), zero, p], rng)
    }

    /// Explicit batched Pauli-error channel.
    fn pauli_error_many<R: rand::Rng + ?Sized>(
        &mut self,
        targets: &[usize],
        p: [C; 3],
        rng: &mut R,
    ) {
        for &q in targets {
            self.pauli_error(q, p.clone(), rng);
        }
    }

    /// Explicit batched `X_ERROR(p)`.
    fn x_error_many<R: rand::Rng + ?Sized>(&mut self, targets: &[usize], p: C, rng: &mut R) {
        for &q in targets {
            self.x_error(q, p.clone(), rng);
        }
    }

    /// Explicit batched `Y_ERROR(p)`.
    fn y_error_many<R: rand::Rng + ?Sized>(&mut self, targets: &[usize], p: C, rng: &mut R) {
        for &q in targets {
            self.y_error(q, p.clone(), rng);
        }
    }

    /// Explicit batched `Z_ERROR(p)`.
    fn z_error_many<R: rand::Rng + ?Sized>(&mut self, targets: &[usize], p: C, rng: &mut R) {
        for &q in targets {
            self.z_error(q, p.clone(), rng);
        }
    }
}

/// Apply the same single-qubit Pauli error channel uniformly to every qubit in
/// the system.
pub trait PauliErrorAll<C: Coefficient> {
    /// Apply the Pauli channel `p = [p_x, p_y, p_z]` to every qubit.
    fn pauli_error_all<R: rand::Rng + ?Sized>(&mut self, p: [C; 3], rng: &mut R);
}

/// Two-qubit Pauli error channel.
pub trait TwoQubitPauliError<C: Coefficient> {
    /// Apply a two-qubit Pauli-error channel to one pair. Probabilities are given
    /// in the order
    /// `{IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ}`.
    fn two_qubit_pauli_error<R: rand::Rng + ?Sized>(
        &mut self,
        qubit0: usize,
        qubit1: usize,
        p: [C; 15],
        rng: &mut R,
    );

    /// Explicit batched two-qubit Pauli-error channel.
    fn two_qubit_pauli_error_many<R: rand::Rng + ?Sized>(
        &mut self,
        pairs: &[(usize, usize)],
        p: [C; 15],
        rng: &mut R,
    ) {
        for &(a, b) in pairs {
            self.two_qubit_pauli_error(a, b, p.clone(), rng);
        }
    }
}

/// Single-qubit depolarizing channel.
pub trait Depolarizing<C: Coefficient> {
    /// Depolarize one qubit with probability `p`.
    fn depolarize1<R: rand::Rng + ?Sized>(&mut self, qubit: usize, p: C, rng: &mut R);

    /// Explicit batched single-qubit depolarizing channel.
    fn depolarize1_many<R: rand::Rng + ?Sized>(&mut self, targets: &[usize], p: C, rng: &mut R) {
        for &q in targets {
            self.depolarize1(q, p.clone(), rng);
        }
    }
}

/// Two-qubit depolarizing channel.
pub trait Depolarizing2<C: Coefficient> {
    /// Depolarize one qubit pair with probability `p`.
    fn depolarize2<R: rand::Rng + ?Sized>(
        &mut self,
        qubit0: usize,
        qubit1: usize,
        p: C,
        rng: &mut R,
    );

    /// Explicit batched two-qubit depolarizing channel.
    fn depolarize2_many<R: rand::Rng + ?Sized>(
        &mut self,
        pairs: &[(usize, usize)],
        p: C,
        rng: &mut R,
    ) {
        for &(a, b) in pairs {
            self.depolarize2(a, b, p.clone(), rng);
        }
    }
}

/// Amplitude-damping channel (single qubit).
pub trait AmplitudeDamping<C: Coefficient> {
    /// Apply amplitude damping with damping parameter `gamma`.
    fn amplitude_damping(&mut self, qubit: usize, gamma: C);
}

/// Single-qubit loss channel — with probability `p`, mark the qubit as lost
/// (see [`LossState`](crate::loss::LossState)).
pub trait LossChannel<C: Coefficient> {
    /// Apply a loss channel to `qubit` with loss probability `p`.
    fn loss_channel<R: rand::Rng + ?Sized>(&mut self, qubit: usize, p: C, rng: &mut R);
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
    fn correlated_loss_channel<R: rand::Rng + ?Sized>(
        &mut self,
        qubit0: usize,
        qubit1: usize,
        p: [C; 3],
        rng: &mut R,
    );
}

/// Reset the loss bit on a qubit — models a re-cooling / re-loading event that
/// brings a previously-lost atom back.
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
    fn asymmetric_loss_channel<R: rand::Rng + ?Sized>(
        &mut self,
        qubit: usize,
        p0: C,
        p1: C,
        rng: &mut R,
    );
}

#[cfg(test)]
mod tests {
    use super::PauliErrorFactors;
    use num::Complex;

    #[test]
    fn deterministic_pauli_errors_have_expected_conjugation_signs() {
        // Inputs are X/Y/Z probabilities; outputs are X/Z/Y eigenvalues.
        for (probabilities, expected) in [
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            ([1.0, 0.0, 0.0], [1.0, -1.0, -1.0]),
            ([0.0, 1.0, 0.0], [-1.0, -1.0, 1.0]),
            ([0.0, 0.0, 1.0], [-1.0, 1.0, -1.0]),
        ] {
            assert_eq!(f64::pauli_error_factors(probabilities), expected);
            assert_eq!(
                Complex::<f64>::pauli_error_factors(probabilities.map(|p| Complex::new(p, 0.0)),),
                expected.map(|factor| Complex::new(factor, 0.0)),
            );
        }
    }
}
