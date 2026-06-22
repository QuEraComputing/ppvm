// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::config::Config;

// FIXME: most channels don't need to own probs, we can just reference them and clean up the code

/// Single-qubit Pauli error channel — apply `X`, `Y`, or `Z` with the
/// three given probabilities.
pub trait PauliError<T: Config> {
    /// Apply a Pauli-error channel to qubit `addr0` with the probability
    /// triple `p = [p_x, p_y, p_z]`.
    fn pauli_error(&mut self, addr0: usize, p: [T::Coeff; 3]);
}

/// Apply the same single-qubit Pauli error channel uniformly to every
/// qubit in the system.
pub trait PauliErrorAll<T: Config> {
    /// Apply `Pauli error` with probabilities `p = [p_x, p_y, p_z]` to
    /// every qubit.
    fn pauli_error_all(&mut self, p: [T::Coeff; 3]);
}

/// Two-qubit Pauli error channel.
pub trait TwoQubitPauliError<T: Config> {
    /// Apply a two-qubit Pauli-error channel to `(addr0, addr1)`.
    /// Probabilities are given in the order:
    /// `{IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ}`.
    fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [T::Coeff; 15]);
}

/// Single-qubit depolarizing channel.
pub trait Depolarizing<T: Config> {
    /// Depolarize qubit `addr0` with probability `p`.
    fn depolarize(&mut self, addr0: usize, p: T::Coeff);
}

/// Two-qubit depolarizing channel.
pub trait Depolarizing2<T: Config> {
    /// Depolarize the pair `(addr0, addr1)` with probability `p`.
    fn depolarize2(&mut self, addr0: usize, addr1: usize, p: T::Coeff);
}

/// Amplitude-damping channel (single qubit).
pub trait AmplitudeDamping<T: Config> {
    /// Apply amplitude damping with damping parameter `gamma`.
    fn amplitude_damping(&mut self, addr0: usize, gamma: T::Coeff);
}

/// Single-qubit loss channel — with probability `p`, mark the qubit as
/// lost (`Pauli::L`).
pub trait LossChannel<T: Config> {
    /// Apply a loss channel to qubit `addr0` with loss probability `p`.
    fn loss_channel(&mut self, addr0: usize, p: T::Coeff);
}

/// Correlated two-qubit loss channel.
pub trait CorrelatedLossChannel<T: Config> {
    /// Apply a correlated loss channel to qubits at `addr0` and `addr1`.
    ///
    /// The three probabilities are:
    /// * `p[0]`: The probability of losing both qubits simultaneously when
    ///   both of them are in the qubit subspace.
    /// * `p[1]`: The probability of losing either one qubit when both of them are
    ///   in the qubit subspace.
    /// * `p[2]`: The probability of losing one qubit when the other one has already
    ///   been lost prior to the channel.
    fn correlated_loss_channel(&mut self, addr0: usize, addr1: usize, p: [T::Coeff; 3]);
}

/// Reset the loss bit on a qubit — used to model a re-cooling /
/// re-loading event that brings a previously-lost atom back.
pub trait ResetLossChannel<T: Config> {
    /// Clear the loss bit at `addr0`.
    fn reset_loss_channel(&mut self, addr0: usize);
}

/// State-dependent ("asymmetric") single-qubit loss channel: a qubit is
/// lost from `|0⟩` with probability `p0` and from `|1⟩` with probability
/// `p1`. Unlike [`LossChannel`], the total loss probability depends on the
/// qubit's populations, so the channel reads the current `⟨Z⟩`.
pub trait AsymmetricLossChannel<T: Config> {
    /// Apply asymmetric loss to qubit `addr0`, with `p0` / `p1` the loss
    /// probabilities from `|0⟩` / `|1⟩`. See the backend impl for the
    /// trajectory approximation used (the survival back-action is omitted).
    fn asymmetric_loss_channel(&mut self, addr0: usize, p0: T::Coeff, p1: T::Coeff);
}
