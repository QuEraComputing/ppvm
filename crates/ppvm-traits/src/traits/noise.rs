// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::config::Config;
use num::Zero;

// FIXME: most channels don't need to own probs, we can just reference them and clean up the code

/// Single-qubit Pauli error channel — apply `X`, `Y`, or `Z` with the
/// three given probabilities.
pub trait PauliError<T: Config> {
    /// Apply a Pauli-error channel `[p_x, p_y, p_z]` to one qubit.
    fn pauli_error(&mut self, addr0: usize, p: [T::Coeff; 3]);

    /// stim `X_ERROR(p)` — apply X with probability `p` to one qubit.
    fn x_error(&mut self, addr0: usize, p: T::Coeff) {
        let zero = T::Coeff::zero();
        self.pauli_error(addr0, [p, zero.clone(), zero])
    }

    /// stim `Y_ERROR(p)` — apply Y with probability `p` to one qubit.
    fn y_error(&mut self, addr0: usize, p: T::Coeff) {
        let zero = T::Coeff::zero();
        self.pauli_error(addr0, [zero.clone(), p, zero])
    }

    /// stim `Z_ERROR(p)` — apply Z with probability `p` to one qubit.
    fn z_error(&mut self, addr0: usize, p: T::Coeff) {
        let zero = T::Coeff::zero();
        self.pauli_error(addr0, [zero.clone(), zero, p])
    }

    /// Explicit batched Pauli-error channel.
    fn pauli_error_many(&mut self, targets: &[usize], p: [T::Coeff; 3]) {
        for &q in targets {
            self.pauli_error(q, p.clone());
        }
    }

    /// Explicit batched `X_ERROR(p)`.
    fn x_error_many(&mut self, targets: &[usize], p: T::Coeff) {
        for &q in targets {
            self.x_error(q, p.clone());
        }
    }

    /// Explicit batched `Y_ERROR(p)`.
    fn y_error_many(&mut self, targets: &[usize], p: T::Coeff) {
        for &q in targets {
            self.y_error(q, p.clone());
        }
    }

    /// Explicit batched `Z_ERROR(p)`.
    fn z_error_many(&mut self, targets: &[usize], p: T::Coeff) {
        for &q in targets {
            self.z_error(q, p.clone());
        }
    }
}

/// Two-qubit Pauli error channel.
pub trait TwoQubitPauliError<T: Config> {
    /// Apply a two-qubit Pauli-error channel to one pair. Probabilities are given in the order:
    /// `{IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ}`.
    fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [T::Coeff; 15]);

    /// Explicit batched two-qubit Pauli-error channel.
    fn two_qubit_pauli_error_many(&mut self, pairs: &[(usize, usize)], p: [T::Coeff; 15]) {
        for &(a, b) in pairs {
            self.two_qubit_pauli_error(a, b, p.clone());
        }
    }
}

/// Single-qubit depolarizing channel.
pub trait Depolarizing<T: Config> {
    /// Depolarize one qubit with probability `p`.
    fn depolarize1(&mut self, addr0: usize, p: T::Coeff);

    /// Explicit batched single-qubit depolarizing channel.
    fn depolarize1_many(&mut self, targets: &[usize], p: T::Coeff) {
        for &q in targets {
            self.depolarize1(q, p.clone());
        }
    }
}

/// Two-qubit depolarizing channel.
pub trait Depolarizing2<T: Config> {
    /// Depolarize one qubit pair with probability `p`.
    fn depolarize2(&mut self, addr0: usize, addr1: usize, p: T::Coeff);

    /// Explicit batched two-qubit depolarizing channel.
    fn depolarize2_many(&mut self, pairs: &[(usize, usize)], p: T::Coeff) {
        for &(a, b) in pairs {
            self.depolarize2(a, b, p.clone());
        }
    }
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
///
/// # The `p[1]` convention (normative)
///
/// This is the one place the parameterization is defined; every backend
/// (`ppvm-pauli-sum`, `ppvm-tableau`'s trajectory, `ppvm-tableau-sum`'s mixture)
/// and every binding cites it rather than restating it. In the paper's notation
/// `p = [p_LL, p_LQ, p_LN]`, and `p[1] = p_LQ` is the probability that a
/// **named** one of the two atoms is lost while the other survives. The two
/// single-loss events are disjoint, so
///
/// * the probability of losing *exactly one* atom is `2·p[1]`, and
/// * the probability that both remain in the qubit subspace — the factor a
///   fully in-subspace observable is scaled by — is `1 − 2·p[1] − p[0]`.
///
/// The channel is completely positive exactly on `p[0], p[1] >= 0`,
/// `p[0] + 2·p[1] <= 1`, `p[2] ∈ [0, 1]`; the tableau backends `debug_assert`
/// that region.
pub trait CorrelatedLossChannel<T: Config> {
    /// Apply a correlated loss channel to qubits at `addr0` and `addr1`.
    ///
    /// The three probabilities are:
    /// * `p[0]`: The probability of losing both qubits simultaneously when
    ///   both of them are in the qubit subspace.
    /// * `p[1]`: The probability of losing a **named** one of the two qubits when
    ///   both of them are in the qubit subspace, so losing *exactly one* has
    ///   probability `2·p[1]` and the both-present survivor is scaled by
    ///   `1 − 2·p[1] − p[0]` (which qubit is lost is 50/50).
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

/// Single-qubit leakage channel — with probability `p0` the qubit is
/// leaked and pinned to `|0⟩`, with probability `p1` leaked and pinned
/// to `|1⟩`. Unlike [`LossChannel`], a leaked qubit remains measurable
/// and returns the pinned computational bit; gates still skip it.
pub trait LeakageChannel<T: Config> {
    /// Apply leakage to qubit `addr0`, with `p0` / `p1` the probabilities
    /// of leaking into `|0⟩` / `|1⟩`.
    fn leakage_channel(&mut self, addr0: usize, p0: T::Coeff, p1: T::Coeff);
}

/// Reset leakage on a qubit — used to model a leakage-reduction event
/// that returns a previously-leaked qubit to the computational subspace
/// in `|0⟩`.
pub trait ResetLeakageChannel<T: Config> {
    /// Clear leakage at `addr0` and reset the qubit to `|0⟩`.
    ///
    /// No-op if the qubit is live. A lost qubit cannot be recovered this way.
    fn reset_leakage_channel(&mut self, addr0: usize);
}
