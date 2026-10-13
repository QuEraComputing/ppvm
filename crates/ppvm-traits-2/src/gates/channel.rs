// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Noise-channel capabilities. `rng` is threaded through each call rather than
//! owned by the state: trajectory backends draw from it, deterministic ones
//! ignore it, and the caller owns seeding and stream sharing.

use crate::arithmetic::Coefficient;

/// A unital single-qubit Pauli error channel `P ↦ λ_P·P`.
pub trait PauliError<C: Coefficient> {
    /// Apply a single-qubit Pauli channel with `X`, `Y`, `Z` probabilities.
    fn pauli_error(&mut self, qubit: usize, probabilities: [C; 3]);

    /// stim `X_ERROR(p)` — apply `X` with probability `p` to one qubit.
    fn x_error(&mut self, qubit: usize, p: C) {
        self.pauli_error(qubit, [p, C::zero(), C::zero()])
    }

    /// stim `Y_ERROR(p)` — apply `Y` with probability `p` to one qubit.
    fn y_error(&mut self, qubit: usize, p: C) {
        self.pauli_error(qubit, [C::zero(), p, C::zero()])
    }

    /// stim `Z_ERROR(p)` — apply `Z` with probability `p` to one qubit.
    fn z_error(&mut self, qubit: usize, p: C) {
        self.pauli_error(qubit, [C::zero(), C::zero(), p])
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
/// (see [`LossState`](crate::loss::LossState)).
pub trait LossChannel<C: Coefficient> {
    /// Apply a loss channel to `qubit` with loss probability `p`.
    fn loss_channel(&mut self, qubit: usize, p: C);
}

/// Correlated two-qubit loss channel. Completely positive exactly on
/// `p[0], p[1] >= 0`, `p[0] + 2·p[1] <= 1`, `p[2] ∈ [0, 1]`.
pub trait CorrelatedLossChannel<C: Coefficient> {
    /// Apply correlated loss `p = [p_LL, p_LQ, p_LN]`: `p[0]` loses both, `p[2]`
    /// the survivor of an earlier loss. `p[1]` loses a **named** one, so exactly
    /// one goes with probability `2·p[1]` and both remain with `1−2·p[1]−p[0]`.
    fn correlated_loss_channel(&mut self, qubit0: usize, qubit1: usize, p: [C; 3]);
}

/// Reset the loss bit on a qubit — models a re-cooling / re-loading event that
/// brings a previously-lost atom back.
pub trait ResetLossChannel {
    /// Clear the loss bit at `qubit`.
    fn reset_loss_channel(&mut self, qubit: usize);
}

/// State-dependent loss: probability `p0` from `|0⟩`, `p1` from `|1⟩`.
/// The total loss probability depends on the populations, so the channel reads `⟨Z⟩`.
pub trait AsymmetricLossChannel<C: Coefficient> {
    /// Apply asymmetric loss to `qubit`, with `p0` / `p1` the loss probabilities
    /// from `|0⟩` / `|1⟩`. See the backend impl for the trajectory approximation
    /// used (the survival back-action is omitted).
    fn asymmetric_loss_channel(&mut self, qubit: usize, p0: C, p1: C);
}
