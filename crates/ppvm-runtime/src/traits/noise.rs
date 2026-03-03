use crate::config::Config;

pub trait PauliError<T: Config> {
    fn pauli_error(&mut self, addr0: usize, p: [T::Coeff; 3]);
}

pub trait PauliErrorAll<T: Config> {
    fn pauli_error_all(&mut self, p: [T::Coeff; 3]);
}

pub trait TwoQubitPauliError<T: Config> {
    /// Probabilities are given in the order:
    /// {IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ}
    fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [T::Coeff; 15]);
}

pub trait Depolarizing<T: Config> {
    fn depolarizing(&mut self, addr0: usize, p: T::Coeff);
}

pub trait AmplitudeDamping<T: Config> {
    fn amplitude_damping(&mut self, addr0: usize, gamma: T::Coeff);
}

pub trait LossChannel<T: Config> {
    fn loss_channel(&mut self, addr0: usize, p: T::Coeff);
}

pub trait CorrelatedLossChannel<T: Config> {
    /// Apply a correlated loss channel to qubits at `addr0` and `addr1`.
    ///
    /// The three probabilities are:
    /// * `p[0]`: The probability of losing both qubits simultaneously when
    ///     both of them are in the qubit subspace.
    /// * `p[1]`: The probability of losing either one qubit when both of them are
    ///     in the qubit subspace.
    /// * `p[2]`: The probability of losing one qubit when the other one has already
    ///     been lost prior to the channel.
    fn correlated_loss_channel(&mut self, addr0: usize, addr1: usize, p: [T::Coeff; 3]);
}

pub trait ResetLossChannel<T: Config> {
    fn reset_loss_channel(&mut self, addr0: usize);
}
