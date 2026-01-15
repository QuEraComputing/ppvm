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
    /// A simple loss channel that reduces the coefficients of all terms by (1 - p)
    /// This is equivalent to multiplying the density matrix by (1 - p), thereby reducing
    /// the trace of the density matrix.
    fn loss_channel(&mut self, addr0: usize, p: T::Coeff);
}
