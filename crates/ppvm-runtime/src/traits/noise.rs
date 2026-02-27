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

pub trait DepolarizingError<T: Config> {
    fn depolarize(&mut self, addr0: usize, p: T::Coeff);
}

pub trait Depolarizing2Error<T: Config> {
    fn depolarize2(&mut self, addr0: usize, addr1: usize, p: T::Coeff);
}

impl<C, T: Config> Depolarizing2Error<T> for C
where
    C: TwoQubitPauliError<T>,
{
    fn depolarize2(&mut self, addr0: usize, addr1: usize, p: T::Coeff) {
        let ps: [T::Coeff; 15] = core::array::from_fn(|_| p.clone() * (1.0 / 15.0));
        self.two_qubit_pauli_error(addr0, addr1, ps);
    }
}

pub trait AmplitudeDamping<T: Config> {
    fn amplitude_damping(&mut self, addr0: usize, gamma: T::Coeff);
}

pub trait LossChannel<T: Config> {
    fn loss_channel(&mut self, addr0: usize, p: T::Coeff);
}

pub trait ResetLossChannel<T: Config> {
    fn reset_loss_channel(&mut self, addr0: usize);
}
