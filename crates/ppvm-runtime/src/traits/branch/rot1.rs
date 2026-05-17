use crate::char::Pauli;
use crate::config::Config;

/// Single-qubit Pauli rotations `exp(-i θ/2 · P)`.
pub trait RotationOne<T: Config> {
    /// Rotate about `axis` (one of `X`, `Y`, `Z`) by angle `theta`.
    fn rotate_1(&mut self, axis: Pauli, addr0: usize, theta: T::Coeff);
    /// `RX(θ)` on qubit `addr0`.
    fn rx(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        self.rotate_1(Pauli::X, addr0, theta.into())
    }
    /// `RY(θ)` on qubit `addr0`.
    fn ry(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        self.rotate_1(Pauli::Y, addr0, theta.into())
    }
    /// `RZ(θ)` on qubit `addr0`.
    fn rz(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        self.rotate_1(Pauli::Z, addr0, theta.into())
    }
}
