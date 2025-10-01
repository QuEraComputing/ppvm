use crate::char::Pauli;
use crate::config::Config;

pub trait RotationOne<T: Config> {
    fn rotate_1(&mut self, axis: Pauli, addr0: usize, theta: T::Value);
    fn rx(&mut self, addr0: usize, theta: T::Value) {
        self.rotate_1(Pauli::X, addr0, theta)
    }
    fn ry(&mut self, addr0: usize, theta: T::Value) {
        self.rotate_1(Pauli::Y, addr0, theta)
    }
    fn rz(&mut self, addr0: usize, theta: T::Value) {
        self.rotate_1(Pauli::Z, addr0, theta)
    }
}
