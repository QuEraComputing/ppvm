use crate::char::Pauli;
use crate::config::Config;

pub trait RotationOneMapInsertClosure<T: Config> {
    fn rotate_1_map_insert_closure(
        k: &T::PauliWordType,
        v: &mut T::Coeff,
        axis: Pauli,
        addr0: usize,
        sin: &T::Coeff,
        cos: &T::Coeff,
    ) -> Option<(T::PauliWordType, T::Coeff)>;
}

pub trait RotationOne<T: Config>: RotationOneMapInsertClosure<T> {
    fn rotate_1(&mut self, axis: Pauli, addr0: usize, theta: T::Coeff);
    fn rx(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        self.rotate_1(Pauli::X, addr0, theta.into())
    }
    fn ry(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        self.rotate_1(Pauli::Y, addr0, theta.into())
    }
    fn rz(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        self.rotate_1(Pauli::Z, addr0, theta.into())
    }
}
