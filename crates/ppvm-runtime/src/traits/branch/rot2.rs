use crate::config::Config;

macro_rules! def_rotation {
    ($name:ident, $x_a:expr, $z_a:expr, $x_b:expr, $z_b:expr) => {
        fn $name(&mut self, a: usize, b: usize, theta: impl Into<T::Coeff>) {
            self.rotate_2([$x_a, $z_a], [$x_b, $z_b], a, b, theta.into())
        }
    };
}

pub trait RotationTwo<T: Config> {
    /// Two-qubit Pauli rotation: `exp(-i * theta/2 * P_a ⊗ P_b)`.
    ///
    /// Each axis is encoded as `[x, z]` bits:
    /// `[0,0]` = I, `[1,0]` = X, `[0,1]` = Z, `[1,1]` = Y.
    fn rotate_2(&mut self, axis_a: [u8; 2], axis_b: [u8; 2], a: usize, b: usize, theta: T::Coeff);
    //                 x, z, x, z
    def_rotation!(rxx, 1, 0, 1, 0);
    def_rotation!(rxy, 1, 0, 1, 1);
    def_rotation!(rxz, 1, 0, 0, 1);

    def_rotation!(ryx, 1, 1, 1, 0);
    def_rotation!(ryy, 1, 1, 1, 1);
    def_rotation!(ryz, 1, 1, 0, 1);

    def_rotation!(rzx, 0, 1, 1, 0);
    def_rotation!(rzy, 0, 1, 1, 1);
    def_rotation!(rzz, 0, 1, 0, 1);
}
