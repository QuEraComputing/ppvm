use crate::config::Config;

macro_rules! def_rotation {
    ($name:ident, $x_a:expr, $z_a:expr, $x_b:expr, $z_b:expr) => {
        fn $name(&mut self, a: usize, b: usize, theta: T::Value) {
            self.rotate_2($x_a, $z_a, $x_b, $z_b, a, b, theta)
        }
    };
}

pub trait RotationTwo<T: Config> {
    fn rotate_2(
        &mut self,
        axis_a_x: u8,
        axis_a_z: u8,
        axis_b_x: u8,
        axis_b_z: u8,
        a: usize,
        b: usize,
        theta: T::Value,
    );
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
