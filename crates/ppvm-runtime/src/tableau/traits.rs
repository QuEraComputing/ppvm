use crate::config::Config;
use crate::traits::Clifford;
use num::complex::Complex;

pub trait TGate<T: Config> {
    fn t(&mut self, addr0: usize);
    fn t_adj(&mut self, addr0: usize);
    fn t_or_t_adj(&mut self, addr0: usize, adjoint: bool);
    fn rz(&mut self, addr0: usize, theta: T::Coeff);
    fn branch_z_with_coefficients(
        &mut self,
        addr0: usize,
        complex_cos: Complex<T::Coeff>,
        complex_sin: Complex<T::Coeff>,
    );
}

pub trait Measure {
    fn measure(&mut self, addr0: usize) -> bool;
}

pub trait CliffordExtensions: Clifford {
    fn s_adj(&mut self, addr0: usize);
    fn sqrt_x(&mut self, addr0: usize) {
        self.h(addr0);
        self.s(addr0);
        self.h(addr0);
    }

    fn sqrt_x_adj(&mut self, addr0: usize) {
        self.h(addr0);
        self.s_adj(addr0);
        self.h(addr0);
    }

    fn sqrt_y(&mut self, addr0: usize) {
        self.s(addr0);
        self.sqrt_x(addr0);
        self.s_adj(addr0);
    }

    fn sqrt_y_adj(&mut self, addr0: usize) {
        self.s_adj(addr0);
        self.sqrt_x_adj(addr0);
        self.s(addr0);
    }
}
