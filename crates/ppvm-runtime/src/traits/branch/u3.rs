use crate::config::Config;

pub trait U3Gate<T: Config> {
    fn u3(&mut self, addr: usize, theta: T::Coeff, phi: T::Coeff, lambda: T::Coeff);
}
