use crate::config::Config;

pub trait CRx<T: Config> {
    fn crx(&mut self, control: usize, target: usize, theta: T::Coeff);
}
