use crate::config::Config;

pub trait TGate<T: Config> {
    fn t(&mut self, addr0: usize);
    fn t_adj(&mut self, addr0: usize);
}
