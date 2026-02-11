pub trait TGate {
    fn t(&mut self, addr0: usize);
    fn t_adj(&mut self, addr0: usize);
    fn t_or_t_adj(&mut self, addr0: usize, adjoint: bool);
}

pub trait Measure {
    fn measure(&mut self, addr0: usize) -> bool;
}
