use crate::traits::Clifford;

pub trait TGate {
    fn t(&mut self, addr0: usize);
    fn t_adj(&mut self, addr0: usize);
    fn t_or_t_adj(&mut self, addr0: usize, adjoint: bool);
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
}
