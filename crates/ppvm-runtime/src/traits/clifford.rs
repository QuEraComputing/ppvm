pub trait Clifford {
    fn x(&mut self, index: usize);
    fn y(&mut self, index: usize);
    fn z(&mut self, index: usize);
    fn h(&mut self, index: usize);
    fn s(&mut self, index: usize);
    fn s_adj(&mut self, index: usize);
    fn cnot(&mut self, control: usize, target: usize);
    fn cz(&mut self, control: usize, target: usize);
}

pub trait CliffordExtensions: Clifford {
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
        self.s(addr0);
        self.sqrt_x_adj(addr0);
        self.s_adj(addr0);
    }
}

impl<C> CliffordExtensions for C where C: Clifford {}
