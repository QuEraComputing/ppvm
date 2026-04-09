pub trait Clifford {
    fn x(&mut self, index: usize);
    fn y(&mut self, index: usize);
    fn z(&mut self, index: usize);
    fn h(&mut self, index: usize);
    fn s(&mut self, index: usize);
    fn cnot(&mut self, control: usize, target: usize);
    fn cz(&mut self, control: usize, target: usize);
}

pub trait CliffordExtensions: Clifford {
    fn s_adj(&mut self, addr0: usize);
    fn sqrt_x(&mut self, addr0: usize);
    fn sqrt_x_adj(&mut self, addr0: usize);
    fn sqrt_y(&mut self, addr0: usize);
    fn sqrt_y_adj(&mut self, addr0: usize);

    fn cy(&mut self, addr0: usize, addr1: usize);
}
