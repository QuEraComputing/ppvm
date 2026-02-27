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
