pub trait Measure {
    fn measure(&mut self, addr0: usize) -> bool;
}
