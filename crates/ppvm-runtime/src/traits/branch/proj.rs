pub trait Projection {
    fn p0(&mut self, pos: usize);
    fn p1(&mut self, pos: usize);
}
