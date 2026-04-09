pub trait Reset {
    fn reset(&mut self, addr0: usize);
}
