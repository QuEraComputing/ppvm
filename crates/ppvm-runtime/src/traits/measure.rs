pub trait Measure {
    fn measure(&mut self, addr0: usize) -> bool;
}

pub trait LossyMeasure {
    fn measure(&mut self, addr0: usize) -> Option<bool>;
}
