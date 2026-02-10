pub trait GeneralizedTableauTGate {
    fn t(&mut self, addr0: usize);
    fn t_adj(&mut self, addr0: usize);
    fn t_or_t_adj(&mut self, addr0: usize, adjoint: bool);
    fn compute_shift_z(&self, addr0: usize) -> usize;
    fn compute_phase_z(&self, addr0: usize, branch_index: usize) -> u8;
}

pub trait Measure {
    fn measure(&mut self, addr0: usize) -> bool;
}
