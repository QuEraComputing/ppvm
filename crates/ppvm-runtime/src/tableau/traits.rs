pub trait GeneralizedTableauTGate {
    fn t(&mut self, addr0: usize);
    fn t_adj(&mut self, addr0: usize);
    fn t_or_t_adj(&mut self, addr0: usize, adjoint: bool);
    fn compute_shift_z(&self, addr0: usize) -> usize;
    fn compute_phase_z(&self, addr0: usize, branch_index: usize) -> u8;
}

pub trait Measure {
    fn measure(&mut self, addr0: usize) -> bool;
    fn find_anticommuting_stabilizer(&self, addr0: usize) -> Option<usize>;
    fn update_tableau_according_to_outcome(&mut self, addr0: usize, q_idx: usize, outcome: bool);
    fn get_deterministic_outcome(&self, addr0: usize) -> bool;
}
