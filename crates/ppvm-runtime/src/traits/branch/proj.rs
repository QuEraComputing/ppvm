/// Projective Z-basis projectors `|0⟩⟨0|` and `|1⟩⟨1|`.
pub trait Projection {
    /// Project qubit `pos` onto `|0⟩`.
    fn p0(&mut self, pos: usize);
    /// Project qubit `pos` onto `|1⟩`.
    fn p1(&mut self, pos: usize);
}
