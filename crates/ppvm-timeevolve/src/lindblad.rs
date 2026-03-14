use ppvm_runtime::prelude::{Config, PhasedPauliWord};

pub enum RateMatrix {
    Vector(Vec<f64>),
    Dense(Vec<Vec<f64>>),
}

impl From<Vec<f64>> for RateMatrix {
    fn from(v: Vec<f64>) -> Self {
        RateMatrix::Vector(v)
    }
}

pub struct CollapseOp<T: Config> {
    #[allow(clippy::type_complexity)]
    pub(crate) terms: Vec<(PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>, f64)>,
    #[allow(dead_code)] // used in LindbladOp::new (Task 3)
    pub(crate) n_qubits: usize,
}

impl<T: Config> CollapseOp<T> {
    pub fn new(n_qubits: usize) -> Self {
        CollapseOp { terms: Vec::new(), n_qubits }
    }

    pub fn push(
        &mut self,
        word: PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>,
        coeff: f64,
    ) {
        self.terms.push((word, coeff));
    }
}

pub struct LindbladOp<T: Config> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: Config> LindbladOp<T> {
    pub fn new(_ops: Vec<CollapseOp<T>>, _rates: RateMatrix) -> Self {
        LindbladOp { _phantom: std::marker::PhantomData }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::prelude::{PauliWord, PhasedPauliWord, config::fxhash::ByteF64};

    type W1 = PauliWord<[u8; 1], fxhash::FxBuildHasher>;
    type PPW1 = PhasedPauliWord<[u8; 1], fxhash::FxBuildHasher, W1>;

    fn ppw(pauli: &str, phase: u8) -> PPW1 {
        PhasedPauliWord::build_from_word(W1::from(pauli), phase)
    }

    #[test]
    fn collapse_op_x_plus_iy() {
        // c = X + iY: X has phase 0 (+1), Y has phase 1 (+i)
        let mut op = CollapseOp::<ByteF64<1>>::new(1);
        op.push(ppw("X", 0), 1.0);
        op.push(ppw("Y", 1), 1.0);
        assert_eq!(op.terms.len(), 2);
        assert_eq!(op.n_qubits, 1);
        assert_eq!(op.terms[0].0.phase, 0);
        assert_eq!(op.terms[1].0.phase, 1);
    }

    #[test]
    fn collapse_op_real_z() {
        let mut op = CollapseOp::<ByteF64<1>>::new(1);
        op.push(ppw("Z", 0), 1.0);
        assert_eq!(op.terms.len(), 1);
        assert_eq!(op.n_qubits, 1);
        assert_eq!(op.terms[0].0.phase, 0);
    }

    #[test]
    fn collapse_op_n_qubits_stored() {
        let op = CollapseOp::<ByteF64<1>>::new(3);
        assert_eq!(op.n_qubits, 3);
    }

    #[test]
    fn rate_matrix_from_vec() {
        let r = RateMatrix::from(vec![1.0, 2.0]);
        match r {
            RateMatrix::Vector(v) => assert_eq!(v, vec![1.0, 2.0]),
            _ => panic!("expected Vector"),
        }
    }

    #[test]
    fn rate_matrix_dense_construction() {
        let r = RateMatrix::Dense(vec![vec![1.0, 0.5], vec![0.5, 1.0]]);
        match r {
            RateMatrix::Dense(m) => {
                assert_eq!(m[0], vec![1.0, 0.5]);
                assert_eq!(m[1], vec![0.5, 1.0]);
            }
            _ => panic!("expected Dense"),
        }
    }
}
