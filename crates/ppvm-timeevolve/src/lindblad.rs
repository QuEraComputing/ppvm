use std::ops::{Mul, MulAssign};

use ppvm_runtime::prelude::{ACMapAddAssign, ACMapIter, Config, PauliSum, PhasedPauliWord};

pub enum RateMatrix {
    Vector(Vec<f64>),
    Dense(Vec<Vec<f64>>),
}

impl From<Vec<f64>> for RateMatrix {
    fn from(v: Vec<f64>) -> Self {
        RateMatrix::Vector(v)
    }
}

fn get_rate(rates: &RateMatrix, i: usize, j: usize) -> f64 {
    match rates {
        RateMatrix::Vector(v) => if i == j { v[i] } else { 0.0 },
        RateMatrix::Dense(m) => m[i][j],
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

// Fields used in apply (Task 5)
#[allow(dead_code)]
pub(crate) struct LindbladTerm<T: Config> {
    pub(crate) left:   PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>,
    pub(crate) right:  PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>,
    pub(crate) a_kl:   PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>,
    pub(crate) weight: f64,
}

pub struct LindbladOp<T: Config> {
    #[allow(dead_code)] // used in apply (Task 5)
    pub(crate) terms: Vec<LindbladTerm<T>>,
}

impl<T: Config> LindbladOp<T>
where
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        Mul<Output = PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>>
        + MulAssign
        + Clone,
    T::PauliWordType: Clone,
{
    pub fn new(ops: Vec<CollapseOp<T>>, rates: RateMatrix) -> Self {
        let mut terms = Vec::new();
        let n = ops.len();

        for i in 0..n {
            for j in 0..n {
                let gamma_ij = get_rate(&rates, i, j);
                if gamma_ij == 0.0 {
                    continue;
                }
                for (sigma_k, r_ik) in &ops[i].terms {
                    for (sigma_l, r_jl) in &ops[j].terms {
                        let weight = gamma_ij * r_ik * r_jl;
                        if weight == 0.0 {
                            continue;
                        }

                        let phi_k_dag = (4 - sigma_k.phase) % 4;
                        let phi_l = sigma_l.phase;

                        let left = PhasedPauliWord::build_from_word(
                            sigma_k.word.clone(),
                            phi_k_dag,
                        );
                        let right = PhasedPauliWord::build_from_word(
                            sigma_l.word.clone(),
                            phi_l,
                        );
                        let a_kl = left.clone() * right.clone();

                        terms.push(LindbladTerm { left, right, a_kl, weight });
                    }
                }
            }
        }

        LindbladOp { terms }
    }
}

/// Returns the real part of i^phase: +1 (phase=0), -1 (phase=2), 0 otherwise.
#[inline]
#[allow(dead_code)] // used in apply and rhs (Task 6)
fn re_phase(phase: u8) -> f64 {
    match phase {
        0 => 1.0,
        2 => -1.0,
        _ => 0.0,
    }
}

impl<T: Config> LindbladOp<T>
where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>,
    T::Coeff: std::ops::AddAssign + Copy + std::ops::Mul<Output = T::Coeff>,
    T::PauliWordType: Clone,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        MulAssign + Clone,
    f64: Into<T::Coeff>,
{
    /// Accumulates `L(P)` into `result`.
    #[allow(dead_code)] // called from rhs (Task 6)
    pub(crate) fn apply(&self, p: &PauliSum<T>, result: &mut PauliSum<T>) {
        for term in &self.terms {
            for (w_a, coeff_a) in p.data().iter() {
                let wa_phased =
                    PhasedPauliWord::<T::Storage, T::BuildHasher, T::PauliWordType>::from(
                        w_a.clone(),
                    );

                // Sandwich: 2 * left * W_a * right
                let mut tmp = term.left.clone();
                tmp *= wa_phased.clone();
                tmp *= term.right.clone();
                let s = re_phase(tmp.phase);
                if s != 0.0 {
                    let c = (2.0 * term.weight * s).into() * *coeff_a;
                    *result += (tmp.word, c);
                }

                // Anticommutator: -(a_kl * W_a + W_a * a_kl)
                let mut t1 = term.a_kl.clone();
                t1 *= wa_phased.clone();
                let s1 = re_phase(t1.phase);
                if s1 != 0.0 {
                    let c = (-term.weight * s1).into() * *coeff_a;
                    *result += (t1.word, c);
                }

                let mut t2 = wa_phased;
                t2 *= term.a_kl.clone();
                let s2 = re_phase(t2.phase);
                if s2 != 0.0 {
                    let c = (-term.weight * s2).into() * *coeff_a;
                    *result += (t2.word, c);
                }
            }
        }
    }
}

#[allow(dead_code)] // called from rhs (Task 6)
/// Accumulates `i[ham, p]` into `result` using real f64 arithmetic.
///
/// For each pair of terms (W_a, h_a) in ham and (W_b, p_b) in p:
///   - Compute tmp = W_a * W_b as a PhasedPauliWord
///   - phase 1 (+i): add -2 * h_a * p_b to tmp.word
///   - phase 3 (-i): add +2 * h_a * p_b to tmp.word
///   - phase 0, 2: skip (commuting pairs cancel)
pub(crate) fn commutator_real<T: Config>(
    ham: &PauliSum<T>,
    p: &PauliSum<T>,
    result: &mut PauliSum<T>,
) where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>,
    T::Coeff: std::ops::AddAssign + Copy + std::ops::Mul<Output = T::Coeff>,
    T::PauliWordType: Clone,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        Mul<Output = PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>> + Clone,
    f64: Into<T::Coeff>,
{
    for (w_a, h_a) in ham.data().iter() {
        for (w_b, p_b) in p.data().iter() {
            let left = PhasedPauliWord::<T::Storage, T::BuildHasher, T::PauliWordType>::from(
                w_a.clone(),
            );
            let right = PhasedPauliWord::<T::Storage, T::BuildHasher, T::PauliWordType>::from(
                w_b.clone(),
            );
            let tmp = left * right;
            let coeff = match tmp.phase {
                1 => (-2.0_f64).into() * (*h_a * *p_b),
                3 => (2.0_f64).into() * (*h_a * *p_b),
                _ => continue,
            };
            *result += (tmp.word, coeff);
        }
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

    fn single_op(pauli: &str, phase: u8) -> CollapseOp<ByteF64<1>> {
        let mut op = CollapseOp::new(1);
        op.push(ppw(pauli, phase), 1.0);
        op
    }

    // ---- Task 3 tests ----

    #[test]
    fn single_real_op_z() {
        // c = Z (phase=0), rate=1.0
        // phi_k=0, phi_k†=0, phi_l=0
        // left.phase = phi_k† = 0, right.phase = 0
        // a_kl = (Z,0)*(Z,0) = (I,0)
        // weight = gamma * r_ik * r_jl = 1.0 * 1.0 * 1.0 = 1.0
        let lop = LindbladOp::new(vec![single_op("Z", 0)], RateMatrix::from(vec![1.0]));
        assert_eq!(lop.terms.len(), 1);
        let t = &lop.terms[0];
        assert_eq!(t.left.word, W1::from("Z"));
        assert_eq!(t.left.phase, 0);
        assert_eq!(t.right.word, W1::from("Z"));
        assert_eq!(t.right.phase, 0);
        assert_eq!(t.a_kl.word, W1::from("I"));
        assert_eq!(t.a_kl.phase, 0);
        assert!((t.weight - 1.0).abs() < 1e-15);
    }

    #[test]
    fn single_imaginary_op_iy() {
        // c = iY (phase=1), rate=1.0
        // phi_k=1, phi_k†=(4-1)%4=3, phi_l=1
        // left.phase = phi_k† = 3, right.phase = phi_l = 1
        // a_kl = (Y,3)*(Y,1): phases 3+1=0; bare YY=I; total (I,0)
        // weight = gamma * 1.0 * 1.0 = 1.0
        let lop = LindbladOp::new(vec![single_op("Y", 1)], RateMatrix::from(vec![1.0]));
        assert_eq!(lop.terms.len(), 1);
        let t = &lop.terms[0];
        assert_eq!(t.left.word, W1::from("Y"));
        assert_eq!(t.left.phase, 3);
        assert_eq!(t.right.word, W1::from("Y"));
        assert_eq!(t.right.phase, 1);
        assert_eq!(t.a_kl.word, W1::from("I"));
        assert_eq!(t.a_kl.phase, 0);
        assert!((t.weight - 1.0).abs() < 1e-15);
    }

    #[test]
    fn two_term_op_x_plus_iy() {
        // c = X + iY has 2 terms, so 2x2 = 4 LindbladTerms
        // none should have weight=0 (gamma=1, r_ik=1, r_jl=1 for all pairs)
        let mut op = CollapseOp::<ByteF64<1>>::new(1);
        op.push(ppw("X", 0), 1.0);
        op.push(ppw("Y", 1), 1.0);
        let lop = LindbladOp::new(vec![op], RateMatrix::from(vec![1.0]));
        assert_eq!(lop.terms.len(), 4);
        assert!(lop.terms.iter().all(|t| t.weight != 0.0));
    }

    #[test]
    fn dense_rate_matrix_off_diagonal() {
        // c1=X, c2=Y, gamma=[[1.0, 0.5],[0.5, 1.0]]
        // 4 (i,j) pairs, each with 1x1 term pair => 4 LindbladTerms
        // off-diagonal (i=0,j=1): weight = gamma_01 * r_ik * r_jl = 0.5 * 1.0 * 1.0 = 0.5
        let ops = vec![single_op("X", 0), single_op("Y", 0)];
        let rates = RateMatrix::Dense(vec![vec![1.0, 0.5], vec![0.5, 1.0]]);
        let lop = LindbladOp::new(ops, rates);
        assert_eq!(lop.terms.len(), 4);
        // The (i=0,j=1) term is index 1 (order: (0,0),(0,1),(1,0),(1,1))
        let off_diag = &lop.terms[1];
        assert!((off_diag.weight - 0.5).abs() < 1e-15);
    }

    // ---- Task 5 tests ----

    fn lindblad_x() -> LindbladOp<ByteF64<1>> {
        LindbladOp::new(vec![single_op("X", 0)], RateMatrix::from(vec![1.0]))
    }

    fn apply_to(lop: &LindbladOp<ByteF64<1>>, word: &str) -> PauliSum<ByteF64<1>> {
        let p = sum1(&[(word, 1.0)]);
        let mut result: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        lop.apply(&p, &mut result);
        result
    }

    #[test]
    fn apply_x_dephasing_lx_is_zero() {
        // c=X, L(X): 2·X·X·X - 2·X = 2X - 2X = 0
        let result = apply_to(&lindblad_x(), "X");
        assert_eq!(get_coeff(&result, "X"), 0.0);
        assert_eq!(get_coeff(&result, "Y"), 0.0);
        assert_eq!(get_coeff(&result, "Z"), 0.0);
    }

    #[test]
    fn apply_x_dephasing_lz_is_minus_4z() {
        // c=X, L(Z): XZX = -Z (XZ=-iY, -iY*X = -i*(-iZ) = -Z)
        // sandwich: 2*weight*1*(-1) = -1 to Z; anticommutator: -(X^2*Z + Z*X^2) = -(Z+Z) = -2Z
        // Total: 2*(-Z) - 2Z = -4Z
        let result = apply_to(&lindblad_x(), "Z");
        assert!((get_coeff(&result, "Z") - (-4.0)).abs() < 1e-15);
        assert_eq!(get_coeff(&result, "X"), 0.0);
    }

    #[test]
    fn apply_x_dephasing_ly_is_minus_4y() {
        // c=X, L(Y) = -4Y (XYX = -Y, analogous)
        let result = apply_to(&lindblad_x(), "Y");
        assert!((get_coeff(&result, "Y") - (-4.0)).abs() < 1e-15);
        assert_eq!(get_coeff(&result, "X"), 0.0);
    }

    #[test]
    fn apply_accumulates() {
        let lop = lindblad_x();
        let p = sum1(&[("Z", 1.0)]);
        let mut result: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        lop.apply(&p, &mut result);
        lop.apply(&p, &mut result);
        assert!((get_coeff(&result, "Z") - (-8.0)).abs() < 1e-15);
    }

    #[test]
    fn apply_lowering_op_lz() {
        // c = X + iY (un-normalised lowering operator), gamma = 1.0
        // Derivation:
        //   c† = X - iY, c†c = (X-iY)(X+iY) = 2I - 2Z
        //   Sandwich: (X-iY)·Z·(X+iY)
        //     Z·(X+iY) = ZX + i·ZY = iY + i·(-iX) = iY + X
        //     (X-iY)·(iY+X): expand term by term:
        //       X·iY = i·XY = i·iZ = -Z
        //       X·X  = I
        //       (-iY)·iY = -i²·YY = I
        //       (-iY)·X  = -i·YX = -i·(-iZ) = -Z
        //     Sum: -Z + I + I - Z = 2I - 2Z
        //   L(Z) = 2·(2I-2Z) - (2I-2Z)·Z - Z·(2I-2Z)
        //        = 4I - 4Z - (2Z - 2I) - (2Z - 2I)
        //        = 4I - 4Z - 2Z + 2I - 2Z + 2I
        //        = 8I - 8Z
        let mut op = CollapseOp::<ByteF64<1>>::new(1);
        op.push(ppw("X", 0), 1.0);
        op.push(ppw("Y", 1), 1.0);
        let lop = LindbladOp::new(vec![op], RateMatrix::from(vec![1.0]));
        let p = sum1(&[("Z", 1.0)]);
        let mut result: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        lop.apply(&p, &mut result);
        assert!((get_coeff(&result, "I") - 8.0).abs() < 1e-14);
        assert!((get_coeff(&result, "Z") - (-8.0)).abs() < 1e-14);
    }

    // ---- Task 4 tests ----

    fn sum1(terms: &[(&str, f64)]) -> PauliSum<ByteF64<1>> {
        let mut s: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        for &(w, c) in terms {
            s += (w, c);
        }
        s
    }

    fn get_coeff(s: &PauliSum<ByteF64<1>>, word: &str) -> f64 {
        use ppvm_runtime::prelude::Trace;
        let w = W1::from(word);
        s.data().trace(&w)
    }

    #[test]
    fn commutator_xx_is_zero() {
        // i[X, X] = 0: XX has phase 0 (commutes)
        let h = sum1(&[("X", 1.0)]);
        let p = sum1(&[("X", 1.0)]);
        let mut result: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        commutator_real(&h, &p, &mut result);
        assert_eq!(get_coeff(&result, "X"), 0.0);
        assert_eq!(get_coeff(&result, "Y"), 0.0);
    }

    #[test]
    fn commutator_zx_is_minus_2y() {
        // i[Z, X]: ZX = +iY (phase 1) → add -2 * 1 * 1 = -2 to Y
        let h = sum1(&[("Z", 1.0)]);
        let p = sum1(&[("X", 1.0)]);
        let mut result: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        commutator_real(&h, &p, &mut result);
        assert!((get_coeff(&result, "Y") - (-2.0)).abs() < 1e-15);
    }

    #[test]
    fn commutator_xz_is_plus_2y() {
        // i[X, Z]: XZ = -iY (phase 3) → add +2 * 1 * 1 = +2 to Y
        let h = sum1(&[("X", 1.0)]);
        let p = sum1(&[("Z", 1.0)]);
        let mut result: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        commutator_real(&h, &p, &mut result);
        assert!((get_coeff(&result, "Y") - 2.0).abs() < 1e-15);
    }

    #[test]
    fn commutator_zy_is_plus_2x() {
        // i[Z, Y]: ZY = -iX (phase 3) → add +2 * 1 * 1 = +2 to X
        let h = sum1(&[("Z", 1.0)]);
        let p = sum1(&[("Y", 1.0)]);
        let mut result: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        commutator_real(&h, &p, &mut result);
        assert!((get_coeff(&result, "X") - 2.0).abs() < 1e-15);
    }

    #[test]
    fn commutator_multiterm_linear() {
        // H = 0.5*Z, P = X + Y
        // i[0.5Z, X] = 0.5 * (-2Y) = -1.0 * Y
        // i[0.5Z, Y] = 0.5 * (+2X) = +1.0 * X
        // Total: X: +1.0, Y: -1.0
        let h = sum1(&[("Z", 0.5)]);
        let p = sum1(&[("X", 1.0), ("Y", 1.0)]);
        let mut result: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        commutator_real(&h, &p, &mut result);
        assert!((get_coeff(&result, "X") - 1.0).abs() < 1e-15);
        assert!((get_coeff(&result, "Y") - (-1.0)).abs() < 1e-15);
    }

    #[test]
    fn commutator_accumulates() {
        // Calling twice should double the result
        let h = sum1(&[("Z", 1.0)]);
        let p = sum1(&[("X", 1.0)]);
        let mut result: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        commutator_real(&h, &p, &mut result);
        commutator_real(&h, &p, &mut result);
        assert!((get_coeff(&result, "Y") - (-4.0)).abs() < 1e-15);
    }

    // ---- Task 2 tests (kept here, same module) ----

    #[test]
    fn collapse_op_x_plus_iy() {
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
