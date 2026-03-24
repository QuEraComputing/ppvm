use std::borrow::Borrow;
use std::ops::{Mul, MulAssign};

use rayon::prelude::*;

use ppvm_runtime::prelude::{
    ACMapAddAssign, ACMapBase, ACMapIter, Config, Pauli, PauliStorage, PauliSum, PauliWord,
    PauliWordTrait, PhasedPauliWord,
};

/// Direction of a ladder operator: Raise (`S₊ = X − iY`) or Lower (`S₋ = X + iY`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LadderDirection {
    Raise,
    Lower,
}

impl LadderDirection {
    /// Swaps Raise ↔ Lower (used when conjugating: `S₊† = S₋`).
    #[inline]
    pub fn flip(&self) -> Self {
        match self {
            LadderDirection::Raise => LadderDirection::Lower,
            LadderDirection::Lower => LadderDirection::Raise,
        }
    }
}

/// A single-qubit ladder operator with a qubit index and direction.
#[derive(Debug, Clone, Copy)]
pub struct LadderOp {
    pub qubit: usize,
    pub direction: LadderDirection,
}

impl LadderOp {
    /// Expands this ladder operator into a two-term `CollapseOp`.
    ///
    /// Lower (`S₋ = X + iY`): X with phase 0, Y with phase 1 (+i).
    /// Raise  (`S₊ = X − iY`): X with phase 0, Y with phase 3 (−i).
    pub fn expand<T: Config>(&self, n_qubits: usize) -> CollapseOp<T>
    where
        T::PauliWordType: PauliWordTrait,
        PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>: From<T::PauliWordType>,
    {
        let identity = T::PauliWordType::new(n_qubits);
        let x_word = identity.set_new(self.qubit, Pauli::X);
        let y_word = identity.set_new(self.qubit, Pauli::Y);
        let y_phase = match self.direction {
            LadderDirection::Lower => 1,
            LadderDirection::Raise => 3,
        };
        let mut op = CollapseOp::new(n_qubits);
        op.push(PhasedPauliWord::build_from_word(x_word, 0), 1.0);
        op.push(PhasedPauliWord::build_from_word(y_word, y_phase), 1.0);
        op
    }
}

/// User-facing input to `LindbladOp::new`: either a generic `CollapseOp` or a
/// single-qubit `LadderOp` (enabling the fast ladder kernel in Task 22+).
pub enum JumpOp<T: Config> {
    Generic(CollapseOp<T>),
    Ladder(LadderOp),
}

impl<T: Config> From<CollapseOp<T>> for JumpOp<T> {
    fn from(op: CollapseOp<T>) -> Self {
        JumpOp::Generic(op)
    }
}

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

                        terms.push(LindbladTerm { left, right, weight });
                    }
                }
            }
        }

        LindbladOp { terms }
    }
}

/// Returns the commutation parity of two Pauli words: 1 if they anticommute, 0 if they commute.
///
/// Computed as `parity((a.xbits & b.zbits) XOR (a.zbits & b.xbits))` over `n_qubits` bits.
/// Accesses `xbits` and `zbits` directly as `BitArray` fields (no trait abstraction).
/// Per-bit indexing via `a.xbits[i]` generates pure scalar `w`-register ARM code; this
/// avoids LLVM auto-vectorizing to NEON instructions, which adds scalar↔vector overhead
/// for the 1–3 byte storage sizes used in practice.
#[inline(always)]
pub(crate) fn comm_parity<A: PauliStorage, S>(
    a: &PauliWord<A, S>,
    b: &PauliWord<A, S>,
    n_qubits: usize,
) -> u8 {
    let mut parity = 0u8;
    for i in 0..n_qubits {
        let ax = a.xbits[i] as u8;
        let az = a.zbits[i] as u8;
        let bx = b.xbits[i] as u8;
        let bz = b.zbits[i] as u8;
        parity ^= (ax & bz) ^ (az & bx);
    }
    parity
}

/// Returns the real part of i^phase: +1 (phase=0), -1 (phase=2), 0 otherwise.
#[inline]
#[allow(dead_code)] // used indirectly via apply, which is called from rhs (used in solve.rs)
fn re_phase(phase: u8) -> f64 {
    match phase {
        0 => 1.0,
        2 => -1.0,
        _ => 0.0,
    }
}

/// Parallel fold/reduce over Lindblad terms for large term counts.
///
/// `#[cold]` + `#[inline(never)]` keeps this entirely out of the hot sequential
/// path. Takes `op: &LindbladOp<T>` (single pointer) to keep the argument count
/// at 3, matching the live registers already held by `rhs_into_par`.
/// Only called when `op.terms.len() >= 200` (PAR_THRESHOLD).
#[cold]
#[inline(never)]
fn apply_par<T: Config>(
    op: &LindbladOp<T>,
    p: &PauliSum<T>,
    result: &mut PauliSum<T>,
) where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>
        + Send
        + Sync,
    T::Coeff: std::ops::AddAssign + Copy + std::ops::Mul<Output = T::Coeff> + Send,
    T::PauliWordType: Clone + Borrow<PauliWord<T::Storage, T::BuildHasher>> + Send + Sync,
    T::BuildHasher: Send + Sync,
    T::Strategy: Send + Sync,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>: MulAssign + Clone,
    f64: Into<T::Coeff>,
{
    let n = p.n_qubits();
    let combined = op.terms
        .par_iter()
        .fold(
            || PauliSum::<T>::builder().n_qubits(n).build(),
            |mut local, term| {
                for (w_a, coeff_a) in p.data().iter() {
                    let wa_phased =
                        PhasedPauliWord::<T::Storage, T::BuildHasher, T::PauliWordType>::from(
                            w_a.clone(),
                        );
                    let p1 = comm_parity(term.left.word.borrow(), wa_phased.word.borrow(), n);
                    let p2 = comm_parity(wa_phased.word.borrow(), term.right.word.borrow(), n);
                    let multiplicity = p1 + p2;
                    if multiplicity > 0 {
                        let mut tmp = term.left.clone();
                        tmp *= wa_phased;
                        tmp *= term.right.clone();
                        let s = re_phase(tmp.phase);
                        if s != 0.0 {
                            let c = (multiplicity as f64 * 2.0 * term.weight * s).into()
                                * *coeff_a;
                            local += (tmp.word, c);
                        }
                    }
                }
                local
            },
        )
        .reduce(
            || PauliSum::<T>::builder().n_qubits(n).build(),
            |mut a, b| {
                for (w, c) in b.data().iter() {
                    a += (w.clone(), *c);
                }
                a
            },
        );
    for (w, c) in combined.data().iter() {
        *result += (w.clone(), *c);
    }
}

impl<T: Config> LindbladOp<T>
where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>
        + Send
        + Sync,
    T::Coeff: std::ops::AddAssign + Copy + std::ops::Mul<Output = T::Coeff> + Send,
    T::PauliWordType: Clone + Borrow<PauliWord<T::Storage, T::BuildHasher>> + Send + Sync,
    T::BuildHasher: Send + Sync,
    T::Strategy: Send + Sync,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        MulAssign + Clone,
    f64: Into<T::Coeff>,
{
    /// Accumulates `L(P)` into `result` sequentially.
    ///
    /// Uses the commutator form (Task 16): contribution is zero when `multiplicity = 0`
    /// (~25% of entries), avoiding MulAssign for those. Expected 1.5 MulAssigns/entry
    /// vs 2.5 in the old sandwich+anticommutator form.
    ///
    /// `#[inline]` allows this function to be inlined into `rhs_into` (which
    /// performs the parallel dispatch), keeping this body free of any cross-path
    /// register pressure from `apply_par`'s calling convention.
    #[inline]
    pub(crate) fn apply(&self, p: &PauliSum<T>, result: &mut PauliSum<T>) {
        let n = p.n_qubits();
        for term in &self.terms {
            for (w_a, coeff_a) in p.data().iter() {
                let wa_phased =
                    PhasedPauliWord::<T::Storage, T::BuildHasher, T::PauliWordType>::from(
                        w_a.clone(),
                    );
                let p1 = comm_parity(term.left.word.borrow(), wa_phased.word.borrow(), n);
                let p2 = comm_parity(wa_phased.word.borrow(), term.right.word.borrow(), n);
                let multiplicity = p1 + p2;
                if multiplicity > 0 {
                    let mut tmp = term.left.clone();
                    tmp *= wa_phased;
                    tmp *= term.right.clone();
                    let s = re_phase(tmp.phase);
                    if s != 0.0 {
                        let c = (multiplicity as f64 * 2.0 * term.weight * s).into() * *coeff_a;
                        *result += (tmp.word, c);
                    }
                }
            }
        }
    }
}

/// Computes `dP/dt = i[ham, P] + L(P)` and returns the result.
///
/// Creates a fresh zero-initialised `PauliSum` using `T::Strategy::default()`, calls
/// `commutator_real` if `ham` is provided, then `lindblad.apply`, then `truncate()`.
pub fn rhs<T: Config>(
    ham: Option<&PauliSum<T>>,
    lindblad: &LindbladOp<T>,
    p: &PauliSum<T>,
) -> PauliSum<T>
where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>
        + ACMapBase
        + Send
        + Sync,
    T::Coeff: std::ops::AddAssign + Copy + std::ops::Mul<Output = T::Coeff> + Send,
    T::PauliWordType: Clone + Borrow<PauliWord<T::Storage, T::BuildHasher>> + Send + Sync,
    T::BuildHasher: Send + Sync,
    T::Strategy: Send + Sync,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        Mul<Output = PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>>
        + MulAssign
        + Clone,
    f64: Into<T::Coeff>,
{
    let mut result = PauliSum::<T>::builder().n_qubits(p.n_qubits()).build();
    rhs_into(ham, lindblad, p, &mut result);
    result
}

/// Parallel implementation of [`rhs_into`] for large Lindblad operators.
///
/// `#[cold]` + `#[inline(never)]` fully isolates all rayon code from the sequential
/// hot path in `rhs_into`. Without this isolation, rayon's atomic operations
/// (acquire/release fences) would be visible to LLVM's optimizer in the inlined
/// `rhs_into` body, preventing memory-access reordering in the sequential loop
/// and causing a ~60 µs regression even when the parallel path is never taken.
#[cold]
#[inline(never)]
fn rhs_into_par<T: Config>(
    ham: Option<&PauliSum<T>>,
    lindblad: &LindbladOp<T>,
    p: &PauliSum<T>,
    result: &mut PauliSum<T>,
) where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>
        + ACMapBase
        + Send
        + Sync,
    T::Coeff: std::ops::AddAssign + Copy + std::ops::Mul<Output = T::Coeff> + Send,
    T::PauliWordType: Clone + Borrow<PauliWord<T::Storage, T::BuildHasher>> + Send + Sync,
    T::BuildHasher: Send + Sync,
    T::Strategy: Send + Sync,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        Mul<Output = PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>>
        + MulAssign
        + Clone,
    f64: Into<T::Coeff>,
{
    result.data_mut().clear();
    if let Some(h) = ham {
        commutator_real(h, p, result);
    }
    apply_par(lindblad, p, result);
    result.truncate();
}

/// In-place version of [`rhs`].
///
/// Clears `result`, computes `dP/dt = i[ham, P] + L(P)` into it, then calls `truncate()`.
/// Retains the allocated capacity of `result`.
pub(crate) fn rhs_into<T: Config>(
    ham: Option<&PauliSum<T>>,
    lindblad: &LindbladOp<T>,
    p: &PauliSum<T>,
    result: &mut PauliSum<T>,
) where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>
        + ACMapBase
        + Send
        + Sync,
    T::Coeff: std::ops::AddAssign + Copy + std::ops::Mul<Output = T::Coeff> + Send,
    T::PauliWordType: Clone + Borrow<PauliWord<T::Storage, T::BuildHasher>> + Send + Sync,
    T::BuildHasher: Send + Sync,
    T::Strategy: Send + Sync,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        Mul<Output = PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>>
        + MulAssign
        + Clone,
    f64: Into<T::Coeff>,
{
    // Empirical crossover: ~200 terms on a 14-core machine. Check before any work
    // so the parallel path is fully isolated in `rhs_into_par` (cold, never-inline).
    // This ensures the sequential hot path below has ZERO rayon code in scope —
    // rayon's atomics would otherwise prevent LLVM from reordering memory accesses
    // in the sequential loop, causing a ~60 µs regression even on the cold branch.
    if lindblad.terms.len() >= 200 {
        return rhs_into_par(ham, lindblad, p, result);
    }
    result.data_mut().clear();
    if let Some(h) = ham {
        commutator_real(h, p, result);
    }
    lindblad.apply(p, result);
    result.truncate();
}

/// Accumulates `i[ham, p]` into `result` using real f64 arithmetic.
#[allow(dead_code)] // called from rhs, which is called from solve.rs (Task 7)
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
        // Hoist left: depends only on w_a, not on w_b.
        let left = PhasedPauliWord::<T::Storage, T::BuildHasher, T::PauliWordType>::from(
            w_a.clone(),
        );
        for (w_b, p_b) in p.data().iter() {
            let right = PhasedPauliWord::<T::Storage, T::BuildHasher, T::PauliWordType>::from(
                w_b.clone(),
            );
            let tmp = left.clone() * right;
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
    use ppvm_runtime::strategy::CoefficientThreshold;

    type W1 = PauliWord<[u8; 1], fxhash::FxBuildHasher>;
    type W3 = PauliWord<[u8; 3], fxhash::FxBuildHasher>;
    type PPW1 = PhasedPauliWord<[u8; 1], fxhash::FxBuildHasher, W1>;

    fn ppw(pauli: &str, phase: u8) -> PPW1 {
        PhasedPauliWord::build_from_word(W1::from(pauli), phase)
    }

    fn single_op(pauli: &str, phase: u8) -> CollapseOp<ByteF64<1>> {
        let mut op = CollapseOp::new(1);
        op.push(ppw(pauli, phase), 1.0);
        op
    }

    // ---- Task 13 tests ----

    #[test]
    fn comm_parity_single_qubit_pairs() {
        // Commuting pairs: parity = 0
        assert_eq!(comm_parity(&W1::from("I"), &W1::from("X"), 1), 0); // IX
        assert_eq!(comm_parity(&W1::from("X"), &W1::from("I"), 1), 0); // XI
        assert_eq!(comm_parity(&W1::from("X"), &W1::from("X"), 1), 0); // XX
        assert_eq!(comm_parity(&W1::from("Y"), &W1::from("Y"), 1), 0); // YY
        assert_eq!(comm_parity(&W1::from("Z"), &W1::from("Z"), 1), 0); // ZZ

        // Anticommuting pairs: parity = 1
        assert_eq!(comm_parity(&W1::from("X"), &W1::from("Y"), 1), 1); // XY
        assert_eq!(comm_parity(&W1::from("X"), &W1::from("Z"), 1), 1); // XZ
        assert_eq!(comm_parity(&W1::from("Y"), &W1::from("Z"), 1), 1); // YZ
        assert_eq!(comm_parity(&W1::from("Y"), &W1::from("X"), 1), 1); // YX
        assert_eq!(comm_parity(&W1::from("Z"), &W1::from("X"), 1), 1); // ZX
        assert_eq!(comm_parity(&W1::from("Z"), &W1::from("Y"), 1), 1); // ZY
    }

    #[test]
    fn comm_parity_multi_qubit() {
        type W2 = PauliWord<[u8; 1], fxhash::FxBuildHasher>;

        // XZ vs ZX: qubit 0 (X,Z)->1, qubit 1 (Z,X)->1; parity = 0 (even number of anticommuting)
        assert_eq!(comm_parity(&W2::from("XZ"), &W2::from("ZX"), 2), 0);

        // XY vs IZ: qubit 0 (X,I)->0, qubit 1 (Y,Z)->1; parity = 1
        assert_eq!(comm_parity(&W2::from("XY"), &W2::from("IZ"), 2), 1);

        // XZ vs XI: qubit 0 (X,X)->0, qubit 1 (Z,I)->0; parity = 0
        assert_eq!(comm_parity(&W2::from("XZ"), &W2::from("XI"), 2), 0);

        // XX vs YI: qubit 0 (X,Y)->1, qubit 1 (X,I)->0; parity = 1
        assert_eq!(comm_parity(&W2::from("XX"), &W2::from("YI"), 2), 1);
    }

    #[test]
    fn imaginary_phase_term_gives_zero_x() {
        // Regression: when left*right has an imaginary phase, re_phase filters it.
        //
        // Use c1=Y (phase=0), c2=Z (phase=0) with γ₁₂=1 (off-diagonal rate matrix).
        // Cross-pair: left={Y,0}, right={Z,0}.
        //   Y*Z = iX (phase=1), so re_phase=0 → no X term from any W_a.
        // Apply to P=Z: verify no X coefficient appears.
        let ops = vec![single_op("Y", 0), single_op("Z", 0)];
        let rates = RateMatrix::Dense(vec![vec![0.0, 1.0], vec![1.0, 0.0]]);
        let lop = LindbladOp::new(ops, rates);

        let p = sum1(&[("Z", 1.0)]);
        let mut result: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        lop.apply(&p, &mut result);

        assert_eq!(get_coeff(&result, "X"), 0.0,
            "imaginary-phase product left*W_a*right must produce zero real coefficient");
    }

    // ---- Task 17 tests ----

    #[test]
    fn comm_parity_n20() {
        // Per-bit parity for n=20 qubits using [u8; 3] storage.

        // Single anticommuting pair at qubit 0: parity = 1
        let x0 = W3::from("XIIIIIIIIIIIIIIIIIII");
        let z0 = W3::from("ZIIIIIIIIIIIIIIIIIII");
        assert_eq!(comm_parity(&x0, &z0, 20), 1);

        // Two anticommuting pairs at qubits 0 and 1: parity = 0 (even count)
        let xz = W3::from("XZIIIIIIIIIIIIIIIIII");
        let zx = W3::from("ZXIIIIIIIIIIIIIIIIII");
        assert_eq!(comm_parity(&xz, &zx, 20), 0);

        // 20 anticommuting pairs: parity = 20 & 1 = 0
        let all_x = W3::from("XXXXXXXXXXXXXXXXXXXX");
        let all_z = W3::from("ZZZZZZZZZZZZZZZZZZZZ");
        assert_eq!(comm_parity(&all_x, &all_z, 20), 0);
    }

    #[test]
    fn comm_parity_zero_padding() {
        // Verify unused high bits (bits 5–7 in [u8; 1] for n=5) don't produce
        // spurious parity. The per-bit loop stops at n_qubits=5, so padding bits
        // are never accessed — this test confirms the loop bound is correct.

        // 5 anticommuting pairs (XXXXX vs ZZZZZ): 5 XOR 1s → parity = 5 & 1 = 1.
        assert_eq!(comm_parity(&W1::from("XXXXX"), &W1::from("ZZZZZ"), 5), 1);

        // 4 anticommuting pairs: parity = 4 & 1 = 0
        assert_eq!(comm_parity(&W1::from("XXXX"), &W1::from("ZZZZ"), 4), 0);
    }

    // ---- Task 16 tests ----

    #[test]
    fn commutator_form_zero() {
        // multiplicity=0 path: both left and right commute with W_a → no contribution.
        //
        // c=Z (phase=0), left={Z,0}, right={Z,0}.
        // W_a = Z: comm_parity(Z,Z)=0 and comm_parity(Z,Z)=0 → multiplicity=0.
        // Apply to P=Z only; expect zero output (L_term(Z) = 0 before summing other
        // contributions; with this single term and W_a=Z, the inner body is skipped).
        let lop = LindbladOp::new(vec![single_op("Z", 0)], RateMatrix::from(vec![1.0]));
        // Verify multiplicity=0 for (left=Z, W_a=Z, right=Z)
        let t = &lop.terms[0];
        let wa = W1::from("Z");
        assert_eq!(comm_parity(&t.left.word, &wa, 1), 0);
        assert_eq!(comm_parity(&wa, &t.right.word, 1), 0);
        // Behavioral check: L(Z) = 2·Z·Z·Z − {Z·Z,Z} = 2Z − 2Z = 0
        let result = apply_to(&lop, "Z");
        assert_eq!(get_coeff(&result, "Z"), 0.0);
        assert_eq!(get_coeff(&result, "X"), 0.0);
        assert_eq!(get_coeff(&result, "Y"), 0.0);
    }

    #[test]
    fn commutator_form_double() {
        // multiplicity=2 path: both parities=1 → coefficient = 4·weight·re_phase·coeff_a.
        //
        // c=X (phase=0), left={X,0}, right={X,0}, weight=1.
        // W_a=Z: comm_parity(X,Z)=1, comm_parity(Z,X)=1 → multiplicity=2.
        // left*W_a*right = X·Z·X: XZ = -iY (phase 3), then (-iY)·X = -i·YX = -i·(-iZ) = -Z.
        // Phase of X·Z·X: 2 (re_phase = -1).
        // Coefficient = 2 * 1 * 1 * (-1) = -4 → Z entry gets -4·coeff_a = -4.
        //
        // Manual derivation: L_X(Z) = 2·X·Z·X − {X·X,Z} = 2·(-Z) − {I,Z} = -2Z − 2Z = -4Z.
        let result = apply_to(&lindblad_x(), "Z");
        // Verify multiplicity=2 for this term
        let t = &lindblad_x().terms[0];
        let wz = W1::from("Z");
        assert_eq!(comm_parity(&t.left.word, &wz, 1) + comm_parity(&wz, &t.right.word, 1), 2);
        assert!((get_coeff(&result, "Z") - (-4.0)).abs() < 1e-15,
            "multiplicity=2 must give coefficient 4*weight*re_phase = -4");
    }

    // ---- Task 3 tests ----

    #[test]
    fn single_real_op_z() {
        // c = Z (phase=0), rate=1.0
        // phi_k=0, phi_k†=0, phi_l=0
        // left.phase = phi_k† = 0, right.phase = 0
        // weight = gamma * r_ik * r_jl = 1.0 * 1.0 * 1.0 = 1.0
        let lop = LindbladOp::new(vec![single_op("Z", 0)], RateMatrix::from(vec![1.0]));
        assert_eq!(lop.terms.len(), 1);
        let t = &lop.terms[0];
        assert_eq!(t.left.word, W1::from("Z"));
        assert_eq!(t.left.phase, 0);
        assert_eq!(t.right.word, W1::from("Z"));
        assert_eq!(t.right.phase, 0);
        assert!((t.weight - 1.0).abs() < 1e-15);
    }

    #[test]
    fn single_imaginary_op_iy() {
        // c = iY (phase=1), rate=1.0
        // phi_k=1, phi_k†=(4-1)%4=3, phi_l=1
        // left.phase = phi_k† = 3, right.phase = phi_l = 1
        // weight = gamma * 1.0 * 1.0 = 1.0
        let lop = LindbladOp::new(vec![single_op("Y", 1)], RateMatrix::from(vec![1.0]));
        assert_eq!(lop.terms.len(), 1);
        let t = &lop.terms[0];
        assert_eq!(t.left.word, W1::from("Y"));
        assert_eq!(t.left.phase, 3);
        assert_eq!(t.right.word, W1::from("Y"));
        assert_eq!(t.right.phase, 1);
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

    // ---- Task 21 tests ----

    #[test]
    fn ladder_direction_flip() {
        assert_eq!(LadderDirection::Lower.flip(), LadderDirection::Raise);
        assert_eq!(LadderDirection::Raise.flip(), LadderDirection::Lower);
    }

    #[test]
    fn ladder_op_expand_lower() {
        let op = LadderOp { qubit: 0, direction: LadderDirection::Lower };
        let expanded = op.expand::<ByteF64<1>>(1);
        assert_eq!(expanded.terms.len(), 2);
        // First term: X with phase 0
        assert_eq!(expanded.terms[0].0.word, W1::from("X"));
        assert_eq!(expanded.terms[0].0.phase, 0);
        assert!((expanded.terms[0].1 - 1.0).abs() < 1e-15);
        // Second term: Y with phase 1 (+i)
        assert_eq!(expanded.terms[1].0.word, W1::from("Y"));
        assert_eq!(expanded.terms[1].0.phase, 1);
        assert!((expanded.terms[1].1 - 1.0).abs() < 1e-15);
    }

    #[test]
    fn ladder_op_expand_raise() {
        let op = LadderOp { qubit: 0, direction: LadderDirection::Raise };
        let expanded = op.expand::<ByteF64<1>>(1);
        assert_eq!(expanded.terms.len(), 2);
        // First term: X with phase 0
        assert_eq!(expanded.terms[0].0.word, W1::from("X"));
        assert_eq!(expanded.terms[0].0.phase, 0);
        // Second term: Y with phase 3 (-i)
        assert_eq!(expanded.terms[1].0.word, W1::from("Y"));
        assert_eq!(expanded.terms[1].0.phase, 3);
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

    // ---- Task 6 tests ----

    fn empty_lindblad() -> LindbladOp<ByteF64<1>> {
        LindbladOp::new(vec![], RateMatrix::from(vec![]))
    }

    #[test]
    fn rhs_pure_hamiltonian() {
        // H = 0.5*Z, no Lindblad, P = X
        // i[0.5Z, X] = 0.5 * (-2Y) = -1.0 * Y
        let h = sum1(&[("Z", 0.5)]);
        let p = sum1(&[("X", 1.0)]);
        let result = rhs(Some(&h), &empty_lindblad(), &p);
        assert!((get_coeff(&result, "Y") - (-1.0)).abs() < 1e-15);
        assert_eq!(get_coeff(&result, "X"), 0.0);
        assert_eq!(get_coeff(&result, "Z"), 0.0);
    }

    #[test]
    fn rhs_pure_lindblad() {
        // c=X, gamma=1, P=Z => L(Z) = -4Z
        let p = sum1(&[("Z", 1.0)]);
        let result = rhs(None, &lindblad_x(), &p);
        assert!((get_coeff(&result, "Z") - (-4.0)).abs() < 1e-15);
        assert_eq!(get_coeff(&result, "X"), 0.0);
    }

    #[test]
    fn rhs_ham_and_lindblad() {
        // H = 0.5*Z, c=X, gamma=1, P=X
        // Hamiltonian: i[0.5Z, X] = -Y; Lindblad: L(X) = 0
        // Result: {Y: -1.0}
        let h = sum1(&[("Z", 0.5)]);
        let p = sum1(&[("X", 1.0)]);
        let result = rhs(Some(&h), &lindblad_x(), &p);
        assert!((get_coeff(&result, "Y") - (-1.0)).abs() < 1e-15);
        assert_eq!(get_coeff(&result, "X"), 0.0);
    }

    #[test]
    fn rhs_no_ham_no_lindblad() {
        // No ham, empty Lindblad, P = X => dP/dt = 0
        let p = sum1(&[("X", 1.0)]);
        let result = rhs(None, &empty_lindblad(), &p);
        assert_eq!(get_coeff(&result, "X"), 0.0);
    }

    #[test]
    fn rhs_truncates_small_terms() {
        // Use CoefficientThreshold strategy with a large threshold (1.0).
        // H = 0.5*Z, P = X gives result Y: -1.0. Since |-1.0| >= 1.0, Y survives.
        // But if we scale P by 1e-13, the Y term (-1e-13) is below the threshold.
        type ThreshConfig = ByteF64<1, CoefficientThreshold>;

        let mut h: PauliSum<ThreshConfig> = PauliSum::builder()
            .n_qubits(1)
            .strategy(CoefficientThreshold(1.0))
            .build();
        h += ("Z", 0.5_f64);
        let mut p: PauliSum<ThreshConfig> = PauliSum::builder()
            .n_qubits(1)
            .strategy(CoefficientThreshold(1.0))
            .build();
        // P = 1e-13 * X: result would be Y: -1e-13, which is below threshold 1.0
        p += ("X", 1e-13_f64);
        let lop = LindbladOp::<ThreshConfig>::new(vec![], RateMatrix::from(vec![]));
        let result = rhs(Some(&h), &lop, &p);
        // Y: -1e-13 * 0.5 * (-2) = -1e-13 should be truncated by threshold 1.0
        use ppvm_runtime::prelude::Trace;
        let w = PauliWord::<[u8; 1], fxhash::FxBuildHasher>::from("Y");
        let y_coeff = result.data().trace(&w);
        assert_eq!(y_coeff, 0.0, "small term should be truncated");
    }

    // ---- Task 18 tests ----

    /// Build the benchmark Lindblad (n=5, dense): 5 lowering ops c_i = X_i + iY_i,
    /// dense 5×5 rate matrix γ_ij = 1/(1+|i−j|). Produces 100 LindbladTerms.
    fn build_benchmark_lindblad() -> LindbladOp<ByteF64<1>> {
        let n = 5usize;
        let ppw5 = |s: &str, ph: u8| -> PhasedPauliWord<[u8; 1], fxhash::FxBuildHasher, W1> {
            PhasedPauliWord::build_from_word(W1::from(s), ph)
        };
        let template = vec!['I'; n];
        let mut ops: Vec<CollapseOp<ByteF64<1>>> = Vec::new();
        for i in 0..n {
            let mut op = CollapseOp::new(n);
            let mut px = template.clone();
            let mut py = template.clone();
            px[i] = 'X';
            py[i] = 'Y';
            let sx: String = px.iter().collect();
            let sy: String = py.iter().collect();
            op.push(ppw5(&sx, 0), 1.0);
            op.push(ppw5(&sy, 1), 1.0);
            ops.push(op);
        }
        let rates: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| 1.0 / (1.0 + (i as f64 - j as f64).abs())).collect())
            .collect();
        LindbladOp::new(ops, RateMatrix::Dense(rates))
    }

    /// Helper: assert determinism and linearity of `lop.apply` for a given state.
    fn check_apply_consistency(lop: &LindbladOp<ByteF64<1>>, p_str: &str, n: usize) {
        use ppvm_runtime::prelude::Trace;

        let mut p: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(n).build();
        p += (p_str, 1.0_f64);

        let mut result1: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(n).build();
        lop.apply(&p, &mut result1);

        // Determinism: second call must agree on every coefficient.
        let mut result2: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(n).build();
        lop.apply(&p, &mut result2);
        for (w, c1) in result1.data().iter() {
            let c2 = result2.data().trace(w);
            assert!(
                (c1 - c2).abs() < 1e-14,
                "non-deterministic apply (n={n}): word {:?}, run1={c1}, run2={c2}", w
            );
        }

        // Linearity: apply(2·P) must equal 2·apply(P).
        let mut p2: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(n).build();
        p2 += (p_str, 2.0_f64);
        let mut result_2p: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(n).build();
        lop.apply(&p2, &mut result_2p);
        for (w, c1) in result1.data().iter() {
            let c_2p = result_2p.data().trace(w);
            assert!(
                (c_2p - 2.0 * c1).abs() < 1e-13,
                "linearity violated (n={n}): word {:?}, 2·apply(P)={}, apply(2P)={}",
                w, 2.0 * c1, c_2p
            );
        }
    }

    #[test]
    fn parallel_matches_sequential() {
        // Sequential path: n=5 benchmark Lindblad, 100 terms < PAR_THRESHOLD=200.
        // Tests that the direct-accumulation path is correct.
        let lop5 = build_benchmark_lindblad();
        assert!(lop5.terms.len() < 200, "expected sequential path");
        check_apply_consistency(&lop5, "ZIIII", 5);

        // Parallel path: n=8 dense lowering Lindblad, 256 terms > PAR_THRESHOLD=200.
        // (8 ops × 2 terms × 8 ops × 2 terms; n=8 fits in [u8;1].)
        // Tests that fold/reduce accumulates identically to sequential.
        let n = 8usize;
        let ppw8 = |s: &str, ph: u8| -> PhasedPauliWord<[u8; 1], fxhash::FxBuildHasher, W1> {
            PhasedPauliWord::build_from_word(W1::from(s), ph)
        };
        let template = vec!['I'; n];
        let mut ops8: Vec<CollapseOp<ByteF64<1>>> = Vec::new();
        for i in 0..n {
            let mut op = CollapseOp::new(n);
            let mut px = template.clone();
            let mut py = template.clone();
            px[i] = 'X';
            py[i] = 'Y';
            let sx: String = px.iter().collect();
            let sy: String = py.iter().collect();
            op.push(ppw8(&sx, 0), 1.0);
            op.push(ppw8(&sy, 1), 1.0);
            ops8.push(op);
        }
        let rates8: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| 1.0 / (1.0 + (i as f64 - j as f64).abs())).collect())
            .collect();
        let lop8 = LindbladOp::new(ops8, RateMatrix::Dense(rates8));
        assert!(lop8.terms.len() > 200, "expected parallel path");
        check_apply_consistency(&lop8, "ZIIIIIII", n);
    }
}
