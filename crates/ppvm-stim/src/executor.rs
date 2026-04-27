use bitvec::view::BitView;
use itertools::Itertools;
use num::Integer;
use num::PrimInt;
use num::{Complex, One, ToPrimitive, Zero};
use num::complex::{Complex64, ComplexFloat};
use std::fmt::Debug;

use ppvm_runtime::prelude::*;
use ppvm_tableau::prelude::*;

use crate::tableau_program::{
    GateKind, Instruction, MeasureKind, NoiseKind, TableauProgram,
};

#[derive(Debug, thiserror::Error, PartialEq)]
#[non_exhaustive]
pub enum ExecError {}

/// Execute a normalized program against a tableau, returning the per-measurement
/// results in circuit order.
pub fn execute<T, I, C>(
    program: &TableauProgram,
    tab: &mut GeneralizedTableau<T, I, C>,
) -> Result<Vec<Option<bool>>, ExecError>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug + Send + Sync,
{
    let mut results = Vec::with_capacity(program.expected_measurement_count);
    execute_slice(&program.instructions, tab, &mut results)?;
    Ok(results)
}

/// Execute many shots, building a fresh tableau per shot via `make_tableau`.
pub fn sample<T, I, C, F>(
    program: &TableauProgram,
    num_shots: usize,
    mut make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One + Zero + Clone + num::Num + ToPrimitive + std::fmt::Debug
        + std::ops::Mul<f64> + PartialOrd<f64> + Send + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>> + From<Complex64>
        + std::ops::MulAssign + std::ops::AddAssign + One + ComplexFloat + Copy,
    I: TableauIndex + Debug + Send + Sync,
    F: FnMut() -> GeneralizedTableau<T, I, C>,
{
    (0..num_shots)
        .map(|_| {
            let mut tab = make_tableau();
            execute(program, &mut tab)
        })
        .collect()
}

fn execute_slice<T, I, C>(
    instructions: &[Instruction],
    tab: &mut GeneralizedTableau<T, I, C>,
    results: &mut Vec<Option<bool>>,
) -> Result<(), ExecError>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One + Zero + Clone + num::Num + ToPrimitive + std::fmt::Debug
        + std::ops::Mul<f64> + PartialOrd<f64> + Send + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>> + From<Complex64>
        + std::ops::MulAssign + std::ops::AddAssign + One + ComplexFloat + Copy,
    I: TableauIndex + Debug + Send + Sync,
{
    for instr in instructions {
        match instr {
            Instruction::Gate { kind, targets, .. } => match kind {
                GateKind::Reset      => targets.iter().for_each(|&q| tab.reset(q)),
                GateKind::X          => targets.iter().for_each(|&q| tab.x(q)),
                GateKind::Y          => targets.iter().for_each(|&q| tab.y(q)),
                GateKind::Z          => targets.iter().for_each(|&q| tab.z(q)),
                GateKind::H          => targets.iter().for_each(|&q| tab.h(q)),
                GateKind::S          => targets.iter().for_each(|&q| tab.s(q)),
                GateKind::SDag       => targets.iter().for_each(|&q| tab.s_adj(q)),
                GateKind::SqrtX      => targets.iter().for_each(|&q| tab.sqrt_x(q)),
                GateKind::SqrtXDag   => targets.iter().for_each(|&q| tab.sqrt_x_adj(q)),
                GateKind::SqrtY      => targets.iter().for_each(|&q| tab.sqrt_y(q)),
                GateKind::SqrtYDag   => targets.iter().for_each(|&q| tab.sqrt_y_adj(q)),
                GateKind::T          => targets.iter().for_each(|&q| tab.t(q)),
                GateKind::TDag       => targets.iter().for_each(|&q| tab.t_adj(q)),
                GateKind::RX { theta } => targets.iter().for_each(|&q| tab.rx(q, *theta)),
                GateKind::RY { theta } => targets.iter().for_each(|&q| tab.ry(q, *theta)),
                GateKind::RZ { theta } => targets.iter().for_each(|&q| tab.rz(q, *theta)),
                GateKind::U3 { theta, phi, lambda } => targets.iter().for_each(|&q| {
                    tab.u3(q, (*theta).into(), (*phi).into(), (*lambda).into())
                }),
                GateKind::CX => targets.chunks_exact(2).for_each(|p| tab.cnot(p[0], p[1])),
                GateKind::CY => targets.chunks_exact(2).for_each(|p| tab.cy(p[0], p[1])),
                GateKind::CZ => targets.chunks_exact(2).for_each(|p| tab.cz(p[0], p[1])),
            },
            Instruction::Noise { kind, targets, args, .. } => match kind {
                NoiseKind::Depolarize1 => {
                    debug_assert_eq!(args.len(), 1);
                    let p = args[0];
                    for &q in targets {
                        tab.depolarize(q, p.into());
                    }
                }
                NoiseKind::Depolarize2 => {
                    debug_assert_eq!(args.len(), 1);
                    let p = args[0];
                    for (a, b) in targets.iter().copied().tuples() {
                        tab.depolarize2(a, b, p.into());
                    }
                }
                NoiseKind::PauliChannel1 => {
                    debug_assert_eq!(args.len(), 3);
                    let ps: [T::Coeff; 3] = [args[0].into(), args[1].into(), args[2].into()];
                    for &q in targets {
                        tab.pauli_error(q, ps.clone());
                    }
                }
                NoiseKind::PauliChannel2 => {
                    debug_assert_eq!(args.len(), 15);
                    let ps: [T::Coeff; 15] = std::array::from_fn(|i| args[i].into());
                    debug_assert!(targets.len().is_even());
                    for (a, b) in targets.iter().copied().tuples() {
                        tab.two_qubit_pauli_error(a, b, ps.clone());
                    }
                }
                NoiseKind::XError => {
                    debug_assert_eq!(args.len(), 1);
                    let ps: [T::Coeff; 3] = [args[0].into(), T::Coeff::zero(), T::Coeff::zero()];
                    for &q in targets {
                        tab.pauli_error(q, ps.clone());
                    }
                }
                NoiseKind::YError => {
                    debug_assert_eq!(args.len(), 1);
                    let ps: [T::Coeff; 3] = [T::Coeff::zero(), args[0].into(), T::Coeff::zero()];
                    for &q in targets {
                        tab.pauli_error(q, ps.clone());
                    }
                }
                NoiseKind::ZError => {
                    debug_assert_eq!(args.len(), 1);
                    let ps: [T::Coeff; 3] = [T::Coeff::zero(), T::Coeff::zero(), args[0].into()];
                    for &q in targets {
                        tab.pauli_error(q, ps.clone());
                    }
                }
                NoiseKind::Loss => {
                    debug_assert_eq!(args.len(), 1);
                    for &q in targets {
                        tab.loss_channel(q, args[0].into());
                    }
                }
                NoiseKind::CorrelatedLoss => {
                    debug_assert_eq!(args.len(), 3);
                    let ps: [T::Coeff; 3] = [args[0].into(), args[1].into(), args[2].into()];
                    for (a, b) in targets.iter().copied().tuples() {
                        tab.correlated_loss_channel(a, b, ps.clone());
                    }
                }
            },
            Instruction::Measure { kind, targets, .. } => match kind {
                MeasureKind::M => {
                    for &q in targets {
                        results.push(tab.measure(q));
                    }
                }
                MeasureKind::MR => {
                    for &q in targets {
                        let outcome = tab.measure(q);
                        if outcome == Some(true) {
                            tab.x(q);
                        }
                        results.push(outcome);
                    }
                }
            },
            Instruction::Annotation { .. } => { /* phase-1 no-op */ }
            Instruction::Repeat { count, body, .. } => {
                for _ in 0..*count {
                    execute_slice(body, tab, results)?;
                }
            }
        }
    }
    Ok(())
}
