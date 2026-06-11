// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::BitView;
use itertools::Itertools;
use num::Integer;
use num::PrimInt;
use num::complex::{Complex64, ComplexFloat};
use num::{Complex, One, ToPrimitive, Zero};
use std::fmt::Debug;

use ppvm_runtime::prelude::*;
use ppvm_tableau::prelude::*;
use stim_parser::ast::{GateName, MeasureName, NoiseName};
use stim_parser::extended::{Axis, ExtendedInstruction, ExtendedProgram, RawPassthrough};

use crate::prepare::{ExecError, prepare};

/// Validate and execute a parsed extended Stim program against a tableau,
/// returning the per-measurement results in circuit order.
pub fn execute<T, I, C>(
    program: &ExtendedProgram,
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
    prepare(program)?;
    let count = program.measurement_count();
    let mut results = Vec::with_capacity(count);
    execute_prepared(&program.instructions, tab, &mut results);
    Ok(results)
}

/// Execute many shots serially, building a fresh tableau for shot `i` via
/// `make_tableau(i)`.
///
/// The shot index lets callers derive a deterministic per-shot seed (e.g.
/// `seed + i`) so results are independent of evaluation order — the same
/// factory then yields identical results from `sample_parallel` (when the `rayon` feature is enabled).
pub fn sample_serial<T, I, C, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
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
    F: Fn(usize) -> GeneralizedTableau<T, I, C>,
{
    prepare(program)?;
    let count = program.measurement_count();
    Ok((0..num_shots)
        .map(|i| {
            let mut tab = make_tableau(i);
            let mut results = Vec::with_capacity(count);
            execute_prepared(&program.instructions, &mut tab, &mut results);
            results
        })
        .collect())
}

/// Execute many shots, building a fresh tableau for shot `i` via
/// `make_tableau(i)`.
///
/// When the `rayon` feature is enabled this dispatches to `sample_parallel`
/// for batches large enough to amortise thread-scheduling overhead, and to
/// [`sample_serial`] otherwise. Without the feature it is always serial. Use
/// [`sample_serial`] / `sample_parallel` directly to force one or the other.
pub fn sample<T, I, C, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
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
    F: Fn(usize) -> GeneralizedTableau<T, I, C> + Sync,
{
    #[cfg(feature = "rayon")]
    {
        // Below ~4 shots per thread, rayon's scheduling overhead outweighs
        // the gain, so stay serial; with a single thread there is no upside.
        let n_threads = rayon::current_num_threads();
        if n_threads <= 1 || num_shots < 4 * n_threads {
            sample_serial(program, num_shots, make_tableau)
        } else {
            sample_parallel(program, num_shots, make_tableau)
        }
    }
    #[cfg(not(feature = "rayon"))]
    {
        sample_serial(program, num_shots, make_tableau)
    }
}

/// Execute many shots in parallel across the global rayon thread pool,
/// building a fresh tableau for shot `i` via `make_tableau(i)`.
///
/// Per-shot seeds derived from `i` make the result independent of how rayon
/// schedules the work, so a seeded factory yields the same shots (in the same
/// order) as [`sample_serial`]. Thread count follows rayon's global pool
/// (set `RAYON_NUM_THREADS` to override).
#[cfg(feature = "rayon")]
pub fn sample_parallel<T, I, C, F>(
    program: &ExtendedProgram,
    num_shots: usize,
    make_tableau: F,
) -> Result<Vec<Vec<Option<bool>>>, ExecError>
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
    F: Fn(usize) -> GeneralizedTableau<T, I, C> + Sync,
{
    use rayon::prelude::*;
    prepare(program)?;
    let count = program.measurement_count();
    Ok((0..num_shots)
        .into_par_iter()
        .map(|i| {
            let mut tab = make_tableau(i);
            let mut results = Vec::with_capacity(count);
            execute_prepared(&program.instructions, &mut tab, &mut results);
            results
        })
        .collect())
}

/// Dispatch a slice of validated instructions onto a tableau, appending
/// measurement bits to `results`. Skip-validates: the caller is responsible
/// for having run [`prepare`] on the originating program. Used by
/// [`execute`] / [`sample`] internally and by the Python `tab.run()` path
/// where `StimProgram` already cached the prepare step.
pub fn execute_prepared<T, I, C>(
    instructions: &[ExtendedInstruction],
    tab: &mut GeneralizedTableau<T, I, C>,
    results: &mut Vec<Option<bool>>,
) where
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
    for instr in instructions {
        match instr {
            ExtendedInstruction::Raw(RawPassthrough::Gate { name, targets, .. }) => match name {
                GateName::Reset | GateName::ResetZ => targets.iter().for_each(|&q| tab.reset(q)),
                GateName::X => targets.iter().for_each(|&q| tab.x(q)),
                GateName::Y => targets.iter().for_each(|&q| tab.y(q)),
                GateName::Z => targets.iter().for_each(|&q| tab.z(q)),
                GateName::H | GateName::HXZ => targets.iter().for_each(|&q| tab.h(q)),
                GateName::S | GateName::SqrtZ => targets.iter().for_each(|&q| tab.s(q)),
                GateName::SDag | GateName::SqrtZDag => {
                    targets.iter().for_each(|&q| tab.s_adj(q));
                }
                GateName::SqrtX => targets.iter().for_each(|&q| tab.sqrt_x(q)),
                GateName::SqrtXDag => targets.iter().for_each(|&q| tab.sqrt_x_adj(q)),
                GateName::SqrtY => targets.iter().for_each(|&q| tab.sqrt_y(q)),
                GateName::SqrtYDag => targets.iter().for_each(|&q| tab.sqrt_y_adj(q)),
                GateName::Identity => {}
                GateName::CX | GateName::ZCX | GateName::CNot => {
                    targets.chunks_exact(2).for_each(|p| tab.cnot(p[0], p[1]));
                }
                GateName::CY | GateName::ZCY => {
                    targets.chunks_exact(2).for_each(|p| tab.cy(p[0], p[1]));
                }
                GateName::CZ | GateName::ZCZ => {
                    targets.chunks_exact(2).for_each(|p| tab.cz(p[0], p[1]));
                }
                GateName::Swap
                | GateName::ISwap
                | GateName::ISwapDag
                | GateName::SqrtXX
                | GateName::SqrtYY
                | GateName::SqrtZZ
                | GateName::CXSwap
                | GateName::SwapCX
                | GateName::XCX
                | GateName::XCY
                | GateName::XCZ
                | GateName::YCX
                | GateName::YCY
                | GateName::YCZ
                | GateName::CXYZ
                | GateName::CZYX
                | GateName::HXY
                | GateName::HYZ => {
                    unreachable!("unsupported gate {name:?} should have been rejected by prepare")
                }
            },
            ExtendedInstruction::T { targets, .. } => targets.iter().for_each(|&q| tab.t(q)),
            ExtendedInstruction::TDag { targets, .. } => {
                targets.iter().for_each(|&q| tab.t_adj(q));
            }
            ExtendedInstruction::Rotation {
                axis,
                theta,
                targets,
                ..
            } => match axis {
                Axis::X => targets.iter().for_each(|&q| tab.rx(q, *theta)),
                Axis::Y => targets.iter().for_each(|&q| tab.ry(q, *theta)),
                Axis::Z => targets.iter().for_each(|&q| tab.rz(q, *theta)),
            },
            ExtendedInstruction::U3 {
                theta,
                phi,
                lambda,
                targets,
                ..
            } => targets
                .iter()
                .for_each(|&q| tab.u3(q, (*theta).into(), (*phi).into(), (*lambda).into())),
            ExtendedInstruction::Raw(RawPassthrough::Noise {
                name,
                targets,
                args,
                ..
            }) => match name {
                NoiseName::Depolarize1 => {
                    debug_assert_eq!(args.len(), 1);
                    let p = args[0];
                    for &q in targets {
                        tab.depolarize(q, p.into());
                    }
                }
                NoiseName::Depolarize2 => {
                    debug_assert_eq!(args.len(), 1);
                    let p = args[0];
                    for (a, b) in targets.iter().copied().tuples() {
                        tab.depolarize2(a, b, p.into());
                    }
                }
                NoiseName::PauliChannel1 => {
                    debug_assert_eq!(args.len(), 3);
                    let ps: [T::Coeff; 3] = [args[0].into(), args[1].into(), args[2].into()];
                    for &q in targets {
                        tab.pauli_error(q, ps.clone());
                    }
                }
                NoiseName::PauliChannel2 => {
                    debug_assert_eq!(args.len(), 15);
                    let ps: [T::Coeff; 15] = std::array::from_fn(|i| args[i].into());
                    debug_assert!(targets.len().is_even());
                    for (a, b) in targets.iter().copied().tuples() {
                        tab.two_qubit_pauli_error(a, b, ps.clone());
                    }
                }
                NoiseName::XError | NoiseName::YError | NoiseName::ZError => {
                    debug_assert_eq!(args.len(), 1);
                    let p: T::Coeff = args[0].into();
                    let zero = T::Coeff::zero();
                    let ps: [T::Coeff; 3] = match name {
                        NoiseName::XError => [p, zero.clone(), zero],
                        NoiseName::YError => [zero.clone(), p, zero],
                        NoiseName::ZError => [zero.clone(), zero, p],
                        _ => unreachable!(),
                    };
                    for &q in targets {
                        tab.pauli_error(q, ps.clone());
                    }
                }
                NoiseName::IError
                | NoiseName::HeraldedErase
                | NoiseName::HeraldedPauliChannel1
                | NoiseName::CorrelatedError
                | NoiseName::ElseCorrelatedError => {
                    unreachable!("unsupported noise {name:?} should have been rejected by prepare")
                }
            },
            ExtendedInstruction::Loss { p, targets, .. } => {
                for &q in targets {
                    tab.loss_channel(q, (*p).into());
                }
            }
            ExtendedInstruction::CorrelatedLoss { ps, targets, .. } => {
                let ps: [T::Coeff; 3] = [ps[0].into(), ps[1].into(), ps[2].into()];
                for &(a, b) in targets {
                    tab.correlated_loss_channel(a, b, ps.clone());
                }
            }
            ExtendedInstruction::Raw(RawPassthrough::Measure {
                name,
                args,
                targets,
                ..
            }) => {
                let noise = args.first().copied().unwrap_or(0.0);
                match name {
                    MeasureName::M | MeasureName::MZ => {
                        for &q in targets {
                            results.push(tab.measure_noisy(q, noise));
                        }
                    }
                    // MR cannot delegate to `measure_noisy` because the reset must use
                    // the *true* outcome — but it shares the readout flip via
                    // `flip_with_prob`, so the RNG-draw shape matches MZ exactly.
                    MeasureName::MR => {
                        for &q in targets {
                            let true_outcome = tab.measure(q);
                            if true_outcome == Some(true) {
                                tab.x(q);
                            }
                            let recorded = true_outcome.map(|b| tab.flip_with_prob(b, noise));
                            results.push(recorded);
                        }
                    }
                    MeasureName::MX
                    | MeasureName::MY
                    | MeasureName::MRX
                    | MeasureName::MRY
                    | MeasureName::MXX
                    | MeasureName::MYY
                    | MeasureName::MZZ
                    | MeasureName::MPP => {
                        unreachable!(
                            "unsupported measure {name:?} should have been rejected by prepare"
                        )
                    }
                }
            }
            ExtendedInstruction::MPad { bits, prob, .. } => {
                let noise = prob.unwrap_or(0.0);
                for &bit in bits {
                    results.push(Some(tab.flip_with_prob(bit, noise)));
                }
            }
            ExtendedInstruction::Raw(RawPassthrough::Annotation { .. }) => { /* no-op */ }
            ExtendedInstruction::Repeat { count, body, .. } => {
                for _ in 0..*count {
                    execute_prepared(body, tab, results);
                }
            }
        }
    }
}
