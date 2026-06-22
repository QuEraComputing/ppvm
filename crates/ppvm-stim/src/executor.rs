// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::BitView;
use itertools::Itertools;
use num::Integer;
use num::PrimInt;
use num::complex::{Complex64, ComplexFloat};
use num::{Complex, One, ToPrimitive, Zero};
use std::fmt::Debug;

use ppvm_pauli_sum::prelude::*;
use ppvm_tableau::measure_all::LossyMeasureAll;
use ppvm_tableau::prelude::*;
use stim_parser::ast::{GateName, MeasureName, NoiseName};
use stim_parser::extended::{Axis, ExtendedInstruction, ExtendedProgram, RawPassthrough};

use crate::validate::{ExecError, validate};

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
        + PartialOrd
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
    validate(program)?;
    let count = program.measurement_count();
    let mut results = Vec::with_capacity(count);
    execute_validated(&program.instructions, tab, &mut results);
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
        + PartialOrd
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
    validate(program)?;
    Ok(sample_serial_validated(
        &program.instructions,
        program.measurement_count(),
        num_shots,
        make_tableau,
    ))
}

/// Like [`sample_serial`] but skips validation — call only when the program
/// has already been validated (e.g. via [`validate`](fn@validate)).
pub fn sample_serial_validated<T, I, C, F>(
    instructions: &[stim_parser::extended::ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
) -> Vec<Vec<Option<bool>>>
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
        + PartialOrd
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
    (0..num_shots)
        .map(|i| {
            let mut tab = make_tableau(i);
            let mut results = Vec::with_capacity(measurement_count);
            execute_validated(instructions, &mut tab, &mut results);
            results
        })
        .collect()
}

/// Execute many shots, building a fresh tableau for shot `i` via
/// `make_tableau(i)`.
///
/// When the `rayon` feature is enabled this dispatches to `sample_parallel`
/// for batches large enough to amortise thread-scheduling overhead, and to
/// [`sample_serial`] otherwise. Without the feature it is always serial. Use
/// [`sample_serial`] / `sample_parallel` directly to force one or the other.
#[cfg(feature = "rayon")]
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
        + PartialOrd
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
    validate(program)?;
    Ok(sample_validated(
        &program.instructions,
        program.measurement_count(),
        num_shots,
        make_tableau,
    ))
}

/// Like [`sample`] but skips validation — call only when the program has
/// already been validated (e.g. via [`validate`](fn@validate)).
///
/// When the `rayon` feature is enabled this dispatches to
/// [`sample_parallel_validated`] for large batches, falling back to
/// [`sample_serial_validated`] for small ones.
#[cfg(feature = "rayon")]
pub fn sample_validated<T, I, C, F>(
    instructions: &[stim_parser::extended::ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
) -> Vec<Vec<Option<bool>>>
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
        + PartialOrd
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
    // Below ~4 shots per thread, rayon's scheduling overhead outweighs
    // the gain, so stay serial; with a single thread there is no upside.
    let n_threads = rayon::current_num_threads();
    if n_threads <= 1 || num_shots < 4 * n_threads {
        sample_serial_validated(instructions, measurement_count, num_shots, make_tableau)
    } else {
        sample_parallel_validated(instructions, measurement_count, num_shots, make_tableau)
    }
}

/// Like [`sample`] but skips validation — call only when the program has
/// already been validated (e.g. via [`validate`](fn@validate)). Without the `rayon`
/// feature this always runs serially.
#[cfg(not(feature = "rayon"))]
pub fn sample_validated<T, I, C, F>(
    instructions: &[stim_parser::extended::ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
) -> Vec<Vec<Option<bool>>>
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
        + PartialOrd
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
    sample_serial_validated(instructions, measurement_count, num_shots, make_tableau)
}

/// Execute many shots, building a fresh tableau for shot `i` via
/// `make_tableau(i)`.
///
/// Without the `rayon` feature execution is always serial. Enable the
/// `rayon` feature to get a version that can dispatch to parallel execution
/// for large batches. Use [`sample_serial`] to force serial execution
/// regardless of the feature.
#[cfg(not(feature = "rayon"))]
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
        + PartialOrd
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
    sample_serial(program, num_shots, make_tableau)
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
        + PartialOrd
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
    validate(program)?;
    Ok(sample_parallel_validated(
        &program.instructions,
        program.measurement_count(),
        num_shots,
        make_tableau,
    ))
}

/// Like [`sample_parallel`] but skips validation — call only when the program
/// has already been validated (e.g. via [`validate`](fn@validate)).
#[cfg(feature = "rayon")]
pub fn sample_parallel_validated<T, I, C, F>(
    instructions: &[stim_parser::extended::ExtendedInstruction],
    measurement_count: usize,
    num_shots: usize,
    make_tableau: F,
) -> Vec<Vec<Option<bool>>>
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
        + PartialOrd
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
    (0..num_shots)
        .into_par_iter()
        .map(|i| {
            let mut tab = make_tableau(i);
            let mut results = Vec::with_capacity(measurement_count);
            execute_validated(instructions, &mut tab, &mut results);
            results
        })
        .collect()
}

/// Dispatch a slice of validated instructions onto a tableau, appending
/// measurement bits to `results`. Skips validation — the caller is responsible
/// for having run [`validate`](fn@validate) on the originating program. Used by
/// [`execute`] / [`sample`] internally and by the Python `tab.run()` path
/// where `StimProgram` already cached the validate step.
pub fn execute_validated<T, I, C>(
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
        + PartialOrd
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
                GateName::X => tab.x_batch(targets),
                GateName::Y => tab.y_batch(targets),
                GateName::Z => tab.z_batch(targets),
                GateName::H | GateName::HXZ => tab.h_batch(targets),
                GateName::S | GateName::SqrtZ => tab.s_batch(targets),
                GateName::SDag | GateName::SqrtZDag => tab.s_adj_batch(targets),
                GateName::SqrtX => tab.sqrt_x_batch(targets),
                GateName::SqrtXDag => tab.sqrt_x_adj_batch(targets),
                GateName::SqrtY => tab.sqrt_y_batch(targets),
                GateName::SqrtYDag => tab.sqrt_y_adj_batch(targets),
                GateName::Identity => {}
                GateName::CX | GateName::ZCX | GateName::CNot => {
                    let mut pairs: Vec<(usize, usize)> = Vec::with_capacity(targets.len() / 2);
                    pairs.extend(targets.chunks_exact(2).map(|p| (p[0], p[1])));
                    tab.cnot_batch(&pairs);
                }
                GateName::CY | GateName::ZCY => {
                    let mut pairs: Vec<(usize, usize)> = Vec::with_capacity(targets.len() / 2);
                    pairs.extend(targets.chunks_exact(2).map(|p| (p[0], p[1])));
                    tab.cy_batch(&pairs);
                }
                GateName::CZ | GateName::ZCZ => {
                    let mut pairs: Vec<(usize, usize)> = Vec::with_capacity(targets.len() / 2);
                    pairs.extend(targets.chunks_exact(2).map(|p| (p[0], p[1])));
                    tab.cz_batch(&pairs);
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
                    unreachable!("unsupported gate {name:?} should have been rejected by validate")
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
                    unreachable!("unsupported noise {name:?} should have been rejected by validate")
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
                        // Noiseless readout: batch the measurement so the shared
                        // scratch is reused across targets, mirroring the batched
                        // gate paths above. With readout noise, keep the per-qubit
                        // loop so each measure/flip pair's RNG draws stay
                        // interleaved exactly as before.
                        if noise > 0.0 {
                            for &q in targets {
                                results.push(tab.measure_noisy(q, noise));
                            }
                        } else {
                            results.extend(tab.measure_batch(targets));
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
                            "unsupported measure {name:?} should have been rejected by validate"
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
                    execute_validated(body, tab, results);
                }
            }
        }
    }
}
