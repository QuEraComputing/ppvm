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
use ppvm_tableau::prelude::*;
use smallvec::SmallVec;
use stim_parser_2::prelude::{
    Axis, ExtendedInstruction, ExtendedProgram, GateName, GateOp, MeasureName, MeasureOp, MppOp,
    NoiseName, NoiseOp, PauliAxis, Target,
};

use crate::validate::{ExecError, validate};

/// Unwrap a plain qubit target. The `validate` pass guarantees that only the
/// control slot of a controlled Pauli may be a `rec[...]`, so every target
/// reached here through a non-control path is a qubit.
#[inline]
fn qubit(t: Target) -> usize {
    t.as_qubit()
        .expect("non-control gate targets are validated as qubits by validate")
}

// Inline capacity for the gate-target buffers below. Sized so typical narrow
// and moderate-width gate instructions stay on the stack; wider broadcasts
// spill to the heap (same as a plain `Vec`), so this is never worse.
const TARGETS_INLINE: usize = 16;

/// Collect plain qubit targets, for the batched single-qubit gate paths.
/// Stack-allocates (no per-instruction heap alloc) for up to [`TARGETS_INLINE`]
/// targets, since `*_many` takes a `&[usize]` and the source is `&[Target]`.
fn qubits(targets: &[Target]) -> SmallVec<[usize; TARGETS_INLINE]> {
    targets.iter().map(|&t| qubit(t)).collect()
}

/// Collect (control, target) qubit pairs, for the batched two-qubit gate paths.
/// Only valid when no target is a measurement record (see [`has_record_control`]).
fn qubit_pairs(targets: &[Target]) -> SmallVec<[(usize, usize); TARGETS_INLINE / 2]> {
    targets
        .chunks_exact(2)
        .map(|p| (qubit(p[0]), qubit(p[1])))
        .collect()
}

/// Whether any target is a measurement-record control `rec[-k]`. Such a gate
/// takes the per-pair feed-forward path rather than the batched fast path.
fn has_record_control(targets: &[Target]) -> bool {
    targets.iter().any(|t| matches!(t, Target::Rec(_)))
}

/// Resolve a measurement-record lookback `rec[-k]` against the running record.
/// `k == 1` is the most recent measurement. A lost-qubit measurement (`None`)
/// resolves to `false` (no feed-forward applied). An out-of-range lookback also
/// resolves to `false`, but `validate` rejects those up front, so on a validated
/// program the only `false`-from-missing case is a lost-qubit measurement.
#[inline]
fn record_bit(record: &[Option<bool>], k: usize) -> bool {
    record
        .len()
        .checked_sub(k)
        .and_then(|i| record.get(i).copied().flatten())
        .unwrap_or(false)
}

/// Measure qubit `q` in the Z basis, reset it to |0>, and return the reported
/// (readout-noise-flipped) outcome — the building block for `MR`/`MRX`/`MRY`.
///
/// The reset must act on the *true* (pre-flip) outcome, so this cannot delegate
/// to `measure_noisy`; instead it applies the same `flip_with_prob` draw and
/// overwrites the record entry with the reported value, so the measurement
/// record matches the returned result exactly (consistent with `M`).
fn measure_reset_z<T, I, C>(
    tab: &mut GeneralizedTableau<T, I, C>,
    q: usize,
    noise: f64,
) -> Option<bool>
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
    let true_outcome = tab.measure(q);
    if true_outcome == Some(true) {
        tab.x(q);
    }
    let recorded = true_outcome.map(|b| tab.flip_with_prob(b, noise));
    tab.overwrite_last_measurement_record(recorded);
    recorded
}

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
    instructions: &[ExtendedInstruction],
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
    instructions: &[ExtendedInstruction],
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
    instructions: &[ExtendedInstruction],
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
    instructions: &[ExtendedInstruction],
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
            ExtendedInstruction::Gate(GateOp { name, targets, .. }) => match name {
                GateName::Reset | GateName::ResetZ => {
                    targets.iter().for_each(|&t| tab.reset(qubit(t)));
                }
                // RX: reset to |+> (Z-basis reset, then rotate into the X basis).
                GateName::ResetX => targets.iter().for_each(|&t| {
                    let q = qubit(t);
                    tab.reset(q);
                    tab.h(q);
                }),
                // RY: reset to |i> (|0> -> |+> -> S|+> = |i>).
                GateName::ResetY => targets.iter().for_each(|&t| {
                    let q = qubit(t);
                    tab.reset(q);
                    tab.h(q);
                    tab.s(q);
                }),
                GateName::X => tab.x_many(&qubits(targets)),
                GateName::Y => tab.y_many(&qubits(targets)),
                GateName::Z => tab.z_many(&qubits(targets)),
                GateName::H | GateName::HXZ => tab.h_many(&qubits(targets)),
                GateName::S | GateName::SqrtZ => tab.s_many(&qubits(targets)),
                GateName::SDag | GateName::SqrtZDag => tab.s_dag_many(&qubits(targets)),
                GateName::SqrtX => tab.sqrt_x_many(&qubits(targets)),
                GateName::SqrtXDag => tab.sqrt_x_dag_many(&qubits(targets)),
                GateName::SqrtY => tab.sqrt_y_many(&qubits(targets)),
                GateName::SqrtYDag => tab.sqrt_y_dag_many(&qubits(targets)),
                GateName::Identity => {}
                // Controlled Paulis. The control slot may be a measurement record
                // `rec[-k]` (classical feed-forward): apply the target Pauli iff
                // the recorded bit is 1, exactly as Stim's `single_cx`/`single_cy`.
                // With no record control present, keep the batched fast path.
                GateName::CX | GateName::ZCX | GateName::CNot => {
                    if has_record_control(targets) {
                        for p in targets.chunks_exact(2) {
                            match p[0] {
                                Target::Qubit(c) => tab.cnot(c, qubit(p[1])),
                                Target::Rec(k) => {
                                    if record_bit(&tab.measurement_record, k) {
                                        tab.x(qubit(p[1]));
                                    }
                                }
                            }
                        }
                    } else {
                        tab.cnot_many(&qubit_pairs(targets));
                    }
                }
                GateName::CY | GateName::ZCY => {
                    if has_record_control(targets) {
                        for p in targets.chunks_exact(2) {
                            match p[0] {
                                Target::Qubit(c) => tab.cy(c, qubit(p[1])),
                                Target::Rec(k) => {
                                    if record_bit(&tab.measurement_record, k) {
                                        tab.y(qubit(p[1]));
                                    }
                                }
                            }
                        }
                    } else {
                        tab.cy_many(&qubit_pairs(targets));
                    }
                }
                GateName::CZ | GateName::ZCZ => {
                    if has_record_control(targets) {
                        for p in targets.chunks_exact(2) {
                            match p[0] {
                                Target::Qubit(c) => tab.cz(c, qubit(p[1])),
                                Target::Rec(k) => {
                                    if record_bit(&tab.measurement_record, k) {
                                        tab.z(qubit(p[1]));
                                    }
                                }
                            }
                        }
                    } else {
                        tab.cz_many(&qubit_pairs(targets));
                    }
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
                GateName::T | GateName::TDag => {
                    unreachable!("T/T_DAG are lowered to ExtendedInstruction::T/TDag by interpret")
                }
            },
            ExtendedInstruction::T { targets, .. } => targets.iter().for_each(|&q| tab.t(q)),
            ExtendedInstruction::TDag { targets, .. } => {
                targets.iter().for_each(|&q| tab.t_dag(q));
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
            ExtendedInstruction::Noise(NoiseOp {
                name,
                targets,
                args,
                ..
            }) => match name {
                NoiseName::Depolarize1 => {
                    debug_assert_eq!(args.len(), 1);
                    let p = args[0];
                    for &q in targets {
                        tab.depolarize1(q, p.into());
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
            ExtendedInstruction::Measure(MeasureOp {
                name,
                args,
                targets,
                ..
            }) => {
                let noise = args.first().copied().unwrap_or(0.0);
                match name {
                    MeasureName::M | MeasureName::MZ => {
                        // Noiseless readout: measure the whole target list at once
                        // so the shared scratch is reused across targets, mirroring
                        // the batched gate paths above. With readout noise, keep the
                        // per-qubit loop so each measure/flip pair's RNG draws stay
                        // interleaved exactly as before.
                        if noise > 0.0 {
                            for &q in targets {
                                results.push(tab.measure_noisy(q, noise));
                            }
                        } else {
                            results.extend(tab.measure_many(targets));
                        }
                    }
                    // MR cannot delegate to `measure_noisy` because the reset must use
                    // the *true* outcome — but it shares the readout flip via
                    // `flip_with_prob`, so the RNG-draw shape matches MZ exactly.
                    MeasureName::MR => {
                        for &q in targets {
                            results.push(measure_reset_z(tab, q, noise));
                        }
                    }
                    // X/Y-basis measurements via Stim's basis-change decomposition:
                    // conjugate so the measured axis maps to Z (`MX = H M H`,
                    // `MY = S_DAG H M H S`), do the Z-basis measurement, then
                    // conjugate back. Reuses `measure_noisy` so the readout-flip /
                    // RNG-draw shape matches `M` exactly.
                    MeasureName::MX => {
                        for &q in targets {
                            tab.h(q);
                            results.push(tab.measure_noisy(q, noise));
                            tab.h(q);
                        }
                    }
                    MeasureName::MY => {
                        for &q in targets {
                            tab.s_dag(q);
                            tab.h(q);
                            results.push(tab.measure_noisy(q, noise));
                            tab.h(q);
                            tab.s(q);
                        }
                    }
                    // MRX/MRY add a reset to the measured eigenstate
                    // (`MRX = H M R H`, `MRY = S_DAG H M R H S`). The reset must
                    // use the *true* outcome, so it mirrors the MR arm rather than
                    // delegating to `measure_noisy`.
                    MeasureName::MRX => {
                        for &q in targets {
                            tab.h(q);
                            results.push(measure_reset_z(tab, q, noise));
                            tab.h(q);
                        }
                    }
                    MeasureName::MRY => {
                        for &q in targets {
                            tab.s_dag(q);
                            tab.h(q);
                            results.push(measure_reset_z(tab, q, noise));
                            tab.h(q);
                            tab.s(q);
                        }
                    }
                    MeasureName::MXX | MeasureName::MYY | MeasureName::MZZ | MeasureName::MPP => {
                        unreachable!(
                            "unsupported measure {name:?} should have been rejected by validate"
                        )
                    }
                }
            }
            ExtendedInstruction::MPad { bits, prob, .. } => {
                let noise = prob.unwrap_or(0.0);
                for &bit in bits {
                    let recorded = Some(tab.flip_with_prob(bit, noise));
                    tab.append_measurement_record(recorded);
                    results.push(recorded);
                }
            }
            // MPP: measure each Pauli product non-destructively via Stim's
            // basis-change + CX-ladder gadget. Each factor's axis is rotated to
            // Z (`X = H · Z · H`, `Y = (S H) · Z · (H S_DAG)`), a CX ladder onto
            // the first qubit maps the product `Z_0 Z_1 ... Z_{m-1}` to a single
            // `Z_0`, that qubit is measured, then the ladder and basis changes
            // are undone so only the product operator is projected.
            ExtendedInstruction::Mpp(MppOp { products, args, .. }) => {
                let noise = args.first().copied().unwrap_or(0.0);
                for product in products {
                    for f in product {
                        match f.axis {
                            PauliAxis::X => tab.h(f.qubit),
                            PauliAxis::Y => {
                                tab.s_dag(f.qubit);
                                tab.h(f.qubit);
                            }
                            PauliAxis::Z => {}
                        }
                    }
                    let q0 = product[0].qubit;
                    for f in &product[1..] {
                        tab.cnot(f.qubit, q0);
                    }
                    results.push(tab.measure_noisy(q0, noise));
                    for f in product[1..].iter().rev() {
                        tab.cnot(f.qubit, q0);
                    }
                    for f in product {
                        match f.axis {
                            PauliAxis::X => tab.h(f.qubit),
                            PauliAxis::Y => {
                                tab.h(f.qubit);
                                tab.s(f.qubit);
                            }
                            PauliAxis::Z => {}
                        }
                    }
                }
            }
            ExtendedInstruction::Annotation(_) => { /* no-op */ }
            ExtendedInstruction::Repeat { count, body, .. } => {
                for _ in 0..*count {
                    execute_validated(body, tab, results);
                }
            }
        }
    }
}
