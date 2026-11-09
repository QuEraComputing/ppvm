// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Custom circuit messages that don't fit vihaco's standard stack protocols.

use smallvec::SmallVec;
use vihaco::PopStack;
use vihaco::word::Word;

/// Qubit batch plus a fixed parameter block; `P = ()` for plain batches.
pub struct QubitBatch<P = ()> {
    pub qubits: SmallVec<[u64; 8]>,
    pub params: P,
}

pub struct TwoQubitPauliErrorMsg {
    pub q0: u64,
    pub q1: u64,
    pub ps: [f64; 15],
}

pub struct TraceMessage {
    pub pauli_pattern: String,
}

/// Fixed-size parameter blocks popped off the operand stack.
pub trait PopParams: Sized {
    fn pop_params<S: PopStack<f64>>(source: &mut S) -> eyre::Result<Self>;
}

impl PopParams for () {
    fn pop_params<S: PopStack<f64>>(_: &mut S) -> eyre::Result<Self> {
        Ok(())
    }
}

impl PopParams for f64 {
    fn pop_params<S: PopStack<f64>>(source: &mut S) -> eyre::Result<Self> {
        source.pop()
    }
}

impl<const N: usize> PopParams for [f64; N] {
    fn pop_params<S: PopStack<f64>>(source: &mut S) -> eyre::Result<Self> {
        let mut params = [0.0; N];
        // Last pushed is on top, so fill from the back to keep push order.
        for p in params.iter_mut().rev() {
            *p = source.pop()?;
        }
        Ok(params)
    }
}

/// Pop a [`QubitBatch`] pushed as `q0 .. q(n-1)`, then `n`, then the params.
///
/// Call this from a composite's `BuildMessage` impl for batched instructions.
pub fn pop_qubit_batch<S, P>(stack: &mut S) -> eyre::Result<QubitBatch<P>>
where
    S: PopStack<Word> + PopStack<f64>,
    P: PopParams,
{
    let params = P::pop_params(stack)?;
    let n = usize::try_from(PopStack::<Word>::pop(stack)?)?;
    let mut qubits = SmallVec::with_capacity(n);
    for _ in 0..n {
        qubits.push(PopStack::<Word>::pop(stack)?);
    }
    qubits.reverse();
    Ok(QubitBatch { qubits, params })
}
