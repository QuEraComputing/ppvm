// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Custom circuit messages that don't fit vihaco's standard stack protocols.

use smallvec::SmallVec;
use vihaco::runtime::{AutomaticMessage, PrepareMessage};
use vihaco::word::Word;
use vihaco::{Composite, Execute, Field, HasOperandStack, PopStack};

/// Qubit batch plus a fixed parameter block; `P = ()` for plain batches.
pub struct QubitBatch<P = ()> {
    pub qubits: SmallVec<[usize; 8]>,
    pub params: P,
}

/// Batch of qubit pairs plus a fixed parameter block; `P = ()` for plain batches.
pub struct QubitPairBatch<P = ()> {
    pub pairs: SmallVec<[(usize, usize); 4]>,
    pub params: P,
}

/// Single two-qubit Pauli channel: pushed as `q0 q1`, then the 15 probabilities.
pub struct TwoQubitPauliErrorMsg {
    pub q0: usize,
    pub q1: usize,
    pub ps: [f64; 15],
}

/// Pop a [`TwoQubitPauliErrorMsg`] pushed as `q0 q1`, then the 15 probabilities.
pub fn pop_two_qubit_pauli_error<S>(stack: &mut S) -> eyre::Result<TwoQubitPauliErrorMsg>
where
    S: PopStack<Word> + PopStack<f64>,
{
    let ps = <[f64; 15]>::pop_params(stack)?;
    let q1 = usize::try_from(PopStack::<Word>::pop(stack)?)?;
    let q0 = usize::try_from(PopStack::<Word>::pop(stack)?)?;
    Ok(TwoQubitPauliErrorMsg { q0, q1, ps })
}

/// Pauli pattern to trace against, pushed as a `const_string` id.
pub struct TraceMessage {
    pub pauli_pattern: String,
}

/// Read access to a program's string table, used to resolve string operands.
pub trait StringTable {
    fn string(&self, id: usize) -> eyre::Result<&str>;
}

impl<I, V, Ty, Info> StringTable for vihaco::module::LocalModule<I, V, Ty, Info> {
    fn string(&self, id: usize) -> eyre::Result<&str> {
        self.strings
            .get(id)
            .map(String::as_str)
            .ok_or_else(|| eyre::eyre!("unknown string id {id}"))
    }
}

/// Pop a [`TraceMessage`] whose pattern is a `const_string` id into `strings`.
pub fn pop_trace<S: PopStack<Word>>(
    stack: &mut S,
    strings: &(impl StringTable + ?Sized),
) -> eyre::Result<TraceMessage> {
    let id = usize::try_from(stack.pop()?)?;
    Ok(TraceMessage {
        pauli_pattern: strings.string(id)?.to_owned(),
    })
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
        let q = usize::try_from(PopStack::<Word>::pop(stack)?)?;
        qubits.push(q);
    }
    qubits.reverse();
    Ok(QubitBatch { qubits, params })
}

/// Pop a [`QubitPairBatch`] pushed as `a0 b0 .. a(n-1) b(n-1)`, then the pair
/// count `n`, then the params.
///
/// Call this from a composite's `BuildMessage` impl for batched instructions.
pub fn pop_qubit_pair_batch<S, P>(stack: &mut S) -> eyre::Result<QubitPairBatch<P>>
where
    S: PopStack<Word> + PopStack<f64>,
    P: PopParams,
{
    let params = P::pop_params(stack)?;
    let n = usize::try_from(PopStack::<Word>::pop(stack)?)?;
    let mut pairs = SmallVec::with_capacity(n);
    for _ in 0..n {
        let b = usize::try_from(PopStack::<Word>::pop(stack)?)?;
        let a = usize::try_from(PopStack::<Word>::pop(stack)?)?;
        pairs.push((a, b));
    }
    pairs.reverse();
    Ok(QubitPairBatch { pairs, params })
}

// TODO: these use vihaco's `#[doc(hidden)]` message traits so the batches build in
// any composite; switch to a public extension point once QuEraComputing/vihaco#169 lands.
macro_rules! automatic_batch_message {
    ($msg:ident, $pop:ident) => {
        impl<P> AutomaticMessage for $msg<P> {}

        impl<C, F, I, Ctx: ?Sized, P> PrepareMessage<C, F, I, Ctx> for $msg<P>
        where
            C: Composite + Field<F> + HasOperandStack<I>,
            C::Child: Execute<I, Message = $msg<P>>,
            C::Source: PopStack<Word> + PopStack<f64>,
            P: PopParams,
        {
            fn prepare(parent: &mut C, _: &I, _: &Ctx) -> eyre::Result<Self> {
                $pop(parent.operand_stack(1))
            }
        }
    };
}

automatic_batch_message!(QubitBatch, pop_qubit_batch);
automatic_batch_message!(QubitPairBatch, pop_qubit_pair_batch);

impl AutomaticMessage for TwoQubitPauliErrorMsg {}

impl<C, F, I, Ctx: ?Sized> PrepareMessage<C, F, I, Ctx> for TwoQubitPauliErrorMsg
where
    C: Composite + Field<F> + HasOperandStack<I>,
    C::Child: Execute<I, Message = TwoQubitPauliErrorMsg>,
    C::Source: PopStack<Word> + PopStack<f64>,
{
    fn prepare(parent: &mut C, _: &I, _: &Ctx) -> eyre::Result<Self> {
        pop_two_qubit_pauli_error(parent.operand_stack(17))
    }
}

impl AutomaticMessage for TraceMessage {}

impl<C, F, I, Ctx> PrepareMessage<C, F, I, Ctx> for TraceMessage
where
    C: Composite + Field<F> + HasOperandStack<I>,
    C::Child: Execute<I, Message = TraceMessage>,
    C::Source: PopStack<Word>,
    Ctx: StringTable + ?Sized,
{
    fn prepare(parent: &mut C, _: &I, strings: &Ctx) -> eyre::Result<Self> {
        pop_trace(parent.operand_stack(1), strings)
    }
}
