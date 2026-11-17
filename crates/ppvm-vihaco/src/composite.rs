// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::component::{Circuit, CircuitOp};
use crate::measurements::{
    CircuitOutcomeEffect, MeasurementObserver, MeasurementResult, TraceObserver,
};
use eyre;
use vihaco::{
    Apply, BuildMessage, Effects, Execute, ExecutionState, HasInstructionSet, PushToStack,
    ReusableMessageRecipe, Step, Type, Value, composite, field, module::LocalModule,
};
use vihaco_stdlib::debugger::Debugger;
use vihaco_stdlib::debugger::{DebuggerState, HasBreakpointInstruction, dialect::debug};

// Device configuration lives in `device_info` (shipped with the circuit
// component so it compiles without the composite). Re-exported here so the
// existing `crate::composite::{…}` paths keep resolving.
pub use crate::device_info::{BackendKind, PPVM_MAGIC, PPVMDeviceInfo};

#[derive(Default)]
#[composite]
pub struct PPVM {
    #[device(0)]
    #[program]
    execution: vihaco_stdlib::ExecutionState<vihaco_stdlib::debugger::DebuggerState>,

    #[device(1)]
    cpu: vihaco_cpu::CPU,

    // NOTE: we need the full import path here for some vihaco macro magic reasons
    #[device(2)]
    circuit: crate::component::Circuit,

    #[device(3)]
    debugger: vihaco_stdlib::debugger::Debugger<u32, Instruction>,

    measurement_record: MeasurementObserver,
    trace_record: TraceObserver,
}

type Execution = ExecutionState<DebuggerState>;
type Instruction = <PPVM as HasInstructionSet>::Runtime;

/// A lowered PPVM program; `extra` carries the header's device configuration.
pub type PPVMModule = LocalModule<Instruction, Value, Type, PPVMDeviceInfo>;

type Program = PPVMModule;

impl<I, M> BuildMessage<field!(PPVM::debugger), I, Program, M> for PPVM
where
    vihaco_stdlib::debugger::Debugger<u32, Instruction>: Execute<I, Message = M>,
    M: ReusableMessageRecipe<Self, field!(PPVM::debugger), I, Program>,
{
    fn build_message(&mut self, instruction: &I, module: &Program) -> eyre::Result<M> {
        M::prepare(self, instruction, module)
    }
}

// A breakpoint suspends the invocation via the debugger's `Break` effect.
impl Step<field!(PPVM::debugger), debug::Breakpoint, Program> for PPVM {
    fn step(&mut self, instruction: &debug::Breakpoint, module: &Program) -> eyre::Result<()> {
        let message = <Self as BuildMessage<field!(PPVM::debugger), debug::Breakpoint, Program>>::build_message(
            self, instruction, module,
        )?;
        let effect = self.debugger.execute(instruction, message)?;
        self.execution.apply(effect);
        Ok(())
    }
}

impl HasBreakpointInstruction for Instruction {
    fn new_breakpoint() -> Self {
        Self::V_debugger__debug__Breakpoint(debug::Breakpoint)
    }
}

impl PPVM {
    /// A fresh machine for the device described by `info`, ready to execute.
    pub fn new(info: &PPVMDeviceInfo, seed: Option<u64>) -> eyre::Result<Self> {
        if info.n_qubits == 0 {
            eyre::bail!("circuit.n_qubits must be declared");
        }
        Ok(Self {
            execution: ExecutionState::default(),
            cpu: vihaco_cpu::CPU::default(),
            circuit: Circuit::new(info, seed)?,
            measurement_record: MeasurementObserver::default(),
            trace_record: TraceObserver::default(),
            debugger: Debugger::default(),
        })
    }

    pub fn measurement_record(&self) -> &[MeasurementResult] {
        &self.measurement_record.record
    }

    pub fn trace_record(&self) -> &[f64] {
        &self.trace_record.record
    }

    /// Consume the machine, returning its measurement and trace records.
    pub fn into_records(self) -> (Vec<MeasurementResult>, Vec<f64>) {
        (self.measurement_record.record, self.trace_record.record)
    }
}

// Every field builds its messages through the message type's recipe.
impl<I, M> BuildMessage<field!(PPVM::execution), I, Program, M> for PPVM
where
    Execution: Execute<I, Message = M>,
    M: ReusableMessageRecipe<Self, field!(PPVM::execution), I, Program>,
{
    fn build_message(&mut self, instruction: &I, module: &Program) -> eyre::Result<M> {
        M::prepare(self, instruction, module)
    }
}

impl<I, M> BuildMessage<field!(PPVM::cpu), I, Program, M> for PPVM
where
    vihaco_cpu::CPU: Execute<I, Message = M>,
    M: ReusableMessageRecipe<Self, field!(PPVM::cpu), I, Program>,
{
    fn build_message(&mut self, instruction: &I, module: &Program) -> eyre::Result<M> {
        M::prepare(self, instruction, module)
    }
}

impl<I, M> BuildMessage<field!(PPVM::circuit), I, Program, M> for PPVM
where
    Circuit: Execute<I, Message = M>,
    M: ReusableMessageRecipe<Self, field!(PPVM::circuit), I, Program>,
{
    fn build_message(&mut self, instruction: &I, module: &Program) -> eyre::Result<M> {
        M::prepare(self, instruction, module)
    }
}

// execution stepping
impl<I> Step<field!(PPVM::execution), I, Program> for PPVM
where
    Self: BuildMessage<field!(PPVM::execution), I, Program>,
    Execution: Execute<I, Effect = (), Error = eyre::Report>,
{
    fn step(&mut self, instruction: &I, module: &Program) -> eyre::Result<()> {
        let message = <Self as BuildMessage<field!(PPVM::execution), I, Program>>::build_message(
            self,
            instruction,
            module,
        )?;
        self.execution.execute(instruction, message)
    }
}

// CPU stepping
impl<I> Step<field!(PPVM::cpu), I, Program> for PPVM
where
    Self: BuildMessage<field!(PPVM::cpu), I, Program>,
    vihaco_cpu::CPU: Execute<I, Error = eyre::Report>,
    Execution: vihaco::Apply<<vihaco_cpu::CPU as Execute<I>>::Effect>,
{
    fn step(&mut self, instruction: &I, module: &Program) -> eyre::Result<()> {
        let message = <Self as BuildMessage<field!(PPVM::cpu), I, Program>>::build_message(
            self,
            instruction,
            module,
        )?;
        let effect = self.cpu.execute(instruction, message)?;
        vihaco::Apply::apply(&mut self.execution, effect);
        Ok(())
    }
}

// handling circuit effects
impl vihaco::Apply<Effects<CircuitOutcomeEffect>> for PPVM {
    fn apply(&mut self, effects: Effects<CircuitOutcomeEffect>) {
        for effect in effects {
            // Programs branch on outcomes: push each as u32 (0/1/2, 2 = Lost),
            // and can consume trace values the same way.
            match &effect {
                CircuitOutcomeEffect::Measurement(m) => {
                    for &outcome in &m.measurement_results {
                        self.execution.apply(PushToStack((outcome as u32,)));
                    }
                }
                CircuitOutcomeEffect::Trace(t) => {
                    self.execution.apply(PushToStack((t.value,)));
                }
            }
            self.record_outcome(effect);
        }
    }
}

/// Record a circuit effect without pushing it onto the operand stack.
pub trait Record<E> {
    fn record(&mut self, effect: E);
}

impl Record<()> for PPVM {
    fn record(&mut self, _effect: ()) {}
}

impl Record<Effects<CircuitOutcomeEffect>> for PPVM {
    fn record(&mut self, effects: Effects<CircuitOutcomeEffect>) {
        for effect in effects {
            self.record_outcome(effect);
        }
    }
}

impl PPVM {
    fn record_outcome(&mut self, effect: CircuitOutcomeEffect) {
        match effect {
            CircuitOutcomeEffect::Measurement(m) => {
                self.measurement_record.record.push(m.measurement_results)
            }
            CircuitOutcomeEffect::Trace(t) => self.trace_record.record.push(t.value),
        }
    }

    /// Execute one circuit instruction on a ready-made message, outside any
    /// running program.
    ///
    /// The operand stack is not touched: outcomes land in the measurement and
    /// trace records only. Qubit indices are not range-checked here; callers
    /// building messages by hand should check them against `n_qubits`.
    pub fn send<I>(&mut self, instruction: &I, message: I::Message) -> eyre::Result<()>
    where
        I: CircuitOp,
        Self: Record<I::Effect>,
    {
        let effect = self.circuit.execute(instruction, message)?;
        self.record(effect);
        Ok(())
    }

    /// Number of qubits on the circuit device.
    pub fn n_qubits(&self) -> usize {
        self.circuit.n_qubits()
    }

    /// Human-readable rendering of the circuit state.
    pub fn state_string(&self) -> String {
        self.circuit.state_string()
    }

    /// Enter `module`'s main without executing anything; `resume` runs it.
    pub fn start(&mut self, module: &Program) -> eyre::Result<()> {
        self.execution.start(module, &[])
    }

    /// Run the started program until it halts, returns, or hits a breakpoint.
    pub fn resume(&mut self, module: &Program) -> eyre::Result<()> {
        vihaco_stdlib::debugger::resume_with_breaks::<field!(PPVM::debugger), _, _, _, _, _>(
            module,
            self,
            vihaco_stdlib::debugger::DebugAction::Continue,
        )
    }

    /// Start `module` and run it to the end, continuing past every breakpoint.
    pub fn run(&mut self, module: &Program) -> eyre::Result<()> {
        self.start(module)?;
        // Breakpoints only pause the interactive debugger; a batch run skips past them.
        self.resume(module)?;
        while self.is_paused() {
            self.resume(module)?;
        }
        Ok(())
    }

    pub fn single_step(&mut self, module: &Program) -> eyre::Result<()> {
        vihaco_stdlib::debugger::resume_with_breaks::<field!(PPVM::debugger), _, _, _, _, _>(
            module,
            self,
            vihaco_stdlib::debugger::DebugAction::Step,
        )
    }

    /// Index of the next instruction to execute.
    pub fn pc(&self) -> u32 {
        let pc = self.execution.pc();
        match self.execution.run_state() {
            // An explicit breakpoint leaves `pc` on itself; resuming skips past it.
            // A gutter breakpoint resumes by running the instruction it displaced.
            vihaco_stdlib::RunState::Other(DebuggerState::Breakpoint { resume_pc })
                if !self.debugger.contains(&pc) =>
            {
                *resume_pc
            }
            _ => pc,
        }
    }

    /// True while a breakpoint has suspended the running program.
    pub fn is_paused(&self) -> bool {
        matches!(
            self.execution.run_state(),
            vihaco_stdlib::RunState::Other(DebuggerState::Breakpoint { .. })
        )
    }

    /// True once the program has finished (or before one was started).
    pub fn is_halted(&self) -> bool {
        self.execution.is_halted()
    }
}

impl vihaco::Apply<()> for PPVM {
    fn apply(&mut self, _effect: ()) {}
}

// Circuit stepping
impl<I> Step<field!(PPVM::circuit), I, Program> for PPVM
where
    Self: BuildMessage<field!(PPVM::circuit), I, Program>,
    Circuit: Execute<I, Error = eyre::Report>,
    Self: vihaco::Apply<<Circuit as Execute<I>>::Effect>,
{
    fn step(&mut self, instruction: &I, module: &Program) -> eyre::Result<()> {
        let message = <Self as BuildMessage<field!(PPVM::circuit), I, Program>>::build_message(
            self,
            instruction,
            module,
        )?;
        let effect = <Self as vihaco::Field<field!(PPVM::circuit)>>::get_mut(self)
            .execute(instruction, message)?;
        vihaco::Apply::apply(self, effect);
        Ok(())
    }
}

#[cfg(test)]
#[allow(dead_code)]
fn assert_dispatch_coverage(
    m: &mut PPVM,
    inst: &Instruction,
    module: &Program,
) -> eyre::Result<()> {
    vihaco::Dispatch::<Program>::dispatch(m, inst, module)
}

impl vihaco::Reset for PPVM {
    fn reset(&mut self) {
        self.cpu.reset();
        self.circuit.reset();
        self.measurement_record.record.clear();
        self.trace_record.record.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::loader::compile_program;
    use crate::measurements::MeasurementOutcome;
    use vihaco::{PopFromStack, PushStack, Reset};
    use vihaco_circuit_isa::circuit;
    use vihaco_circuit_isa::messages::TraceMessage;

    /// Wrap header declarations and a straight-line `main` into SST v1 source.
    /// `strings` are appended to the string table after `main` and the file name.
    fn program(header: &str, strings: &[&str], body: &[&str]) -> String {
        let strings: Vec<String> = ["main", "test.sst"]
            .iter()
            .chain(strings)
            .map(|s| format!("{s:?}"))
            .collect();
        format!(
            "sst v1
.section(root):
    .section(ppvm):
    .header(ppvm):
        {header}
        metadata {{
            strings: [{strings}],
            functions: [
                {{
                    name: 0,
                    signature: {{ params: [], ret: [] }},
                    local_count: 0,
                    start_address: 0,
                    end_address: {end},
                    file: 1,
                }},
            ],
            labels: [],
            constants: [],
            main_function: 0,
            source_symbols: [],
            extra: {{}},
        }};
    .header(ppvm).
    .text(ppvm):
        fn @main() {{
            {body}
        }}
    .text(ppvm).
    .section(ppvm).
.section(root).
",
            strings = strings.join(", "),
            end = body.len(),
            body = body.join("\n            "),
        )
    }

    /// Compile `source` and run it to the end on a fresh machine.
    fn run_source(source: &str) -> eyre::Result<PPVM> {
        let module = compile_program(source)?;
        let mut machine = PPVM::new(&module.extra, None)?;
        machine.run(&module)?;
        Ok(machine)
    }

    fn info(n_qubits: usize) -> PPVMDeviceInfo {
        PPVMDeviceInfo {
            n_qubits,
            ..Default::default()
        }
    }

    fn paulisum_info(n_qubits: usize, observable: Option<&str>) -> PPVMDeviceInfo {
        PPVMDeviceInfo {
            n_qubits,
            backend: BackendKind::PauliSum,
            observable: observable.map(str::to_owned),
            ..Default::default()
        }
    }

    /// Bell pair on q0, q1, then measure both.
    const BELL: &[&str] = &[
        "execution::program::const_u64 0",
        "circuit::circuit::h",
        "execution::program::const_u64 0",
        "execution::program::const_u64 1",
        "circuit::circuit::cnot",
        "execution::program::const_u64 0",
        "circuit::circuit::measure",
        "execution::program::const_u64 1",
        "circuit::circuit::measure",
        "execution::program::ret 0",
    ];

    /// Pause, then measure q0 (|0> -> 0 unless a gate is injected at the pause).
    const BREAK_THEN_MEASURE: &[&str] = &[
        "debugger::debug::breakpoint",
        "execution::program::const_u64 0",
        "circuit::circuit::measure",
        "execution::program::ret 0",
    ];

    /// The outcomes recorded so far, one per single-qubit measurement.
    fn outcomes(machine: &PPVM) -> Vec<MeasurementOutcome> {
        machine
            .measurement_record()
            .iter()
            .flatten()
            .copied()
            .collect()
    }

    /// Start `source` and run it until its first breakpoint.
    fn paused_at_breakpoint(source: &str) -> eyre::Result<(PPVM, Program)> {
        let module = compile_program(source)?;
        let mut machine = PPVM::new(&module.extra, None)?;
        machine.start(&module)?;
        machine.resume(&module)?;
        assert!(
            machine.is_paused(),
            "program should pause at its breakpoint"
        );
        Ok((machine, module))
    }

    // ─── Construction ─────────────────────────────────────────────────────

    #[test]
    fn new_rejects_zero_qubits() {
        assert!(PPVM::new(&info(0), None).is_err());
    }

    #[test]
    fn new_fails_when_n_qubits_undeclared() -> eyre::Result<()> {
        let module = compile_program(&program("", &[], &["execution::program::ret 0"]))?;
        let err = PPVM::new(&module.extra, None)
            .err()
            .expect("construction should fail");
        assert!(err.to_string().contains("circuit.n_qubits"), "err: {err}");
        Ok(())
    }

    #[test]
    fn state_string_renders_the_circuit() -> eyre::Result<()> {
        let machine = PPVM::new(&info(2), None)?;
        let rendered = machine.state_string();
        assert!(
            !rendered.is_empty(),
            "state_string should render the tableau"
        );
        assert_eq!(rendered, machine.circuit.state_string());
        Ok(())
    }

    #[test]
    fn paulisum_rejects_missing_observable() {
        let err = PPVM::new(&paulisum_info(1, None), None)
            .err()
            .expect("construction should fail");
        assert!(err.to_string().contains("observable"), "err: {err}");
    }

    #[test]
    fn paulisum_rejects_mismatched_observable_length() {
        // Three letters but only two qubits.
        let err = PPVM::new(&paulisum_info(2, Some("ZZZ")), None)
            .err()
            .expect("construction should fail");
        assert!(err.to_string().contains("invalid Pauli-sum"), "err: {err}");
    }

    // ─── Compiling and running source ─────────────────────────────────────

    #[test]
    fn compile_program_reads_device_info_and_code() -> eyre::Result<()> {
        let source = program(
            "circuit.n_qubits 2;\n        circuit.coefficient_threshold 1e-8;",
            &[],
            &[
                "execution::program::const_u64 0",
                "circuit::circuit::h",
                "execution::program::ret 0",
            ],
        );
        let module = compile_program(&source)?;
        assert_eq!(module.extra.n_qubits, 2);
        assert_eq!(module.extra.coefficient_threshold, 1e-8);
        assert_eq!(module.code.len(), 3);
        Ok(())
    }

    #[test]
    fn compile_program_reports_parse_errors() {
        let source = program(
            "circuit.n_qubits 2;",
            &[],
            &[
                "circuit::circuit::not_a_real_gate",
                "execution::program::ret 0",
            ],
        );
        assert!(compile_program(&source).is_err());
    }

    #[test]
    fn missing_operand_is_an_error() {
        // `h` with an empty stack has no qubit to pop; this must surface as `Err`.
        let source = program(
            "circuit.n_qubits 1;",
            &[],
            &["circuit::circuit::h", "execution::program::ret 0"],
        );
        assert!(run_source(&source).is_err());
    }

    #[test]
    fn bell_outcomes_agree() -> eyre::Result<()> {
        let machine = run_source(&program("circuit.n_qubits 2;", &[], BELL))?;
        let outcomes = outcomes(&machine);
        assert_eq!(outcomes.len(), 2);
        assert_eq!(outcomes[0], outcomes[1], "Bell outcomes are correlated");
        Ok(())
    }

    #[test]
    fn ghz_outcomes_agree() -> eyre::Result<()> {
        // 5-qubit GHZ: H on q0, CNOT down the chain, then measure every qubit.
        let mut body = vec![
            "execution::program::const_u64 0".to_owned(),
            "circuit::circuit::h".to_owned(),
        ];
        for i in 0..4 {
            body.push(format!("execution::program::const_u64 {i}"));
            body.push(format!("execution::program::const_u64 {}", i + 1));
            body.push("circuit::circuit::cnot".to_owned());
        }
        for q in 0..5 {
            body.push(format!("execution::program::const_u64 {q}"));
            body.push("circuit::circuit::measure".to_owned());
        }
        body.push("execution::program::ret 0".to_owned());
        let body: Vec<&str> = body.iter().map(String::as_str).collect();

        let header = "circuit.n_qubits 5;\n        circuit.coefficient_threshold 1e-10;";
        let machine = run_source(&program(header, &[], &body))?;
        let outcomes = outcomes(&machine);
        assert_eq!(outcomes.len(), 5);
        assert!(
            outcomes.iter().all(|&o| o == outcomes[0]),
            "GHZ outcomes all agree: {outcomes:?}"
        );
        Ok(())
    }

    #[test]
    fn reset_clears_the_records() -> eyre::Result<()> {
        let mut machine = run_source(&program("circuit.n_qubits 2;", &[], BELL))?;
        assert_eq!(machine.measurement_record().len(), 2);
        machine.reset();
        assert!(machine.measurement_record().is_empty());
        assert!(machine.trace_record().is_empty());
        Ok(())
    }

    #[test]
    fn two_qubit_operands_pop_in_push_order() -> eyre::Result<()> {
        // Operands are pushed in argument order (q0, q1, then floats), so every
        // two-qubit gate reads q0 as the first value pushed.
        let module = compile_program(&program(
            "circuit.n_qubits 8;",
            &[],
            &["execution::program::ret 0"],
        ))?;
        let mut m = PPVM::new(&module.extra, None)?;
        m.start(&module)?;

        m.execution.push(2_u64);
        m.execution.push(5_u64);
        let PopFromStack(cnot) = <PPVM as BuildMessage<
            field!(PPVM::circuit),
            circuit::CNOT,
            Program,
        >>::build_message(&mut m, &circuit::CNOT, &module)?;
        assert_eq!(cnot, (2, 5));

        m.execution.push(2_u64);
        m.execution.push(5_u64);
        m.execution.push(0.3_f64);
        let PopFromStack(rxx) = <PPVM as BuildMessage<
            field!(PPVM::circuit),
            circuit::RXX,
            Program,
        >>::build_message(&mut m, &circuit::RXX, &module)?;
        assert_eq!(rxx, (2, 5, 0.3));

        m.execution.push(2_u64);
        m.execution.push(5_u64);
        for p in [0.1_f64, 0.2, 0.3] {
            m.execution.push(p);
        }
        let PopFromStack(loss) = <PPVM as BuildMessage<
            field!(PPVM::circuit),
            circuit::CorrelatedLoss,
            Program,
        >>::build_message(
            &mut m, &circuit::CorrelatedLoss, &module
        )?;
        assert_eq!(loss, (2, 5, 0.1, 0.2, 0.3));

        // Leakage takes one qubit and two floats.
        m.execution.push(2_u64);
        m.execution.push(0.1_f64);
        m.execution.push(0.2_f64);
        let PopFromStack(leakage) = <PPVM as BuildMessage<
            field!(PPVM::circuit),
            circuit::Leakage,
            Program,
        >>::build_message(&mut m, &circuit::Leakage, &module)?;
        assert_eq!(leakage, (2, 0.1, 0.2));
        assert!(m.execution.stack().is_empty(), "every operand was consumed");
        Ok(())
    }

    // ─── Breakpoints ──────────────────────────────────────────────────────

    #[test]
    fn run_ignores_breakpoints() -> eyre::Result<()> {
        // A breakpoint between the two Bell measurements must not cut the run short.
        let mut body = BELL.to_vec();
        body.insert(7, "debugger::debug::breakpoint");
        let machine = run_source(&program("circuit.n_qubits 2;", &[], &body))?;
        let outcomes = outcomes(&machine);
        assert_eq!(outcomes.len(), 2);
        assert_eq!(outcomes[0], outcomes[1]);
        Ok(())
    }

    #[test]
    fn single_step_advances_past_breakpoint() -> eyre::Result<()> {
        let source = program("circuit.n_qubits 1;", &[], BREAK_THEN_MEASURE);
        let (mut m, module) = paused_at_breakpoint(&source)?;
        let pc_at_break = m.pc();
        m.single_step(&module)?;
        assert!(!m.is_paused(), "stepping must leave the breakpoint");
        assert!(m.pc() > pc_at_break, "pc must advance");
        Ok(())
    }

    // ─── Sending instructions directly ────────────────────────────────────

    #[test]
    fn send_persists_state_across_calls() -> eyre::Result<()> {
        let mut m = PPVM::new(&info(1), None)?;
        m.send(&circuit::X, PopFromStack((0,)))?;
        assert!(m.measurement_record().is_empty());
        // The X from the first call persists, so q0 measures 1.
        m.send(&circuit::Measure, PopFromStack((0,)))?;
        assert_eq!(outcomes(&m), [MeasurementOutcome::One]);
        Ok(())
    }

    #[test]
    fn send_measure_does_not_grow_the_stack() -> eyre::Result<()> {
        // Only running bytecode consumes outcomes from the stack, so a sent
        // measurement records its outcome without pushing it.
        let mut m = PPVM::new(&info(1), None)?;
        m.send(&circuit::Measure, PopFromStack((0,)))?;
        assert!(m.execution.stack().is_empty());
        assert_eq!(m.measurement_record().len(), 1);
        Ok(())
    }

    #[test]
    fn send_at_breakpoint_preserves_position() -> eyre::Result<()> {
        let source = program("circuit.n_qubits 1;", &[], BREAK_THEN_MEASURE);
        let (mut m, module) = paused_at_breakpoint(&source)?;
        let pc = m.pc();
        let depth = m.execution.stack().len();

        m.send(&circuit::X, PopFromStack((0,)))?;
        assert_eq!(m.pc(), pc, "pc must be preserved");
        assert_eq!(m.execution.stack().len(), depth, "stack must be unchanged");
        assert!(m.is_paused());

        // The injected X took effect: resuming measures |1>.
        m.resume(&module)?;
        assert!(m.is_halted());
        assert_eq!(outcomes(&m), [MeasurementOutcome::One]);
        Ok(())
    }

    #[test]
    fn failed_send_at_breakpoint_preserves_position() -> eyre::Result<()> {
        let source = program("circuit.n_qubits 1;", &[], BREAK_THEN_MEASURE);
        let (mut m, _) = paused_at_breakpoint(&source)?;
        let pc = m.pc();
        let depth = m.execution.stack().len();

        let bad = TraceMessage {
            pauli_pattern: "not a pattern".to_owned(),
        };
        assert!(m.send(&circuit::Trace, bad).is_err());
        assert_eq!(m.pc(), pc, "pc must be preserved on error");
        assert_eq!(m.execution.stack().len(), depth);
        assert!(m.is_paused());
        Ok(())
    }

    // ─── Backends: truncate, trace, reset ─────────────────────────────────

    #[test]
    fn tableau_truncate_is_silent_no_op() -> eyre::Result<()> {
        let source = program(
            "circuit.n_qubits 1;",
            &[],
            &["circuit::circuit::truncate", "execution::program::ret 0"],
        );
        let machine = run_source(&source)?;
        assert!(machine.measurement_record().is_empty());
        assert!(machine.trace_record().is_empty());
        Ok(())
    }

    #[test]
    fn paulisum_truncate_runs_without_error() -> eyre::Result<()> {
        let header =
            "circuit.n_qubits 1;\n        circuit.backend paulisum;\n        circuit.observable Z;";
        let source = program(
            header,
            &[],
            &["circuit::circuit::truncate", "execution::program::ret 0"],
        );
        run_source(&source)?;
        Ok(())
    }

    /// Run a trace of `pattern` on `header`'s device and return the traces.
    fn trace(header: &str, pattern: &str) -> eyre::Result<Vec<f64>> {
        let push = format!("execution::program::const_string {pattern:?}");
        let source = program(
            header,
            &[pattern],
            &[
                &push,
                "circuit::circuit::trace",
                "execution::program::ret 0",
            ],
        );
        Ok(run_source(&source)?.trace_record().to_vec())
    }

    #[test]
    fn tableau_trace_emits_expectation_on_zero_state() -> eyre::Result<()> {
        // <0|Z|0> = 1.
        let traces = trace("circuit.n_qubits 1;", "Z0")?;
        assert_eq!(traces.len(), 1);
        assert!((traces[0] - 1.0).abs() < 1e-12, "got {}", traces[0]);
        Ok(())
    }

    #[test]
    fn paulisum_trace_populates_trace_record() -> eyre::Result<()> {
        // The seeded observable is Z; tracing Z0 picks up that one unit term.
        let header =
            "circuit.n_qubits 1;\n        circuit.backend paulisum;\n        circuit.observable Z;";
        assert_eq!(trace(header, "Z0")?, [1.0]);
        Ok(())
    }

    #[test]
    fn paulisum_multi_term_observable_seeds_all_terms() -> eyre::Result<()> {
        // 1.0*ZZ + 0.5*XX: a pattern matching both words sums to 1.5.
        let header = "circuit.n_qubits 2;\n        circuit.backend paulisum;\n        circuit.observable 1.0*ZZ+0.5*XX;";
        assert_eq!(trace(header, "[XZ]0[XZ]1")?, [1.5]);
        Ok(())
    }

    /// H(0) changes the propagated observable; reset must restore the seeded one.
    fn assert_reset_restores_seeded_observable(backend: &str) -> eyre::Result<()> {
        let header = format!(
            "circuit.n_qubits 1;\n        circuit.backend {backend};\n        circuit.observable Z;"
        );
        let source = program(
            &header,
            &[],
            &[
                "execution::program::const_u64 0",
                "circuit::circuit::h",
                "execution::program::ret 0",
            ],
        );
        let module = compile_program(&source)?;
        let mut machine = PPVM::new(&module.extra, None)?;
        let seeded = machine.state_string();
        machine.run(&module)?;
        assert_ne!(
            machine.state_string(),
            seeded,
            "H(0) should change the state"
        );
        machine.reset();
        assert_eq!(machine.state_string(), seeded);
        Ok(())
    }

    #[test]
    fn paulisum_reset_restores_seeded_observable() -> eyre::Result<()> {
        assert_reset_restores_seeded_observable("paulisum")
    }

    #[test]
    fn lossy_paulisum_reset_restores_seeded_observable() -> eyre::Result<()> {
        assert_reset_restores_seeded_observable("lossy_paulisum")
    }
}
