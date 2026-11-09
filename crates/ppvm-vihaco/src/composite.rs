// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::component::Circuit;
use crate::measurements::{
    CircuitOutcomeEffect, MeasurementObserver, MeasurementResult, TraceObserver,
};
use eyre;
use vihaco::{
    BuildMessage, Effects, Execute, ExecutionState, HasInstructionSet, PushToStack, Step, Type,
    Value, composite, field, module::LocalModule,
};

// Device configuration lives in `device_info` (shipped with the circuit
// component so it compiles without the composite). Re-exported here so the
// existing `crate::composite::{…}` paths keep resolving.
pub use crate::device_info::{BackendKind, PPVM_MAGIC, PPVMDeviceInfo};

#[derive(Default)]
#[composite]
pub struct PPVM {
    #[device(0)]
    #[program]
    execution: vihaco_stdlib::ExecutionState,

    #[device(1)]
    cpu: vihaco_cpu::CPU,

    // NOTE: we need the full import path here for some vihaco macro magic reasons
    #[device(2)]
    circuit: crate::component::Circuit,

    measurement_record: MeasurementObserver,

    trace_record: TraceObserver,
}

type Instruction = <PPVM as HasInstructionSet>::Runtime;

/// A lowered PPVM program; `extra` carries the header's device configuration.
pub type PPVMModule = LocalModule<Instruction, Value, Type, PPVMDeviceInfo>;

type Program = PPVMModule;

impl PPVM {
    /// A fresh machine for the device described by `info`, ready to execute.
    pub fn new(info: &PPVMDeviceInfo, seed: Option<u64>) -> eyre::Result<Self> {
        if info.n_qubits == 0 {
            eyre::bail!("circuit.n_qubits must be declared");
        }
        Ok(Self {
            execution: ExecutionState::new(),
            cpu: vihaco_cpu::CPU::default(),
            circuit: Circuit::new(info, seed)?,
            measurement_record: MeasurementObserver::default(),
            trace_record: TraceObserver::default(),
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

// Batched circuit instructions pop their operands through custom messages.
// TODO: per-composite builders conflict with vihaco's blanket `BuildMessage` (E0119)
// since the messages are foreign here; the ISA crate implements them for now instead.
// impl<I, P: PopParams> BuildMessage<field!(PPVM::circuit), I, Program, QubitBatch<P>> for PPVM
// where
//     Circuit: Execute<I, Message = QubitBatch<P>>,
// {
//     fn build_message(&mut self, _: &I, _: &Program) -> eyre::Result<QubitBatch<P>> {
//         pop_qubit_batch(&mut self.execution)
//     }
// }

// impl<I, P: PopParams> BuildMessage<field!(PPVM::circuit), I, Program, QubitPairBatch<P>> for PPVM
// where
//     Circuit: Execute<I, Message = QubitPairBatch<P>>,
// {
//     fn build_message(&mut self, _: &I, _: &Program) -> eyre::Result<QubitPairBatch<P>> {
//         pop_qubit_pair_batch(&mut self.execution)
//     }
// }

// execution stepping
impl<I> Step<field!(PPVM::execution), I, Program> for PPVM
where
    Self: BuildMessage<field!(PPVM::execution), I, Program>,
    ExecutionState: Execute<I, Effect = (), Error = eyre::Report>,
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
    ExecutionState: vihaco::Apply<<vihaco_cpu::CPU as Execute<I>>::Effect>,
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
            match effect {
                CircuitOutcomeEffect::Measurement(m) => {
                    // Programs branch on outcomes: push each as u32 (0/1/2, 2 = Lost).
                    for &outcome in &m.measurement_results {
                        self.execution.apply(PushToStack((outcome as u32,)));
                    }
                    self.measurement_record.record.push(m.measurement_results)
                }
                CircuitOutcomeEffect::Trace(t) => {
                    // Programs can consume the expectation value, as with outcomes.
                    self.execution.apply(PushToStack((t.value,)));
                    self.trace_record.record.push(t.value);
                }
            }
        }
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

// #[cfg(test)]
// mod tests {
//     use vihaco::{Type, Value, module::LocalModule};

//     use super::*;

//     #[test]
//     fn test_run_ppvm() -> eyre::Result<()> {
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();

//         module.extra.n_qubits = 2;

//         /*
//         cpu::cpu.const u64, 0
//         circuit::circuit.h
//         */
//         let zero = PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//             Type::U64,
//             Value::U64(0),
//         ));
//         let one = PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//             Type::U64,
//             Value::U64(1),
//         ));
//         module.code.push(zero.clone());
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::H));

//         /*
//         cpu::cpu.const u64, 0
//         circuit::circuit.t
//         */
//         module.code.push(zero.clone());
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::T));

//         /*
//         cpu::cpu.const u64, 0
//         cpu::cpu.const u64, 1
//         circuit::circuit.cnot
//         */
//         module.code.push(zero.clone());
//         module.code.push(one.clone());
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::CNOT));

//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;

//         for _ in 0..module.code.len() {
//             machine.step_once()?;
//             assert!(machine.cpu.stack().len() <= 2);
//         }

//         let num_coefficients = match &machine.circuit {
//             Circuit::Tableau(TableauCircuit::Bits64(ex)) => {
//                 println!("{}", ex.tab);
//                 ex.tab.coefficients.len()
//             }
//             Circuit::Tableau(TableauCircuit::Bits128(ex)) => {
//                 println!("{}", ex.tab);
//                 ex.tab.coefficients.len()
//             }
//             Circuit::Tableau(TableauCircuit::Bits256(ex)) => {
//                 println!("{}", ex.tab);
//                 ex.tab.coefficients.len()
//             }
//             Circuit::Tableau(TableauCircuit::Bits512(ex)) => {
//                 println!("{}", ex.tab);
//                 ex.tab.coefficients.len()
//             }
//             Circuit::Tableau(TableauCircuit::Bits1024(ex)) => {
//                 println!("{}", ex.tab);
//                 ex.tab.coefficients.len()
//             }
//             Circuit::Tableau(TableauCircuit::Bits2048(ex)) => {
//                 println!("{}", ex.tab);
//                 ex.tab.coefficients.len()
//             }
//             Circuit::PauliSum(_) | Circuit::LossyPauliSum(_) => {
//                 panic!("test expects the default Tableau backend, got a PauliSum variant");
//             }
//         };

//         assert_eq!(num_coefficients, 2);
//         Ok(())
//     }

//     #[test]
//     fn test_device_decl() -> eyre::Result<()> {
//         // Equivalent .sst source:
//         //
//         //     device circuit.n_qubits 5;
//         //     device circuit.coefficient_threshold 1e-10;
//         //
//         //     fn @main() { ...5-qubit GHZ + 5 measurements... }

//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();

//         module.extra.n_qubits = 5;
//         module.extra.coefficient_threshold = 1e-10;

//         // 5-qubit GHZ: H on q0, then CNOT(q_i, q_{i+1}) for i = 0..4.
//         /*
//         cpu::cpu.const u64, 0
//         circuit::circuit.h
//         */
//         module
//             .code
//             .push(PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                 Type::U64,
//                 Value::U64(0),
//             )));
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::H));

//         for i in 0..4u64 {
//             /*
//             cpu::cpu.const u64, i
//             cpu::cpu.const u64, i+1
//             circuit::circuit.cnot
//             */
//             module
//                 .code
//                 .push(PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                     Type::U64,
//                     Value::U64(i),
//                 )));
//             module
//                 .code
//                 .push(PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                     Type::U64,
//                     Value::U64(i + 1),
//                 )));
//             module
//                 .code
//                 .push(PPVMInstruction::Circuit(CircuitInstruction::CNOT));
//         }

//         // Measure all 5 qubits.
//         for q in 0..5u64 {
//             /*
//             cpu::cpu.const u64, q
//             circuit::circuit.measure
//             */
//             module
//                 .code
//                 .push(PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                     Type::U64,
//                     Value::U64(q),
//                 )));
//             module
//                 .code
//                 .push(PPVMInstruction::Circuit(CircuitInstruction::Measure));
//         }

//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;

//         for _ in 0..module.code.len() {
//             machine.step_once()?;
//         }

//         assert_eq!(machine.measurement_record().len(), 5);
//         Ok(())
//     }

//     // ─── Incremental execution (REPL) ─────────────────────────────────────

//     #[test]
//     fn execute_single_instruction_persists_state_across_calls() -> eyre::Result<()> {
//         use crate::measurements::MeasurementOutcome;

//         // A 1-qubit device with no code; the REPL builds up instructions
//         // incrementally, one command at a time, rather than loading a program.
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 1;

//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;

//         // First command: X on q0 (|0> -> |1>).
//         let x = [
//             PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                 Type::U64,
//                 Value::U64(0),
//             )),
//             PPVMInstruction::Circuit(CircuitInstruction::X),
//         ];
//         machine.execute_single_instruction(&x)?;
//         // No measurement yet.
//         assert!(machine.measurement_record().is_empty());

//         // Second command: measure q0. The X from the first command must persist,
//         // so the outcome is deterministically |1>.
//         let measure = [
//             PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                 Type::U64,
//                 Value::U64(0),
//             )),
//             PPVMInstruction::Circuit(CircuitInstruction::Measure),
//         ];
//         machine.execute_single_instruction(&measure)?;

//         let record = machine.measurement_record();
//         assert_eq!(record.len(), 1);
//         assert_eq!(record[0].as_slice(), [MeasurementOutcome::One]);
//         Ok(())
//     }

//     #[test]
//     fn execute_single_instruction_propagates_engine_errors() -> eyre::Result<()> {
//         // The REPL relies on engine errors surfacing as `Err` (so it can print
//         // them and keep looping) rather than panicking. A circuit with no qubit
//         // operand on the stack is one such propagating error.
//         //
//         // NOTE: an out-of-range qubit index (>= n_qubits) currently *panics* in
//         // the tableau rather than erroring, so the REPL command layer must
//         // bounds-check qubit indices before calling `execute`.
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 1;

//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;

//         // `circuit::circuit.h` with nothing on the stack: `pop_u64` fails.
//         let missing_operand = [PPVMInstruction::Circuit(CircuitInstruction::H)];
//         assert!(
//             machine
//                 .execute_single_instruction(&missing_operand)
//                 .is_err()
//         );
//         Ok(())
//     }

//     #[test]
//     fn state_string_renders_a_small_device() -> eyre::Result<()> {
//         let source = "device circuit.n_qubits 2;\nfn @main() { cpu::cpu.ret 0 }\n";
//         let mut machine = PPVM::default();
//         machine.load_program(source)?;
//         machine.init()?;

//         let rendered = machine.state_string();
//         assert!(
//             !rendered.is_empty(),
//             "state_string should render the tableau"
//         );
//         // PPVM delegates to the circuit.
//         assert_eq!(rendered, machine.circuit.state_string());
//         Ok(())
//     }

//     #[test]
//     fn resolve_circuit_pops_operands_in_reverse_of_push_order() -> eyre::Result<()> {
//         // Convention: operands are pushed in argument order (q0, q1, then any
//         // floats) and popped in reverse. So every two-qubit circuit must read q0 as
//         // the first operand pushed, consistently, with or without trailing
//         // floats. (CNOT already obeyed this; the float-carrying arms did not.)
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 8;
//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;

//         // CNOT: push q0=2, q1=5.
//         machine.cpu.stack_push(Value::U64(2));
//         machine.cpu.stack_push(Value::U64(5));
//         assert_eq!(
//             machine.resolve_circuit(&CircuitInstruction::CNOT)?,
//             CircuitMessage::TwoQubit(2, 5)
//         );

//         // RXX: push q0=2, q1=5, theta — same qubit order as CNOT.
//         machine.cpu.stack_push(Value::U64(2));
//         machine.cpu.stack_push(Value::U64(5));
//         machine.cpu.stack_push(Value::F64(0.3));
//         assert_eq!(
//             machine.resolve_circuit(&CircuitInstruction::RXX)?,
//             CircuitMessage::TwoQubitAndFloat(2, 5, 0.3)
//         );

//         // CorrelatedLoss: push q0=2, q1=5, p0, p1, p2.
//         machine.cpu.stack_push(Value::U64(2));
//         machine.cpu.stack_push(Value::U64(5));
//         machine.cpu.stack_push(Value::F64(0.1));
//         machine.cpu.stack_push(Value::F64(0.2));
//         machine.cpu.stack_push(Value::F64(0.3));
//         assert_eq!(
//             machine.resolve_circuit(&CircuitInstruction::CorrelatedLoss)?,
//             CircuitMessage::TwoQubitAndFloatArr3(2, 5, [0.1, 0.2, 0.3])
//         );

//         // Leakage: push q=2, p0, p1 — two floats, unlike Loss.
//         machine.cpu.stack_push(Value::U64(2));
//         machine.cpu.stack_push(Value::F64(0.1));
//         machine.cpu.stack_push(Value::F64(0.2));
//         assert_eq!(
//             machine.resolve_circuit(&CircuitInstruction::Leakage)?,
//             CircuitMessage::QubitAndTwoFloats(2, 0.1, 0.2)
//         );
//         Ok(())
//     }

//     #[test]
//     fn with_qubits_builds_an_initialized_device() -> eyre::Result<()> {
//         use crate::measurements::MeasurementOutcome;

//         let mut m = PPVM::with_qubits(2)?;
//         // The device is ready to take instructions immediately.
//         m.apply_circuit_instruction(CircuitInstruction::X, &[0], &[])?;
//         m.apply_circuit_instruction(CircuitInstruction::Measure, &[0], &[])?;
//         let record = m.measurement_record();
//         assert_eq!(record.len(), 1);
//         assert_eq!(record[0].as_slice(), [MeasurementOutcome::One]);
//         Ok(())
//     }

//     #[test]
//     fn with_qubits_zero_is_an_error() {
//         assert!(PPVM::with_qubits(0).is_err());
//     }

//     #[test]
//     fn apply_circuit_instruction_bounds_checks_qubits() -> eyre::Result<()> {
//         // q1 is out of range on a 1-qubit device: this must error rather than
//         // panic in the tableau.
//         let mut m = PPVM::with_qubits(1)?;
//         let err = m
//             .apply_circuit_instruction(CircuitInstruction::X, &[1], &[])
//             .unwrap_err();
//         assert!(err.to_string().contains("out of range"), "got: {err}");
//         Ok(())
//     }

//     // ─── Parser-driven entry points ───────────────────────────────────────

//     #[test]
//     fn load_program_populates_device_info_and_code() -> eyre::Result<()> {
//         let source = "device circuit.n_qubits 2;\n\
//                       device circuit.coefficient_threshold 1e-8;\n\
//                       fn @main() {\n\
//                           cpu::cpu.const u64, 0\n\
//                           circuit::circuit.h\n\
//                           cpu::cpu.ret 0\n\
//                       }\n";
//         let mut machine = PPVM::default();
//         machine.load_program(source)?;
//         assert_eq!(machine.loader.module.extra.n_qubits, 2);
//         assert_eq!(machine.loader.module.extra.coefficient_threshold, 1e-8);
//         // cpu::cpu.const u64, 0 / circuit::circuit.h / ret = 3
//         assert_eq!(machine.loader.module.code.len(), 3);
//         Ok(())
//     }

//     #[test]
//     fn run_program_executes_bell_circuit() -> eyre::Result<()> {
//         let source = "device circuit.n_qubits 2;\n\
//                       fn @main() {\n\
//                           cpu::cpu.const u64, 0\n\
//                           circuit::circuit.h\n\
//                           cpu::cpu.const u64, 0\n\
//                           cpu::cpu.const u64, 1\n\
//                           circuit::circuit.cnot\n\
//                           cpu::cpu.const u64, 0\n\
//                           circuit::circuit.measure\n\
//                           cpu::cpu.const u64, 1\n\
//                           circuit::circuit.measure\n\
//                           cpu::cpu.ret 0\n\
//                       }\n";
//         let mut machine = PPVM::default();
//         machine.run_program(source)?;
//         let record = machine.measurement_record();
//         assert_eq!(record.len(), 2);
//         Ok(())
//     }

//     #[test]
//     fn reset_clears_the_measurement_record() -> eyre::Result<()> {
//         use vihaco::Reset;

//         let source = "device circuit.n_qubits 2;\n\
//                       fn @main() {\n\
//                           cpu::cpu.const u64, 0\n\
//                           circuit::circuit.h\n\
//                           cpu::cpu.const u64, 0\n\
//                           cpu::cpu.const u64, 1\n\
//                           circuit::circuit.cnot\n\
//                           cpu::cpu.const u64, 0\n\
//                           circuit::circuit.measure\n\
//                           cpu::cpu.const u64, 1\n\
//                           circuit::circuit.measure\n\
//                           cpu::cpu.ret 0\n\
//                       }\n";
//         let mut machine = PPVM::default();
//         machine.run_program(source)?;
//         assert_eq!(machine.measurement_record().len(), 2);

//         // Resetting the machine must discard the recorded measurements, so a
//         // subsequent run does not see stale results leaking in from before.
//         machine.reset();
//         assert!(
//             machine.measurement_record().is_empty(),
//             "reset must clear the measurement record"
//         );

//         Ok(())
//     }

//     #[test]
//     fn init_fails_when_n_qubits_undeclared() -> eyre::Result<()> {
//         let source = "fn @main() { cpu::cpu.ret 0 }\n";
//         let mut machine = PPVM::default();
//         machine.load_program(source)?;
//         let err = machine.init().unwrap_err();
//         assert!(err.to_string().contains("circuit.n_qubits"), "err: {err}");
//         Ok(())
//     }

//     #[test]
//     fn run_program_reports_parse_errors() {
//         let source = "device circuit.n_qubits 2;\n\
//                       fn @main() {\n\
//                           circuit::circuit.not_a_real_gate\n\
//                           cpu::cpu.ret 0\n\
//                       }\n";
//         let mut machine = PPVM::default();
//         let err = machine.run_program(source).unwrap_err();
//         assert!(
//             err.to_string().contains("parsing failed")
//                 || err.to_string().contains("parsing functions failed")
//                 || err.to_string().contains("unhandled raw form"),
//             "err: {err}"
//         );
//     }

//     // ─── Breakpoints ──────────────────────────────────────────────────────

//     /// Bell circuit with a `cpu::cpu.breakpoint` between the two measurements.
//     const BREAKPOINT_PROGRAM: &str = "device circuit.n_qubits 2;\n\
//                                       fn @main() {\n\
//                                           cpu::cpu.const u64, 0\n\
//                                           circuit::circuit.h\n\
//                                           cpu::cpu.const u64, 0\n\
//                                           cpu::cpu.const u64, 1\n\
//                                           circuit::circuit.cnot\n\
//                                           cpu::cpu.const u64, 0\n\
//                                           circuit::circuit.measure\n\
//                                           cpu::cpu.breakpoint\n\
//                                           cpu::cpu.const u64, 1\n\
//                                           circuit::circuit.measure\n\
//                                           cpu::cpu.ret 0\n\
//                                       }\n";

//     #[test]
//     fn run_ignores_breakpoints() -> eyre::Result<()> {
//         // A batch run must execute straight through the cpu::cpu.breakpoint and record
//         // both measurements, exactly as if it weren't there.
//         let mut machine = PPVM::default();
//         machine.run_program(BREAKPOINT_PROGRAM)?;
//         assert_eq!(machine.measurement_record().len(), 2);
//         Ok(())
//     }

//     #[test]
//     fn step_once_advances_past_breakpoint() -> eyre::Result<()> {
//         let mut machine = PPVM::default();
//         machine.load_program(BREAKPOINT_PROGRAM)?;
//         machine.init()?;

//         // Step until the cpu::cpu.breakpoint pauses us.
//         let mut outcome = StepOutcome::Continue;
//         for _ in 0..machine.loader.module.code.len() {
//             outcome = machine.step_once()?;
//             if outcome == StepOutcome::Breakpoint {
//                 break;
//             }
//         }
//         assert_eq!(
//             outcome,
//             StepOutcome::Breakpoint,
//             "cpu::cpu.breakpoint should pause"
//         );
//         let pc_at_break = machine.current_pc();

//         // Stepping again must make progress (advance the pc) rather than
//         // re-hitting the same cpu::cpu.breakpoint instruction.
//         let next = machine.step_once()?;
//         assert_ne!(
//             next,
//             StepOutcome::Breakpoint,
//             "must move past the cpu::cpu.breakpoint"
//         );
//         assert!(machine.current_pc() > pc_at_break, "pc must advance");

//         Ok(())
//     }

//     #[test]
//     fn paulisum_truncate_runs_without_error() -> eyre::Result<()> {
//         // Smoke test: a `circuit::circuit.truncate` reaches the PauliSum executor's
//         // Truncate arm and calls `state.truncate()`. Task 8 makes the
//         // observable mandatory for PauliSum init, so seed `Z` here.
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 1;
//         module.extra.backend = BackendKind::PauliSum;
//         module.extra.observable = Some("Z".to_string());
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::Truncate));

//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;
//         machine.step_once()?;
//         Ok(())
//     }

//     #[test]
//     fn paulisum_trace_populates_trace_record() -> eyre::Result<()> {
//         // End-to-end Trace pipeline: with the observable `Z` seeded (Task 8),
//         // tracing the `Z0` pattern picks up that one term with coefficient
//         // 1.0, so the trace should be exactly 1.0.
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 1;
//         module.extra.backend = BackendKind::PauliSum;
//         module.extra.observable = Some("Z".to_string());
//         // `Z0` matches a Z on qubit 0; the parser requires position anchors.
//         module.strings.push("Z0".to_string());

//         module
//             .code
//             .push(PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                 Type::String,
//                 Value::String(0),
//             )));
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::Trace));

//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;
//         for _ in 0..module.code.len() {
//             machine.step_once()?;
//         }

//         assert_eq!(machine.trace_record(), vec![1.0]);
//         Ok(())
//     }

//     #[test]
//     fn paulisum_multi_term_observable_seeds_all_terms() -> eyre::Result<()> {
//         // Task 11: a sum-valued observable seeds every term. With
//         // `"ZZ + 0.5*XX"` the state holds `1.0 * ZZ + 0.5 * XX`; tracing
//         // `[XZ]0[XZ]1` matches both words and returns 1.0 + 0.5 = 1.5.
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 2;
//         module.extra.backend = BackendKind::PauliSum;
//         module.extra.observable = Some("ZZ + 0.5*XX".to_string());
//         module.strings.push("[XZ]0[XZ]1".to_string());

//         module
//             .code
//             .push(PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                 Type::String,
//                 Value::String(0),
//             )));
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::Trace));

//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;
//         for _ in 0..module.code.len() {
//             machine.step_once()?;
//         }

//         assert_eq!(machine.trace_record(), vec![1.5]);
//         Ok(())
//     }

//     #[test]
//     fn paulisum_init_rejects_missing_observable() {
//         // Task 8 requires `device circuit.observable` for PauliSum / Lossy.
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 1;
//         module.extra.backend = BackendKind::PauliSum;

//         let mut machine = PPVM::default();
//         machine.load(&module).unwrap();
//         let err = machine.init().unwrap_err();
//         assert!(
//             err.to_string().contains("observable"),
//             "expected observable-related error, got: {err}"
//         );
//     }

//     #[test]
//     fn paulisum_init_rejects_mismatched_observable_length() {
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 2;
//         module.extra.backend = BackendKind::PauliSum;
//         // Three letters but only two qubits — should error.
//         module.extra.observable = Some("ZZZ".to_string());

//         let mut machine = PPVM::default();
//         machine.load(&module).unwrap();
//         let err = machine.init().unwrap_err();
//         assert!(
//             err.to_string().contains("invalid Pauli-sum"),
//             "expected parser rejection, got: {err}"
//         );
//     }

//     #[test]
//     fn tableau_truncate_is_silent_no_op() -> eyre::Result<()> {
//         // Task 9: `circuit::circuit.truncate` on the default Tableau backend should run
//         // without error — the tableau prunes via coefficient_threshold during
//         // every gate, so the explicit Truncate instruction has nothing to do.
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 1;
//         // backend defaults to Tableau; no observable needed.
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::Truncate));

//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;
//         machine.step_once()?;
//         // No observer effects emitted by Truncate.
//         assert!(machine.measurement_record().is_empty());
//         assert!(machine.trace_record().is_empty());
//         Ok(())
//     }

//     #[test]
//     fn apply_circuit_instruction_preserves_pc_and_code_len() -> eyre::Result<()> {
//         use crate::measurements::MeasurementOutcome;

//         // cpu::cpu.breakpoint; then measure q0. Step to the cpu::cpu.breakpoint, inject X, resume.
//         let src = "device circuit.n_qubits 1;\n\
//                    fn @main() { cpu::cpu.breakpoint\n cpu::cpu.const u64, 0\n circuit::circuit.measure\n cpu::cpu.ret 0 }\n";
//         let mut m = PPVM::default();
//         m.load_program(src)?;
//         m.init()?;

//         // Run until the cpu::cpu.breakpoint pauses us.
//         loop {
//             if m.step_once()? == StepOutcome::Breakpoint {
//                 break;
//             }
//         }
//         let pc = m.current_pc();
//         let len = m.loader.module.code.len();
//         let depth = m.cpu.stack().len();

//         // Inject X on q0 while "paused".
//         m.apply_circuit_instruction(CircuitInstruction::X, &[0], &[])?;

//         // The debugger's position must be untouched.
//         assert_eq!(m.current_pc(), pc, "pc must be preserved");
//         assert_eq!(
//             m.loader.module.code.len(),
//             len,
//             "appended op must be truncated back"
//         );
//         assert_eq!(
//             m.cpu.stack().len(),
//             depth,
//             "operand stack must be left unchanged"
//         );

//         // And the X took effect: resuming the program measures |1>.
//         while !matches!(m.step_once()?, StepOutcome::Return | StepOutcome::Halt) {}
//         let rec = m.measurement_record();
//         assert_eq!(rec.len(), 1);
//         assert_eq!(rec[0].as_slice(), [MeasurementOutcome::One]);
//         Ok(())
//     }

//     #[test]
//     fn apply_circuit_instruction_preserves_pc_and_code_on_error() -> eyre::Result<()> {
//         // A program stepped to a known pc; an injected gate that errors mid-block
//         // must NOT corrupt the code vector or the program counter.
//         let src = "device circuit.n_qubits 1;\n\
//                    fn @main() { cpu::cpu.breakpoint\n cpu::cpu.const u64, 0\n circuit::circuit.measure\n cpu::cpu.ret 0 }\n";
//         let mut m = PPVM::default();
//         m.load_program(src)?;
//         m.init()?;
//         loop {
//             if m.step_once()? == StepOutcome::Breakpoint {
//                 break;
//             }
//         }
//         let pc = m.current_pc();
//         let len = m.loader.module.code.len();
//         let depth = m.cpu.stack().len();

//         // RX with no float param errors during execution: resolve_circuit pops
//         // the f64 first and finds the qubit value (U64) instead, returning Err.
//         let err = m.apply_circuit_instruction(CircuitInstruction::RX, &[0], &[]);
//         assert!(err.is_err(), "expected the injected gate to error");

//         // The failed injection must leave the debugger untouched.
//         assert_eq!(m.current_pc(), pc, "pc must be preserved on error");
//         assert_eq!(
//             m.loader.module.code.len(),
//             len,
//             "code must be truncated back on error"
//         );
//         assert_eq!(
//             m.cpu.stack().len(),
//             depth,
//             "operand stack must be left unchanged on error"
//         );
//         Ok(())
//     }

//     #[test]
//     fn apply_circuit_instruction_measurement_does_not_grow_the_stack() -> eyre::Result<()> {
//         // `circuit::circuit.measure` pushes its outcome onto the CPU operand stack for
//         // bytecode to consume. An injected measurement has no such consumer, so
//         // that push must be rolled back: otherwise a paused program resumes with
//         // a stray operand, and a REPL session's stack grows without bound.
//         let mut m = PPVM::with_qubits(1)?;
//         let depth = m.cpu.stack().len();
//         m.apply_circuit_instruction(CircuitInstruction::Measure, &[0], &[])?;
//         assert_eq!(
//             m.cpu.stack().len(),
//             depth,
//             "injected measurement must not leave its outcome on the stack"
//         );
//         // The outcome still lands in the measurement record.
//         assert_eq!(m.measurement_record().len(), 1);
//         Ok(())
//     }

//     #[test]
//     fn tableau_trace_emits_expectation_on_zero_state() {
//         // Task 16: `circuit::circuit.trace` on the Tableau backend now computes
//         // Σ_{P matches pat} ⟨ψ|P|ψ⟩ via `GeneralizedTableau::trace`. On the
//         // freshly-initialized |0⟩ state, pattern `Z0` matches the single
//         // Pauli Z and ⟨0|Z|0⟩ = 1, so the trace_record gets one entry: 1.0.
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 1;
//         module.strings.push("Z0".to_string());
//         module
//             .code
//             .push(PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                 Type::String,
//                 Value::String(0),
//             )));
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::Trace));

//         let mut machine = PPVM::default();
//         machine.load(&module).unwrap();
//         machine.init().unwrap();
//         machine.step_once().unwrap(); // const.string
//         machine.step_once().unwrap(); // circuit::circuit.trace
//         let trace = machine.trace_record();
//         assert_eq!(trace.len(), 1);
//         assert!(
//             (trace[0] - 1.0).abs() < 1e-12,
//             "expected ⟨0|Z|0⟩ = 1.0, got {}",
//             trace[0]
//         );
//     }

//     #[test]
//     fn paulisum_reset_restores_seeded_observable() -> eyre::Result<()> {
//         use vihaco::Reset;

//         // Seed the observable `Z` (PauliSum backend), then apply H(0), which
//         // conjugates Z -> X in the Heisenberg picture and changes the state.
//         // reset() must rebuild the state from the seeded observable.
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 1;
//         module.extra.backend = BackendKind::PauliSum;
//         module.extra.observable = Some("Z".to_string());

//         module
//             .code
//             .push(PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                 Type::U64,
//                 Value::U64(0),
//             )));
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::H));

//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;

//         let seeded = machine.state_string();
//         for _ in 0..module.code.len() {
//             machine.step_once()?;
//         }
//         assert_ne!(
//             machine.state_string(),
//             seeded,
//             "H(0) should have changed the propagated observable"
//         );

//         machine.reset();
//         assert_eq!(
//             machine.state_string(),
//             seeded,
//             "reset must rebuild the state from the seeded observable"
//         );
//         Ok(())
//     }

//     #[test]
//     fn lossy_paulisum_reset_restores_seeded_observable() -> eyre::Result<()> {
//         use vihaco::Reset;

//         // Same as `paulisum_reset_restores_seeded_observable`, but through the
//         // LossyPauliSum dispatch path.
//         let mut module: LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo> =
//             LocalModule::default();
//         module.extra.n_qubits = 1;
//         module.extra.backend = BackendKind::LossyPauliSum;
//         module.extra.observable = Some("Z".to_string());

//         module
//             .code
//             .push(PPVMInstruction::Cpu(vihaco_cpu::RuntimeInstruction::Const(
//                 Type::U64,
//                 Value::U64(0),
//             )));
//         module
//             .code
//             .push(PPVMInstruction::Circuit(CircuitInstruction::H));

//         let mut machine = PPVM::default();
//         machine.load(&module)?;
//         machine.init()?;

//         let seeded = machine.state_string();
//         for _ in 0..module.code.len() {
//             machine.step_once()?;
//         }
//         assert_ne!(
//             machine.state_string(),
//             seeded,
//             "H(0) should have changed the propagated observable"
//         );

//         machine.reset();
//         assert_eq!(
//             machine.state_string(),
//             seeded,
//             "reset must rebuild the state from the seeded observable"
//         );
//         Ok(())
//     }
// }
