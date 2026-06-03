use vihaco::frame::Frame;
use vihaco::machine::StackFrame;
use vihaco::observer::stdio::{StdoutEffect, StdoutObserver};
use vihaco::traits::{GetProgramGlobal, ProgramCounter, StackMemory};
use vihaco::{Effects, Observe, ProgramLoader, Value, composite, observe};
use vihaco_cpu::{CPU, CPUMessage};

/// Re-exported so consumers (e.g. the CLI debugger) can match on step results
/// without depending on `vihaco-cpu` directly.
pub use vihaco_cpu::StepOutcome;

use crate::component::{Circuit, CircuitEffect};
use crate::instruction::CircuitInstruction;
use crate::measurements::{MeasurementEffect, MeasurementObserver, MeasurementResult};
use crate::message::CircuitMessage;

pub const PPVM_MAGIC: u32 = 0x5050564D;

#[derive(Debug, Clone, PartialEq)]
pub struct PPVMDeviceInfo {
    pub magic: u32,
    pub n_qubits: usize,
    pub coefficient_threshold: f64,
}

impl Default for PPVMDeviceInfo {
    fn default() -> Self {
        Self {
            magic: PPVM_MAGIC,
            n_qubits: 0,
            coefficient_threshold: 1e-10,
        }
    }
}

pub type Instruction = PPVMInstruction;

#[composite]
#[derive(Default)]
pub struct PPVM {
    #[program]
    loader: ProgramLoader<PPVMInstruction, PPVMDeviceInfo>,

    #[device(0x00, resolve_with = resolve_cpu)]
    cpu: CPU,

    #[device(0x01, resolve_with = resolve_circuit)]
    circuit: Circuit,

    stdout: StdoutObserver,

    measurement_record: MeasurementObserver,
}

#[derive(Debug, Clone)]
pub enum PPVMEffect {
    Step(StepOutcome),
    Stdout(StdoutEffect),
    Circuit(CircuitEffect),
    Measurement(MeasurementEffect),
}

#[observe(vihaco::observer::stdio::StdoutEffect, effect = PPVMEffect)]
impl PPVM {
    fn observe_stdout_effect(&mut self, effect: &StdoutEffect) -> eyre::Result<Effects<()>> {
        Observe::<StdoutEffect>::observe(&mut self.stdout, effect)
    }
}

#[observe(CircuitEffect, effect = PPVMEffect)]
impl PPVM {
    fn observe_circuit_effect(
        &mut self,
        effect: &CircuitEffect,
    ) -> eyre::Result<Effects<MeasurementEffect>> {
        Observe::<CircuitEffect>::observe(&mut self.circuit, effect)
    }
}

#[observe(MeasurementEffect, effect = PPVMEffect)]
impl PPVM {
    fn observe_measurement_effect(
        &mut self,
        effect: &MeasurementEffect,
    ) -> eyre::Result<Effects<()>> {
        Observe::<MeasurementEffect>::observe(&mut self.measurement_record, effect)
    }
}

impl From<StdoutEffect> for PPVMEffect {
    fn from(value: StdoutEffect) -> Self {
        Self::Stdout(value)
    }
}

impl From<MeasurementEffect> for PPVMEffect {
    fn from(value: MeasurementEffect) -> Self {
        Self::Measurement(value)
    }
}

impl std::fmt::Display for PPVMInstruction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PPVMInstruction::Cpu(inst) => inst.fmt(f),
            PPVMInstruction::Circuit(inst) => inst.fmt(f),
        }
    }
}

impl PartialEq for PPVMInstruction {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (PPVMInstruction::Cpu(a), PPVMInstruction::Cpu(b)) => a == b,
            (PPVMInstruction::Circuit(a), PPVMInstruction::Circuit(b)) => a == b,
            _ => false,
        }
    }
}

impl From<vihaco_cpu::Instruction> for PPVMInstruction {
    fn from(value: vihaco_cpu::Instruction) -> Self {
        Self::Cpu(value)
    }
}

impl From<CircuitInstruction> for PPVMInstruction {
    fn from(value: CircuitInstruction) -> Self {
        Self::Circuit(value)
    }
}

impl PPVM {
    fn resolve_cpu(&mut self, inst: &vihaco_cpu::Instruction) -> eyre::Result<CPUMessage> {
        match inst {
            vihaco_cpu::Instruction::IndirectCall => {
                let function_id: u32 = self.cpu.stack_top()?.get_function_ref()?;
                let function = self.loader.get_function(function_id as usize)?;
                Ok(CPUMessage::FunctionInfo {
                    arity: function.signature.params.len() as u32,
                    start_address: function.start_address,
                })
            }
            vihaco_cpu::Instruction::Print => {
                let value = *self.cpu.stack_top()?;
                match value {
                    vihaco::Value::String(addr) => {
                        let string = self.loader.get_string(addr as usize)?.clone();
                        Ok(CPUMessage::Print(string))
                    }
                    value => Ok(CPUMessage::Print(value.to_string())),
                }
            }
            _ => Ok(CPUMessage::None),
        }
    }

    fn resolve_circuit(&mut self, inst: &CircuitInstruction) -> eyre::Result<CircuitMessage> {
        use crate::instruction::CircuitInstruction::*;
        match inst {
            X | Y | Z | H | S | SAdj | SqrtX | SqrtY | SqrtXAdj | SqrtYAdj | T | TAdj | Measure
            | Reset => {
                let q = self.pop_qubit()?;
                Ok(CircuitMessage::Qubit(q))
            }
            CNOT | CZ => {
                let q1 = self.pop_qubit()?;
                let q0 = self.pop_qubit()?;
                Ok(CircuitMessage::TwoQubit(q0, q1))
            }
            RX | RY | RZ | Depolarize | Loss => {
                let theta = self.pop_f64()?;
                let q = self.pop_qubit()?;
                Ok(CircuitMessage::QubitAndFloat(q, theta))
            }
            RXX | RYY | RZZ | Depolarize2 => {
                let theta = self.pop_f64()?;
                let q0 = self.pop_qubit()?;
                let q1 = self.pop_qubit()?;
                Ok(CircuitMessage::TwoQubitAndFloat(q0, q1, theta))
            }
            U3 => {
                let lam = self.pop_f64()?;
                let phi = self.pop_f64()?;
                let theta = self.pop_f64()?;
                let q = self.pop_qubit()?;
                Ok(CircuitMessage::QubitU3(q, theta, phi, lam))
            }

            // TODO: pop actual float arrays?
            PauliError => {
                let pz = self.pop_f64()?;
                let py = self.pop_f64()?;
                let px = self.pop_f64()?;
                let q = self.pop_qubit()?;
                Ok(CircuitMessage::QubitAndFloatArr3(q, [px, py, pz]))
            }
            CorrelatedLoss => {
                let p2 = self.pop_f64()?;
                let p1 = self.pop_f64()?;
                let p0 = self.pop_f64()?;
                let q0 = self.pop_qubit()?;
                let q1 = self.pop_qubit()?;
                Ok(CircuitMessage::TwoQubitAndFloatArr3(q0, q1, [p0, p1, p2]))
            }
            TwoQubitPauliError => {
                let mut ps = [0.0; 15];
                for p in ps.iter_mut() {
                    *p = self.pop_f64()?;
                }
                let q0 = self.pop_qubit()?;
                let q1 = self.pop_qubit()?;
                Ok(CircuitMessage::TwoQubitAndFloatArr15(q0, q1, ps))
            }
        }
    }

    fn pop_qubit(&mut self) -> eyre::Result<usize> {
        match self.cpu.stack_pop()? {
            vihaco::Value::U32(v) => Ok(v as usize),
            vihaco::Value::U64(v) => usize::try_from(v).map_err(Into::into),
            vihaco::Value::I64(v) => usize::try_from(v).map_err(Into::into),
            v => Err(eyre::eyre!("Expected qubit address, got {:?}", v)),
        }
    }

    fn pop_f64(&mut self) -> eyre::Result<f64> {
        match self.cpu.stack_pop()? {
            vihaco::Value::F64(v) => Ok(v),
            v => Err(eyre::eyre!("Expected f64 argument, got {:?}", v)),
        }
    }

    pub fn init(&mut self) -> eyre::Result<()> {
        self.init_inner(None)
    }

    /// Like [`PPVM::init`], but seed the circuit RNG deterministically so the
    /// run is reproducible.
    pub fn init_with_seed(&mut self, seed: u64) -> eyre::Result<()> {
        self.init_inner(Some(seed))
    }

    fn init_inner(&mut self, seed: Option<u64>) -> eyre::Result<()> {
        let info = &self.loader.module.extra;
        if info.n_qubits == 0 {
            return Err(eyre::eyre!("device circuit.n_qubits must be declared"));
        }
        self.circuit = match seed {
            Some(seed) => Circuit::new_with_seed(info.n_qubits, info.coefficient_threshold, seed),
            None => Circuit::new(info.n_qubits, info.coefficient_threshold),
        };

        // push entry frame
        self.cpu.push_frame(Frame {
            base: 0,
            span: (0, 0, 0),
            function: None,
            ret_pc: 0,
        });

        Ok(())
    }

    pub fn load(
        &mut self,
        module: &vihaco::module::Module<Instruction, vihaco::Value, vihaco::Type, PPVMDeviceInfo>,
    ) -> eyre::Result<()> {
        self.loader.module = module.clone();
        Ok(())
    }

    pub fn step_once(&mut self) -> eyre::Result<StepOutcome> {
        let inst = self.peek_instruction()?.clone();
        let effects = self.execute_effects(inst)?;
        self.continue_effects(effects)
    }

    /// Program counter: the index of the next instruction to execute.
    pub fn current_pc(&self) -> u32 {
        self.loader.pc()
    }

    /// The next instruction to execute, or `None` once execution has run off
    /// the end of the code. Intended for debuggers/inspection.
    pub fn current_instruction(&self) -> Option<PPVMInstruction> {
        self.peek_instruction().ok().cloned()
    }

    fn execute_effects(&mut self, inst: Instruction) -> eyre::Result<Effects<PPVMEffect>> {
        log::debug!("exec inst: {:?}, stack: {:?}", inst, self.cpu.stack());
        match inst {
            PPVMInstruction::Cpu(cpu_inst) => {
                let msg = self.resolve_cpu(&cpu_inst)?;
                self.cpu.set_current_pc(self.loader.pc());
                let stdout_effect = match (&cpu_inst, &msg) {
                    (vihaco_cpu::Instruction::Print, vihaco_cpu::CPUMessage::Print(text)) => {
                        Some(PPVMEffect::Stdout(StdoutEffect(text.clone())))
                    }
                    _ => None,
                };
                let outcome = vihaco::expect_exactly_one_effect(
                    vihaco::GeneratedComponent::execute_generated(&mut self.cpu, cpu_inst, msg)?,
                )?;
                // Advance past a breakpoint as well, so the debugger that paused
                // on it doesn't re-hit the same instruction on the next step.
                if matches!(outcome, StepOutcome::Continue | StepOutcome::Breakpoint) {
                    if let Some(target) = self.cpu.take_pending_pc() {
                        *self.loader.pc_mut() = target;
                    } else {
                        *self.loader.pc_mut() += 1;
                    }
                }
                let mut effects = Effects::one(PPVMEffect::Step(outcome));
                if let Some(stdout_effect) = stdout_effect {
                    effects = effects.append(stdout_effect);
                }
                Ok(effects)
            }
            PPVMInstruction::Circuit(inst) => {
                let msg = self.resolve_circuit(&inst)?;
                let circuit_effects = <Circuit as vihaco::GeneratedComponent>::execute_generated(
                    &mut self.circuit,
                    inst,
                    msg,
                )?;
                *self.loader.pc_mut() += 1;
                let mut effects = Effects::one(PPVMEffect::Step(StepOutcome::Continue));
                for measurement_effect in circuit_effects {
                    effects = effects.append(PPVMEffect::Measurement(measurement_effect));
                }
                Ok(effects)
            }
        }
    }

    fn continue_effects(&mut self, effects: Effects<PPVMEffect>) -> eyre::Result<StepOutcome> {
        let mut step_outcome = None;
        for effect in effects {
            match effect {
                PPVMEffect::Step(outcome) => {
                    if step_outcome.replace(outcome).is_some() {
                        return Err(eyre::eyre!(
                            "expected exactly one PPVM step effect, got multiple"
                        ));
                    }
                }
                effect => self.continue_observer_effect(effect)?,
            }
        }

        step_outcome.ok_or_else(|| eyre::eyre!("expected exactly one PPVM step effect, got 0"))
    }

    fn continue_observer_effect(&mut self, effect: PPVMEffect) -> eyre::Result<()> {
        match effect {
            PPVMEffect::Stdout(effect) => {
                let follow_ups = Observe::<StdoutEffect>::observe(self, &effect)?;
                self.continue_observer_effects(follow_ups)
            }
            PPVMEffect::Circuit(effect) => {
                let follow_ups = Observe::<CircuitEffect>::observe(self, &effect)?;
                self.continue_observer_effects(follow_ups)
            }
            PPVMEffect::Measurement(effect) => {
                let follow_ups = Observe::<MeasurementEffect>::observe(self, &effect)?;
                // NOTE: push measurements to stack; two booleans: outcome, is_lost
                for outcome in effect.measurement_results {
                    let m = Value::U32(outcome as u32);
                    self.cpu.stack_push(m);
                }
                self.continue_observer_effects(follow_ups)
            }
            PPVMEffect::Step(_) => Err(eyre::eyre!(
                "unexpected Step effect while continuing PPVM observer follow-ups"
            )),
        }
    }

    fn continue_observer_effects(&mut self, effects: Effects<PPVMEffect>) -> eyre::Result<()> {
        for effect in effects {
            self.continue_observer_effect(effect)?;
        }
        Ok(())
    }

    pub fn run(&mut self) -> eyre::Result<StepOutcome> {
        self.run_with_seed(None)
    }

    /// Like [`PPVM::run`], but seed the circuit RNG deterministically when
    /// `seed` is `Some`, making the run reproducible.
    pub fn run_with_seed(&mut self, seed: Option<u64>) -> eyre::Result<StepOutcome> {
        match seed {
            Some(seed) => self.init_with_seed(seed)?,
            None => self.init()?,
        }

        loop {
            // Breakpoints only pause the interactive debugger; a batch run
            // skips straight past them.
            match self.step_once()? {
                StepOutcome::Continue | StepOutcome::Breakpoint => continue,
                action => return Ok(action),
            }
        }
    }

    pub fn stdout(&self) -> &[u8] {
        self.stdout.output()
    }

    pub fn measurement_record(&self) -> Vec<MeasurementResult> {
        self.measurement_record.record.clone()
    }

    pub fn load_program(&mut self, program: &str) -> eyre::Result<()> {
        let module = crate::compile_program(program)?;
        self.load(&module)?;
        Ok(())
    }

    /// Load from a file, auto-detecting the format: if it starts with the PPVM
    /// magic it is loaded as `.ssb` bytecode, otherwise it is parsed as `.sst`
    /// source text. A magic match commits to the bytecode path — a corrupt
    /// `.ssb` errors rather than silently falling back to the text parser.
    pub fn load_file(&mut self, path: &str) -> eyre::Result<()> {
        let bytes = std::fs::read(path)?;
        if crate::bytecode::is_bytecode(&bytes) {
            self.load_bytecode(&bytes)
        } else {
            self.load_program(std::str::from_utf8(&bytes)?)
        }
    }

    /// Load a module from an in-memory `.ssb` byte stream.
    pub fn load_bytecode(&mut self, bytes: &[u8]) -> eyre::Result<()> {
        let module = crate::bytecode::module_from_bytes(bytes)?;
        self.load(&module)
    }

    /// Read a `.ssb` file and load the module it contains.
    pub fn load_bytecode_file(&mut self, path: &str) -> eyre::Result<()> {
        let bytes = std::fs::read(path)?;
        self.load_bytecode(&bytes)
    }

    pub fn run_program(&mut self, program: &str) -> eyre::Result<()> {
        self.load_program(program)?;
        self.run()?;
        Ok(())
    }

    pub fn run_file(&mut self, path: &str) -> eyre::Result<()> {
        self.load_file(path)?;
        self.run()?;
        Ok(())
    }
}

impl vihaco::Reset for PPVM {
    fn reset(&mut self) {
        self.cpu.reset();
        self.circuit.reset();
        self.loader.pc = 0;
        self.measurement_record.record.clear();
    }
}

#[cfg(test)]
mod tests {
    use vihaco::{Type, Value, module::Module};

    use super::*;

    #[test]
    fn test_run_ppvm() -> eyre::Result<()> {
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();

        module.extra.n_qubits = 2;

        /*
        const.u64 0
        gate h
        */
        let zero = PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(Value::U64(0)));
        let one = PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(Value::U64(1)));
        module.code.push(zero.clone());
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::H));

        /*
        const.u64 0
        gate t
        */

        module.code.push(zero.clone());
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::T));

        /*
        const.u64 0
        const.u64 1
        gate cnot
        */
        module.code.push(zero.clone());
        module.code.push(one.clone());
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::CNOT));

        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;

        for _ in 0..module.code.len() {
            machine.step_once()?;
            assert!(machine.cpu.stack().len() <= 2);
        }

        let num_coefficients = match &machine.circuit {
            Circuit::Bits64(ex) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::Bits128(ex) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::Bits256(ex) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::Bits512(ex) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::Bits1024(ex) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::Bits2048(ex) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
        };

        assert_eq!(num_coefficients, 2);
        Ok(())
    }

    #[test]
    fn test_device_decl() -> eyre::Result<()> {
        // Equivalent .sst source:
        //
        //     device circuit.n_qubits 5;
        //     device circuit.coefficient_threshold 1e-10;
        //
        //     fn @main() { ...5-qubit GHZ + 5 measurements... }

        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();

        module.extra.n_qubits = 5;
        module.extra.coefficient_threshold = 1e-10;

        // 5-qubit GHZ: H on q0, then CNOT(q_i, q_{i+1}) for i = 0..4.
        /*
        const.u64 0
        gate h
        */
        module
            .code
            .push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                Value::U64(0),
            )));
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::H));

        for i in 0..4u64 {
            /*
            const.u64 i
            const.u64 i+1
            gate cnot
            */
            module
                .code
                .push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                    Value::U64(i),
                )));
            module
                .code
                .push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                    Value::U64(i + 1),
                )));
            module
                .code
                .push(PPVMInstruction::Circuit(CircuitInstruction::CNOT));
        }

        // Measure all 5 qubits.
        for q in 0..5u64 {
            /*
            const.u64 q
            gate measure
            */
            module
                .code
                .push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                    Value::U64(q),
                )));
            module
                .code
                .push(PPVMInstruction::Circuit(CircuitInstruction::Measure));
        }

        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;

        for _ in 0..module.code.len() {
            machine.step_once()?;
        }

        assert_eq!(machine.measurement_record().len(), 5);
        Ok(())
    }

    // ─── Parser-driven entry points ───────────────────────────────────────

    #[test]
    fn load_program_populates_device_info_and_code() -> eyre::Result<()> {
        let source = "device circuit.n_qubits 2;\n\
                      device circuit.coefficient_threshold 1e-8;\n\
                      fn @main() {\n\
                          const.u64 0\n\
                          gate h\n\
                          ret\n\
                      }\n";
        let mut machine = PPVM::default();
        machine.load_program(source)?;
        assert_eq!(machine.loader.module.extra.n_qubits, 2);
        assert_eq!(machine.loader.module.extra.coefficient_threshold, 1e-8);
        // const.u64 0 / gate h / ret = 3
        assert_eq!(machine.loader.module.code.len(), 3);
        Ok(())
    }

    #[test]
    fn run_program_executes_bell_circuit() -> eyre::Result<()> {
        let source = "device circuit.n_qubits 2;\n\
                      fn @main() {\n\
                          const.u64 0\n\
                          gate h\n\
                          const.u64 0\n\
                          const.u64 1\n\
                          gate cnot\n\
                          const.u64 0\n\
                          gate measure\n\
                          const.u64 1\n\
                          gate measure\n\
                          ret\n\
                      }\n";
        let mut machine = PPVM::default();
        machine.run_program(source)?;
        let record = machine.measurement_record();
        assert_eq!(record.len(), 2);
        Ok(())
    }

    #[test]
    fn reset_clears_the_measurement_record() -> eyre::Result<()> {
        use vihaco::Reset;

        let source = "device circuit.n_qubits 2;\n\
                      fn @main() {\n\
                          const.u64 0\n\
                          gate h\n\
                          const.u64 0\n\
                          const.u64 1\n\
                          gate cnot\n\
                          const.u64 0\n\
                          gate measure\n\
                          const.u64 1\n\
                          gate measure\n\
                          ret\n\
                      }\n";
        let mut machine = PPVM::default();
        machine.run_program(source)?;
        assert_eq!(machine.measurement_record().len(), 2);

        // Resetting the machine must discard the recorded measurements, so a
        // subsequent run does not see stale results leaking in from before.
        machine.reset();
        assert!(
            machine.measurement_record().is_empty(),
            "reset must clear the measurement record"
        );

        Ok(())
    }

    #[test]
    fn init_fails_when_n_qubits_undeclared() -> eyre::Result<()> {
        let source = "fn @main() { ret }\n";
        let mut machine = PPVM::default();
        machine.load_program(source)?;
        let err = machine.init().unwrap_err();
        assert!(err.to_string().contains("circuit.n_qubits"), "err: {err}");
        Ok(())
    }

    #[test]
    fn run_program_reports_parse_errors() {
        let source = "device circuit.n_qubits 2;\n\
                      fn @main() {\n\
                          gate not_a_real_gate\n\
                          ret\n\
                      }\n";
        let mut machine = PPVM::default();
        let err = machine.run_program(source).unwrap_err();
        assert!(
            err.to_string().contains("parsing failed")
                || err.to_string().contains("unhandled raw form"),
            "err: {err}"
        );
    }

    // ─── Breakpoints ──────────────────────────────────────────────────────

    /// Bell circuit with a `breakpoint` between the two measurements.
    const BREAKPOINT_PROGRAM: &str = "device circuit.n_qubits 2;\n\
                                      fn @main() {\n\
                                          const.u64 0\n\
                                          gate h\n\
                                          const.u64 0\n\
                                          const.u64 1\n\
                                          gate cnot\n\
                                          const.u64 0\n\
                                          gate measure\n\
                                          breakpoint\n\
                                          const.u64 1\n\
                                          gate measure\n\
                                          ret\n\
                                      }\n";

    #[test]
    fn run_ignores_breakpoints() -> eyre::Result<()> {
        // A batch run must execute straight through the breakpoint and record
        // both measurements, exactly as if it weren't there.
        let mut machine = PPVM::default();
        machine.run_program(BREAKPOINT_PROGRAM)?;
        assert_eq!(machine.measurement_record().len(), 2);
        Ok(())
    }

    #[test]
    fn step_once_advances_past_breakpoint() -> eyre::Result<()> {
        let mut machine = PPVM::default();
        machine.load_program(BREAKPOINT_PROGRAM)?;
        machine.init()?;

        // Step until the breakpoint pauses us.
        let mut outcome = StepOutcome::Continue;
        for _ in 0..machine.loader.module.code.len() {
            outcome = machine.step_once()?;
            if outcome == StepOutcome::Breakpoint {
                break;
            }
        }
        assert_eq!(outcome, StepOutcome::Breakpoint, "breakpoint should pause");
        let pc_at_break = machine.current_pc();

        // Stepping again must make progress (advance the pc) rather than
        // re-hitting the same breakpoint instruction.
        let next = machine.step_once()?;
        assert_ne!(
            next,
            StepOutcome::Breakpoint,
            "must move past the breakpoint"
        );
        assert!(machine.current_pc() > pc_at_break, "pc must advance");

        Ok(())
    }
}
