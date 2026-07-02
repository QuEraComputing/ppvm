// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use eyre::{Result, eyre};
use vihaco::frame::Frame;
use vihaco::machine::StackFrame;
use vihaco::observer::stdio::{StdoutEffect, StdoutObserver};
use vihaco::traits::{GetProgramGlobal, ProgramCounter, StackMemory};
use vihaco::{Effects, Observe, ProgramLoader, Value, composite, observe};
use vihaco_cpu::{CPU, CPUMessage};

/// Re-exported so consumers (e.g. the CLI debugger) can match on step results
/// without depending on `vihaco-cpu` directly.
pub use vihaco_cpu::StepOutcome;

use crate::component::Circuit;
#[cfg(test)]
use crate::component::TableauCircuit;
use crate::measurements::{
    CircuitOutcomeEffect, MeasurementEffect, MeasurementObserver, MeasurementResult, TraceEffect,
    TraceObserver,
};
use vihaco_circuit_isa::{CircuitEffect, CircuitInstruction, CircuitMessage};

pub const PPVM_MAGIC: u32 = 0x5050564D;

/// Which execution backend the circuit runs on. Selected via the
/// `device circuit.backend` header; defaults to `Tableau` so existing
/// programs that don't declare a backend keep working.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, vihaco_parser::Parse)]
pub enum BackendKind {
    #[default]
    Tableau,
    PauliSum,
    #[token = "lossy_paulisum"]
    LossyPauliSum,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PPVMDeviceInfo {
    pub magic: u32,
    pub n_qubits: usize,
    pub coefficient_threshold: f64,
    pub backend: BackendKind,
    pub observable: Option<String>,
    pub max_pauli_weight: Option<usize>,
}

impl Default for PPVMDeviceInfo {
    fn default() -> Self {
        Self {
            magic: PPVM_MAGIC,
            n_qubits: 0,
            coefficient_threshold: 1e-10,
            backend: BackendKind::default(),
            observable: None,
            max_pauli_weight: None,
        }
    }
}

pub type Instruction = PPVMInstruction;

#[composite]
#[derive(Default)]
pub struct PPVM {
    #[program]
    loader: ProgramLoader<PPVMInstruction, PPVMDeviceInfo>,

    #[device(0x00)]
    cpu: CPU,

    #[device(0x01)]
    circuit: Circuit,

    stdout: StdoutObserver,

    measurement_record: MeasurementObserver,

    trace_record: TraceObserver,
}

#[derive(Debug, Clone)]
pub enum PPVMEffect {
    Step(StepOutcome),
    Stdout(StdoutEffect),
    Circuit(Box<CircuitEffect>),
    Measurement(MeasurementEffect),
    Trace(TraceEffect),
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
    ) -> eyre::Result<Effects<CircuitOutcomeEffect>> {
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

#[observe(TraceEffect, effect = PPVMEffect)]
impl PPVM {
    fn observe_trace_effect(&mut self, effect: &TraceEffect) -> eyre::Result<Effects<()>> {
        Observe::<TraceEffect>::observe(&mut self.trace_record, effect)
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

impl From<TraceEffect> for PPVMEffect {
    fn from(value: TraceEffect) -> Self {
        Self::Trace(value)
    }
}

impl From<CircuitOutcomeEffect> for PPVMEffect {
    fn from(value: CircuitOutcomeEffect) -> Self {
        match value {
            CircuitOutcomeEffect::Measurement(m) => Self::Measurement(m),
            CircuitOutcomeEffect::Trace(t) => Self::Trace(t),
        }
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
        use CircuitInstruction::*;
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
                let q1 = self.pop_qubit()?;
                let q0 = self.pop_qubit()?;
                Ok(CircuitMessage::TwoQubitAndFloat(q0, q1, theta))
            }
            R => {
                let theta = self.pop_f64()?;
                let axis_angle = self.pop_f64()?;
                let q = self.pop_qubit()?;
                Ok(CircuitMessage::QubitAndTwoFloats(q, axis_angle, theta))
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
                let q1 = self.pop_qubit()?;
                let q0 = self.pop_qubit()?;
                Ok(CircuitMessage::TwoQubitAndFloatArr3(q0, q1, [p0, p1, p2]))
            }
            TwoQubitPauliError => {
                let mut ps = [0.0; 15];
                for p in ps.iter_mut().rev() {
                    *p = self.pop_f64()?;
                }
                let q1 = self.pop_qubit()?;
                let q0 = self.pop_qubit()?;
                Ok(CircuitMessage::TwoQubitAndFloatArr15(q0, q1, ps))
            }
            Trace => {
                let s = self.pop_string()?;
                Ok(CircuitMessage::PauliPatternStr(s))
            }
            Truncate => Ok(CircuitMessage::None),
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

    /// Pop a `Value::String(addr)`, look up the addr in the module's string
    /// table, and return the owned string. Used by `Trace` to resolve its
    /// Pauli-pattern operand before the executor sees the message.
    fn pop_string(&mut self) -> eyre::Result<String> {
        match self.cpu.stack_pop()? {
            vihaco::Value::String(addr) => Ok(self.loader.get_string(addr as usize)?.clone()),
            v => Err(eyre::eyre!("Expected string operand, got {:?}", v)),
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
        self.circuit = match (info.backend, seed) {
            (BackendKind::Tableau, None) => Circuit::tableau(info),
            (BackendKind::Tableau, Some(seed)) => Circuit::tableau_with_seed(info, seed),
            (BackendKind::PauliSum, _) => {
                let terms = parse_observable_terms(info)?;
                Circuit::paulisum(info, &terms)
            }
            (BackendKind::LossyPauliSum, _) => {
                let terms = parse_observable_terms(info)?;
                Circuit::lossy_paulisum(info, &terms)
            }
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

    /// Append one REPL command's lowered VM ops and run just that block against
    /// the persistent state, advancing the pc through it. NOTE: an out-of-range
    /// qubit index *panics* in the tableau rather than erroring, so callers must
    /// bounds-check qubit operands against `n_qubits` first.
    pub fn execute_single_instruction(&mut self, instrs: &[PPVMInstruction]) -> eyre::Result<()> {
        let start = self.loader.module.code.len() as u32;
        self.loader.module.code.extend_from_slice(instrs);
        *self.loader.pc_mut() = start;
        for _ in 0..instrs.len() {
            self.step_once()?;
        }
        Ok(())
    }

    /// Render the current circuit state (tableau / Pauli sum) for the REPL's
    /// `show` command. Delegates to the circuit's size-specific tableau.
    pub fn state_string(&self) -> String {
        self.circuit.state_string()
    }

    /// Build a fresh, initialized `n_qubits`-qubit device with no code. The
    /// REPL's `device` command uses this to (re)create the machine. Errors if
    /// `n_qubits` is zero (a device must have at least one qubit).
    pub fn with_qubits(n_qubits: usize) -> eyre::Result<Self> {
        let mut machine = Self::default();
        let mut module = vihaco::module::Module::<
            PPVMInstruction,
            Value,
            vihaco::Type,
            PPVMDeviceInfo,
        >::default();
        module.extra.n_qubits = n_qubits;
        machine.load(&module)?;
        machine.init()?;
        Ok(machine)
    }

    /// Lower a single circuit instruction — qubit operands first, then float
    /// params, per the push-in-order / pop-in-reverse convention — and execute
    /// it against the persistent state. Qubit indices are bounds-checked against
    /// the device size first, because the tableau panics (rather than erroring)
    /// on an out-of-range qubit.
    pub fn apply_circuit_instruction(
        &mut self,
        inst: CircuitInstruction,
        qubits: &[usize],
        params: &[f64],
    ) -> eyre::Result<()> {
        let n_qubits = self.loader.module.extra.n_qubits;
        for &q in qubits {
            if q >= n_qubits {
                eyre::bail!("qubit {q} out of range for {n_qubits}-qubit device");
            }
        }

        let mut instrs = Vec::with_capacity(qubits.len() + params.len() + 1);
        for &q in qubits {
            instrs.push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                Value::U64(q as u64),
            )));
        }
        for &p in params {
            instrs.push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                Value::F64(p),
            )));
        }
        instrs.push(PPVMInstruction::Circuit(inst));

        // Run the appended block, then roll back pc, code, and operand stack (on
        // every path, including error) so the injected op is transparent to a
        // paused program — only its tableau/measurement effects persist.
        let saved_pc = self.loader.pc();
        let saved_len = self.loader.module.code.len();
        let saved_stack = self.cpu.stack_len();
        let result = self.execute_single_instruction(&instrs);
        self.loader.module.code.truncate(saved_len);
        *self.loader.pc_mut() = saved_pc;
        self.cpu.stack_mut().truncate(saved_stack);
        result
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
                for outcome in circuit_effects {
                    effects = effects.append(PPVMEffect::from(outcome));
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
            PPVMEffect::Trace(effect) => {
                let value = effect.value;
                let follow_ups = Observe::<TraceEffect>::observe(self, &effect)?;
                // Mirror the measurement wiring: append to the trace record
                // (via the observer above) AND push the value onto the CPU
                // stack so user bytecode can consume it. Plan Task 7.
                self.cpu.stack_push(Value::F64(value));
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

    /// Per-trace values collected by `Trace` instructions during the run.
    /// Parallel to [`PPVM::measurement_record`]: one f64 per `Trace` executed,
    /// in execution order.
    pub fn trace_record(&self) -> Vec<f64> {
        self.trace_record.record.clone()
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
        self.trace_record.record.clear();
    }
}

/// Parse the `device circuit.observable` header into Pauli-sum terms ready to
/// seed a `PauliSum` / `LossyPauliSum` state. Single-Pauli observables from
/// Phase 2 keep working as the degenerate one-term case; multi-term sums like
/// `"1.0*ZZ + 0.5*XX"` are handled by [`parse_pauli_sum_terms`].
fn parse_observable_terms(info: &PPVMDeviceInfo) -> Result<Vec<(String, f64)>> {
    let observable = info.observable.as_deref().ok_or_else(|| {
        eyre!(
            "the {:?} backend requires `device circuit.observable` to be set",
            info.backend
        )
    })?;
    crate::observable::parse_pauli_sum_terms(observable, info.n_qubits)
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
        circuit.h
        */
        let zero = PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(Value::U64(0)));
        let one = PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(Value::U64(1)));
        module.code.push(zero.clone());
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::H));

        /*
        const.u64 0
        circuit.t
        */

        module.code.push(zero.clone());
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::T));

        /*
        const.u64 0
        const.u64 1
        circuit.cnot
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
            Circuit::Tableau(TableauCircuit::Bits64(ex)) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::Tableau(TableauCircuit::Bits128(ex)) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::Tableau(TableauCircuit::Bits256(ex)) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::Tableau(TableauCircuit::Bits512(ex)) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::Tableau(TableauCircuit::Bits1024(ex)) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::Tableau(TableauCircuit::Bits2048(ex)) => {
                println!("{}", ex.tab);
                ex.tab.coefficients.len()
            }
            Circuit::PauliSum(_) | Circuit::LossyPauliSum(_) => {
                panic!("test expects the default Tableau backend, got a PauliSum variant");
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
        circuit.h
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
            circuit.cnot
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
            circuit.measure
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

    // ─── Incremental execution (REPL) ─────────────────────────────────────

    #[test]
    fn execute_single_instruction_persists_state_across_calls() -> eyre::Result<()> {
        use crate::measurements::MeasurementOutcome;

        // A 1-qubit device with no code; the REPL builds up instructions
        // incrementally, one command at a time, rather than loading a program.
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 1;

        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;

        // First command: X on q0 (|0> -> |1>).
        let x = [
            PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(Value::U64(0))),
            PPVMInstruction::Circuit(CircuitInstruction::X),
        ];
        machine.execute_single_instruction(&x)?;
        // No measurement yet.
        assert!(machine.measurement_record().is_empty());

        // Second command: measure q0. The X from the first command must persist,
        // so the outcome is deterministically |1>.
        let measure = [
            PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(Value::U64(0))),
            PPVMInstruction::Circuit(CircuitInstruction::Measure),
        ];
        machine.execute_single_instruction(&measure)?;

        let record = machine.measurement_record();
        assert_eq!(record.len(), 1);
        assert_eq!(record[0].as_slice(), [MeasurementOutcome::One]);
        Ok(())
    }

    #[test]
    fn execute_single_instruction_propagates_engine_errors() -> eyre::Result<()> {
        // The REPL relies on engine errors surfacing as `Err` (so it can print
        // them and keep looping) rather than panicking. A circuit with no qubit
        // operand on the stack is one such propagating error.
        //
        // NOTE: an out-of-range qubit index (>= n_qubits) currently *panics* in
        // the tableau rather than erroring, so the REPL command layer must
        // bounds-check qubit indices before calling `execute`.
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 1;

        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;

        // `circuit.h` with nothing on the stack: `pop_qubit` fails.
        let missing_operand = [PPVMInstruction::Circuit(CircuitInstruction::H)];
        assert!(
            machine
                .execute_single_instruction(&missing_operand)
                .is_err()
        );
        Ok(())
    }

    #[test]
    fn state_string_renders_a_small_device() -> eyre::Result<()> {
        let source = "device circuit.n_qubits 2;\nfn @main() { ret }\n";
        let mut machine = PPVM::default();
        machine.load_program(source)?;
        machine.init()?;

        let rendered = machine.state_string();
        assert!(
            !rendered.is_empty(),
            "state_string should render the tableau"
        );
        // PPVM delegates to the circuit.
        assert_eq!(rendered, machine.circuit.state_string());
        Ok(())
    }

    #[test]
    fn resolve_circuit_pops_operands_in_reverse_of_push_order() -> eyre::Result<()> {
        // Convention: operands are pushed in argument order (q0, q1, then any
        // floats) and popped in reverse. So every two-qubit circuit must read q0 as
        // the first operand pushed, consistently, with or without trailing
        // floats. (CNOT already obeyed this; the float-carrying arms did not.)
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 8;
        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;

        // CNOT: push q0=2, q1=5.
        machine.cpu.stack_push(Value::U32(2));
        machine.cpu.stack_push(Value::U32(5));
        assert_eq!(
            machine.resolve_circuit(&CircuitInstruction::CNOT)?,
            CircuitMessage::TwoQubit(2, 5)
        );

        // RXX: push q0=2, q1=5, theta — same qubit order as CNOT.
        machine.cpu.stack_push(Value::U32(2));
        machine.cpu.stack_push(Value::U32(5));
        machine.cpu.stack_push(Value::F64(0.3));
        assert_eq!(
            machine.resolve_circuit(&CircuitInstruction::RXX)?,
            CircuitMessage::TwoQubitAndFloat(2, 5, 0.3)
        );

        // CorrelatedLoss: push q0=2, q1=5, p0, p1, p2.
        machine.cpu.stack_push(Value::U32(2));
        machine.cpu.stack_push(Value::U32(5));
        machine.cpu.stack_push(Value::F64(0.1));
        machine.cpu.stack_push(Value::F64(0.2));
        machine.cpu.stack_push(Value::F64(0.3));
        assert_eq!(
            machine.resolve_circuit(&CircuitInstruction::CorrelatedLoss)?,
            CircuitMessage::TwoQubitAndFloatArr3(2, 5, [0.1, 0.2, 0.3])
        );
        Ok(())
    }

    #[test]
    fn with_qubits_builds_an_initialized_device() -> eyre::Result<()> {
        use crate::measurements::MeasurementOutcome;

        let mut m = PPVM::with_qubits(2)?;
        // The device is ready to take instructions immediately.
        m.apply_circuit_instruction(CircuitInstruction::X, &[0], &[])?;
        m.apply_circuit_instruction(CircuitInstruction::Measure, &[0], &[])?;
        let record = m.measurement_record();
        assert_eq!(record.len(), 1);
        assert_eq!(record[0].as_slice(), [MeasurementOutcome::One]);
        Ok(())
    }

    #[test]
    fn with_qubits_zero_is_an_error() {
        assert!(PPVM::with_qubits(0).is_err());
    }

    #[test]
    fn apply_circuit_instruction_bounds_checks_qubits() -> eyre::Result<()> {
        // q1 is out of range on a 1-qubit device: this must error rather than
        // panic in the tableau.
        let mut m = PPVM::with_qubits(1)?;
        let err = m
            .apply_circuit_instruction(CircuitInstruction::X, &[1], &[])
            .unwrap_err();
        assert!(err.to_string().contains("out of range"), "got: {err}");
        Ok(())
    }

    // ─── Parser-driven entry points ───────────────────────────────────────

    #[test]
    fn load_program_populates_device_info_and_code() -> eyre::Result<()> {
        let source = "device circuit.n_qubits 2;\n\
                      device circuit.coefficient_threshold 1e-8;\n\
                      fn @main() {\n\
                          const.u64 0\n\
                          circuit.h\n\
                          ret\n\
                      }\n";
        let mut machine = PPVM::default();
        machine.load_program(source)?;
        assert_eq!(machine.loader.module.extra.n_qubits, 2);
        assert_eq!(machine.loader.module.extra.coefficient_threshold, 1e-8);
        // const.u64 0 / circuit.h / ret = 3
        assert_eq!(machine.loader.module.code.len(), 3);
        Ok(())
    }

    #[test]
    fn run_program_executes_bell_circuit() -> eyre::Result<()> {
        let source = "device circuit.n_qubits 2;\n\
                      fn @main() {\n\
                          const.u64 0\n\
                          circuit.h\n\
                          const.u64 0\n\
                          const.u64 1\n\
                          circuit.cnot\n\
                          const.u64 0\n\
                          circuit.measure\n\
                          const.u64 1\n\
                          circuit.measure\n\
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
                          circuit.h\n\
                          const.u64 0\n\
                          const.u64 1\n\
                          circuit.cnot\n\
                          const.u64 0\n\
                          circuit.measure\n\
                          const.u64 1\n\
                          circuit.measure\n\
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
                          circuit.not_a_real_gate\n\
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
                                          circuit.h\n\
                                          const.u64 0\n\
                                          const.u64 1\n\
                                          circuit.cnot\n\
                                          const.u64 0\n\
                                          circuit.measure\n\
                                          breakpoint\n\
                                          const.u64 1\n\
                                          circuit.measure\n\
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

    #[test]
    fn paulisum_truncate_runs_without_error() -> eyre::Result<()> {
        // Smoke test: a `circuit.truncate` reaches the PauliSum executor's
        // Truncate arm and calls `state.truncate()`. Task 8 makes the
        // observable mandatory for PauliSum init, so seed `Z` here.
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 1;
        module.extra.backend = BackendKind::PauliSum;
        module.extra.observable = Some("Z".to_string());
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::Truncate));

        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;
        machine.step_once()?;
        Ok(())
    }

    #[test]
    fn paulisum_trace_populates_trace_record() -> eyre::Result<()> {
        // End-to-end Trace pipeline: with the observable `Z` seeded (Task 8),
        // tracing the `Z0` pattern picks up that one term with coefficient
        // 1.0, so the trace should be exactly 1.0.
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 1;
        module.extra.backend = BackendKind::PauliSum;
        module.extra.observable = Some("Z".to_string());
        // `Z0` matches a Z on qubit 0; the parser requires position anchors.
        module.strings.push("Z0".to_string());

        module
            .code
            .push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                Value::String(0),
            )));
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::Trace));

        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;
        for _ in 0..module.code.len() {
            machine.step_once()?;
        }

        assert_eq!(machine.trace_record(), vec![1.0]);
        Ok(())
    }

    #[test]
    fn paulisum_multi_term_observable_seeds_all_terms() -> eyre::Result<()> {
        // Task 11: a sum-valued observable seeds every term. With
        // `"ZZ + 0.5*XX"` the state holds `1.0 * ZZ + 0.5 * XX`; tracing
        // `[XZ]0[XZ]1` matches both words and returns 1.0 + 0.5 = 1.5.
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 2;
        module.extra.backend = BackendKind::PauliSum;
        module.extra.observable = Some("ZZ + 0.5*XX".to_string());
        module.strings.push("[XZ]0[XZ]1".to_string());

        module
            .code
            .push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                Value::String(0),
            )));
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::Trace));

        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;
        for _ in 0..module.code.len() {
            machine.step_once()?;
        }

        assert_eq!(machine.trace_record(), vec![1.5]);
        Ok(())
    }

    #[test]
    fn paulisum_init_rejects_missing_observable() {
        // Task 8 requires `device circuit.observable` for PauliSum / Lossy.
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 1;
        module.extra.backend = BackendKind::PauliSum;

        let mut machine = PPVM::default();
        machine.load(&module).unwrap();
        let err = machine.init().unwrap_err();
        assert!(
            err.to_string().contains("observable"),
            "expected observable-related error, got: {err}"
        );
    }

    #[test]
    fn paulisum_init_rejects_mismatched_observable_length() {
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 2;
        module.extra.backend = BackendKind::PauliSum;
        // Three letters but only two qubits — should error.
        module.extra.observable = Some("ZZZ".to_string());

        let mut machine = PPVM::default();
        machine.load(&module).unwrap();
        let err = machine.init().unwrap_err();
        assert!(
            err.to_string().contains("invalid Pauli-sum"),
            "expected parser rejection, got: {err}"
        );
    }

    #[test]
    fn tableau_truncate_is_silent_no_op() -> eyre::Result<()> {
        // Task 9: `circuit.truncate` on the default Tableau backend should run
        // without error — the tableau prunes via coefficient_threshold during
        // every gate, so the explicit Truncate instruction has nothing to do.
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 1;
        // backend defaults to Tableau; no observable needed.
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::Truncate));

        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;
        machine.step_once()?;
        // No observer effects emitted by Truncate.
        assert!(machine.measurement_record().is_empty());
        assert!(machine.trace_record().is_empty());
        Ok(())
    }

    #[test]
    fn apply_circuit_instruction_preserves_pc_and_code_len() -> eyre::Result<()> {
        use crate::measurements::MeasurementOutcome;

        // breakpoint; then measure q0. Step to the breakpoint, inject X, resume.
        let src = "device circuit.n_qubits 1;\n\
                   fn @main() { breakpoint\n const.u64 0\n circuit.measure\n ret }\n";
        let mut m = PPVM::default();
        m.load_program(src)?;
        m.init()?;

        // Run until the breakpoint pauses us.
        loop {
            if m.step_once()? == StepOutcome::Breakpoint {
                break;
            }
        }
        let pc = m.current_pc();
        let len = m.loader.module.code.len();
        let depth = m.cpu.stack().len();

        // Inject X on q0 while "paused".
        m.apply_circuit_instruction(CircuitInstruction::X, &[0], &[])?;

        // The debugger's position must be untouched.
        assert_eq!(m.current_pc(), pc, "pc must be preserved");
        assert_eq!(
            m.loader.module.code.len(),
            len,
            "appended op must be truncated back"
        );
        assert_eq!(
            m.cpu.stack().len(),
            depth,
            "operand stack must be left unchanged"
        );

        // And the X took effect: resuming the program measures |1>.
        while !matches!(m.step_once()?, StepOutcome::Return | StepOutcome::Halt) {}
        let rec = m.measurement_record();
        assert_eq!(rec.len(), 1);
        assert_eq!(rec[0].as_slice(), [MeasurementOutcome::One]);
        Ok(())
    }

    #[test]
    fn apply_circuit_instruction_preserves_pc_and_code_on_error() -> eyre::Result<()> {
        // A program stepped to a known pc; an injected gate that errors mid-block
        // must NOT corrupt the code vector or the program counter.
        let src = "device circuit.n_qubits 1;\n\
                   fn @main() { breakpoint\n const.u64 0\n circuit.measure\n ret }\n";
        let mut m = PPVM::default();
        m.load_program(src)?;
        m.init()?;
        loop {
            if m.step_once()? == StepOutcome::Breakpoint {
                break;
            }
        }
        let pc = m.current_pc();
        let len = m.loader.module.code.len();
        let depth = m.cpu.stack().len();

        // RX with no float param errors during execution: resolve_circuit pops
        // the f64 first and finds the qubit value (U64) instead, returning Err.
        let err = m.apply_circuit_instruction(CircuitInstruction::RX, &[0], &[]);
        assert!(err.is_err(), "expected the injected gate to error");

        // The failed injection must leave the debugger untouched.
        assert_eq!(m.current_pc(), pc, "pc must be preserved on error");
        assert_eq!(
            m.loader.module.code.len(),
            len,
            "code must be truncated back on error"
        );
        assert_eq!(
            m.cpu.stack().len(),
            depth,
            "operand stack must be left unchanged on error"
        );
        Ok(())
    }

    #[test]
    fn apply_circuit_instruction_measurement_does_not_grow_the_stack() -> eyre::Result<()> {
        // `circuit.measure` pushes its outcome onto the CPU operand stack for
        // bytecode to consume. An injected measurement has no such consumer, so
        // that push must be rolled back: otherwise a paused program resumes with
        // a stray operand, and a REPL session's stack grows without bound.
        let mut m = PPVM::with_qubits(1)?;
        let depth = m.cpu.stack().len();
        m.apply_circuit_instruction(CircuitInstruction::Measure, &[0], &[])?;
        assert_eq!(
            m.cpu.stack().len(),
            depth,
            "injected measurement must not leave its outcome on the stack"
        );
        // The outcome still lands in the measurement record.
        assert_eq!(m.measurement_record().len(), 1);
        Ok(())
    }

    #[test]
    fn tableau_trace_emits_expectation_on_zero_state() {
        // Task 16: `circuit.trace` on the Tableau backend now computes
        // Σ_{P matches pat} ⟨ψ|P|ψ⟩ via `GeneralizedTableau::trace`. On the
        // freshly-initialized |0⟩ state, pattern `Z0` matches the single
        // Pauli Z and ⟨0|Z|0⟩ = 1, so the trace_record gets one entry: 1.0.
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 1;
        module.strings.push("Z0".to_string());
        module
            .code
            .push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                Value::String(0),
            )));
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::Trace));

        let mut machine = PPVM::default();
        machine.load(&module).unwrap();
        machine.init().unwrap();
        machine.step_once().unwrap(); // const.string
        machine.step_once().unwrap(); // circuit.trace
        let trace = machine.trace_record();
        assert_eq!(trace.len(), 1);
        assert!(
            (trace[0] - 1.0).abs() < 1e-12,
            "expected ⟨0|Z|0⟩ = 1.0, got {}",
            trace[0]
        );
    }

    #[test]
    fn paulisum_reset_restores_seeded_observable() -> eyre::Result<()> {
        use vihaco::Reset;

        // Seed the observable `Z` (PauliSum backend), then apply H(0), which
        // conjugates Z -> X in the Heisenberg picture and changes the state.
        // reset() must rebuild the state from the seeded observable.
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 1;
        module.extra.backend = BackendKind::PauliSum;
        module.extra.observable = Some("Z".to_string());

        module
            .code
            .push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                Value::U64(0),
            )));
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::H));

        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;

        let seeded = machine.state_string();
        for _ in 0..module.code.len() {
            machine.step_once()?;
        }
        assert_ne!(
            machine.state_string(),
            seeded,
            "H(0) should have changed the propagated observable"
        );

        machine.reset();
        assert_eq!(
            machine.state_string(),
            seeded,
            "reset must rebuild the state from the seeded observable"
        );
        Ok(())
    }

    #[test]
    fn lossy_paulisum_reset_restores_seeded_observable() -> eyre::Result<()> {
        use vihaco::Reset;

        // Same as `paulisum_reset_restores_seeded_observable`, but through the
        // LossyPauliSum dispatch path.
        let mut module: Module<PPVMInstruction, Value, Type, PPVMDeviceInfo> = Module::default();
        module.extra.n_qubits = 1;
        module.extra.backend = BackendKind::LossyPauliSum;
        module.extra.observable = Some("Z".to_string());

        module
            .code
            .push(PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(
                Value::U64(0),
            )));
        module
            .code
            .push(PPVMInstruction::Circuit(CircuitInstruction::H));

        let mut machine = PPVM::default();
        machine.load(&module)?;
        machine.init()?;

        let seeded = machine.state_string();
        for _ in 0..module.code.len() {
            machine.step_once()?;
        }
        assert_ne!(
            machine.state_string(),
            seeded,
            "H(0) should have changed the propagated observable"
        );

        machine.reset();
        assert_eq!(
            machine.state_string(),
            seeded,
            "reset must rebuild the state from the seeded observable"
        );
        Ok(())
    }
}
