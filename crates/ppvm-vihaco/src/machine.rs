use num::complex::Complex64;
use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_tableau::data::GeneralizedTableau;
use vihaco::machine::Machine;
use vihaco::observer::stdio::{StdoutEvent, StdoutObserver};
use vihaco::traits::{GetProgramGlobal, ProgramCounter, StackMemory};
use vihaco::{GeneratedMachine, ProgramLoader};
use vihaco_cpu::{CPU, CPUMessage, StepOutcome};

use crate::message::CircuitMessage;
use crate::prelude::{Circuit, CircuitInstruction};

pub type Instruction = PPVM128Instruction;

#[derive(vihaco::Machine)]
pub struct PPVM128 {
    #[program]
    loader: ProgramLoader<PPVM128Instruction>,

    #[device(0x00, resolve_with = resolve_cpu, custom_parser)]
    cpu: CPU,

    #[device(0x01, resolve_with = resolve_circuit)]
    circuit: Circuit<Byte8F64<2>, u128, Vec<(Complex64, u128)>>,

    #[observe(StdoutEvent)]
    stdout: StdoutObserver,
}

impl From<vihaco_cpu::Instruction> for PPVM128Instruction {
    fn from(value: vihaco_cpu::Instruction) -> Self {
        Self::Cpu(value)
    }
}

impl From<CircuitInstruction> for PPVM128Instruction {
    fn from(value: CircuitInstruction) -> Self {
        Self::Circuit(value)
    }
}

impl PPVM128 {
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
            vihaco_cpu::Instruction::Const(v) => {
                self.cpu.stack_push(*v);
                Ok(CPUMessage::None)
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
                todo!()
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
}

impl vihaco::Reset for PPVM128 {
    fn reset(&mut self) {
        self.cpu.reset();
        self.circuit.reset();
        self.loader.pc = 0;
    }
}

impl Machine<Instruction> for PPVM128 {
    type MachineStepResult = StepOutcome;

    fn init(&mut self) -> eyre::Result<()> {
        self.circuit = Circuit {
            tab: GeneralizedTableau::new(10, 1e-10),
        };
        Ok(())
    }

    fn load(
        &mut self,
        module: &vihaco::module::Module<
            Instruction,
            vihaco::Value,
            vihaco::Type,
            vihaco::module::NoInfo,
        >,
    ) -> eyre::Result<()> {
        self.loader.module = module.clone();
        Ok(())
    }

    fn step(&mut self) -> eyre::Result<StepOutcome> {
        let mut ctx = vihaco::ExecContext::new(0);
        let inst = self.peek_instruction()?.clone();
        let outcome = match inst {
            PPVM128Instruction::Cpu(cpu_inst) => {
                let msg = self.resolve_cpu(&cpu_inst)?;
                vihaco::GeneratedComponent::execute_generated(
                    &mut self.cpu,
                    cpu_inst,
                    msg,
                    &mut ctx,
                )?
            }
            PPVM128Instruction::Circuit(circuit_inst) => {
                let msg = self.resolve_circuit(&circuit_inst)?;
                vihaco::GeneratedComponent::execute_generated(
                    &mut self.circuit,
                    circuit_inst,
                    msg,
                    &mut ctx,
                )?;
                StepOutcome::Continue
            }
        };

        if outcome == StepOutcome::Continue {
            if let Some(target) = self.cpu.take_pending_pc() {
                *self.pc_mut() = target;
            } else {
                *self.pc_mut() += 1;
            }
        }

        for event in ctx.into_events() {
            let _ = <Self as GeneratedMachine>::deliver_any(self, event.as_ref());
        }

        Ok(outcome)
    }

    fn run(&mut self) -> eyre::Result<()> {
        self.init()?;
        loop {
            match Machine::step(self)? {
                StepOutcome::Continue => continue,
                StepOutcome::Breakpoint | StepOutcome::Halt => break,
                StepOutcome::Return => return Ok(()),
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use vihaco::{Type, Value, module::Module};

    use super::*;

    #[test]
    fn test_run_ppvm() {
        let mut module: Module<PPVM128Instruction, Value, Type> = Module::default();

        /*
        const.u64 0
        gate h
        */
        let zero = PPVM128Instruction::Cpu(vihaco_cpu::Instruction::Const(Value::U64(0)));
        let one = PPVM128Instruction::Cpu(vihaco_cpu::Instruction::Const(Value::U64(1)));
        module.code.push(zero.clone());
        module
            .code
            .push(PPVM128Instruction::Circuit(CircuitInstruction::H));

        /*
        const.u64 0
        gate t
        */

        module.code.push(zero.clone());
        module
            .code
            .push(PPVM128Instruction::Circuit(CircuitInstruction::T));

        /*
        const.u64 0
        const.u64 1
        gate cnot
        */
        module.code.push(zero.clone());
        module.code.push(one.clone());
        module
            .code
            .push(PPVM128Instruction::Circuit(CircuitInstruction::CNOT));

        let mut machine = PPVM128 {
            loader: ProgramLoader::default(),
            cpu: CPU::default(),
            circuit: Circuit {
                tab: GeneralizedTableau::new(2, 1e-10),
            },
            stdout: StdoutObserver::default(),
        };

        println!("{:?}", module.code);

        machine.load(&module).unwrap();

        for _ in 0..module.code.len() {
            machine.step().unwrap();
        }

        println!("{}", machine.circuit.tab);
        assert_eq!(machine.circuit.tab.coefficients.len(), 2);
    }
}
