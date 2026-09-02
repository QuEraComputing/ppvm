// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashMap;

use chumsky::{IterParser, Parser};
use vihaco::module::{FunctionInfo, LabelInfo, LocalModule, Parameter, Signature};
use vihaco::syntax::{Param, ParsedFunction, ParsedModule, Resolve, skip};
use vihaco::{Parse, Type, Value};
use vihaco_circuit_isa::{CircuitInstruction, CircuitSurfaceInstruction};
use vihaco_cpu::{RuntimeInstruction as CpuRuntime, SurfaceInstruction as CpuSurface};
use vihaco_cpu::{SurfaceType, SurfaceValue};
use vihaco_parser::{BareToken, Ident, QuotedString};

use crate::composite::ppvm_module as ppvm;
use crate::composite::{BackendKind, PPVMDeviceInfo, PPVMInstruction};

#[derive(Debug, Clone, PartialEq, vihaco_parser_derive::Parse)]
#[syntax_class(value)]
pub enum BackendKindSyntax {
    #[pattern = "`tableau`"]
    Tableau,
    #[pattern = "`paulisum`"]
    PauliSum,
    #[pattern = "`lossy_paulisum`"]
    LossyPauliSum,
}

impl From<BackendKindSyntax> for BackendKind {
    fn from(value: BackendKindSyntax) -> Self {
        match value {
            BackendKindSyntax::Tableau => Self::Tableau,
            BackendKindSyntax::PauliSum => Self::PauliSum,
            BackendKindSyntax::LossyPauliSum => Self::LossyPauliSum,
        }
    }
}

#[derive(Debug, Clone, PartialEq, vihaco_parser_derive::Parse)]
#[syntax_class(metadata, head = "device")]
pub enum PPVMHeader {
    #[pattern = "circuit.n_qubits $0"]
    NumQubits(usize),
    #[pattern = "circuit.coefficient_threshold $0"]
    CoefficientThreshold(f64),
    #[pattern = "circuit.backend $0"]
    Backend(BackendKindSyntax),
    #[pattern = "circuit.observable $0"]
    Observable(Ident),
    #[pattern = "circuit.max_pauli_weight $0"]
    MaxPauliWeight(usize),
}

#[derive(Debug, Default)]
pub struct PPVMResolver {
    strings: Vec<String>,
}

impl PPVMResolver {
    pub fn new() -> Self {
        Self::default()
    }

    fn apply_header(info: &mut PPVMDeviceInfo, header: PPVMHeader) -> eyre::Result<()> {
        match header {
            PPVMHeader::NumQubits(n) => info.n_qubits = n,
            PPVMHeader::CoefficientThreshold(t) => info.coefficient_threshold = t,
            PPVMHeader::Backend(b) => info.backend = b.into(),
            PPVMHeader::Observable(s) => info.observable = Some(s.as_str().to_owned()),
            PPVMHeader::MaxPauliWeight(w) => info.max_pauli_weight = Some(w),
        }
        Ok(())
    }

    fn intern(&mut self, value: &str) -> u32 {
        if let Some(index) = self.strings.iter().position(|item| item == value) {
            return index as u32;
        }
        let index = self.strings.len() as u32;
        self.strings.push(value.to_owned());
        index
    }

    fn runtime_type(ty: SurfaceType) -> Type {
        match ty {
            SurfaceType::Undefined => Type::Undefined,
            SurfaceType::String => Type::String,
            SurfaceType::Bool => Type::Bool,
            SurfaceType::I64 => Type::I64,
            SurfaceType::U32 => Type::U32,
            SurfaceType::U64 => Type::U64,
            SurfaceType::F64 => Type::F64,
            SurfaceType::FunctionRef => Type::FunctionRef,
            SurfaceType::HeapRef => Type::HeapRef,
        }
    }

    fn lower_value(&mut self, ty: SurfaceType, value: SurfaceValue) -> eyre::Result<Value> {
        let text = match value {
            SurfaceValue::Quoted(QuotedString(value)) => {
                return Ok(Value::String(self.intern(&value)));
            }
            SurfaceValue::Bare(BareToken(value)) => value,
        };
        Ok(match ty {
            SurfaceType::Undefined => Value::Undefined,
            SurfaceType::String => Value::String(self.intern(&text)),
            SurfaceType::Bool => Value::Bool(text.parse()?),
            SurfaceType::I64 => Value::I64(text.parse()?),
            SurfaceType::U32 => Value::U32(text.parse()?),
            SurfaceType::U64 => Value::U64(text.parse()?),
            SurfaceType::F64 => Value::F64(text.parse()?),
            SurfaceType::FunctionRef => Value::FunctionRef(text.parse()?),
            SurfaceType::HeapRef => Value::HeapRef(text.parse()?),
        })
    }

    fn lower_cpu(
        &mut self,
        instruction: CpuSurface,
        labels: &HashMap<String, u32>,
        functions: &HashMap<String, u32>,
    ) -> eyre::Result<Option<CpuRuntime>> {
        use CpuRuntime as R;
        use CpuSurface as S;
        let runtime = match instruction {
            S::Span(a, b, c) => R::Span(a, b, c),
            S::Label(_) => return Ok(None),
            S::FunctionStart => R::FunctionStart,
            S::FunctionEnd => R::FunctionEnd,
            S::Breakpoint => R::Breakpoint,
            S::Branch(name) => R::Branch(
                *labels
                    .get(name.as_str())
                    .ok_or_else(|| eyre::eyre!("undefined label `@{}`", name.as_str()))?,
            ),
            S::ConditionalBranch(a, b) => R::ConditionalBranch(
                *labels
                    .get(a.as_str())
                    .ok_or_else(|| eyre::eyre!("undefined label `@{}`", a.as_str()))?,
                *labels
                    .get(b.as_str())
                    .ok_or_else(|| eyre::eyre!("undefined label `@{}`", b.as_str()))?,
            ),
            S::Return(n) => R::Return(n),
            S::IndirectCall => R::IndirectCall,
            S::Call(arity, name) => R::Call(
                arity,
                *functions
                    .get(name.as_str())
                    .ok_or_else(|| eyre::eyre!("undefined function `@{}`", name.as_str()))?,
            ),
            S::Halt => R::Halt,
            S::Print => R::Print,
            S::Load(ty, slot) => R::Load(Self::runtime_type(ty), slot),
            S::Store(ty, slot) => R::Store(Self::runtime_type(ty), slot),
            S::Dup => R::Dup,
            S::HeapAlloc(n) => R::HeapAlloc(n),
            S::GetItem => R::GetItem,
            S::HeapDealloc => R::HeapDealloc,
            S::Const(ty, value) => {
                let value = self.lower_value(ty, value)?;
                R::Const(Self::runtime_type(ty), value)
            }
            S::Add(ty) => R::Add(Self::runtime_type(ty)),
            S::Sub(ty) => R::Sub(Self::runtime_type(ty)),
            S::Mul(ty) => R::Mul(Self::runtime_type(ty)),
            S::Div(ty) => R::Div(Self::runtime_type(ty)),
            S::Rem(ty) => R::Rem(Self::runtime_type(ty)),
            S::Neg(ty) => R::Neg(Self::runtime_type(ty)),
            S::Shl(ty) => R::Shl(Self::runtime_type(ty)),
            S::Shr(ty) => R::Shr(Self::runtime_type(ty)),
            S::Rol(ty) => R::Rol(Self::runtime_type(ty)),
            S::Ror(ty) => R::Ror(Self::runtime_type(ty)),
            S::BitAnd(ty) => R::BitAnd(Self::runtime_type(ty)),
            S::BitOr(ty) => R::BitOr(Self::runtime_type(ty)),
            S::BitXor(ty) => R::BitXor(Self::runtime_type(ty)),
            S::Not => R::Not,
            S::And => R::And,
            S::Or => R::Or,
            S::Xor => R::Xor,
            S::Eq(ty) => R::Eq(Self::runtime_type(ty)),
            S::Ne(ty) => R::Ne(Self::runtime_type(ty)),
            S::Lt(ty) => R::Lt(Self::runtime_type(ty)),
            S::Gt(ty) => R::Gt(Self::runtime_type(ty)),
            S::Le(ty) => R::Le(Self::runtime_type(ty)),
            S::Ge(ty) => R::Ge(Self::runtime_type(ty)),
        };
        Ok(Some(runtime))
    }

    fn lower_circuit(instruction: CircuitSurfaceInstruction) -> CircuitInstruction {
        use CircuitInstruction as R;
        use CircuitSurfaceInstruction as S;
        match instruction {
            S::TwoQubitPauliError => R::TwoQubitPauliError,
            S::Truncate => R::Truncate,
            S::Trace => R::Trace,
            S::X => R::X,
            S::Y => R::Y,
            S::Z => R::Z,
            S::H => R::H,
            S::SqrtXAdj => R::SqrtXAdj,
            S::SqrtX => R::SqrtX,
            S::SqrtYAdj => R::SqrtYAdj,
            S::SqrtY => R::SqrtY,
            S::SAdj => R::SAdj,
            S::S => R::S,
            S::CNOT => R::CNOT,
            S::CZ => R::CZ,
            S::TAdj => R::TAdj,
            S::T => R::T,
            S::RXX => R::RXX,
            S::RYY => R::RYY,
            S::RZZ => R::RZZ,
            S::RX => R::RX,
            S::RY => R::RY,
            S::RZ => R::RZ,
            S::U3 => R::U3,
            S::Measure => R::Measure,
            S::Reset => R::Reset,
            S::R => R::R,
            S::Loss => R::Loss,
            S::CorrelatedLoss => R::CorrelatedLoss,
            S::PauliError => R::PauliError,
            S::Depolarize2 => R::Depolarize2,
            S::Depolarize => R::Depolarize,
        }
    }
}

impl Resolve<ppvm::syntax::Instruction, SurfaceType, Vec<PPVMHeader>> for PPVMResolver {
    type Module = LocalModule<PPVMInstruction, Value, Type, PPVMDeviceInfo>;

    fn resolve_module(
        &mut self,
        parsed: ParsedModule<ppvm::syntax::Instruction, SurfaceType, Vec<PPVMHeader>>,
    ) -> eyre::Result<Self::Module> {
        let mut module = LocalModule::default();
        for header in parsed.header {
            Self::apply_header(&mut module.extra, header)?;
        }

        let mut functions = HashMap::new();
        let mut function_address = 0u32;
        for function in &parsed.functions {
            functions.insert(function.name.as_str().to_owned(), function_address);
            function_address += function
                .body
                .iter()
                .filter(|instruction| {
                    !matches!(
                        instruction,
                        ppvm::syntax::Instruction::Cpu(CpuSurface::Label(_))
                    )
                })
                .count() as u32;
        }

        let mut labels = HashMap::new();
        let mut address = 0u32;
        for function in &parsed.functions {
            for instruction in &function.body {
                if let ppvm::syntax::Instruction::Cpu(CpuSurface::Label(name)) = instruction {
                    labels.insert(name.as_str().to_owned(), address);
                    module.labels.push(LabelInfo {
                        address,
                        name: self.intern(name.as_str()),
                    });
                } else {
                    address += 1;
                }
            }
        }

        for function in parsed.functions {
            let start_address = module.code.len() as u32;
            for instruction in function.body {
                let runtime = match instruction {
                    ppvm::syntax::Instruction::Cpu(cpu) => self
                        .lower_cpu(cpu, &labels, &functions)?
                        .map(PPVMInstruction::Cpu),
                    ppvm::syntax::Instruction::Circuit(circuit) => {
                        Some(PPVMInstruction::Circuit(Self::lower_circuit(circuit)))
                    }
                };
                if let Some(runtime) = runtime {
                    module.code.push(runtime);
                }
            }
            let end_address = module.code.len() as u32;
            module.functions.push(FunctionInfo {
                name: self.intern(function.name.as_str()),
                signature: Signature {
                    params: function
                        .params
                        .into_iter()
                        .map(|Param { name, ty }| Parameter {
                            name: self.intern(name.as_str()),
                            ty: Self::runtime_type(ty),
                        })
                        .collect(),
                    ret: function
                        .return_ty
                        .map(Self::runtime_type)
                        .into_iter()
                        .collect(),
                },
                local_count: 0,
                start_address,
                end_address,
                file: 0,
            });
        }
        module.main_function = functions.get("main").copied();
        module.strings = std::mem::take(&mut self.strings);
        Ok(module)
    }
}

pub fn parse_headers(source: &str) -> eyre::Result<(Vec<PPVMHeader>, String)> {
    let mut headers = Vec::new();
    let mut body = String::new();
    for line in source.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with("device ") {
            let header = trimmed.strip_suffix(';').unwrap_or(trimmed);
            let parsed = PPVMHeader::parser()
                .parse(header)
                .into_result()
                .map_err(|errors| eyre::eyre!("invalid device header `{header}`: {errors:?}"))?;
            headers.push(parsed);
        } else {
            body.push_str(line);
            body.push('\n');
        }
    }
    Ok((headers, body))
}

pub fn parse_functions(
    source: &str,
) -> eyre::Result<Vec<ParsedFunction<ppvm::syntax::Instruction, SurfaceType>>> {
    skip()
        .ignore_then(
            ParsedFunction::<ppvm::syntax::Instruction, SurfaceType>::parser()
                .repeated()
                .collect::<Vec<_>>(),
        )
        .then_ignore(skip())
        .parse(source)
        .into_result()
        .map_err(|errors| eyre::eyre!("parsing functions failed: {errors:?}"))
}
