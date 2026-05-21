use std::collections::HashMap;

use chumsky::{Parser, error::Simple, extra};
use vihaco::{
    Type, Value,
    module::Module,
    syntax::{BodyItem, RawForm, RawOperand, Resolve},
};
use vihaco_parser_core::Parse;

use crate::{
    composite::{PPVMDeviceInfo, PPVMInstruction},
    instruction::CircuitInstruction,
};

#[derive(Debug, Clone, PartialEq, vihaco_parser::Parse)]
#[head = "device "]
pub enum PPVMHeader {
    #[token = "circuit.n_qubits"]
    NumQubits(usize),

    #[token = "circuit.coefficient_threshold"]
    CoefficientThrehsold(f64),
}

#[derive(Debug, Default)]
pub struct PPVMResolver {
    strings: Vec<String>,
}

impl PPVMResolver {
    pub fn new() -> Self {
        Self::default()
    }

    // fn intern(&mut self, s: &str) -> u32 {
    //     if let Some(idx) = self.strings.iter().position(|existing| existing == s) {
    //         return idx as u32;
    //     }
    //     let idx = self.strings.len() as u32;
    //     self.strings.push(s.to_string());
    //     idx
    // }

    fn apply_header(info: &mut PPVMDeviceInfo, header: PPVMHeader) -> eyre::Result<()> {
        match header {
            PPVMHeader::NumQubits(n) => {
                info.n_qubits = n;
            }
            PPVMHeader::CoefficientThrehsold(t) => {
                info.coefficient_threshold = t;
            }
        }
        Ok(())
    }

    fn lower_raw(&mut self, raw: RawForm) -> eyre::Result<Vec<PPVMInstruction>> {
        match raw.mnemonic.as_str() {
            "ret" => {
                require_no_operands(&raw)?;
                Ok(vec![vihaco_cpu::Instruction::Return(0).into()])
            }
            other => Err(eyre::eyre!(
                "PPVMResolver: unhandled raw form `{other}` (operands: {:?})",
                raw.operands
            )),
        }
    }
}

impl Resolve<PPVMInstruction, PPVMHeader> for PPVMResolver {
    type Module = Module<PPVMInstruction, Value, Type, PPVMDeviceInfo>;
    fn resolve_module(
        &mut self,
        parsed: vihaco::syntax::ParsedModule<PPVMInstruction, PPVMHeader>,
    ) -> eyre::Result<Self::Module> {
        let mut info = PPVMDeviceInfo::default();
        for header in parsed.headers {
            Self::apply_header(&mut info, header)?;
        }

        let mut code: Vec<PPVMInstruction> = Vec::new();
        let mut labels: HashMap<String, u32> = HashMap::new();
        let mut patches: Vec<(usize, BranchPatch)> = Vec::new();
        for function in parsed.functions {
            for item in function.body {
                match item {
                    BodyItem::Direct(inst) => code.push(inst),
                    BodyItem::Raw(raw) => {
                        if let Some(name) = raw_as_label(&raw) {
                            if labels.insert(name.clone(), code.len() as u32).is_some() {
                                return Err(eyre::eyre!("duplicate label `@{name}`"));
                            }
                            continue;
                        }
                        if let Some(patch) = raw_as_branch(&raw) {
                            let idx = code.len();
                            code.push(patch.placeholder());
                            patches.push((idx, patch));
                            continue;
                        }
                        code.extend(self.lower_raw(raw)?);
                    }
                }
            }
        }
        for (idx, patch) in patches {
            patch.apply(&mut code, idx, &labels)?;
        }

        let mut module = Module::default();
        module.code = code;
        module.strings = std::mem::take(&mut self.strings);
        module.extra = info;
        Ok(module)
    }
}

type E<'src> = extra::Err<Simple<'src, char>>;

impl<'src> Parse<'src> for PPVMInstruction {
    fn parser() -> impl Parser<'src, &'src str, Self, E<'src>> {
        use chumsky::prelude::*;

        let cpu = <vihaco_cpu::Instruction as Parse>::parser().map(PPVMInstruction::Cpu);

        // Reuse the derived parser for all CircuitInstruction variants;
        // just gate it behind the `gate ` keyword.
        let circuit = just("gate")
            .then(text::whitespace().at_least(1))
            .ignore_then(<CircuitInstruction as Parse>::parser())
            .map(PPVMInstruction::Circuit);

        // Try `gate ...` first so CPU doesn't see "gate" as an identifier.
        choice((circuit, cpu))
    }
}

// ---- Everything below is 1:1 copy from Acamar with Acamar -> PPVM renaming ----

fn require_no_operands(raw: &RawForm) -> eyre::Result<()> {
    if !raw.operands.is_empty() {
        return Err(eyre::eyre!(
            "`{}` takes no operands, got {}",
            raw.mnemonic,
            raw.operands.len()
        ));
    }
    Ok(())
}

/// A deferred branch whose target(s) couldn't be resolved at lowering time
/// because the label may appear later in the function body. Patched in a
/// second pass once all labels are known.
#[derive(Debug)]
enum BranchPatch {
    /// `br @target` — fills the `u32` in `cpu::Instruction::Branch`.
    Unconditional(String),
    /// `br @t, @f` / `cond_br @t, @f` — fills both `u32`s in
    /// `cpu::Instruction::ConditionalBranch`.
    Conditional(String, String),
}

/// `@name:` → `Some("name")`. Body parser already emits `@entry:` as a single
/// raw mnemonic with no operands, so the check is purely on the mnemonic
/// shape.
fn raw_as_label(raw: &RawForm) -> Option<String> {
    if !raw.operands.is_empty() {
        return None;
    }
    let m = raw.mnemonic.as_str();
    let stripped = m.strip_prefix('@')?.strip_suffix(':')?;
    if stripped.is_empty() {
        return None;
    }
    Some(stripped.to_string())
}

/// `br @t` / `br @t, @f` / `cond_br @t, @f`.
fn raw_as_branch(raw: &RawForm) -> Option<BranchPatch> {
    let symbols: Vec<&str> = raw
        .operands
        .iter()
        .map(|op| match op {
            RawOperand::Symbol(s) => Some(s.as_str()),
            _ => None,
        })
        .collect::<Option<Vec<_>>>()?;

    match (raw.mnemonic.as_str(), symbols.as_slice()) {
        ("br", [t]) => Some(BranchPatch::Unconditional((*t).to_string())),
        ("br", [t, f]) | ("cond_br", [t, f]) => {
            Some(BranchPatch::Conditional((*t).to_string(), (*f).to_string()))
        }
        _ => None,
    }
}

impl BranchPatch {
    fn placeholder(&self) -> PPVMInstruction {
        match self {
            BranchPatch::Unconditional(_) => vihaco_cpu::Instruction::Branch(u32::MAX).into(),
            BranchPatch::Conditional(_, _) => {
                vihaco_cpu::Instruction::ConditionalBranch(u32::MAX, u32::MAX).into()
            }
        }
    }

    fn apply(
        self,
        code: &mut [PPVMInstruction],
        idx: usize,
        labels: &HashMap<String, u32>,
    ) -> eyre::Result<()> {
        let lookup = |name: &str| {
            labels
                .get(name)
                .copied()
                .ok_or_else(|| eyre::eyre!("undefined label `@{name}`"))
        };
        let resolved = match self {
            BranchPatch::Unconditional(t) => vihaco_cpu::Instruction::Branch(lookup(&t)?).into(),
            BranchPatch::Conditional(t, f) => {
                vihaco_cpu::Instruction::ConditionalBranch(lookup(&t)?, lookup(&f)?).into()
            }
        };
        code[idx] = resolved;
        Ok(())
    }
}
