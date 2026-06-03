// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashMap;

use chumsky::{Parser, error::Simple, extra};
use vihaco::{
    Type, Value,
    module::Module,
    syntax::{BodyItem, RawForm, RawOperand, Resolve},
};
use vihaco_circuit_isa::CircuitInstruction;
use vihaco_parser_core::Parse;

use crate::composite::{PPVMDeviceInfo, PPVMInstruction};

#[derive(Debug, Clone, PartialEq, vihaco_parser::Parse)]
#[head = "device "]
pub enum PPVMHeader {
    #[token = "circuit.n_qubits"]
    #[delimiters(open = "", close = "", separator = "")]
    NumQubits(usize),

    #[token = "circuit.coefficient_threshold"]
    #[delimiters(open = "", close = "", separator = "")]
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
                let keep = match raw.operands.as_slice() {
                    [] => 0u32,
                    [RawOperand::UInt(n)] => u32::try_from(*n)
                        .map_err(|_| eyre::eyre!("`ret` keep count {n} does not fit in u32"))?,
                    other => {
                        return Err(eyre::eyre!(
                            "`ret` takes 0 or 1 unsigned int operands, got {other:?}"
                        ));
                    }
                };
                Ok(vec![vihaco_cpu::Instruction::Return(keep).into()])
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
        let mut branch_patches: Vec<(usize, BranchPatch)> = Vec::new();
        let mut call_patches: Vec<(usize, CallPatch)> = Vec::new();
        for function in parsed.functions {
            if labels
                .insert(function.name.clone(), code.len() as u32)
                .is_some()
            {
                return Err(eyre::eyre!("duplicate function name `@{}`", function.name));
            }
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
                            branch_patches.push((idx, patch));
                            continue;
                        }
                        if let Some(patch) = raw_as_call(&raw)? {
                            let idx = code.len();
                            code.push(patch.placeholder());
                            call_patches.push((idx, patch));
                            continue;
                        }
                        code.extend(self.lower_raw(raw)?);
                    }
                }
            }
        }
        for (idx, patch) in branch_patches {
            patch.apply(&mut code, idx, &labels)?;
        }
        for (idx, patch) in call_patches {
            patch.apply(&mut code, idx, &labels)?;
        }

        let mut module = Module::default();
        module.code = code;
        module.strings = std::mem::take(&mut self.strings);
        module.extra = info;
        Ok(module)
    }
}

type Err<'src> = extra::Err<Simple<'src, char>>;

impl<'src> Parse<'src> for PPVMInstruction {
    fn parser() -> impl Parser<'src, &'src str, Self, Err<'src>> {
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

/// `call <arity>, @target` — symbolic target resolved in a second pass against
/// the same label table that holds branch targets and function entry points.
#[derive(Debug)]
struct CallPatch {
    arity: u32,
    target: String,
}

/// `call <arity>, @target` → `Some(CallPatch)`. Returns `Ok(None)` for any
/// other mnemonic so the resolver can fall through to `lower_raw`.
fn raw_as_call(raw: &RawForm) -> eyre::Result<Option<CallPatch>> {
    if raw.mnemonic != "call" {
        return Ok(None);
    }
    match raw.operands.as_slice() {
        [RawOperand::UInt(arity), RawOperand::Symbol(target)] => {
            let arity = u32::try_from(*arity)
                .map_err(|_| eyre::eyre!("`call` arity {arity} does not fit in u32"))?;
            Ok(Some(CallPatch {
                arity,
                target: target.clone(),
            }))
        }
        other => Err(eyre::eyre!(
            "`call` expects `<arity:uint>, @<target>`, got operands {other:?}"
        )),
    }
}

impl CallPatch {
    fn placeholder(&self) -> PPVMInstruction {
        vihaco_cpu::Instruction::Call(self.arity, u32::MAX).into()
    }

    fn apply(
        self,
        code: &mut [PPVMInstruction],
        idx: usize,
        labels: &HashMap<String, u32>,
    ) -> eyre::Result<()> {
        let target = labels
            .get(&self.target)
            .copied()
            .ok_or_else(|| eyre::eyre!("undefined function `@{}`", self.target))?;
        code[idx] = vihaco_cpu::Instruction::Call(self.arity, target).into();
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vihaco::syntax::ParsedModule;

    fn parse_module(source: &str) -> ParsedModule<PPVMInstruction, PPVMHeader> {
        ParsedModule::<PPVMInstruction, PPVMHeader>::parser()
            .parse(source)
            .into_result()
            .unwrap_or_else(|e| panic!("parse failed: {e:?}"))
    }

    fn raw(mnemonic: &str, operands: Vec<RawOperand>) -> RawForm {
        RawForm {
            mnemonic: mnemonic.to_string(),
            operands,
        }
    }

    // ─── Header parsing ───────────────────────────────────────────────────

    #[test]
    fn header_parses_n_qubits() {
        let got = <PPVMHeader as Parse>::parser()
            .parse("device circuit.n_qubits 5")
            .into_result()
            .unwrap_or_else(|e| panic!("parse failed: {e:?}"));
        assert_eq!(got, PPVMHeader::NumQubits(5));
    }

    #[test]
    fn header_parses_coefficient_threshold() {
        let got = <PPVMHeader as Parse>::parser()
            .parse("device circuit.coefficient_threshold 1e-10")
            .into_result()
            .unwrap_or_else(|e| panic!("parse failed: {e:?}"));
        assert_eq!(got, PPVMHeader::CoefficientThrehsold(1e-10));
    }

    #[test]
    fn header_n_qubits_rejects_extra_operand() {
        // The variant has exactly one field, so the parser consumes one
        // integer. Wrapped in a full module, a second integer must trip the
        // module-level parser.
        let result = ParsedModule::<PPVMInstruction, PPVMHeader>::parser()
            .parse(
                "device circuit.n_qubits 5 6;\n\
                 fn @main() { ret }\n",
            )
            .into_result();
        assert!(result.is_err(), "expected parse error, got {result:?}");
    }

    #[test]
    fn header_coefficient_threshold_rejects_extra_operand() {
        let result = ParsedModule::<PPVMInstruction, PPVMHeader>::parser()
            .parse(
                "device circuit.coefficient_threshold 1e-10 0.5;\n\
                 fn @main() { ret }\n",
            )
            .into_result();
        assert!(result.is_err(), "expected parse error, got {result:?}");
    }

    #[test]
    fn apply_header_sets_n_qubits() {
        let mut info = PPVMDeviceInfo::default();
        PPVMResolver::apply_header(&mut info, PPVMHeader::NumQubits(7)).unwrap();
        assert_eq!(info.n_qubits, 7);
    }

    #[test]
    fn apply_header_sets_coefficient_threshold() {
        let mut info = PPVMDeviceInfo::default();
        PPVMResolver::apply_header(&mut info, PPVMHeader::CoefficientThrehsold(5e-6)).unwrap();
        assert_eq!(info.coefficient_threshold, 5e-6);
    }

    // ─── PPVMInstruction parser dispatch ──────────────────────────────────

    #[test]
    fn ppvm_instruction_parses_cpu_const() {
        let got = <PPVMInstruction as Parse>::parser()
            .parse("const.u64 7")
            .into_result()
            .unwrap();
        assert!(matches!(
            got,
            PPVMInstruction::Cpu(vihaco_cpu::Instruction::Const(Value::U64(7)))
        ));
    }

    #[test]
    fn ppvm_instruction_parses_gate_h() {
        let got = <PPVMInstruction as Parse>::parser()
            .parse("gate h")
            .into_result()
            .unwrap();
        assert!(matches!(
            got,
            PPVMInstruction::Circuit(CircuitInstruction::H)
        ));
    }

    #[test]
    fn ppvm_instruction_parses_gate_cnot() {
        let got = <PPVMInstruction as Parse>::parser()
            .parse("gate cnot")
            .into_result()
            .unwrap();
        assert!(matches!(
            got,
            PPVMInstruction::Circuit(CircuitInstruction::CNOT)
        ));
    }

    #[test]
    fn ppvm_instruction_parses_gate_measure() {
        let got = <PPVMInstruction as Parse>::parser()
            .parse("gate measure")
            .into_result()
            .unwrap();
        assert!(matches!(
            got,
            PPVMInstruction::Circuit(CircuitInstruction::Measure)
        ));
    }

    #[test]
    fn ppvm_instruction_parses_gate_rx() {
        let got = <PPVMInstruction as Parse>::parser()
            .parse("gate rx")
            .into_result()
            .unwrap();
        assert!(matches!(
            got,
            PPVMInstruction::Circuit(CircuitInstruction::RX)
        ));
    }

    #[test]
    fn ppvm_instruction_rejects_bare_circuit_token_without_gate_prefix() {
        // `h` on its own must not parse as Circuit(H) — only `gate h` does.
        // Without `gate `, the CPU parser is tried, which should reject
        // `h` (not a CPU mnemonic).
        let result = <PPVMInstruction as Parse>::parser()
            .parse("h")
            .into_result();
        assert!(result.is_err(), "expected parse error, got {result:?}");
    }

    // ─── lower_raw ────────────────────────────────────────────────────────

    #[test]
    fn lower_raw_ret_emits_return_zero() {
        let mut r = PPVMResolver::new();
        let out = r.lower_raw(raw("ret", vec![])).unwrap();
        assert_eq!(out.len(), 1);
        assert!(matches!(
            out[0],
            PPVMInstruction::Cpu(vihaco_cpu::Instruction::Return(0))
        ));
    }

    #[test]
    fn lower_raw_ret_with_uint_operand_emits_return_n() {
        let mut r = PPVMResolver::new();
        let out = r.lower_raw(raw("ret", vec![RawOperand::UInt(2)])).unwrap();
        assert_eq!(out.len(), 1);
        assert!(matches!(
            out[0],
            PPVMInstruction::Cpu(vihaco_cpu::Instruction::Return(2))
        ));
    }

    #[test]
    fn lower_raw_ret_with_non_uint_operand_errors() {
        let mut r = PPVMResolver::new();
        let err = r
            .lower_raw(raw("ret", vec![RawOperand::Symbol("foo".into())]))
            .unwrap_err();
        assert!(
            err.to_string().contains("`ret` takes 0 or 1 unsigned int"),
            "err: {err}"
        );
    }

    #[test]
    fn lower_raw_unknown_mnemonic_errors() {
        let mut r = PPVMResolver::new();
        let err = r.lower_raw(raw("nope", vec![])).unwrap_err();
        assert!(err.to_string().contains("unhandled raw form"), "err: {err}");
    }

    // ─── End-to-end resolver behaviour ────────────────────────────────────

    #[test]
    fn resolver_populates_device_info_from_headers() {
        let parsed = parse_module(
            "device circuit.n_qubits 3;\n\
             device circuit.coefficient_threshold 1e-8;\n\
             fn @main() { ret }\n",
        );
        let m = PPVMResolver::new().resolve_module(parsed).unwrap();
        assert_eq!(m.extra.n_qubits, 3);
        assert_eq!(m.extra.coefficient_threshold, 1e-8);
    }

    #[test]
    fn resolver_lowers_simple_bell_body() {
        // Smoke test the whole pipeline on a tiny bell-like body.
        let parsed = parse_module(
            "device circuit.n_qubits 2;\n\
             fn @main() {\n\
                 const.u64 0\n\
                 gate h\n\
                 const.u64 0\n\
                 const.u64 1\n\
                 gate cnot\n\
                 ret\n\
             }\n",
        );
        let m = PPVMResolver::new().resolve_module(parsed).unwrap();
        // const.u64 0 / gate h / const.u64 0 / const.u64 1 / gate cnot / ret
        assert_eq!(m.code.len(), 6);
        assert!(matches!(
            m.code[1],
            PPVMInstruction::Circuit(CircuitInstruction::H)
        ));
        assert!(matches!(
            m.code[4],
            PPVMInstruction::Circuit(CircuitInstruction::CNOT)
        ));
        assert!(matches!(
            m.code[5],
            PPVMInstruction::Cpu(vihaco_cpu::Instruction::Return(0))
        ));
    }

    #[test]
    fn resolver_resolves_forward_branch_targets() {
        let parsed = parse_module(
            "fn @main() {\n\
                 @loop:\n\
                     br @done\n\
                 @done:\n\
                     ret\n\
             }\n",
        );
        let m = PPVMResolver::new().resolve_module(parsed).unwrap();
        assert!(matches!(
            m.code[0],
            PPVMInstruction::Cpu(vihaco_cpu::Instruction::Branch(1))
        ));
    }

    #[test]
    fn resolver_resolves_conditional_branch_with_two_targets() {
        let parsed = parse_module(
            "fn @main() {\n\
                 @head:\n\
                     br @head, @exit\n\
                 @exit:\n\
                     ret\n\
             }\n",
        );
        let m = PPVMResolver::new().resolve_module(parsed).unwrap();
        assert!(matches!(
            m.code[0],
            PPVMInstruction::Cpu(vihaco_cpu::Instruction::ConditionalBranch(0, 1))
        ));
    }

    #[test]
    fn resolver_rejects_undefined_branch_target() {
        let parsed = parse_module(
            "fn @main() {\n\
                 br @missing\n\
                 ret\n\
             }\n",
        );
        let err = PPVMResolver::new().resolve_module(parsed).unwrap_err();
        assert!(
            err.to_string().contains("undefined label `@missing`"),
            "err: {err}"
        );
    }

    #[test]
    fn resolver_rejects_duplicate_label() {
        let parsed = parse_module(
            "fn @main() {\n\
                 @same:\n\
                     ret\n\
                 @same:\n\
                     ret\n\
             }\n",
        );
        let err = PPVMResolver::new().resolve_module(parsed).unwrap_err();
        assert!(
            err.to_string().contains("duplicate label `@same`"),
            "err: {err}"
        );
    }
}
