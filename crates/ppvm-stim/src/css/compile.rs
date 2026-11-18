// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashMap;

use stim_parser::prelude::{
    ExtendedInstruction, ExtendedProgram, GateName, GateOp, LineMap, MeasureName, MeasureOp,
    NoiseName, NoiseOp, Span, Target,
};

use super::{CssError, CssProgram};
use crate::executor::qubit;
use crate::validate::{ExecError, validate};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Basis {
    Z,
    X,
}

/// One step of the frame program. `log1m_*` fields hold `ln(1 - p)`.
#[derive(Debug, Clone)]
pub(crate) enum Op {
    Reset {
        q: usize,
        basis: Basis,
    },
    ResetState(usize),
    ResetLoss(usize),
    ResetLeakage(usize),
    Measure {
        q: usize,
        basis: Basis,
        reset: bool,
        log1m_flip: f64,
    },
    Pad {
        bit: bool,
        log1m_flip: f64,
    },
    H(Vec<usize>),
    Pauli {
        q: usize,
        x: bool,
        z: bool,
    },
    Feedback {
        lookback: usize,
        q: usize,
        x: bool,
        z: bool,
    },
    Cx(Vec<[usize; 2]>),
    Cz(Vec<[usize; 2]>),
    Pauli1(Box<Channel<usize>>),
    Pauli2(Box<Channel<[usize; 2]>>),
    Loss(Box<Channel<usize>>),
    Leakage(Box<Channel<usize>>),
    CorrelatedLoss(Box<PairLoss>),
    Repeat(Box<Repeat>),
}

#[derive(Debug, Clone)]
pub(crate) struct Repeat {
    pub count: u64,
    pub body: Vec<Op>,
}

/// Events at rate `p` per target and shot; each event's outcome is drawn from `cdf`.
#[derive(Debug, Clone)]
pub(crate) struct Channel<T> {
    pub targets: Vec<T>,
    pub log1m_p: f64,
    pub cdf: Vec<f64>,
}

impl<T> Channel<T> {
    pub fn pick(&self, u: f64) -> usize {
        self.cdf
            .iter()
            .position(|&c| u < c)
            .unwrap_or(self.cdf.len() - 1)
    }
}

/// `I_ERROR[correlated_loss](p0, p1, p2)`: events at rate `p = max(p0 + 2 p1, p2)`,
/// thinned per shot with the thresholds below (in units of `p`).
#[derive(Debug, Clone)]
pub(crate) struct PairLoss {
    pub pairs: Vec<[usize; 2]>,
    pub log1m_p: f64,
    pub both: f64,
    pub first: f64,
    pub second: f64,
    pub other: f64,
}

impl CssProgram {
    /// Compile a Stim program for the CSS frame sampler.
    ///
    /// Fails with [`CssError::Exec`] if the program does not pass [`validate`] or uses an
    /// instruction outside the CSS gate set, and with [`CssError::NotCss`] if no choice of
    /// Hadamard-conjugated qubit segments turns every `CX` / `CZ` into a CSS gate.
    pub fn compile(program: &ExtendedProgram) -> Result<Self, CssError> {
        validate(program)?;
        let num_qubits = program.num_qubits();
        let mut labels = Labels::default();
        let mut compiler = Compiler {
            line_map: &program.line_map,
            ops: Vec::new(),
            segment: (0..num_qubits).map(|_| labels.add()).collect(),
            hadamard: vec![false; num_qubits],
            labels,
        };
        compiler.block(&program.instructions)?;
        Ok(Self {
            ops: compiler.ops,
            num_qubits,
            num_measurements: program.measurement_count(),
        })
    }
}

/// Union-find over the Hadamard labels of qubit segments (a segment starts at a reset or a
/// measurement), with the parity of each label relative to its root.
#[derive(Default)]
struct Labels {
    parent: Vec<usize>,
    parity: Vec<bool>,
    size: Vec<usize>,
}

impl Labels {
    fn add(&mut self) -> usize {
        self.parent.push(self.parent.len());
        self.parity.push(false);
        self.size.push(1);
        self.parent.len() - 1
    }

    fn find(&mut self, a: usize) -> (usize, bool) {
        let parent = self.parent[a];
        if parent == a {
            return (a, false);
        }
        let (root, parity) = self.find(parent);
        self.parent[a] = root;
        self.parity[a] ^= parity;
        (root, self.parity[a])
    }

    /// Require `label(a) ^ label(b) == parity`; false if that contradicts earlier requirements.
    fn relate(&mut self, a: usize, b: usize, parity: bool) -> bool {
        let ((ra, pa), (rb, pb)) = (self.find(a), self.find(b));
        if ra == rb {
            return pa ^ pb == parity;
        }
        let (small, large) = if self.size[ra] < self.size[rb] {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.parent[small] = large;
        self.parity[small] = pa ^ pb ^ parity;
        self.size[large] += self.size[small];
        true
    }
}

struct Compiler<'a> {
    line_map: &'a LineMap,
    ops: Vec<Op>,
    labels: Labels,
    segment: Vec<usize>,
    hadamard: Vec<bool>,
}

impl Compiler<'_> {
    fn block(&mut self, instructions: &[ExtendedInstruction]) -> Result<(), CssError> {
        instructions.iter().try_for_each(|i| self.instruction(i))
    }

    fn instruction(&mut self, instruction: &ExtendedInstruction) -> Result<(), CssError> {
        use ExtendedInstruction as I;
        match instruction {
            I::Gate(op) => self.gate(op),
            I::Noise(op) => self.noise(op),
            I::Measure(op) => self.measure(op),
            I::Annotation(_) => Ok(()),
            I::Loss { p, targets, span } => {
                if let Some(c) = self.channel(targets, &[*p], "I_ERROR[loss]", *span)? {
                    self.ops.push(Op::Loss(c));
                }
                Ok(())
            }
            I::Leakage {
                p0,
                p1,
                targets,
                span,
            } => {
                if let Some(c) = self.channel(targets, &[*p0, *p1], "I_ERROR[leakage]", *span)? {
                    self.ops.push(Op::Leakage(c));
                }
                Ok(())
            }
            I::CorrelatedLoss { ps, targets, span } => self.correlated_loss(*ps, targets, *span),
            I::ResetLoss { targets, .. } => {
                targets
                    .iter()
                    .for_each(|&q| self.ops.push(Op::ResetLoss(q)));
                Ok(())
            }
            I::ResetLeakage { targets, .. } => {
                targets
                    .iter()
                    .for_each(|&q| self.ops.push(Op::ResetLeakage(q)));
                Ok(())
            }
            I::ResetState { targets, .. } => {
                for &q in targets {
                    self.ops.push(Op::ResetState(q));
                    self.new_segment(q);
                }
                Ok(())
            }
            I::MPad { prob, bits, .. } => {
                let log1m_flip = log1m(prob.unwrap_or(0.0));
                bits.iter()
                    .for_each(|&bit| self.ops.push(Op::Pad { bit, log1m_flip }));
                Ok(())
            }
            I::Repeat { count, body, .. } => self.repeat(*count, body),
            I::Mpp(op) => Err(self.unsupported("MPP", op.span)),
            I::T { span, .. } => Err(self.unsupported("T", *span)),
            I::TDag { span, .. } => Err(self.unsupported("T_DAG", *span)),
            I::Rotation { axis, span, .. } => Err(self.unsupported(&format!("R_{axis:?}"), *span)),
            I::U3 { span, .. } => Err(self.unsupported("U3", *span)),
        }
    }

    fn gate(&mut self, op: &GateOp) -> Result<(), CssError> {
        use GateName::*;
        let qubits = op.targets.iter().map(|&t| qubit(t));
        match op.name {
            Reset | ResetZ => qubits.for_each(|q| self.reset(q, Basis::Z)),
            ResetX => qubits.for_each(|q| self.reset(q, Basis::X)),
            X => qubits.for_each(|q| self.pauli(q, true, false)),
            Y => qubits.for_each(|q| self.pauli(q, true, true)),
            Z => qubits.for_each(|q| self.pauli(q, false, true)),
            H | HXZ => self.hadamards(qubits.collect()),
            // SQRT_Y = H Z and SQRT_Y_DAG = H X: a Hadamard after a Pauli.
            SqrtY => qubits.for_each(|q| {
                self.pauli(q, false, true);
                self.hadamards(vec![q]);
            }),
            SqrtYDag => qubits.for_each(|q| {
                self.pauli(q, true, false);
                self.hadamards(vec![q]);
            }),
            Identity => {}
            CX | ZCX | CNot => return self.controlled(op, true, false),
            CY | ZCY => return self.controlled(op, true, true),
            CZ | ZCZ => return self.controlled(op, false, true),
            _ => return Err(self.unsupported(op.name.canonical_name(), op.span)),
        }
        Ok(())
    }

    /// A controlled Pauli: classical feedback from a record, or a `CX` / `CZ` between qubits.
    /// Consecutive qubit pairs become one op.
    fn controlled(&mut self, op: &GateOp, x: bool, z: bool) -> Result<(), CssError> {
        let mut pairs = Vec::new();
        for pair in op.targets.as_chunks::<2>().0 {
            let q = qubit(pair[1]);
            match pair[0] {
                Target::Rec(lookback) => {
                    self.push_pairs(std::mem::take(&mut pairs), z);
                    self.ops.push(Op::Feedback { lookback, q, x, z });
                }
                Target::Qubit(_) if x && z => {
                    return Err(self.unsupported(op.name.canonical_name(), op.span));
                }
                Target::Qubit(c) => {
                    // A CX needs equal Hadamard labels on both qubits, a CZ different ones.
                    let parity = self.hadamard[c] ^ self.hadamard[q] ^ z;
                    if !self.labels.relate(self.segment[c], self.segment[q], parity) {
                        return Err(CssError::NotCss {
                            line: op.span.line(self.line_map),
                        });
                    }
                    pairs.push([c, q]);
                }
            }
        }
        self.push_pairs(pairs, z);
        Ok(())
    }

    fn push_pairs(&mut self, pairs: Vec<[usize; 2]>, z: bool) {
        if !pairs.is_empty() {
            self.ops.push(if z { Op::Cz(pairs) } else { Op::Cx(pairs) });
        }
    }

    fn noise(&mut self, op: &NoiseOp) -> Result<(), CssError> {
        use NoiseName::*;
        let a = &op.args;
        match op.name {
            XError => self.pauli1(op, [a[0], 0.0, 0.0]),
            YError => self.pauli1(op, [0.0, a[0], 0.0]),
            ZError => self.pauli1(op, [0.0, 0.0, a[0]]),
            Depolarize1 => self.pauli1(op, [a[0] / 3.0; 3]),
            PauliChannel1 => self.pauli1(op, [a[0], a[1], a[2]]),
            Depolarize2 => self.pauli2(op, &[a[0] / 15.0; 15]),
            PauliChannel2 => self.pauli2(op, a),
            IError
            | HeraldedErase
            | HeraldedPauliChannel1
            | CorrelatedError
            | ElseCorrelatedError => Err(self.unsupported(op.name.canonical_name(), op.span)),
        }
    }

    fn pauli1(&mut self, op: &NoiseOp, probs: [f64; 3]) -> Result<(), CssError> {
        if let Some(c) = self.channel(&op.targets, &probs, op.name.canonical_name(), op.span)? {
            self.ops.push(Op::Pauli1(c));
        }
        Ok(())
    }

    fn pauli2(&mut self, op: &NoiseOp, probs: &[f64]) -> Result<(), CssError> {
        let pairs = op.targets.as_chunks::<2>().0;
        if let Some(c) = self.channel(pairs, probs, op.name.canonical_name(), op.span)? {
            self.ops.push(Op::Pauli2(c));
        }
        Ok(())
    }

    fn measure(&mut self, op: &MeasureOp) -> Result<(), CssError> {
        let (basis, reset) = match op.name {
            MeasureName::M | MeasureName::MZ => (Basis::Z, false),
            MeasureName::MR => (Basis::Z, true),
            MeasureName::MX => (Basis::X, false),
            MeasureName::MRX => (Basis::X, true),
            _ => return Err(self.unsupported(op.name.canonical_name(), op.span)),
        };
        let log1m_flip = log1m(op.args.first().copied().unwrap_or(0.0));
        for &q in &op.targets {
            self.ops.push(Op::Measure {
                q,
                basis,
                reset,
                log1m_flip,
            });
            self.new_segment(q);
        }
        Ok(())
    }

    fn correlated_loss(
        &mut self,
        [p0, p1, p2]: [f64; 3],
        targets: &[(usize, usize)],
        span: Span,
    ) -> Result<(), CssError> {
        let name = "I_ERROR[correlated_loss]";
        for p in [p0, p1, p2, p0 + 2.0 * p1] {
            self.check(p, name, span)?;
        }
        let p = (p0 + 2.0 * p1).max(p2).min(1.0);
        if p > 0.0 && !targets.is_empty() {
            self.ops.push(Op::CorrelatedLoss(Box::new(PairLoss {
                pairs: targets.iter().map(|&(a, b)| [a, b]).collect(),
                log1m_p: log1m(p),
                both: p0 / p,
                first: (p0 + p1) / p,
                second: (p0 + 2.0 * p1) / p,
                other: p2 / p,
            })));
        }
        Ok(())
    }

    fn channel<T: Copy>(
        &self,
        targets: &[T],
        probs: &[f64],
        name: &str,
        span: Span,
    ) -> Result<Option<Box<Channel<T>>>, CssError> {
        for &p in probs {
            self.check(p, name, span)?;
        }
        let total: f64 = probs.iter().sum();
        self.check(total, name, span)?;
        if total == 0.0 || targets.is_empty() {
            return Ok(None);
        }
        let cdf = probs
            .iter()
            .scan(0.0, |acc, p| {
                *acc += p;
                Some(*acc / total)
            })
            .collect();
        Ok(Some(Box::new(Channel {
            targets: targets.to_vec(),
            log1m_p: log1m(total.min(1.0)),
            cdf,
        })))
    }

    fn check(&self, p: f64, name: &str, span: Span) -> Result<(), CssError> {
        if (0.0..=1.0 + 1e-12).contains(&p) {
            Ok(())
        } else {
            Err(ExecError::InvalidProbability {
                name: name.to_string(),
                line: span.line(self.line_map),
                value: p,
            }
            .into())
        }
    }

    fn unsupported(&self, name: &str, span: Span) -> CssError {
        ExecError::Unsupported {
            name: name.to_string(),
            line: span.line(self.line_map),
        }
        .into()
    }

    /// Keep the ops of the first iteration and check the labels of the others. Once the labels of
    /// the live segments repeat, every further period adds the same constraints, so whole periods
    /// are skipped and only the remainder is checked.
    fn repeat(&mut self, count: u64, body: &[ExtendedInstruction]) -> Result<(), CssError> {
        let start = self.ops.len();
        let mut seen = HashMap::new();
        let mut iteration = 0;
        while iteration < count {
            if let Some(first) = seen.insert(self.label_state(), iteration) {
                iteration = count - (count - iteration) % (iteration - first);
                seen.clear();
                continue;
            }
            let len = self.ops.len();
            self.block(body)?;
            if iteration > 0 {
                self.ops.truncate(len);
            }
            iteration += 1;
        }
        let body = self.ops.split_off(start);
        self.ops.push(Op::Repeat(Box::new(Repeat { count, body })));
        Ok(())
    }

    /// Hadamard parity and label class of each qubit's live segment, with classes numbered in
    /// order of first appearance and labels relative to the class's first qubit.
    fn label_state(&mut self) -> Vec<(usize, bool, bool)> {
        let mut classes = HashMap::new();
        let mut state = Vec::with_capacity(self.segment.len());
        for (&segment, &hadamard) in self.segment.iter().zip(&self.hadamard) {
            let (root, parity) = self.labels.find(segment);
            let next = classes.len();
            let &mut (class, first) = classes.entry(root).or_insert((next, parity));
            state.push((class, parity ^ first, hadamard));
        }
        state
    }

    fn reset(&mut self, q: usize, basis: Basis) {
        self.ops.push(Op::Reset { q, basis });
        self.new_segment(q);
    }

    fn pauli(&mut self, q: usize, x: bool, z: bool) {
        self.ops.push(Op::Pauli { q, x, z });
    }

    fn hadamards(&mut self, qubits: Vec<usize>) {
        for &q in &qubits {
            self.hadamard[q] ^= true;
        }
        self.ops.push(Op::H(qubits));
    }

    fn new_segment(&mut self, q: usize) {
        self.segment[q] = self.labels.add();
        self.hadamard[q] = false;
    }
}

fn log1m(p: f64) -> f64 {
    (-p).ln_1p()
}
