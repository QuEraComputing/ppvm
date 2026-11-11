// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Benchmark and profiling harness: runs a Stim circuit's Clifford gates,
//! measurements and resets on `InverseTableau`, one shot at a time with a fresh
//! tableau per shot. Noise and annotations are skipped (noise sampling is the
//! circuit runner's job, not the tableau's).
//!
//! ```text
//! cargo run --release -p ppvm-tableau-2 --example run_stim -- <circuit.stim> [shots]
//! ```
//!
//! The clifft-bench surface code (d = 7, 7 rounds) is fetched by
//! `scripts/fetch_clifft_bench.sh` into `target/clifft-bench/`.

use std::time::Instant;

use ppvm_tableau_2::InverseTableau;
use ppvm_traits_2::Clifford;
use stim_parser::instructions::{GateName, MeasureName};
use stim_parser::prelude::{Instruction, Target};

/// One gate, or one measurement instruction over `targets[at..at + len]`
/// (kept together so they run as one `measure_many`).
#[derive(Clone, Copy)]
enum Op {
    One(GateName, usize),
    Two(GateName, usize, usize),
    Measure { at: usize, len: usize },
    Reset(usize),
    MeasureReset { at: usize, len: usize },
}

fn qubits(targets: &[Target]) -> Vec<usize> {
    targets
        .iter()
        .map(|t| t.as_qubit().expect("only qubit targets are supported"))
        .collect()
}

/// Flatten the program into ops, unrolling `REPEAT` blocks. Measurement
/// targets go to `targets`.
fn lower(instructions: &[Instruction], ops: &mut Vec<Op>, targets: &mut Vec<usize>) {
    for instruction in instructions {
        match instruction {
            Instruction::Gate(g) => {
                let qs = qubits(&g.targets);
                match g.name {
                    GateName::Reset | GateName::ResetZ => {
                        ops.extend(qs.iter().map(|&q| Op::Reset(q)))
                    }
                    GateName::CX | GateName::CZ | GateName::CY => {
                        ops.extend(qs.chunks(2).map(|p| Op::Two(g.name, p[0], p[1])))
                    }
                    name => ops.extend(qs.iter().map(|&q| Op::One(name, q))),
                }
            }
            Instruction::Measure(m) => {
                let (at, len) = (targets.len(), m.targets.len());
                targets.extend_from_slice(&m.targets);
                ops.push(match m.name {
                    MeasureName::M | MeasureName::MZ => Op::Measure { at, len },
                    MeasureName::MR => Op::MeasureReset { at, len },
                    other => panic!("unsupported measurement {other:?}"),
                });
            }
            Instruction::Repeat { count, body, .. } => {
                let mut block = Vec::new();
                lower(body, &mut block, targets);
                for _ in 0..*count {
                    ops.extend_from_slice(&block);
                }
            }
            Instruction::Noise(_) | Instruction::Annotation(_) => {}
            other => panic!("unsupported instruction {other:?}"),
        }
    }
}

fn apply_one(t: &mut InverseTableau, name: GateName, q: usize) {
    match name {
        GateName::X => t.x(q),
        GateName::Y => t.y(q),
        GateName::Z => t.z(q),
        GateName::H | GateName::HXZ => t.h(q),
        GateName::S | GateName::SqrtZ => t.s(q),
        GateName::SDag | GateName::SqrtZDag => t.s_dag(q),
        GateName::SqrtX => t.sqrt_x(q),
        GateName::SqrtXDag => t.sqrt_x_dag(q),
        GateName::SqrtY => t.sqrt_y(q),
        GateName::SqrtYDag => t.sqrt_y_dag(q),
        other => panic!("unsupported gate {other:?}"),
    }
}

fn apply_two(t: &mut InverseTableau, name: GateName, a: usize, b: usize) {
    match name {
        GateName::CX => t.cnot(a, b),
        GateName::CZ => t.cz(a, b),
        _ => t.cy(a, b),
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: run_stim <circuit.stim> [shots]");
    let shots: usize = args.next().map_or(1000, |s| s.parse().expect("shots"));
    let program = stim_parser::parse(&std::fs::read_to_string(&path).expect("read circuit"))
        .unwrap_or_else(|e| panic!("parse error: {e:?}"));
    let (mut ops, mut targets) = (Vec::new(), Vec::new());
    lower(&program.instructions, &mut ops, &mut targets);
    let n = 1 + ops
        .iter()
        .flat_map(|op| match *op {
            Op::Two(_, a, b) => vec![a, b],
            Op::One(_, q) | Op::Reset(q) => vec![q],
            Op::Measure { at, len } | Op::MeasureReset { at, len } => {
                targets[at..at + len].to_vec()
            }
        })
        .max()
        .expect("empty circuit");
    let mut outcomes = vec![false; targets.len()];

    let mut seed = 0x9e37_79b9_7f4a_7c15u64;
    let mut coin = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed & 1 == 1
    };
    let (mut parity, mut random) = (false, 0usize);
    let start = Instant::now();
    for _ in 0..shots {
        let mut t = InverseTableau::new(n);
        for &op in &ops {
            let mut flip = || {
                random += 1;
                coin()
            };
            match op {
                Op::One(name, q) => apply_one(&mut t, name, q),
                Op::Two(name, a, b) => apply_two(&mut t, name, a, b),
                Op::Measure { at, len } => {
                    let out = &mut outcomes[..len];
                    t.measure_many(&targets[at..at + len], &mut flip, out);
                    parity ^= out.iter().fold(false, |p, &m| p ^ m);
                }
                Op::Reset(q) => t.reset(q, &mut flip),
                Op::MeasureReset { at, len } => {
                    let out = &mut outcomes[..len];
                    t.measure_many(&targets[at..at + len], &mut flip, out);
                    for (&q, &m) in targets[at..at + len].iter().zip(out.iter()) {
                        parity ^= m;
                        if m {
                            t.x(q);
                        }
                    }
                }
            }
        }
    }
    let elapsed = start.elapsed().as_secs_f64();
    let count = |reset: bool| -> usize {
        ops.iter()
            .map(|op| match *op {
                Op::Measure { len, .. } | Op::MeasureReset { len, .. } => len,
                Op::Reset(_) => usize::from(reset),
                _ => 0,
            })
            .sum()
    };
    let measurements = count(false);
    println!(
        "{path}: {n} qubits, {} ops ({measurements} measurements) per shot, {shots} shots",
        ops.len()
    );
    println!(
        "{:.1} µs/shot, {:.0} shots/s, {:.1}% of measurements and resets random (parity {parity})",
        elapsed / shots as f64 * 1e6,
        shots as f64 / elapsed,
        100.0 * random as f64 / (shots * count(true)) as f64,
    );
}
