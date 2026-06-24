// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! General STIM-program benchmark: parse a `.stim` source once, then time the
//! tableau [`execute`] cost on a freshly built [`GeneralizedTableau`] per
//! iteration. Parsing stays outside the timed loop — the optimization target is
//! the tableau branching workload, not the parser.
//!
//! Lives in `ppvm-stim` (not `ppvm-tableau`) because the STIM parser only
//! exists here: `ppvm-stim` depends on `ppvm-tableau`, so the parser cannot be
//! reached from the tableau crate. This matches the sibling `tableau-msd-stim`
//! bench.
//!
//! The headline input is `cultivation_d5.stim`, the magic-state cultivation
//! circuit: it exercises the T-gate branching path (native `T`/`T_DAG`
//! rotations split the tableau into a superposition of branches), which is the
//! expensive, branchy workload we want to track. The remaining inputs form a
//! small representative basket (repetition code, GHZ, Bell pair, classically
//! controlled feedback) to keep the bench honest as a general-program runner.
//!
//! The tableau type is `GeneralizedTableau<ByteFxHashF64<8>, usize>`, matching
//! `tests/cultivation.rs` so correctness is identical to the e2e test.

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_stim::{execute, parse_extended};
use ppvm_tableau::prelude::*;
use stim_parser::prelude::{ExtendedInstruction, ExtendedProgram};

type Tab = GeneralizedTableau<ByteFxHashF64<8>, usize>;

/// Tableau qubit width needed to run `program`: the largest qubit index touched
/// by any instruction, plus one.
///
/// The tableau is fixed-width (set at construction; it does not auto-grow), and
/// the `ExtendedProgram` AST exposes no qubit-count accessor — only
/// [`ExtendedProgram::measurement_count`], which counts recorded bits, not
/// qubits. So we walk the instruction tree and take `max index + 1`, recursing
/// through `REPEAT` bodies. `rec[-k]` record-lookback targets are skipped
/// (they index the measurement record, not a qubit). A program with no qubit
/// targets needs zero qubits.
fn required_qubits(program: &ExtendedProgram) -> usize {
    fn max_in_slice(instructions: &[ExtendedInstruction]) -> Option<usize> {
        let mut max: Option<usize> = None;
        let mut bump = |q: usize| max = Some(max.map_or(q, |m: usize| m.max(q)));

        for instr in instructions {
            match instr {
                // Gate operands carry the qubit / `rec[-k]` distinction:
                // skip record lookbacks, count only `Target::Qubit`.
                ExtendedInstruction::Gate(op) => {
                    for t in &op.targets {
                        if let Some(q) = t.as_qubit() {
                            bump(q);
                        }
                    }
                }
                // Noise / Measure / Annotation operands are plain qubit indices.
                ExtendedInstruction::Noise(op) => {
                    for &q in &op.targets {
                        bump(q);
                    }
                }
                ExtendedInstruction::Measure(op) => {
                    for &q in &op.targets {
                        bump(q);
                    }
                }
                ExtendedInstruction::Annotation(op) => {
                    for &q in &op.targets {
                        bump(q);
                    }
                }
                ExtendedInstruction::T { targets, .. }
                | ExtendedInstruction::TDag { targets, .. }
                | ExtendedInstruction::Rotation { targets, .. }
                | ExtendedInstruction::U3 { targets, .. }
                | ExtendedInstruction::Loss { targets, .. } => {
                    for &q in targets {
                        bump(q);
                    }
                }
                ExtendedInstruction::CorrelatedLoss { targets, .. } => {
                    for &(a, b) in targets {
                        bump(a);
                        bump(b);
                    }
                }
                ExtendedInstruction::Mpp(op) => {
                    for product in &op.products {
                        for factor in product {
                            bump(factor.qubit);
                        }
                    }
                }
                // `MPad` records bits without touching a qubit.
                ExtendedInstruction::MPad { .. } => {}
                ExtendedInstruction::Repeat { body, .. } => {
                    if let Some(q) = max_in_slice(body) {
                        bump(q);
                    }
                }
            }
        }
        max
    }

    max_in_slice(&program.instructions).map_or(0, |m| m + 1)
}

/// Parse `src` once (panicking with context on failure), then benchmark
/// [`execute`] on a freshly constructed tableau per iteration. The fresh
/// tableau is built in the `iter_batched_ref` setup closure, so the timed
/// routine measures `execute` including the per-run tableau state — matching
/// the sibling `tableau-msd-stim` bench.
fn bench_stim_program(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    name: &str,
    src: &str,
) {
    let prog =
        parse_extended(src).unwrap_or_else(|e| panic!("parse_extended failed for {name}: {e}"));
    let n_qubits = required_qubits(&prog);

    group.bench_function(name, |b| {
        b.iter_batched_ref(
            || GeneralizedTableau::new(n_qubits, 1e-10),
            |tab: &mut Tab| execute(&prog, tab).expect("execute"),
            criterion::BatchSize::SmallInput,
        );
    });
}

pub fn stim_circuits_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("stim-circuits");

    // Headline: the branchy magic-state cultivation circuit (T/T_DAG, MPP).
    bench_stim_program(
        &mut group,
        "cultivation_d5",
        include_str!("../tests/data/cultivation_d5.stim"),
    );

    // Representative basket of runnable corpus programs.
    bench_stim_program(
        &mut group,
        "repetition_code_d3_r3",
        include_str!("../tests/data/repetition_code_d3_r3.stim"),
    );
    bench_stim_program(&mut group, "ghz", include_str!("../tests/data/ghz.stim"));
    bench_stim_program(
        &mut group,
        "bell_pair",
        include_str!("../tests/data/bell_pair.stim"),
    );
    bench_stim_program(
        &mut group,
        "feedback_cx",
        include_str!("../tests/data/feedback_cx.stim"),
    );

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50);
    targets = stim_circuits_benchmarks
}
criterion_main!(benches);
