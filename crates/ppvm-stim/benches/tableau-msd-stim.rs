// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_stim::{execute, parse_extended};
use ppvm_tableau::prelude::*;
use stim_parser::extended::ExtendedProgram;

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

fn msd_stim_string() -> String {
    let qubits_per_code_block = 17;
    let n_qubits = qubits_per_code_block * 5;

    let qubit_addrs: Vec<usize> = (0..n_qubits).collect();
    let ql: Vec<&[usize]> = qubit_addrs.chunks_exact(qubits_per_code_block).collect();

    let mut lines: Vec<String> = Vec::new();

    for q in ql.iter() {
        let encoding_qubit = q[7];
        lines.push(format!("H {encoding_qubit}"));
        lines.push(format!("S[T] {encoding_qubit}"));
        encode_stim(&mut lines, q);
    }

    for i in [0, 1, 4] {
        lines.push(fmt_gate("SQRT_X", ql[i]));
    }

    lines.push(fmt_cz_pairs(ql[0], ql[1]));
    lines.push(fmt_cz_pairs(ql[2], ql[3]));

    lines.push(fmt_gate("SQRT_Y", ql[0]));
    lines.push(fmt_gate("SQRT_Y", ql[3]));

    lines.push(fmt_cz_pairs(ql[0], ql[2]));
    lines.push(fmt_cz_pairs(ql[3], ql[4]));

    lines.push(fmt_gate("SQRT_X_DAG", ql[0]));

    lines.push(fmt_cz_pairs(ql[0], ql[4]));
    lines.push(fmt_cz_pairs(ql[1], ql[3]));

    for q in ql.iter().take(5) {
        lines.push(fmt_gate("SQRT_X_DAG", q));
    }

    let all_qubits: Vec<usize> = (0..n_qubits).collect();
    lines.push(fmt_gate("M", &all_qubits));

    lines.join("\n")
}

fn encode_stim(lines: &mut Vec<String>, q: &[usize]) {
    assert_eq!(q.len(), 17);

    lines.push(fmt_gate_indices(
        "SQRT_Y",
        q,
        &[0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16],
    ));

    lines.push(fmt_cz_index_pairs(
        q,
        &[[1, 3], [7, 10], [12, 14], [13, 16]],
    ));
    lines.push(fmt_gate_indices("SQRT_Y_DAG", q, &[7, 16]));

    lines.push(fmt_cz_index_pairs(
        q,
        &[[4, 7], [8, 10], [11, 14], [15, 16]],
    ));
    lines.push(fmt_gate_indices("SQRT_Y_DAG", q, &[4, 10, 14, 16]));

    lines.push(fmt_cz_index_pairs(
        q,
        &[[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]],
    ));
    lines.push(fmt_gate_indices("SQRT_Y", q, &[3, 6, 9, 10, 12, 13]));

    lines.push(fmt_cz_index_pairs(
        q,
        &[[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]],
    ));
    lines.push(fmt_gate_indices(
        "SQRT_Y",
        q,
        &[1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14],
    ));

    lines.push(fmt_cz_index_pairs(
        q,
        &[[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]],
    ));
    lines.push(fmt_gate_indices("SQRT_Y_DAG", q, &[0, 2, 5, 6, 8, 10, 12]));
}

fn fmt_gate(gate: &str, qubits: &[usize]) -> String {
    let targets: Vec<String> = qubits.iter().map(|q| q.to_string()).collect();
    format!("{gate} {}", targets.join(" "))
}

fn fmt_gate_indices(gate: &str, q: &[usize], indices: &[usize]) -> String {
    fmt_gate(gate, &indices.iter().map(|&i| q[i]).collect::<Vec<_>>())
}

fn fmt_cz_pairs(controls: &[usize], targets: &[usize]) -> String {
    let pairs: Vec<String> = controls
        .iter()
        .zip(targets)
        .map(|(c, t)| format!("{c} {t}"))
        .collect();
    format!("CZ {}", pairs.join(" "))
}

fn fmt_cz_index_pairs(q: &[usize], pairs: &[[usize; 2]]) -> String {
    let targets: Vec<String> = pairs
        .iter()
        .map(|[i, j]| format!("{} {}", q[*i], q[*j]))
        .collect();
    format!("CZ {}", targets.join(" "))
}

fn msd_stim_func(prog: &ExtendedProgram) {
    let n_qubits = 17 * 5;
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);
    execute(prog, &mut tab).expect("execute");
}

pub fn benchmark_suite_msd_stim(c: &mut Criterion, name: impl AsRef<str>) {
    let circuit = msd_stim_string();
    let prog = parse_extended(&circuit).expect("parse_extended");

    let mut group = c.benchmark_group(name.as_ref());
    group.bench_function("msd-stim-0", |b| {
        b.iter_batched_ref(
            || {},
            |_| msd_stim_func(&prog),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

pub fn msd_stim_benchmarks(c: &mut Criterion) {
    benchmark_suite_msd_stim(c, "msd-stim");
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50);
    targets = msd_stim_benchmarks
}
criterion_main!(benches);
