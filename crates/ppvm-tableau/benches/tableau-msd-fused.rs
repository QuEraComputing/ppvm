// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::{measure_all::LossyMeasureAll, prelude::*};

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

fn msd_func_fused<const MEASURE: bool>() -> (String, Tab) {
    let qubits_per_code_block = 17;
    let n_qubits = qubits_per_code_block * 5;
    debug_assert!(
        n_qubits < 8 * 11,
        "Make sure to update the bytes to match the qubit number"
    );
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);

    let qubit_addrs: Vec<usize> = (0..n_qubits).collect();
    let ql: Vec<&[usize]> = qubit_addrs.chunks_exact(qubits_per_code_block).collect();

    debug_assert_eq!(ql.len(), 5);

    for q in ql.iter() {
        let encoding_qubit = if q.len() == 7 { q[6] } else { q[7] };
        tab.h(encoding_qubit);
        tab.t(encoding_qubit);
        encode_fused(&mut tab, q);
    }

    // sqrt_x on blocks 0, 1, 4 — batched per block
    tab.sqrt_x_many(ql[0]);
    tab.sqrt_x_many(ql[1]);
    tab.sqrt_x_many(ql[4]);

    // ql[0] x ql[1]: pairs (0,17)...(16,33) — all in word 0
    tab.cz_block_pairs(0, 17, 17);

    // ql[2] x ql[3]: pairs (34,51)...(50,67)
    tab.cz_block_pairs(34, 17, 13); // (34,51)...(46,63) in word 0
    // (47,64)...(50,67): controls word 0 bits 47-50, targets word 1 bits 0-3
    tab.cz_block_pairs_cross_word(0, 47, 1, 0, 4);

    // sqrt_y on ql[0] and ql[3]
    tab.sqrt_y_many(ql[0]);
    tab.sqrt_y_many(ql[3]);

    // ql[0] x ql[2]: pairs (0,34)...(16,50) — all in word 0
    tab.cz_block_pairs(0, 34, 17);

    // ql[3] x ql[4]: (51,68)...(67,84)
    // (51,68)...(63,80): controls word 0 bits 51-63, targets word 1 bits 4-16
    tab.cz_block_pairs_cross_word(0, 51, 1, 4, 13);
    tab.cz_block_pairs(64, 17, 4); // (64,81)...(67,84) both in word 1

    tab.sqrt_x_dag_many(ql[0]);

    // ql[0] x ql[4]: (0,68)...(16,84)
    // controls word 0 bits 0-16, targets word 1 bits 4-20
    tab.cz_block_pairs_cross_word(0, 0, 1, 4, 17);

    // ql[1] x ql[3]: (17,51)...(33,67)
    tab.cz_block_pairs(17, 34, 13); // (17,51)...(29,63) in word 0
    // (30,64)...(33,67): controls word 0 bits 30-33, targets word 1 bits 0-3
    tab.cz_block_pairs_cross_word(0, 30, 1, 0, 4);

    // sqrt_x_dag on all blocks
    for block in ql.iter().take(5) {
        tab.sqrt_x_dag_many(block);
    }

    if MEASURE {
        let bit_string: String = tab
            .measure_all()
            .into_iter()
            .map(|outcome| if outcome.unwrap() { '1' } else { '0' })
            .collect();
        (bit_string, tab)
    } else {
        ("".to_owned(), tab)
    }
}

fn encode_fused(tab: &mut Tab, qubits: &[usize]) {
    if qubits.len() != 7 && qubits.len() != 17 {
        panic!("Unsupported number of qubits {}", qubits.len());
    }

    if qubits.len() == 7 {
        tab.sqrt_y_dag_many(&[
            qubits[0], qubits[1], qubits[2], qubits[3], qubits[4], qubits[5],
        ]);

        tab.cz_many(&[
            (qubits[1], qubits[2]),
            (qubits[3], qubits[4]),
            (qubits[5], qubits[6]),
        ]);

        tab.sqrt_y(qubits[6]);

        tab.cz_many(&[
            (qubits[0], qubits[3]),
            (qubits[2], qubits[5]),
            (qubits[4], qubits[6]),
        ]);

        tab.sqrt_y_many(&[qubits[2], qubits[3], qubits[4], qubits[5], qubits[6]]);

        tab.cz_many(&[
            (qubits[0], qubits[1]),
            (qubits[2], qubits[3]),
            (qubits[4], qubits[5]),
        ]);

        tab.sqrt_y_many(&[qubits[1], qubits[2], qubits[4]]);

        return;
    }

    // NOTE: len == 17 here
    tab.sqrt_y_many(&[
        qubits[0], qubits[1], qubits[2], qubits[3], qubits[4], qubits[5], qubits[6], qubits[8],
        qubits[9], qubits[10], qubits[11], qubits[12], qubits[13], qubits[14], qubits[15],
        qubits[16],
    ]);

    tab.cz_many(&[
        (qubits[1], qubits[3]),
        (qubits[7], qubits[10]),
        (qubits[12], qubits[14]),
        (qubits[13], qubits[16]),
    ]);

    tab.sqrt_y_dag_many(&[qubits[7], qubits[16]]);

    tab.cz_many(&[
        (qubits[4], qubits[7]),
        (qubits[8], qubits[10]),
        (qubits[11], qubits[14]),
        (qubits[15], qubits[16]),
    ]);

    tab.sqrt_y_dag_many(&[qubits[4], qubits[10], qubits[14], qubits[16]]);

    tab.cz_many(&[
        (qubits[2], qubits[4]),
        (qubits[6], qubits[8]),
        (qubits[7], qubits[9]),
        (qubits[10], qubits[13]),
        (qubits[14], qubits[16]),
    ]);

    tab.sqrt_y_many(&[
        qubits[3], qubits[6], qubits[9], qubits[10], qubits[12], qubits[13],
    ]);

    tab.cz_many(&[
        (qubits[0], qubits[2]),
        (qubits[3], qubits[6]),
        (qubits[5], qubits[8]),
        (qubits[10], qubits[12]),
        (qubits[11], qubits[13]),
    ]);

    tab.sqrt_y_many(&[
        qubits[1], qubits[2], qubits[3], qubits[4], qubits[6], qubits[7], qubits[8], qubits[9],
        qubits[11], qubits[12], qubits[14],
    ]);

    tab.cz_many(&[
        (qubits[0], qubits[1]),
        (qubits[2], qubits[3]),
        (qubits[4], qubits[5]),
        (qubits[6], qubits[7]),
        (qubits[8], qubits[9]),
        (qubits[12], qubits[15]),
    ]);

    tab.sqrt_y_dag_many(&[
        qubits[0], qubits[2], qubits[5], qubits[6], qubits[8], qubits[10], qubits[12],
    ]);
}

pub fn benchmark_suite_msd_fused(c: &mut Criterion, name: impl AsRef<str>) {
    let mut group = c.benchmark_group(name.as_ref());
    group.bench_function("msd-fused-0", |b| {
        b.iter_batched_ref(
            || {},
            |_| {
                msd_func_fused::<true>();
            },
            criterion::BatchSize::SmallInput,
        );
    });

    // propagate up to measurements
    let (_, tab) = msd_func_fused::<false>();
    group.bench_function("msd-fused-sample", |b| {
        b.iter_batched_ref(
            || {},
            |_| {
                let mut tab_sample = tab.fork(None);
                let bit_string: String = (0..85)
                    .map(|i| tab_sample.measure(i))
                    .map(|outcome| if outcome.unwrap() { '1' } else { '0' })
                    .collect();
                std::hint::black_box(bit_string);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

pub fn msd_fused_benchmarks(c: &mut Criterion) {
    benchmark_suite_msd_fused(c, "msd-fused");
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(50);
    targets = msd_fused_benchmarks
}
criterion_main!(benches);
