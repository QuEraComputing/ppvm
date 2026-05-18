use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;
use ppvm_tableau_sum::{data::GeneralizedTableauSum, sampler::Sampler};

type GTabSum = GeneralizedTableauSum<Byte8F64<2>, u128>;

fn encode(tab: &mut GTabSum, qubits: &[usize]) {
    if qubits.len() == 17 {
        for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16] {
            tab.sqrt_y(qubits[i]);
        }
        for [i, j] in [[1, 3], [7, 10], [12, 14], [13, 16]] {
            tab.cz(qubits[i], qubits[j]);
        }
        for i in [7, 16] {
            tab.sqrt_y_adj(qubits[i]);
        }
        for [i, j] in [[4, 7], [8, 10], [11, 14], [15, 16]] {
            tab.cz(qubits[i], qubits[j]);
        }
        for i in [4, 10, 14, 16] {
            tab.sqrt_y_adj(qubits[i]);
        }
        for [i, j] in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]] {
            tab.cz(qubits[i], qubits[j]);
        }
        for i in [3, 6, 9, 10, 12, 13] {
            tab.sqrt_y(qubits[i]);
        }
        for [i, j] in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]] {
            tab.cz(qubits[i], qubits[j]);
        }
        for i in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14] {
            tab.sqrt_y(qubits[i]);
        }
        for [i, j] in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]] {
            tab.cz(qubits[i], qubits[j]);
        }
        for i in [0, 2, 5, 6, 8, 10, 12] {
            tab.sqrt_y_adj(qubits[i]);
        }
    }
}

fn msd() -> Sampler<Byte8F64<2>, u128> {
    let n_qubits = 85;

    let mut tab: GTabSum = GeneralizedTableauSum::new(n_qubits, 1e-10, 1e-8);
    let qubit_addrs: Vec<usize> = (0..n_qubits).collect();
    let ql: Vec<&[usize]> = qubit_addrs.chunks_exact(17).collect();

    // Phase 1: Encoding (H + T + encode per block)
    for q in ql.iter() {
        let encoding_qubit = q[7];
        tab.h(encoding_qubit);
        tab.t(encoding_qubit);
        encode(&mut tab, q);
    }

    // Phase 2: Middle gates (sqrt_x, cz, sqrt_y, sqrt_x_adj layers)
    for i in [0, 1, 4] {
        for q in ql[i] {
            tab.sqrt_x(*q);
        }
    }
    for (control, target) in ql[0].iter().zip(ql[1]) {
        tab.cz(*control, *target);
    }
    for (control, target) in ql[2].iter().zip(ql[3]) {
        tab.cz(*control, *target);
    }
    for q in ql[0] {
        tab.sqrt_y(*q);
    }
    for q in ql[3] {
        tab.sqrt_y(*q);
    }
    for (control, target) in ql[0].iter().zip(ql[2]) {
        tab.cz(*control, *target);
    }
    for (control, target) in ql[3].iter().zip(ql[4]) {
        tab.cz(*control, *target);
    }
    for q in ql[0] {
        tab.sqrt_x_adj(*q);
    }
    for (control, target) in ql[0].iter().zip(ql[4]) {
        tab.cz(*control, *target);
    }
    for (control, target) in ql[1].iter().zip(ql[3]) {
        tab.cz(*control, *target);
    }
    for block in ql.iter().take(5) {
        for q in *block {
            tab.sqrt_x_adj(*q);
        }
    }

    tab.sampler()
}

pub fn benchmark_suite_msd(c: &mut Criterion, name: impl AsRef<str>) {
    let mut group = c.benchmark_group(name.as_ref());
    let mut sampler = msd();
    group.bench_function("msd-sampler", |b| {
        b.iter_batched_ref(
            || {},
            |_| sampler.sample(),
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

pub fn msd_benchmarks(c: &mut Criterion) {
    benchmark_suite_msd(c, "msd-sampler");
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(100);
    targets = msd_benchmarks
}
criterion_main!(benches);
