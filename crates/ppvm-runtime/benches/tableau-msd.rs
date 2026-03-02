use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_runtime::config::indexmap::ByteFxHashF64;
use ppvm_runtime::prelude::*;
use ppvm_runtime::tableau::CliffordExtensions;
use rayon::current_num_threads;

type Tab = GeneralizedTableau<ByteFxHashF64<11>, u128>;

fn msd_func() -> String {
    // from Rafael:  https://www.notion.so/Simulating-85-qubit-MSD-circuit-using-stabilizer-rank-decomposition-and-pyzx-288f86eeff3c802fb262ef1cfa69dfae?source=copy_link#28df86eeff3c80bfa087ed15bcf49b77
    let qubits_per_code_block = 17;
    let n_qubits = qubits_per_code_block * 5;
    debug_assert!(
        n_qubits < 8 * 11,
        "Make sure to update the bytes to match the qubit number"
    );
    let mut tab: Tab = GeneralizedTableau::new(n_qubits, 1e-10);

    let qubit_addrs: Vec<usize> = (0..n_qubits).collect();

    // split qubits in 5 groups of n qubits
    let ql: Vec<&[usize]> = qubit_addrs.chunks_exact(qubits_per_code_block).collect(); //np.array_split(qubits, 5)

    debug_assert_eq!(ql.len(), 5);

    for q in ql.iter() {
        let encoding_qubit = if q.len() == 7 { q[6] } else { q[7] };
        tab.h(encoding_qubit);
        tab.t(encoding_qubit);
        encode(&mut tab, q);
    }

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

    for i in 0..5 {
        for q in ql[i] {
            tab.sqrt_x_adj(*q);
        }
    }

    let bit_string: String = (0..n_qubits)
        .map(|i| tab.measure(i))
        .map(|outcome| if outcome { '1' } else { '0' })
        .collect();

    bit_string
}

fn encode(tab: &mut Tab, qubits: &[usize]) {
    if qubits.len() != 7 && qubits.len() != 17 {
        panic!("Unsupported number of qubits {}", qubits.len());
    }

    if qubits.len() == 7 {
        for (idx, q) in qubits.iter().enumerate() {
            if idx == 6 {
                continue;
            }

            tab.sqrt_y_adj(*q);
        }

        tab.cz(qubits[1], qubits[2]);
        tab.cz(qubits[3], qubits[4]);
        tab.cz(qubits[5], qubits[6]);

        tab.sqrt_y(qubits[6]);

        tab.cz(qubits[0], qubits[3]);
        tab.cz(qubits[2], qubits[5]);
        tab.cz(qubits[4], qubits[6]);

        for (idx, q) in qubits.iter().enumerate() {
            if idx < 2 {
                continue;
            }
            tab.sqrt_y(*q);
        }

        tab.cz(qubits[0], qubits[1]);
        tab.cz(qubits[2], qubits[3]);
        tab.cz(qubits[4], qubits[5]);

        tab.sqrt_y(qubits[1]);
        tab.sqrt_y(qubits[2]);
        tab.sqrt_y(qubits[4]);

        return;
    }

    // NOTE: len == 17 here
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

pub fn benchmark_suite_msd(c: &mut Criterion, name: impl AsRef<str>) {
    let mut group = c.benchmark_group(name.as_ref());
    group.bench_function("msd-0", |b| {
        b.iter_batched_ref(
            || {},
            |_| {
                msd_func();
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.finish();
}

pub fn msd_benchmarks(c: &mut Criterion) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build_global()
        .unwrap();
    println!("Using {} threads", current_num_threads());
    benchmark_suite_msd(c, "msd");
}

criterion_group!(benches, msd_benchmarks);
criterion_main!(benches);
