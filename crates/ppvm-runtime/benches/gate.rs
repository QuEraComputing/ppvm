use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_runtime::prelude::*;
use rayon::current_num_threads;

pub fn benchmark_suite<T: Config>(c: &mut Criterion, name: impl AsRef<str>) {
    let mut group = c.benchmark_group(name.as_ref());
    let n_qubits = 12;
    let mut state: PauliSum<T> = PauliSum::builder()
        .n_qubits(n_qubits)
        .capacity(1 << 20)
        .build();
    let mut term = PauliWord::new(n_qubits);
    term.set(0, Pauli::Z);
    term.set(1, Pauli::Z);
    state += (term.clone(), T::Coeff::from(1.0));
    for _ in 0..2 {
        for i in 0..n_qubits {
            state.rz(i, 1.1);
            state.ry(i, 2.1);
            state.rz(i, 1.1);
        }
        for i in 0..n_qubits {
            state.cnot(i, (i + 1) % n_qubits);
        }
    }
    println!("Initial state has {} terms", state.len());

    group.bench_function("x", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.x(0);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.bench_function("y", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.y(0);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.bench_function("z", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.z(0);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.bench_function("h", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.h(0);
            },
            criterion::BatchSize::SmallInput,
        );
    });
    group.bench_function("cnot", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.cnot(0, 1);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.bench_function("cz", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.cz(0, 1);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.bench_function("rx", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.rx(1, 0.5);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.bench_function("rz", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.rz(1, 0.5);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.bench_function("ry", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.ry(1, 0.5);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.bench_function("rxx", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.rxx(1, 2, 0.5);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.bench_function("ryy", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.ryy(1, 2, 0.5);
            },
            criterion::BatchSize::SmallInput,
        );
    });

    group.bench_function("rzz", |b| {
        b.iter_batched_ref(
            || state.clone(),
            |state| {
                state.rzz(1, 2, 0.5);
            },
            criterion::BatchSize::SmallInput,
        );
    });
    group.finish();
}

pub fn micro_benchmark(c: &mut Criterion) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build_global()
        .unwrap();
    println!("Using {} threads", current_num_threads());
    benchmark_suite::<config::gxhash::ByteF64<2>>(c, "ByteF64GxHashMap<2>");
    benchmark_suite::<config::fxhash::ByteF64<2>>(c, "ByteF64FxHashMap<2>");
    benchmark_suite::<config::dashmap::ByteFxHashF64<2>>(c, "ByteF64FxDashMap<2>");
    benchmark_suite::<config::dashmap::ByteGxHashF64<2>>(c, "ByteF64GxDashMap<2>");
    benchmark_suite::<config::indexmap::ByteFxHashF64<2>>(c, "ByteF64FxIndexMap<2>");
    benchmark_suite::<config::indexmap::ByteGxHashF64<2>>(c, "ByteF64GxIndexMap<2>");
}

criterion_group!(benches, micro_benchmark);
criterion_main!(benches);
