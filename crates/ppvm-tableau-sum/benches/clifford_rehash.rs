//! Benchmark the cost of a Clifford burst followed by a noise call.
//!
//! Used to estimate the overhead of moving from today's lazy-fingerprint
//! scheme to an eager scheme that re-fingerprints every entry on every
//! Clifford gate (the cost model implied by storing tableaus as HashMap
//! keys). Run once against the current code, then re-run after modifying
//! the Clifford gate macro in `gates/clifford.rs` to eagerly recompute
//! `entry_fingerprints` instead of clearing them.

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_runtime::traits::{Clifford, CliffordExtensions, Depolarizing};
use ppvm_tableau_sum::data::GeneralizedTableauSum;

type GTabSum = GeneralizedTableauSum<Byte8F64<2>, u128>;

/// Build a state with enough entries that per-entry rehash cost is visible.
fn build_state(n_qubits: usize) -> GTabSum {
    let mut tab: GTabSum = GeneralizedTableauSum::new_with_seed(n_qubits, 1e-12, 1e-8, 42);
    tab.entries.reserve(1024);
    tab.entry_fingerprints.reserve(1024);
    tab.h(0);
    tab.cnot(0, 1);
    tab.cnot(1, 2);
    tab.depolarize(0, 0.5);
    tab.depolarize(1, 0.5);
    tab.depolarize(2, 0.5);
    tab.depolarize(3, 0.5);
    tab
}

/// A long Clifford run between noise events. Every gate here mutates every
/// entry and invalidates its cached fingerprint under today's scheme.
fn apply_clifford_burst(tab: &mut GTabSum) {
    tab.h(4);
    tab.cnot(4, 5);
    tab.s(5);
    tab.cz(0, 5);
    tab.cnot(2, 3);
    tab.sqrt_x(6);
    tab.cnot(6, 7);
    tab.h(3);
    tab.cz(1, 4);
    tab.s_adj(2);
    tab.cnot(7, 0);
    tab.sqrt_y(5);
    tab.x(3);
    tab.y(6);
    tab.cnot(0, 4);
    tab.h(7);
}

fn clifford_rehash_benchmark(c: &mut Criterion) {
    let n_qubits = 8;
    let mut group = c.benchmark_group("clifford-rehash");
    group.bench_function("burst-then-noise", |b| {
        b.iter_batched(
            || build_state(n_qubits),
            |mut tab| {
                apply_clifford_burst(&mut tab);
                // Trailing noise forces fingerprint recomputation under
                // today's lazy scheme. An eager-rehash variant shifts
                // that cost from here into the burst.
                tab.depolarize(4, 0.3);
                // Return so `tab` is dropped outside the timing window.
                tab
            },
            criterion::BatchSize::SmallInput,
        );
    });
    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
        .sample_size(100);
    targets = clifford_rehash_benchmark
}
criterion_main!(benches);
