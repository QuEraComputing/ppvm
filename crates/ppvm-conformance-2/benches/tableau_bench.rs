// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Same-build tableau perf gate: OLD `ppvm-tableau` vs NEW `ppvm-tableau-2`,
//! both engines in **one** binary so the new/old ratio cancels the code-layout
//! bias (Mytkowicz et al., "Producing Wrong Data Without Doing Anything
//! Obviously Wrong") that makes cross-build absolutes unreliable.
//!
//! Follows the pattern of `benches/pauli_sum_integration.rs`.
//!
//! # What is the gate and what is diagnostic
//!
//! The HEADLINE metric is the group `tableau-integration`: whole real workloads
//! (the 85-qubit MSD circuit, naive and fused; the rot2 brickwork; the fused
//! T-gate circuit; the CNOT-chain scaling sweep; the 85-qubit measure-all sweep;
//! the 4000-shot noisy sampler; the branch-coalesce regimes), each reported as a
//! TOTAL wall-clock ratio. A deep circuit is the only thing that can see
//! cumulative per-gate costs — allocation churn, a dropped double-buffer, a batch
//! trait silently taking its loop default — because a tight one-gate `iter` loop
//! lets the allocator recycle one warm page.
//!
//! `tableau-micro` is **diagnostic only**: it attributes a movement in the
//! headline number to a specific kernel. A regression is never read off it.
//!
//! # Fair configuration (apples to apples)
//!
//! Every pair holds the algebraic config identical on both sides — same storage
//! width, same coefficient type, same (`FxHash`) hasher. The matched aliases live
//! in `ppvm_conformance_2::tableau`:
//!
//! * `OldWide`/`NewWide` — `[usize; 2]` storage, `u128` index (MSD, fused-T,
//!   branch-coalesce, scaling);
//! * `OldNarrow`/`NewNarrow` — `[u8; 8]` storage, `usize` index (rot2 brickwork,
//!   noisy sampler).
//!
//! Both sides also run the *same* workload source: every circuit below is a
//! single generic function over the `Driver` trait, instantiated twice.
//!
//! # Scaling-sweep note
//!
//! The old `tableau-scaling-{96,128}` benches used a `usize` index, where
//! `compute_decomposition`'s `1 << i` overflows for `i ≥ 64` (it panics in debug
//! and silently masks the shift in release). Both engines reproduce that
//! verbatim, so the sweep here is run on the `u128` index at every `n` — the same
//! config on both sides, and the only one that is well-defined above 64 qubits.
//!
//! Run:
//! ```bash
//! cargo bench -p ppvm-conformance-2 --bench tableau_bench
//! ```

use std::hint::black_box;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

use ppvm_conformance_2::tableau::*;

/// Longer than criterion's default so the CIs are tight enough to read a ratio.
fn integration_config() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(5))
        .sample_size(20)
}

// ===========================================================================
// 1. MSD-85q — naive and fused (integration baselines #1 and #2)
// ===========================================================================

fn bench_msd(c: &mut Criterion) {
    let mut g = c.benchmark_group("tableau-integration/msd-85q");
    g.warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(6))
        .sample_size(20);

    g.bench_function("naive/old", |b| {
        b.iter(|| black_box(msd_bitstring::<OldWide>(None)))
    });
    g.bench_function("naive/new", |b| {
        b.iter(|| black_box(msd_bitstring::<NewWide>(None)))
    });
    g.bench_function("fused/old", |b| {
        b.iter(|| black_box(msd_bitstring_fused::<OldWide>(None)))
    });
    g.bench_function("fused/new", |b| {
        b.iter(|| black_box(msd_bitstring_fused::<NewWide>(None)))
    });
    g.finish();
}

// ===========================================================================
// 2. rot2 brickwork (integration baseline #3)
// ===========================================================================

fn bench_rot2(c: &mut Criterion) {
    let mut g = c.benchmark_group("tableau-integration/rot2-brickwork");
    g.warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(4))
        .sample_size(20);

    for &(n, layers) in &[(8usize, 4usize), (10, 4), (12, 3)] {
        let m = rot2_brickwork::<NewNarrow>(n, layers).n_coeffs();
        g.bench_with_input(
            BenchmarkId::new("old", format!("n{n}_l{layers}_m{m}")),
            &(n, layers),
            |b, &(n, layers)| {
                b.iter(|| black_box(rot2_brickwork::<OldNarrow>(n, layers).n_coeffs()))
            },
        );
        g.bench_with_input(
            BenchmarkId::new("new", format!("n{n}_l{layers}_m{m}")),
            &(n, layers),
            |b, &(n, layers)| {
                b.iter(|| black_box(rot2_brickwork::<NewNarrow>(n, layers).n_coeffs()))
            },
        );
    }
    g.finish();
}

// ===========================================================================
// 3. fused T-gate circuit (integration baseline #4)
// ===========================================================================

fn bench_fused_tgate(c: &mut Criterion) {
    let mut g = c.benchmark_group("tableau-integration/fused-tgate");
    g.warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(6))
        .sample_size(20);

    for n_tgates in [8usize, 12, 16] {
        let old_setup: OldWide = fused_tgate_setup(n_tgates);
        let new_setup: NewWide = fused_tgate_setup(n_tgates);
        g.bench_with_input(
            BenchmarkId::new("old", format!("{n_tgates}t-85q")),
            &n_tgates,
            |b, &k| {
                b.iter_batched_ref(
                    || old_setup.fork(Some(0)),
                    |tab| black_box(fused_tgate_body(tab, k)),
                    criterion::BatchSize::LargeInput,
                )
            },
        );
        g.bench_with_input(
            BenchmarkId::new("new", format!("{n_tgates}t-85q")),
            &n_tgates,
            |b, &k| {
                b.iter_batched_ref(
                    || new_setup.fork(Some(0)),
                    |tab| black_box(fused_tgate_body(tab, k)),
                    criterion::BatchSize::LargeInput,
                )
            },
        );
    }
    g.finish();
}

// ===========================================================================
// 4. CNOT-chain scaling sweep + the n-qubit measurement sweep (baseline #5)
// ===========================================================================

fn bench_scaling(c: &mut Criterion) {
    let mut g = c.benchmark_group("tableau-integration/scaling");
    g.warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(4))
        .sample_size(30);

    for n in [32usize, 64, 96, 128] {
        let old_base: OldWide = Driver::new_seeded(n, 1e-10, 1);
        let new_base: NewWide = Driver::new_seeded(n, 1e-10, 1);
        g.bench_with_input(BenchmarkId::new("gates/old", n), &n, |b, _| {
            b.iter_batched_ref(
                || old_base.fork(None),
                |tab| black_box(scaling_circuit(tab)),
                criterion::BatchSize::SmallInput,
            )
        });
        g.bench_with_input(BenchmarkId::new("gates/new", n), &n, |b, _| {
            b.iter_batched_ref(
                || new_base.fork(None),
                |tab| black_box(scaling_circuit(tab)),
                criterion::BatchSize::SmallInput,
            )
        });

        // The whole n-qubit measurement sweep on a prepared state.
        let mut old_prepared: OldWide = Driver::new_seeded(n, 1e-10, 1);
        let mut new_prepared: NewWide = Driver::new_seeded(n, 1e-10, 1);
        scaling_prepare(&mut old_prepared);
        scaling_prepare(&mut new_prepared);
        g.bench_with_input(BenchmarkId::new("measure-sweep/old", n), &n, |b, _| {
            b.iter_batched_ref(
                || old_prepared.fork(None),
                |tab| {
                    for i in 0..n {
                        black_box(tab.measure(i));
                    }
                },
                criterion::BatchSize::SmallInput,
            )
        });
        g.bench_with_input(BenchmarkId::new("measure-sweep/new", n), &n, |b, _| {
            b.iter_batched_ref(
                || new_prepared.fork(None),
                |tab| {
                    for i in 0..n {
                        black_box(tab.measure(i));
                    }
                },
                criterion::BatchSize::SmallInput,
            )
        });
    }
    g.finish();
}

// ===========================================================================
// 5. measure_all / measure_many on the prepared MSD state (baseline #6)
// ===========================================================================

fn bench_measure_all(c: &mut Criterion) {
    let mut g = c.benchmark_group("tableau-integration/measure-all-msd");
    g.warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(4))
        .sample_size(30);

    let old_state: OldWide = msd_state(Some(3));
    let new_state: NewWide = msd_state(Some(3));
    let all: Vec<usize> = (0..MSD_QUBITS).collect();

    g.bench_function("measure_all/old", |b| {
        b.iter_batched_ref(
            || old_state.fork(None),
            |tab| black_box(tab.measure_all()),
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("measure_all/new", |b| {
        b.iter_batched_ref(
            || new_state.fork(None),
            |tab| black_box(tab.measure_all()),
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("measure_many/old", |b| {
        b.iter_batched_ref(
            || old_state.fork(None),
            |tab| black_box(tab.measure_many(&all)),
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("measure_many/new", |b| {
        b.iter_batched_ref(
            || new_state.fork(None),
            |tab| black_box(tab.measure_many(&all)),
            criterion::BatchSize::LargeInput,
        )
    });
    // The naive per-qubit loop: the gap to `measure_all` IS the scratch-reuse win.
    g.bench_function("measure_loop/old", |b| {
        b.iter_batched_ref(
            || old_state.fork(None),
            |tab| {
                for i in 0..MSD_QUBITS {
                    black_box(tab.measure(i));
                }
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("measure_loop/new", |b| {
        b.iter_batched_ref(
            || new_state.fork(None),
            |tab| {
                for i in 0..MSD_QUBITS {
                    black_box(tab.measure(i));
                }
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.finish();
}

// ===========================================================================
// 6. Noisy-Clifford shot averaging (integration baseline #7)
// ===========================================================================

fn bench_noisy_shots(c: &mut Criterion) {
    let mut g = c.benchmark_group("tableau-integration/noisy-shots");
    g.warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(4))
        .sample_size(20);

    const SHOTS: u64 = 4000;
    g.throughput(Throughput::Elements(SHOTS));
    g.bench_function("old", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for shot in 0..SHOTS {
                acc += noisy_shot::<OldNarrow>(shot);
            }
            black_box(acc)
        })
    });
    g.bench_function("new", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for shot in 0..SHOTS {
                acc += noisy_shot::<NewNarrow>(shot);
            }
            black_box(acc)
        })
    });
    g.finish();
}

// ===========================================================================
// 7. Branch-coalesce scaling (integration baseline #8)
// ===========================================================================

fn bench_branch_coalesce(c: &mut Criterion) {
    let mut g = c.benchmark_group("tableau-integration/branch-coalesce");
    g.warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_secs(3))
        .sample_size(20);

    // j straddles the packed-path cutoff (m ≤ 65535): 2^16 = 65536 falls back to
    // the generic `(I, u32)` sort, so a slope change there is the packing's
    // contribution.
    for j in [2usize, 5, 8, 11, 14, 16] {
        let m = 1u64 << j;
        g.throughput(Throughput::Elements(m));

        // doubling: the benched T flips a FRESH index bit ⇒ output 2m, no merges.
        let mut old_d: OldWide = branch_grow(j);
        let mut new_d: NewWide = branch_grow(j);
        old_d.h(j);
        new_d.h(j);
        g.bench_with_input(BenchmarkId::new("doubling/old", m), &j, |b, &j| {
            b.iter_batched_ref(
                || old_d.fork(None),
                |tab| tab.t(j),
                criterion::BatchSize::LargeInput,
            )
        });
        g.bench_with_input(BenchmarkId::new("doubling/new", m), &j, |b, &j| {
            b.iter_batched_ref(
                || new_d.fork(None),
                |tab| tab.t(j),
                criterion::BatchSize::LargeInput,
            )
        });

        // merge: the benched T reuses an already-branched bit ⇒ output m, all merges.
        let old_m: OldWide = branch_grow(j);
        let new_m: NewWide = branch_grow(j);
        g.bench_with_input(BenchmarkId::new("merge/old", m), &j, |b, _| {
            b.iter_batched_ref(
                || old_m.fork(None),
                |tab| tab.t(0),
                criterion::BatchSize::LargeInput,
            )
        });
        g.bench_with_input(BenchmarkId::new("merge/new", m), &j, |b, _| {
            b.iter_batched_ref(
                || new_m.fork(None),
                |tab| tab.t(0),
                criterion::BatchSize::LargeInput,
            )
        });
    }
    g.finish();
}

// ===========================================================================
// 8. Diagnostic microbenches (attribution only — never the gate)
// ===========================================================================

/// An 85-qubit state whose `measure(0)` takes the **case-a** branch, with a
/// support of exactly `2^n_branch` amplitudes: `h(0)` makes `Z₀` anticommute
/// with a stabilizer, and the `h;t` pairs on qubits `1..=n_branch` grow the
/// amplitude vector without touching the `Z₀` dichotomy.
fn case_a_state<D: Driver>(n_branch: usize) -> D {
    let mut t: D = Driver::new_seeded(85, 1e-10, 7);
    t.h(0);
    for q in 1..=n_branch {
        t.h(q);
        t.t(q);
    }
    t
}

fn bench_micro(c: &mut Criterion) {
    let mut g = c.benchmark_group("tableau-micro");
    g.warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_secs(2))
        .sample_size(50);

    const N: usize = 85;
    let old_base: OldWide = Driver::new_seeded(N, 1e-10, 1);
    let new_base: NewWide = Driver::new_seeded(N, 1e-10, 1);
    let block: Vec<usize> = (0..17).collect();

    macro_rules! micro {
        ($name:literal, $old:expr, $new:expr) => {
            g.bench_function(concat!($name, "/old"), |b| {
                b.iter_batched_ref(
                    || old_base.fork(None),
                    $old,
                    criterion::BatchSize::SmallInput,
                )
            });
            g.bench_function(concat!($name, "/new"), |b| {
                b.iter_batched_ref(
                    || new_base.fork(None),
                    $new,
                    criterion::BatchSize::SmallInput,
                )
            });
        };
    }

    micro!("h", |t: &mut OldWide| t.h(0), |t: &mut NewWide| t.h(0));
    micro!("s", |t: &mut OldWide| t.s(0), |t: &mut NewWide| t.s(0));
    micro!(
        "sqrt_y",
        |t: &mut OldWide| t.sqrt_y(0),
        |t: &mut NewWide| t.sqrt_y(0)
    );
    micro!(
        "cnot",
        |t: &mut OldWide| t.cnot(0, 70),
        |t: &mut NewWide| t.cnot(0, 70)
    );
    micro!("cz", |t: &mut OldWide| t.cz(0, 70), |t: &mut NewWide| t
        .cz(0, 70));
    micro!("t-gate", |t: &mut OldWide| t.t(0), |t: &mut NewWide| t.t(0));
    micro!(
        "measure",
        |t: &mut OldWide| {
            black_box(t.measure(0));
        },
        |t: &mut NewWide| {
            black_box(t.measure(0));
        }
    );
    // The fused batch kernels vs their per-qubit expansions: this pair is what
    // detects a batch trait that silently took its loop default.
    micro!(
        "sqrt_y_many17",
        |t: &mut OldWide| t.sqrt_y_many(&block),
        |t: &mut NewWide| t.sqrt_y_many(&block)
    );
    micro!(
        "sqrt_y_loop17",
        |t: &mut OldWide| {
            for &q in &block {
                t.sqrt_y(q)
            }
        },
        |t: &mut NewWide| {
            for &q in &block {
                t.sqrt_y(q)
            }
        }
    );
    micro!(
        "cz_block17",
        |t: &mut OldWide| t.cz_block(0, 17, 17),
        |t: &mut NewWide| t.cz_block(0, 17, 17)
    );
    micro!(
        "cz_loop17",
        |t: &mut OldWide| {
            for i in 0..17 {
                t.cz(i, i + 17)
            }
        },
        |t: &mut NewWide| {
            for i in 0..17 {
                t.cz(i, i + 17)
            }
        }
    );

    // ATTRIBUTION for the `measure_loop` ratio. OLD builds a fresh
    // `MeasureScratch` on the stack per `measure` call; NEW keeps one per
    // tableau, so it constructs (and heap-allocates) at most ONE for a whole
    // sweep and reuses it. The new scratch additionally owns the five
    // sort-merge working Vecs (`by_idx`/`shifted`/`a`/`bt`/`merged`), so it is
    // a bigger object to build.
    //
    // This pair therefore measures 85 construct+drop pairs, i.e. what OLD pays
    // per 85-qubit sweep and what NEW pays per 85 *tableaux*. Its ~3.7 ratio is
    // NOT a per-sweep cost on the new side and must not be read as one — the
    // real workload is `msd_sweep_loop` (0.79). It is kept only to bound the
    // one-off construction term visible in `msd_measure_single`.
    g.bench_function("scratch_new_x85/old", |b| {
        b.iter(|| {
            for _ in 0..85 {
                black_box(ppvm_tableau::measure::MeasureScratch::<u128, f64>::new());
            }
        })
    });
    g.bench_function("scratch_new_x85/new", |b| {
        b.iter(|| {
            for _ in 0..85 {
                black_box(ppvm_tableau_2::MeasureScratch::<u128>::new());
            }
        })
    });

    // ATTRIBUTION, second variable: the same single `measure` and a
    // `z_expectation` on the BRANCHY 85-qubit MSD state. `z_expectation` runs
    // `compute_decomposition` (the O(n²) frame walk, i.e. `Row::mul_assign` over
    // the anticommuting generators) plus the overlap and NOTHING else — no
    // scratch, no allocation, no projection, no RNG draw. So if the per-qubit
    // `measure` ratio and the `z_expectation` ratio move together, the cost is
    // in the frame walk, not in the scratch.
    let old_msd: OldWide = msd_state(Some(3));
    let new_msd: NewWide = msd_state(Some(3));
    g.bench_function("msd_measure_single/old", |b| {
        b.iter_batched_ref(
            || old_msd.fork(None),
            |t| black_box(t.measure(0)),
            criterion::BatchSize::SmallInput,
        )
    });
    g.bench_function("msd_measure_single/new", |b| {
        b.iter_batched_ref(
            || new_msd.fork(None),
            |t| black_box(t.measure(0)),
            criterion::BatchSize::SmallInput,
        )
    });
    g.bench_function("msd_z_expectation/old", |b| {
        b.iter(|| black_box(old_msd.z_expectation(0)))
    });
    g.bench_function("msd_z_expectation/new", |b| {
        b.iter(|| black_box(new_msd.z_expectation(0)))
    });

    // ATTRIBUTION, case-a fixed overhead vs per-element work. The `measure`
    // micro above is CASE B (a fresh frame has `Z₀` as a stabilizer), so it says
    // nothing about the sort-merge branch. These two isolate case a at two
    // support sizes on an otherwise identical frame:
    //
    //   * `case_a_m1`  — one `h(0)`, so `Z₀` anticommutes but the support is a
    //     single amplitude. Everything measured is FIXED overhead: the five
    //     working-Vec allocations, the scratch, the two length-1 sorts, and the
    //     `O(n)` frame projection.
    //   * `case_a_m32` — five `h;t` pairs first, so the same kernel runs over 32
    //     amplitudes. `(case_a_m32 − case_a_m1)` is the per-element work.
    //
    // A ratio that is high on `m1` and flat on `(m32 − m1)` is a fixed-overhead
    // regression; the reverse is an inner-loop regression.
    //
    // MEASURED (this machine, same build, after the `Row::site_probe` hoist):
    // m1 = 311/283 ns (old/new), m32 = 697/694. So `(m32 − m1)` is 386 old vs
    // 411 new — the per-element sort-merge inner loop is at parity-ish — while
    // the FIXED case-a cost now favours new by ~28 ns, which is the frame
    // half (`frame_project`, 0.64) paying for itself. Allocation COUNTS are
    // identical (6 per case-a measurement on both); `MeasureScratch::new` is
    // 3 ns (`scratch_new_x85`, and new constructs at most ONE per tableau), and
    // moving the five working Vecs into the scratch is a strict win wherever
    // the scratch is reused (`msd_sweep_*`, 0.79).
    let old_a1: OldWide = case_a_state(1);
    let new_a1: NewWide = case_a_state(1);
    let old_a32: OldWide = case_a_state(5);
    let new_a32: NewWide = case_a_state(5);
    g.bench_function("case_a_m1/old", |b| {
        b.iter_batched_ref(
            || old_a1.fork(None),
            |t| black_box(t.measure(0)),
            criterion::BatchSize::SmallInput,
        )
    });
    g.bench_function("case_a_m1/new", |b| {
        b.iter_batched_ref(
            || new_a1.fork(None),
            |t| black_box(t.measure(0)),
            criterion::BatchSize::SmallInput,
        )
    });
    g.bench_function("case_a_m32/old", |b| {
        b.iter_batched_ref(
            || old_a32.fork(None),
            |t| black_box(t.measure(0)),
            criterion::BatchSize::SmallInput,
        )
    });
    g.bench_function("case_a_m32/new", |b| {
        b.iter_batched_ref(
            || new_a32.fork(None),
            |t| black_box(t.measure(0)),
            criterion::BatchSize::SmallInput,
        )
    });

    // ATTRIBUTION, the case-a FRAME projection on its own. `Measure` on the bare
    // `Tableau` is exactly `find_z_anticommuting_stabilizer` + one RNG bool +
    // `update_tableau_according_to_outcome` — the `O(n)` row sweep case a runs
    // and case b does not. No amplitudes, no scratch, no allocation. This is the
    // one case-a step `msd_z_expectation` (which stops after
    // `compute_decomposition`) cannot see.
    g.bench_function("frame_project/old", |b| {
        b.iter_batched_ref(
            || old_a1.tableau.clone(),
            |t| black_box(ppvm_traits::Measure::measure(t, 0)),
            criterion::BatchSize::SmallInput,
        )
    });
    g.bench_function("frame_project/new", |b| {
        b.iter_batched_ref(
            || new_a1.tableau.clone(),
            |t| black_box(ppvm_traits_2::Measure::measure(t, 0)),
            criterion::BatchSize::SmallInput,
        )
    });

    // ATTRIBUTION, third variable: NEW's `measure_all` vs `measure_many(all)`
    // vs `measure_many_with_scratch(all, shared)` — all three sweep the same 85
    // qubits with ONE scratch, so any gap between them is not algorithmic.
    let all85: Vec<usize> = (0..MSD_QUBITS).collect();
    g.bench_function("msd_sweep_all/new", |b| {
        b.iter_batched_ref(
            || new_msd.fork(None),
            |t| black_box(ppvm_tableau_2::GeneralizedTableau::measure_all(t)),
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("msd_sweep_many/new", |b| {
        b.iter_batched_ref(
            || new_msd.fork(None),
            |t| black_box(ppvm_traits_2::Measure::measure_many(t, &all85)),
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("msd_sweep_many_scratch/new", |b| {
        b.iter_batched_ref(
            || {
                (
                    new_msd.fork(None),
                    ppvm_tableau_2::MeasureScratch::<u128>::new(),
                )
            },
            |(t, sc)| black_box(t.measure_many_with_scratch(&all85, sc)),
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("msd_sweep_loop/old", |b| {
        b.iter_batched_ref(
            || old_msd.fork(None),
            |t| {
                for i in 0..MSD_QUBITS {
                    black_box(t.measure(i));
                }
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("msd_sweep_loop/new", |b| {
        b.iter_batched_ref(
            || new_msd.fork(None),
            |t| {
                for i in 0..MSD_QUBITS {
                    black_box(t.measure(i));
                }
            },
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("msd_sweep_all/old", |b| {
        b.iter_batched_ref(
            || old_msd.fork(None),
            |t| black_box(ppvm_tableau::measure_all::LossyMeasureAll::measure_all(t)),
            criterion::BatchSize::LargeInput,
        )
    });
    g.bench_function("msd_sweep_all_scratch/new", |b| {
        b.iter_batched_ref(
            || {
                (
                    new_msd.fork(None),
                    ppvm_tableau_2::MeasureScratch::<u128>::new(),
                )
            },
            |(t, sc)| black_box(t.measure_all_with_scratch(sc)),
            criterion::BatchSize::LargeInput,
        )
    });

    // Construction + fork: every bench's setup and the sampler's hot path.
    g.bench_function("construct/old", |b| {
        b.iter(|| black_box(<OldWide as Driver>::new_seeded(N, 1e-10, 1).n_qubits()))
    });
    g.bench_function("construct/new", |b| {
        b.iter(|| black_box(<NewWide as Driver>::new_seeded(N, 1e-10, 1).n_qubits()))
    });
    g.bench_function("fork/old", |b| {
        b.iter(|| black_box(old_base.fork(Some(0)).n_qubits()))
    });
    g.bench_function("fork/new", |b| {
        b.iter(|| black_box(new_base.fork(Some(0)).n_qubits()))
    });

    g.finish();
}

// ===========================================================================
// 9. ATTRIBUTION for `scaling/measure-sweep` (a resolved regression — kept as
//    the standing guard on `compute_decomposition`)
// ===========================================================================

/// This group holds ONE variable at a time to find out which half of a
/// measurement a `scaling/measure-sweep` movement lives in, rather than
/// guessing. It is what found the ~1.11–1.15 regression this sweep used to
/// carry, and it stays as the regression guard for the kernel it blamed.
///
/// The workload is `n` per-qubit `measure` calls on the CNOT-chain state, whose
/// amplitude support is a constant 4–8 branches. A measurement is two separable
/// halves:
///
/// * the FRAME half — `compute_decomposition` (`O(n)`) plus, in case a, the
///   `O(n)` row sweep of `update_tableau_according_to_outcome`;
/// * the AMPLITUDE half — the case-b in-place overlap or the case-a sort-merge.
///
/// `frame_sweep` runs the identical `n`-qubit sweep on the BARE `Tableau`, i.e.
/// the frame half with the amplitude half deleted on both engines. Since the
/// support is tiny and `n`-independent while the ratio is flat in `n`, the frame
/// half is the only term that can carry a proportional cost — and this measures
/// it directly instead of inferring it.
///
/// Reading it: if `frame_sweep` ≈ the full `measure-sweep` ratio, a movement is
/// in the frame; if `frame_sweep` ≈ 1.0, it is in the amplitude path and the
/// frame is exonerated.
///
/// # HISTORY — the regression this group found, and its fix
///
/// The sweep once ran at 1.11–1.15 at every `n`. The single-variable A/B here
/// exonerated the row projection (`frame_sweep` 1.002) and the per-call scratch
/// handling (`many_sweep` still 1.16 with ONE scratch for the whole sweep), and
/// pinned `decomp_only` — `compute_decomposition` alone — as the cause: at
/// `n = 32` it was 2.33 µs of the 2.54 µs sweep, i.e. 92 % of it.
///
/// A `samply` profile then put ~81 % of the sweep's self time inside that one
/// function on BOTH engines, and the instruction histogram put ~40 % of it in
/// the two anticommutation scans. The scans were re-indexed through `bitvec`'s
/// per-bit addressing once per row; hoisting the site to a `(word, x-probe,
/// z-probe)` triple (`Row::site_probe`) turns the per-row test into two `AND`s
/// and an `XOR` on raw words, cutting the loop body from ~12 instructions to
/// ~7. The same hoist went into `find_z_anticommuting_stabilizer`,
/// `get_deterministic_outcome` and `update_tableau_according_to_outcome`.
///
/// | probe | n=32 | n=128 |
/// |:--|--:|--:|
/// | `frame_sweep` (row projection) | 0.97 | 0.98 |
/// | `many_sweep` (one scratch for the whole sweep) | 0.89 | 0.93 |
/// | `decomp_only` (`compute_decomposition` alone) | 0.82 | 0.86 |
/// | full `measure-sweep` | 0.91 | 0.90 |
///
/// The change is codegen only — same visit order, same multiplication order,
/// same accumulated phase — so the differential suite pins that nothing moved.
///
/// `decomp_sweep` (the `z_expectation` sweep) is kept as context but is NOT the
/// attribution probe: `z_expectation` builds a case-a `FxHashMap` over the whole
/// support, which swamps the decomposition term — and its absolute numbers moved
/// 42.9 → 47.3 µs on the OLD side across a relink with no change to that code,
/// a textbook code-layout swing. Read `decomp_only` instead.
fn bench_measure_sweep_attribution(c: &mut Criterion) {
    let mut g = c.benchmark_group("tableau-attrib/measure-sweep");
    g.warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(4))
        .sample_size(30);

    for n in [32usize, 128] {
        let mut old_prepared: OldWide = Driver::new_seeded(n, 1e-10, 1);
        let mut new_prepared: NewWide = Driver::new_seeded(n, 1e-10, 1);
        scaling_prepare(&mut old_prepared);
        scaling_prepare(&mut new_prepared);

        // FRAME half only: the same sweep on the bare frame. No amplitudes, no
        // scratch, no allocation — `Measure for Tableau` is decomposition + (in
        // case a) one RNG bool + the row projection.
        g.bench_with_input(BenchmarkId::new("frame_sweep/old", n), &n, |b, _| {
            b.iter_batched_ref(
                || old_prepared.tableau.clone(),
                |t| {
                    for i in 0..n {
                        black_box(ppvm_traits::Measure::measure(t, i));
                    }
                },
                criterion::BatchSize::SmallInput,
            )
        });
        g.bench_with_input(BenchmarkId::new("frame_sweep/new", n), &n, |b, _| {
            b.iter_batched_ref(
                || new_prepared.tableau.clone(),
                |t| {
                    for i in 0..n {
                        black_box(ppvm_traits_2::Measure::measure(t, i));
                    }
                },
                criterion::BatchSize::SmallInput,
            )
        });

        // FRAME half, decomposition only: `z_expectation` stops after
        // `compute_decomposition` and never projects, so
        // `(frame_sweep − decomp_sweep)` is the projection's share.
        g.bench_with_input(BenchmarkId::new("decomp_sweep/old", n), &n, |b, _| {
            b.iter(|| {
                let mut acc = 0.0;
                for i in 0..n {
                    acc += old_prepared.z_expectation(i);
                }
                black_box(acc)
            })
        });
        g.bench_with_input(BenchmarkId::new("decomp_sweep/new", n), &n, |b, _| {
            b.iter(|| {
                let mut acc = 0.0;
                for i in 0..n {
                    acc += new_prepared.z_expectation(i);
                }
                black_box(acc)
            })
        });

        // The remaining variable: per-call scratch handling. `measure_many` runs
        // the SAME `n` measurements in the same order with ONE scratch held for
        // the whole sweep, whereas the per-qubit `measure` loop above goes
        // through `with_scratch` once per call. Everything else — decomposition,
        // case-a/case-b dispatch, arithmetic, RNG order, record — is identical
        // by construction (contract 7). So `(measure loop − measure_many)` on
        // one engine is exactly that engine's per-call scratch overhead, and
        // comparing the two engines' gaps isolates it as a single variable.
        // The dominant per-measurement term, benched with NOTHING else attached.
        // `decomp_sweep` above goes through `z_expectation`, whose case-a
        // `FxHashMap` over the whole support swamps the decomposition and hides
        // its contribution; this calls the `O(n)` kernel directly on both
        // engines, same frame, same qubit order.
        g.bench_with_input(BenchmarkId::new("decomp_only/old", n), &n, |b, _| {
            b.iter(|| {
                for i in 0..n {
                    black_box(old_prepared.compute_decomposition(i, ppvm_traits::Pauli::Z));
                }
            })
        });
        g.bench_with_input(BenchmarkId::new("decomp_only/new", n), &n, |b, _| {
            b.iter(|| {
                for i in 0..n {
                    black_box(new_prepared.compute_decomposition(i, ppvm_traits_2::Pauli::Z));
                }
            })
        });

        let all: Vec<usize> = (0..n).collect();
        g.bench_with_input(BenchmarkId::new("many_sweep/old", n), &n, |b, _| {
            b.iter_batched_ref(
                || old_prepared.fork(None),
                |t| black_box(t.measure_many(&all)),
                criterion::BatchSize::SmallInput,
            )
        });
        g.bench_with_input(BenchmarkId::new("many_sweep/new", n), &n, |b, _| {
            b.iter_batched_ref(
                || new_prepared.fork(None),
                |t| black_box(t.measure_many(&all)),
                criterion::BatchSize::SmallInput,
            )
        });
    }
    g.finish();
}

criterion_group! {
    name = benches;
    config = integration_config();
    targets = bench_msd, bench_rot2, bench_fused_tgate, bench_scaling,
              bench_measure_all, bench_noisy_shots, bench_branch_coalesce,
              bench_micro, bench_measure_sweep_attribution
}
criterion_main!(benches);
