// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Pure-Clifford scaling comparison: a GHZ chain (`H 0`, `CX i i+1`, `M *`) on
//! both engines, split into the gate sweep and the measurement sweep so the
//! per-measurement transpose can be attributed separately from the gate
//! kernels.

use ppvm_conformance_2::tableau::{Driver, NewWide, OldWide};
use std::time::{Duration, Instant};

fn ghz_gates<D: Driver>(tab: &mut D, n: usize) {
    tab.h(0);
    for i in 0..n - 1 {
        tab.cnot(i, i + 1);
    }
}

/// Time `body` over enough repetitions to accumulate ~`budget`, returning the
/// mean per-repetition duration.
fn time<D: Driver, F: FnMut(&mut D)>(base: &D, budget: Duration, mut body: F) -> Duration {
    let mut reps = 0u32;
    let start = Instant::now();
    while start.elapsed() < budget {
        for _ in 0..16 {
            let mut tab = base.fork(Some(1));
            body(&mut tab);
            reps += 1;
        }
    }
    start.elapsed() / reps
}

fn row<D: Driver>(n: usize, budget: Duration) -> (Duration, Duration, Duration, Duration) {
    let base: D = Driver::new_seeded(n, 1e-10, 1);
    let gates = time(&base, budget, |tab: &mut D| {
        ghz_gates(tab, n);
        std::hint::black_box(());
    });

    let mut prepared: D = Driver::new_seeded(n, 1e-10, 1);
    ghz_gates(&mut prepared, n);
    let first = time(&prepared, budget, |tab: &mut D| {
        std::hint::black_box(tab.measure(0));
    });
    let measure = time(&prepared, budget, |tab: &mut D| {
        for i in 0..n {
            std::hint::black_box(tab.measure(i));
        }
    });

    let whole = time(&base, budget, |tab: &mut D| {
        ghz_gates(tab, n);
        for i in 0..n {
            std::hint::black_box(tab.measure(i));
        }
    });
    (gates, first, measure, whole)
}

fn main() {
    let budget = Duration::from_millis(
        std::env::args()
            .nth(1)
            .and_then(|s| s.parse().ok())
            .unwrap_or(300),
    );

    println!(
        "{:>6}  {:>24}  {:>24}  {:>24}  {:>24}",
        "n", "gates old/new", "1st measure old/new", "measure-all old/new", "whole old/new"
    );
    for n in [2usize, 4, 8, 16, 32, 64, 128] {
        let (go, fo, mo, wo) = row::<OldWide>(n, budget);
        let (gn, fnn, mn, wn) = row::<NewWide>(n, budget);
        let ratio = |old: Duration, new: Duration| old.as_secs_f64() / new.as_secs_f64();
        println!(
            "{n:>6}  {:>9.2?} {:>9.2?} {:>4.2}x  {:>9.2?} {:>9.2?} {:>4.2}x  {:>9.2?} {:>9.2?} {:>4.2}x  {:>9.2?} {:>9.2?} {:>4.2}x",
            go,
            gn,
            ratio(go, gn),
            fo,
            fnn,
            ratio(fo, fnn),
            mo,
            mn,
            ratio(mo, mn),
            wo,
            wn,
            ratio(wo, wn),
        );
    }
}
