// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! TFIM Trotter runtime vs qubit count, swept across the same storage-tier
//! ladder the Python bindings dispatch on, for both `fxhash` and `gxhash`.
//!
//! This reproduces the "bump then converge" scaling the PR investigated: as
//! qubits grow within a storage tier the map fills up, and with `fxhash` the
//! cached-hash low bits cluster hashbrown's buckets at high fill, so runtime
//! balloons toward the top of a tier and then drops when the next (wider,
//! lower-relative-fill) tier kicks in. `gxhash` avalanches well even on short
//! keys, so it should stay smooth.
//!
//! Circuit / parameters mirror `benches/trotter.rs` and the Python
//! `test_trotter.py` (so the numbers line up with PauliPropagation.jl):
//! TFIM, h = 1, dt = 0.1, truncation 1e-6, depolarizing noise 1e-4 applied as
//! `pauli_error([1e-4/4; 3])`. `J` and the Trotter `STEPS` are env-tunable;
//! the defaults (J = 1/8, 10 steps) match the committed benchmark but keep the
//! state small (~1e3 terms), so all three hashers are flat. The cliff is a
//! *high-fill* effect — drive the state large with `J=1.0 STEPS=20` to
//! reproduce the bump (this is what `benchmarks/plot_tfim_sweep.py` plots).
//!
//! Three ppvm series are emitted per qubit count: `fxhash_nofold` (a no-fold
//! `HashFinalize` newtype reproducing pre-PR fxhash), `fxhash` (folded, this
//! PR), and `gxhash`.
//!
//! Storage tier for `n` qubits = `2^k` bytes, the smallest power-of-two byte
//! width with `2^k > ceil(n / 8)` — exactly Python's `_init_ppvm_interface`
//! dispatch.
//!
//! Usage (gxhash needs AES):
//! ```bash
//! RUSTFLAGS="-C target-feature=+aes" J=1.0 STEPS=20 \
//!   QUBITS="8,16,24,32,40,48,56,64,80,96,112,122" ITERS=2 \
//!   cargo run --release -p ppvm-pauli-sum --example trotter_qubit_sweep > sweep.csv
//! ```

use std::hash::BuildHasher;
use std::time::Instant;

use ppvm_pauli_sum::prelude::*;
use ppvm_pauli_sum::strategy::CoefficientThreshold;
use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::traits::{HashFinalize, PauliWordTrait};

/// `FxBuildHasher` with the identity `HashFinalize` — i.e. the pre-PR
/// behavior with **no** high-bit fold. Used to reproduce the bucket-clustering
/// "bump" that the fold (and gxhash) remove.
#[derive(Clone, Default)]
struct RawFx(fxhash::FxBuildHasher);

impl BuildHasher for RawFx {
    type Hasher = <fxhash::FxBuildHasher as BuildHasher>::Hasher;
    #[inline]
    fn build_hasher(&self) -> Self::Hasher {
        self.0.build_hasher()
    }
}

// Identity finalize == no fold == the cliff is back.
impl HashFinalize for RawFx {}

const H: f64 = 1.0;
const DT: f64 = 0.1 / H;
const NOISE: [f64; 3] = [1e-4 / 4.0; 3];

/// Tunable circuit parameters (defaults mirror `trotter.rs` / PP.jl).
#[derive(Clone, Copy)]
struct Params {
    steps: usize,
    theta_x: f64,
    theta_zz: f64,
    min_abs_coeff: f64,
}

impl Params {
    fn from_env() -> Self {
        let env_f64 = |k: &str, d: f64| {
            std::env::var(k)
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(d)
        };
        let j = env_f64("J", 1.0 / 8.0) * H;
        let steps = std::env::var("STEPS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(10);
        Params {
            steps,
            theta_x: DT * H,
            theta_zz: DT * j,
            min_abs_coeff: env_f64("MIN_ABS_COEFF", 1e-6),
        }
    }
}

/// Smallest power-of-two byte width with `2^k > ceil(n/8)` — Python's tier.
fn storage_bytes(n: usize) -> usize {
    let need = n.div_ceil(8);
    let mut k = 0;
    while (1usize << k) <= need {
        k += 1;
    }
    1usize << k
}

/// Run the TFIM Trotter circuit `iters` times, return (min seconds, |state|).
fn run_trotter<C, T>(n: usize, iters: usize, p: Params) -> (f64, usize)
where
    C: Config<Coeff = f64, Strategy = CoefficientThreshold, PauliWordType = T>,
    T: PauliWordTrait,
    for<'a> &'a str: Into<T>,
{
    let strat = CoefficientThreshold(p.min_abs_coeff);
    let mut seed: PauliSum<C> = PauliSum::builder()
        .n_qubits(n)
        .strategy(strat)
        .capacity(n * n)
        .build();
    // initial observable: sum_i Z_i
    for i in 0..n {
        let term: String = (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect();
        seed += (term.as_str(), 1.0);
    }

    let mut best = f64::INFINITY;
    let mut final_len = 0usize;
    for _ in 0..iters {
        let mut state = seed.clone();
        let t0 = Instant::now();
        for _ in 0..p.steps {
            for i in 0..n {
                state.rx(i, p.theta_x);
                state.truncate();
                state.pauli_error(i, NOISE);
                state.truncate();
            }
            for i in 0..n - 1 {
                state.rzz(i, i + 1, p.theta_zz);
                state.truncate();
                state.pauli_error(i, NOISE);
                state.truncate();
                state.pauli_error(i + 1, NOISE);
                state.truncate();
            }
        }
        best = best.min(t0.elapsed().as_secs_f64());
        final_len = state.len();
    }
    (best, final_len)
}

/// Dispatch `n` to the right storage tier for a given indexmap config family.
macro_rules! sweep {
    ($family:ident, $n:expr, $iters:expr, $p:expr) => {
        match storage_bytes($n) {
            2 => {
                run_trotter::<config::indexmap::$family<2, CoefficientThreshold>, _>($n, $iters, $p)
            }
            4 => {
                run_trotter::<config::indexmap::$family<4, CoefficientThreshold>, _>($n, $iters, $p)
            }
            8 => {
                run_trotter::<config::indexmap::$family<8, CoefficientThreshold>, _>($n, $iters, $p)
            }
            16 => run_trotter::<config::indexmap::$family<16, CoefficientThreshold>, _>(
                $n, $iters, $p,
            ),
            32 => run_trotter::<config::indexmap::$family<32, CoefficientThreshold>, _>(
                $n, $iters, $p,
            ),
            b => panic!("no config for {b}-byte storage (n={})", $n),
        }
    };
}

/// Same dispatch but with the no-fold `RawFx` word hasher (pre-PR fxhash).
macro_rules! sweep_rawfx {
    ($n:expr, $iters:expr, $p:expr) => {{
        type Cfg<const N: usize> =
            config::indexmap::ByteFxHash<N, f64, CoefficientThreshold, PauliWord<[u8; N], RawFx>>;
        match storage_bytes($n) {
            2 => run_trotter::<Cfg<2>, _>($n, $iters, $p),
            4 => run_trotter::<Cfg<4>, _>($n, $iters, $p),
            8 => run_trotter::<Cfg<8>, _>($n, $iters, $p),
            16 => run_trotter::<Cfg<16>, _>($n, $iters, $p),
            32 => run_trotter::<Cfg<32>, _>($n, $iters, $p),
            b => panic!("no config for {b}-byte storage (n={})", $n),
        }
    }};
}

fn main() {
    let qubits: Vec<usize> = std::env::var("QUBITS")
        .ok()
        .map(|s| s.split(',').filter_map(|t| t.trim().parse().ok()).collect())
        .unwrap_or_else(|| {
            vec![
                8, 16, 24, 32, 40, 44, 48, 52, 56, 60, 64, 72, 80, 88, 96, 104, 112, 120, 122,
            ]
        });
    let iters: usize = std::env::var("ITERS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);

    let p = Params::from_env();
    eprintln!(
        "params: steps={} theta_x={:.4} theta_zz={:.4} min_abs_coeff={:.0e}",
        p.steps, p.theta_x, p.theta_zz, p.min_abs_coeff
    );

    // CSV to stdout; progress to stderr so a redirect keeps clean data.
    println!("qubits,hasher,bytes,time_s,terms");
    for &n in &qubits {
        let bytes = storage_bytes(n);
        let it = if n > 72 { iters.min(2) } else { iters };

        let (raw_t, raw_len) = sweep_rawfx!(n, it, p);
        println!("{n},fxhash_nofold,{bytes},{raw_t:.6},{raw_len}");
        eprintln!("n={n:3} bytes={bytes:2} fxhash_nofold {raw_t:8.4}s ({raw_len} terms)");

        let (fx_t, fx_len) = sweep!(ByteFxHashF64, n, it, p);
        println!("{n},fxhash,{bytes},{fx_t:.6},{fx_len}");
        eprintln!("n={n:3} bytes={bytes:2} fxhash        {fx_t:8.4}s ({fx_len} terms)");

        let (gx_t, gx_len) = sweep!(ByteGxHashF64, n, it, p);
        println!("{n},gxhash,{bytes},{gx_t:.6},{gx_len}");
        eprintln!("n={n:3} bytes={bytes:2} gxhash        {gx_t:8.4}s ({gx_len} terms)");
    }
}
