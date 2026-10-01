// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! CPU baseline: TFIM / Heisenberg Trotter benchmarks on ppvm's `PauliSum`.
//!
//! ```bash
//! RAYON_NUM_THREADS=1 cargo run --release -p ppvm-cupauliprop --example trotter_cpu -- \
//!     --model heisenberg --out results_cpu_heisenberg.json
//! ```

mod common;

use clap::Parser;
use common::{BenchArgs, bench_size, write_results};
use ppvm_cupauliprop::CpuPauliSum;
use ppvm_cupauliprop::cpu::storage_bytes;
use serde_json::{Value, json};

#[derive(Parser)]
#[command(about = "TFIM / Heisenberg Trotter benchmark on the ppvm CPU PauliSum")]
struct Cli {
    #[command(flatten)]
    bench: BenchArgs,
}

fn run<const B: usize>(args: &BenchArgs, size: usize) -> Value {
    bench_size(
        args,
        size,
        |n, terms| CpuPauliSum::<B>::new(n, terms, args.cutoff),
        |_| {},
    )
}

fn main() {
    let args = Cli::parse().bench;
    let mut results = Vec::new();
    for size in args.sizes() {
        let row = match storage_bytes(args.n_qubits(size)) {
            2 => run::<2>(&args, size),
            4 => run::<4>(&args, size),
            8 => run::<8>(&args, size),
            16 => run::<16>(&args, size),
            32 => run::<32>(&args, size),
            64 => run::<64>(&args, size),
            128 => run::<128>(&args, size),
            256 => run::<256>(&args, size),
            b => panic!("unsupported storage width {b} bytes"),
        };
        results.push(row);
    }
    write_results(&args, "ppvm (rust)", json!({ "name": "cpu" }), results);
}
