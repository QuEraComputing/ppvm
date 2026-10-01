// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! GPU run: TFIM / Heisenberg Trotter benchmarks on `CudaPauliSum`.
//!
//! ```bash
//! CUPAULIPROP_LIB_DIR=... cargo run --release -p ppvm-cupauliprop --features cuda \
//!     --example trotter_gpu -- --model heisenberg --out results_gpu_heisenberg.json
//! ```
//!
//! `LD_LIBRARY_PATH` must contain the cuPauliProp and CUDA runtime library
//! directories.

mod common;

use clap::Parser;
use common::{BenchArgs, bench_size, write_results};
use ppvm_cupauliprop::cuda::{CallTimes, device_memory, library_version};
use ppvm_cupauliprop::{CudaOptions, CudaPauliSum};
use serde_json::json;

#[derive(Parser)]
#[command(about = "TFIM / Heisenberg Trotter benchmark on the cuPauliProp GPU backend")]
struct Cli {
    #[command(flatten)]
    bench: BenchArgs,
    /// Term capacity per expansion [default: derived from free GPU memory].
    #[arg(long)]
    capacity: Option<usize>,
    /// Workspace size in MiB [default: derived from free GPU memory].
    #[arg(long)]
    workspace_mib: Option<usize>,
    /// Print host time per apply phase (view, prepare, workspace, compute).
    #[arg(long)]
    timings: bool,
}

/// Per-apply mean host time of each phase, for the last round.
fn print_call_times(t: &CallTimes) {
    let per = |d: std::time::Duration| d.as_secs_f64() * 1e6 / t.applies.max(1) as f64;
    eprintln!(
        "    per apply ({} applies): view {:.1} µs  prepare {:.1} µs  workspace {:.1} µs  compute {:.1} µs",
        t.applies,
        per(t.view),
        per(t.prepare),
        per(t.workspace),
        per(t.compute)
    );
}

fn main() {
    let cli = Cli::parse();
    let args = &cli.bench;
    let options = CudaOptions {
        capacity: cli.capacity,
        workspace_bytes: cli.workspace_mib.map(|m| m << 20),
    };
    let (free, total) = device_memory();
    eprintln!(
        "cuPauliProp {} | GPU memory free {} / {} MiB",
        library_version(),
        free >> 20,
        total >> 20
    );

    let results = args
        .sizes()
        .into_iter()
        .map(|size| {
            let (row, state) = bench_size(
                args,
                size,
                |n, terms| CudaPauliSum::new(n, terms, args.cutoff, options),
                CudaPauliSum::synchronize,
            );
            if cli.timings {
                print_call_times(&state.call_times());
            }
            row
        })
        .collect();
    let info = json!({
        "name": "cupauliprop",
        "cupauliprop_version": library_version(),
        "gpu_memory_total_bytes": total,
        "capacity": cli.capacity,
        "workspace_mib": cli.workspace_mib,
    });
    write_results(args, "ppvm-cupauliprop", info, results);
}
