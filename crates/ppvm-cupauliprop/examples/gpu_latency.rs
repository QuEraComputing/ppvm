// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Per-call GPU round-trip latency on this host, to separate cuPauliProp's
//! own per-apply cost from scheduling delays (e.g. a GPU shared with other
//! processes).
//!
//! 1. raw CUDA: empty sync, 8-byte memset + sync, 8-byte device→host copy;
//! 2. `CudaPauliSum::rx` on a 1-term state, split into apply phases.
//!
//! ```bash
//! cargo run --release -p ppvm-cupauliprop --features cuda --example gpu_latency
//! ```

use std::ffi::c_void;
use std::ptr;
use std::time::Instant;

use ppvm_cupauliprop::{CudaOptions, CudaPauliSum};
use ppvm_cupauliprop_sys::*;

const REPS: usize = 200;

fn check(err: cudaError_t) {
    assert_eq!(err, 0, "CUDA call failed with error {err}");
}

/// Print min and median of `REPS` timed calls of `f`, in µs.
fn measure(label: &str, mut f: impl FnMut()) {
    let mut us: Vec<f64> = (0..REPS)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed().as_secs_f64() * 1e6
        })
        .collect();
    us.sort_by(f64::total_cmp);
    println!(
        "{label:<32} min {:>9.1} µs   median {:>9.1} µs",
        us[0],
        us[REPS / 2]
    );
}

fn main() {
    let mut dev: *mut c_void = ptr::null_mut();
    unsafe {
        check(cudaMalloc(&mut dev, 8));
        check(cudaDeviceSynchronize());
    }
    let mut host = 0u64;

    measure("cudaDeviceSynchronize (empty)", || unsafe {
        check(cudaDeviceSynchronize());
    });
    measure("cudaMemset 8 B + sync", || unsafe {
        check(cudaMemset(dev, 0, 8));
        check(cudaDeviceSynchronize());
    });
    measure("cudaMemcpy 8 B device->host", || unsafe {
        check(cudaMemcpy(
            (&raw mut host).cast(),
            dev,
            8,
            CUDA_MEMCPY_DEVICE_TO_HOST,
        ));
    });

    let options = CudaOptions {
        capacity: Some(1 << 16),
        workspace_bytes: Some(64 << 20),
    };
    let mut state = CudaPauliSum::new(4, &["ZIII".to_owned()], 1e-6, options);
    measure("CudaPauliSum::rx, 1 term", || state.rx(1, 0.1));
    let t = state.call_times();
    let per = |d: std::time::Duration| d.as_secs_f64() * 1e6 / t.applies as f64;
    println!(
        "  mean per apply: view {:.1} µs  prepare {:.1} µs  workspace {:.1} µs  compute {:.1} µs",
        per(t.view),
        per(t.prepare),
        per(t.workspace),
        per(t.compute)
    );

    unsafe { check(cudaFree(dev)) };
}
