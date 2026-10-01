// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! GPU equivalence tests. Need a CUDA device:
//! `cargo test -p ppvm-cupauliprop --features cuda -- --ignored`.
#![cfg(feature = "cuda")]

use ppvm_cupauliprop::circuits::{heisenberg, lattice_edges, sum_z_terms, tfim};
use ppvm_cupauliprop::{CpuPauliSum, CudaOptions, CudaPauliSum, Propagator};

// Small fixed buffers so tests don't reserve the whole device.
const OPTIONS: CudaOptions = CudaOptions {
    capacity: Some(1 << 20),
    workspace_bytes: Some(256 << 20),
};

fn single_z(n: usize, q: usize) -> Vec<String> {
    vec![(0..n).map(|i| if i == q { 'Z' } else { 'I' }).collect()]
}

#[test]
#[ignore = "needs a CUDA GPU"]
fn rx_on_z_matches_analytic() {
    // R†ZR with R = exp(-iθX/2) = cos θ Z ± sin θ Y; only Z overlaps |0⟩.
    // Qubit 70 checks the second packed word.
    for q in [0, 3, 70] {
        let mut s = CudaPauliSum::new(80, &single_z(80, q), 0.0, OPTIONS);
        s.rx(q, 0.3);
        assert!(
            (s.overlap_with_zero() - 0.3f64.cos()).abs() < 1e-12,
            "q={q}"
        );
        assert_eq!(s.len(), 2);
    }
}

#[test]
#[ignore = "needs a CUDA GPU"]
fn rotation_on_other_qubit_is_identity() {
    let mut s = CudaPauliSum::new(8, &single_z(8, 0), 0.0, OPTIONS);
    s.rx(1, 0.3);
    assert!((s.overlap_with_zero() - 1.0).abs() < 1e-14);
    assert_eq!(s.len(), 1);
}

#[test]
#[ignore = "needs a CUDA GPU"]
fn pauli_error_damps_z() {
    // Heisenberg Pauli channel on Z: (1 - 2(px + py)) Z.
    let mut s = CudaPauliSum::new(4, &single_z(4, 2), 0.0, OPTIONS);
    s.pauli_error(2, [0.01, 0.02, 0.03]);
    assert!((s.overlap_with_zero() - (1.0 - 2.0 * (0.01 + 0.02))).abs() < 1e-14);
}

fn assert_close(gpu: &CudaPauliSum, cpu: &impl Propagator) {
    let (g, c) = (gpu.overlap_with_zero(), cpu.overlap_with_zero());
    assert!(
        (g - c).abs() <= 1e-10 * c.abs().max(1.0),
        "gpu {g} vs cpu {c}"
    );
    let (gl, cl) = (gpu.len() as f64, cpu.len() as f64);
    assert!((gl - cl).abs() <= 1e-3 * cl, "terms gpu {gl} vs cpu {cl}");
}

#[test]
#[ignore = "needs a CUDA GPU"]
fn tfim_matches_cpu() {
    let (n, cutoff, noise) = (8, 1e-6, [1e-4 / 4.0; 3]);
    let terms = sum_z_terms(n);
    let mut cpu = CpuPauliSum::<2>::new(n, &terms, cutoff);
    let mut gpu = CudaPauliSum::new(n, &terms, cutoff, OPTIONS);
    tfim(&mut cpu, n, 4, 0.1, 0.1, noise);
    tfim(&mut gpu, n, 4, 0.1, 0.1, noise);
    assert_close(&gpu, &cpu);
}

#[test]
#[ignore = "needs a CUDA GPU"]
fn heisenberg_matches_cpu() {
    let (l, cutoff, noise) = (3, 1e-6, [1e-4 / 4.0; 3]);
    let terms = sum_z_terms(l * l);
    let edges = lattice_edges(l);
    let mut cpu = CpuPauliSum::<4>::new(l * l, &terms, cutoff);
    let mut gpu = CudaPauliSum::new(l * l, &terms, cutoff, OPTIONS);
    heisenberg(&mut cpu, &edges, 2, 0.1, noise);
    heisenberg(&mut gpu, &edges, 2, 0.1, noise);
    assert_close(&gpu, &cpu);
}
