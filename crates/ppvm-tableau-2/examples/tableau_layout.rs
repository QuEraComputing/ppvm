// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Native worker for `benchmarks/tableau-layout/run.py`.
//!
//! The path import compiles the production storage and kernels without exposing
//! private internals in the library API. `ppvm-native-kernel` excludes the public
//! frame's inverse-sign cache; `ppvm-public` includes that cache and hash updates.
//! Each CSV checksum describes one sweep on the original fixture, outside timing.

#[allow(dead_code, unused_imports)]
#[path = "../src/storage/mod.rs"]
mod storage;

use std::hint::black_box;
use std::path::Path;
use std::time::{Duration, Instant};

use ppvm_tableau_2::Tableau;
use ppvm_traits_2::{Clifford, Pauli, StabilizerFrame};
use storage::{HALVES, Half, Orientation, Plane, TableauData, blocks};

const FNV_OFFSET: u64 = 14_695_981_039_346_656_037;
const FNV_PRIME: u64 = 1_099_511_628_211;

#[derive(Clone, Copy)]
enum Operation {
    Comm,
    Mul,
    H,
    S,
    Cnot,
    Circuit,
    TransposeRoundtrip,
}

impl Operation {
    fn name(self) -> &'static str {
        match self {
            Self::Comm => "comm",
            Self::Mul => "mul",
            Self::H => "h",
            Self::S => "s",
            Self::Cnot => "cnot",
            Self::Circuit => "circuit",
            Self::TransposeRoundtrip => "transpose_roundtrip",
        }
    }

    fn count(self, n: usize) -> usize {
        match self {
            Self::Comm | Self::Mul => 2 * n,
            Self::Circuit => 3 * n,
            Self::TransposeRoundtrip => 1,
            _ => n,
        }
    }
}

#[derive(Clone, Copy)]
enum Layout {
    Row,
    Column,
    TransposeBatch,
}

impl Layout {
    fn name(self) -> &'static str {
        match self {
            Self::Row => "row",
            Self::Column => "column",
            Self::TransposeBatch => "column-transpose-batch",
        }
    }
}

#[derive(Clone)]
struct Native {
    data: TableauData,
    source_x: Vec<u64>,
    source_z: Vec<u64>,
    target_x: Vec<u64>,
    target_z: Vec<u64>,
    commutation_sum: u64,
}

impl Native {
    fn from_fixture(path: &Path) -> Self {
        let text = std::fs::read_to_string(path).expect("read fixture");
        let mut lines = text.lines();
        let n: usize = lines.next().expect("fixture size").parse().expect("size");
        assert!(
            n >= 2,
            "the cyclic two-qubit sweep requires at least two qubits"
        );
        let mut data = TableauData::identity(n);
        data.transpose_quadrants();
        for generator in 0..2 * n {
            let mut fields = lines.next().expect("fixture row").split_whitespace();
            let phase: u8 = fields
                .next()
                .expect("phase")
                .parse()
                .expect("phase integer");
            assert!(phase < 4);
            let (half, row) = Half::split(generator, n);
            data.set_phase_of(half, row, phase);
            for plane in [Plane::X, Plane::Z] {
                let bits = fields.next().expect("bitstring").as_bytes();
                assert_eq!(bits.len(), n);
                let words = data.major_mut(half, plane, row);
                words.fill(0);
                for (qubit, &value) in bits.iter().enumerate() {
                    assert!(value == b'0' || value == b'1');
                    TableauData::set_bit(words, qubit, value == b'1');
                }
            }
            assert!(fields.next().is_none());
        }
        assert!(lines.all(|line| line.trim().is_empty()));
        data.transpose_quadrants();
        let stride = data.stride();
        Self {
            data,
            source_x: vec![0; stride],
            source_z: vec![0; stride],
            target_x: vec![0; stride],
            target_z: vec![0; stride],
            commutation_sum: 0,
        }
    }

    fn commutation(&mut self, dst: usize) -> u8 {
        let n = self.data.n_qubits();
        let (src_half, src_i) = Half::split((dst + 1) % (2 * n), n);
        let (dst_half, dst_i) = Half::split(dst, n);
        let count = if self.data.orientation() == Orientation::RowMajor {
            blocks::and_count(
                self.data.major(dst_half, Plane::X, dst_i),
                self.data.major(src_half, Plane::Z, src_i),
            ) + blocks::and_count(
                self.data.major(dst_half, Plane::Z, dst_i),
                self.data.major(src_half, Plane::X, src_i),
            )
        } else {
            self.data
                .gather_row(src_half, src_i, &mut self.source_x, &mut self.source_z);
            self.data
                .gather_row(dst_half, dst_i, &mut self.target_x, &mut self.target_z);
            blocks::and_count(&self.target_x, &self.source_z)
                + blocks::and_count(&self.target_z, &self.source_x)
        };
        (count & 1) as u8
    }

    fn rows(&mut self, operation: Operation) {
        let n = self.data.n_qubits();
        self.commutation_sum = 0;
        for dst in 0..2 * n {
            if matches!(operation, Operation::Comm) {
                self.commutation_sum += u64::from(self.commutation(dst));
                continue;
            }
            let src = (dst + 1) % (2 * n);
            let (src_half, src_i) = Half::split(src, n);
            let (dst_half, dst_i) = Half::split(dst, n);
            if self.data.orientation() == Orientation::RowMajor {
                self.source_x
                    .copy_from_slice(self.data.major(src_half, Plane::X, src_i));
                self.source_z
                    .copy_from_slice(self.data.major(src_half, Plane::Z, src_i));
                let phase = self.data.phase_of(src_half, src_i);
                self.data
                    .multiply_row_by(dst_half, dst_i, &self.source_x, &self.source_z, phase);
            } else {
                self.data
                    .gather_row(src_half, src_i, &mut self.source_x, &mut self.source_z);
                self.data
                    .gather_row(dst_half, dst_i, &mut self.target_x, &mut self.target_z);
                let delta = blocks::row_multiply(
                    &mut self.target_x,
                    &mut self.target_z,
                    &self.source_x,
                    &self.source_z,
                );
                let phase = (self.data.phase_of(dst_half, dst_i)
                    + self.data.phase_of(src_half, src_i)
                    + delta)
                    % 4;
                self.data.set_phase_of(dst_half, dst_i, phase);
                for q in 0..n {
                    TableauData::set_bit(
                        self.data.major_mut(dst_half, Plane::X, q),
                        dst_i,
                        TableauData::bit(&self.target_x, q),
                    );
                    TableauData::set_bit(
                        self.data.major_mut(dst_half, Plane::Z, q),
                        dst_i,
                        TableauData::bit(&self.target_z, q),
                    );
                }
            }
        }
    }

    fn gate(&mut self, operation: Operation, q: usize) {
        for half in HALVES {
            if matches!(operation, Operation::Cnot) {
                let target = (q + 1) % self.data.n_qubits();
                let (xc, zc, xt, zt, phase) = self.data.gate2_mut(half, q, target);
                blocks::cnot(xc, zc, xt, zt, phase);
            } else {
                let (x, z, phase) = self.data.gate1_mut(half, q);
                match operation {
                    Operation::H => blocks::h(x, z, phase),
                    Operation::S => blocks::s(x, z, phase),
                    _ => unreachable!(),
                }
            }
        }
    }

    fn sweep(&mut self, layout: Layout, operation: Operation) {
        if matches!(operation, Operation::TransposeRoundtrip) {
            self.data.transpose_quadrants();
            self.data.transpose_quadrants();
        } else if matches!(operation, Operation::Comm | Operation::Mul) {
            if matches!(layout, Layout::TransposeBatch) {
                self.data.transpose_quadrants();
            }
            self.rows(operation);
            if matches!(layout, Layout::TransposeBatch) {
                self.data.transpose_quadrants();
            }
        } else {
            for q in 0..self.data.n_qubits() {
                if matches!(operation, Operation::Circuit) {
                    for gate in [Operation::H, Operation::S, Operation::Cnot] {
                        self.gate(gate, q);
                    }
                } else {
                    self.gate(operation, q);
                }
            }
        }
    }

    fn checksum(&self, layout: Layout, operation: Operation) -> u64 {
        if matches!(operation, Operation::Comm) {
            // Keep the FNV work outside timing, using the same row kernel to
            // expose every result to the harness's independent scalar oracle.
            let n = self.data.n_qubits();
            let mut data = self.clone();
            if matches!(layout, Layout::TransposeBatch) {
                data.data.transpose_quadrants();
            }
            let mut hash = FNV_OFFSET;
            let mut sum = 0;
            for dst in 0..2 * n {
                let anticommutes = data.commutation(dst);
                sum += u64::from(anticommutes);
                hash = fnv(hash, anticommutes);
            }
            assert_eq!(
                self.commutation_sum, sum,
                "commutation checksum differs from consumed result"
            );
            return hash;
        }
        matrix_checksum(
            self.data.n_qubits(),
            |row| self.data.phase(row),
            |row, q| u8::from(self.data.x_bit(row, q)) | (u8::from(self.data.z_bit(row, q)) << 1),
        )
    }
}

fn fnv(hash: u64, byte: u8) -> u64 {
    (hash ^ u64::from(byte)).wrapping_mul(FNV_PRIME)
}

fn matrix_checksum(
    n: usize,
    phase: impl Fn(usize) -> u8,
    cell: impl Fn(usize, usize) -> u8,
) -> u64 {
    let mut hash = FNV_OFFSET;
    for row in 0..2 * n {
        hash = fnv(hash, phase(row));
        for q in 0..n {
            hash = fnv(hash, cell(row, q));
        }
    }
    hash
}

fn public_checksum(data: &Tableau) -> u64 {
    matrix_checksum(
        data.n_qubits(),
        |row| data.row_phase(row),
        |row, q| match data.row_site(row, q) {
            Pauli::I => 0,
            Pauli::X => 1,
            Pauli::Z => 2,
            Pauli::Y => 3,
        },
    )
}

fn public_sweep(data: &mut Tableau, operation: Operation) {
    let n = data.n_qubits();
    if matches!(operation, Operation::Mul) {
        for dst in 0..2 * n {
            data.row_multiply((dst + 1) % (2 * n), dst);
        }
        return;
    }
    for q in 0..n {
        match operation {
            Operation::H => data.h(q),
            Operation::S => data.s(q),
            Operation::Cnot => data.cnot(q, (q + 1) % n),
            Operation::Circuit => {
                data.h(q);
                data.s(q);
                data.cnot(q, (q + 1) % n);
            }
            _ => unreachable!(),
        }
    }
}

/// Initialize through the public API when the orchestrator supplied preparation gates.
fn public_fixture(path: &Path, native: &Native) -> Option<Tableau> {
    let gates = std::fs::read_to_string(path.with_extension("gates")).ok()?;
    let mut data = Tableau::new(native.data.n_qubits());
    for line in gates.lines().filter(|line| !line.trim().is_empty()) {
        let fields: Vec<_> = line.split_whitespace().collect();
        let q: usize = fields[1].parse().expect("gate target");
        match fields[0] {
            "h" => data.h(q),
            "s" => data.s(q),
            "cnot" => data.cnot(q, fields[2].parse().expect("CNOT target")),
            _ => panic!("unknown preparation gate"),
        }
    }
    assert_eq!(
        public_checksum(&data),
        native.checksum(Layout::Column, Operation::H),
        "preparation gates differ from fixture"
    );
    Some(data)
}

struct Timing {
    n: usize,
    samples: usize,
    minimum: Duration,
}

impl Timing {
    fn measure<T: Clone>(
        &self,
        labels: (&str, &str),
        operation: Operation,
        original: &T,
        sweep: impl Fn(&mut T),
        checksum: impl Fn(&T) -> u64,
    ) {
        let (implementation, layout) = labels;
        let mut check = original.clone();
        sweep(&mut check);
        let digest = checksum(&check);
        black_box(check);

        let time = |iterations: u64| {
            let mut data = original.clone();
            let start = Instant::now();
            for _ in 0..iterations {
                sweep(black_box(&mut data));
            }
            let elapsed = start.elapsed();
            black_box(data);
            elapsed
        };
        let mut iterations = 1u64;
        while time(iterations) < self.minimum {
            iterations = iterations
                .checked_mul(2)
                .expect("timer calibration overflow");
        }
        for sample in 1..=self.samples {
            let elapsed = time(iterations);
            println!(
                "{implementation},{layout},64,{},{},{iterations},{sample},{},{},{digest:016x}",
                self.n,
                operation.name(),
                elapsed.as_nanos(),
                operation.count(self.n)
            );
        }
    }
}

fn benchmark_fixture(path: &Path, samples: usize, minimum: Duration) {
    let native = Native::from_fixture(path);
    let public = public_fixture(path, &native);
    let timing = Timing {
        n: native.data.n_qubits(),
        samples,
        minimum,
    };
    assert!(timing.samples > 0 && !timing.minimum.is_zero());
    for layout in [Layout::Row, Layout::Column, Layout::TransposeBatch] {
        let mut original = native.clone();
        if matches!(layout, Layout::Row) {
            original.data.transpose_quadrants();
        }
        let operations: &[Operation] = if matches!(layout, Layout::Column) {
            &[
                Operation::Comm,
                Operation::Mul,
                Operation::H,
                Operation::S,
                Operation::Cnot,
                Operation::Circuit,
                Operation::TransposeRoundtrip,
            ]
        } else {
            &[Operation::Comm, Operation::Mul]
        };
        for &operation in operations {
            timing.measure(
                ("ppvm-native-kernel", layout.name()),
                operation,
                &original,
                |data| data.sweep(layout, operation),
                |data| data.checksum(layout, operation),
            );
        }
    }
    if let Some(original) = public {
        for operation in [
            Operation::Mul,
            Operation::H,
            Operation::S,
            Operation::Cnot,
            Operation::Circuit,
        ] {
            timing.measure(
                ("ppvm-public", "column"),
                operation,
                &original,
                |data| public_sweep(data, operation),
                public_checksum,
            );
        }
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(
        args.len(),
        4,
        "usage: tableau_layout FIXTURE_OR_DIRECTORY SAMPLES MIN_MS"
    );
    let path = Path::new(&args[1]);
    let samples = args[2].parse().expect("sample count");
    let minimum =
        Duration::from_secs_f64(args[3].parse::<f64>().expect("minimum milliseconds") / 1000.0);
    let mut paths = if path.is_dir() {
        std::fs::read_dir(path)
            .expect("fixture directory")
            .map(|entry| entry.expect("fixture entry").path())
            .filter(|path| path.extension().is_some_and(|extension| extension == "txt"))
            .collect::<Vec<_>>()
    } else {
        vec![path.to_path_buf()]
    };
    paths.sort();
    assert!(!paths.is_empty(), "no fixtures found");
    println!(
        "implementation,layout,word_bits,n,operation,iterations,sample,elapsed_ns,operations_per_iteration,checksum"
    );
    for path in paths {
        benchmark_fixture(&path, samples, minimum);
    }
}
