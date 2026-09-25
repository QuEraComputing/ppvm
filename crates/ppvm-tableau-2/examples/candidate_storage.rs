// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Standalone storage experiment, not a replacement simulator backend.
//!
//! Run with `<fixture file or directory> <samples> <minimum milliseconds>`.
//! Each fixture contains `n`, then `2n` lines of `phase x-bits z-bits`, with
//! qubit zero first. Phases mean `i^phase` times a tensor of I/X/Y/Z.
//! Both packing axes use the minimum whole-word padding and two packed phase
//! planes. There are no inverse caches, quadrant splits, or alignment promises
//! beyond the element type. Setup and conversion are outside measured regions.

use std::hint::black_box;
use std::ops::{BitAnd, BitOr, BitXor, Not, Shl};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

trait Word:
    Copy
    + Default
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + BitXor<Output = Self>
    + Not<Output = Self>
    + Shl<usize, Output = Self>
    + PartialEq
{
    const BITS: usize;
    const ONE: Self;
    fn count(self) -> u32;
}

macro_rules! words {
    ($($t:ty),*) => {$(impl Word for $t {
        const BITS: usize = <$t>::BITS as usize;
        const ONE: Self = 1;
        fn count(self) -> u32 { self.count_ones() }
    })*};
}
words!(u8, u16, u32, u64, u128);

#[derive(Clone, Debug, PartialEq, Eq)]
struct Row {
    phase: u8,
    sites: Vec<u8>,
}

#[derive(Clone)]
struct Packed<W, const COLUMNS: bool> {
    n: usize,
    stride: usize,
    x: Vec<W>,
    z: Vec<W>,
    lo: Vec<W>,
    hi: Vec<W>,
}

fn bit<W: Word>(words: &[W], i: usize) -> bool {
    words[i / W::BITS] & (W::ONE << (i % W::BITS)) != W::default()
}

fn set_bit<W: Word>(words: &mut [W], i: usize, value: bool) {
    let mask = W::ONE << (i % W::BITS);
    let word = &mut words[i / W::BITS];
    *word = (*word & !mask) | if value { mask } else { W::default() };
}

impl<W: Word, const COLUMNS: bool> Packed<W, COLUMNS> {
    fn new(rows: &[Row]) -> Self {
        let n = rows[0].sites.len();
        let stride = (if COLUMNS { 2 * n } else { n }).div_ceil(W::BITS);
        let len = (if COLUMNS { n } else { 2 * n }) * stride;
        let mut data = Self {
            n,
            stride,
            x: vec![W::default(); len],
            z: vec![W::default(); len],
            lo: vec![W::default(); (2 * n).div_ceil(W::BITS)],
            hi: vec![W::default(); (2 * n).div_ceil(W::BITS)],
        };
        for (r, row) in rows.iter().enumerate() {
            data.set_phase(r, row.phase);
            for (q, &site) in row.sites.iter().enumerate() {
                data.set_site(r, q, site);
            }
        }
        data
    }

    fn index(&self, r: usize, q: usize) -> usize {
        if COLUMNS {
            q * self.stride * W::BITS + r
        } else {
            r * self.stride * W::BITS + q
        }
    }

    fn site(&self, r: usize, q: usize) -> u8 {
        let i = self.index(r, q);
        u8::from(bit(&self.x, i)) | (u8::from(bit(&self.z, i)) << 1)
    }

    fn set_site(&mut self, r: usize, q: usize, site: u8) {
        let i = self.index(r, q);
        set_bit(&mut self.x, i, site & 1 != 0);
        set_bit(&mut self.z, i, site & 2 != 0);
    }

    fn phase(&self, r: usize) -> u8 {
        u8::from(bit(&self.lo, r)) | (u8::from(bit(&self.hi, r)) << 1)
    }

    fn set_phase(&mut self, r: usize, phase: u8) {
        set_bit(&mut self.lo, r, phase & 1 != 0);
        set_bit(&mut self.hi, r, phase & 2 != 0);
    }

    fn rows(&self) -> Vec<Row> {
        (0..2 * self.n)
            .map(|r| Row {
                phase: self.phase(r),
                sites: (0..self.n).map(|q| self.site(r, q)).collect(),
            })
            .collect()
    }

    fn anticommutes(&self, a: usize, b: usize) -> bool {
        if COLUMNS {
            let mut parity = 0;
            for q in 0..self.n {
                let (a, b) = (self.site(a, q), self.site(b, q));
                parity ^= ((a & 1) & (b >> 1)) ^ ((a >> 1) & (b & 1));
            }
            parity != 0
        } else {
            let mut parity = W::default();
            for i in 0..self.stride {
                let (a, b) = (a * self.stride + i, b * self.stride + i);
                parity = parity ^ (self.x[a] & self.z[b]) ^ (self.z[a] & self.x[b]);
            }
            parity.count() & 1 != 0
        }
    }

    // dst <- dst * src, including the phase of anticommuting products.
    fn multiply(&mut self, dst: usize, src: usize) {
        let mut phase = u64::from(self.phase(dst) + self.phase(src));
        if COLUMNS {
            for q in 0..self.n {
                let (a, b) = (self.site(dst, q), self.site(src, q));
                phase += u64::from(PRODUCT_PHASE[a as usize][b as usize]);
                self.set_site(dst, q, a ^ b);
            }
        } else {
            for i in 0..self.stride {
                let (dst, src) = (dst * self.stride + i, src * self.stride + i);
                let (a, b, c, d) = (self.x[dst], self.z[dst], self.x[src], self.z[src]);
                let sign = (a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d);
                let imag = (a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d);
                phase += u64::from(2 * sign.count() + imag.count());
                self.x[dst] = a ^ c;
                self.z[dst] = b ^ d;
            }
        }
        self.set_phase(dst, (phase & 3) as u8);
    }

    fn gate(&mut self, operation: &str, q: usize) {
        let target = (q + 1) % self.n;
        if COLUMNS {
            for i in 0..self.stride {
                let a = q * self.stride + i;
                let (x, z) = (self.x[a], self.z[a]);
                let sign = match operation {
                    "h" => {
                        self.x[a] = z;
                        self.z[a] = x;
                        x & z
                    }
                    "s" => {
                        self.z[a] = z ^ x;
                        x & z
                    }
                    "cnot" => {
                        let b = target * self.stride + i;
                        let (xt, zt) = (self.x[b], self.z[b]);
                        self.x[b] = xt ^ x;
                        self.z[a] = z ^ zt;
                        x & zt & !(xt ^ z)
                    }
                    _ => unreachable!(),
                };
                self.hi[i] = self.hi[i] ^ sign;
            }
        } else {
            for r in 0..2 * self.n {
                let a = self.site(r, q);
                let (x, z) = (a & 1, a >> 1);
                let (site, sign) = match operation {
                    "h" => ((x << 1) | z, x & z),
                    "s" => (x | ((z ^ x) << 1), x & z),
                    "cnot" => {
                        let b = self.site(r, target);
                        let (xt, zt) = (b & 1, b >> 1);
                        self.set_site(r, target, (xt ^ x) | (zt << 1));
                        (x | ((z ^ zt) << 1), x & zt & !(xt ^ z))
                    }
                    _ => unreachable!(),
                };
                self.set_site(r, q, site);
                self.set_phase(r, self.phase(r) ^ (2 * sign));
            }
        }
    }

    fn sweep(&mut self, operation: &str) -> u64 {
        match operation {
            "comm" => (0..2 * self.n)
                .map(|dst| u64::from(black_box(self.anticommutes(dst, (dst + 1) % (2 * self.n)))))
                .sum(),
            "mul" => {
                for dst in 0..2 * self.n {
                    self.multiply(dst, (dst + 1) % (2 * self.n));
                }
                0
            }
            "circuit" => {
                for q in 0..self.n {
                    for gate in ["h", "s", "cnot"] {
                        self.gate(gate, q);
                    }
                }
                0
            }
            gate => {
                for q in 0..self.n {
                    self.gate(gate, q);
                }
                0
            }
        }
    }
}

// Scalar Pauli algebra oracle: site codes are I=0, X=1, Z=2, Y=3.
const PRODUCT_PHASE: [[u8; 4]; 4] = [[0, 0, 0, 0], [0, 0, 3, 1], [0, 1, 0, 3], [0, 3, 1, 0]];
const CNOT: [(u8, u8, u8); 16] = [
    (0, 0, 0),
    (0, 1, 0),
    (2, 2, 0),
    (2, 3, 0),
    (1, 1, 0),
    (1, 0, 0),
    (3, 3, 2),
    (3, 2, 0),
    (2, 0, 0),
    (2, 1, 0),
    (0, 2, 0),
    (0, 3, 0),
    (3, 1, 0),
    (3, 0, 0),
    (1, 3, 0),
    (1, 2, 2),
];
const OPERATIONS: [&str; 6] = ["comm", "mul", "h", "s", "cnot", "circuit"];

fn scalar_gate(rows: &mut [Row], operation: &str, q: usize) {
    let target = (q + 1) % rows[0].sites.len();
    for row in rows {
        let site = row.sites[q] as usize;
        let (next, phase) = match operation {
            "h" => [(0, 0), (2, 0), (1, 0), (3, 2)][site],
            "s" => [(0, 0), (3, 0), (2, 0), (1, 2)][site],
            "cnot" => {
                let (control, target_site, phase) = CNOT[4 * site + row.sites[target] as usize];
                row.sites[target] = target_site;
                (control, phase)
            }
            _ => unreachable!(),
        };
        row.sites[q] = next;
        row.phase = (row.phase + phase) & 3;
    }
}

fn scalar_sweep(rows: &mut [Row], operation: &str) {
    if operation == "mul" {
        for dst in 0..rows.len() {
            let src = (dst + 1) % rows.len();
            let mut phase = u64::from(rows[dst].phase + rows[src].phase);
            for q in 0..rows[dst].sites.len() {
                let (a, b) = (rows[dst].sites[q], rows[src].sites[q]);
                phase += u64::from(PRODUCT_PHASE[a as usize][b as usize]);
                rows[dst].sites[q] = a ^ b;
            }
            rows[dst].phase = (phase & 3) as u8;
        }
    } else {
        for q in 0..rows[0].sites.len() {
            if operation == "circuit" {
                for gate in ["h", "s", "cnot"] {
                    scalar_gate(rows, gate, q);
                }
            } else {
                scalar_gate(rows, operation, q);
            }
        }
    }
}

fn hash(bytes: impl IntoIterator<Item = u8>) -> u64 {
    bytes.into_iter().fold(14695981039346656037, |h, byte| {
        (h ^ u64::from(byte)).wrapping_mul(1099511628211)
    })
}

fn verify<W: Word, const COLUMNS: bool>(fixture: &[Row], operation: &str) -> u64 {
    let mut packed = Packed::<W, COLUMNS>::new(fixture);
    assert_eq!(packed.rows(), fixture);
    if operation == "comm" {
        let answers: Vec<_> = (0..fixture.len())
            .map(|dst| {
                let src = (dst + 1) % fixture.len();
                let expected = fixture[dst]
                    .sites
                    .iter()
                    .zip(&fixture[src].sites)
                    .filter(|&(a, b)| a != b && *a != 0 && *b != 0)
                    .count()
                    % 2
                    != 0;
                let actual = packed.anticommutes(dst, src);
                assert_eq!(actual, expected);
                u8::from(actual)
            })
            .collect();
        hash(answers)
    } else {
        let mut expected = fixture.to_vec();
        scalar_sweep(&mut expected, operation);
        packed.sweep(operation);
        let actual = packed.rows();
        assert_eq!(
            actual,
            expected,
            "{operation}, {} bits, columns={COLUMNS}",
            W::BITS
        );
        hash(
            actual
                .iter()
                .flat_map(|row| std::iter::once(row.phase).chain(row.sites.iter().copied())),
        )
    }
}

fn benchmark<W: Word, const COLUMNS: bool>(fixture: &[Row], samples: usize, minimum: Duration) {
    let original = Packed::<W, COLUMNS>::new(fixture);
    let layout = if COLUMNS {
        "generator_bits"
    } else {
        "qubit_bits"
    };
    let n = original.n;
    for operation in OPERATIONS {
        let checksum = verify::<W, COLUMNS>(fixture, operation);
        let operations = match operation {
            "comm" | "mul" => 2 * n,
            "circuit" => 3 * n,
            _ => n,
        };
        let timed = |iterations| {
            let mut data = original.clone();
            let start = Instant::now();
            for _ in 0..iterations {
                black_box(black_box(&mut data).sweep(black_box(operation)));
            }
            let elapsed = start.elapsed();
            black_box(data);
            elapsed
        };
        // Also warms each monomorphized kernel before recording samples.
        let mut iterations = 1usize;
        while timed(iterations) < minimum {
            iterations = iterations.checked_mul(2).expect("iteration count overflow");
        }
        for sample in 0..samples {
            let elapsed = timed(iterations).as_nanos();
            println!(
                "candidate,{layout},{},{n},{operation},{iterations},{sample},{elapsed},{operations},{checksum:016x}",
                W::BITS
            );
        }
    }
}

fn read_fixture(path: &Path) -> Vec<Row> {
    let text = std::fs::read_to_string(path).expect("read fixture");
    let mut fields = text.split_whitespace();
    let n: usize = fields.next().expect("n").parse().expect("integer n");
    assert!(n >= 2, "fixtures must have at least two qubits");
    let rows = (0..2 * n)
        .map(|_| {
            let phase: u8 = fields
                .next()
                .expect("phase")
                .parse()
                .expect("integer phase");
            assert!(phase < 4);
            let x = fields.next().expect("X bits");
            let z = fields.next().expect("Z bits");
            assert_eq!(x.len(), n);
            assert_eq!(z.len(), n);
            let sites = x
                .bytes()
                .zip(z.bytes())
                .map(|(x, z)| {
                    assert!(matches!(x, b'0' | b'1') && matches!(z, b'0' | b'1'));
                    (x - b'0') | ((z - b'0') << 1)
                })
                .collect();
            Row { phase, sites }
        })
        .collect();
    assert!(fields.next().is_none(), "unexpected fixture content");
    rows
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(
        args.len(),
        4,
        "usage: candidate_storage FIXTURE_OR_DIRECTORY SAMPLES MIN_MS"
    );
    let input = PathBuf::from(&args[1]);
    let samples: usize = args[2].parse().expect("integer samples");
    let milliseconds: f64 = args[3].parse().expect("minimum milliseconds");
    assert!(samples > 0 && milliseconds.is_finite() && milliseconds > 0.0);
    let minimum = Duration::from_secs_f64(milliseconds / 1000.0);
    let mut files = if input.is_dir() {
        std::fs::read_dir(input)
            .expect("read fixture directory")
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.extension().is_some_and(|extension| extension == "txt"))
            .collect()
    } else {
        vec![input]
    };
    files.sort();
    assert!(!files.is_empty(), "no fixtures");
    println!(
        "implementation,layout,word_bits,n,operation,iterations,sample,elapsed_ns,operations_per_iteration,checksum"
    );
    for path in files {
        let rows = read_fixture(&path);
        macro_rules! run { ($($t:ty),*) => {$(
            benchmark::<$t, false>(&rows, samples, minimum);
            benchmark::<$t, true>(&rows, samples, minimum);
        )*}; }
        run!(u8, u16, u32, u64, u128);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_sweeps_match_pauli_algebra_at_word_boundaries() {
        for n in [
            2, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 255, 256, 257,
        ] {
            let fixture: Vec<_> = (0..2 * n)
                .map(|r| Row {
                    phase: (r % 4) as u8,
                    sites: (0..n)
                        .map(|q| ((r * 17 + q * 3 + r * q / 3) % 4) as u8)
                        .collect(),
                })
                .collect();
            for operation in OPERATIONS {
                let expected = verify::<u64, false>(&fixture, operation);
                macro_rules! check { ($($t:ty),*) => {$(
                    assert_eq!(verify::<$t, false>(&fixture, operation), expected);
                    assert_eq!(verify::<$t, true>(&fixture, operation), expected);
                )*}; }
                check!(u8, u16, u32, u64, u128);
            }
        }
    }
}
