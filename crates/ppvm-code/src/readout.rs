// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::{Code, Row};

/// A parity of single-qubit measurement bits, XORed with a fixed offset.
#[derive(Debug, PartialEq, Eq)]
pub struct Parity {
    pub qubits: Vec<usize>,
    pub offset: bool,
}

impl Row {
    fn incompatible(&self, i: usize, axis: u8) -> bool {
        let x = self.bit(i);
        let z = self.bit(self.n + i);
        match axis {
            b'X' => z,
            b'Z' => x,
            b'Y' => x ^ z,
            _ => unreachable!("validated measurement axis"),
        }
    }

    fn parity(&self) -> Parity {
        Parity {
            qubits: (0..self.n)
                .filter(|&i| self.bit(i) || self.bit(self.n + i))
                .collect(),
            offset: self.negative,
        }
    }
}

impl Code {
    /// Derive deterministic terminal checks and logical parities. Axes are in
    /// physical-qubit order. Logical axes are in logical order; I skips a slot.
    /// Each resulting detector is zero on an ideal shot from the codespace.
    pub fn terminal_readout(
        &self,
        axes: &[u8],
        logical_axes: &[u8],
    ) -> Result<(Vec<Parity>, Vec<Parity>), String> {
        if axes.len() != self.n || !axes.iter().all(|b| b"XYZ".contains(b)) {
            return Err("supply one X, Y or Z measurement axis for every physical qubit".into());
        }
        if logical_axes.len() != self.k() || !logical_axes.iter().all(|b| b"IXYZ".contains(b)) {
            return Err("supply one X, Y, Z or I (skip) for every logical qubit".into());
        }
        // Eliminate components that anticommute with a local readout axis,
        // carrying the signed physical row alongside the projected pivot.
        let mut pivots: Vec<(usize, Row)> = Vec::new();
        let mut detectors = Vec::new();
        for check in self.independent_stabilizers() {
            let mut row = check.clone();
            reduce(&mut row, axes, &pivots);
            if let Some(i) = axes
                .iter()
                .enumerate()
                .position(|(i, &a)| row.incompatible(i, a))
            {
                let pos = pivots.partition_point(|(c, _)| *c < i);
                pivots.insert(pos, (i, row));
            } else {
                detectors.push(row.parity());
            }
        }

        let mut observables = Vec::new();
        for (i, (&axis, (x, z))) in logical_axes.iter().zip(self.logical_pairs()).enumerate() {
            let mut row = match axis {
                b'X' => x.clone(),
                b'Z' => z.clone(),
                b'Y' => {
                    // Logical Y = i X Z in the signed logical convention.
                    let mut y = x.clone();
                    let phase = ppvm_pauli_kernels::row_product(&mut y.x, &mut y.z, &z.x, &z.z);
                    y.negative ^= z.negative ^ ((phase + 1) % 4 == 2);
                    y
                }
                _ => continue,
            };
            reduce(&mut row, axes, &pivots);
            if axes
                .iter()
                .enumerate()
                .any(|(i, &a)| row.incompatible(i, a))
            {
                return Err(format!(
                    "logical {} at index {i} is inaccessible in this readout",
                    char::from(axis)
                ));
            }
            observables.push(row.parity());
        }
        Ok((detectors, observables))
    }
}

fn reduce(row: &mut Row, axes: &[u8], pivots: &[(usize, Row)]) {
    for (i, p) in pivots {
        if row.incompatible(*i, axes[*i]) {
            row.multiply_check(p);
        }
    }
}
