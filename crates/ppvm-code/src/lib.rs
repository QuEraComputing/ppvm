// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Qubit stabilizer codes defined by signed Hermitian Pauli generators.
mod readout;
pub use readout::Parity;

use ppvm_pauli_kernels::row_product;
use std::{fmt, sync::OnceLock};

/// A signed Hermitian Pauli row. Site zero is the leftmost letter.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Row {
    n: usize,
    x: Vec<u64>,
    z: Vec<u64>,
    negative: bool,
}

impl Row {
    fn identity(n: usize) -> Self {
        let words = n.div_ceil(64);
        Self {
            n,
            x: vec![0; words],
            z: vec![0; words],
            negative: false,
        }
    }

    /// Build a row from equal-length binary X/Z planes; (1,1) is Hermitian Y.
    pub fn from_bits(xs: &[bool], zs: &[bool], negative: bool) -> Result<Self, String> {
        if xs.len() != zs.len() {
            return Err("X and Z rows must have equal width".into());
        }
        let mut row = Self::identity(xs.len());
        row.negative = negative;
        for (i, (&x, &z)) in xs.iter().zip(zs).enumerate() {
            row.x[i / 64] |= u64::from(x) << (i % 64);
            row.z[i / 64] |= u64::from(z) << (i % 64);
        }
        Ok(row)
    }

    fn bit(&self, column: usize) -> bool {
        let (plane, i) = if column < self.n {
            (&self.x, column)
        } else {
            (&self.z, column - self.n)
        };
        plane[i / 64] >> (i % 64) & 1 != 0
    }

    fn set(&mut self, column: usize) {
        let (plane, i) = if column < self.n {
            (&mut self.x, column)
        } else {
            (&mut self.z, column - self.n)
        };
        plane[i / 64] |= 1 << (i % 64);
    }

    // Call only on equal-width rows.
    fn anticommutes(&self, other: &Self) -> bool {
        let mut parity = 0;
        for i in 0..self.x.len() {
            parity ^= (self.x[i] & other.z[i]) ^ (self.z[i] & other.x[i]);
        }
        parity.count_ones() % 2 == 1
    }

    fn multiply_check(&mut self, other: &Self) {
        let phase = row_product(&mut self.x, &mut self.z, &other.x, &other.z);
        debug_assert!(phase == 0 || phase == 2);
        self.negative ^= other.negative ^ (phase == 2);
    }

    // Only binary representatives matter when choosing a new logical basis.
    fn add_support(&mut self, other: &Self) {
        for i in 0..self.x.len() {
            self.x[i] ^= other.x[i];
            self.z[i] ^= other.z[i];
        }
    }
}

impl std::str::FromStr for Row {
    type Err = String;
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let letters = text.strip_prefix(['+', '-']).unwrap_or(text);
        let mut row = Row::identity(letters.len());
        row.negative = text.starts_with('-');
        for (i, letter) in letters.bytes().enumerate() {
            match letter {
                b'I' | b'_' => (),
                b'X' => row.set(i),
                b'Z' => row.set(row.n + i),
                b'Y' => {
                    row.set(i);
                    row.set(row.n + i);
                }
                _ => return Err(format!("invalid signed Pauli {text:?} at site {i}")),
            }
        }
        Ok(row)
    }
}

impl fmt::Display for Row {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(if self.negative { "-" } else { "+" })?;
        for i in 0..self.n {
            f.write_str(match (self.bit(i), self.bit(self.n + i)) {
                (false, false) => "I",
                (true, false) => "X",
                (false, true) => "Z",
                (true, true) => "Y",
            })?;
        }
        Ok(())
    }
}

/// The joint +1 eigenspace of the supplied signed generators.
#[derive(Clone, Debug)]
pub struct Code {
    n: usize,
    checks: Vec<Row>,
    independent: Vec<usize>,
    // A reduced basis keyed by leading column. Also used to find the centralizer.
    pivots: Vec<(usize, Row)>,
    logicals: OnceLock<Vec<(Row, Row)>>,
}

impl Code {
    /// Preserve all input checks in order, including consistent redundant checks.
    pub fn new(n: usize, checks: Vec<Row>) -> Result<Self, String> {
        if n > usize::MAX / 2 {
            return Err("tableau width is too large".into());
        }
        let mut pivots: Vec<(usize, Row)> = Vec::new();
        let mut independent = Vec::new();
        for (i, check) in checks.iter().enumerate() {
            if check.n != n {
                return Err(format!("row {i} has width {}, expected {n}", check.n));
            }
            for &j in &independent {
                if check.anticommutes(&checks[j]) {
                    return Err(format!("rows {j} and {i} anticommute"));
                }
            }
            let mut reduced = check.clone();
            for (column, basis) in &pivots {
                if reduced.bit(*column) {
                    reduced.multiply_check(basis);
                }
            }
            if let Some(column) = (0..2 * n).find(|&c| reduced.bit(c)) {
                for (_, basis) in &mut pivots {
                    if basis.bit(column) {
                        basis.multiply_check(&reduced);
                    }
                }
                let pos = pivots.partition_point(|(c, _)| *c < column);
                pivots.insert(pos, (column, reduced));
                independent.push(i);
            } else if reduced.negative {
                return Err(format!(
                    "row {i} is inconsistent: the generators produce -I"
                ));
            }
        }
        Ok(Self {
            n,
            checks,
            independent,
            pivots,
            logicals: OnceLock::new(),
        })
    }

    pub fn n(&self) -> usize {
        self.n
    }
    pub fn k(&self) -> usize {
        self.n - self.pivots.len()
    }
    /// Test signed membership; e.g. a code stabilized by -ZZ excludes +ZZ.
    pub fn contains_stabilizer(&self, op: &Row) -> Result<bool, String> {
        if op.n != self.n {
            return Err("Pauli width must equal code width".into());
        }
        if self.pivots.iter().any(|(_, s)| op.anticommutes(s)) {
            return Ok(false);
        }
        let mut row = op.clone();
        for (pivot, basis) in &self.pivots {
            if row.bit(*pivot) {
                row.multiply_check(basis);
            }
        }
        Ok(!row.negative && row.x.iter().chain(&row.z).all(|&w| w == 0))
    }

    pub fn stabilizers(&self) -> &[Row] {
        &self.checks
    }
    pub fn independent_stabilizers(&self) -> impl Iterator<Item = &Row> {
        self.independent.iter().map(|&i| &self.checks[i])
    }

    /// Deterministic logical X/Z pairs in the centralizer modulo the stabilizers.
    pub fn logical_pairs(&self) -> &[(Row, Row)] {
        self.logicals.get_or_init(|| self.derive_logicals())
    }

    /// Select a logical convention; commutation modulo the checks guarantees
    /// that all 2k representatives are independent modulo the stabilizer group.
    pub fn with_logicals(mut self, xs: Vec<Row>, zs: Vec<Row>) -> Result<Self, String> {
        if xs.len() != self.k() || zs.len() != self.k() {
            return Err(format!("expected {} logical X/Z pairs", self.k()));
        }
        for (i, row) in xs.iter().chain(&zs).enumerate() {
            if row.n != self.n {
                return Err(format!("logical row {i} has incorrect width"));
            }
            if self.pivots.iter().any(|(_, s)| s.anticommutes(row)) {
                return Err(format!("logical row {i} anticommutes with a stabilizer"));
            }
        }
        for i in 0..xs.len() {
            for j in 0..xs.len() {
                if xs[i].anticommutes(&xs[j])
                    || zs[i].anticommutes(&zs[j])
                    || xs[i].anticommutes(&zs[j]) != (i == j)
                {
                    return Err(format!(
                        "invalid commutation between logical operators at indices {i} and {j}"
                    ));
                }
            }
        }
        self.logicals = OnceLock::from(xs.into_iter().zip(zs).collect::<Vec<_>>());
        Ok(self)
    }

    fn derive_logicals(&self) -> Vec<(Row, Row)> {
        if self.k() == 0 {
            return Vec::new();
        }
        let n = self.n;
        // Null vectors of [X|Z] are centralizing Paulis after swapping X and Z.
        let swap = |c: usize| if c < n { c + n } else { c - n };
        let mut space = Vec::new();
        for free in 0..2 * n {
            if self.pivots.binary_search_by_key(&free, |(p, _)| *p).is_ok() {
                continue;
            }
            let mut row = Row::identity(n);
            row.set(swap(free));
            for (pivot, equation) in &self.pivots {
                if equation.bit(free) {
                    row.set(swap(*pivot));
                }
            }
            space.push(row);
        }
        // Symplectic Gram-Schmidt. Unpaired rows are in the stabilizer span.
        let mut pairs = Vec::with_capacity(self.k());
        while let Some(x) = space.pop() {
            let Some(i) = space.iter().position(|r| x.anticommutes(r)) else {
                continue;
            };
            let z = space.remove(i);
            for r in &mut space {
                let toggle_z = r.anticommutes(&x);
                let toggle_x = r.anticommutes(&z);
                if toggle_z {
                    r.add_support(&z);
                }
                if toggle_x {
                    r.add_support(&x);
                }
            }
            pairs.push((x, z));
        }
        debug_assert_eq!(pairs.len(), self.k());
        pairs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn code(checks: &[&str]) -> Result<Code, String> {
        Code::new(
            checks[0].len() - 1,
            checks.iter().map(|s| s.parse()).collect::<Result<_, _>>()?,
        )
    }

    #[test]
    fn signed_dependencies() {
        let good = code(&["+XX", "+ZZ", "-YY", "+II"]).unwrap();
        assert_eq!((good.n(), good.k(), good.stabilizers().len()), (2, 0, 4));
        assert!(code(&["+XX", "+ZZ", "+YY"]).unwrap_err().contains("-I"));
        assert!(code(&["+II", "-II"]).is_err());
        assert!(code(&["+XZ", "+ZZ"]).is_err());
    }

    #[test]
    fn perfect_five_logicals() {
        let c = code(&["+XZZXI", "+IXZZX", "+XIXZZ", "+ZXIXZ"]).unwrap();
        let pairs = c.logical_pairs();
        assert_eq!(pairs.len(), 1);
        let (x, z) = &pairs[0];
        assert!(x.anticommutes(z));
        assert!(
            c.checks
                .iter()
                .all(|s| !s.anticommutes(x) && !s.anticommutes(z))
        );
    }

    #[test]
    fn logical_pairs_across_storage_words() {
        let n = 65;
        let mut g = Row::identity(n);
        g.set(0);
        g.set(n - 1);
        let c = Code::new(n, vec![g]).unwrap();
        let (xs, zs) = c.logical_pairs().iter().cloned().unzip();
        assert_eq!(c.k(), 64);
        c.with_logicals(xs, zs).unwrap();
    }
}
