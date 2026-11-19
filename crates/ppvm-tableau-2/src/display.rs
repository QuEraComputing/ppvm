// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `Display` for the simulator and its parts, built from their public
//! accessors.

use std::fmt::{Display, Formatter, Result};

use crate::{
    ActiveDormant, ActiveState, DenseState, GeneralizedTableau, InverseTableau, RowRef,
    VirtualState,
};

/// A row as a signed Pauli string over `n` qubits, qubit 0 first: `-XZI`.
fn write_row(f: &mut Formatter<'_>, row: RowRef<'_>, n: usize) -> Result {
    write!(f, "{}", if row.negative { '-' } else { '+' })?;
    for q in 0..n {
        let c = match (row.x_bit(q), row.z_bit(q)) {
            (false, false) => 'I',
            (true, false) => 'X',
            (true, true) => 'Y',
            (false, true) => 'Z',
        };
        write!(f, "{c}")?;
    }
    Ok(())
}

/// One line per qubit: the rows `C⁻¹ X_q C` and `C⁻¹ Z_q C`.
impl Display for InverseTableau {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let n = self.num_qubits();
        writeln!(
            f,
            "InverseTableau ({n} qubits, rows C⁻¹ X_q C and C⁻¹ Z_q C):"
        )?;
        for q in 0..n {
            write!(f, "  {q}: ")?;
            write_row(f, self.ix(q), n)?;
            write!(f, "  ")?;
            write_row(f, self.iz(q), n)?;
            writeln!(f)?;
        }
        Ok(())
    }
}

/// The nonzero amplitudes, axis 0 as the rightmost bit.
impl Display for DenseState {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let k = self.num_axes();
        for (alpha, c) in self.amplitudes().iter().enumerate() {
            if c.norm_sqr() > 0.0 {
                writeln!(f, "  |{alpha:0k$b}⟩: {c}")?;
            }
        }
        Ok(())
    }
}

/// The virtual qubit of each axis, then `φ` over the axes.
impl<S: ActiveState + Display> Display for ActiveDormant<S> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        writeln!(
            f,
            "Active virtual qubits (axis 0 first): {:?}",
            self.active()
        )?;
        write!(f, "{}", self.phi())
    }
}

impl<V: VirtualState + Display> Display for GeneralizedTableau<V> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let n = self.num_qubits();
        writeln!(f, "GeneralizedTableau ({n} qubits):")?;
        write!(f, "{}", self.tableau())?;
        write!(f, "{}", self.state())?;
        let lost: Vec<usize> = (0..n).filter(|&q| self.is_lost(q)).collect();
        let leaked: Vec<usize> = (0..n).filter(|&q| self.is_leaked(q)).collect();
        writeln!(f, "Lost: {lost:?}, leaked: {leaked:?}")
    }
}
