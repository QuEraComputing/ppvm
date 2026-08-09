// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Human-readable rendering of the frame and the generalized tableau.
//!
//! Ported verbatim from `ppvm-tableau/src/display.rs` — byte-for-byte the same
//! output, including the blank line [`GeneralizedTableau`] emits after the
//! nested frame (old wrote `writeln!(f, "{}", self.tableau)` and the frame's own
//! rendering already ends in a newline). `Display` is a user-facing surface, so
//! under the prime directive the layout is reproduced rather than improved.
//!
//! Row rendering matches `ppvm-pauli-word`'s `Display for PhasedPauliWord`: a
//! `+` / `+i` / `-` / `-i` phase prefix followed by one `I`/`X`/`Y`/`Z` letter
//! per site.

use std::fmt::{self, Display};

use ppvm_traits_2::Pauli;

use crate::data::{Bitstring, GeneralizedTableau, Row, RowStorage, Tableau};

/// `+` / `+i` / `-` / `-i` for the `ℤ/4` phase convention `0: +1, 1: +i,
/// 2: −1, 3: −i`.
const PHASE_PREFIX: [&str; 4] = ["+", "+i", "-", "-i"];

/// One row as `<phase-prefix><Pauli letters>`, e.g. `-iXYZI`.
///
/// A free function rather than a `Display` impl on [`Row`] because a row does
/// not carry its own width — the frame supplies it.
fn fmt_row<A: RowStorage>(
    f: &mut fmt::Formatter<'_>,
    row: &Row<A>,
    n_qubits: usize,
) -> fmt::Result {
    debug_assert!(row.phase < 4, "Invalid phase value: {}", row.phase);
    f.write_str(PHASE_PREFIX[(row.phase % 4) as usize])?;
    for i in 0..n_qubits {
        f.write_str(match row.get(i) {
            Pauli::I => "I",
            Pauli::X => "X",
            Pauli::Y => "Y",
            Pauli::Z => "Z",
        })?;
    }
    Ok(())
}

impl<A: RowStorage, H> Display for Tableau<A, H> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Tableau ({} qubits):", self.n_qubits)?;
        writeln!(f, "  Destabilizers: [")?;
        for row in self.destabilizers().iter() {
            f.write_str("    ")?;
            fmt_row(f, row, self.n_qubits)?;
            f.write_str("\n")?;
        }
        writeln!(f, "  ]")?;
        writeln!(f, "  Stabilizers: [")?;
        for row in self.stabilizers().iter() {
            f.write_str("    ")?;
            fmt_row(f, row, self.n_qubits)?;
            f.write_str("\n")?;
        }
        writeln!(f, "  ]")?;
        Ok(())
    }
}

impl<A: RowStorage, I: Bitstring + Display, H> Display for GeneralizedTableau<A, I, H> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Generalized Tableau ({} qubits):", self.tableau.n_qubits)?;
        writeln!(f, "  Tableau:")?;
        writeln!(f, "{}", self.tableau)?;
        writeln!(f, "  Coefficients:")?;
        for &(coeff, idx) in self.coefficients.iter() {
            writeln!(f, "    Index {idx}: {coeff}")?;
        }
        writeln!(f, "  Is Lost: [")?;
        for (i, &lost) in self.is_lost.iter().enumerate() {
            writeln!(f, "    Qubit {i}: {lost}")?;
        }
        writeln!(f, "  ]")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    #[test]
    fn tableau_display_matches_the_old_layout() {
        let mut tab: Tableau = Tableau::new(2);
        tab.h(0);
        assert_eq!(
            tab.to_string(),
            "Tableau (2 qubits):\n  \
             Destabilizers: [\n    \
             +ZI\n    \
             +IX\n  \
             ]\n  \
             Stabilizers: [\n    \
             +XI\n    \
             +IZ\n  \
             ]\n"
        );
    }

    #[test]
    fn generalized_tableau_display_renders_frame_coefficients_and_loss() {
        let tab: GeneralizedTableau = GeneralizedTableau::new(1, 1e-12);
        let rendered = tab.to_string();
        assert!(rendered.starts_with("Generalized Tableau (1 qubits):\n  Tableau:\n"));
        // The nested frame ends in a newline and `writeln!` adds another, so a
        // blank line separates the frame from the coefficients — as in old.
        assert!(rendered.contains("  ]\n\n  Coefficients:\n    Index 0: 1+0i\n"));
        assert!(rendered.ends_with("  Is Lost: [\n    Qubit 0: false\n  ]\n"));
    }
}
