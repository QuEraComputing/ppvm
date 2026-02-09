use super::GeneralizedTableau;
use super::Tableau;
use super::sparsevec::SparseVector;
use crate::config::Config;
use num::complex::Complex;
use std::fmt::Display;

impl<const N: usize, T: Config> Display for Tableau<N, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Tableau ({} qubits):", N)?;
        writeln!(f, "  Destabilizers: [")?;
        for stab in self.destabilizers.iter() {
            writeln!(f, "    {}", stab)?;
        }
        writeln!(f, "  ]")?;
        writeln!(f, "  Stabilizers: [")?;
        for stab in self.stabilizers.iter() {
            writeln!(f, "    {}", stab)?;
        }
        writeln!(f, "  ]")?;
        Ok(())
    }
}

impl<const N: usize, T: Config, C: SparseVector<Complex<T::Coeff>>> Display
    for GeneralizedTableau<N, T, C>
where
    Complex<T::Coeff>: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Generalized Tableau ({} qubits):", N)?;
        writeln!(f, "  Tableau:")?;
        writeln!(f, "{}", self.tableau)?;
        writeln!(f, "  Coefficients:")?;
        for (coeff, idx) in self.coefficients.clone().into_iter() {
            writeln!(f, "    Index {}: {}", idx, coeff)?;
        }
        writeln!(f, "  Is Lost: [")?;
        for (i, &lost) in self.is_lost.iter().enumerate() {
            writeln!(f, "    Qubit {}: {}", i, lost)?;
        }
        writeln!(f, "  ]")?;
        Ok(())
    }
}
