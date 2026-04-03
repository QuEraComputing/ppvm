use super::GeneralizedTableau;
use super::Tableau;
use super::sparsevec::SparseVector;
use crate::config::Config;
use crate::tableau::tableau_index::TableauIndex;
use num::complex::Complex;
use std::fmt::Display;
use std::ops::{BitAnd, Shl};

impl<T: Config> Display for Tableau<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Tableau ({} qubits):", self.n_qubits)?;
        writeln!(f, "  Destabilizers: [")?;
        for stab in self.destabilizers().iter() {
            writeln!(f, "    {}", stab)?;
        }
        writeln!(f, "  ]")?;
        writeln!(f, "  Stabilizers: [")?;
        for stab in self.stabilizers().iter() {
            writeln!(f, "    {}", stab)?;
        }
        writeln!(f, "  ]")?;
        Ok(())
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> Display for GeneralizedTableau<T, I, C>
where
    Complex<T::Coeff>: Display,
    <T as Config>::Coeff: num::Num,
    I: TableauIndex + Display,
    <I as BitAnd<<I as Shl<usize>>::Output>>::Output: PartialEq<I>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Generalized Tableau ({} qubits):", self.tableau.n_qubits)?;
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
