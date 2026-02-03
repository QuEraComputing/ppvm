use super::Tableau;
use crate::config::Config;
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
