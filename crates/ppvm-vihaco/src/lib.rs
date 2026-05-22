use crate::composite::PPVM;

pub mod component;
pub mod composite;
pub mod instruction;
pub mod measurements;
pub mod message;
mod syntax;

pub fn run_file(path: &str) -> eyre::Result<PPVM> {
    let mut machine = PPVM::default();
    machine.run_file(path)?;
    Ok(machine)
}

pub fn run_program(program: &str) -> eyre::Result<PPVM> {
    let mut machine = PPVM::default();
    machine.run_program(program)?;
    Ok(machine)
}

pub mod prelude {
    pub use crate::component::Circuit;
    pub use crate::composite::PPVM;
}
