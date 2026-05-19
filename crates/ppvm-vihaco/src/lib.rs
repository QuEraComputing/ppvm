pub mod component;
pub mod instruction;
// pub mod machine;
pub mod message;

pub mod prelude {
    pub use crate::component::{Circuit, CircuitEffect};
    pub use crate::instruction::CircuitInstruction;
    pub use crate::message::CircuitMessage;
}
