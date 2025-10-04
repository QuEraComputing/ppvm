pub mod char;
pub mod config;
mod map;
pub mod pattern;
pub mod phase;
pub mod strategy;
pub mod sum;
pub mod traits;
pub mod word;

pub mod prelude {
    pub use crate::char::Pauli;
    pub use crate::config;
    pub use crate::config::Config;
    pub use crate::pattern::PauliPattern;
    pub use crate::phase::PhasedPauliWord;
    pub use crate::sum::PauliSum;
    pub use crate::traits::*;
    pub use crate::word::PauliWord;
}
