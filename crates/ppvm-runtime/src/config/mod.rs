use crate::traits::{ACMap, Coefficient, PauliStorage};

pub trait Config {
    type Storage: PauliStorage;
    type Coeff: Coefficient;
    type Map: ACMap<Self::Storage, Self::Coeff>;
}

pub mod gxhash;

#[cfg(feature = "dashmap")]
pub mod dashmap;

#[cfg(feature = "indexmap")]
pub mod indexmap;

#[cfg(feature = "fxhash")]
pub mod fxhash;
