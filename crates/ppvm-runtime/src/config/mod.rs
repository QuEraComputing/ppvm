use std::hash::BuildHasher;

use crate::traits::{ACMap, Coefficient, PauliStorage};

pub trait Config {
    type Storage: PauliStorage;
    type Coeff: Coefficient;
    type BuildHasher: BuildHasher + Clone + Default;
    type Map: ACMap<Self::Storage, Self::Coeff, Self::BuildHasher>;
}

pub mod fxhash;

#[cfg(feature = "dashmap")]
pub mod dashmap;

#[cfg(feature = "indexmap")]
pub mod indexmap;

#[cfg(feature = "gxhash")]
pub mod gxhash;
