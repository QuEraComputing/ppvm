use std::hash::BuildHasher;

use crate::traits::{ACMap, Coefficient, PauliStorage, Strategy};

pub trait Config: Clone {
    type Storage: PauliStorage;
    type Coeff: Coefficient;
    type Strategy: Strategy;
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
