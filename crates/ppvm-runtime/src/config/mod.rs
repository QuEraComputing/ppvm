use std::hash::BuildHasher;

use crate::traits::{ACMap, Coefficient, PauliStorage, PauliWordTrait, Strategy};

pub trait Config: Clone {
    type Storage: PauliStorage;
    type Coeff: Coefficient;
    type Strategy: Strategy;
    type BuildHasher: BuildHasher + Clone + Default;
    type PauliWordType: PauliWordTrait;
    type Map: ACMap<Self::Storage, Self::Coeff, Self::BuildHasher, Self::PauliWordType>;
}

pub mod fxhash;

#[cfg(feature = "dashmap")]
pub mod dashmap;

#[cfg(feature = "indexmap")]
pub mod indexmap;

#[cfg(feature = "gxhash")]
pub mod gxhash;
