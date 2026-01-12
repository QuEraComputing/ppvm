use std::hash::BuildHasher;

use ppvm_runtime::traits::Coefficient;

use crate::map::TableauMap;

pub trait Config: Clone {
    type Coeff: Coefficient + 'static;
    type BuildHasher: BuildHasher + Clone + Default;
    type Map: TableauMap<Self::Coeff, Self::BuildHasher>;
}

pub mod fxhash;

#[cfg(feature = "indexmap")]
pub mod indexmap;

#[cfg(feature = "gxhash")]
pub mod gxhash;
