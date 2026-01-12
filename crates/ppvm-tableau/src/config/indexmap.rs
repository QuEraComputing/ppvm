use std::marker::PhantomData;

use num::complex::Complex64;
use ppvm_runtime::traits::Coefficient;

use crate::{config::Config, Tableau};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Index<C: Coefficient + 'static>(PhantomData<C>);

impl<C: Coefficient + 'static> Config for Index<C> {
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type Map = indexmap::IndexMap<Tableau, C, Self::BuildHasher>;
}

pub type IndexComplex = Index<Complex64>;
