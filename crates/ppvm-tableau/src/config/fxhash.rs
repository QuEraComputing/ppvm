use std::collections::HashMap;
use std::marker::PhantomData;

use num::complex::Complex64;
use ppvm_runtime::traits::Coefficient;

use crate::{config::Config, Tableau};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Fx<C: Coefficient + 'static>(PhantomData<C>);

impl<C: Coefficient + 'static> Config for Fx<C> {
    type Coeff = C;
    type BuildHasher = fxhash::FxBuildHasher;
    type Map = HashMap<Tableau, C, Self::BuildHasher>;
}

pub type FxComplex = Fx<Complex64>;
