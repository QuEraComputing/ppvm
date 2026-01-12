use std::collections::HashMap;
use std::marker::PhantomData;

use num::complex::Complex64;
use ppvm_runtime::traits::Coefficient;

use crate::{config::Config, Tableau};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Gx<C: Coefficient + 'static>(PhantomData<C>);

impl<C: Coefficient + 'static> Config for Gx<C> {
    type Coeff = C;
    type BuildHasher = gxhash::GxBuildHasher;
    type Map = HashMap<Tableau, C, Self::BuildHasher>;
}

pub type GxComplex = Gx<Complex64>;
