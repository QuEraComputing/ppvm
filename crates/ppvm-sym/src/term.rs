use std::{collections::BTreeMap, hash::Hash};

use fxhash::FxHashMap;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Item {
    Sin(u32),
    Cos(u32),
}

/// <coff> sin^m cos^n
///
/// note: the order of the variables matters, we always sort them in ascending order
/// so that we can have a canonical representation of the terms.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Prod {
    pub(crate) sin: BTreeMap<u32, u32>,
    pub(crate) cos: BTreeMap<u32, u32>,
}

impl Prod {
    pub fn new() -> Self {
        Self {
            sin: BTreeMap::new(),
            cos: BTreeMap::new(),
        }
    }

    pub fn sin(id: u32) -> Self {
        let mut p = Self::new();
        p.sin.insert(id, 1);
        p
    }

    pub fn sin_pow(&self) -> usize {
        self.sin.len()
    }

    pub fn cos_pow(&self) -> usize {
        self.cos.len()
    }

    pub fn cos(id: u32) -> Self {
        let mut p = Self::new();
        p.cos.insert(id, 1);
        p
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Sum {
    pub(crate) c0: f64,
    pub(crate) terms: FxHashMap<Prod, f64>,
    pub(crate) max: usize, // max sin pow
}

impl Sum {
    pub fn new(max: usize) -> Self {
        Self {
            c0: 0.0,
            terms: FxHashMap::default(),
            max,
        }
    }
}
