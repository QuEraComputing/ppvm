use std::{collections::BTreeMap, hash::Hash};

use fxhash::FxHashMap;

/// <coff> sin^m cos^n
///
/// note: the order of the variables matters, we always sort them in ascending order
/// so that we can have a canonical representation of the terms.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Prod {
    pub(crate) sin: BTreeMap<u32, u32>,
    pub(crate) cos: BTreeMap<u32, u32>,
    pub(crate) sin_pow: usize,
    pub(crate) cos_pow: usize,
    /// phase factor mod 4, encoded as:
    /// |  | sign | imag |
    /// |--|------|------|
    /// |+1|    0 |    0 |
    /// |+i|    0 |    1 |
    /// |-1|    1 |    0 |
    /// |-i|    1 |    1 |
    pub(crate) phase: u8,
}

impl Prod {
    pub fn new() -> Self {
        Self {
            sin: BTreeMap::new(),
            cos: BTreeMap::new(),
            sin_pow: 0,
            cos_pow: 0,
            phase: 0,
        }
    }

    pub fn add_phase(&mut self, phase: u8) {
        self.phase = (self.phase + phase) % 4;
    }

    pub fn sin(id: u32) -> Self {
        let mut p = Self::new();
        p.sin.insert(id, 1);
        p.sin_pow = 1;
        p
    }

    pub fn pow(&self) -> usize {
        self.sin_pow + self.cos_pow
    }

    pub fn sin_pow(&self) -> usize {
        self.sin_pow
    }

    pub fn cos_pow(&self) -> usize {
        self.cos_pow
    }

    pub fn cos(id: u32) -> Self {
        let mut p = Self::new();
        p.cos.insert(id, 1);
        p.cos_pow = 1;
        p
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Sum {
    pub(crate) c0: f64,
    pub(crate) terms: FxHashMap<Prod, f64>,
}

impl Sum {
    pub fn new() -> Self {
        Self {
            c0: 0.0,
            terms: FxHashMap::default(),
        }
    }

    pub fn add_const(&mut self, c: f64, min_eps: f64) {
        if c.abs() < min_eps {
            return;
        }
        self.c0 += c;
    }

    pub fn add_term(&mut self, p: Prod, coeff: f64, max: usize, min_eps: f64) {
        if p.sin_pow() > max || coeff.abs() < min_eps {
            return;
        }

        if p.pow() == 0 {
            self.c0 += coeff;
            return;
        }
        *self.terms.entry(p).or_insert(0.0) += coeff;
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Inner {
    Sum(Sum),
    One(Prod, f64),
    Var(u32),
    Const(f64),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Term {
    pub(crate) inner: Inner,
    pub(crate) max_sin: usize, // max sin power
    pub(crate) min_eps: f64,   // min coefficient to keep
}

impl Term {
    pub fn set_max_sin(&mut self, max: usize) {
        self.max_sin = max;
    }

    pub fn set_min_eps(&mut self, eps: f64) {
        self.min_eps = eps;
    }

    pub fn var(u: u32) -> Self {
        Self {
            inner: Inner::Var(u),
            max_sin: usize::MAX,
            min_eps: f64::EPSILON,
        }
    }

    pub fn sin(mut self) -> Self {
        match &mut self.inner {
            Inner::Var(u) => {
                self.inner = Inner::One(Prod::sin(*u), 1.0);
            }
            Inner::Const(c) => {
                *c = (*c).sin();
            }
            _ => {
                panic!("only variable or constant can be input of sin");
            }
        }
        return self;
    }

    pub fn cos(mut self) -> Self {
        match &mut self.inner {
            Inner::Var(u) => {
                self.inner = Inner::One(Prod::cos(*u), 1.0);
            }
            Inner::Const(c) => {
                *c = (*c).cos();
            }
            _ => {
                panic!("only variable or constant can be input of cos");
            }
        }
        return self;
    }

    pub fn from_f64(c: f64) -> Self {
        Self {
            inner: Inner::Const(c),
            max_sin: usize::MAX,
            min_eps: f64::EPSILON,
        }
    }
}
