use ppvm_runtime::traits::{Clifford, ComplexCoefficient};
use rayon::prelude::*;

use crate::config::Config;
use crate::map::{
    TableauMap, TableauMapAddAssign, TableauMapBase, TableauMapDrain, TableauMapGet,
};
use crate::Tableau;

#[derive(Clone, Debug, PartialEq)]
pub struct TableauSum<T: Config> {
    n_qubits: usize,
    pub(crate) map: T::Map,
    par_threshold: usize,
}

const PAR_THRESHOLD: usize = 256;

impl<T: Config> TableauSum<T>
where
    T::Map: TableauMap<T::Coeff, T::BuildHasher> + Send,
{
    pub fn new(n_qubits: usize) -> Self {
        let mut map = T::Map::with_capacity(1);
        map.add_assign(Tableau::new(n_qubits), T::Coeff::from(1.0));
        Self {
            n_qubits,
            map,
            par_threshold: PAR_THRESHOLD,
        }
    }

    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.map.len() == 0
    }

    pub fn with_parallel_threshold(mut self, threshold: usize) -> Self {
        self.par_threshold = threshold;
        self
    }

    pub fn set_parallel_threshold(&mut self, threshold: usize) {
        self.par_threshold = threshold;
    }

    pub fn coeff(&self, tableau: &Tableau) -> Option<T::Coeff> {
        self.map.get(tableau).cloned()
    }

    pub fn x(&mut self, qubit: usize) {
        self.map_clifford(|t| t.x(qubit));
    }

    pub fn y(&mut self, qubit: usize) {
        self.map_clifford(|t| t.y(qubit));
    }

    pub fn z(&mut self, qubit: usize) {
        self.map_clifford(|t| t.z(qubit));
    }

    pub fn h(&mut self, qubit: usize) {
        self.map_clifford(|t| t.h(qubit));
    }

    pub fn s(&mut self, qubit: usize) {
        self.map_clifford(|t| t.s(qubit));
    }

    pub fn t(&mut self, qubit: usize)
    where
        T::Coeff: ComplexCoefficient,
    {
        let (a, b) = t_coeffs::<T::Coeff>();
        if self.map.len() >= self.par_threshold {
            let entries: Vec<(Tableau, T::Coeff)> = self.map.drain().collect();
            let next = entries
                .into_par_iter()
                .fold(
                    || T::Map::with_capacity(8),
                    |mut acc, (tableau, coeff)| {
                        let left = tableau.clone();
                        let mut right = tableau;
                        right.z(qubit);
                        add_coeff::<T>(&mut acc, left, coeff.clone() * a.clone());
                        add_coeff::<T>(&mut acc, right, coeff * b.clone());
                        acc
                    },
                )
                .reduce(|| T::Map::with_capacity(0), merge_maps::<T>);
            self.map = next;
        } else {
            let mut next = T::Map::with_capacity(self.map.len() * 2 + 1);
            for (tableau, coeff) in self.map.drain() {
                let left = tableau.clone();
                let mut right = tableau;
                right.z(qubit);
                add_coeff::<T>(&mut next, left, coeff.clone() * a.clone());
                add_coeff::<T>(&mut next, right, coeff * b.clone());
            }
            self.map = next;
        }
    }

    fn map_clifford<F>(&mut self, f: F)
    where
        F: Fn(&mut Tableau) + Sync + Send,
    {
        if self.map.len() >= self.par_threshold {
            let entries: Vec<(Tableau, T::Coeff)> = self.map.drain().collect();
            let next = entries
                .into_par_iter()
                .fold(
                    || T::Map::with_capacity(8),
                    |mut acc, (mut tableau, coeff)| {
                        f(&mut tableau);
                        add_coeff::<T>(&mut acc, tableau, coeff);
                        acc
                    },
                )
                .reduce(|| T::Map::with_capacity(0), merge_maps::<T>);
            self.map = next;
        } else {
            let mut next = T::Map::with_capacity(self.map.len());
            for (mut tableau, coeff) in self.map.drain() {
                f(&mut tableau);
                add_coeff::<T>(&mut next, tableau, coeff);
            }
            self.map = next;
        }
    }
}

impl<T: Config> Clifford for TableauSum<T>
where
    T::Map: TableauMap<T::Coeff, T::BuildHasher> + Send,
{
    fn x(&mut self, index: usize) {
        TableauSum::<T>::x(self, index);
    }

    fn y(&mut self, index: usize) {
        TableauSum::<T>::y(self, index);
    }

    fn z(&mut self, index: usize) {
        TableauSum::<T>::z(self, index);
    }

    fn h(&mut self, index: usize) {
        TableauSum::<T>::h(self, index);
    }

    fn s(&mut self, index: usize) {
        TableauSum::<T>::s(self, index);
    }

    fn cnot(&mut self, _control: usize, _target: usize) {
        unimplemented!("CNOT not yet implemented for TableauSum");
    }

    fn cz(&mut self, _control: usize, _target: usize) {
        unimplemented!("CZ not yet implemented for TableauSum");
    }
}

fn add_coeff<T: Config>(map: &mut T::Map, tableau: Tableau, coeff: T::Coeff)
where
    T::Map: TableauMapAddAssign<T::Coeff, T::BuildHasher>,
{
    map.add_assign(tableau, coeff);
}

fn merge_maps<T: Config>(mut left: T::Map, mut right: T::Map) -> T::Map
where
    T::Map: TableauMapAddAssign<T::Coeff, T::BuildHasher> + TableauMapDrain<T::Coeff>,
{
    for (tableau, coeff) in right.drain() {
        add_coeff::<T>(&mut left, tableau, coeff);
    }
    left
}

pub(crate) fn t_coeffs<C: ComplexCoefficient>() -> (C, C) {
    let angle = std::f64::consts::FRAC_PI_4;
    let cos = angle.cos();
    let sin = angle.sin();
    let w = C::from(cos) + C::from(sin).mul_phase(1);
    let half = C::from(0.5);
    let one = C::from(1.0);
    let a = (one.clone() + w.clone()) * half.clone();
    let b = (one - w) * half;
    (a, b)
}
