// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use num::complex::Complex64;
use ppvm_traits_2::{
    Clifford, CliffordBatch, CliffordExtensions, CliffordExtensionsBatch, Pauli, RotXY,
    RotationOne, RotationTwo, TGate, U3Gate,
};

use super::GeneralizedTableauMixture;
use crate::{Bitstring, RowStorage};

macro_rules! one_qubit {
    ($name:ident) => {
        fn $name(&mut self, qubit: usize) {
            for (tab, _) in &mut self.entries {
                tab.$name(qubit);
            }
            self.mark_dirty();
        }
    };
}

macro_rules! two_qubit {
    ($name:ident) => {
        fn $name(&mut self, a: usize, b: usize) {
            for (tab, _) in &mut self.entries {
                tab.$name(a, b);
            }
            self.mark_dirty();
        }
    };
}

impl<A, I, H> Clifford for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    one_qubit!(x);
    one_qubit!(y);
    one_qubit!(z);
    one_qubit!(h);
    one_qubit!(s);
    two_qubit!(cnot);
    two_qubit!(cz);
}

impl<A, I, H> CliffordExtensions for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    one_qubit!(s_dag);
    one_qubit!(sqrt_x);
    one_qubit!(sqrt_x_dag);
    one_qubit!(sqrt_y);
    one_qubit!(sqrt_y_dag);
    two_qubit!(cy);
}

impl<A, I, H> CliffordBatch for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
}

impl<A, I, H> CliffordExtensionsBatch for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
}

impl<A, I, H> TGate for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    one_qubit!(t);
    one_qubit!(t_dag);
}

impl<A, I, H> RotationOne<Complex64, f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn rotate_1(&mut self, axis: Pauli, qubit: usize, theta: f64) {
        for (tab, _) in &mut self.entries {
            tab.rotate_1(axis, qubit, theta);
        }
        self.mark_dirty();
    }
}

impl<A, I, H> RotationTwo<Complex64, f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn rotate_2(&mut self, axis_a: [u8; 2], axis_b: [u8; 2], a: usize, b: usize, theta: f64) {
        for (tab, _) in &mut self.entries {
            tab.rotate_2(axis_a, axis_b, a, b, theta);
        }
        self.mark_dirty();
    }
}

impl<A, I, H> RotXY<Complex64, f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn r(&mut self, qubit: usize, axis_angle: f64, theta: f64) {
        for (tab, _) in &mut self.entries {
            tab.r(qubit, axis_angle, theta);
        }
        self.mark_dirty();
    }
}

impl<A, I, H> U3Gate<Complex64, f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn u3(&mut self, qubit: usize, theta: f64, phi: f64, lambda: f64) {
        for (tab, _) in &mut self.entries {
            tab.u3(qubit, theta, phi, lambda);
        }
        self.mark_dirty();
    }
}
