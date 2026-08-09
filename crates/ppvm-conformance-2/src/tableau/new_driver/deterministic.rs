// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::complex::Complex64;
use ppvm_tableau_2::{Bitstring, GeneralizedTableau, RowStorage};
use ppvm_traits_2::{
    Clifford, CliffordBatch, CliffordExtensions, CliffordExtensionsBatch, Pauli, RotXY,
    RotationOne, RotationTwo, TGate, U3Gate,
};

use super::super::NewDriver;

type Driven<A, I, H> = NewDriver<GeneralizedTableau<A, I, H>>;

impl<A: RowStorage, I: Bitstring, H> Clifford for Driven<A, I, H> {
    fn x(&mut self, q: usize) {
        self.tab.x(q);
    }
    fn y(&mut self, q: usize) {
        self.tab.y(q);
    }
    fn z(&mut self, q: usize) {
        self.tab.z(q);
    }
    fn h(&mut self, q: usize) {
        self.tab.h(q);
    }
    fn s(&mut self, q: usize) {
        self.tab.s(q);
    }
    fn cnot(&mut self, control: usize, target: usize) {
        self.tab.cnot(control, target);
    }
    fn cz(&mut self, a: usize, b: usize) {
        self.tab.cz(a, b);
    }
}

impl<A: RowStorage, I: Bitstring, H> CliffordExtensions for Driven<A, I, H> {
    fn s_dag(&mut self, q: usize) {
        self.tab.s_dag(q);
    }
    fn sqrt_x(&mut self, q: usize) {
        self.tab.sqrt_x(q);
    }
    fn sqrt_x_dag(&mut self, q: usize) {
        self.tab.sqrt_x_dag(q);
    }
    fn sqrt_y(&mut self, q: usize) {
        self.tab.sqrt_y(q);
    }
    fn sqrt_y_dag(&mut self, q: usize) {
        self.tab.sqrt_y_dag(q);
    }
    fn cy(&mut self, control: usize, target: usize) {
        self.tab.cy(control, target);
    }
}

impl<A: RowStorage, I: Bitstring, H> CliffordBatch for Driven<A, I, H> {
    fn x_many(&mut self, q: &[usize]) {
        self.tab.x_many(q);
    }
    fn y_many(&mut self, q: &[usize]) {
        self.tab.y_many(q);
    }
    fn z_many(&mut self, q: &[usize]) {
        self.tab.z_many(q);
    }
    fn h_many(&mut self, q: &[usize]) {
        self.tab.h_many(q);
    }
    fn s_many(&mut self, q: &[usize]) {
        self.tab.s_many(q);
    }
    fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
        self.tab.cnot_many(pairs);
    }
    fn cz_many(&mut self, pairs: &[(usize, usize)]) {
        self.tab.cz_many(pairs);
    }
}

impl<A: RowStorage, I: Bitstring, H> CliffordExtensionsBatch for Driven<A, I, H> {
    fn s_dag_many(&mut self, q: &[usize]) {
        self.tab.s_dag_many(q);
    }
    fn sqrt_x_many(&mut self, q: &[usize]) {
        self.tab.sqrt_x_many(q);
    }
    fn sqrt_x_dag_many(&mut self, q: &[usize]) {
        self.tab.sqrt_x_dag_many(q);
    }
    fn sqrt_y_many(&mut self, q: &[usize]) {
        self.tab.sqrt_y_many(q);
    }
    fn sqrt_y_dag_many(&mut self, q: &[usize]) {
        self.tab.sqrt_y_dag_many(q);
    }
    fn cy_many(&mut self, pairs: &[(usize, usize)]) {
        self.tab.cy_many(pairs);
    }
}

impl<A: RowStorage, I: Bitstring, H> TGate for Driven<A, I, H> {
    fn t(&mut self, q: usize) {
        self.tab.t(q);
    }
    fn t_dag(&mut self, q: usize) {
        self.tab.t_dag(q);
    }
}

impl<A: RowStorage, I: Bitstring, H> RotationOne<Complex64, f64> for Driven<A, I, H> {
    fn rotate_1(&mut self, axis: Pauli, q: usize, theta: f64) {
        self.tab.rotate_1(axis, q, theta);
    }
}

impl<A: RowStorage, I: Bitstring, H> RotationTwo<Complex64, f64> for Driven<A, I, H> {
    fn rotate_2(&mut self, a_axis: [u8; 2], b_axis: [u8; 2], a: usize, b: usize, theta: f64) {
        self.tab.rotate_2(a_axis, b_axis, a, b, theta);
    }
}

impl<A: RowStorage, I: Bitstring, H> RotXY<Complex64, f64> for Driven<A, I, H> {
    fn r(&mut self, q: usize, axis_angle: f64, theta: f64) {
        self.tab.r(q, axis_angle, theta);
    }
}

impl<A: RowStorage, I: Bitstring, H> U3Gate<Complex64, f64> for Driven<A, I, H> {
    fn u3(&mut self, q: usize, theta: f64, phi: f64, lambda: f64) {
        self.tab.u3(q, theta, phi, lambda);
    }
}
