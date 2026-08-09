// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod sealed {
    pub trait Sealed {}

    #[cfg(feature = "legacy")]
    impl<T, I, C> Sealed for ppvm_tableau_legacy::prelude::GeneralizedTableau<T, I, C>
    where
        T: ppvm_pauli_sum_legacy::prelude::Config,
        <<T as ppvm_pauli_sum_legacy::prelude::Config>::Storage as bitvec::view::BitView>::Store:
            num::PrimInt,
        C: ppvm_tableau_legacy::prelude::SparseVector<num::Complex<T::Coeff>, I> + std::fmt::Debug,
        T::Coeff: num::One
            + num::Zero
            + Clone
            + num::Num
            + num::ToPrimitive
            + std::fmt::Debug
            + std::ops::Mul<f64>
            + PartialOrd<f64>
            + PartialOrd
            + Send
            + Sync,
        num::Complex<T::Coeff>: std::ops::Mul<Output = num::Complex<T::Coeff>>
            + From<num::complex::Complex64>
            + std::ops::MulAssign
            + std::ops::AddAssign
            + num::One
            + num::complex::ComplexFloat
            + Copy,
        I: ppvm_tableau_legacy::prelude::TableauIndex + std::fmt::Debug + Send + Sync,
    {
    }

    #[cfg(feature = "traits-2")]
    impl<A, I, H> Sealed for ppvm_tableau_2::GeneralizedTableau<A, I, H>
    where
        A: ppvm_tableau_2::RowStorage,
        I: ppvm_tableau_2::Bitstring,
    {
    }
}

/// Semantic surface required by the Stim instruction executor.
///
/// This trait is sealed: downstream crates consume it through the generic
/// execution APIs, while backend implementations remain controlled here.
pub trait StimTableau: sealed::Sealed {
    #[doc(hidden)]
    type Config;
    #[doc(hidden)]
    type Index;
    #[doc(hidden)]
    type Coefficients;

    fn reset<R: rand::Rng + ?Sized>(&mut self, q: usize, rng: &mut R);
    fn x(&mut self, q: usize);
    fn y(&mut self, q: usize);
    fn z(&mut self, q: usize);
    fn h(&mut self, q: usize);
    fn s(&mut self, q: usize);
    fn s_dag(&mut self, q: usize);
    fn sqrt_x(&mut self, q: usize);
    fn sqrt_x_dag(&mut self, q: usize);
    fn sqrt_y(&mut self, q: usize);
    fn sqrt_y_dag(&mut self, q: usize);
    fn cnot(&mut self, a: usize, b: usize);
    fn cy(&mut self, a: usize, b: usize);
    fn cz(&mut self, a: usize, b: usize);
    fn x_many(&mut self, q: &[usize]);
    fn y_many(&mut self, q: &[usize]);
    fn z_many(&mut self, q: &[usize]);
    fn h_many(&mut self, q: &[usize]);
    fn s_many(&mut self, q: &[usize]);
    fn s_dag_many(&mut self, q: &[usize]);
    fn sqrt_x_many(&mut self, q: &[usize]);
    fn sqrt_x_dag_many(&mut self, q: &[usize]);
    fn sqrt_y_many(&mut self, q: &[usize]);
    fn sqrt_y_dag_many(&mut self, q: &[usize]);
    fn cnot_many(&mut self, q: &[(usize, usize)]);
    fn cy_many(&mut self, q: &[(usize, usize)]);
    fn cz_many(&mut self, q: &[(usize, usize)]);
    fn t(&mut self, q: usize);
    fn t_dag(&mut self, q: usize);
    fn rx(&mut self, q: usize, theta: f64);
    fn ry(&mut self, q: usize, theta: f64);
    fn rz(&mut self, q: usize, theta: f64);
    fn u3(&mut self, q: usize, theta: f64, phi: f64, lambda: f64);
    fn depolarize1<R: rand::Rng + ?Sized>(&mut self, q: usize, p: f64, rng: &mut R);
    fn depolarize2<R: rand::Rng + ?Sized>(&mut self, a: usize, b: usize, p: f64, rng: &mut R);
    fn pauli_error<R: rand::Rng + ?Sized>(&mut self, q: usize, p: [f64; 3], rng: &mut R);
    fn two_qubit_pauli_error<R: rand::Rng + ?Sized>(
        &mut self,
        a: usize,
        b: usize,
        p: [f64; 15],
        rng: &mut R,
    );
    fn loss_channel<R: rand::Rng + ?Sized>(&mut self, q: usize, p: f64, rng: &mut R);
    fn correlated_loss_channel<R: rand::Rng + ?Sized>(
        &mut self,
        a: usize,
        b: usize,
        p: [f64; 3],
        rng: &mut R,
    );
    fn measure<R: rand::Rng + ?Sized>(&mut self, q: usize, rng: &mut R) -> Option<bool>;
    fn measure_many<R: rand::Rng + ?Sized>(
        &mut self,
        q: &[usize],
        rng: &mut R,
    ) -> Vec<Option<bool>>;
    fn measure_noisy<R: rand::Rng + ?Sized>(
        &mut self,
        q: usize,
        p: f64,
        rng: &mut R,
    ) -> Option<bool>;
    fn flip_with_prob<R: rand::Rng + ?Sized>(&mut self, bit: bool, p: f64, rng: &mut R) -> bool;
    fn measurement_record(&self) -> &[Option<bool>];
    fn append_measurement_record(&mut self, result: Option<bool>);
    fn overwrite_last_measurement_record(&mut self, result: Option<bool>);
}

#[doc(hidden)]
pub struct SelectedBackend;

#[doc(hidden)]
pub trait TableauType<C, I, S> {
    type Tableau: StimTableau<Config = C, Index = I, Coefficients = S>;
}

#[cfg(feature = "legacy")]
impl<C, I, S> TableauType<C, I, S> for SelectedBackend
where
    C: ppvm_pauli_sum_legacy::prelude::Config,
    S: ppvm_tableau_legacy::prelude::SparseVector<num::Complex<C::Coeff>, I>,
    ppvm_tableau_legacy::prelude::GeneralizedTableau<C, I, S>:
        StimTableau<Config = C, Index = I, Coefficients = S>,
{
    type Tableau = ppvm_tableau_legacy::prelude::GeneralizedTableau<C, I, S>;
}

#[cfg(feature = "traits-2")]
impl<C, I, S> TableauType<C, I, S> for SelectedBackend
where
    C: ppvm_tableau_2::RowStorage,
    I: ppvm_tableau_2::Bitstring,
    ppvm_tableau_2::GeneralizedTableau<C, I, S>:
        StimTableau<Config = C, Index = I, Coefficients = S>,
{
    type Tableau = ppvm_tableau_2::GeneralizedTableau<C, I, S>;
}

#[cfg(feature = "legacy")]
mod legacy;
#[cfg(feature = "traits-2")]
mod traits_2;
