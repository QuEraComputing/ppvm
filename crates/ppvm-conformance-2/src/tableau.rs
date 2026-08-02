// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Phase-4 tableau differential harness: matched OLD (`ppvm-tableau`) and NEW
//! (`ppvm-tableau-2`) engines driven through **one** circuit description, so a
//! differential test or a same-build benchmark can never accidentally run two
//! different gate sequences.
//!
//! The two engines have incompatible trait towers (`ppvm-traits` vs
//! `ppvm-traits-2`) and incompatible generic parameters (`Config` vs a bare
//! storage array), so the shared surface is expressed as the object-safe-ish
//! [`Driver`] trait below and implemented once per engine. Every workload in
//! this module (MSD-85q naive and fused, the rot2 brickwork, the fused-T-gate
//! circuit, the CNOT-chain scaling sweep, …) is generic over [`Driver`] and is
//! therefore *literally* the same code on both sides.
//!
//! # Matched configurations
//!
//! Correctness is storage-independent, but a perf ratio is not: `BitArray<u64>`
//! and `BitArray<[u8; 8]>` differ by a few percent in codegen, which would fold
//! into an engine-to-engine ratio. Each alias pair below therefore pins the
//! **same** storage width and the same coefficient/index types on both sides:
//!
//! | workload | old | new | storage |
//! |:--|:--|:--|:--|
//! | MSD-85q, fused-T | `GeneralizedTableau<Byte8F64<2>, u128>` | `GeneralizedTableau<[usize; 2], u128>` | `[usize; 2]` |
//! | rot2 brickwork | `GeneralizedTableau<ByteFxHashF64<8>, usize>` | `GeneralizedTableau<[u8; 8], usize>` | `[u8; 8]` |
//! | scaling sweep | `GeneralizedTableau<Byte8F64<2>, usize>` | `GeneralizedTableau<[usize; 2], usize>` | `[usize; 2]` |
//!
//! (`Byte8F64<N>` is `[usize; N]`-backed; `ByteFxHashF64<N>` is `[u8; N]`-backed.
//! Both use `FxHash`, which is also the new crate's default `H`.)

use num::complex::Complex64;

// --- OLD engine ------------------------------------------------------------
use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_tableau::data::GeneralizedTableau as OldGeneralizedTableau;
use ppvm_tableau::measure_all::LossyMeasureAll;

// --- NEW engine ------------------------------------------------------------
use ppvm_tableau_2::GeneralizedTableau as NewGeneralizedTableau;

// Every gate call below is written in fully-qualified (UFCS) form: the two trait
// towers export same-named traits (`Clifford`, `TGate`, …) *and* the concrete
// types carry same-named inherent methods, so importing either set into scope
// makes the calls ambiguous. UFCS names the exact trait on the exact side.

/// OLD 85-qubit MSD / fused-T configuration (`[usize; 2]` storage, `u128` index).
pub type OldWide = OldGeneralizedTableau<Byte8F64<2>, u128>;
/// NEW 85-qubit MSD / fused-T configuration, storage-matched to [`OldWide`].
pub type NewWide = NewGeneralizedTableau<[usize; 2], u128>;

/// OLD rot2-brickwork configuration (`[u8; 8]` storage, `usize` index).
pub type OldNarrow = OldGeneralizedTableau<ByteFxHashF64<8>, usize>;
/// NEW rot2-brickwork configuration, storage-matched to [`OldNarrow`].
pub type NewNarrow = NewGeneralizedTableau<[u8; 8], usize>;

/// OLD scaling-sweep configuration (`[usize; 2]` storage, `usize` index).
pub type OldScaling = OldGeneralizedTableau<Byte8F64<2>, usize>;
/// NEW scaling-sweep configuration, storage-matched to [`OldScaling`].
pub type NewScaling = NewGeneralizedTableau<[usize; 2], usize>;

/// A `(x-plane, z-plane, phase)` snapshot of one tableau row, normalized to
/// `Vec<u64>` planes so the OLD `[usize; N]` / `[u8; N]` rows and the NEW ones
/// compare element-for-element regardless of the storage element type.
pub type RowSnapshot = (Vec<u64>, Vec<u64>, u8);

/// The gate/measurement surface both engines share, so a workload can be written
/// **once** and replayed identically on each.
///
/// Deliberately narrow: only what the integration baseline's workloads need,
/// plus the observation hooks (`rows`, `coeffs`, `record`) the differential
/// assertions read. Everything returns plain data (`Vec<u64>`, `u128`, …) so a
/// test can compare across the two type towers.
pub trait Driver {
    /// A fresh `|0…0⟩` state with a deterministic RNG seed.
    fn new_seeded(n_qubits: usize, threshold: f64, seed: u64) -> Self;
    /// A fresh `|0…0⟩` state seeded from OS entropy.
    fn new_entropy(n_qubits: usize, threshold: f64) -> Self;
    /// Clone the state and reseed the RNG (`None` ⇒ OS entropy).
    fn fork(&self, seed: Option<u64>) -> Self;

    /// Number of qubits.
    fn n_qubits(&self) -> usize;

    // --- Clifford ---------------------------------------------------------
    /// Pauli `X`.
    fn x(&mut self, q: usize);
    /// Pauli `Y`.
    fn y(&mut self, q: usize);
    /// Pauli `Z`.
    fn z(&mut self, q: usize);
    /// Hadamard.
    fn h(&mut self, q: usize);
    /// Phase gate `S`.
    fn s(&mut self, q: usize);
    /// `S†`.
    fn s_dag(&mut self, q: usize);
    /// `√X`.
    fn sqrt_x(&mut self, q: usize);
    /// `√X†`.
    fn sqrt_x_dag(&mut self, q: usize);
    /// `√Y`.
    fn sqrt_y(&mut self, q: usize);
    /// `√Y†`.
    fn sqrt_y_dag(&mut self, q: usize);
    /// CNOT.
    fn cnot(&mut self, c: usize, t: usize);
    /// CZ.
    fn cz(&mut self, a: usize, b: usize);
    /// CY.
    fn cy(&mut self, c: usize, t: usize);
    /// `zcx` (X-controlled X).
    fn zcx(&mut self, c: usize, t: usize);
    /// `zcy`.
    fn zcy(&mut self, c: usize, t: usize);
    /// `zcz`.
    fn zcz(&mut self, a: usize, b: usize);

    // --- batched Clifford --------------------------------------------------
    /// Batched `H`.
    fn h_many(&mut self, qs: &[usize]);
    /// Batched `X`.
    fn x_many(&mut self, qs: &[usize]);
    /// Batched `S`.
    fn s_many(&mut self, qs: &[usize]);
    /// Batched `√X`.
    fn sqrt_x_many(&mut self, qs: &[usize]);
    /// Batched `√X†`.
    fn sqrt_x_dag_many(&mut self, qs: &[usize]);
    /// Batched `√Y`.
    fn sqrt_y_many(&mut self, qs: &[usize]);
    /// Batched `√Y†`.
    fn sqrt_y_dag_many(&mut self, qs: &[usize]);
    /// Batched CZ.
    fn cz_many(&mut self, pairs: &[(usize, usize)]);
    /// Batched CNOT.
    fn cnot_many(&mut self, pairs: &[(usize, usize)]);
    /// Fused constant-offset CZ block.
    fn cz_block(&mut self, control_base: usize, target_base: usize, count: usize);
    /// Fused same-word CZ pairs.
    fn cz_block_pairs(&mut self, base: usize, offset: usize, count: usize);

    // --- non-Clifford ------------------------------------------------------
    /// `T`.
    fn t(&mut self, q: usize);
    /// `T†`.
    fn t_dag(&mut self, q: usize);
    /// `RX(θ)`.
    fn rx(&mut self, q: usize, theta: f64);
    /// `RY(θ)`.
    fn ry(&mut self, q: usize, theta: f64);
    /// `RZ(θ)`.
    fn rz(&mut self, q: usize, theta: f64);
    /// `RXX(θ)`.
    fn rxx(&mut self, a: usize, b: usize, theta: f64);
    /// `RYY(θ)`.
    fn ryy(&mut self, a: usize, b: usize, theta: f64);
    /// `RZZ(θ)`.
    fn rzz(&mut self, a: usize, b: usize, theta: f64);

    // --- measurement / channels -------------------------------------------
    /// Measure one qubit in the Z basis.
    fn measure(&mut self, q: usize) -> Option<bool>;
    /// Measure every qubit, ascending.
    fn measure_all(&mut self) -> Vec<Option<bool>>;
    /// Measure the given qubits, in the caller's order.
    fn measure_many(&mut self, qs: &[usize]) -> Vec<Option<bool>>;
    /// Measure with readout noise.
    fn measure_noisy(&mut self, q: usize, flip_prob: f64) -> Option<bool>;
    /// Measure-and-reset to `|0⟩`.
    fn reset(&mut self, q: usize);
    /// Single-qubit depolarizing channel.
    fn depolarize1(&mut self, q: usize, p: f64);
    /// Two-qubit depolarizing channel.
    fn depolarize2(&mut self, a: usize, b: usize, p: f64);
    /// Photon-loss channel.
    fn loss_channel(&mut self, q: usize, p: f64);
    /// State-dependent loss channel.
    fn asymmetric_loss_channel(&mut self, q: usize, p0: f64, p1: f64);
    /// Correlated two-qubit loss channel.
    fn correlated_loss_channel(&mut self, a: usize, b: usize, p: [f64; 3]);
    /// Clear the loss flag on a qubit.
    fn reset_loss_channel(&mut self, q: usize);

    // --- observation -------------------------------------------------------
    /// All `2n` rows as `(x-plane, z-plane, phase)`, destabilizers first.
    fn rows(&self) -> Vec<RowSnapshot>;
    /// The amplitude vector **in stored order**, as `(index, value)`.
    fn coeffs(&self) -> Vec<(u128, Complex64)>;
    /// The amplitude vector sorted ascending by index.
    fn coeffs_sorted(&self) -> Vec<(u128, Complex64)> {
        let mut v = self.coeffs();
        v.sort_by_key(|e| e.0);
        v
    }
    /// Number of stored amplitudes.
    fn n_coeffs(&self) -> usize;
    /// The measurement record so far.
    fn record(&self) -> Vec<Option<bool>>;
    /// Per-qubit loss flags.
    fn lost(&self) -> Vec<bool>;
    /// `⟨Z_q⟩`, non-destructively.
    fn z_expectation(&self, q: usize) -> f64;
    /// `⟨ψ|P|ψ⟩` for the Pauli string `word` (e.g. `"ZZ"`).
    fn expectation_str(&self, word: &str) -> f64;
    /// The next `f64` from the state's RNG (drains the stream — test-only, used
    /// to prove two runs left the stream in the same place).
    fn peek_rng_f64(&mut self) -> f64;
}

// ===========================================================================
// OLD impl
// ===========================================================================

macro_rules! impl_old_driver {
    ($ty:ty, $cfg:ty, $idx:ty, $store:ty) => {
        impl Driver for $ty {
            fn new_seeded(n_qubits: usize, threshold: f64, seed: u64) -> Self {
                <$ty>::new_with_seed(n_qubits, threshold, seed)
            }
            fn new_entropy(n_qubits: usize, threshold: f64) -> Self {
                <$ty>::new(n_qubits, threshold)
            }
            fn fork(&self, seed: Option<u64>) -> Self {
                OldGeneralizedTableau::fork(self, seed)
            }
            fn n_qubits(&self) -> usize {
                self.tableau.n_qubits
            }

            fn x(&mut self, q: usize) {
                ppvm_traits::traits::Clifford::x(self, q)
            }
            fn y(&mut self, q: usize) {
                ppvm_traits::traits::Clifford::y(self, q)
            }
            fn z(&mut self, q: usize) {
                ppvm_traits::traits::Clifford::z(self, q)
            }
            fn h(&mut self, q: usize) {
                ppvm_traits::traits::Clifford::h(self, q)
            }
            fn s(&mut self, q: usize) {
                ppvm_traits::traits::Clifford::s(self, q)
            }
            fn s_dag(&mut self, q: usize) {
                ppvm_traits::traits::CliffordExtensions::s_dag(self, q)
            }
            fn sqrt_x(&mut self, q: usize) {
                ppvm_traits::traits::CliffordExtensions::sqrt_x(self, q)
            }
            fn sqrt_x_dag(&mut self, q: usize) {
                ppvm_traits::traits::CliffordExtensions::sqrt_x_dag(self, q)
            }
            fn sqrt_y(&mut self, q: usize) {
                ppvm_traits::traits::CliffordExtensions::sqrt_y(self, q)
            }
            fn sqrt_y_dag(&mut self, q: usize) {
                ppvm_traits::traits::CliffordExtensions::sqrt_y_dag(self, q)
            }
            fn cnot(&mut self, c: usize, t: usize) {
                ppvm_traits::traits::Clifford::cnot(self, c, t)
            }
            fn cz(&mut self, a: usize, b: usize) {
                ppvm_traits::traits::Clifford::cz(self, a, b)
            }
            fn cy(&mut self, c: usize, t: usize) {
                ppvm_traits::traits::CliffordExtensions::cy(self, c, t)
            }
            fn zcx(&mut self, c: usize, t: usize) {
                ppvm_traits::traits::Clifford::zcx(self, c, t)
            }
            fn zcy(&mut self, c: usize, t: usize) {
                ppvm_traits::traits::CliffordExtensions::zcy(self, c, t)
            }
            fn zcz(&mut self, a: usize, b: usize) {
                ppvm_traits::traits::Clifford::zcz(self, a, b)
            }

            fn h_many(&mut self, qs: &[usize]) {
                ppvm_traits::traits::CliffordBatch::h_many(self, qs)
            }
            fn x_many(&mut self, qs: &[usize]) {
                ppvm_traits::traits::CliffordBatch::x_many(self, qs)
            }
            fn s_many(&mut self, qs: &[usize]) {
                ppvm_traits::traits::CliffordBatch::s_many(self, qs)
            }
            fn sqrt_x_many(&mut self, qs: &[usize]) {
                ppvm_traits::traits::CliffordExtensionsBatch::sqrt_x_many(self, qs)
            }
            fn sqrt_x_dag_many(&mut self, qs: &[usize]) {
                ppvm_traits::traits::CliffordExtensionsBatch::sqrt_x_dag_many(self, qs)
            }
            fn sqrt_y_many(&mut self, qs: &[usize]) {
                ppvm_traits::traits::CliffordExtensionsBatch::sqrt_y_many(self, qs)
            }
            fn sqrt_y_dag_many(&mut self, qs: &[usize]) {
                ppvm_traits::traits::CliffordExtensionsBatch::sqrt_y_dag_many(self, qs)
            }
            fn cz_many(&mut self, pairs: &[(usize, usize)]) {
                ppvm_traits::traits::CliffordBatch::cz_many(self, pairs)
            }
            fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
                ppvm_traits::traits::CliffordBatch::cnot_many(self, pairs)
            }
            fn cz_block(&mut self, control_base: usize, target_base: usize, count: usize) {
                OldGeneralizedTableau::cz_block(self, control_base, target_base, count)
            }
            fn cz_block_pairs(&mut self, base: usize, offset: usize, count: usize) {
                OldGeneralizedTableau::cz_block_pairs(self, base, offset, count)
            }

            fn t(&mut self, q: usize) {
                ppvm_traits::traits::TGate::t(self, q)
            }
            fn t_dag(&mut self, q: usize) {
                ppvm_traits::traits::TGate::t_dag(self, q)
            }
            fn rx(&mut self, q: usize, theta: f64) {
                ppvm_traits::traits::RotationOne::rx(self, q, theta)
            }
            fn ry(&mut self, q: usize, theta: f64) {
                ppvm_traits::traits::RotationOne::ry(self, q, theta)
            }
            fn rz(&mut self, q: usize, theta: f64) {
                ppvm_traits::traits::RotationOne::rz(self, q, theta)
            }
            fn rxx(&mut self, a: usize, b: usize, theta: f64) {
                ppvm_traits::traits::RotationTwo::rxx(self, a, b, theta)
            }
            fn ryy(&mut self, a: usize, b: usize, theta: f64) {
                ppvm_traits::traits::RotationTwo::ryy(self, a, b, theta)
            }
            fn rzz(&mut self, a: usize, b: usize, theta: f64) {
                ppvm_traits::traits::RotationTwo::rzz(self, a, b, theta)
            }

            fn measure(&mut self, q: usize) -> Option<bool> {
                ppvm_traits::traits::LossyMeasure::measure(self, q)
            }
            fn measure_all(&mut self) -> Vec<Option<bool>> {
                <Self as LossyMeasureAll>::measure_all(self)
            }
            fn measure_many(&mut self, qs: &[usize]) -> Vec<Option<bool>> {
                ppvm_traits::traits::LossyMeasure::measure_many(self, qs)
            }
            fn measure_noisy(&mut self, q: usize, flip_prob: f64) -> Option<bool> {
                OldGeneralizedTableau::measure_noisy(self, q, flip_prob)
            }
            fn reset(&mut self, q: usize) {
                ppvm_traits::traits::Reset::reset(self, q)
            }
            fn depolarize1(&mut self, q: usize, p: f64) {
                ppvm_traits::traits::Depolarizing::<$cfg>::depolarize1(self, q, p)
            }
            fn depolarize2(&mut self, a: usize, b: usize, p: f64) {
                ppvm_traits::traits::Depolarizing2::<$cfg>::depolarize2(self, a, b, p)
            }
            fn loss_channel(&mut self, q: usize, p: f64) {
                ppvm_traits::traits::LossChannel::<$cfg>::loss_channel(self, q, p)
            }
            fn asymmetric_loss_channel(&mut self, q: usize, p0: f64, p1: f64) {
                ppvm_traits::traits::AsymmetricLossChannel::<$cfg>::asymmetric_loss_channel(
                    self, q, p0, p1,
                )
            }
            fn correlated_loss_channel(&mut self, a: usize, b: usize, p: [f64; 3]) {
                ppvm_traits::traits::CorrelatedLossChannel::<$cfg>::correlated_loss_channel(
                    self, a, b, p,
                )
            }
            fn reset_loss_channel(&mut self, q: usize) {
                ppvm_traits::traits::ResetLossChannel::<$cfg>::reset_loss_channel(self, q)
            }

            fn rows(&self) -> Vec<RowSnapshot> {
                self.tableau
                    .data
                    .iter()
                    .map(|r| {
                        (
                            r.word.xbits.data.iter().map(|&w| w as u64).collect(),
                            r.word.zbits.data.iter().map(|&w| w as u64).collect(),
                            r.phase,
                        )
                    })
                    .collect()
            }
            fn coeffs(&self) -> Vec<(u128, Complex64)> {
                self.coefficients
                    .iter()
                    .map(|&(c, i)| (i as u128, c))
                    .collect()
            }
            fn n_coeffs(&self) -> usize {
                self.coefficients.len()
            }
            fn record(&self) -> Vec<Option<bool>> {
                self.measurement_record.clone()
            }
            fn lost(&self) -> Vec<bool> {
                self.is_lost.clone()
            }
            fn z_expectation(&self, q: usize) -> f64 {
                OldGeneralizedTableau::z_expectation(self, q)
            }
            fn expectation_str(&self, word: &str) -> f64 {
                let w: ppvm_pauli_word::word::PauliWord<$store> = word.into();
                OldGeneralizedTableau::expectation(self, &w)
            }
            fn peek_rng_f64(&mut self) -> f64 {
                OldGeneralizedTableau::bernoulli(self, 0.5) as u8 as f64
            }
        }
    };
}

impl_old_driver!(OldWide, Byte8F64<2>, u128, [usize; 2]);
impl_old_driver!(OldNarrow, ByteFxHashF64<8>, usize, [u8; 8]);
impl_old_driver!(OldScaling, Byte8F64<2>, usize, [usize; 2]);

// ===========================================================================
// NEW impl
// ===========================================================================

macro_rules! impl_new_driver {
    ($ty:ty, $store:ty, $idx:ty) => {
        impl Driver for $ty {
            fn new_seeded(n_qubits: usize, threshold: f64, seed: u64) -> Self {
                <$ty>::new_with_seed(n_qubits, threshold, seed)
            }
            fn new_entropy(n_qubits: usize, threshold: f64) -> Self {
                <$ty>::new(n_qubits, threshold)
            }
            fn fork(&self, seed: Option<u64>) -> Self {
                NewGeneralizedTableau::fork(self, seed)
            }
            fn n_qubits(&self) -> usize {
                NewGeneralizedTableau::n_qubits(self)
            }

            fn x(&mut self, q: usize) {
                ppvm_traits_2::Clifford::x(self, q)
            }
            fn y(&mut self, q: usize) {
                ppvm_traits_2::Clifford::y(self, q)
            }
            fn z(&mut self, q: usize) {
                ppvm_traits_2::Clifford::z(self, q)
            }
            fn h(&mut self, q: usize) {
                ppvm_traits_2::Clifford::h(self, q)
            }
            fn s(&mut self, q: usize) {
                ppvm_traits_2::Clifford::s(self, q)
            }
            fn s_dag(&mut self, q: usize) {
                ppvm_traits_2::CliffordExtensions::s_dag(self, q)
            }
            fn sqrt_x(&mut self, q: usize) {
                ppvm_traits_2::CliffordExtensions::sqrt_x(self, q)
            }
            fn sqrt_x_dag(&mut self, q: usize) {
                ppvm_traits_2::CliffordExtensions::sqrt_x_dag(self, q)
            }
            fn sqrt_y(&mut self, q: usize) {
                ppvm_traits_2::CliffordExtensions::sqrt_y(self, q)
            }
            fn sqrt_y_dag(&mut self, q: usize) {
                ppvm_traits_2::CliffordExtensions::sqrt_y_dag(self, q)
            }
            fn cnot(&mut self, c: usize, t: usize) {
                ppvm_traits_2::Clifford::cnot(self, c, t)
            }
            fn cz(&mut self, a: usize, b: usize) {
                ppvm_traits_2::Clifford::cz(self, a, b)
            }
            fn cy(&mut self, c: usize, t: usize) {
                ppvm_traits_2::CliffordExtensions::cy(self, c, t)
            }
            fn zcx(&mut self, c: usize, t: usize) {
                ppvm_traits_2::Clifford::zcx(self, c, t)
            }
            fn zcy(&mut self, c: usize, t: usize) {
                ppvm_traits_2::CliffordExtensions::zcy(self, c, t)
            }
            fn zcz(&mut self, a: usize, b: usize) {
                ppvm_traits_2::Clifford::zcz(self, a, b)
            }

            fn h_many(&mut self, qs: &[usize]) {
                ppvm_traits_2::CliffordBatch::h_many(self, qs)
            }
            fn x_many(&mut self, qs: &[usize]) {
                ppvm_traits_2::CliffordBatch::x_many(self, qs)
            }
            fn s_many(&mut self, qs: &[usize]) {
                ppvm_traits_2::CliffordBatch::s_many(self, qs)
            }
            fn sqrt_x_many(&mut self, qs: &[usize]) {
                ppvm_traits_2::CliffordExtensionsBatch::sqrt_x_many(self, qs)
            }
            fn sqrt_x_dag_many(&mut self, qs: &[usize]) {
                ppvm_traits_2::CliffordExtensionsBatch::sqrt_x_dag_many(self, qs)
            }
            fn sqrt_y_many(&mut self, qs: &[usize]) {
                ppvm_traits_2::CliffordExtensionsBatch::sqrt_y_many(self, qs)
            }
            fn sqrt_y_dag_many(&mut self, qs: &[usize]) {
                ppvm_traits_2::CliffordExtensionsBatch::sqrt_y_dag_many(self, qs)
            }
            fn cz_many(&mut self, pairs: &[(usize, usize)]) {
                ppvm_traits_2::CliffordBatch::cz_many(self, pairs)
            }
            fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
                ppvm_traits_2::CliffordBatch::cnot_many(self, pairs)
            }
            fn cz_block(&mut self, control_base: usize, target_base: usize, count: usize) {
                NewGeneralizedTableau::cz_block(self, control_base, target_base, count)
            }
            fn cz_block_pairs(&mut self, base: usize, offset: usize, count: usize) {
                NewGeneralizedTableau::cz_block_pairs(self, base, offset, count)
            }

            fn t(&mut self, q: usize) {
                ppvm_traits_2::TGate::t(self, q)
            }
            fn t_dag(&mut self, q: usize) {
                ppvm_traits_2::TGate::t_dag(self, q)
            }
            fn rx(&mut self, q: usize, theta: f64) {
                ppvm_traits_2::RotationOne::<Complex64, f64>::rx(self, q, theta)
            }
            fn ry(&mut self, q: usize, theta: f64) {
                ppvm_traits_2::RotationOne::<Complex64, f64>::ry(self, q, theta)
            }
            fn rz(&mut self, q: usize, theta: f64) {
                ppvm_traits_2::RotationOne::<Complex64, f64>::rz(self, q, theta)
            }
            fn rxx(&mut self, a: usize, b: usize, theta: f64) {
                ppvm_traits_2::RotationTwo::<Complex64, f64>::rxx(self, a, b, theta)
            }
            fn ryy(&mut self, a: usize, b: usize, theta: f64) {
                ppvm_traits_2::RotationTwo::<Complex64, f64>::ryy(self, a, b, theta)
            }
            fn rzz(&mut self, a: usize, b: usize, theta: f64) {
                ppvm_traits_2::RotationTwo::<Complex64, f64>::rzz(self, a, b, theta)
            }

            fn measure(&mut self, q: usize) -> Option<bool> {
                ppvm_traits_2::Measure::measure(self, q)
            }
            fn measure_all(&mut self) -> Vec<Option<bool>> {
                NewGeneralizedTableau::measure_all(self)
            }
            fn measure_many(&mut self, qs: &[usize]) -> Vec<Option<bool>> {
                ppvm_traits_2::Measure::measure_many(self, qs)
            }
            fn measure_noisy(&mut self, q: usize, flip_prob: f64) -> Option<bool> {
                NewGeneralizedTableau::measure_noisy(self, q, flip_prob)
            }
            fn reset(&mut self, q: usize) {
                ppvm_traits_2::Reset::reset(self, q)
            }
            fn depolarize1(&mut self, q: usize, p: f64) {
                ppvm_traits_2::Depolarizing::<f64>::depolarize1(self, q, p)
            }
            fn depolarize2(&mut self, a: usize, b: usize, p: f64) {
                ppvm_traits_2::Depolarizing2::<f64>::depolarize2(self, a, b, p)
            }
            fn loss_channel(&mut self, q: usize, p: f64) {
                ppvm_traits_2::LossChannel::<f64>::loss_channel(self, q, p)
            }
            fn asymmetric_loss_channel(&mut self, q: usize, p0: f64, p1: f64) {
                ppvm_traits_2::AsymmetricLossChannel::<f64>::asymmetric_loss_channel(
                    self, q, p0, p1,
                )
            }
            fn correlated_loss_channel(&mut self, a: usize, b: usize, p: [f64; 3]) {
                ppvm_traits_2::CorrelatedLossChannel::<f64>::correlated_loss_channel(self, a, b, p)
            }
            fn reset_loss_channel(&mut self, q: usize) {
                ppvm_traits_2::ResetLossChannel::reset_loss_channel(self, q)
            }

            fn rows(&self) -> Vec<RowSnapshot> {
                self.tableau
                    .rows()
                    .map(|(x, z, p)| {
                        (
                            x.iter().map(|&w| w as u64).collect(),
                            z.iter().map(|&w| w as u64).collect(),
                            p,
                        )
                    })
                    .collect()
            }
            fn coeffs(&self) -> Vec<(u128, Complex64)> {
                self.coefficients
                    .iter()
                    .map(|&(c, i)| (i as u128, c))
                    .collect()
            }
            fn n_coeffs(&self) -> usize {
                self.coefficients.len()
            }
            fn record(&self) -> Vec<Option<bool>> {
                self.measurement_record.clone()
            }
            fn lost(&self) -> Vec<bool> {
                self.is_lost.clone()
            }
            fn z_expectation(&self, q: usize) -> f64 {
                NewGeneralizedTableau::z_expectation(self, q)
            }
            fn expectation_str(&self, word: &str) -> f64 {
                let w: ppvm_pauli_word_2::PauliWord<$store> = word.into();
                NewGeneralizedTableau::expectation(self, &w)
            }
            fn peek_rng_f64(&mut self) -> f64 {
                NewGeneralizedTableau::bernoulli(self, 0.5) as u8 as f64
            }
        }
    };
}

impl_new_driver!(NewWide, [usize; 2], u128);
impl_new_driver!(NewNarrow, [u8; 8], usize);
impl_new_driver!(NewScaling, [usize; 2], usize);

// ===========================================================================
// Integration-baseline workloads — written once, replayed on both engines
// ===========================================================================

/// The 17-qubit MSD encoder (`encode` in `benches/tableau-msd.rs`), gate for
/// gate.
pub fn msd_encode<D: Driver>(tab: &mut D, q: &[usize]) {
    assert_eq!(q.len(), 17, "only the 17-qubit code block is used");
    for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16] {
        tab.sqrt_y(q[i]);
    }
    for [i, j] in [[1, 3], [7, 10], [12, 14], [13, 16]] {
        tab.cz(q[i], q[j]);
    }
    for i in [7, 16] {
        tab.sqrt_y_dag(q[i]);
    }
    for [i, j] in [[4, 7], [8, 10], [11, 14], [15, 16]] {
        tab.cz(q[i], q[j]);
    }
    for i in [4, 10, 14, 16] {
        tab.sqrt_y_dag(q[i]);
    }
    for [i, j] in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]] {
        tab.cz(q[i], q[j]);
    }
    for i in [3, 6, 9, 10, 12, 13] {
        tab.sqrt_y(q[i]);
    }
    for [i, j] in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]] {
        tab.cz(q[i], q[j]);
    }
    for i in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14] {
        tab.sqrt_y(q[i]);
    }
    for [i, j] in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]] {
        tab.cz(q[i], q[j]);
    }
    for i in [0, 2, 5, 6, 8, 10, 12] {
        tab.sqrt_y_dag(q[i]);
    }
}

/// The batched ("fused") 17-qubit MSD encoder from
/// `benches/tableau-msd-fused.rs`.
pub fn msd_encode_fused<D: Driver>(tab: &mut D, q: &[usize]) {
    assert_eq!(q.len(), 17);
    let pick = |ids: &[usize]| -> Vec<usize> { ids.iter().map(|&i| q[i]).collect() };
    let pairs = |ps: &[[usize; 2]]| -> Vec<(usize, usize)> {
        ps.iter().map(|&[i, j]| (q[i], q[j])).collect()
    };
    tab.sqrt_y_many(&pick(&[
        0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    ]));
    tab.cz_many(&pairs(&[[1, 3], [7, 10], [12, 14], [13, 16]]));
    tab.sqrt_y_dag_many(&pick(&[7, 16]));
    tab.cz_many(&pairs(&[[4, 7], [8, 10], [11, 14], [15, 16]]));
    tab.sqrt_y_dag_many(&pick(&[4, 10, 14, 16]));
    tab.cz_many(&pairs(&[[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]]));
    tab.sqrt_y_many(&pick(&[3, 6, 9, 10, 12, 13]));
    tab.cz_many(&pairs(&[[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]]));
    tab.sqrt_y_many(&pick(&[1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14]));
    tab.cz_many(&pairs(&[[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]]));
    tab.sqrt_y_dag_many(&pick(&[0, 2, 5, 6, 8, 10, 12]));
}

/// Number of qubits in one MSD code block.
pub const MSD_BLOCK: usize = 17;
/// Total qubits in the 5-block MSD workload.
pub const MSD_QUBITS: usize = MSD_BLOCK * 5;
/// The MSD workload's truncation threshold.
pub const MSD_THRESHOLD: f64 = 1e-10;

/// Build the naive (per-gate) 85-qubit MSD state, **without** measuring.
///
/// Integration baseline #1's Clifford + T portion, gate for gate from
/// `crates/ppvm-tableau/benches/tableau-msd.rs::msd_func`.
pub fn msd_state<D: Driver>(seed: Option<u64>) -> D {
    let n = MSD_QUBITS;
    let mut tab: D = match seed {
        Some(s) => D::new_seeded(n, MSD_THRESHOLD, s),
        None => D::new_entropy(n, MSD_THRESHOLD),
    };
    let addrs: Vec<usize> = (0..n).collect();
    let ql: Vec<&[usize]> = addrs.chunks_exact(MSD_BLOCK).collect();

    for q in ql.iter() {
        let enc = q[7];
        tab.h(enc);
        tab.t(enc);
        msd_encode(&mut tab, q);
    }
    for i in [0, 1, 4] {
        for &q in ql[i] {
            tab.sqrt_x(q);
        }
    }
    for (&c, &t) in ql[0].iter().zip(ql[1]) {
        tab.cz(c, t);
    }
    for (&c, &t) in ql[2].iter().zip(ql[3]) {
        tab.cz(c, t);
    }
    for &q in ql[0] {
        tab.sqrt_y(q);
    }
    for &q in ql[3] {
        tab.sqrt_y(q);
    }
    for (&c, &t) in ql[0].iter().zip(ql[2]) {
        tab.cz(c, t);
    }
    for (&c, &t) in ql[3].iter().zip(ql[4]) {
        tab.cz(c, t);
    }
    for &q in ql[0] {
        tab.sqrt_x_dag(q);
    }
    for (&c, &t) in ql[0].iter().zip(ql[4]) {
        tab.cz(c, t);
    }
    for (&c, &t) in ql[1].iter().zip(ql[3]) {
        tab.cz(c, t);
    }
    for block in ql.iter().take(5) {
        for &q in *block {
            tab.sqrt_x_dag(q);
        }
    }
    tab
}

/// Full naive MSD workload: build, then `measure(i)` for `i` in `0..85`,
/// collected into an 85-character bitstring (`msd_func`).
pub fn msd_bitstring<D: Driver>(seed: Option<u64>) -> String {
    let mut tab: D = msd_state(seed);
    (0..MSD_QUBITS)
        .map(|i| if tab.measure(i).unwrap() { '1' } else { '0' })
        .collect()
}

/// Build the **fused** 85-qubit MSD state (batched gates + `cz_block`), without
/// measuring — integration baseline #2's Clifford portion.
pub fn msd_state_fused<D: Driver>(seed: Option<u64>) -> D {
    let n = MSD_QUBITS;
    let mut tab: D = match seed {
        Some(s) => D::new_seeded(n, MSD_THRESHOLD, s),
        None => D::new_entropy(n, MSD_THRESHOLD),
    };
    let addrs: Vec<usize> = (0..n).collect();
    let ql: Vec<&[usize]> = addrs.chunks_exact(MSD_BLOCK).collect();

    for q in ql.iter() {
        let enc = q[7];
        tab.h(enc);
        tab.t(enc);
        msd_encode_fused(&mut tab, q);
    }
    tab.sqrt_x_many(ql[0]);
    tab.sqrt_x_many(ql[1]);
    tab.sqrt_x_many(ql[4]);

    tab.cz_block(ql[0][0], ql[1][0], MSD_BLOCK);
    tab.cz_block(ql[2][0], ql[3][0], MSD_BLOCK);

    tab.sqrt_y_many(ql[0]);
    tab.sqrt_y_many(ql[3]);

    tab.cz_block(ql[0][0], ql[2][0], MSD_BLOCK);
    tab.cz_block(ql[3][0], ql[4][0], MSD_BLOCK);

    tab.sqrt_x_dag_many(ql[0]);

    tab.cz_block(ql[0][0], ql[4][0], MSD_BLOCK);
    tab.cz_block(ql[1][0], ql[3][0], MSD_BLOCK);

    for block in ql.iter().take(5) {
        tab.sqrt_x_dag_many(block);
    }
    tab
}

/// Full fused MSD workload: build, then `measure_all()`, as an 85-character
/// bitstring (`msd_func_fused`).
pub fn msd_bitstring_fused<D: Driver>(seed: Option<u64>) -> String {
    let mut tab: D = msd_state_fused(seed);
    tab.measure_all()
        .into_iter()
        .map(|o| if o.unwrap() { '1' } else { '0' })
        .collect()
}

/// Integration baseline #3: the branchy two-qubit-rotation brickwork
/// (`benches/rot2-apply.rs::rot2_brickwork`).
pub fn rot2_brickwork<D: Driver>(n: usize, layers: usize) -> D {
    use std::f64::consts::PI;
    let mut tab: D = D::new_seeded(n, 1e-10, 1);
    for q in (0..n).step_by(2) {
        tab.h(q);
    }
    for layer in 0..layers {
        for a in (0..n.saturating_sub(1)).step_by(2) {
            tab.rxx(a, a + 1, 0.3 * PI);
            tab.ryy(a, a + 1, 0.4 * PI);
        }
        for a in (1..n.saturating_sub(1)).step_by(2) {
            tab.rzz(a, a + 1, 0.25 * PI);
            tab.rxx(a, a + 1, 0.15 * PI);
        }
        if layer % 2 == 0 {
            for q in (1..n).step_by(2) {
                tab.h(q);
            }
        }
    }
    tab
}

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

/// Fork `tab` over `seeds` independent RNG streams, measure every qubit, and
/// fold every outcome bit into an FNV-1a digest — the old crate's
/// `tests/apply_path.rs` golden-master instrument, reproduced exactly.
pub fn measure_record_digest<D: Driver>(tab: &D, n: usize, seeds: u64) -> u64 {
    let mut h = FNV_OFFSET;
    for seed in 0..seeds {
        let mut forked = tab.fork(Some(seed));
        for q in 0..n {
            let bit = forked.measure(q).expect("no lost qubits in this circuit");
            h ^= bit as u64;
            h = h.wrapping_mul(FNV_PRIME);
        }
    }
    h
}

/// Integration baseline #4 setup: batched Cliffords + `H` on the T targets.
pub fn fused_tgate_setup<D: Driver>(n_tgates: usize) -> D {
    let n_qubits = MSD_QUBITS;
    let mut tab: D = D::new_entropy(n_qubits, MSD_THRESHOLD);
    let block1: Vec<usize> = (0..17).collect();
    let block2: Vec<usize> = (17..34).collect();
    tab.sqrt_y_many(&block1);
    tab.sqrt_x_many(&block2);
    tab.cz_block_pairs(0, 17, 17);
    for i in 0..n_tgates {
        tab.h(i);
    }
    tab
}

/// Integration baseline #4 timed body: the T layer, two batched Clifford
/// layers, then all 85 measurements.
pub fn fused_tgate_body<D: Driver>(tab: &mut D, n_tgates: usize) -> Vec<Option<bool>> {
    let block1: Vec<usize> = (0..17).collect();
    let block2: Vec<usize> = (17..34).collect();
    for i in 0..n_tgates {
        tab.t(i);
    }
    tab.sqrt_x_dag_many(&block1);
    tab.sqrt_y_dag_many(&block2);
    (0..MSD_QUBITS).map(|i| tab.measure(i)).collect()
}

/// Integration baseline #5: the CNOT-chain + T-gate circuit, then `measure(0)`.
pub fn scaling_circuit<D: Driver>(tab: &mut D) -> Option<bool> {
    let n = tab.n_qubits();
    tab.h(0);
    tab.t(0);
    for i in 0..n - 1 {
        tab.cnot(i, i + 1);
    }
    tab.t(n - 1);
    tab.t(n - 2);
    tab.measure(0)
}

/// Integration baseline #5 without the trailing measurement (so the frame can be
/// snapshotted, or the whole `n`-qubit measurement sweep timed separately).
pub fn scaling_prepare<D: Driver>(tab: &mut D) {
    let n = tab.n_qubits();
    tab.h(0);
    tab.t(0);
    for i in 0..n - 1 {
        tab.cnot(i, i + 1);
    }
    tab.t(n - 1);
    tab.t(n - 2);
}

/// Integration baseline #7: one noisy-Clifford shot on 2 qubits, returning
/// `⟨ZZ⟩`.
pub fn noisy_shot<D: Driver>(shot: u64) -> f64 {
    let mut tab: D = D::new_seeded(2, 1e-12, shot);
    tab.h(0);
    tab.depolarize1(0, 0.05);
    tab.cnot(0, 1);
    tab.depolarize1(0, 0.05);
    tab.depolarize1(1, 0.05);
    tab.expectation_str("ZZ")
}

/// Integration baseline #8: grow the amplitude vector to exactly `2^j` branches
/// with `j` `H`+`T` pairs on an 80-qubit, threshold-0 state.
pub fn branch_grow<D: Driver>(j: usize) -> D {
    let mut tab: D = D::new_seeded(80, 0.0, 12345);
    for i in 0..j {
        tab.h(i);
        tab.t(i);
    }
    tab
}
