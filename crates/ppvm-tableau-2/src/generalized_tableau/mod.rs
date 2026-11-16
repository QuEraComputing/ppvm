// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`GeneralizedTableau`]: the simulator, `|ψ⟩ = C |χ⟩`.
//!
//! - `gates`: Clifford, T, rotation and `U3` gates, and reset.
//! - `measure`: Z measurement.
//! - `loss`: loss and leakage channels.
//! - `noise`: Pauli and depolarizing channels.

mod gates;
mod loss;
mod measure;
mod noise;
mod status;

use rand::SeedableRng;
use rand::rngs::SmallRng;

use num_complex::Complex64;
use ppvm_traits_2::{Clifford, Pauli};

use crate::{InverseTableau, VirtualState};
use status::Status;

/// The state `|ψ⟩ = C |χ⟩`: the Clifford `C` as an [`InverseTableau`], and the
/// virtual state `|χ⟩` in the representation `V`.
///
/// Physical gates `C ← G·C` touch only the tableau. Measurements ask `V` for
/// the outcome weights, sample here, and let `V` collapse. Lost and leaked
/// qubits are tracked here; the tableau and `V` never see them.
#[derive(Clone, Debug)]
pub struct GeneralizedTableau<V: VirtualState> {
    tableau: InverseTableau,
    state: V,
    status: Status,
    rng: SmallRng,
}

impl<V: VirtualState> GeneralizedTableau<V> {
    /// `|0…0⟩` on `n` qubits, with the RNG seeded from OS entropy.
    pub fn new(n: usize) -> Self {
        Self::with_rng(n, rand::make_rng())
    }

    /// Same as [`new`](Self::new), but seed the RNG deterministically.
    pub fn new_with_seed(n: usize, seed: u64) -> Self {
        Self::with_rng(n, SmallRng::seed_from_u64(seed))
    }

    fn with_rng(n: usize, rng: SmallRng) -> Self {
        Self {
            tableau: InverseTableau::new(n),
            state: V::new(n),
            status: Status::new(n),
            rng,
        }
    }

    pub fn num_qubits(&self) -> usize {
        self.tableau.num_qubits()
    }

    pub fn tableau(&self) -> &InverseTableau {
        &self.tableau
    }

    pub fn state(&self) -> &V {
        &self.state
    }

    /// Whether `q` is lost: gates skip it and measuring it gives `None`.
    pub fn is_lost(&self, q: usize) -> bool {
        self.status.is_lost(q)
    }

    /// Whether `q` is leaked: gates skip it and measuring it gives the pinned bit.
    pub fn is_leaked(&self, q: usize) -> bool {
        self.status.is_leaked(q)
    }

    /// `|ψ⟩ ← (a·I + b·∏ Z_q) |ψ⟩` over `qubits`, which must be distinct.
    /// Inactive qubits are dropped, as gates skip them.
    fn rotate_z(&mut self, qubits: &[usize], a: Complex64, b: Complex64) {
        let live: Vec<usize>;
        let qubits = if self.status.any() {
            live = qubits
                .iter()
                .copied()
                .filter(|&q| !self.status.is_inactive(q))
                .collect();
            if live.is_empty() {
                return;
            }
            &live[..]
        } else {
            qubits
        };
        let mut view = self.tableau.virtual_view();
        self.state.rotate(&mut view, qubits, a, b);
    }

    /// `exp(−i θ/2 ∏ Z_q) = cos(θ/2)·I − i sin(θ/2)·∏ Z_q`.
    fn rz_product(&mut self, qubits: &[usize], theta: f64) {
        let (s, c) = (theta / 2.0).sin_cos();
        self.rotate_z(qubits, Complex64::new(c, 0.0), Complex64::new(0.0, -s));
    }

    /// Map `axis` on `q` to `Z`: apply `W` with `W·axis·W† = Z`.
    fn axis_to_z(&mut self, axis: Pauli, q: usize) {
        match axis {
            Pauli::X => self.h(q),
            Pauli::Y => {
                self.s_dag(q);
                self.h(q);
            }
            Pauli::Z | Pauli::I => {}
        }
    }

    /// Undo [`axis_to_z`](Self::axis_to_z): apply `W†`.
    fn z_to_axis(&mut self, axis: Pauli, q: usize) {
        match axis {
            Pauli::X => self.h(q),
            Pauli::Y => {
                self.h(q);
                self.s(q);
            }
            Pauli::Z | Pauli::I => {}
        }
    }
}
