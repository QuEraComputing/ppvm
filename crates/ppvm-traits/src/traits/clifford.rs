// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::traits::PauliWordTrait;

/// The minimal Clifford gate set: the single-qubit Paulis (`X`, `Y`, `Z`),
/// Hadamard (`H`), phase gate (`S`), and the two entangling Cliffords
/// `CNOT` and `CZ`.
///
/// Implemented by `PauliSum`, by every tableau type, and — via the
/// blanket impl in this module — by every [`PauliWordTrait`]
/// implementor.
///
/// # Examples
///
/// Build the GHZ-preparation circuit on a `PauliSum` (Heisenberg picture,
/// so gates are applied in reverse). `PauliSum` lives in the downstream
/// `ppvm-pauli-sum` crate, so this example is `ignore`d here:
///
/// ```ignore
/// use ppvm_pauli_sum::prelude::*;
///
/// let mut state: PauliSum<config::indexmap::ByteFxHashF64<1>> =
///     PauliSum::builder().n_qubits(2).build();
/// state += ("ZZ", 1.0);
/// state.cnot([0, 1]);
/// state.h(0);
/// assert_eq!(state.len(), 1);
/// ```
pub trait Clifford {
    /// Apply Pauli `X` to each target.
    fn x(&mut self, targets: impl crate::traits::Targets);
    /// Apply Pauli `Y` to each target.
    fn y(&mut self, targets: impl crate::traits::Targets);
    /// Apply Pauli `Z` to each target.
    fn z(&mut self, targets: impl crate::traits::Targets);
    /// Apply Hadamard `H` to each target.
    fn h(&mut self, targets: impl crate::traits::Targets);
    /// Apply phase gate `S` to each target.
    fn s(&mut self, targets: impl crate::traits::Targets);
    /// Apply `CNOT` to each consecutive `(control, target)` pair.
    fn cnot(&mut self, targets: impl crate::traits::Targets);
    /// Apply `CZ` to each consecutive pair.
    fn cz(&mut self, targets: impl crate::traits::Targets);

    /// stim alias for [`cnot`](Clifford::cnot).
    fn cx(&mut self, targets: impl crate::traits::Targets) {
        self.cnot(targets)
    }
    /// stim alias for [`cnot`](Clifford::cnot).
    fn zcx(&mut self, targets: impl crate::traits::Targets) {
        self.cnot(targets)
    }
    /// stim alias for [`cz`](Clifford::cz).
    fn zcz(&mut self, targets: impl crate::traits::Targets) {
        self.cz(targets)
    }
}

/// Additional Clifford gates beyond the minimal set: `S†`, `√X`, `√X†`,
/// `√Y`, `√Y†`, and `CY`.
pub trait CliffordExtensions: Clifford {
    /// Apply `S†` to each target.
    fn s_dag(&mut self, targets: impl crate::traits::Targets);
    /// Apply `√X` to each target.
    fn sqrt_x(&mut self, targets: impl crate::traits::Targets);
    /// Apply `(√X)†` to each target.
    fn sqrt_x_dag(&mut self, targets: impl crate::traits::Targets);
    /// Apply `√Y` to each target.
    fn sqrt_y(&mut self, targets: impl crate::traits::Targets);
    /// Apply `(√Y)†` to each target.
    fn sqrt_y_dag(&mut self, targets: impl crate::traits::Targets);
    /// Apply `CY` to each consecutive pair.
    fn cy(&mut self, targets: impl crate::traits::Targets);
    /// stim alias for [`cy`](CliffordExtensions::cy).
    fn zcy(&mut self, targets: impl crate::traits::Targets) {
        self.cy(targets)
    }
}

// === Blanket Clifford impl for PauliWordTrait ===
//
// Any type implementing [`PauliWordTrait`] automatically gets word-level
// Clifford gate behavior. X/Y/Z are no-ops at the word level (they only
// affect phase, tracked separately in `PhasedPauliWord`). H, S, CNOT, CZ
// transform the bit representation. All gates honor loss bits via
// `get_lbit`, which returns `false` for non-lossy types.
//
// New PauliWordTrait implementors get this behavior for free.
//
// Breaking change for downstream crates: this blanket impl conflicts with
// any external `impl Clifford for MyWord where MyWord: PauliWordTrait`.
// Downstreams that need custom Clifford semantics must not implement
// `PauliWordTrait` on that type — see `PauliWordTrait`'s docs and the
// `PhasedPauliWord` impl in `crate::phase::clifford` for the pattern.

impl<T: PauliWordTrait> Clifford for T {
    #[inline]
    fn x(&mut self, _targets: impl crate::traits::Targets) {
        // X * I * X = I    00 -> 00, 0
        // X * X * X = X    10 -> 10, 0
        // X * Z * X = -Z   01 -> 01, 1
        // X * Y * X = -Y   11 -> 11, 1
        // word-level no-op: phase tracked at PhasedPauliWord level
    }

    #[inline]
    fn y(&mut self, _targets: impl crate::traits::Targets) {
        // word-level no-op
    }

    #[inline]
    fn z(&mut self, _targets: impl crate::traits::Targets) {
        // word-level no-op
    }

    #[inline]
    fn h(&mut self, targets: impl crate::traits::Targets) {
        // H * I * H = I, H * X * H = Z, H * Z * H = X, H * Y * H = -Y
        for index in targets.each() {
            if self.get_lbit(index) {
                continue;
            }
            let ix = self.get_xbit(index);
            let iz = self.get_zbit(index);
            self.set_xbit(index, iz);
            self.set_zbit(index, ix);
            self.rehash();
        }
    }

    #[inline]
    fn s(&mut self, targets: impl crate::traits::Targets) {
        // S * I * S = I, S * X * S = Y, S * Z * S = Z, S * Y * S = -X
        for index in targets.each() {
            if self.get_lbit(index) {
                continue;
            }
            let z = self.get_xbit(index) ^ self.get_zbit(index);
            self.set_zbit(index, z);
            self.rehash();
        }
    }

    #[inline]
    fn cnot(&mut self, targets: impl crate::traits::Targets) {
        for (control, target) in targets.pairs() {
            if self.get_lbit(control) || self.get_lbit(target) {
                continue;
            }
            let control_z = self.get_zbit(target) ^ self.get_zbit(control);
            let target_x = self.get_xbit(control) ^ self.get_xbit(target);
            self.set_zbit(control, control_z);
            self.set_xbit(target, target_x);
            self.rehash();
        }
    }

    #[inline]
    fn cz(&mut self, targets: impl crate::traits::Targets) {
        for (control, target) in targets.pairs() {
            if self.get_lbit(control) || self.get_lbit(target) {
                continue;
            }
            // flip the control z if target x is 1
            let control_z = self.get_zbit(control) ^ self.get_xbit(target);
            self.set_zbit(control, control_z);
            // flip the target z if control x is 1
            let target_z = self.get_zbit(target) ^ self.get_xbit(control);
            self.set_zbit(target, target_z);
            self.rehash();
        }
    }
}

/// Batched Clifford gates: apply the same gate to many qubits in one call.
///
/// Default implementations loop over the corresponding single-qubit
/// (or single-pair) method on [`Clifford`]. Types like the stabilizer
/// `Tableau` override these methods with a fused inner-loop or bitmask
/// implementation. Types that don't need specialization can implement
/// this trait with an empty `impl` to use the defaults.
pub trait CliffordBatch: Clifford {
    /// Apply Pauli `X` to every qubit in `indices`.
    fn x_batch(&mut self, indices: &[usize]) {
        for &q in indices {
            self.x(q);
        }
    }
    /// Apply Pauli `Y` to every qubit in `indices`.
    fn y_batch(&mut self, indices: &[usize]) {
        for &q in indices {
            self.y(q);
        }
    }
    /// Apply Pauli `Z` to every qubit in `indices`.
    fn z_batch(&mut self, indices: &[usize]) {
        for &q in indices {
            self.z(q);
        }
    }
    /// Apply Hadamard `H` to every qubit in `indices`.
    fn h_batch(&mut self, indices: &[usize]) {
        for &q in indices {
            self.h(q);
        }
    }
    /// Apply phase gate `S` to every qubit in `indices`.
    fn s_batch(&mut self, indices: &[usize]) {
        for &q in indices {
            self.s(q);
        }
    }
    /// Apply `CNOT` to every `(control, target)` pair.
    fn cnot_batch(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cnot([c, t]);
        }
    }
    /// Apply `CZ` to every `(control, target)` pair.
    fn cz_batch(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cz([c, t]);
        }
    }
}

/// Batched form of [`CliffordExtensions`].
pub trait CliffordExtensionsBatch: CliffordExtensions + CliffordBatch {
    /// Apply `S†` to every qubit in `indices`.
    fn s_adj_batch(&mut self, indices: &[usize]) {
        for &q in indices {
            self.s_dag(q);
        }
    }
    /// Apply `√X` to every qubit in `indices`.
    fn sqrt_x_batch(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_x(q);
        }
    }
    /// Apply `(√X)†` to every qubit in `indices`.
    fn sqrt_x_adj_batch(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_x_dag(q);
        }
    }
    /// Apply `√Y` to every qubit in `indices`.
    fn sqrt_y_batch(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_y(q);
        }
    }
    /// Apply `(√Y)†` to every qubit in `indices`.
    fn sqrt_y_adj_batch(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_y_dag(q);
        }
    }
    /// Apply `CY` to every `(control, target)` pair.
    fn cy_batch(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cy([c, t]);
        }
    }
}

impl<T: PauliWordTrait> CliffordExtensions for T {
    // |    Gate    |  X  |  Y  |  Z  |
    // |:----------:|:---:|:---:|:---:|
    // |     s      | -Y  |  X  |  Z  |
    // |   s_adj    |  Y  | -X  |  Z  |
    // |   sqrt_x   |  X  | -Z  |  Y  |
    // | sqrt*x*adj |  X  |  Z  | -Y  |
    // |   sqrt_y   |  Z  |  Y  | -X  |
    // | sqrt*y*adj | -Z  |  Y  |  X  |

    #[inline]
    fn s_dag(&mut self, targets: impl crate::traits::Targets) {
        // s_dag has the same bit mapping as s (only phases differ)
        for addr0 in targets.each() {
            self.s(addr0);
        }
    }

    #[inline]
    fn sqrt_x(&mut self, targets: impl crate::traits::Targets) {
        for addr0 in targets.each() {
            if self.get_lbit(addr0) {
                continue;
            }
            let x = self.get_xbit(addr0);
            let z = self.get_zbit(addr0);
            self.set_xbit(addr0, x ^ z);
            self.rehash();
        }
    }

    #[inline]
    fn sqrt_x_dag(&mut self, targets: impl crate::traits::Targets) {
        for addr0 in targets.each() {
            if self.get_lbit(addr0) {
                continue;
            }
            let x = self.get_xbit(addr0);
            let z = self.get_zbit(addr0);
            self.set_xbit(addr0, x ^ z);
            self.rehash();
        }
    }

    #[inline]
    fn sqrt_y(&mut self, targets: impl crate::traits::Targets) {
        for addr0 in targets.each() {
            if self.get_lbit(addr0) {
                continue;
            }
            let x = self.get_xbit(addr0);
            let z = self.get_zbit(addr0);
            self.set_xbit(addr0, z);
            self.set_zbit(addr0, x);
            self.rehash();
        }
    }

    #[inline]
    fn sqrt_y_dag(&mut self, targets: impl crate::traits::Targets) {
        for addr0 in targets.each() {
            if self.get_lbit(addr0) {
                continue;
            }
            let x = self.get_xbit(addr0);
            let z = self.get_zbit(addr0);
            self.set_xbit(addr0, z);
            self.set_zbit(addr0, x);
            self.rehash();
        }
    }

    // | CY  |  I  |  X  |  Y  |  Z  |
    // |:---:|:---:|:---:|:---:|:---:|
    // |  I  | II  | ZX  | IY  | ZZ  |
    // |  X  | XY  | -YZ | XI  | YX  |
    // |  Y  | YY  | XZ  | YI  | -XX |
    // |  Z  | ZI  | IX  | ZY  | IZ  |
    #[inline]
    fn cy(&mut self, targets: impl crate::traits::Targets) {
        for (addr0, addr1) in targets.pairs() {
            if self.get_lbit(addr0) || self.get_lbit(addr1) {
                continue;
            }
            let xc = self.get_xbit(addr0);
            let zc = self.get_zbit(addr0);
            let xt = self.get_xbit(addr1);
            let zt = self.get_zbit(addr1);
            self.set_zbit(addr0, zc ^ xt ^ zt);
            self.set_xbit(addr1, xt ^ xc);
            self.set_zbit(addr1, zt ^ xc);
            self.rehash();
        }
    }
}
