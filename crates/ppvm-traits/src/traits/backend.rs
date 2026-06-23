// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Device-portable Pauli backend abstraction.
//!
//! This is the **de-biased** replacement for the closure-driven [`ACMap`] trait
//! family. See `docs/design/device-portable-pauli-backend.md` for the full
//! rationale. The short version:
//!
//! * [`ACMap`] expresses every transformation as a host Rust closure
//!   (`Fn(&W, &V) -> …`) and iterates by host reference (`(&W, &V)`). Neither can
//!   cross to a CUDA kernel / device-resident map.
//! * [`PauliBackend`] instead exposes a **fixed, reified vocabulary** of Pauli
//!   propagation operations (gate identity + qubit indices + scalar coeffs). Each
//!   backend implements them with whatever it likes — a CPU loop or a GPU kernel —
//!   so the *same* `PauliSum<T>` can be driven by a Rust `HashMap` or a cuco
//!   device map, giving `CuPauliSum = PauliSum<CudaConfig>`.
//!
//! [`ACMap`]: crate::traits::ACMap
//!
//! ## Status
//!
//! Scaffold for `refactor/device-portable-pauli-map`. The trait shape and the
//! Clifford-1q slice of the [`reference`] CPU backend are implemented and tested;
//! the branching ops (rotations) and 2q Cliffords are marked `ITERATE` for
//! follow-up. Nothing in the existing crate is wired to this yet (additive).

use crate::char::Pauli;
use crate::traits::Coefficient;

/// Single-qubit Clifford gates, reified so a backend can dispatch on them
/// (CPU branch / GPU kernel selection) instead of receiving a host closure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Clifford1q {
    /// Pauli `X` (word-level no-op; conjugation flips the sign of Z/Y terms).
    X,
    /// Pauli `Y` (flips the sign of X/Z terms).
    Y,
    /// Pauli `Z` (flips the sign of X/Y terms).
    Z,
    /// Hadamard.
    H,
    /// Phase gate `S`.
    S,
    /// `S†`.
    Sdag,
    /// `√X`.
    SqrtX,
    /// `(√X)†`.
    SqrtXdag,
    /// `√Y`.
    SqrtY,
    /// `(√Y)†`.
    SqrtYdag,
}

/// Two-qubit Clifford gates, reified.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Clifford2q {
    /// Controlled-NOT.
    Cnot,
    /// Controlled-Z.
    Cz,
    /// Controlled-Y.
    Cy,
}

/// A Pauli word packed into two bitmasks (one bit per qubit): `x` carries the
/// X-component, `z` the Z-component, so `(x,z) = (0,0)=I (1,0)=X (0,1)=Z (1,1)=Y`.
///
/// This is the canonical host-interop / export form — deliberately the same
/// encoding the CUDA backend uses on the device. `(u64, u64)` covers
/// `n_qubits <= 64`; wider words are a follow-up (see the design doc).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct PackedWord {
    /// X-component bitmask.
    pub x: u64,
    /// Z-component bitmask.
    pub z: u64,
}

impl PackedWord {
    /// `true` if qubit `q` carries an X-component.
    #[inline]
    pub fn xbit(&self, q: usize) -> bool {
        (self.x >> q) & 1 == 1
    }
    /// `true` if qubit `q` carries a Z-component.
    #[inline]
    pub fn zbit(&self, q: usize) -> bool {
        (self.z >> q) & 1 == 1
    }
    /// Set/clear qubit `q`'s X-component.
    #[inline]
    pub fn set_xbit(&mut self, q: usize, v: bool) {
        self.x = (self.x & !(1u64 << q)) | ((v as u64) << q);
    }
    /// Set/clear qubit `q`'s Z-component.
    #[inline]
    pub fn set_zbit(&mut self, q: usize, v: bool) {
        self.z = (self.z & !(1u64 << q)) | ((v as u64) << q);
    }
}

/// The device-portable backend a [`PauliSum`](../../../ppvm_pauli_sum/sum/struct.PauliSum.html)
/// can be built on. Implemented by CPU maps and by the cuco device map.
///
/// Every method is closure-free and word-type-free in the hot path: operations
/// are described by data (gate + indices + scalars), and the backend owns its key
/// representation. Host interop is confined to [`export`](Self::export) and the
/// constructors.
pub trait PauliBackend: Sized {
    /// Coefficient type stored per Pauli term.
    type Coeff: Coefficient;

    // ---- container primitives (closure-free) -------------------------------

    /// Empty backend over `n_qubits` qubits with reserved `capacity`.
    fn with_capacity(n_qubits: usize, capacity: usize) -> Self;
    /// Number of qubits the words span.
    fn n_qubits(&self) -> usize;
    /// Number of stored (live) terms.
    fn len(&self) -> usize;
    /// `true` if there are no terms.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Remove every term, retaining capacity.
    fn clear(&mut self);
    /// Multiply every coefficient by `factor`.
    fn scale_all(&mut self, factor: Self::Coeff);
    /// Drain `other` into `self`, summing coefficients on key collision.
    fn merge(&mut self, other: &mut Self);

    // ---- reified gate operations -------------------------------------------

    /// Apply a single-qubit Clifford (Heisenberg conjugation `P ↦ U† P U`).
    fn apply_clifford_1q(&mut self, gate: Clifford1q, q: usize);
    /// Apply a two-qubit Clifford to the `(a, b)` pair.
    fn apply_clifford_2q(&mut self, gate: Clifford2q, a: usize, b: usize);
    /// Apply a single-qubit Pauli rotation about `axis` by an angle whose
    /// `sin`/`cos` are precomputed by the caller. May branch each term into two.
    fn apply_rotation_1q(&mut self, axis: Pauli, q: usize, sin: Self::Coeff, cos: Self::Coeff);
    /// Apply a two-qubit Pauli rotation; axes are `[x_bit, z_bit]` Pauli codes.
    fn apply_rotation_2q(
        &mut self,
        axis_a: [u8; 2],
        axis_b: [u8; 2],
        a: usize,
        b: usize,
        sin: Self::Coeff,
        cos: Self::Coeff,
    );

    // ---- reductions / truncation / interop ---------------------------------

    /// Drop every term whose coefficient magnitude is `< eps` (SPD truncation,
    /// applied after merge).
    fn truncate_abs(&mut self, eps: f64);
    /// Coefficient overlap `Σ_k self[k] · other[k]` over matching keys.
    fn overlap(&self, other: &Self) -> Self::Coeff;
    /// Copy live terms to the host in the canonical packed encoding.
    fn export(&self) -> Vec<(PackedWord, Self::Coeff)>;
}

pub mod reference {
    //! Illustrative CPU reference backend over packed `u64` keys.
    //!
    //! This exists to **validate the [`PauliBackend`] shape** and to give the
    //! migration a correctness oracle — it is intentionally over `PackedWord`
    //! keys (mirroring the CUDA device encoding), not the production
    //! `HashMap<PauliWord, C>` path, and it pins `Coeff = f64`. The Clifford-1q
    //! family is fully implemented and tested; branching ops are `ITERATE`.

    use super::*;
    use std::collections::HashMap;

    /// Reference `PauliBackend` over `HashMap<PackedWord, f64>`.
    #[derive(Clone, Debug, Default)]
    pub struct HashMapBackend {
        map: HashMap<PackedWord, f64>,
        n_qubits: usize,
    }

    impl HashMapBackend {
        /// Build directly from `(word, coeff)` terms (test/interop convenience).
        pub fn from_terms(n_qubits: usize, terms: &[(PackedWord, f64)]) -> Self {
            let mut b = Self::with_capacity(n_qubits, terms.len());
            for &(k, v) in terms {
                *b.map.entry(k).or_insert(0.0) += v;
            }
            b
        }
    }

    /// For X/Y/Z: does conjugation flip the sign of the term carrying Pauli at
    /// `q`? (Mirrors `ppvm-pauli-sum/src/sum/clifford.rs`.)
    #[inline]
    fn sign_flip_pauli(gate: Clifford1q, k: &PackedWord, q: usize) -> bool {
        let (x, z) = (k.xbit(q), k.zbit(q));
        match gate {
            Clifford1q::X => z,     // flips Z, Y
            Clifford1q::Y => x ^ z, // flips X, Z
            Clifford1q::Z => x,     // flips X, Y
            _ => false,
        }
    }

    /// For the key-changing single-qubit Cliffords: rewrite `k` in place and
    /// return whether the coefficient sign flips. (Mirrors clifford.rs exactly.)
    #[inline]
    fn transform_clifford(gate: Clifford1q, k: &mut PackedWord, q: usize) -> bool {
        let (x, z) = (k.xbit(q), k.zbit(q));
        match gate {
            Clifford1q::H => {
                k.set_xbit(q, z);
                k.set_zbit(q, x);
                x & z // flip for Y
            }
            Clifford1q::S => {
                k.set_zbit(q, x ^ z);
                x & !z // flip for X
            }
            Clifford1q::Sdag => {
                k.set_zbit(q, x ^ z);
                x & z // flip for Y
            }
            Clifford1q::SqrtX => {
                k.set_xbit(q, x ^ z);
                x & z // flip for Y
            }
            Clifford1q::SqrtXdag => {
                k.set_xbit(q, x ^ z);
                !x & z // flip for Z
            }
            Clifford1q::SqrtY => {
                k.set_xbit(q, z);
                k.set_zbit(q, x);
                !x & z // flip for Z
            }
            Clifford1q::SqrtYdag => {
                k.set_xbit(q, z);
                k.set_zbit(q, x);
                x & !z // flip for X
            }
            // X/Y/Z don't change the key; handled by `sign_flip_pauli`.
            Clifford1q::X | Clifford1q::Y | Clifford1q::Z => false,
        }
    }

    impl PauliBackend for HashMapBackend {
        type Coeff = f64;

        fn with_capacity(n_qubits: usize, capacity: usize) -> Self {
            Self {
                map: HashMap::with_capacity(capacity),
                n_qubits,
            }
        }
        fn n_qubits(&self) -> usize {
            self.n_qubits
        }
        fn len(&self) -> usize {
            self.map.len()
        }
        fn clear(&mut self) {
            self.map.clear();
        }
        fn scale_all(&mut self, factor: f64) {
            for v in self.map.values_mut() {
                *v *= factor;
            }
        }
        fn merge(&mut self, other: &mut Self) {
            for (k, v) in other.map.drain() {
                *self.map.entry(k).or_insert(0.0) += v;
            }
        }

        fn apply_clifford_1q(&mut self, gate: Clifford1q, q: usize) {
            match gate {
                // Sign-only gates: keys are unchanged, mutate values in place.
                Clifford1q::X | Clifford1q::Y | Clifford1q::Z => {
                    for (k, v) in self.map.iter_mut() {
                        if sign_flip_pauli(gate, k, q) {
                            *v = -*v;
                        }
                    }
                }
                // Key-changing gates: rebuild, accumulating on collision.
                _ => {
                    let old = std::mem::take(&mut self.map);
                    self.map.reserve(old.len());
                    for (mut k, v) in old {
                        let v = if transform_clifford(gate, &mut k, q) {
                            -v
                        } else {
                            v
                        };
                        *self.map.entry(k).or_insert(0.0) += v;
                    }
                }
            }
        }

        fn apply_clifford_2q(&mut self, _gate: Clifford2q, _a: usize, _b: usize) {
            unimplemented!("ITERATE: 2q Clifford (cnot/cz/cy) on the packed-key reference backend");
        }
        fn apply_rotation_1q(&mut self, _axis: Pauli, _q: usize, _sin: f64, _cos: f64) {
            unimplemented!(
                "ITERATE: 1q rotation — branching (≤2 terms) + merge; reuse levi_civita"
            );
        }
        fn apply_rotation_2q(
            &mut self,
            _axis_a: [u8; 2],
            _axis_b: [u8; 2],
            _a: usize,
            _b: usize,
            _sin: f64,
            _cos: f64,
        ) {
            unimplemented!("ITERATE: 2q rotation — branching + merge; reuse comm_2");
        }

        fn truncate_abs(&mut self, eps: f64) {
            self.map.retain(|_, v| v.abs() >= eps);
        }
        fn overlap(&self, other: &Self) -> f64 {
            // Iterate the smaller map; look up the other.
            let (small, big) = if self.map.len() <= other.map.len() {
                (&self.map, &other.map)
            } else {
                (&other.map, &self.map)
            };
            small
                .iter()
                .filter_map(|(k, v)| big.get(k).map(|w| v * w))
                .sum()
        }
        fn export(&self) -> Vec<(PackedWord, f64)> {
            self.map.iter().map(|(&k, &v)| (k, v)).collect()
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        // Pauli encoding helpers for 1-qubit words.
        const I: PackedWord = PackedWord { x: 0, z: 0 };
        const X: PackedWord = PackedWord { x: 1, z: 0 };
        const Z: PackedWord = PackedWord { x: 0, z: 1 };
        const Y: PackedWord = PackedWord { x: 1, z: 1 };

        fn one(word: PackedWord, c: f64) -> HashMapBackend {
            HashMapBackend::from_terms(1, &[(word, c)])
        }

        fn only_term(b: &HashMapBackend) -> (PackedWord, f64) {
            assert_eq!(b.len(), 1, "expected exactly one term");
            let mut t = b.export();
            t.pop().unwrap()
        }

        #[test]
        fn pauli_sign_flips_match_conjugation() {
            // X conjugation flips Z and Y, leaves I and X.
            for (w, flips) in [(I, false), (X, false), (Z, true), (Y, true)] {
                let mut b = one(w, 3.0);
                b.apply_clifford_1q(Clifford1q::X, 0);
                assert_eq!(only_term(&b), (w, if flips { -3.0 } else { 3.0 }));
            }
        }

        #[test]
        fn hadamard_swaps_x_and_z() {
            // H: X→Z, Z→X (no sign), Y→-Y, I→I.
            let mut b = one(X, 1.0);
            b.apply_clifford_1q(Clifford1q::H, 0);
            assert_eq!(only_term(&b), (Z, 1.0));

            let mut b = one(Y, 1.0);
            b.apply_clifford_1q(Clifford1q::H, 0);
            assert_eq!(only_term(&b), (Y, -1.0)); // sign flip for Y
        }

        #[test]
        fn s_gate_maps_x_to_y() {
            // S: X→-Y? clifford.rs: S flips sign for X (xbit set, zbit clear),
            // and sets z = x^z so X(1,0)→(1,1)=Y. So X → -Y.
            let mut b = one(X, 1.0);
            b.apply_clifford_1q(Clifford1q::S, 0);
            assert_eq!(only_term(&b), (Y, -1.0));
            // Z is unchanged by S (z=x^z with x=0).
            let mut b = one(Z, 2.0);
            b.apply_clifford_1q(Clifford1q::S, 0);
            assert_eq!(only_term(&b), (Z, 2.0));
        }

        #[test]
        fn key_changing_gate_merges_collisions() {
            // X term and Z term; after H they become Z and X respectively — no
            // collision. Use H on {X:1, Z:1} → {Z:1, X:1}: still two terms.
            let mut b = HashMapBackend::from_terms(1, &[(X, 1.0), (Z, 1.0)]);
            b.apply_clifford_1q(Clifford1q::H, 0);
            assert_eq!(b.len(), 2);
            // Now construct a real collision: {Z:1, X:1} then H maps both to the
            // other's slot; values are distinct so still 2. A genuine merge:
            // start {X:2, Z:3}; apply H twice (identity) → back to {X:2, Z:3}.
            b.apply_clifford_1q(Clifford1q::H, 0);
            let mut got = b.export();
            got.sort_by_key(|(k, _)| (k.x, k.z));
            // Sorted by (x,z): Z=(0,1) precedes X=(1,0).
            assert_eq!(got, vec![(Z, 1.0), (X, 1.0)]);
        }

        #[test]
        fn scale_truncate_overlap() {
            let mut b = HashMapBackend::from_terms(1, &[(X, 2.0), (Z, 0.001)]);
            b.scale_all(2.0);
            b.truncate_abs(0.01); // drops the 0.002 Z term
            assert_eq!(b.len(), 1);
            assert_eq!(only_term(&b), (X, 4.0));

            let a = HashMapBackend::from_terms(1, &[(X, 2.0), (Z, 3.0)]);
            let c = HashMapBackend::from_terms(1, &[(X, 5.0), (Y, 7.0)]);
            assert_eq!(a.overlap(&c), 2.0 * 5.0); // only X matches
        }

        #[test]
        fn merge_accumulates() {
            let mut a = HashMapBackend::from_terms(1, &[(X, 1.0), (Z, 2.0)]);
            let mut b = HashMapBackend::from_terms(1, &[(Z, 10.0), (Y, 100.0)]);
            a.merge(&mut b);
            assert_eq!(a.len(), 3);
            assert!(b.is_empty());
            let mut got = a.export();
            got.sort_by_key(|(k, _)| (k.x, k.z));
            // Sorted by (x,z): Z=(0,1), X=(1,0), Y=(1,1).
            assert_eq!(got, vec![(Z, 12.0), (X, 1.0), (Y, 100.0)]);
        }
    }
}
