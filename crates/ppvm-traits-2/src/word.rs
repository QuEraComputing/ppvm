// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// The common **read-only** concept for an indexed algebraic monomial.
pub trait Word {
    /// The operator alphabet at one index.
    type Site;

    /// Number of sites (for a dense Pauli word, the qubit width).
    fn n_sites(&self) -> usize;

    /// Read the site at `index`.
    fn get(&self, index: usize) -> Self::Site;

    /// Number of non-identity factors in the site alphabet.
    /// Representations without explicit identities may have `weight() == n_sites()`.
    fn weight(&self) -> usize;

    /// Iterate the sites in index order.
    fn iter(&self) -> impl Iterator<Item = Self::Site>;
}

/// Mutable single-vector X/Z access — a point of `GF(2)^{2n}`. Composite
/// defaults compose the four scalar accessors; a packed backend may override
/// any of them to refresh metadata once per call instead of once per bit.
pub trait PauliBits: Word {
    /// Read the X bit at index `i`.
    fn x_bit(&self, i: usize) -> bool;

    /// Read the Z bit at index `i`.
    fn z_bit(&self, i: usize) -> bool;

    /// Set the X bit at index `i` (refreshes any structural-hash cache).
    fn set_x_bit(&mut self, i: usize, v: bool);

    /// Set the Z bit at index `i` (refreshes any structural-hash cache).
    fn set_z_bit(&mut self, i: usize, v: bool);

    /// Set both bit planes at one site.
    #[inline(always)]
    fn set_xz_bits(&mut self, i: usize, x: bool, z: bool) {
        self.set_x_bit(i, x);
        self.set_z_bit(i, z);
    }
    /// Set both bit planes at two sites.
    #[inline(always)]
    fn set_xz_bits2(&mut self, i: usize, xi: bool, zi: bool, j: usize, xj: bool, zj: bool) {
        self.set_xz_bits(i, xi, zi);
        self.set_xz_bits(j, xj, zj);
    }
    /// Set one X bit and one Z bit together, as required by CNOT.
    #[inline(always)]
    fn set_x_bit_and_z_bit(&mut self, x_i: usize, x: bool, z_i: usize, z: bool) {
        self.set_x_bit(x_i, x);
        self.set_z_bit(z_i, z);
    }
    /// Set two Z bits together, as required by CZ, leaving X bits unchanged.
    #[inline(always)]
    fn set_z_bit_pair(&mut self, i: usize, zi: bool, j: usize, zj: bool) {
        self.set_z_bit(i, zi);
        self.set_z_bit(j, zj);
    }

    /// Packed local Pauli code: `0=I, 1=X, 2=Z, 3=Y`.
    #[inline(always)]
    fn pauli_code(&self, i: usize) -> u8 {
        (self.x_bit(i) as u8) | ((self.z_bit(i) as u8) << 1)
    }

    /// Copy this word and toggle selected X/Z bits at one site.
    #[inline]
    fn toggled_bits(&self, i: usize, toggle_x: bool, toggle_z: bool) -> Self
    where
        Self: Sized + Clone,
    {
        self.clone().into_toggled_bits(i, [toggle_x, toggle_z])
    }

    /// Copy this word once and toggle selected X/Z bits at two sites.
    /// Each mask is `[toggle_x, toggle_z]`; only one copy is made.
    #[inline]
    fn toggled_bits2(&self, i: usize, toggle_i: [bool; 2], j: usize, toggle_j: [bool; 2]) -> Self
    where
        Self: Sized + Clone,
    {
        self.clone().into_toggled_bits2(i, toggle_i, j, toggle_j)
    }

    /// Consume this word and toggle the bits selected by `[toggle_x, toggle_z]`
    /// at one site. The in-place primitive every variant above is built from.
    #[inline]
    fn into_toggled_bits(mut self, i: usize, toggle: [bool; 2]) -> Self
    where
        Self: Sized,
    {
        if toggle[0] {
            let x = self.x_bit(i);
            self.set_x_bit(i, !x);
        }
        if toggle[1] {
            let z = self.z_bit(i);
            self.set_z_bit(i, !z);
        }
        self
    }

    /// Consume this word and toggle two sites using `[toggle_x, toggle_z]` masks.
    #[inline]
    fn into_toggled_bits2(
        self,
        i: usize,
        toggle_i: [bool; 2],
        j: usize,
        toggle_j: [bool; 2],
    ) -> Self
    where
        Self: Sized,
    {
        self.into_toggled_bits(i, toggle_i)
            .into_toggled_bits(j, toggle_j)
    }

    /// Whether this word anticommutes with `pauli = (x_bit, z_bit)` at site `i`.
    /// Computes the local symplectic form `x_P·z_Q ⊕ z_P·x_Q`.
    #[inline]
    fn anticommutes_at(&self, i: usize, pauli: (bool, bool)) -> bool {
        (self.x_bit(i) & pauli.1) ^ (self.z_bit(i) & pauli.0)
    }
}
