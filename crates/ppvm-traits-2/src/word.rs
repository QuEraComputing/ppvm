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

    /// Number of non-identity factors according to the concrete site alphabet.
    ///
    /// A Pauli-motivated read (the `MaxPauliWeight` policy needs it); an ordered
    /// representation that stores no explicit identities may have
    /// `weight() == n_sites()`.
    fn weight(&self) -> usize;

    /// Iterate the sites in index order.
    fn iter(&self) -> impl Iterator<Item = Self::Site>;
}

// Mutable single-vector X/Z access — a point of `GF(2)^{2n}`.
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
    ///
    /// The default composes the scalar setters. Packed words whose setters
    /// refresh structural metadata eagerly may override this to refresh once.
    #[inline(always)]
    fn set_xz_bits(&mut self, i: usize, x: bool, z: bool) {
        self.set_x_bit(i, x);
        self.set_z_bit(i, z);
    }
    /// Set both packed bit planes at two sites.
    ///
    /// The default composes the one-site setter. Packed words with eager
    /// structural metadata can override this to refresh exactly once after the
    /// four writes.
    #[inline(always)]
    fn set_xz_bits2(&mut self, i: usize, xi: bool, zi: bool, j: usize, xj: bool, zj: bool) {
        self.set_xz_bits(i, xi, zi);
        self.set_xz_bits(j, xj, zj);
    }
    /// Set one X-plane bit and one Z-plane bit as one structural mutation.
    ///
    /// `CNOT` updates exactly this pair (`x_target`, `z_control`). Packed words
    /// with eager structural metadata can override this to refresh once without
    /// also reading or rewriting the two unchanged companion bits.
    #[inline(always)]
    fn set_x_bit_and_z_bit(&mut self, x_i: usize, x: bool, z_i: usize, z: bool) {
        self.set_x_bit(x_i, x);
        self.set_z_bit(z_i, z);
    }
    /// Set two Z-plane bits as one structural mutation.
    ///
    /// `CZ` updates exactly this pair (`z_a`, `z_b`) and leaves both X bits
    /// alone — the counterpart of [`set_x_bit_and_z_bit`](PauliBits::set_x_bit_and_z_bit)
    /// for `CNOT`. Routing it through the four-bit
    /// [`set_xz_bits2`](PauliBits::set_xz_bits2) instead makes a packed word read
    /// and rewrite the two unchanged X bits on every gate. The default composes
    /// the scalar setters.
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

    /// A copy of this word with the X and/or Z bit at `i` toggled — the
    /// **rotation-branch key builder** (`iGP` from a diagonal `P`).
    ///
    /// Provided as clone-then-flip so every `PauliBits` implementer gets a branch
    /// builder for free and the rotation/branching kernels can be generic over the
    /// word type (the ordinary and the lossy key run the *same* kernel — see
    /// `ppvm-pauli-sum-2`'s rotation and loss modules). `PauliWord` overrides it
    /// with a direct plane copy that computes the digest exactly once, skipping
    /// the redundant refresh a `clone` + `set_*_bit` pair performs.
    fn toggled_bits(&self, i: usize, toggle_x: bool, toggle_z: bool) -> Self
    where
        Self: Sized + Clone,
    {
        let mut out = self.clone();
        if toggle_x {
            let b = out.x_bit(i);
            out.set_x_bit(i, !b);
        }
        if toggle_z {
            let b = out.z_bit(i);
            out.set_z_bit(i, !b);
        }
        out
    }

    /// A copy of this word with the X and/or Z bits at **two** sites toggled —
    /// the **two-qubit** rotation-branch key builder (`iG_aG_b·P`).
    ///
    /// Chaining [`toggled_bits`](PauliBits::toggled_bits) twice would build two
    /// whole words per produced branch term, and for a packed word each of those is
    /// a full copy of *both* bit planes plus a word rebuild — so the intermediate
    /// is pure waste on the hot path of `rzz`/`rxx`/`ryy`/`rotate_2` (old built
    /// **one** `k.clone()` and then wrote four bits into it,
    /// `ppvm-pauli-sum/src/sum/rot2.rs`). The two-site entry point makes the single
    /// copy the *only* copy, and it scales with the storage tier: at `[u8; 32]` the
    /// chained form moved 64 redundant bytes per branch term.
    ///
    /// The default is clone-then-flip (one clone, up to four bit writes), which is
    /// already old's shape; `PauliWord` overrides it with a direct plane copy that
    /// computes the digest exactly once, as it does for the single-site form.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn toggled_bits2(
        &self,
        i: usize,
        toggle_x_i: bool,
        toggle_z_i: bool,
        j: usize,
        toggle_x_j: bool,
        toggle_z_j: bool,
    ) -> Self
    where
        Self: Sized + Clone,
    {
        let mut out = self.clone();
        if toggle_x_i {
            let b = out.x_bit(i);
            out.set_x_bit(i, !b);
        }
        if toggle_z_i {
            let b = out.z_bit(i);
            out.set_z_bit(i, !b);
        }
        if toggle_x_j {
            let b = out.x_bit(j);
            out.set_x_bit(j, !b);
        }
        if toggle_z_j {
            let b = out.z_bit(j);
            out.set_z_bit(j, !b);
        }
        out
    }

    /// Consume a word and toggle X/Z bits at two sites in place.
    ///
    /// Re-keying kernels already own their key, so this avoids the structural
    /// copy required by [`toggled_bits2`](PauliBits::toggled_bits2). Packed words
    /// may override it to defer derived-hash refresh until all writes complete.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn into_toggled_bits2(
        mut self,
        i: usize,
        toggle_x_i: bool,
        toggle_z_i: bool,
        j: usize,
        toggle_x_j: bool,
        toggle_z_j: bool,
    ) -> Self
    where
        Self: Sized,
    {
        if toggle_x_i {
            let b = self.x_bit(i);
            self.set_x_bit(i, !b);
        }
        if toggle_z_i {
            let b = self.z_bit(i);
            self.set_z_bit(i, !b);
        }
        if toggle_x_j {
            let b = self.x_bit(j);
            self.set_x_bit(j, !b);
        }
        if toggle_z_j {
            let b = self.z_bit(j);
            self.set_z_bit(j, !b);
        }
        self
    }

    /// Whether this word anticommutes with the single-qubit Pauli
    /// `pauli = (x_bit, z_bit)` at index `i`, i.e. whether the symplectic form
    /// `ω(P, Q) = x_P·z_Q ⊕ z_P·x_Q` is `1` there.
    ///
    /// The old `PauliWordTrait::anticommutes_at`
    /// (`ppvm-traits/src/traits/word_trait.rs`), reproduced verbatim as a
    /// provided method — it is the pivot test the tableau measurement search
    /// runs (`ppvm-tableau/src/data.rs`), and it is derivable from the two bit
    /// reads, so it needs no new required method.
    #[inline]
    fn anticommutes_at(&self, i: usize, pauli: (bool, bool)) -> bool {
        (self.x_bit(i) & pauli.1) ^ (self.z_bit(i) & pauli.0)
    }
}
