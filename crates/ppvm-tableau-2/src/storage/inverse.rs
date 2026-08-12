// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The inverse tableau — which is *almost* free, and the part that is not.
//!
//! Stim's `TableauSimulator` stores `inv_state`, the inverse of the Clifford its
//! frame represents, because what measurement needs is a **row of the inverse**:
//! for `U` with `sᵢ = U ZᵢU†` and `dᵢ = U XᵢU†`,
//!
//! ```text
//! ω(P, sᵢ) = ω(U†PU, Zᵢ) = x-bit i of U†PU
//! ω(P, dᵢ) = ω(U†PU, Xᵢ) = z-bit i of U†PU
//! ```
//!
//! so the two anticommutation masks a measurement needs are the two bit planes
//! of one inverse row, and the `ℤ/4` phase the generalized algorithm needs is
//! that row's sign plus an `O(1)` correction
//! ([`Tableau::decomposition_phase`](crate::data::Tableau)).
//!
//! # The bits are already here
//!
//! Read those identities the other way round: the `x`-bit at site `j` of the
//! inverse row `ix_q = U†X_qU` is `ω(X_q, s_j)`, i.e. the **Z bit of stabilizer
//! `j` at qubit `q`** — and "bit `j` of the stabilizer-Z quadrant's major `q`" is
//! one contiguous major of the canonical column-major arena. All four planes of a
//! qubit's two inverse rows are forward majors:
//!
//! | inverse row | its X plane | its Z plane |
//! |:--|:--|:--|
//! | `ix_q = U†X_qU` | `major(Stab, Z, q)` | `major(Destab, Z, q)` |
//! | `iz_q = U†Z_qU` | `major(Stab, X, q)` | `major(Destab, X, q)` |
//!
//! This is the same fact that makes Stim's `Tableau::inverse()` a quadrant
//! transpose. Here the forward frame is *already* stored transposed relative to
//! its rows, so the canonical orientation is the inverse tableau held
//! row-contiguously — exactly Stim's `inv_state` layout. **No second bit matrix
//! is allocated and none is maintained**: every forward gate updates the
//! inverse's bits as a side effect of updating its own columns.
//!
//! # What is not free
//!
//! Signs. `U†X_qU` is a conjugate of a Hermitian Pauli, so it is Hermitian, so
//! its `ℤ/4` phase is `0` or `2` — one bit per row, `2n` bits in total, in
//! [`InverseSigns`]. Those the gates have to maintain, which is what
//! [`crate::inverse`] does.
//!
//! # Validity
//!
//! [`TableauData::inverse_valid`] is the escape hatch for frame mutations whose
//! inverse rule is not (yet) implemented: the public `SymplecticColumns` /
//! `PhaseTrack` primitives, which move bits and signs *independently* and so do
//! not correspond to any single Clifford, and `row_multiply`, which can leave a
//! generator non-Hermitian. When the flag is clear every reader falls back to the
//! row fold it replaced, so a missing rule costs speed and never correctness.

use super::{Half, Plane, TableauData, blocks};

/// Which of a qubit's two inverse rows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum InvRow {
    /// `ix_q = U† X_q U`.
    X,
    /// `iz_q = U† Z_q U`.
    Z,
}

/// The `2n` inverse-row signs, plus the working planes their updates need.
///
/// Excluded from [`TableauData`]'s equality and hashing, which compare the arena
/// alone: these are a *derived* cache of the forward frame, so two frames with
/// identical quadrants and phases have identical inverse signs whenever both are
/// valid — and comparing them would otherwise make a frame unequal to itself
/// across an invalidation.
#[derive(Debug)]
pub(crate) struct InverseSigns {
    /// Bit `q` is set when `U†X_qU` carries phase `2`.
    x: Vec<u64>,
    /// Bit `q` is set when `U†Z_qU` carries phase `2`.
    z: Vec<u64>,
    /// Working planes, `stride` words each, grown on first use and reused for
    /// the life of the frame: the `CY` rule's intermediate row and the
    /// projection's pivot / site / selector planes.
    scratch: Vec<u64>,
    /// Whether the signs describe the current frame.
    valid: bool,
}

impl InverseSigns {
    /// All-`+1` signs, which is the identity frame's inverse (`U = I`).
    pub(crate) fn identity(stride: usize) -> Self {
        Self {
            x: vec![0; stride],
            z: vec![0; stride],
            scratch: Vec::new(),
            valid: true,
        }
    }

    /// Restore the identity signs in place.
    pub(crate) fn reset(&mut self) {
        self.x.fill(0);
        self.z.fill(0);
        self.valid = true;
    }
}

/// Cloning drops the working planes: they are scratch, so a fork that inherits
/// them would pay `n/8` bytes and an allocation for buffers the clone has not
/// used yet. The next update re-grows them ([`TableauData::take_inv_scratch`]).
impl Clone for InverseSigns {
    fn clone(&self) -> Self {
        Self {
            x: self.x.clone(),
            z: self.z.clone(),
            scratch: Vec::new(),
            valid: self.valid,
        }
    }
}

impl TableauData {
    /// Whether the inverse-row signs currently describe this frame.
    #[inline]
    pub(crate) fn inverse_valid(&self) -> bool {
        self.inverse.valid
    }

    /// Abandon the inverse signs; readers fall back to the row fold.
    #[inline]
    pub(crate) fn invalidate_inverse(&mut self) {
        self.inverse.valid = false;
    }

    /// Declare the signs current again, after a rebuild.
    #[cfg(test)]
    #[inline]
    pub(crate) fn revalidate_inverse(&mut self) {
        self.inverse.valid = true;
    }

    /// Detach `planes` zeroed working planes of `stride` words each.
    ///
    /// Detached rather than borrowed because every user also reads the arena,
    /// which lives in the same struct; the pair is
    /// [`Self::restore_inv_scratch`], which must run before the next update.
    pub(crate) fn take_inv_scratch(&mut self, planes: usize) -> Vec<u64> {
        let want = planes * self.stride;
        let mut scratch = std::mem::take(&mut self.inverse.scratch);
        if scratch.len() < want {
            scratch.resize(want, 0);
        }
        scratch[..want].fill(0);
        scratch
    }

    /// Hand the working planes back for the next update to reuse.
    #[inline]
    pub(crate) fn restore_inv_scratch(&mut self, scratch: Vec<u64>) {
        self.inverse.scratch = scratch;
    }

    /// One family's sign plane, bit-indexed by qubit.
    ///
    /// The `ph` argument of a [`blocks`] kernel driven over inverse rows — see
    /// [`Tableau::project_inverse`](crate::data::Tableau).
    #[inline]
    pub(crate) fn inv_sign_plane_mut(&mut self, row: InvRow) -> &mut [u64] {
        match row {
            InvRow::X => &mut self.inverse.x,
            InvRow::Z => &mut self.inverse.z,
        }
    }

    /// The `ℤ/4` phase of inverse row `row` at qubit `q` — `0` or `2`.
    #[inline]
    pub(crate) fn inv_sign(&self, row: InvRow, q: usize) -> u8 {
        let plane = match row {
            InvRow::X => &self.inverse.x,
            InvRow::Z => &self.inverse.z,
        };
        u8::from(Self::bit(plane, q)) << 1
    }

    /// Overwrite [`Self::inv_sign`]. `phase` must be even — an inverse row is
    /// Hermitian, so its phase is a pure sign.
    #[inline]
    pub(crate) fn set_inv_sign(&mut self, row: InvRow, q: usize, phase: u8) {
        debug_assert_eq!(
            phase % 2,
            0,
            "inverse row {row:?} at qubit {q} is Hermitian, so its phase cannot be imaginary"
        );
        let plane = match row {
            InvRow::X => &mut self.inverse.x,
            InvRow::Z => &mut self.inverse.z,
        };
        Self::set_bit(plane, q, phase & 2 == 2);
    }

    /// The two bit planes of one inverse row, both site-indexed and both
    /// contiguous — see the module note's table.
    ///
    /// Column-major only: in the transposed orientation these majors hold a
    /// generator's bits, not an inverse row's.
    #[inline]
    pub(crate) fn inv_planes(&self, row: InvRow, q: usize) -> (&[u64], &[u64]) {
        debug_assert_eq!(self.orientation(), super::Orientation::ColumnMajor);
        let plane = match row {
            InvRow::X => Plane::Z,
            InvRow::Z => Plane::X,
        };
        (
            self.major(Half::Stab, plane, q),
            self.major(Half::Destab, plane, q),
        )
    }

    /// The `ℤ/4` phase of `U†PU` for a single-site Pauli `P` at qubit `q`.
    ///
    /// `X` and `Z` are stored signs. `Y_q = i·X_q·Z_q`, so its image is
    /// `i·ix_q·iz_q` and the phase picks up that product's `g`-rule — one pass
    /// over the four planes above, no allocation and no writes.
    pub(crate) fn inv_phase_of(&self, q: usize, pauli: ppvm_traits_2::Pauli) -> u8 {
        use ppvm_traits_2::Pauli;
        match pauli {
            Pauli::I => 0,
            Pauli::X => self.inv_sign(InvRow::X, q),
            Pauli::Z => self.inv_sign(InvRow::Z, q),
            Pauli::Y => (1 + self.inv_pair_phase((InvRow::X, q), (InvRow::Z, q))) % 4,
        }
    }

    /// The `ℤ/4` phase of the product of two inverse rows, in the given order:
    /// `phase(ia) + phase(ib) + g`, with `g` the Aaronson–Gottesman rule over
    /// their bits.
    pub(crate) fn inv_pair_phase(&self, a: (InvRow, usize), b: (InvRow, usize)) -> u8 {
        let (ax, az) = self.inv_planes(a.0, a.1);
        let (bx, bz) = self.inv_planes(b.0, b.1);
        let g = blocks::row_multiply_phase(ax, az, bx, bz);
        (self.inv_sign(a.0, a.1) + self.inv_sign(b.0, b.1) + g) % 4
    }

    /// The `ℤ/4` phase of the ordered product of three inverse rows.
    ///
    /// Only `CY` needs this (`Φ(X_c·Y_t) = ix_c · i·ix_t·iz_t`), and it is the one
    /// update that has to materialize an intermediate row — hence the scratch
    /// planes in [`InverseSigns`].
    pub(crate) fn inv_triple_phase(
        &mut self,
        a: (InvRow, usize),
        b: (InvRow, usize),
        c: (InvRow, usize),
    ) -> u8 {
        let head = self.inv_pair_phase(a, b);
        let mut scratch = self.take_inv_scratch(2);
        let (sx, rest) = scratch.split_at_mut(self.stride);
        let sz = &mut rest[..self.stride];
        {
            let (ax, az) = self.inv_planes(a.0, a.1);
            sx.copy_from_slice(ax);
            sz.copy_from_slice(az);
        }
        let (bx, bz) = self.inv_planes(b.0, b.1);
        blocks::row_multiply(sx, sz, bx, bz);
        let (cx, cz) = self.inv_planes(c.0, c.1);
        let g = blocks::row_multiply_phase(sx, sz, cx, cz);
        let phase = (head + self.inv_sign(c.0, c.1) + g) % 4;
        self.restore_inv_scratch(scratch);
        phase
    }
}
