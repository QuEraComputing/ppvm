// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`TableauData`]: the frame's physical storage — one aligned, contiguous,
//! bit-packed allocation divided by computed offsets.
//!
//! Design: `docs/design/tableau-data-structure.md` §"Physical storage",
//! §"Column-major orientation", §"Temporary transposition", §"Loss ownership".
//!
//! # Why this replaced `Vec<Row<A>>`
//!
//! The shipped-then-replaced representation was one `BitArray<A>` pair per
//! generator, `A` a *compile-time* width. Two costs followed, and together they
//! were the whole distance to Stim:
//!
//! 1. Every Clifford gate walked all `2n` generators to touch one machine word
//!    per plane per row — `O(n)` **strided** work for a one-qubit gate, where
//!    Stim's equivalent is a contiguous SIMD sweep.
//! 2. The stride was `size_of::<Row<A>>()`, set by `A` rather than by `n`. An
//!    85-qubit frame stored at `U2048` spent 87 KB carrying 3.7 KB of live bits,
//!    and the gate loop touched all of it.
//!
//! Here the *generator* dimension is contiguous for a fixed qubit, so `h(q)` is
//! two `2n`-bit sweeps over words that are adjacent in memory, and the stride is
//! `n.div_ceil(64)` rounded to a block — a runtime quantity, so there is no
//! compile-time qubit cap and no padding beyond one block.
//!
//! # Four square quadrants
//!
//! The logical object is a `2n × n` X matrix and a `2n × n` Z matrix. They are
//! stored as **four square `n × n` quadrants** — destabilizer-X,
//! destabilizer-Z, stabilizer-X, stabilizer-Z — because a rectangular bit-matrix
//! transpose is not an in-place permutation while a square one is, and the
//! `TransposedTableau` guard needs the
//! orientation change to be in-place. Stim splits its tableau for exactly this
//! reason (`do_transpose_quadrants`).
//!
//! Each quadrant is `n` majors of [`TableauData::stride`] words. In the
//! canonical [`Orientation::ColumnMajor`] the major index is the **qubit** and
//! the bit index within a major is the **generator**; under the guard those swap.
//!
//! # The `ℤ/4` phase is two bit planes
//!
//! A generator's phase is `ℤ/4` (`0: +1, 1: +i, 2: −1, 3: −i`), but every
//! Clifford gate only ever adds `2` — it flips the *sign*. Storing the phase as
//! one `u8` per generator would make that flip a scalar loop over `2n` bytes and
//! give back everything the column-major layout just bought. So the phase is
//! split into a low plane (the `i` bit) and a high plane (the sign bit), both
//! bit-packed and generator-indexed: a Clifford sign update is
//! `phase_hi ^= <predicate>`, one contiguous sweep, and the full `ℤ/4` addition
//! a row multiply needs is the two-plane ripple `carry = lo & rhs_lo;
//! lo ^= rhs_lo; hi ^= rhs_hi ^ carry`.
//!
//! # Padding is zero, always
//!
//! Words past `n.div_ceil(64)` in a major, and bits past `n` in the last word,
//! are held at zero. Every kernel here is elementwise and every phase predicate
//! is guarded by an `AND` against an X or Z plane, so zero padding stays zero;
//! [`transpose`] gathers through a zero-filled scratch for the same reason. That
//! invariant is what lets equality and hashing compare the whole arena in bulk
//! instead of masking each ragged edge.

use std::hash::{Hash, Hasher};

pub(crate) mod blocks;
pub(crate) mod inverse;
pub(crate) mod transpose;

pub(crate) use inverse::{InvRow, InverseSigns};
pub(crate) use transpose::BITS_PER_WORD;

/// Words per allocation block. The arena is a `Vec<Block>`, so every region
/// start is 32-byte aligned and every major is a whole number of blocks.
pub(crate) const WORDS_PER_BLOCK: usize = 4;

/// A 32-byte aligned run of four machine words.
///
/// Allocating `Vec<Block>` rather than `Vec<u64>` is how the arena gets its
/// alignment without a custom allocator: `Vec` aligns to its element type.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(C, align(32))]
struct Block([u64; WORDS_PER_BLOCK]);

/// Which physical arrangement the X/Z quadrants are currently in.
///
/// [`ColumnMajor`](Orientation::ColumnMajor) is canonical and is the only
/// orientation a public method may return in; see [`TransposedTableau`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Orientation {
    /// Major index is the qubit, bit index within a major is the generator.
    /// Gates and the measurement pivot scan are contiguous.
    ColumnMajor,
    /// Major index is the generator, bit index is the qubit. Row multiplication
    /// and elimination are contiguous.
    RowMajor,
}

/// Which half of the frame a generator belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Half {
    /// Generators `0..n`.
    Destab = 0,
    /// Generators `n..2n`.
    Stab = 1,
}

/// Which symplectic bit plane.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Plane {
    /// The X plane.
    X = 0,
    /// The Z plane.
    Z = 1,
}

/// The five disjoint plane slices a two-qubit Clifford kernel sweeps:
/// `(X_a, Z_a, X_b, Z_b, phase-hi)` for one half of the frame.
pub(crate) type Gate2Slices<'a> = (
    &'a mut [u64],
    &'a mut [u64],
    &'a mut [u64],
    &'a mut [u64],
    &'a mut [u64],
);

/// Both halves of the frame, in the order every sweep visits them.
pub(crate) const HALVES: [Half; 2] = [Half::Destab, Half::Stab];

impl Half {
    /// The half a global generator index `0..2n` lives in, and its index within
    /// that half.
    #[inline]
    pub(crate) fn split(generator: usize, n_qubits: usize) -> (Self, usize) {
        if generator < n_qubits {
            (Half::Destab, generator)
        } else {
            (Half::Stab, generator - n_qubits)
        }
    }
}

// ─── TableauData ──────────────────────────────────────────────────────────

/// The frame's bits: four square X/Z quadrants, two phase planes, one loss
/// plane, in a single aligned allocation.
#[derive(Clone)]
pub(crate) struct TableauData {
    /// The whole arena. Region starts are computed, never stored.
    blocks: Vec<Block>,
    n_qubits: usize,
    /// Words per major, `n.div_ceil(64)` rounded up to a whole block.
    stride: usize,
    orientation: Orientation,
    /// Signs of the inverse tableau's rows — a derived cache whose bits live in
    /// the quadrants above. See [`inverse`]; excluded from equality and hashing.
    inverse: InverseSigns,
}

impl TableauData {
    /// The identity frame on `n_qubits`: destabilizer `i` is `X_i`, stabilizer
    /// `i` is `Z_i`, all phases `+1`, nothing lost.
    pub(crate) fn identity(n_qubits: usize) -> Self {
        let stride = Self::stride_for(n_qubits);
        let words = Self::total_words(n_qubits, stride);
        let mut data = Self {
            blocks: vec![Block([0; WORDS_PER_BLOCK]); words / WORDS_PER_BLOCK],
            n_qubits,
            stride,
            orientation: Orientation::ColumnMajor,
            inverse: InverseSigns::identity(stride),
        };
        data.write_identity();
        data
    }

    /// Restore the identity frame without reallocating.
    pub(crate) fn reset_to_identity(&mut self) {
        self.blocks.fill(Block([0; WORDS_PER_BLOCK]));
        self.orientation = Orientation::ColumnMajor;
        self.inverse.reset();
        self.write_identity();
    }

    /// Set the `2n` diagonal bits. In either orientation the identity frame is
    /// symmetric — bit `(q, q)` of the destabilizer-X and stabilizer-Z
    /// quadrants — so this does not consult [`Self::orientation`].
    fn write_identity(&mut self) {
        let (n, stride) = (self.n_qubits, self.stride);
        for q in 0..n {
            let dx = Self::major_start(n, stride, Half::Destab, Plane::X) + q * stride;
            let sz = Self::major_start(n, stride, Half::Stab, Plane::Z) + q * stride;
            let (word, bit) = (q / BITS_PER_WORD, q % BITS_PER_WORD);
            let words = self.words_mut();
            words[dx + word] |= 1u64 << bit;
            words[sz + word] |= 1u64 << bit;
        }
    }

    /// Words per major for an `n`-qubit frame: enough for `n` bits, rounded up
    /// to a whole 32-byte block, and never zero.
    #[inline]
    fn stride_for(n_qubits: usize) -> usize {
        n_qubits
            .div_ceil(BITS_PER_WORD)
            .next_multiple_of(WORDS_PER_BLOCK)
            .max(WORDS_PER_BLOCK)
    }

    /// Total arena size: four `n × stride` quadrants, two 2-major phase planes,
    /// one 1-major loss plane.
    #[inline]
    fn total_words(n_qubits: usize, stride: usize) -> usize {
        4 * n_qubits * stride + 5 * stride
    }

    /// Word offset of the first major of a quadrant.
    #[inline]
    fn major_start(n_qubits: usize, stride: usize, half: Half, plane: Plane) -> usize {
        (half as usize * 2 + plane as usize) * n_qubits * stride
    }

    /// Word offset of a phase-plane major. `hi` selects the sign plane.
    #[inline]
    fn phase_start(n_qubits: usize, stride: usize, half: Half, hi: bool) -> usize {
        4 * n_qubits * stride + (usize::from(hi) * 2 + half as usize) * stride
    }

    /// Word offset of the loss plane.
    #[allow(dead_code)]
    #[inline]
    fn loss_start(n_qubits: usize, stride: usize) -> usize {
        4 * n_qubits * stride + 4 * stride
    }

    /// Number of qubits the frame is over.
    #[inline]
    pub(crate) fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Words per major.
    #[inline]
    pub(crate) fn stride(&self) -> usize {
        self.stride
    }

    /// The current physical arrangement.
    #[inline]
    pub(crate) fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// The arena as machine words.
    #[inline]
    fn words(&self) -> &[u64] {
        // SAFETY: `Block` is `#[repr(C, align(32))]` around `[u64;
        // WORDS_PER_BLOCK]`, so it is exactly that array's size with no padding
        // and the same alignment requirement satisfied. A `&[Block]` of length
        // `k` therefore covers `k * WORDS_PER_BLOCK` initialised, aligned `u64`s
        // with the same provenance and lifetime.
        unsafe {
            std::slice::from_raw_parts(
                self.blocks.as_ptr().cast::<u64>(),
                self.blocks.len() * WORDS_PER_BLOCK,
            )
        }
    }

    /// The arena as mutable machine words.
    #[inline]
    fn words_mut(&mut self) -> &mut [u64] {
        // SAFETY: as [`Self::words`], and `&mut self` makes the borrow unique.
        unsafe {
            std::slice::from_raw_parts_mut(
                self.blocks.as_mut_ptr().cast::<u64>(),
                self.blocks.len() * WORDS_PER_BLOCK,
            )
        }
    }

    /// The word range of major `i` of a quadrant.
    #[inline]
    fn major_range(&self, half: Half, plane: Plane, i: usize) -> std::ops::Range<usize> {
        debug_assert!(i < self.n_qubits);
        let start = Self::major_start(self.n_qubits, self.stride, half, plane) + i * self.stride;
        start..start + self.stride
    }

    /// Major `i` of a quadrant. In [`Orientation::ColumnMajor`] `i` is a qubit
    /// and the bits are generator-indexed; in [`Orientation::RowMajor`] the
    /// reverse.
    #[inline]
    pub(crate) fn major(&self, half: Half, plane: Plane, i: usize) -> &[u64] {
        &self.words()[self.major_range(half, plane, i)]
    }

    /// Mutable [`Self::major`].
    #[inline]
    pub(crate) fn major_mut(&mut self, half: Half, plane: Plane, i: usize) -> &mut [u64] {
        let range = self.major_range(half, plane, i);
        &mut self.words_mut()[range]
    }

    /// The word range of a phase-plane major.
    #[inline]
    fn phase_range(&self, half: Half, hi: bool) -> std::ops::Range<usize> {
        let start = Self::phase_start(self.n_qubits, self.stride, half, hi);
        start..start + self.stride
    }

    /// A phase plane for one half, bit-indexed by generator-within-half.
    /// `hi` selects the sign bit (`+2`), `!hi` the `i` bit (`+1`).
    #[inline]
    pub(crate) fn phase_plane(&self, half: Half, hi: bool) -> &[u64] {
        &self.words()[self.phase_range(half, hi)]
    }

    /// Both phase planes of one half, `(low, high)`, mutably.
    #[inline]
    pub(crate) fn phase_planes_mut(&mut self, half: Half) -> (&mut [u64], &mut [u64]) {
        let ranges = [self.phase_range(half, false), self.phase_range(half, true)];
        let [lo, hi] = self
            .words_mut()
            .get_disjoint_mut(ranges)
            .expect("the low and high phase planes are disjoint regions");
        (lo, hi)
    }

    /// Mutable [`Self::phase_plane`].
    #[inline]
    pub(crate) fn phase_plane_mut(&mut self, half: Half, hi: bool) -> &mut [u64] {
        let range = self.phase_range(half, hi);
        &mut self.words_mut()[range]
    }

    /// Gather one plane's bits at qubit `addr0` over all generators of `half`,
    /// into a generator-indexed bit vector.
    ///
    /// Orientation-aware, and the reason it has to be: the measurement path may
    /// run under a [`Orientation::RowMajor`] guard held by an outer caller (one
    /// transpose for a whole `measure_all`), and in that orientation this is a
    /// *strided* read — one bit from each of `n` majors. That is `O(n)` scalar
    /// work against the `O(n²/64)` the elimination costs either way, so it is
    /// the cheap side of the trade; in the canonical orientation it is a
    /// `memcpy` of one contiguous column.
    pub(crate) fn gather_column(&self, half: Half, plane: Plane, addr0: usize, out: &mut [u64]) {
        debug_assert_eq!(out.len(), self.stride);
        match self.orientation {
            Orientation::ColumnMajor => out.copy_from_slice(self.major(half, plane, addr0)),
            Orientation::RowMajor => {
                // The `addr0` bit of `n` consecutive majors. Walking the majors
                // by stride and accumulating a whole output word before storing
                // it keeps this to one load and one shift per generator, where
                // going through `major`/`set_bit` re-derives the quadrant base
                // and re-checks two bounds on every bit.
                let words = self.words();
                let base = Self::major_start(self.n_qubits, self.stride, half, plane)
                    + addr0 / BITS_PER_WORD;
                let mask = 1u64 << (addr0 % BITS_PER_WORD);
                let mut src = base;
                for chunk in out.iter_mut() {
                    let mut acc = 0u64;
                    let end =
                        (src + self.stride * BITS_PER_WORD).min(base + self.stride * self.n_qubits);
                    let mut bit = 0;
                    while src < end {
                        acc |= u64::from(words[src] & mask != 0) << bit;
                        src += self.stride;
                        bit += 1;
                    }
                    *chunk = acc;
                }
            }
        }
    }

    /// Materialize generator `i` of `half` as a qubit-indexed `(x, z)` pair.
    ///
    /// The transpose of [`Self::gather_column`]: `O(n)` strided bit reads in the
    /// canonical orientation, a `memcpy` under a row guard. Folding `k` selected
    /// generators this way costs `k·n` bit reads against the guard's
    /// `O(n²/64)` word ops for the transpose *pair*, so it is the cheaper route
    /// whenever `k` is well below `n/16` — which is every measurement in a
    /// local code, where the decomposition weight does not scale with `n`.
    pub(crate) fn gather_row(&self, half: Half, i: usize, out_x: &mut [u64], out_z: &mut [u64]) {
        debug_assert_eq!(out_x.len(), self.stride);
        match self.orientation {
            Orientation::RowMajor => {
                out_x.copy_from_slice(self.major(half, Plane::X, i));
                out_z.copy_from_slice(self.major(half, Plane::Z, i));
            }
            Orientation::ColumnMajor => {
                // The mirror image of [`Self::gather_column`]'s strided branch,
                // and word-oriented for the same reason: bit `i` of `n`
                // consecutive majors, one output word at a time.
                let words = self.words();
                let (word, mask) = (i / BITS_PER_WORD, 1u64 << (i % BITS_PER_WORD));
                for (plane, out) in [(Plane::X, out_x), (Plane::Z, out_z)] {
                    let base = Self::major_start(self.n_qubits, self.stride, half, plane) + word;
                    let limit = base + self.stride * self.n_qubits;
                    let mut src = base;
                    for chunk in out.iter_mut() {
                        let mut acc = 0u64;
                        let end = (src + self.stride * BITS_PER_WORD).min(limit);
                        let mut bit = 0;
                        while src < end {
                            acc |= u64::from(words[src] & mask != 0) << bit;
                            src += self.stride;
                            bit += 1;
                        }
                        *chunk = acc;
                    }
                }
            }
        }
    }

    /// The four X/Z quadrants as one contiguous byte range.
    ///
    /// Zero padding makes this a faithful digest of the frame's bits: two frames
    /// with the same logical content have byte-identical quadrants.
    pub(crate) fn xz_bytes(&self) -> &[u8] {
        debug_assert_eq!(self.orientation, Orientation::ColumnMajor);
        words_as_bytes(&self.words()[..4 * self.n_qubits * self.stride])
    }

    /// The per-qubit loss plane. Never transposed — it is a vector, not a
    /// matrix, and its index is the qubit in both orientations.
    ///
    /// **Reserved, not yet wired.** `tableau-data-structure.md` §"Loss
    /// ownership" puts the loss plane on `Tableau`, but `GeneralizedTableau`
    /// still owns `is_lost: Vec<bool>` and ~99 call sites read it. The region
    /// is allocated and transposition-safe so that move is a mechanical
    /// follow-up rather than another layout change.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn loss_plane(&self) -> &[u64] {
        let start = Self::loss_start(self.n_qubits, self.stride);
        &self.words()[start..start + self.stride]
    }

    /// Mutable [`Self::loss_plane`].
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn loss_plane_mut(&mut self) -> &mut [u64] {
        let start = Self::loss_start(self.n_qubits, self.stride);
        let stride = self.stride;
        &mut self.words_mut()[start..start + stride]
    }

    // ─── Disjoint borrows for the gate kernels ────────────────────────────

    /// The `(X, Z, phase-hi)` slices a one-qubit Clifford sweeps, for one half.
    ///
    /// Column-major only: `q` indexes a qubit column.
    #[inline]
    pub(crate) fn gate1_mut(
        &mut self,
        half: Half,
        q: usize,
    ) -> (&mut [u64], &mut [u64], &mut [u64]) {
        debug_assert_eq!(self.orientation, Orientation::ColumnMajor);
        let ranges = [
            self.major_range(half, Plane::X, q),
            self.major_range(half, Plane::Z, q),
            self.phase_range(half, true),
        ];
        let [x, z, ph] = self
            .words_mut()
            .get_disjoint_mut(ranges)
            .expect("quadrant and phase-plane regions are disjoint by construction");
        (x, z, ph)
    }

    /// The `(X_a, Z_a, X_b, Z_b, phase-hi)` slices a two-qubit Clifford sweeps,
    /// for one half. Requires `a != b`.
    #[inline]
    pub(crate) fn gate2_mut(&mut self, half: Half, a: usize, b: usize) -> Gate2Slices<'_> {
        debug_assert_eq!(self.orientation, Orientation::ColumnMajor);
        debug_assert_ne!(a, b, "two-qubit gate needs distinct qubits");
        let ranges = [
            self.major_range(half, Plane::X, a),
            self.major_range(half, Plane::Z, a),
            self.major_range(half, Plane::X, b),
            self.major_range(half, Plane::Z, b),
            self.phase_range(half, true),
        ];
        let [xa, za, xb, zb, ph] = self
            .words_mut()
            .get_disjoint_mut(ranges)
            .expect("distinct qubit columns and the phase plane are disjoint");
        (xa, za, xb, zb, ph)
    }

    // ─── Logical bit access ───────────────────────────────────────────────

    /// Read bit `i` of `words`.
    #[inline]
    pub(crate) fn bit(words: &[u64], i: usize) -> bool {
        words[i / BITS_PER_WORD] >> (i % BITS_PER_WORD) & 1 == 1
    }

    /// Write bit `i` of `words`.
    #[inline]
    pub(crate) fn set_bit(words: &mut [u64], i: usize, value: bool) {
        let mask = 1u64 << (i % BITS_PER_WORD);
        let word = &mut words[i / BITS_PER_WORD];
        if value {
            *word |= mask;
        } else {
            *word &= !mask;
        }
    }

    /// The `(major, bit)` pair addressing logical cell `(generator, qubit)` in
    /// the current orientation.
    #[inline]
    fn address(&self, generator: usize, qubit: usize) -> (Half, usize, usize) {
        let (half, g) = Half::split(generator, self.n_qubits);
        match self.orientation {
            Orientation::ColumnMajor => (half, qubit, g),
            Orientation::RowMajor => (half, g, qubit),
        }
    }

    /// The X bit of `(generator, qubit)`.
    #[inline]
    pub(crate) fn x_bit(&self, generator: usize, qubit: usize) -> bool {
        let (half, major, bit) = self.address(generator, qubit);
        Self::bit(self.major(half, Plane::X, major), bit)
    }

    /// The Z bit of `(generator, qubit)`.
    #[inline]
    pub(crate) fn z_bit(&self, generator: usize, qubit: usize) -> bool {
        let (half, major, bit) = self.address(generator, qubit);
        Self::bit(self.major(half, Plane::Z, major), bit)
    }

    /// Write the X and Z bits of `(generator, qubit)`. Test-only: the frame's
    /// own algorithms write whole columns or whole generators.
    #[cfg(test)]
    #[inline]
    pub(crate) fn set_xz_bit(&mut self, generator: usize, qubit: usize, x: bool, z: bool) {
        let (half, major, bit) = self.address(generator, qubit);
        Self::set_bit(self.major_mut(half, Plane::X, major), bit, x);
        Self::set_bit(self.major_mut(half, Plane::Z, major), bit, z);
    }

    /// The `ℤ/4` phase of generator `i` within `half`.
    #[inline]
    pub(crate) fn phase_of(&self, half: Half, i: usize) -> u8 {
        u8::from(Self::bit(self.phase_plane(half, false), i))
            | (u8::from(Self::bit(self.phase_plane(half, true), i)) << 1)
    }

    /// Overwrite the `ℤ/4` phase of generator `i` within `half`.
    #[inline]
    pub(crate) fn set_phase_of(&mut self, half: Half, i: usize, phase: u8) {
        Self::set_bit(self.phase_plane_mut(half, false), i, phase & 1 == 1);
        Self::set_bit(self.phase_plane_mut(half, true), i, phase & 2 == 2);
    }

    /// The `ℤ/4` phase of a global generator index `0..2n`.
    #[inline]
    pub(crate) fn phase(&self, generator: usize) -> u8 {
        let (half, g) = Half::split(generator, self.n_qubits);
        self.phase_of(half, g)
    }

    /// Overwrite the `ℤ/4` phase of a global generator index `0..2n`.
    /// Test-only; the frame writes phases per half.
    #[cfg(test)]
    #[inline]
    pub(crate) fn set_phase(&mut self, generator: usize, phase: u8) {
        let (half, g) = Half::split(generator, self.n_qubits);
        self.set_phase_of(half, g, phase);
    }

    /// Multiply the phased Pauli `(src_x, src_z, src_phase)` into generator `i`
    /// of `half` using the Aaronson–Gottesman `g`-rule.
    ///
    /// Row-major only: the generator's bits must be contiguous. This is the
    /// [`StabilizerFrame::row_multiply`](ppvm_traits_2::StabilizerFrame::row_multiply)
    /// primitive, and the elimination step of the measurement projection.
    pub(crate) fn multiply_row_by(
        &mut self,
        half: Half,
        i: usize,
        src_x: &[u64],
        src_z: &[u64],
        src_phase: u8,
    ) {
        debug_assert_eq!(self.orientation, Orientation::RowMajor);
        let ranges = [
            self.major_range(half, Plane::X, i),
            self.major_range(half, Plane::Z, i),
        ];
        let [dst_x, dst_z] = self
            .words_mut()
            .get_disjoint_mut(ranges)
            .expect("the X and Z quadrants are disjoint regions");
        let g = blocks::row_multiply(dst_x, dst_z, src_x, src_z);
        let phase = (self.phase_of(half, i) + g + src_phase) % 4;
        self.set_phase_of(half, i, phase);
    }

    // ─── Loss plane ───────────────────────────────────────────────────────

    /// Whether qubit `q` has been lost. See [`Self::loss_plane`] — reserved.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn is_lost(&self, q: usize) -> bool {
        Self::bit(self.loss_plane(), q)
    }

    /// Mark or clear qubit `q`'s loss bit. See [`Self::loss_plane`] — reserved.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn set_lost(&mut self, q: usize, value: bool) {
        Self::set_bit(self.loss_plane_mut(), q, value);
    }

    // ─── Orientation ──────────────────────────────────────────────────────

    /// Transpose all four quadrants and flip [`Self::orientation`].
    ///
    /// A physical reordering of the same logical bits: `x_bit`, `z_bit`,
    /// `phase` and `is_lost` all read the same values afterwards.
    pub(crate) fn transpose_quadrants(&mut self) {
        let (n, stride) = (self.n_qubits, self.stride);
        for half in [Half::Destab, Half::Stab] {
            for plane in [Plane::X, Plane::Z] {
                let start = Self::major_start(n, stride, half, plane);
                let words = &mut self.words_mut()[start..start + n * stride];
                transpose::transpose_square(words, stride, n);
            }
        }
        self.orientation = match self.orientation {
            Orientation::ColumnMajor => Orientation::RowMajor,
            Orientation::RowMajor => Orientation::ColumnMajor,
        };
    }
}

/// Bulk comparison of the whole arena.
///
/// Sound because padding is held at zero (module note) and because two frames
/// with the same logical content necessarily share `n_qubits`, hence `stride`,
/// hence arena length. Orientation is *not* excluded: the caller cannot observe
/// a transposed frame, since [`TransposedTableau`] holds `&mut` for its whole
/// lifetime and the borrow checker forbids a shared read through it.
impl PartialEq for TableauData {
    fn eq(&self, other: &Self) -> bool {
        debug_assert_eq!(self.orientation, Orientation::ColumnMajor);
        debug_assert_eq!(other.orientation, Orientation::ColumnMajor);
        self.n_qubits == other.n_qubits && self.blocks == other.blocks
    }
}

impl Eq for TableauData {}

impl Hash for TableauData {
    fn hash<Hr: Hasher>(&self, state: &mut Hr) {
        debug_assert_eq!(self.orientation, Orientation::ColumnMajor);
        self.n_qubits.hash(state);
        state.write(words_as_bytes(self.words()));
    }
}

/// Reinterpret words as bytes for a byte-oriented `Hasher::write`.
#[inline]
fn words_as_bytes(words: &[u64]) -> &[u8] {
    // SAFETY: `u64` has no padding or invalid bit patterns, so any initialised
    // `[u64]` is a valid `[u8]` of eight times the length; the returned slice
    // borrows from `words` and is read-only.
    unsafe { std::slice::from_raw_parts(words.as_ptr().cast::<u8>(), std::mem::size_of_val(words)) }
}

impl std::fmt::Debug for TableauData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TableauData")
            .field("n_qubits", &self.n_qubits)
            .field("stride", &self.stride)
            .field("orientation", &self.orientation)
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_frame_is_the_diagonal() {
        let data = TableauData::identity(70);
        for g in 0..140 {
            for q in 0..70 {
                let (want_x, want_z) = match (g < 70, g % 70 == q) {
                    (true, true) => (true, false),
                    (false, true) => (false, true),
                    _ => (false, false),
                };
                assert_eq!(data.x_bit(g, q), want_x, "x({g}, {q})");
                assert_eq!(data.z_bit(g, q), want_z, "z({g}, {q})");
            }
            assert_eq!(data.phase(g), 0);
        }
    }

    #[test]
    fn arena_is_block_aligned() {
        for n in [1usize, 63, 64, 65, 200] {
            let data = TableauData::identity(n);
            assert_eq!(data.words().as_ptr() as usize % 32, 0, "n={n}");
            assert!(data.stride().is_multiple_of(WORDS_PER_BLOCK), "n={n}");
        }
    }

    #[test]
    fn transpose_round_trips_and_preserves_logical_bits() {
        for n in [1usize, 7, 64, 65, 100] {
            let mut data = TableauData::identity(n);
            // An asymmetric pattern: a symmetric one could survive a no-op.
            for g in 0..2 * n {
                for q in 0..n {
                    data.set_xz_bit(g, q, (g * 5 + q * 3) % 7 == 0, (g + 2 * q) % 5 == 0);
                }
                data.set_phase(g, (g % 4) as u8);
            }
            let before = data.clone();

            data.transpose_quadrants();
            assert_eq!(data.orientation(), Orientation::RowMajor);
            for g in 0..2 * n {
                for q in 0..n {
                    assert_eq!(data.x_bit(g, q), before.x_bit(g, q), "n={n} x({g}, {q})");
                    assert_eq!(data.z_bit(g, q), before.z_bit(g, q), "n={n} z({g}, {q})");
                }
                assert_eq!(data.phase(g), before.phase(g), "n={n} phase({g})");
            }

            data.transpose_quadrants();
            assert_eq!(data.orientation(), Orientation::ColumnMajor);
            assert_eq!(data, before, "n={n}: transpose is not an involution");
        }
    }

    #[test]
    fn padding_never_affects_equality_or_hashing() {
        use std::hash::BuildHasher;
        // n = 65 leaves 63 padding bits in the second word and two padding
        // words in the block; the identity frame must hash the same as one
        // reached by a round trip through the transposed orientation.
        let a = TableauData::identity(65);
        let mut b = TableauData::identity(65);
        b.transpose_quadrants();
        b.transpose_quadrants();
        assert_eq!(a, b);
        let hasher = fxhash::FxBuildHasher::default();
        assert_eq!(hasher.hash_one(&a), hasher.hash_one(&b));
    }

    #[test]
    fn loss_plane_survives_transposition() {
        let mut data = TableauData::identity(70);
        data.set_lost(3, true);
        data.set_lost(69, true);
        data.transpose_quadrants();
        assert!(data.is_lost(3) && data.is_lost(69) && !data.is_lost(4));
        data.transpose_quadrants();
        assert!(data.is_lost(3) && data.is_lost(69) && !data.is_lost(4));
    }
}
