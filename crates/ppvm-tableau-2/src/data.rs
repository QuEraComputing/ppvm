// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The tableau representation: the packed [`Row`], the Clifford frame
//! [`Tableau`], the amplitude store [`Amplitudes`], and the coefficient-aware
//! [`GeneralizedTableau`].
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits"
//! (the frame is `SymplecticColumns` + `PhaseTrack` + `StabilizerFrame`),
//! §"A third instantiation: the generalized tableau" (frame + amplitude
//! vector), §"Tableau indexability" (the frame is itself `Indexable`), and
//! `word-data-structures.md` §"`PauliWord` packed representation" (the packed
//! X/Z plane layout the rows reuse).
//!
//! Ported from `ppvm-tableau/src/{data,sparsevec,tableau_index}.rs`. The packed
//! plane layout, the `2n`-row ordering, the sort-merge branch coalesce and every
//! magnitude-cutoff comparison are reproduced verbatim so the observable
//! behaviour — including the *order* of the amplitude vector, which is public —
//! is byte-for-byte the old crate's.

use std::fmt::Debug;
use std::hash::{BuildHasher, Hash, Hasher};
use std::marker::PhantomData;
use std::ops::{BitAnd, BitOrAssign, BitXor, Shl};
use std::sync::atomic::{AtomicU64, Ordering};

use num::PrimInt;
use num::complex::Complex64;
use ppvm_traits_2::{Indexable, Pauli, Scale, Support};

use crate::storage::{BITS_PER_WORD, HALVES, Half, Orientation, Plane, TableauData, blocks};

/// Bit-string index type addressing one branch of the amplitude vector.
///
/// Blanket-implemented for every primitive (and `bnum`-style) unsigned integer
/// supporting the required bit operations. Pick `usize` for ≤ 64 qubits, `u128`
/// for ≤ 128, `bnum::types::U256`/`U512` beyond. Ported from
/// `ppvm-tableau/src/tableau_index.rs` (`TableauIndex`); renamed to `Bitstring`
/// because in the `-2` tower the amplitude key is the *bitstring* of the
/// graded algebra `C[Bitstring]`, not the tableau
/// (`traits-2-configuration-and-hashing.md` §"A third instantiation").
pub trait Bitstring:
    PartialEq
    + Eq
    + Hash
    + Copy
    + Debug
    + Send
    + Sync
    + From<u8>
    + Shl<usize, Output = Self>
    + BitOrAssign<Self>
    + BitAnd<Self, Output = Self>
    + BitXor<Output = Self>
    + PrimInt
{
}

impl<I> Bitstring for I where
    I: PartialEq
        + Eq
        + Hash
        + Copy
        + Debug
        + Send
        + Sync
        + From<u8>
        + Shl<usize, Output = Self>
        + BitOrAssign<Self>
        + BitAnd<Self, Output = Self>
        + BitXor<Output = Self>
        + PrimInt
{
}

/// `i^φ` for `φ ∈ ℤ/4` — the phase convention `0: +1, 1: +i, 2: −1, 3: −i`.
pub(crate) const COMPLEX_PHASE_CONVERSION: [Complex64; 4] = [
    Complex64 { re: 1.0, im: 0.0 },
    Complex64 { re: 0.0, im: 1.0 },
    Complex64 { re: -1.0, im: 0.0 },
    Complex64 { re: 0.0, im: -1.0 },
];

/// The [`Tableau::hash_cache`] sentinel meaning "not yet computed".
///
/// A frame whose finalized digest happens to equal this value is simply
/// recomputed on every `key_hash()` — the result is still the true digest, only
/// the caching is skipped, so this is a 1-in-2⁶⁴ perf non-event with **no**
/// correctness consequence.
pub(crate) const HASH_UNCACHED: u64 = 0;

/// Symplectic inner product of two bitstrings — the count of shared set bits,
/// used in the stabilizer phase calculations.
#[inline]
pub fn symplectic_inner<I: Bitstring>(alpha: I, beta: I) -> u32 {
    (alpha & beta).count_ones()
}

/// Widen a generator-indexed bit plane into a branch index.
///
/// The frame's masks (`destab_anticomm_bits`, `stab_anticomm_bits`, the
/// odd-phase destabilizer mask) are produced as contiguous `u64` planes but
/// consumed as `I`, which may be a `bnum` big integer wider than a machine word.
/// Set bits are OR'd in one at a time — exactly the `mask |= one << i` the
/// replaced per-row scans did, so a wide `I` pays no more than it did before.
#[inline]
pub(crate) fn bits_to_index<I: Bitstring>(words: &[u64], n_bits: usize) -> I {
    let mut acc = I::zero();
    let one = I::one();
    for (w, &word) in words.iter().enumerate() {
        if word == 0 {
            continue;
        }
        let base = w * BITS_PER_WORD;
        let mut rest = word;
        while rest != 0 {
            let i = base + rest.trailing_zeros() as usize;
            if i >= n_bits {
                break;
            }
            acc |= one << i;
            rest &= rest - 1;
        }
    }
    acc
}

/// The per-coefficient phase contribution, given the destabilizer
/// anticommutation bits and the pre-hoisted odd-phase destabilizer mask.
///
/// A pure function of its arguments (two popcounts, no memory access): the mask
/// is computed **once** per gate/measurement by
/// [`GeneralizedTableau::odd_phase_destabilizer_mask`] and folded in here, which
/// is what turns the `O(n·m)` per-coefficient destabilizer walk into `O(n + m)`
/// per gate.
///
/// # Why this formula
///
/// It is not a convention: `lean/PPVM/Tableau/BranchPhase.lean` *derives* it.
/// Modelling the amplitude basis as `|j⟩ = ∏_l d_l^{j_l}|ψ₀⟩`,
/// `frameOp_eq_shiftOp` proves that the frame-conjugated Pauli
/// `i^{phase_decomp}·D_L·S_G` acting on amplitudes is exactly
/// `phase_decomp + 2·⟨destab_anticomm, j⟩ + 2·popcount(j ∧ stab_anticomm ∧ mask)`
/// — the `⟨destab, ·⟩` term read at the **original** index, and the mask term
/// coming from `d_l² = (−1)^{phase_l}` on doubly-selected generators.
/// `selfInverse_branchPhase` then discharges the `SelfInverse` hypothesis that
/// every case-a/case-b theorem in `lean/PPVM/Tableau/Projection.lean` is stated
/// under, given the `ℤ/2` frame identity
/// `phase_decomp + ⟨destab, stab⟩ + popcount(stab ∧ mask) ≡ 0` — which
/// `frameOp_involutive_iff` shows is just `M² = I`, and which
/// `crates/ppvm-conformance-2/tests/tableau_lean.rs` checks on every
/// decomposition of a random Clifford+`T` sweep.
#[inline]
pub(crate) fn compute_phase_with_mask_static<I: Bitstring>(
    destab_anticomm_bits: I,
    basis_index: I,
    stab_anticomm_bits: I,
    odd_phase_mask: I,
) -> u8 {
    let mut phase = (2 * symplectic_inner(destab_anticomm_bits, basis_index) as u8) % 4;
    let active = basis_index & stab_anticomm_bits;
    let parity = (active & odd_phase_mask).count_ones() % 2;
    phase = (phase + 2 * parity as u8) % 4;
    phase
}

// ─── ScratchRow ───────────────────────────────────────────────────────────

/// A heap-allocated phased Pauli word used as an *accumulator* by the frame's
/// row algebra — the residual `p_word` of
/// [`GeneralizedTableau::compute_decomposition`], the product built by
/// [`Tableau::get_deterministic_outcome`], and the pivot copy the measurement
/// projection multiplies into its neighbours.
///
/// It is not the tableau's storage. The frame's own generators live in
/// [`TableauData`]'s four square quadrants; this is the one place a *single*
/// generator has to exist on its own, and it uses the same
/// `(x-words, z-words, ℤ/4 phase)` shape so [`blocks::row_multiply`] serves both.
///
/// # Why not `Row<A>` any more
///
/// The replaced `Row<A>` was `Copy` over a compile-time-sized `BitArray<A>`,
/// which is what forced the whole tableau to be a `Vec<Row<A>>` — and with it
/// the compile-time qubit cap and the padding blow-up documented in
/// [`crate::storage`]. A scratch row is allocated once per measurement scratch
/// and reused, so it does not need to be `Copy`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ScratchRow {
    /// X bits, qubit-indexed. Length is the frame's stride; padding stays zero.
    pub(crate) x: Vec<u64>,
    /// Z bits, qubit-indexed.
    pub(crate) z: Vec<u64>,
    /// Phase in `ℤ/4` with the convention `0: +1, 1: +i, 2: −1, 3: −i`.
    pub(crate) phase: u8,
}

impl ScratchRow {
    /// The identity row `+I…I`, sized for a frame of the given stride.
    #[inline]
    pub(crate) fn zeroed(stride: usize) -> Self {
        Self {
            x: vec![0; stride],
            z: vec![0; stride],
            phase: 0,
        }
    }

    /// Write the Pauli at site `i`.
    ///
    /// Note the encoding is `X^x Z^z`, **not** the Hermitian Pauli: `Pauli::Y`
    /// sets both bits and leaves the phase at zero, so the row denotes
    /// `X Z = −i Y`. That is the replaced `Row::set`'s convention verbatim, and
    /// every phase the decomposition returns is stated relative to it.
    #[inline]
    pub(crate) fn set(&mut self, i: usize, pauli: Pauli) {
        let (x, z) = match pauli {
            Pauli::I => (false, false),
            Pauli::X => (true, false),
            Pauli::Y => (true, true),
            Pauli::Z => (false, true),
        };
        TableauData::set_bit(&mut self.x, i, x);
        TableauData::set_bit(&mut self.z, i, z);
    }

    /// `phase += delta (mod 4)`.
    #[inline]
    pub(crate) fn add_phase(&mut self, delta: u8) {
        self.phase = (self.phase + delta) % 4;
    }

    /// Multiply the generator at `major` of `data` into `self` using the
    /// Aaronson–Gottesman `g`-rule.
    ///
    /// `data` must be in [`Orientation::RowMajor`], where a generator's bits are
    /// contiguous; `major` is the generator's index *within* `half`. This is the
    /// replaced `Row::mul_assign` — same predicates, same `ℤ/4` accumulation —
    /// reading the multiplicand straight out of the arena instead of from a
    /// copied row.
    #[inline]
    pub(crate) fn mul_generator(&mut self, data: &TableauData, half: Half, major: usize) {
        debug_assert_eq!(data.orientation(), Orientation::RowMajor);
        let src_x = data.major(half, Plane::X, major);
        let src_z = data.major(half, Plane::Z, major);
        let g = blocks::row_multiply(&mut self.x, &mut self.z, src_x, src_z);
        self.add_phase(g);
        self.add_phase(data.phase_of(half, major));
    }
}

// ─── Tableau ──────────────────────────────────────────────────────────────

/// A `2n`-generator stabilizer / destabilizer frame.
///
/// Generators `0..n` are the destabilizers, `n..2n` the stabilizers. The frame
/// is a genuine symplectic basis: the `2n` generators satisfy `ω(dᵢ, sⱼ) = δᵢⱼ`,
/// are linearly independent, start as such and stay such under every Clifford
/// generator — machine-checked in `lean/PPVM/Tableau/Frame.lean`
/// (`IsSymplecticFrame`, `frame_linearIndependent`, `isSymplecticFrame_identity`,
/// `isSymplecticFrame_hAct`/`sAct`/`cnotAct`/`czAct`).
///
/// # Storage
///
/// The bits live in [`TableauData`]: one aligned contiguous allocation holding
/// four square `n × n` X/Z quadrants, two bit-packed `ℤ/4` phase planes and the
/// per-qubit loss plane. The canonical orientation is **column-major** — the
/// generator dimension is contiguous for a fixed qubit — so a one-qubit Clifford
/// is two `n`-bit sweeps rather than a walk over all `2n` generators. See
/// [`crate::storage`] for why, and for the padding invariant equality and
/// hashing rest on.
///
/// Operations that need the opposite orientation (row multiplication, the
/// measurement projection, the decomposition's phase fold) take a
/// [`TransposedTableau`] guard, which transposes on construction and restores
/// the canonical orientation on drop.
///
/// Design: `tableau-data-structure.md`; `traits-2-configuration-and-hashing.md`
/// §"Pauli algebra traits" and §"Tableau indexability" — the tableau may itself
/// key a classical mixture, so it implements [`Indexable`] directly, owning the
/// lazy digest cache behind the contract (which fixes the digest *value*, not
/// the mechanism).
///
/// # Cache representation
///
/// The design sketches the lazy digest as a `OnceLock<u64>`; this uses the
/// sentinel [`AtomicU64`] that `ppvm-pauli-word-2` settled on. Same contract
/// (lazy, interior-mutable, `Send + Sync`, same finalized value), measurably
/// cheaper to *invalidate* — which every Clifford gate must do.
pub struct Tableau<H = fxhash::FxBuildHasher> {
    pub(crate) data: TableauData,
    /// How many [`TransposedTableau`] guards are currently held. Zero outside a
    /// guard, which is the only state any public method can be entered in.
    /// Excluded from equality, hashing and cloning — it is a borrow-lifetime
    /// fact, not part of the frame's identity.
    pub(crate) transpose_depth: usize,
    /// Lazy structural digest (Design: §"Lazy hashing and interior mutability").
    /// Holds [`HASH_UNCACHED`] until [`Indexable::key_hash`] first populates it;
    /// every structural mutation resets it through `&mut self`.
    pub(crate) hash_cache: AtomicU64,
    pub(crate) _hasher: PhantomData<fn() -> H>,
}

impl<H> Tableau<H> {
    /// Construct a fresh frame initialised to `|0…0⟩`.
    pub fn new(n_qubits: usize) -> Self {
        Self {
            data: TableauData::identity(n_qubits),
            transpose_depth: 0,
            hash_cache: AtomicU64::new(HASH_UNCACHED),
            _hasher: PhantomData,
        }
    }

    /// Restore the identity frame.
    pub fn reset_all(&mut self) {
        self.data.reset_to_identity();
        self.invalidate_hash();
    }

    /// Number of qubits.
    #[inline]
    pub fn n_qubits(&self) -> usize {
        self.data.n_qubits()
    }

    /// Clear the lazy digest after a structural mutation.
    ///
    /// Exclusive `&mut self` access — a plain store, no atomic RMW; the next
    /// `key_hash()` recomputes.
    #[inline]
    pub(crate) fn invalidate_hash(&mut self) {
        *self.hash_cache.get_mut() = HASH_UNCACHED;
    }

    /// The `(x-bits, z-bits, phase)` triple of every generator, destabilizers
    /// first, each bit vector qubit-indexed.
    ///
    /// The differential/snapshot view the tests reach for. Materializing a
    /// generator out of the column-major arena costs `O(n)` bit gathers, so this
    /// is a debugging and test surface, not a hot path — the frame's own
    /// algorithms read columns, or take the [`TransposedTableau`] guard.
    pub fn rows(&self) -> impl Iterator<Item = (Vec<u64>, Vec<u64>, u8)> + '_ {
        (0..2 * self.n_qubits()).map(|g| self.row(g))
    }

    /// The stabilizer generators' `(x-bits, z-bits, phase)` triples.
    pub fn stabilizer_rows(&self) -> impl Iterator<Item = (Vec<u64>, Vec<u64>, u8)> + '_ {
        let n = self.n_qubits();
        (n..2 * n).map(|g| self.row(g))
    }

    /// The destabilizer generators' `(x-bits, z-bits, phase)` triples.
    pub fn destabilizer_rows(&self) -> impl Iterator<Item = (Vec<u64>, Vec<u64>, u8)> + '_ {
        (0..self.n_qubits()).map(|g| self.row(g))
    }

    /// Materialize one generator as `(x-bits, z-bits, phase)`.
    pub fn row(&self, generator: usize) -> (Vec<u64>, Vec<u64>, u8) {
        let stride = self.data.stride();
        let mut x = vec![0u64; stride];
        let mut z = vec![0u64; stride];
        for q in 0..self.n_qubits() {
            TableauData::set_bit(&mut x, q, self.data.x_bit(generator, q));
            TableauData::set_bit(&mut z, q, self.data.z_bit(generator, q));
        }
        (x, z, self.data.phase(generator))
    }

    /// The Pauli at `(generator, qubit)`, destabilizers first.
    pub fn row_site(&self, generator: usize, qubit: usize) -> Pauli {
        match (
            self.data.x_bit(generator, qubit),
            self.data.z_bit(generator, qubit),
        ) {
            (false, false) => Pauli::I,
            (true, false) => Pauli::X,
            (false, true) => Pauli::Z,
            (true, true) => Pauli::Y,
        }
    }

    /// The `ℤ/4` phase of a generator, destabilizers first.
    #[inline]
    pub fn row_phase(&self, generator: usize) -> u8 {
        self.data.phase(generator)
    }

    /// Enter the row-major orientation, or nest inside an enclosing guard.
    ///
    /// The manual counterpart of [`TransposedTableau`], for a caller that cannot
    /// hold the guard's borrow — see `measure.rs`'s `RowGuard`. Pair every call
    /// with exactly one [`Self::exit_row_major`], from a `Drop`.
    #[inline]
    pub(crate) fn enter_row_major(&mut self) {
        if self.transpose_depth == 0 {
            self.data.transpose_quadrants();
        }
        self.transpose_depth += 1;
    }

    /// Leave the row-major orientation, restoring column-major at depth zero.
    #[inline]
    pub(crate) fn exit_row_major(&mut self) {
        self.transpose_depth -= 1;
        if self.transpose_depth == 0 {
            self.data.transpose_quadrants();
        }
    }

    /// The frame's X/Z bits as one contiguous byte range, for fingerprinting.
    #[inline]
    pub(crate) fn xz_bytes(&self) -> &[u8] {
        self.data.xz_bytes()
    }

    /// First stabilizer anticommuting with `Z_addr0`, if any.
    ///
    /// `ω(Z_q, sᵢ) = x_{sᵢ}[q]`, so this is the lowest set bit of the
    /// stabilizer-X quadrant's column `addr0` — one contiguous scan in the
    /// canonical orientation, where the replaced layout probed one bit in each
    /// of `n` separately addressed rows.
    pub(crate) fn find_z_anticommuting_stabilizer(&self, addr0: usize) -> Option<usize> {
        if self.data.orientation() == Orientation::ColumnMajor {
            let column = self.data.major(Half::Stab, Plane::X, addr0);
            return blocks::first_set(column).filter(|&i| i < self.n_qubits());
        }
        // Under a row guard the same predicate is a strided probe. `O(n)` scalar
        // reads, against the `O(n²/64)` the caller is under the guard for.
        (0..self.n_qubits())
            .find(|&i| TableauData::bit(self.data.major(Half::Stab, Plane::X, i), addr0))
    }

    /// The anticommutation mask of a single-site Pauli against one half of the
    /// frame: bit `i` is `ω(P at addr0, gᵢ)`.
    ///
    /// `ω(P, g) = x_P·z_g ⊕ z_P·x_g`, so with `P` supported on one site this is
    /// one of the two columns at `addr0`, or their `XOR` for `Y` — contiguous in
    /// the canonical orientation, where the replaced code probed the same site
    /// in each of `n` separately addressed rows.
    pub(crate) fn anticommutation_column(
        &self,
        half: Half,
        addr0: usize,
        pauli: Pauli,
    ) -> Vec<u64> {
        let mut out = vec![0u64; self.data.stride()];
        match pauli {
            Pauli::I => {}
            Pauli::X => self.data.gather_column(half, Plane::Z, addr0, &mut out),
            Pauli::Z => self.data.gather_column(half, Plane::X, addr0, &mut out),
            Pauli::Y => {
                let mut z = vec![0u64; self.data.stride()];
                self.data.gather_column(half, Plane::X, addr0, &mut out);
                self.data.gather_column(half, Plane::Z, addr0, &mut z);
                for (o, &zw) in out.iter_mut().zip(z.iter()) {
                    *o ^= zw;
                }
            }
        }
        out
    }

    /// The deterministic (case-b) measurement outcome for `Z_addr0`.
    ///
    /// `±Z_addr0` is a stabilizer; it is recovered as the product of the
    /// stabilizers whose destabilizer partner anticommutes with `Z_addr0`. The
    /// product must be real — the debug assert pins that, exactly as before.
    ///
    /// Takes the row guard: the product is a fold of whole generators, which is
    /// contiguous only in [`Orientation::RowMajor`].
    pub(crate) fn get_deterministic_outcome(&mut self, addr0: usize) -> bool {
        let n = self.n_qubits();
        // The selector is a *column* read, so it is taken before transposing.
        let selector = self.data.major(Half::Destab, Plane::X, addr0).to_vec();
        let stride = self.data.stride();

        let mut result = ScratchRow::zeroed(stride);
        let guard = TransposedTableau::new(self);
        for i in 0..n {
            if TableauData::bit(&selector, i) {
                result.mul_generator(guard.data(), Half::Stab, i);
            }
        }
        drop(guard);

        debug_assert!(
            result.phase == 0 || result.phase == 2,
            "Measurement result cannot be imaginary!"
        );
        result.phase >= 2
    }

    /// Project the frame onto the sampled case-a outcome (Aaronson–Gottesman
    /// row reduction against the pivot `q_idx`).
    ///
    /// This is the crate's **only** non-unitary frame mutation, so it is outside
    /// `IsSymplecticFrame.map` (which needs an `ω`-isometry). That it still maps
    /// a symplectic basis to a symplectic basis — the fact `canonicalize`'s
    /// no-op and every subsequent `compute_decomposition` rest on — is
    /// machine-checked as `isSymplecticFrame_projectFrame` in
    /// `lean/PPVM/Tableau/Frame.lean` (`projectFrame` is this sweep;
    /// `rowUpdate_eq_ite` is the `xbits[addr0]` conditional multiply).
    ///
    /// The two selectors are column reads and the eliminations are row
    /// multiplies, so the column reads are snapshotted first and the rest runs
    /// under the [`TransposedTableau`] guard — the same split Stim makes when it
    /// wraps `collapse_qubit_z` in `TableauTransposedRaii`.
    pub(crate) fn update_tableau_according_to_outcome(
        &mut self,
        addr0: usize,
        q_idx: usize,
        outcome: bool,
    ) {
        let n = self.n_qubits();
        self.invalidate_hash();
        let stride = self.data.stride();
        let mut stab_selector = vec![0u64; stride];
        let mut destab_selector = vec![0u64; stride];
        self.data
            .gather_column(Half::Stab, Plane::X, addr0, &mut stab_selector);
        self.data
            .gather_column(Half::Destab, Plane::X, addr0, &mut destab_selector);

        let mut guard = TransposedTableau::new(self);
        let data = guard.data_mut();

        // The pivot generator, copied once before the loop rewrites its
        // neighbours — the replaced code's `let g_q = stabilizers[q_idx];`.
        let mut pivot = ScratchRow::zeroed(stride);
        pivot
            .x
            .copy_from_slice(data.major(Half::Stab, Plane::X, q_idx));
        pivot
            .z
            .copy_from_slice(data.major(Half::Stab, Plane::Z, q_idx));
        pivot.phase = data.phase_of(Half::Stab, q_idx);

        for i in 0..n {
            if i == q_idx {
                continue;
            }
            if TableauData::bit(&stab_selector, i) {
                data.multiply_row_by(Half::Stab, i, &pivot.x, &pivot.z, pivot.phase);
            }
            if TableauData::bit(&destab_selector, i) {
                data.multiply_row_by(Half::Destab, i, &pivot.x, &pivot.z, pivot.phase);
            }
        }

        // The pivot becomes the new destabilizer; the new stabilizer is the
        // measured `±Z_addr0`.
        data.major_mut(Half::Destab, Plane::X, q_idx)
            .copy_from_slice(&pivot.x);
        data.major_mut(Half::Destab, Plane::Z, q_idx)
            .copy_from_slice(&pivot.z);
        data.set_phase_of(Half::Destab, q_idx, pivot.phase);

        data.major_mut(Half::Stab, Plane::X, q_idx).fill(0);
        let stab_z = data.major_mut(Half::Stab, Plane::Z, q_idx);
        stab_z.fill(0);
        TableauData::set_bit(stab_z, addr0, true);
        data.set_phase_of(Half::Stab, q_idx, if outcome { 2 } else { 0 });
    }
}

// ─── TransposedTableau ────────────────────────────────────────────────────

/// A frame temporarily held in [`Orientation::RowMajor`], where a generator's
/// bits are contiguous.
///
/// Construction transposes the four quadrants; [`Drop`] transposes them back, so
/// a public method always returns with the frame in its canonical column-major
/// orientation — including on unwind. The guard borrows `&mut Tableau` for its
/// whole lifetime, which is also what makes "hashing only ever observes the
/// canonical orientation" a *borrow-checker* fact rather than a convention: no
/// shared `&self` read can coexist with the guard.
///
/// The direct analogue of Stim's `TableauTransposedRaii`. Row multiplication and
/// elimination want the opposite grain from gates, and paying one transpose for
/// a whole batch of them beats fighting the layout on every row.
/// # Re-entrant
///
/// Guards nest. Only the outermost one transposes; inner ones just bump a depth
/// counter. That is what makes `measure_all` affordable: the transpose is
/// `O(n²/64)` and so is the elimination it enables, so paying it per
/// measurement would be a constant-factor disaster (measured at ~27× on a
/// 1889-qubit `measure_all`). One guard around the whole sweep amortizes it over
/// `n` measurements, exactly as Stim wraps a run of `collapse_qubit_z` calls in
/// a single `TableauTransposedRaii`. Everything the inner code reads is
/// orientation-aware, so it does not care which guard it is under.
pub(crate) struct TransposedTableau<'a, H> {
    tableau: &'a mut Tableau<H>,
}

impl<'a, H> TransposedTableau<'a, H> {
    /// Transpose into row-major and hold the frame there, unless an enclosing
    /// guard already has.
    #[inline]
    pub(crate) fn new(tableau: &'a mut Tableau<H>) -> Self {
        if tableau.transpose_depth == 0 {
            debug_assert_eq!(tableau.data.orientation(), Orientation::ColumnMajor);
            tableau.data.transpose_quadrants();
        }
        tableau.transpose_depth += 1;
        Self { tableau }
    }

    /// The row-major arena.
    #[inline]
    pub(crate) fn data(&self) -> &TableauData {
        &self.tableau.data
    }

    /// The row-major arena, mutably.
    #[inline]
    pub(crate) fn data_mut(&mut self) -> &mut TableauData {
        &mut self.tableau.data
    }
}

impl<H> Drop for TransposedTableau<'_, H> {
    #[inline]
    fn drop(&mut self) {
        self.tableau.transpose_depth -= 1;
        if self.tableau.transpose_depth == 0 {
            self.tableau.data.transpose_quadrants();
            debug_assert_eq!(
                self.tableau.data.orientation(),
                Orientation::ColumnMajor,
                "the outermost guard must restore the canonical orientation"
            );
        }
    }
}

/// Hand-written so the digest algorithm `H` — a private representation
/// parameter that is never a runtime value — does not have to be `Clone`.
impl<H> Clone for Tableau<H> {
    fn clone(&self) -> Self {
        Self {
            data: self.data.clone(),
            transpose_depth: 0,
            hash_cache: AtomicU64::new(self.hash_cache.load(Ordering::Relaxed)),
            _hasher: PhantomData,
        }
    }
}

/// Hand-written for the same reason as [`Clone`]; the digest cache is omitted
/// because it is not part of the frame's identity.
impl<H> Debug for Tableau<H> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Tableau")
            .field("n_qubits", &self.n_qubits())
            .field("data", &self.data)
            .finish()
    }
}

impl<H> PartialEq for Tableau<H> {
    /// Structural: width and bits. The digest cache is not part of the frame's
    /// identity. Compares the whole arena in bulk, which is sound because
    /// padding is held at zero (see [`crate::storage`]).
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl<H> Eq for Tableau<H> {}

impl<H: BuildHasher + Default> Hash for Tableau<H> {
    /// Per the [`Indexable`] contract: exactly `write_u64(self.key_hash())`.
    #[inline]
    fn hash<S: Hasher>(&self, state: &mut S) {
        state.write_u64(self.key_hash());
    }
}

impl<H: BuildHasher + Default> Indexable for Tableau<H> {
    /// The finalized structural digest of the frame.
    ///
    /// Design: §"Tableau indexability" — a tableau may key a classical mixture,
    /// so it is `Indexable` in its own right, owning the cache behind the
    /// contract (which fixes the *value*, not the mechanism). The raw digest of
    /// `H` is passed through a `splitmix64` finalizer so both the low bits (the
    /// hashbrown bucket) and the top 7 (the control tag) avalanche, which is the
    /// property the pass-through storage contract needs.
    ///
    /// The digest is taken over the arena's canonical ranges rather than
    /// generator by generator: zero padding makes the bulk read equivalent, and
    /// it is the whole reason the frame no longer has to materialize `2n` rows
    /// to hash itself.
    fn key_hash(&self) -> u64 {
        let cached = self.hash_cache.load(Ordering::Relaxed);
        if cached != HASH_UNCACHED {
            return cached;
        }
        let mut hasher = H::default().build_hasher();
        hasher.write_usize(self.n_qubits());
        self.data.hash(&mut hasher);
        let digest = finalize(hasher.finish());
        self.hash_cache.store(digest, Ordering::Relaxed);
        digest
    }
}

/// `splitmix64`'s finalizer — an avalanche-quality mix of a raw digest.
#[inline]
fn finalize(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

// ─── Amplitudes ───────────────────────────────────────────────────────────

/// The sparse amplitude vector `C[Bitstring]` of a [`GeneralizedTableau`].
///
/// Design: `traits-2-configuration-and-hashing.md` §"A third instantiation: the
/// generalized tableau" — the amplitudes are the same graded algebra as
/// `PauliSum`, over bitstring keys, with a `Vec` backend.
///
/// # Order is observable
///
/// The entry sequence is public behaviour: a branching gate leaves the vector in
/// ascending index order (the sort-merge output), a case-b measurement preserves
/// the pre-existing order (in-place `retain`), a case-a measurement leaves it
/// ascending, and `rotate_2` leaves it in its merge map's iteration order. Every
/// later fold (`normalize`, the overlaps) sums in vector order, so the order
/// also pins the float rounding. The tuple layout is `(value, index)` — the old
/// crate's, kept so a ported caller reads the same shape.
#[derive(Clone, Debug, PartialEq)]
pub struct Amplitudes<I> {
    entries: Vec<(Complex64, I)>,
}

impl<I: Bitstring> Default for Amplitudes<I> {
    fn default() -> Self {
        Self::new()
    }
}

impl<I: Bitstring> Amplitudes<I> {
    /// An empty amplitude vector.
    #[inline]
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// The `|0…0⟩` amplitude vector `{0 ↦ 1 + 0i}`.
    #[inline]
    pub fn unit() -> Self {
        Self {
            entries: vec![(Complex64::new(1.0, 0.0), I::zero())],
        }
    }

    /// The stored entries, in order, as `(value, index)` pairs.
    #[inline]
    pub fn entries(&self) -> &[(Complex64, I)] {
        &self.entries
    }

    /// Number of stored entries.
    #[inline]
    #[allow(clippy::len_without_is_empty)] // `is_empty` is right below.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether no entries are stored.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Borrow the stored entries without consuming the vector.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, (Complex64, I)> {
        self.entries.iter()
    }

    /// Insert without checking whether `index` is already present.
    ///
    /// A bare `push`: every caller on the hot path is provably collision-free
    /// (the XOR relabel is a bijection; the sort-merge emits each key once).
    #[inline]
    pub fn unsafe_insert(&mut self, index: I, value: Complex64) {
        self.entries.push((value, index));
    }

    /// Add `value` into the entry at `index`, creating it if absent
    /// (linear scan — a cold path only).
    pub fn add_or_insert(&mut self, index: I, value: Complex64) {
        for (v, i) in self.entries.iter_mut() {
            if *i == index {
                *v += value;
                return;
            }
        }
        self.entries.push((value, index));
    }

    /// The value at `index`, or zero if absent (linear scan — a cold path only).
    pub fn get(&self, index: &I) -> Complex64 {
        for (v, i) in self.entries.iter() {
            if i == index {
                return *v;
            }
        }
        Complex64::new(0.0, 0.0)
    }

    /// Multiply every entry's value by `factor`.
    pub fn mul_by(&mut self, factor: Complex64) {
        for (v, _) in self.entries.iter_mut() {
            *v *= factor;
        }
    }

    /// Multiply the value at `index` by `factor`. No-op if absent.
    pub fn mul_element_by(&mut self, index: I, factor: Complex64) {
        for (v, i) in self.entries.iter_mut() {
            if *i == index {
                *v *= factor;
                return;
            }
        }
    }

    /// Drop entries whose magnitude is at most `|cutoff|`.
    pub fn trim(&mut self, cutoff: Complex64) {
        let cutoff_sq = cutoff.norm_sqr();
        self.entries
            .retain(|(element, _)| element.norm_sqr() > cutoff_sq);
    }

    /// Drop entries failing the predicate.
    pub fn retain_entries(&mut self, f: impl FnMut(&(Complex64, I)) -> bool) {
        self.entries.retain(f);
    }

    /// Reserve capacity for at least `additional` more entries.
    #[inline]
    pub fn reserve(&mut self, additional: usize) {
        self.entries.reserve(additional);
    }

    /// L2-normalize in place.
    ///
    /// # Panics
    ///
    /// Panics with `"Zero norm encountered during normalization"` on a zero-norm
    /// vector. Returns early (leaving the entries byte-identical) when the norm
    /// is exactly `1`.
    pub fn normalize(&mut self) {
        // `re*re + im*im` directly; `abs()*abs()` would compute a `hypot` only to
        // square it again.
        let norm: f64 = self
            .entries
            .iter()
            .fold(0.0, |acc, (v, _)| acc + v.re * v.re + v.im * v.im);

        if norm == 0.0 {
            panic!("Zero norm encountered during normalization");
        }
        if norm == 1.0 {
            return;
        }
        let inv_norm_sqrt = 1.0 / norm.sqrt();
        let scale = Complex64::new(inv_norm_sqrt, 0.0);
        for (v, _) in self.entries.iter_mut() {
            *v *= scale;
        }
    }

    /// Take the entries out, leaving an empty vector behind.
    #[inline]
    pub(crate) fn take(&mut self) -> Vec<(Complex64, I)> {
        std::mem::take(&mut self.entries)
    }
}

impl<I: Bitstring> IntoIterator for Amplitudes<I> {
    type Item = (Complex64, I);
    type IntoIter = std::vec::IntoIter<(Complex64, I)>;

    fn into_iter(self) -> Self::IntoIter {
        self.entries.into_iter()
    }
}

/// L0 of the graded algebra (`ppvm-traits-2` §"The map is a graded algebra over
/// `C[K]`"): the amplitude vector is a finitely-supported `Bitstring ⇀ ℂ`.
impl<I: Bitstring> Support for Amplitudes<I> {
    type Key = I;
    type Coeff = Complex64;

    #[inline]
    fn len(&self) -> usize {
        self.entries.len()
    }

    #[inline]
    fn get(&self, key: &I) -> Option<Complex64> {
        self.entries.iter().find(|(_, i)| i == key).map(|(v, _)| *v)
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = (I, Complex64)> {
        self.entries.iter().map(|&(v, i)| (i, v))
    }
}

/// L2, the `ℂ`-module action.
impl<I: Bitstring> Scale for Amplitudes<I> {
    #[inline]
    fn scale(&mut self, s: &Complex64) {
        self.mul_by(*s);
    }
}

/// The one non-algebraic map operation (`ppvm-traits-2::Retain`). Dropping
/// supported terms breaks module exactness, which is why it sits outside the
/// graded algebra.
impl<I: Bitstring> ppvm_traits_2::Retain<I, Complex64> for Amplitudes<I> {
    #[inline]
    fn retain(&mut self, keep: impl Fn(&I, &Complex64) -> bool) {
        self.entries.retain(|(v, i)| keep(i, v));
    }
}

// ─── GeneralizedTableau ───────────────────────────────────────────────────

/// A [`Tableau`] extended with a sparse amplitude vector, so non-Clifford gates
/// are representable.
///
/// The state is `U|c⟩`: a Clifford **frame** `U` (the [`Tableau`]) plus an
/// amplitude vector over bitstrings. A Clifford gate updates the frame only and
/// leaves the amplitudes fixed; a non-Clifford gate branches the amplitudes and
/// leaves the frame fixed (`lean/PPVM/Tableau/Bitstring.lean`: the XOR relabel
/// `idx ^ stab_anticomm_bits` is a bijection). Measurement is the one operation
/// that couples them.
///
/// Design: `traits-2-configuration-and-hashing.md` §"A third instantiation: the
/// generalized tableau".
///
/// # Examples
///
/// ```
/// use ppvm_tableau_2::GeneralizedTableau;
/// use ppvm_traits_2::{Clifford, Measure};
/// use rand::SeedableRng;
///
/// let mut tab: GeneralizedTableau = GeneralizedTableau::new(2, 1e-12);
/// let mut rng = rand::rngs::SmallRng::seed_from_u64(0);
/// tab.h(0);
/// tab.cnot(0, 1);
/// assert_eq!(tab.measure(0, &mut rng), tab.measure(1, &mut rng));
/// ```
pub struct GeneralizedTableau<I = usize, H = fxhash::FxBuildHasher> {
    /// The underlying Clifford frame.
    pub tableau: Tableau<H>,
    /// The sparse amplitude vector indexed by bitstrings.
    pub coefficients: Amplitudes<I>,
    /// Per-qubit loss flags.
    pub is_lost: Vec<bool>,
    /// Coefficient-magnitude threshold below which branches are dropped.
    pub coefficient_threshold: f64,
    /// Ordered log of every measurement performed (mirrors stim's record).
    pub measurement_record: Vec<Option<bool>>,
    /// Measurement working buffers, kept allocated between calls.
    ///
    /// Private and non-observable: the cached
    /// [`odd_phase_mask`](crate::MeasureScratch::odd_phase_mask) is reset to
    /// `None` every time a measurement entry point picks this up, so the *only*
    /// thing carried across calls is the `Vec` capacities — old constructed a
    /// fresh `MeasureScratch` per entry point and got the same `None`. See
    /// [`Measure::measure`](ppvm_traits_2::Measure::measure) for why.
    ///
    /// Boxed and optional so a tableau that is never measured pays **nothing**:
    /// `new` writes a null pointer, `clone`/`fork` copy a null pointer, and the
    /// measurement entry points move one pointer in and out rather than the
    /// whole multi-`Vec` struct. The 4000-shot sampler workload builds a fresh
    /// tableau per shot, so an inline field would tax every construction.
    ///
    /// Trade-off: a tableau that *has* been measured retains buffers sized to
    /// the largest support it ever projected (peak footprint is unchanged —
    /// those buffers existed transiently anyway — but the steady-state footprint
    /// grows). [`Self::reset_all`] and [`Self::fork`] both drop them.
    pub(crate) scratch: Option<Box<crate::measure::MeasureScratch<I>>>,
}

/// Hand-written so the digest algorithm `H` need not be `Clone`; `fork` and
/// every bench's `iter_batched_ref` setup lean on this staying cheap (`2n` rows
/// plus the support).
impl<I: Clone, H> Clone for GeneralizedTableau<I, H> {
    fn clone(&self) -> Self {
        Self {
            tableau: self.tableau.clone(),
            coefficients: self.coefficients.clone(),
            is_lost: self.is_lost.clone(),
            coefficient_threshold: self.coefficient_threshold,
            measurement_record: self.measurement_record.clone(),
            // A *fresh* (i.e. absent) working set, not a copy: the buffers are
            // logically empty between calls, so cloning them would only pay for
            // allocations the clone has not earned — and `fork` is on every
            // bench's setup path and the sampler's hot path.
            scratch: None,
        }
    }
}

impl<I: Debug, H> Debug for GeneralizedTableau<I, H> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GeneralizedTableau")
            .field("tableau", &self.tableau)
            .field("coefficients", &self.coefficients)
            .field("is_lost", &self.is_lost)
            .field("coefficient_threshold", &self.coefficient_threshold)
            .field("measurement_record", &self.measurement_record)
            .finish()
    }
}

impl<I: Bitstring, H> GeneralizedTableau<I, H> {
    /// Construct a generalized tableau in the `|0…0⟩` state.
    ///
    /// Branches whose coefficient magnitude falls at or below
    /// `coefficient_threshold` are dropped **during gate application** — unlike
    /// `PauliSum`, where truncation is caller-driven, the tableau's gates
    /// auto-truncate and there is no `truncate()` entry point.
    pub fn new(n_qubits: usize, coefficient_threshold: f64) -> Self {
        Self {
            tableau: Tableau::new(n_qubits),
            coefficients: Amplitudes::unit(),
            is_lost: vec![false; n_qubits],
            coefficient_threshold,
            measurement_record: Vec::new(),
            scratch: None,
        }
    }

    /// Restore the `|0…0⟩` state: fresh frame rows, amplitudes `{0 ↦ 1}`, all
    /// loss flags cleared, empty measurement record.
    pub fn reset_all(&mut self) {
        self.tableau.reset_all();
        self.coefficients = Amplitudes::unit();
        for l in self.is_lost.iter_mut() {
            *l = false;
        }
        self.measurement_record.clear();
        self.scratch = None;
    }

    /// Clone the whole logical state, including the measurement record.
    ///
    /// Randomness is supplied by the caller and is therefore never duplicated
    /// by a fork.
    pub fn fork(&self) -> Self {
        self.clone()
    }

    /// Number of qubits.
    #[inline]
    pub fn n_qubits(&self) -> usize {
        self.tableau.n_qubits()
    }

    /// All measurement outcomes recorded so far, in order.
    #[inline]
    pub fn current_measurement_record(&self) -> &[Option<bool>] {
        &self.measurement_record
    }

    /// Append an externally defined measurement result to the record (stim
    /// `MPAD`).
    pub fn append_measurement_record(&mut self, result: Option<bool>) {
        self.measurement_record.push(result);
    }

    /// Replace the most recent measurement record entry.
    pub fn overwrite_last_measurement_record(&mut self, result: Option<bool>) {
        if let Some(last) = self.measurement_record.last_mut() {
            *last = result;
        }
    }

    /// Decompose a single-qubit Pauli into stabilizer/destabilizer products.
    ///
    /// Any Pauli can be written `P_addr0 = i^φ · ∏ dₖ^γₖ · ∏ sₗ^λₗ` with
    /// `γₖ = 1` iff `P` anticommutes with `sₖ` and `λₗ = 1` iff it anticommutes
    /// with `dₗ` (Lemma 5 of T. J. Yoder, 2012). `O(n²)`.
    ///
    /// Returns `(phase, stab_anticomm_bits, destab_anticomm_bits)`;
    /// `stab_anticomm_bits` is the index shift a branch applies (`β` in Eq. (4)
    /// of the SOFT paper). The decomposition rests on the measurement dichotomy
    /// machine-checked in `lean/PPVM/Tableau/Frame.lean` (`measurement_dichotomy`,
    /// `measure_deterministic_iff_xfree`).
    ///
    /// Correctness of the two bitmasks is `frame_coordinate_expansion` in the
    /// same file: for a symplectic frame, `v = Σᵢ ω(v, sᵢ)·dᵢ + Σᵢ ω(v, dᵢ)·sᵢ`,
    /// so the anticommutation tests accumulated here *are* `v`'s coordinates in
    /// the frame basis (`frame_surjective` supplies the spanning half). That is
    /// why multiplying those generators into `p_word` provably cancels it down
    /// to the identity and the returned `p_word.phase` is meaningful — the
    /// residual is never silently a non-identity Pauli.
    ///
    /// The `ℤ/4` phase on top of those bits is covered by
    /// `lean/PPVM/Tableau/BranchPhase.lean`. The `add_phase(8 - 2 * phase)`
    /// step is "multiply the generator in, divide its phase squared out":
    /// `destabAction_sq` proves `g² = (−1)^{phase(g)}` and
    /// `add_phase_eight_sub` that `8 − 2·ph = 2·ph` in `ℤ/4`. The fixed
    /// two-loop visit order (**all stabilizers first, then all destabilizers**)
    /// is a genuine convention rather than a free choice:
    /// `stab_destab_commute_sign` shows the opposite order shifts the returned
    /// phase by `2·⟨destab_anticomm, stab_anticomm⟩`, and `frameOp_eq_shiftOp`
    /// confirms this order is the one the downstream per-coefficient formula
    /// ([`compute_phase_with_mask_static`]) is stated for.
    ///
    /// # Why this takes `&mut self`
    ///
    /// The two anticommutation masks are *column* reads, taken in the canonical
    /// column-major orientation where they are contiguous. The `ℤ/4` residual,
    /// by contrast, is a fold of whole generators, contiguous only in
    /// [`Orientation::RowMajor`] — so the multiply loop runs under the
    /// [`TransposedTableau`] guard, and the guard needs unique access. Nothing
    /// logically mutates: the frame is byte-identical on return.
    pub fn compute_decomposition(&mut self, addr0: usize, pauli: Pauli) -> (u8, I, I) {
        debug_assert_ne!(pauli, Pauli::I);
        let n = self.n_qubits();
        let stride = self.tableau.data.stride();

        // `ω(P, g) = x_g[addr0]·z_P ⊕ z_g[addr0]·x_P`, so each half's
        // anticommutation mask is one contiguous column of the arena. The
        // replaced code probed the same site on all `2n` separately addressed
        // rows; the *values*, the visit order and the accumulated phase below
        // are unchanged.
        let destab_anticomm = self
            .tableau
            .anticommutation_column(Half::Destab, addr0, pauli);
        let stab_anticomm = self
            .tableau
            .anticommutation_column(Half::Stab, addr0, pauli);
        let destab_anticomm_bits = bits_to_index::<I>(&destab_anticomm, n);
        let stab_anticomm_bits = bits_to_index::<I>(&stab_anticomm, n);

        let mut p_word = ScratchRow::zeroed(stride);
        p_word.set(addr0, pauli);

        let guard = TransposedTableau::new(&mut self.tableau);
        let data = guard.data();
        for i in 0..n {
            if TableauData::bit(&destab_anticomm, i) {
                // The stabilizer is its own inverse up to its phase; rather than
                // inverting we multiply and divide out the phase squared.
                let phase = data.phase_of(Half::Stab, i);
                p_word.mul_generator(data, Half::Stab, i);
                p_word.add_phase(8 - 2 * phase);
            }
        }
        for i in 0..n {
            if TableauData::bit(&stab_anticomm, i) {
                let phase = data.phase_of(Half::Destab, i);
                p_word.mul_generator(data, Half::Destab, i);
                p_word.add_phase(8 - 2 * phase);
            }
        }
        drop(guard);

        (p_word.phase, stab_anticomm_bits, destab_anticomm_bits)
    }

    /// Multi-qubit generalization of [`Self::compute_decomposition`]: conjugate
    /// an arbitrary Pauli word through the frame.
    ///
    /// Calls [`Self::compute_decomposition`] per non-identity site and multiplies
    /// the single-qubit conjugates in canonical form `i^φ X^x Z^z`, picking up
    /// the `(−1)^{popcount(z_running & x_new)}` cross-phase from
    /// `Z^{z_a} X^{x_b} = (−1)^{z_a·x_b} X^{x_b} Z^{z_a}`.
    ///
    /// That the fold reproduces the *genuine* conjugate of the whole word — not
    /// merely a self-consistent convention — is machine-checked in
    /// `lean/PPVM/Pauli/Word.lean`: `phaseExpN_eq_canon` / `Canon.toG_mul`
    /// reconcile this canonical `X^x Z^z` form with the `g(x,z) = i^{x·z} X^x Z^z`
    /// normalization that `phaseExpN` — and hence the real `2ⁿ×2ⁿ` matrix product
    /// of `PauliMatrix.tensorPauli_mul` — is stated for, and
    /// `crossPhase_cocycle` / `Canon.foldl_eq_prod` show this left-fold against a
    /// *running* Z-mask equals the ordered product of the per-site conjugates.
    /// There is no single-qubit oracle for this: the cross term vanishes
    /// identically at weight 1.
    pub(crate) fn compute_decomposition_word<W>(&mut self, word: &W) -> (u8, I, I)
    where
        W: ppvm_traits_2::Word<Site = Pauli>,
    {
        let mut phase = 0u8;
        let mut stab_anticomm = I::zero();
        let mut destab_anticomm = I::zero();
        for q in 0..self.n_qubits() {
            let p_q = word.get(q);
            if p_q == Pauli::I {
                continue;
            }
            let (q_phase, q_stab, q_destab) = self.compute_decomposition(q, p_q);
            let cross = 2 * (symplectic_inner(destab_anticomm, q_stab) as u8 % 2);
            phase = (phase + q_phase + cross) % 4;
            stab_anticomm = stab_anticomm ^ q_stab;
            destab_anticomm = destab_anticomm ^ q_destab;
        }
        (phase, stab_anticomm, destab_anticomm)
    }

    /// Bitmask whose bit `i` is set iff destabilizer `i` has an odd (imaginary)
    /// phase.
    ///
    /// Computed **once** per gate/measurement and then folded into the
    /// per-coefficient [`compute_phase_with_mask_static`], which is what keeps
    /// the branching inner loop `O(n + m)` rather than `O(n·m)`.
    ///
    /// The parity bit of a generator's `ℤ/4` phase has its own bit plane
    /// ([`crate::storage`]), so the "odd phase" predicate over the whole
    /// destabilizer half is that plane read verbatim — no per-generator walk.
    pub fn odd_phase_destabilizer_mask(&self) -> I {
        let n = self.n_qubits();
        bits_to_index::<I>(self.tableau.data.phase_plane(Half::Destab, false), n)
    }

    /// The replaced per-generator walk, kept as a test oracle for the bit-plane
    /// read above.
    #[cfg(test)]
    pub(crate) fn odd_phase_destabilizer_mask_by_walk(&self) -> I {
        let mut mask = I::zero();
        let one = I::one();
        for i in 0..self.n_qubits() {
            if self.tableau.row_phase(i) % 2 != 0 {
                mask |= one << i;
            }
        }
        mask
    }

    /// The per-coefficient phase, walking the destabilizers directly. Only used
    /// as a test oracle for the mask-hoisted [`compute_phase_with_mask_static`].
    #[cfg(test)]
    pub(crate) fn compute_phase(
        &self,
        destab_anticomm_bits: I,
        basis_index: I,
        stab_anticomm_bits: I,
    ) -> u8 {
        let one = I::one();
        let zero = I::zero();
        let mut phase = (2 * symplectic_inner(destab_anticomm_bits, basis_index) as u8) % 4;
        let active = basis_index & stab_anticomm_bits;
        for i in 0..self.n_qubits() {
            if active & (one << i) == zero {
                continue;
            }
            if self.tableau.row_phase(i) % 2 != 0 {
                phase = (phase + 2) % 4;
            }
        }
        phase
    }

    /// Branch every amplitude by a single-qubit Pauli: `c ↦ f·c` on the original
    /// index and `i^φ·g·c` on `idx ^ stab_anticomm_bits`.
    ///
    /// This is the shared kernel behind `t`/`t_dag`/`rotate_1` (hence
    /// `rx`/`ry`/`rz`/`r`/`u3`).
    ///
    /// # Truncation
    ///
    /// The magnitude cutoff is applied **inline** while writing the new vector:
    /// an entry survives iff `|c|² > threshold²`, i.e. **strictly** greater and
    /// **absolute**. That boundary differs from `ppvm-pauli-sum-2`'s
    /// `CoefficientThreshold` (`magnitude() >= threshold`), a genuine mismatch
    /// machine-checked in `lean/PPVM/Algebra/Truncation.lean` (`cutoff_mismatch`);
    /// the tableau's rule is reproduced verbatim.
    ///
    /// # Algorithm
    ///
    /// Sort-merge coalesce: build the non-branch stream `nb`, the branch values
    /// `brv`, and a `u64`-packed `(key << 16) | build_pos` key stream in one
    /// pass, `sort_unstable` the packed `u64` array (half the data movement of
    /// `(I, u32)` elements), then 2-way merge straight into the destination with
    /// the cutoff applied inline — no intermediate output `Vec`, no map, no
    /// per-key hashing. The `nb_sorted` pre-check skips the sort entirely in the
    /// common case (the previous gate left the vector sorted). Falls back to a
    /// generic `(I, u32)` sort when there are more than `0xFFFF` coefficients or
    /// any branch key needs 47+ bits.
    pub(crate) fn branch_with_coefficients(
        &mut self,
        addr0: usize,
        pauli: Pauli,
        coefficient_factor: Complex64,
        branch_factor: Complex64,
    ) {
        if self.is_lost[addr0] {
            return;
        }

        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(addr0, pauli);

        let odd_phase_mask = self.odd_phase_destabilizer_mask();
        let old_coefficients = self.coefficients.take();
        let n_coefficients = old_coefficients.len();

        let cutoff_sq = self.coefficient_threshold * self.coefficient_threshold;

        #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
        if n_coefficients >= 16_384 {
            // Match the old large-support FxHashMap fold exactly: branch before
            // non-branch insertion fixes collision identity, floating-point
            // addition order, and final iteration order. Computing each pair
            // immediately avoids the old intermediate parallel `pairs` vector.
            let mut map: fxhash::FxHashMap<I, Complex64> =
                fxhash::FxHashMap::with_capacity_and_hasher(2 * n_coefficients, Default::default());
            for (coeff, idx) in old_coefficients {
                let branch_idx = idx ^ stab_anticomm_bits;
                let bpc = compute_phase_with_mask_static(
                    destab_anticomm_bits,
                    idx,
                    stab_anticomm_bits,
                    odd_phase_mask,
                );
                let branch_phase = (bpc + phase_decomp) % 4;
                let phase_factor = COMPLEX_PHASE_CONVERSION[branch_phase as usize];
                let branch_coeff = phase_factor * coeff * branch_factor;
                let nonbranch_coeff = coeff * coefficient_factor;
                *map.entry(branch_idx).or_insert(Complex64::new(0.0, 0.0)) += branch_coeff;
                *map.entry(idx).or_insert(Complex64::new(0.0, 0.0)) += nonbranch_coeff;
            }
            self.coefficients.reserve(map.len());
            for (idx, coeff) in map {
                if coeff.norm_sqr() > cutoff_sq {
                    self.coefficients.unsafe_insert(idx, coeff);
                }
            }
            return;
        }

        let mut nb: Vec<(I, Complex64)> = Vec::with_capacity(n_coefficients);
        let mut brv: Vec<Complex64> = Vec::with_capacity(n_coefficients);
        let mut packed: Vec<u64> = Vec::with_capacity(n_coefficients);
        let mut packable = n_coefficients <= 0xFFFF;
        let mut nb_sorted = true;
        let mut prev: Option<I> = None;
        for (pos, (coeff, idx)) in (0_u32..).zip(old_coefficients) {
            debug_assert!(
                !(coeff.re == 0.0 && coeff.im == 0.0),
                "Coefficient should not be zero"
            );
            let branch_index = idx ^ stab_anticomm_bits;
            let bpc = compute_phase_with_mask_static(
                destab_anticomm_bits,
                idx,
                stab_anticomm_bits,
                odd_phase_mask,
            );
            let branch_phase = (bpc + phase_decomp) % 4;
            let pf = COMPLEX_PHASE_CONVERSION[branch_phase as usize];
            brv.push(pf * coeff * branch_factor);
            match <u64 as num::NumCast>::from(branch_index) {
                Some(k) if k < (1u64 << 47) => packed.push((k << 16) | (pos as u64)),
                _ => {
                    packable = false;
                    packed.push(pos as u64);
                }
            }
            nb.push((idx, coeff * coefficient_factor));
            if let Some(p) = prev
                && idx < p
            {
                nb_sorted = false;
            }
            prev = Some(idx);
        }

        self.coefficients.reserve(nb.len() + brv.len());
        let mut i = 0;
        if packable {
            if !nb_sorted {
                nb.sort_unstable_by_key(|a| a.0);
            }
            packed.sort_unstable();
            // Decode the 47-bit key byte-by-byte from `I::from(u8)` — `NumCast`
            // panics for `bnum` types.
            let decode_key = |w: u64| -> I {
                let k = w >> 16; // k < 2^47; 6 bytes suffice
                let mut v = I::zero();
                for b in 0..6usize {
                    let byte = ((k >> (b * 8)) & 0xFF) as u8;
                    v |= <I as From<u8>>::from(byte) << (b * 8);
                }
                v
            };
            let mut j = 0;
            while i < nb.len() && j < packed.len() {
                let bp = (packed[j] & 0xFFFF) as usize;
                let bk = decode_key(packed[j]);
                match nb[i].0.cmp(&bk) {
                    std::cmp::Ordering::Less => {
                        if nb[i].1.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(nb[i].0, nb[i].1);
                        }
                        i += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        let v = brv[bp];
                        if v.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(bk, v);
                        }
                        j += 1;
                    }
                    std::cmp::Ordering::Equal => {
                        let mut sv = nb[i].1;
                        sv += brv[bp];
                        if sv.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(nb[i].0, sv);
                        }
                        i += 1;
                        j += 1;
                    }
                }
            }
            while j < packed.len() {
                let bp = (packed[j] & 0xFFFF) as usize;
                let bk = decode_key(packed[j]);
                let v = brv[bp];
                if v.norm_sqr() > cutoff_sq {
                    self.coefficients.unsafe_insert(bk, v);
                }
                j += 1;
            }
        } else {
            // Fallback for wide index types, large keys (≥ 2^47), or > 65535
            // coefficients. `nb` is still in build order here; reconstruct the
            // branch keys from it *before* sorting `nb`, so build-position `p`
            // still indexes `brv[p]`.
            let mut brk: Vec<(I, u32)> = (0_u32..)
                .zip(nb.iter())
                .map(|(p, e)| (e.0 ^ stab_anticomm_bits, p))
                .collect();
            if !nb_sorted {
                nb.sort_unstable_by_key(|a| a.0);
            }
            brk.sort_unstable_by_key(|a| a.0);
            let mut j = 0;
            while i < nb.len() && j < brk.len() {
                let (bk, bp) = brk[j];
                match nb[i].0.cmp(&bk) {
                    std::cmp::Ordering::Less => {
                        if nb[i].1.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(nb[i].0, nb[i].1);
                        }
                        i += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        let v = brv[bp as usize];
                        if v.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(bk, v);
                        }
                        j += 1;
                    }
                    std::cmp::Ordering::Equal => {
                        let mut sv = nb[i].1;
                        sv += brv[bp as usize];
                        if sv.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(nb[i].0, sv);
                        }
                        i += 1;
                        j += 1;
                    }
                }
            }
            while j < brk.len() {
                let (bk, bp) = brk[j];
                let v = brv[bp as usize];
                if v.norm_sqr() > cutoff_sq {
                    self.coefficients.unsafe_insert(bk, v);
                }
                j += 1;
            }
        }
        while i < nb.len() {
            if nb[i].1.norm_sqr() > cutoff_sq {
                self.coefficients.unsafe_insert(nb[i].0, nb[i].1);
            }
            i += 1;
        }
    }

    /// Relabel `coefficients` by applying the single-qubit Pauli `pauli` at
    /// `addr0` — the bijective half of the branch (`idx ↦ idx ^
    /// stab_anticomm_bits`).
    ///
    /// XOR by a fixed constant is a bijection, so distinct inputs always map to
    /// distinct outputs (`lean/PPVM/Tableau/Bitstring.lean`): unlike the T-gate
    /// split there are no collisions at all, so a per-index coalesce could never
    /// merge anything and the flat `Vec` relabel is exact. Only `rotate_2` calls
    /// this. Applies the same inline absolute cutoff as
    /// [`Self::branch_with_coefficients`].
    pub(crate) fn compute_coefficients_after_pauli_apply(
        &mut self,
        coefficients: &mut Amplitudes<I>,
        addr0: usize,
        pauli: Pauli,
    ) {
        if self.is_lost[addr0] {
            return;
        }

        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(addr0, pauli);

        let odd_phase_mask = self.odd_phase_destabilizer_mask();
        let old_coefficients = coefficients.take();
        let n_coefficients = old_coefficients.len();

        let mut new_coefficients: Vec<(I, Complex64)> = Vec::with_capacity(n_coefficients);
        for (coeff, idx) in old_coefficients {
            debug_assert!(
                !(coeff.re == 0.0 && coeff.im == 0.0),
                "Coefficient should not be zero"
            );
            let branch_index = idx ^ stab_anticomm_bits;
            let bpc = compute_phase_with_mask_static(
                destab_anticomm_bits,
                idx,
                stab_anticomm_bits,
                odd_phase_mask,
            );
            let branch_phase = (bpc + phase_decomp) % 4;
            let phase_factor = COMPLEX_PHASE_CONVERSION[branch_phase as usize];
            new_coefficients.push((branch_index, phase_factor * coeff));
        }

        let cutoff_sq = self.coefficient_threshold * self.coefficient_threshold;
        coefficients.reserve(new_coefficients.len());
        for (idx, coeff) in new_coefficients {
            if coeff.norm_sqr() > cutoff_sq {
                coefficients.unsafe_insert(idx, coeff);
            }
        }
    }
}

// ─── Fused CZ blocks ──────────────────────────────────────────────────────

/// The word width the `cz_block` family segments on.
///
/// In the replaced layout this was `size_of::<A::Store>() * 8` — the compile-time
/// storage word — and it decided both *where* a run was split and, for
/// overlapping runs, *what the answer was*. The frame is runtime-sized now, so
/// the constant is pinned at 64, which is what every shipped configuration
/// (`[usize; K]` storage on a 64-bit target) already used.
const CZ_BLOCK_WORD: usize = BITS_PER_WORD;

impl<H> Tableau<H> {
    /// Apply CZ to `count` pairs at a constant offset — `(base + i, base +
    /// offset + i)` — **simultaneously**, all reads taken before any write.
    ///
    /// The single phase parity replaces `count` sequential `ℤ/4` updates. That
    /// is sound **only** because the pairs have pairwise-disjoint supports:
    /// `lean/PPVM/Tableau/Batch.lean` proves `czSeq_phase` under exactly that
    /// hypothesis, and `czSeq_phase_needs_disjoint` exhibits two overlapping
    /// pairs on which the batched parity and the per-pair loop **disagree** (a
    /// pair's sign reads a `z`-bit an earlier pair already rewrote).
    ///
    /// # The overlapping case is reproduced, not fixed
    ///
    /// When `offset < count` a qubit is both a control and a target. The
    /// replaced word kernel combined the two z-deltas with `|`, not `^`
    /// (`((x >> offset) & mask_c) | ((x << offset) & mask_t)`), which is
    /// `G-060` in `docs/lean-gap.md`. The delta pass below reproduces that `|`
    /// exactly. Adjacent-pair brickwork is precisely this case, so "fixing" it
    /// here would silently change every such circuit's output while the ledger
    /// is still adjudicating it.
    pub fn cz_block_pairs(&mut self, base: usize, offset: usize, count: usize) {
        if count == 0 {
            return;
        }
        self.invalidate_hash();
        let stride = self.data.stride();
        let mut delta = vec![0u64; stride];

        for half in HALVES {
            // Phase pass: every predicate reads pre-update `z`, so it runs to
            // completion before the delta pass touches a single z-word.
            for k in 0..count {
                let (c, t) = (base + k, base + offset + k);
                if c == t {
                    continue;
                }
                let (xc, zc, xt, zt, ph) = self.data.gate2_mut(half, c, t);
                for i in 0..ph.len() {
                    ph[i] ^= xc[i] & xt[i] & (zc[i] ^ zt[i]);
                }
            }
            if offset == 0 {
                // Degenerate: control and target coincide, so `x[u] | x[u]`.
                for u in base..base + count {
                    delta.copy_from_slice(self.data.major(half, Plane::X, u));
                    let z = self.data.major_mut(half, Plane::Z, u);
                    for (zw, &dw) in z.iter_mut().zip(delta.iter()) {
                        *zw ^= dw;
                    }
                }
                continue;
            }

            // Delta pass: `z[u] ^= ctrl_delta(u) | tgt_delta(u)`, where a qubit
            // acting as a control picks up its partner's x-column and one acting
            // as a target picks up its control's. `x` is never written, so both
            // reads are of pre-update values regardless of visit order.
            for u in base..base + offset + count {
                let is_control = u < base + count;
                let is_target = u >= base + offset;
                if !is_control && !is_target {
                    continue;
                }
                delta.fill(0);
                if is_control {
                    for (dw, &xw) in
                        delta
                            .iter_mut()
                            .zip(self.data.major(half, Plane::X, u + offset))
                    {
                        *dw |= xw;
                    }
                }
                if is_target {
                    for (dw, &xw) in
                        delta
                            .iter_mut()
                            .zip(self.data.major(half, Plane::X, u - offset))
                    {
                        *dw |= xw;
                    }
                }
                let z = self.data.major_mut(half, Plane::Z, u);
                for (zw, &dw) in z.iter_mut().zip(delta.iter()) {
                    *zw ^= dw;
                }
            }
        }
    }

    /// Apply CZ to `count` pairs whose controls and targets live in *different*
    /// `CZ_BLOCK_WORD`-sized index words.
    ///
    /// Distinct words mean the control set and the target set are disjoint, so
    /// every qubit appears in exactly one pair and the simultaneous form and the
    /// sequential loop coincide — `czSeq_phase`'s disjointness hypothesis holds.
    #[inline]
    pub fn cz_block_pairs_cross_word(
        &mut self,
        word_c: usize,
        base_bit_c: usize,
        word_t: usize,
        base_bit_t: usize,
        count: usize,
    ) {
        if count == 0 {
            return;
        }
        debug_assert!(base_bit_c + count <= CZ_BLOCK_WORD);
        debug_assert!(base_bit_t + count <= CZ_BLOCK_WORD);
        debug_assert_ne!(word_c, word_t);
        for k in 0..count {
            let c = word_c * CZ_BLOCK_WORD + base_bit_c + k;
            let t = word_t * CZ_BLOCK_WORD + base_bit_t + k;
            ppvm_traits_2::Clifford::cz(self, c, t);
        }
    }
}

impl<I: Bitstring, H> GeneralizedTableau<I, H> {
    /// Loss-aware [`Tableau::cz_block_pairs`]: falls back to a per-pair `cz`
    /// loop (skipping lost pairs) when any qubit in the range is lost.
    pub fn cz_block_pairs(&mut self, base: usize, offset: usize, count: usize) {
        let any_lost =
            (0..count).any(|i| self.is_lost[base + i] || self.is_lost[base + offset + i]);
        if !any_lost {
            self.tableau.cz_block_pairs(base, offset, count);
        } else {
            for i in 0..count {
                let c = base + i;
                let t = base + offset + i;
                if !self.is_lost[c] && !self.is_lost[t] {
                    ppvm_traits_2::Clifford::cz(&mut self.tableau, c, t);
                }
            }
        }
    }

    /// Loss-aware [`Tableau::cz_block_pairs_cross_word`].
    pub fn cz_block_pairs_cross_word(
        &mut self,
        word_c: usize,
        base_bit_c: usize,
        word_t: usize,
        base_bit_t: usize,
        count: usize,
    ) {
        let any_lost = (0..count).any(|i| {
            let c = word_c * CZ_BLOCK_WORD + base_bit_c + i;
            let t = word_t * CZ_BLOCK_WORD + base_bit_t + i;
            self.is_lost[c] || self.is_lost[t]
        });
        if !any_lost {
            self.tableau
                .cz_block_pairs_cross_word(word_c, base_bit_c, word_t, base_bit_t, count);
        } else {
            for i in 0..count {
                let c = word_c * CZ_BLOCK_WORD + base_bit_c + i;
                let t = word_t * CZ_BLOCK_WORD + base_bit_t + i;
                if !self.is_lost[c] && !self.is_lost[t] {
                    ppvm_traits_2::Clifford::cz(&mut self.tableau, c, t);
                }
            }
        }
    }

    /// Apply CZ to `count` pairs `(control_base + i, target_base + i)`.
    ///
    /// The high-level entry point for a fused block of CZs: splits the run at
    /// `CZ_BLOCK_WORD` boundaries and dispatches each segment to
    /// [`Self::cz_block_pairs`] (same word) or
    /// [`Self::cz_block_pairs_cross_word`] (straddling two). CZ is symmetric, so
    /// the two bases may be passed in either order.
    ///
    /// The segmentation is kept even though the layout no longer has a
    /// "storage word": it is observable. Where a run is split decides which
    /// pairs land in the same simultaneous batch, and for an overlapping run
    /// that changes the answer (see [`Tableau::cz_block_pairs`]).
    pub fn cz_block(&mut self, control_base: usize, target_base: usize, count: usize) {
        if count == 0 {
            return;
        }
        // `cz_block_pairs` needs a non-negative offset; CZ is symmetric.
        let (lo, hi) = if control_base <= target_base {
            (control_base, target_base)
        } else {
            (target_base, control_base)
        };
        let mut i = 0;
        while i < count {
            let (c, t) = (lo + i, hi + i);
            let (wc, bc) = (c / CZ_BLOCK_WORD, c % CZ_BLOCK_WORD);
            let (wt, bt) = (t / CZ_BLOCK_WORD, t % CZ_BLOCK_WORD);
            // Longest run before either index crosses into the next word.
            let run = (CZ_BLOCK_WORD - bc).min(CZ_BLOCK_WORD - bt).min(count - i);
            if wc == wt {
                self.cz_block_pairs(c, t - c, run);
            } else {
                self.cz_block_pairs_cross_word(wc, bc, wt, bt, run);
            }
            i += run;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::{Clifford, RotationOne};

    type TestTableau = GeneralizedTableau<usize>;

    fn snapshot<H>(tab: &Tableau<H>) -> Vec<(Vec<u64>, Vec<u64>, u8)> {
        tab.rows().collect()
    }

    // ─── compute_phase vs the mask-hoisted static form ────────────────

    #[test]
    fn compute_phase_z_on_single_qubit_plus_state() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);

        let (phase_decomp, stab, destab) = tab.compute_decomposition(0, Pauli::Z);
        assert_eq!(phase_decomp + tab.compute_phase(destab, 0, stab), 0);
        assert_eq!(phase_decomp + tab.compute_phase(destab, 1, stab), 0);
    }

    #[test]
    fn compute_phase_z_on_y_stabilizer() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.s(0);

        let (phase_decomp, stab, destab) = tab.compute_decomposition(0, Pauli::Z);
        assert_eq!(phase_decomp + tab.compute_phase(destab, 0, stab), 0);
        assert_eq!(phase_decomp + tab.compute_phase(destab, 1, stab), 0);
    }

    #[test]
    fn compute_phase_z_on_hsh_state() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.s(0);
        tab.tableau.h(0);

        let (phase_decomp, stab, destab) = tab.compute_decomposition(0, Pauli::Z);
        assert_eq!((phase_decomp + tab.compute_phase(destab, 0, stab)) % 4, 1);
        assert_eq!((phase_decomp + tab.compute_phase(destab, 1, stab)) % 4, 3);
    }

    /// The hoisted mask must reproduce the per-coefficient destabilizer walk
    /// exactly — that equality is what licenses the `O(n + m)` inner loop.
    #[test]
    fn mask_hoisted_phase_matches_destabilizer_walk() {
        let mut tab: TestTableau = GeneralizedTableau::new(4, 1e-12);
        tab.tableau.h(0);
        tab.tableau.s(0);
        tab.tableau.cnot(0, 1);
        tab.tableau.h(2);
        tab.tableau.s(2);
        tab.tableau.cz(2, 3);

        let (_, stab, destab) = tab.compute_decomposition(1, Pauli::Z);
        let mask = tab.odd_phase_destabilizer_mask();
        for idx in 0usize..16 {
            assert_eq!(
                tab.compute_phase(destab, idx, stab),
                compute_phase_with_mask_static(destab, idx, stab, mask),
                "phase mismatch at index {idx}"
            );
        }
    }

    /// The odd-phase mask is now the low phase plane read verbatim rather than a
    /// per-generator walk. This pins the two against each other on a frame whose
    /// destabilizers carry every `ℤ/4` residue.
    #[test]
    fn odd_phase_mask_bit_plane_matches_the_generator_walk() {
        for n in [1usize, 4, 65, 70] {
            let mut tab: GeneralizedTableau<u128> = GeneralizedTableau::new(n, 1e-12);
            for q in 0..n {
                tab.tableau.h(q);
                if q % 2 == 0 {
                    tab.tableau.s(q);
                }
                if q % 3 == 0 && q + 1 < n {
                    tab.tableau.cnot(q, q + 1);
                }
            }
            assert_eq!(
                tab.odd_phase_destabilizer_mask(),
                tab.odd_phase_destabilizer_mask_by_walk(),
                "n = {n}"
            );
        }
    }

    // ─── cz_block family ──────────────────────────────────────────────

    #[test]
    fn cz_block_pairs_matches_individual() {
        let n = 8;
        let mut tab1: Tableau = Tableau::new(n);
        tab1.h(0);
        tab1.h(3);
        tab1.s(1);
        let mut tab2 = tab1.clone();

        for i in 0..4 {
            tab1.cz(i, 4 + i);
        }
        tab2.cz_block_pairs(0, 4, 4);

        assert_eq!(snapshot(&tab1), snapshot(&tab2));
    }

    #[test]
    fn cz_block_pairs_offset_17() {
        // The MSD ladder shape: (0,17)..(16,33), all in one u64 word.
        let n = 34;
        let mut tab1: Tableau<[u64; 2]> = Tableau::new(n);
        for i in 0..n {
            tab1.h(i);
        }
        let mut tab2 = tab1.clone();

        for i in 0..17 {
            tab1.cz(i, 17 + i);
        }
        tab2.cz_block_pairs(0, 17, 17);

        assert_eq!(snapshot(&tab1), snapshot(&tab2));
    }

    #[test]
    fn cz_block_zero_count_is_a_noop() {
        let mut tab: Tableau = Tableau::new(8);
        tab.h(0);
        let before = snapshot(&tab);
        tab.cz_block_pairs(0, 4, 0);
        assert_eq!(before, snapshot(&tab));
    }

    #[test]
    fn cz_block_matches_individual_across_word_boundary() {
        // (34,51)..(46,63) sits in word 0; (47,64)..(50,67) straddles.
        let n = 85;
        let mut tab1: GeneralizedTableau<u128> = GeneralizedTableau::new(n, 1e-12);
        for i in 0..n {
            tab1.h(i);
        }
        let mut tab2 = tab1.clone();

        let (control_base, target_base, count) = (34, 51, 17);
        for i in 0..count {
            tab1.cz(control_base + i, target_base + i);
        }
        tab2.cz_block(control_base, target_base, count);
        assert_eq!(snapshot(&tab1.tableau), snapshot(&tab2.tableau));

        // CZ is symmetric: reversed bases must agree.
        let mut tab3: GeneralizedTableau<u128> = GeneralizedTableau::new(n, 1e-12);
        for i in 0..n {
            tab3.h(i);
        }
        tab3.cz_block(target_base, control_base, count);
        assert_eq!(snapshot(&tab1.tableau), snapshot(&tab3.tableau));
    }

    #[test]
    fn cz_block_pairs_with_loss_falls_back_per_pair() {
        let n = 8;
        let mut tab1: TestTableau = GeneralizedTableau::new(n, 1e-12);
        for i in 0..n {
            tab1.h(i);
        }
        tab1.is_lost[2] = true;
        let mut tab2 = tab1.clone();

        for i in 0..4 {
            let (c, t) = (i, 4 + i);
            if !tab1.is_lost[c] && !tab1.is_lost[t] {
                tab1.tableau.cz(c, t);
            }
        }
        tab2.cz_block_pairs(0, 4, 4);

        assert_eq!(snapshot(&tab1.tableau), snapshot(&tab2.tableau));
    }

    // ─── reset_all ────────────────────────────────────────────────────

    #[test]
    fn reset_all_restores_fresh_state() {
        let mut tab: TestTableau = GeneralizedTableau::new(3, 1e-12);
        let fresh: TestTableau = GeneralizedTableau::new(3, 1e-12);

        tab.h(0);
        tab.cnot(0, 1);
        tab.ry(2, 0.7); // non-Clifford: branches the amplitude vector
        assert!(tab.coefficients.len() > 1);

        tab.reset_all();

        assert_eq!(snapshot(&tab.tableau), snapshot(&fresh.tableau));
        assert_eq!(tab.coefficients.entries(), fresh.coefficients.entries());
    }

    #[test]
    fn reset_all_clears_record_and_loss() {
        let mut tab: TestTableau = GeneralizedTableau::new(3, 1e-12);
        tab.append_measurement_record(Some(true));
        tab.append_measurement_record(None);
        tab.is_lost[0] = true;
        tab.is_lost[2] = true;

        tab.reset_all();

        assert!(tab.current_measurement_record().is_empty());
        assert!(tab.is_lost.iter().all(|&lost| !lost));
    }

    // ─── Indexable ────────────────────────────────────────────────────

    /// The digest is a function of the frame alone: equal frames agree, a gate
    /// changes it, and the cached value survives a clone.
    #[test]
    fn key_hash_is_structural_and_cache_transparent() {
        let mut a: Tableau = Tableau::new(4);
        let mut b: Tableau = Tableau::new(4);
        assert_eq!(a.key_hash(), b.key_hash());

        a.h(0);
        assert_ne!(a.key_hash(), b.key_hash());
        b.h(0);
        assert_eq!(a.key_hash(), b.key_hash());

        let cloned = a.clone();
        assert_eq!(cloned.key_hash(), a.key_hash());
        assert_eq!(a, cloned);
    }

    /// A digest recomputed after invalidation must equal a fresh tableau's.
    #[test]
    fn key_hash_invalidates_on_mutation() {
        let mut a: Tableau = Tableau::new(3);
        let _ = a.key_hash(); // populate the cache
        a.cnot(0, 1);
        let recomputed = a.key_hash();

        let mut b: Tableau = Tableau::new(3);
        b.cnot(0, 1);
        assert_eq!(recomputed, b.key_hash());
    }

    // ─── Amplitudes ───────────────────────────────────────────────────

    #[test]
    fn amplitudes_normalize_and_trim() {
        let mut v: Amplitudes<usize> = Amplitudes::new();
        v.unsafe_insert(0, Complex64::new(3.0, 0.0));
        v.unsafe_insert(1, Complex64::new(4.0, 0.0));
        v.normalize();
        assert!((v.get(&0).re - 0.6).abs() < 1e-12);
        assert!((v.get(&1).re - 0.8).abs() < 1e-12);

        v.trim(Complex64::new(0.7, 0.0));
        assert_eq!(v.len(), 1);
        assert_eq!(Support::len(&v), 1);
    }

    #[test]
    #[should_panic(expected = "Zero norm encountered during normalization")]
    fn amplitudes_normalize_panics_on_zero_norm() {
        let mut v: Amplitudes<usize> = Amplitudes::new();
        v.unsafe_insert(0, Complex64::new(0.0, 0.0));
        v.normalize();
    }

    /// `normalize` returns early at norm exactly `1`, leaving the entries
    /// byte-identical.
    #[test]
    fn amplitudes_normalize_is_a_noop_at_unit_norm() {
        let mut v: Amplitudes<usize> = Amplitudes::unit();
        let before = v.entries().to_vec();
        v.normalize();
        assert_eq!(before, v.entries());
    }

    #[test]
    fn amplitudes_scale_is_the_module_action() {
        let mut v: Amplitudes<usize> = Amplitudes::unit();
        Scale::scale(&mut v, &Complex64::new(0.0, 2.0));
        assert_eq!(v.get(&0), Complex64::new(0.0, 2.0));
    }
}
