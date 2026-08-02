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

use bitvec::array::BitArray;
use bitvec::view::BitView;
use num::complex::Complex64;
use num::{One, PrimInt, Zero};
use ppvm_pauli_word_2::PauliStorage;
use ppvm_traits_2::{Indexable, Pauli, Scale, Support};
use rand::SeedableRng;
use rand::rngs::SmallRng;

/// Backing storage for a tableau row's packed X/Z bit planes.
///
/// Reuses [`PauliStorage`] from `ppvm-pauli-word-2` (the same blob bound the
/// packed word uses — `word-data-structures.md` §"`PauliWord` packed
/// representation") and adds the `PrimInt` bound on the raw machine word, which
/// is what lets every Clifford gate hoist `index / bits` and touch **one** word
/// per plane per row instead of going through `bitvec`'s bounds-checked per-bit
/// addressing (the old crate's `<T::Storage as BitView>::Store: PrimInt` clause).
pub trait RowStorage: PauliStorage + BitView<Store: PrimInt> {}

impl<A> RowStorage for A where A: PauliStorage + BitView<Store: PrimInt> {}

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

// ─── Row ──────────────────────────────────────────────────────────────────

/// One tableau row: a phased Pauli word over packed X/Z bit planes with a
/// `ℤ/4` phase.
///
/// # Why not `Phased<PauliWord>`
///
/// A row *is* a phased Pauli word mathematically, but the shipped
/// `ppvm-pauli-word-2::PauliWord` (a) carries a lazy `AtomicU64` hash cache and
/// is therefore **not `Copy`**, and (b) keeps its packed planes `pub(crate)`, so
/// no downstream crate can reach the raw machine words. Both are fatal here: the
/// tableau copies rows on every measurement projection and `2n` times per
/// construction (`let g_q = stabilizers[q_idx];`), and every Clifford gate needs
/// one raw-word read/write per plane per row. The row therefore owns its own
/// hash-free, `Copy` packed representation; the digest that
/// [`Indexable`](ppvm_traits_2::Indexable) exposes is computed over the *whole
/// tableau*, never cached per row.
///
/// Design: `word-data-structures.md` §"Logical Pauli model" for the `(x, z)`
/// encoding; the `ℤ/4` phase and the Aaronson–Gottesman `g`-rule product are the
/// "extension part" of `traits-2-configuration-and-hashing.md` §"Pauli algebra
/// traits".
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Row<A: RowStorage> {
    /// X-bit plane (unused high bits are permanently `0`).
    pub(crate) xbits: BitArray<A>,
    /// Z-bit plane (unused high bits are permanently `0`).
    pub(crate) zbits: BitArray<A>,
    /// Phase in `ℤ/4` with the convention `0: +1, 1: +i, 2: −1, 3: −i`.
    pub(crate) phase: u8,
}

impl<A: RowStorage> Row<A> {
    /// The identity row `+I…I`.
    #[inline]
    pub(crate) fn new(n_qubits: usize) -> Self {
        debug_assert!(
            n_qubits <= 8 * std::mem::size_of::<A>(),
            "n_qubits {n_qubits} exceeds the {}-bit backing storage",
            8 * std::mem::size_of::<A>(),
        );
        Self {
            xbits: BitArray::ZERO,
            zbits: BitArray::ZERO,
            phase: 0,
        }
    }

    /// Write the Pauli at site `i`.
    #[inline]
    pub(crate) fn set(&mut self, i: usize, pauli: Pauli) {
        let (x, z) = match pauli {
            Pauli::I => (false, false),
            Pauli::X => (true, false),
            Pauli::Y => (true, true),
            Pauli::Z => (false, true),
        };
        self.xbits.set(i, x);
        self.zbits.set(i, z);
    }

    /// Read the Pauli at site `i`.
    #[inline]
    pub(crate) fn get(&self, i: usize) -> Pauli {
        match (self.xbits[i], self.zbits[i]) {
            (false, false) => Pauli::I,
            (true, false) => Pauli::X,
            (false, true) => Pauli::Z,
            (true, true) => Pauli::Y,
        }
    }

    /// Resolve site `i` to the `(word index, single-bit mask)` pair the
    /// `_masked` accessors take.
    ///
    /// Hoisting this out of a row loop is the same move every Clifford gate
    /// makes (module note in `clifford.rs`): the frame's `O(n)` scans —
    /// `compute_decomposition`, `find_z_anticommuting_stabilizer`,
    /// `get_deterministic_outcome`, `update_tableau_according_to_outcome` — all
    /// probe the **same** site on all `2n` rows, so the `i / bits`, `i % bits`
    /// and shift belong outside the loop, leaving one raw-word `AND` per row
    /// instead of `bitvec`'s bounds-checked per-bit addressing.
    #[inline]
    pub(crate) fn site_word(i: usize) -> (usize, <A as BitView>::Store) {
        let bits = std::mem::size_of::<<A as BitView>::Store>() * 8;
        let one = <A as BitView>::Store::one();
        (i / bits, one << (i % bits))
    }

    /// Resolve site `i` **and** a fixed probe Pauli `pauli = (x_bit, z_bit)`
    /// into the `(word index, x-probe, z-probe)` triple
    /// [`Self::anticommutes_at_probe`] takes.
    ///
    /// The two probes are the site mask gated by the *opposite* component of
    /// the probe Pauli, because `ω(P, Q) = x_P·z_Q ⊕ z_P·x_Q`. Folding the
    /// gating in here turns the per-row test into two `AND`s and an `XOR` on
    /// whole words with no per-row `bool` materialization.
    #[inline]
    pub(crate) fn site_probe(
        i: usize,
        pauli: (bool, bool),
    ) -> (usize, <A as BitView>::Store, <A as BitView>::Store) {
        let (wi, mask) = Self::site_word(i);
        let zero = <A as BitView>::Store::zero();
        (
            wi,
            if pauli.1 { mask } else { zero },
            if pauli.0 { mask } else { zero },
        )
    }

    /// Whether this row anticommutes with the probe Pauli pre-resolved by
    /// [`Self::site_probe`] — `ω(P, Q) = x_P·z_Q ⊕ z_P·x_Q`, the same value as
    /// [`PauliBits::anticommutes_at`](ppvm_traits_2::PauliBits::anticommutes_at)
    /// would give at that site.
    #[inline]
    pub(crate) fn anticommutes_at_probe(
        &self,
        wi: usize,
        x_probe: <A as BitView>::Store,
        z_probe: <A as BitView>::Store,
    ) -> bool {
        let xw = self.xbits.data.as_raw_slice()[wi];
        let zw = self.zbits.data.as_raw_slice()[wi];
        ((xw & x_probe) ^ (zw & z_probe)) != <A as BitView>::Store::zero()
    }

    /// The X bit at a site pre-resolved by [`Self::site_word`].
    #[inline]
    pub(crate) fn x_at_masked(&self, wi: usize, mask: <A as BitView>::Store) -> bool {
        (self.xbits.data.as_raw_slice()[wi] & mask) != <A as BitView>::Store::zero()
    }

    /// `phase += delta (mod 4)`.
    #[inline]
    pub(crate) fn add_phase(&mut self, delta: u8) {
        self.phase = (self.phase + delta) % 4;
    }

    /// Multiply `rhs` into `self` using the Aaronson–Gottesman `g`-rule.
    ///
    /// Ported word-for-word from `ppvm-pauli-word/src/phase/mul.rs`
    /// (`MulAssign for PhasedPauliWord`): the `sign`/`imag` masks are popcounted
    /// per machine word, so the phase is recovered without a per-site table
    /// lookup. This is the `row_multiply` primitive of
    /// [`StabilizerFrame`](ppvm_traits_2::StabilizerFrame).
    #[inline]
    pub(crate) fn mul_assign(&mut self, rhs: &Self) {
        let mut sign_count = 0u32;
        let mut imag_count = 0u32;
        let lhs_x = &mut self.xbits.data;
        let lhs_z = &mut self.zbits.data;
        let rhs_x = &rhs.xbits.data;
        let rhs_z = &rhs.zbits.data;
        for i in 0..lhs_x.as_raw_slice().len() {
            let a = lhs_x.as_raw_slice()[i];
            let b = lhs_z.as_raw_slice()[i];
            let c = rhs_x.as_raw_slice()[i];
            let d = rhs_z.as_raw_slice()[i];
            let sign = (a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d);
            let imag = (a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d);
            sign_count += sign.count_ones();
            imag_count += imag.count_ones();
            lhs_x.as_raw_mut_slice()[i] = a ^ c;
            lhs_z.as_raw_mut_slice()[i] = b ^ d;
        }
        self.add_phase(((2 * sign_count + imag_count) % 4) as u8);
        self.add_phase(rhs.phase);
    }
}

// ─── Tableau ──────────────────────────────────────────────────────────────

/// A `2n`-row stabilizer / destabilizer frame.
///
/// Rows `0..n` hold the destabilizers, rows `n..2n` the stabilizers. Each row is
/// a phased Pauli word tracking its X/Z bits and a `ℤ/4` phase. The frame is a
/// genuine symplectic basis: the `2n` generators satisfy `ω(dᵢ, sⱼ) = δᵢⱼ`, are
/// linearly independent, start as such and stay such under every Clifford
/// generator — machine-checked in `lean/PPVM/Tableau/Frame.lean`
/// (`IsSymplecticFrame`, `frame_linearIndependent`, `isSymplecticFrame_identity`,
/// `isSymplecticFrame_hAct`/`sAct`/`cnotAct`/`czAct`).
///
/// Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits" and
/// §"Tableau indexability" — the tableau may itself key a classical mixture, so
/// it implements [`Indexable`] directly, owning the lazy digest cache behind the
/// contract (which fixes the digest *value*, not the mechanism).
///
/// # Cache representation
///
/// The design sketches the lazy cache as an `OnceLock<u64>`; this uses the
/// sentinel [`AtomicU64`] that `ppvm-pauli-word-2` already settled on, for the
/// same measured reason. Every Clifford gate must invalidate, and
/// `OnceLock::take()` resets `Once`'s state word as well as the payload:
/// A/B-ing the two representations on the 85-qubit MSD workload (same build,
/// interleaved, 5 separate processes) moved the new/old ratio from ~1.04 to
/// ~1.01. The contract is unchanged — lazy, interior-mutable, `Send + Sync`, and
/// the same finalized digest.
///
/// # Examples
///
/// ```
/// use ppvm_tableau_2::Tableau;
/// use ppvm_traits_2::Clifford;
///
/// let mut tab: Tableau = Tableau::new(2);
/// tab.h(0);
/// tab.cnot(0, 1);
/// assert_eq!(tab.n_qubits(), 2);
/// assert_eq!(tab.stabilizer_rows().count(), 2);
/// ```
pub struct Tableau<A: RowStorage = u64, H = fxhash::FxBuildHasher> {
    pub(crate) n_qubits: usize,
    /// Destabilizers in `0..n`, stabilizers in `n..2n`.
    pub(crate) data: Vec<Row<A>>,
    pub(crate) rng: SmallRng,
    /// Lazy structural digest (Design: §"Lazy hashing and interior mutability").
    /// Holds [`HASH_UNCACHED`] until [`Indexable::key_hash`] first populates it;
    /// every structural mutation resets it through `&mut self`.
    pub(crate) hash_cache: AtomicU64,
    pub(crate) _hasher: PhantomData<fn() -> H>,
}

impl<A: RowStorage, H> Tableau<A, H> {
    fn new_data(n_qubits: usize) -> Vec<Row<A>> {
        let mut data: Vec<Row<A>> = Vec::with_capacity(2 * n_qubits);
        let pw_cache = Row::<A>::new(n_qubits);
        for i in 0..n_qubits {
            // destabilizer
            let mut pw = pw_cache;
            pw.set(i, Pauli::X);
            data.push(pw);
        }
        for i in 0..n_qubits {
            // stabilizer
            let mut pw = pw_cache;
            pw.set(i, Pauli::Z);
            data.push(pw);
        }
        data
    }

    /// Construct a fresh frame initialised to `|0…0⟩`, OS-seeded.
    pub fn new(n_qubits: usize) -> Self {
        Self {
            n_qubits,
            data: Self::new_data(n_qubits),
            rng: rand::make_rng(),
            hash_cache: AtomicU64::new(HASH_UNCACHED),
            _hasher: PhantomData,
        }
    }

    /// Same as [`Tableau::new`], but seed the RNG deterministically.
    pub fn new_with_seed(n_qubits: usize, seed: u64) -> Self {
        let mut t = Self::new(n_qubits);
        t.rng = SmallRng::seed_from_u64(seed);
        t
    }

    /// Restore the identity frame. Does **not** reseed the RNG.
    pub fn reset_all(&mut self) {
        self.data = Self::new_data(self.n_qubits);
        self.invalidate_hash();
    }

    /// Number of qubits.
    #[inline]
    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Clear the lazy digest after a structural mutation.
    ///
    /// Exclusive `&mut self` access — a plain store, no atomic RMW; the next
    /// `key_hash()` recomputes.
    #[inline]
    pub(crate) fn invalidate_hash(&mut self) {
        *self.hash_cache.get_mut() = HASH_UNCACHED;
    }

    /// The `(x-plane, z-plane, phase)` triple of every row, destabilizers first.
    ///
    /// The differential/snapshot view the old crate's tests reach for by
    /// touching `tab.data[i].word.xbits.data`; exposed as a read-only projection
    /// so the packed planes themselves stay private.
    pub fn rows(&self) -> impl Iterator<Item = (A, A, u8)> + '_ {
        self.data
            .iter()
            .map(|r| (r.xbits.data, r.zbits.data, r.phase))
    }

    /// The stabilizer rows' `(x-plane, z-plane, phase)` triples.
    pub fn stabilizer_rows(&self) -> impl Iterator<Item = (A, A, u8)> + '_ {
        self.data[self.n_qubits..]
            .iter()
            .map(|r| (r.xbits.data, r.zbits.data, r.phase))
    }

    /// The destabilizer rows' `(x-plane, z-plane, phase)` triples.
    pub fn destabilizer_rows(&self) -> impl Iterator<Item = (A, A, u8)> + '_ {
        self.data[..self.n_qubits]
            .iter()
            .map(|r| (r.xbits.data, r.zbits.data, r.phase))
    }

    /// The Pauli at `(row, qubit)`, destabilizers first.
    pub fn row_site(&self, row: usize, qubit: usize) -> Pauli {
        self.data[row].get(qubit)
    }

    #[inline]
    pub(crate) fn stabilizers(&self) -> &[Row<A>] {
        &self.data[self.n_qubits..]
    }

    #[inline]
    pub(crate) fn destabilizers(&self) -> &[Row<A>] {
        &self.data[..self.n_qubits]
    }

    /// First stabilizer anticommuting with `Z_addr0`, if any.
    pub(crate) fn find_z_anticommuting_stabilizer(&self, addr0: usize) -> Option<usize> {
        let (wi, mask) = Row::<A>::site_word(addr0);
        self.stabilizers()
            .iter()
            .position(|stab| stab.x_at_masked(wi, mask))
    }

    /// The deterministic (case-b) measurement outcome for `Z_addr0`.
    ///
    /// `±Z_addr0` is a stabilizer; it is recovered as the product of the
    /// stabilizers whose destabilizer partner anticommutes with `Z_addr0`. The
    /// product must be real — the debug assert pins that, exactly as old.
    pub(crate) fn get_deterministic_outcome(&self, addr0: usize) -> bool {
        let destabilizers = self.destabilizers();
        let stabilizers = self.stabilizers();
        let n = self.n_qubits;
        let mut result = Row::<A>::new(n);
        // Same site on every row, so the word index and mask are hoisted and the
        // stabilizer partner is zipped rather than re-indexed (both slices are
        // exactly `n` long, so the pairing is the `i`-th of each).
        let (wi, mask) = Row::<A>::site_word(addr0);
        for (destab, stab) in destabilizers.iter().zip(stabilizers.iter()) {
            if destab.x_at_masked(wi, mask) {
                result.mul_assign(stab);
            }
        }

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
    pub(crate) fn update_tableau_according_to_outcome(
        &mut self,
        addr0: usize,
        q_idx: usize,
        outcome: bool,
    ) {
        let n = self.n_qubits;
        self.invalidate_hash();
        let (destabilizers, stabilizers) = self.data.split_at_mut(n);

        // Copy g_q once before the loop (a register/memcpy copy — this is why
        // rows must stay `Copy`).
        let g_q = stabilizers[q_idx];

        // One site, all `2n` rows: hoist the word index and mask (module note in
        // `clifford.rs`) and walk the two halves in lockstep.
        let (wi, mask) = Row::<A>::site_word(addr0);
        for (i, (destab, stab)) in destabilizers
            .iter_mut()
            .zip(stabilizers.iter_mut())
            .enumerate()
        {
            if i == q_idx {
                continue;
            }
            if stab.x_at_masked(wi, mask) {
                stab.mul_assign(&g_q);
            }
            if destab.x_at_masked(wi, mask) {
                destab.mul_assign(&g_q);
            }
        }

        destabilizers[q_idx] = g_q;

        let stab_q = &mut stabilizers[q_idx];
        stab_q.xbits = BitArray::ZERO;
        stab_q.zbits = BitArray::ZERO;
        stab_q.zbits.set(addr0, true);
        stab_q.phase = if outcome { 2 } else { 0 };
    }
}

/// Hand-written so the digest algorithm `H` — a private representation
/// parameter that is never a runtime value — does not have to be `Clone`.
impl<A: RowStorage, H> Clone for Tableau<A, H> {
    fn clone(&self) -> Self {
        Self {
            n_qubits: self.n_qubits,
            data: self.data.clone(),
            rng: self.rng.clone(),
            hash_cache: AtomicU64::new(self.hash_cache.load(Ordering::Relaxed)),
            _hasher: PhantomData,
        }
    }
}

/// Hand-written for the same reason as [`Clone`]; the RNG state and the digest
/// cache are omitted (neither is part of the frame's identity).
impl<A: RowStorage, H> Debug for Tableau<A, H> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Tableau")
            .field("n_qubits", &self.n_qubits)
            .field("data", &self.data)
            .finish_non_exhaustive()
    }
}

impl<A: RowStorage, H> PartialEq for Tableau<A, H> {
    /// Structural: width and rows. The RNG and the digest cache are not part of
    /// the frame's identity.
    fn eq(&self, other: &Self) -> bool {
        self.n_qubits == other.n_qubits && self.data == other.data
    }
}

impl<A: RowStorage, H> Eq for Tableau<A, H> {}

impl<A: RowStorage, H: BuildHasher + Default> Hash for Tableau<A, H> {
    /// Per the [`Indexable`] contract: exactly `write_u64(self.key_hash())`.
    #[inline]
    fn hash<S: Hasher>(&self, state: &mut S) {
        state.write_u64(self.key_hash());
    }
}

impl<A: RowStorage, H: BuildHasher + Default> Indexable for Tableau<A, H> {
    /// The finalized structural digest of the frame.
    ///
    /// Design: §"Tableau indexability" — a tableau may key a classical mixture,
    /// so it is `Indexable` in its own right, owning the cache behind the
    /// contract (which fixes the *value*, not the mechanism). The raw digest of
    /// `H` is passed through a `splitmix64` finalizer so both the low bits (the
    /// hashbrown bucket) and the top 7 (the control tag) avalanche, which is the
    /// property the pass-through storage contract needs.
    fn key_hash(&self) -> u64 {
        let cached = self.hash_cache.load(Ordering::Relaxed);
        if cached != HASH_UNCACHED {
            return cached;
        }
        let mut hasher = H::default().build_hasher();
        hasher.write_usize(self.n_qubits);
        for row in &self.data {
            row.xbits.data.hash(&mut hasher);
            row.zbits.data.hash(&mut hasher);
            hasher.write_u8(row.phase);
        }
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
///
/// let mut tab: GeneralizedTableau = GeneralizedTableau::new_with_seed(2, 1e-12, 0);
/// tab.h(0);
/// tab.cnot(0, 1);
/// assert_eq!(tab.measure(0), tab.measure(1));
/// ```
pub struct GeneralizedTableau<A: RowStorage = u64, I = usize, H = fxhash::FxBuildHasher> {
    /// The underlying Clifford frame.
    pub tableau: Tableau<A, H>,
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
impl<A: RowStorage, I: Clone, H> Clone for GeneralizedTableau<A, I, H> {
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

impl<A: RowStorage, I: Debug, H> Debug for GeneralizedTableau<A, I, H> {
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

impl<A: RowStorage, I: Bitstring, H> GeneralizedTableau<A, I, H> {
    /// Construct a generalized tableau in the `|0…0⟩` state, OS-seeded.
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

    /// Same as [`GeneralizedTableau::new`], but seed the RNG deterministically.
    pub fn new_with_seed(n_qubits: usize, coefficient_threshold: f64, seed: u64) -> Self {
        let mut s = Self::new(n_qubits, coefficient_threshold);
        s.tableau.rng = SmallRng::seed_from_u64(seed);
        s
    }

    /// Restore the `|0…0⟩` state: fresh frame rows, amplitudes `{0 ↦ 1}`, all
    /// loss flags cleared, empty measurement record. Does **not** reseed the RNG.
    pub fn reset_all(&mut self) {
        self.tableau.reset_all();
        self.coefficients = Amplitudes::unit();
        for l in self.is_lost.iter_mut() {
            *l = false;
        }
        self.measurement_record.clear();
        self.scratch = None;
    }

    /// Clone the whole state (including the measurement record) and reseed the
    /// RNG, producing an independent trajectory. `None` seeds from OS entropy.
    pub fn fork(&self, seed: Option<u64>) -> Self {
        let mut cloned = self.clone();
        cloned.tableau.rng = match seed {
            Some(s) => SmallRng::seed_from_u64(s),
            None => rand::make_rng(),
        };
        cloned
    }

    /// Number of qubits.
    #[inline]
    pub fn n_qubits(&self) -> usize {
        self.tableau.n_qubits
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
    pub fn compute_decomposition(&self, addr0: usize, pauli: Pauli) -> (u8, I, I) {
        let n = self.n_qubits();

        let mut p_word = Row::<A>::new(n);
        p_word.set(addr0, pauli);

        let mut destab_anticomm_bits = I::zero();
        let mut stab_anticomm_bits = I::zero();

        debug_assert_ne!(pauli, Pauli::I);
        let pauli_bits = match pauli {
            Pauli::I => (false, false),
            Pauli::X => (true, false),
            Pauli::Y => (true, true),
            Pauli::Z => (false, true),
        };

        let stabilizers = self.tableau.stabilizers();
        let destabilizers = self.tableau.destabilizers();
        let one = I::one();

        // Both scans probe site `addr0` on all `2n` rows, so the word index and
        // the single-bit mask are hoisted once (the same move every Clifford
        // gate makes — module note in `clifford.rs`) and the `i`-th partner row
        // comes from a `zip` rather than a bounds-checked re-index. Purely a
        // codegen change: the visit order, the multiplication order and the
        // accumulated phase are exactly the old crate's.
        let (wi, x_probe, z_probe) = Row::<A>::site_probe(addr0, pauli_bits);

        for (i, (stab, destab)) in stabilizers.iter().zip(destabilizers.iter()).enumerate() {
            if !destab.anticommutes_at_probe(wi, x_probe, z_probe) {
                continue;
            }
            destab_anticomm_bits |= one << i;
            // The stabilizer is its own inverse up to its phase; rather than
            // inverting we multiply and divide out the phase squared.
            p_word.mul_assign(stab);
            p_word.add_phase(8 - 2 * stab.phase);
        }

        for (i, (stab, destab)) in stabilizers.iter().zip(destabilizers.iter()).enumerate() {
            if !stab.anticommutes_at_probe(wi, x_probe, z_probe) {
                continue;
            }
            stab_anticomm_bits |= one << i;
            p_word.mul_assign(destab);
            p_word.add_phase(8 - 2 * destab.phase);
        }

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
    pub(crate) fn compute_decomposition_word<W>(&self, word: &W) -> (u8, I, I)
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
    /// Computed **once** per gate/measurement (an `O(n)` scan) and then folded
    /// into the per-coefficient [`compute_phase_with_mask_static`], which is what
    /// keeps the branching inner loop `O(n + m)` rather than `O(n·m)`.
    pub fn odd_phase_destabilizer_mask(&self) -> I {
        let mut mask = I::zero();
        let one = I::one();
        for (i, destab) in self.tableau.destabilizers().iter().enumerate() {
            if destab.phase % 2 != 0 {
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
        for (i, destab) in self.tableau.destabilizers().iter().enumerate() {
            if active & (one << i) == zero {
                continue;
            }
            if destab.phase % 2 != 0 {
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
        &self,
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

impl<A: RowStorage, H> Tableau<A, H> {
    /// Apply CZ to `count` pairs at a constant offset — `(base + i, base +
    /// offset + i)` — as a **single** shift+XOR word operation per row.
    ///
    /// All pairs must live in the same storage word. Replaces `count`
    /// full `2n`-row sweeps with one.
    ///
    /// The single `count_ones() & 1` phase parity replaces `count` sequential
    /// `ℤ/4` updates. That is sound **only** because the pairs have
    /// pairwise-disjoint supports: `lean/PPVM/Tableau/Batch.lean` proves
    /// `czSeq_phase` under exactly that hypothesis, and
    /// `czSeq_phase_needs_disjoint` exhibits two overlapping pairs on which the
    /// batched parity and the per-pair loop **disagree** (a pair's sign reads a
    /// `z`-bit an earlier pair already rewrote).
    ///
    /// # Panics
    ///
    /// Debug-asserts that all bits are within the same word.
    #[inline]
    pub fn cz_block_pairs(&mut self, base: usize, offset: usize, count: usize) {
        if count == 0 {
            return;
        }
        self.invalidate_hash();
        let bits_per_word = std::mem::size_of::<<A as BitView>::Store>() * 8;
        let base_bit = base % bits_per_word;
        let word_idx = base / bits_per_word;

        debug_assert_eq!(
            (base + offset + count - 1) / bits_per_word,
            word_idx,
            "All CZ pairs must be in the same word"
        );

        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let count_mask = if count >= bits_per_word {
            !zero
        } else {
            (one << count) - one
        };
        let mask_c = count_mask << base_bit;
        let mask_t = count_mask << (base_bit + offset);

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let x = xp[word_idx];
            let z = zp[word_idx];

            // Phase must use the original z, before the update.
            let xc = (x >> base_bit) & count_mask;
            let xt = (x >> (base_bit + offset)) & count_mask;
            let zc = (z >> base_bit) & count_mask;
            let zt = (z >> (base_bit + offset)) & count_mask;
            let phase_bits = xc & xt & (zc ^ zt);
            pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;

            // z[c] ^= x[t], z[t] ^= x[c]
            let z_delta = ((x >> offset) & mask_c) | ((x << offset) & mask_t);
            zp[word_idx] = z ^ z_delta;
        });
    }

    /// Apply CZ to `count` pairs whose controls and targets live in *different*
    /// storage words.
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
        self.invalidate_hash();
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let bits_per_word = std::mem::size_of::<<A as BitView>::Store>() * 8;

        debug_assert!(base_bit_c + count <= bits_per_word);
        debug_assert!(base_bit_t + count <= bits_per_word);
        debug_assert_ne!(word_c, word_t);

        let count_mask = if count >= bits_per_word {
            !zero
        } else {
            (one << count) - one
        };

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();

            let xc = (xp[word_c] >> base_bit_c) & count_mask;
            let xt = (xp[word_t] >> base_bit_t) & count_mask;
            let zc = (zp[word_c] >> base_bit_c) & count_mask;
            let zt = (zp[word_t] >> base_bit_t) & count_mask;

            let phase_bits = xc & xt & (zc ^ zt);
            pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;

            zp[word_c] = zp[word_c] ^ (xt << base_bit_c);
            zp[word_t] = zp[word_t] ^ (xc << base_bit_t);
        });
    }
}

impl<A: RowStorage, I: Bitstring, H> GeneralizedTableau<A, I, H> {
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
        let bits_per_word = std::mem::size_of::<<A as BitView>::Store>() * 8;
        let any_lost = (0..count).any(|i| {
            let c = word_c * bits_per_word + base_bit_c + i;
            let t = word_t * bits_per_word + base_bit_t + i;
            self.is_lost[c] || self.is_lost[t]
        });
        if !any_lost {
            self.tableau
                .cz_block_pairs_cross_word(word_c, base_bit_c, word_t, base_bit_t, count);
        } else {
            for i in 0..count {
                let c = word_c * bits_per_word + base_bit_c + i;
                let t = word_t * bits_per_word + base_bit_t + i;
                if !self.is_lost[c] && !self.is_lost[t] {
                    ppvm_traits_2::Clifford::cz(&mut self.tableau, c, t);
                }
            }
        }
    }

    /// Apply CZ to `count` pairs `(control_base + i, target_base + i)`.
    ///
    /// The high-level entry point for a fused block of CZs: splits the run at
    /// storage-word boundaries internally and dispatches each segment to
    /// [`Self::cz_block_pairs`] (same word) or
    /// [`Self::cz_block_pairs_cross_word`] (straddling two). CZ is symmetric, so
    /// the two bases may be passed in either order.
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
        let bits_per_word = std::mem::size_of::<<A as BitView>::Store>() * 8;
        let mut i = 0;
        while i < count {
            let (c, t) = (lo + i, hi + i);
            let (wc, bc) = (c / bits_per_word, c % bits_per_word);
            let (wt, bt) = (t / bits_per_word, t % bits_per_word);
            // Longest run before either index crosses into the next word.
            let run = (bits_per_word - bc).min(bits_per_word - bt).min(count - i);
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

    type TestTableau = GeneralizedTableau<u64, usize>;

    fn snapshot<A: RowStorage, H>(tab: &Tableau<A, H>) -> Vec<(A, A, u8)> {
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
        let mut tab1: GeneralizedTableau<[u64; 2], u128> = GeneralizedTableau::new(n, 1e-12);
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
        let mut tab3: GeneralizedTableau<[u64; 2], u128> = GeneralizedTableau::new(n, 1e-12);
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
        let mut b: Tableau = Tableau::new_with_seed(4, 99);
        assert_eq!(
            a.key_hash(),
            b.key_hash(),
            "the RNG is not part of identity"
        );

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
