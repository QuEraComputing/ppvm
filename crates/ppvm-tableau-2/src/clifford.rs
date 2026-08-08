// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Clifford row operations on the frame: the symplectic/phase primitives, the
//! **fused** [`Clifford`]/[`CliffordExtensions`] impls, the frame primitives,
//! and the batched mask sweeps.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits".
//!
//! # Why `Tableau` opts *out* of [`BlanketClifford`](ppvm_traits_2::BlanketClifford)
//!
//! The design lists `Tableau` among the blanket's implementers on the reasoning
//! that its phase primitive is "SIMD-wide", i.e. free. That does not hold for
//! this representation. The blanket runs the phase primitive and the column
//! primitives as *separate* steps, and for a tableau **each primitive is a full
//! sweep over all `2n` rows**: `h` would cost two sweeps instead of one and
//! `cnot`/`cz` three instead of one. Every Clifford touches all `2n` rows (the
//! MSD workload runs ~600 of them at `n = 85`), so taking the blanket is a 2–3×
//! regression on the single hottest loop in the crate.
//!
//! `Tableau` therefore does exactly what `Phased<W>` does for the same reason
//! (design §"The `BlanketClifford` marker and the fused phased override"): it
//! implements [`SymplecticColumns`] and [`PhaseTrack`] so the primitive contract
//! is honoured and auditable, stays out of the `BlanketClifford` marker, and
//! supplies its own read-once fused `impl Clifford` / `impl CliffordExtensions`
//! below. Those fused bodies are ported bit-for-bit from
//! `ppvm-tableau/src/gates/clifford.rs`, so the *values* are the blanket's; only
//! the redundant row sweeps are gone.
//!
//! Every gate hoists `bits`, `wi = index / bits`, `off = index % bits` and the
//! mask **outside** the row loop and then does one raw-word read (and at most one
//! write) per plane per row — never `bitvec`'s bounds-checked per-bit indexing
//! inside the loop.
//!
//! Each generator is an `Sp(2n, 2)` isometry (machine-checked per generator in
//! `lean/PPVM/Pauli/Symplectic.lean`: `hAct_isometry`, `sAct_isometry`,
//! `cnotAct_isometry`, `czAct_isometry`) with the sign action of
//! `lean/PPVM/Pauli/Conjugation.lean` (`conjH_sign`, `conjS_sign`, …), and it
//! maps a symplectic frame to a symplectic frame
//! (`lean/PPVM/Tableau/Frame.lean`, `IsSymplecticFrame.map`).

use bitvec::view::BitView;
use num::{One, PrimInt, Zero};
use ppvm_traits_2::{
    Clifford, CliffordBatch, CliffordExtensions, CliffordExtensionsBatch, PhaseTrack,
    StabilizerFrame, SymplecticColumns,
};
use smallvec::{SmallVec, smallvec};

use crate::data::{Bitstring, GeneralizedTableau, RowStorage, Tableau};

/// Per-word bitmask buffer used by the batched Clifford gates. Stack-allocates
/// for up to 8 storage words; spills to the heap beyond, so there is no hard
/// qubit cap.
type MaskBuf<A> = SmallVec<[<A as BitView>::Store; 8]>;

/// The width in bits of one raw storage word.
macro_rules! bits_per_word {
    ($a:ty) => {
        std::mem::size_of::<<$a as BitView>::Store>() * 8
    };
}

// ─── Sp / phase primitives ────────────────────────────────────────────────

/// The `Sp`-part of conjugation: whole-column bit-plane algebra over all `2n`
/// rows. `Tableau` deliberately does **not** opt into
/// [`BlanketClifford`](ppvm_traits_2::BlanketClifford) (see the module note); the
/// primitives are supplied so the design's audited contract holds and so a
/// mixture algebra keyed on `Tableau` can reach them.
impl<A: RowStorage, H> SymplecticColumns for Tableau<A, H> {
    #[inline]
    fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    #[inline]
    fn swap_xz(&mut self, q: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (q / bits, q % bits);
        let one = <A as BitView>::Store::one();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            xp[wi] = (xw & !mask) | (zw & mask);
            zp[wi] = (zw & !mask) | (xw & mask);
        });
    }

    #[inline]
    fn xor_z_from_x(&mut self, q: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (q / bits, q % bits);
        let one = <A as BitView>::Store::one();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xw = pw.xbits.data.as_raw_slice()[wi];
            let zp = pw.zbits.data.as_raw_mut_slice();
            zp[wi] = zp[wi] ^ (xw & mask);
        });
    }

    #[inline]
    fn xor_x_col(&mut self, ctrl: usize, tgt: usize) {
        let bits = bits_per_word!(A);
        let (wc, sc) = (ctrl / bits, ctrl % bits);
        let (wt, st) = (tgt / bits, tgt % bits);
        let one = <A as BitView>::Store::one();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let xa = (xp[wc] >> sc) & one;
            xp[wt] = xp[wt] ^ (xa << st);
        });
    }

    #[inline]
    fn xor_z_col(&mut self, tgt: usize, ctrl: usize) {
        let bits = bits_per_word!(A);
        let (wc, sc) = (ctrl / bits, ctrl % bits);
        let (wt, st) = (tgt / bits, tgt % bits);
        let one = <A as BitView>::Store::one();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let zp = pw.zbits.data.as_raw_mut_slice();
            let zb = (zp[wt] >> st) & one;
            zp[wc] = zp[wc] ^ (zb << sc);
        });
    }

    #[inline]
    fn cz_bits(&mut self, a: usize, b: usize) {
        let bits = bits_per_word!(A);
        let (wc, sc) = (a / bits, a % bits);
        let (wt, st) = (b / bits, b % bits);
        let one = <A as BitView>::Store::one();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_slice();
            let xa = (xp[wc] >> sc) & one;
            let xb = (xp[wt] >> st) & one;
            let zp = pw.zbits.data.as_raw_mut_slice();
            zp[wc] = zp[wc] ^ (xb << sc);
            zp[wt] = zp[wt] ^ (xa << st);
        });
    }
}

/// The extension part: the tableau's phase algebra is a `ℤ₂` sign (the `+2`
/// bit of the row's `ℤ₄` phase) plus the Aaronson–Gottesman `g` rule, which
/// lives behind [`StabilizerFrame::row_multiply`].
impl<A: RowStorage, H> PhaseTrack for Tableau<A, H> {
    #[inline]
    fn flip_phase_where_xz(&mut self, q: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (q / bits, q % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xw = pw.xbits.data.as_raw_slice()[wi];
            let zw = pw.zbits.data.as_raw_slice()[wi];
            pw.phase ^= ((((xw & zw) & mask) != zero) as u8) << 1;
        });
    }

    /// The backward `S` sign rule: flip where `x & z`, i.e. the same predicate
    /// as [`Self::flip_phase_where_xz`] (`S` and `H` share it; they differ only
    /// in the bit map).
    #[inline]
    fn s_phase(&mut self, q: usize) {
        self.flip_phase_where_xz(q);
    }

    #[inline]
    fn cnot_phase(&mut self, ctrl: usize, tgt: usize) {
        let bits = bits_per_word!(A);
        let (wc, sc) = (ctrl / bits, ctrl % bits);
        let (wt, st) = (tgt / bits, tgt % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_slice();
            let zp = pw.zbits.data.as_raw_slice();
            let xa = (xp[wc] >> sc) & one;
            let za = (zp[wc] >> sc) & one;
            let xb = (xp[wt] >> st) & one;
            let zb = (zp[wt] >> st) & one;
            let phase_flip = (xa & zb) & (xb ^ za ^ one);
            pw.phase ^= (((phase_flip & one) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn cz_phase(&mut self, a: usize, b: usize) {
        let bits = bits_per_word!(A);
        let (wc, sc) = (a / bits, a % bits);
        let (wt, st) = (b / bits, b % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_slice();
            let zp = pw.zbits.data.as_raw_slice();
            let xa = (xp[wc] >> sc) & one;
            let za = (zp[wc] >> sc) & one;
            let xb = (xp[wt] >> st) & one;
            let zb = (zp[wt] >> st) & one;
            let phase_flip = (xa & xb) & (za ^ zb);
            pw.phase ^= (((phase_flip & one) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn x_phase(&mut self, q: usize) {
        Clifford::x(self, q);
    }

    #[inline]
    fn y_phase(&mut self, q: usize) {
        Clifford::y(self, q);
    }

    #[inline]
    fn z_phase(&mut self, q: usize) {
        Clifford::z(self, q);
    }
}

/// The role-exclusive frame primitives.
///
/// Design: §"Pauli algebra traits" (`StabilizerFrame`). The frame is a genuine
/// symplectic basis and stays one under every Clifford generator
/// (`lean/PPVM/Tableau/Frame.lean`: `IsSymplecticFrame`,
/// `frame_linearIndependent`, `isSymplecticFrame_identity`,
/// `isSymplecticFrame_hAct`/`sAct`/`cnotAct`/`czAct`); the pivot search rests on
/// `measurement_dichotomy` / `measure_deterministic_iff_xfree`.
impl<A: RowStorage, H> StabilizerFrame for Tableau<A, H> {
    /// Index (within the stabilizer half) of the first generator anticommuting
    /// with `Z_qubit`, or `None` when the outcome is deterministic.
    #[inline]
    fn anticommuting_pivot(&self, qubit: usize) -> Option<usize> {
        self.find_z_anticommuting_stabilizer(qubit)
    }

    /// Multiply generator `src` into `dst` (row indices over all `2n` rows,
    /// destabilizers first) using the Aaronson–Gottesman `g`-rule.
    #[inline]
    fn row_multiply(&mut self, src: usize, dst: usize) {
        assert_ne!(src, dst, "row_multiply needs distinct rows");
        self.invalidate_hash();
        let s = self.data[src];
        self.data[dst].mul_assign(&s);
    }

    /// No-op: this representation never leaves canonical form.
    ///
    /// Every gate is an `Sp(2n, 2)` isometry on the frame
    /// (`isSymplecticFrame_*` in `lean/PPVM/Tableau/Frame.lean`) and the
    /// measurement projection restores the destabilizer/stabilizer pairing in
    /// place (`update_tableau_according_to_outcome`), so the `2n` rows are a
    /// symplectic basis after every public operation and there is nothing left
    /// to restore.
    ///
    /// The projection half is **not** covered by `IsSymplecticFrame.map` (that
    /// lemma needs an `ω`-isometry, which the projection is not — it overwrites
    /// two rows). It is machine-checked separately as
    /// `isSymplecticFrame_projectFrame` in `lean/PPVM/Tableau/Frame.lean`, whose
    /// `projectFrame` is exactly this crate's row sweep and whose
    /// `rowUpdate_eq_ite` is the `if row.xbits[addr0] { row.mul_assign(&g_q) }`
    /// conditional. Without that theorem this no-op would be an unproved claim
    /// carrying every downstream `compute_decomposition`.
    ///
    /// Kept as an explicit no-op rather than dropped so a caller written against
    /// the trait is portable to a backend that *does* defer canonicalization.
    #[inline]
    fn canonicalize(&mut self) {}
}

// ─── Fused Clifford ───────────────────────────────────────────────────────

/// Single source of truth for the per-gate Clifford phase/bit logic: every
/// method loops over the rows operating directly on the packed planes (raw
/// integer slices plus a hoisted word index and mask), never through `bitvec`'s
/// bounds-checked single-bit indexing inside the per-row loop.
///
/// Every caller — a bare `Tableau`, a `GeneralizedTableau` (which delegates here)
/// and the fused batch paths — runs through this one implementation, so there is
/// no parallel copy that can silently diverge.
impl<A: RowStorage, H> Clifford for Tableau<A, H> {
    #[inline]
    fn x(&mut self, qubit: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (qubit / bits, qubit % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let zw = pw.zbits.data.as_raw_slice()[wi];
            pw.phase ^= (((zw & mask) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn y(&mut self, qubit: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (qubit / bits, qubit % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xw = pw.xbits.data.as_raw_slice()[wi];
            let zw = pw.zbits.data.as_raw_slice()[wi];
            pw.phase ^= ((((xw ^ zw) & mask) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn z(&mut self, qubit: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (qubit / bits, qubit % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xw = pw.xbits.data.as_raw_slice()[wi];
            pw.phase ^= (((xw & mask) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn h(&mut self, qubit: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (qubit / bits, qubit % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            let xb = xw & mask;
            let zb = zw & mask;
            xp[wi] = (xw & !mask) | zb;
            zp[wi] = (zw & !mask) | xb;
            pw.phase ^= (((xb & zb) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn s(&mut self, qubit: usize) {
        // NOTE: S is the only Clifford where forward and backward propagation
        // differ (it is non-Hermitian); only the phase rule differs.
        let bits = bits_per_word!(A);
        let (wi, off) = (qubit / bits, qubit % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xw = pw.xbits.data.as_raw_slice()[wi];
            let zp = pw.zbits.data.as_raw_mut_slice();
            let zw = zp[wi];
            pw.phase ^= ((((xw & zw) & mask) != zero) as u8) << 1;
            zp[wi] = zw ^ (xw & mask);
        });
    }

    #[inline]
    fn cnot(&mut self, control: usize, target: usize) {
        let bits = bits_per_word!(A);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let (wc, sc) = (control / bits, control % bits);
        let (wt, st) = (target / bits, target % bits);
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let xa = (xp[wc] >> sc) & one;
            let za = (zp[wc] >> sc) & one;
            let xb = (xp[wt] >> st) & one;
            let zb = (zp[wt] >> st) & one;
            let phase_flip = (xa & zb) & (xb ^ za ^ one);
            pw.phase ^= (((phase_flip & one) != zero) as u8) << 1;
            zp[wc] = zp[wc] ^ (zb << sc);
            xp[wt] = xp[wt] ^ (xa << st);
        });
    }

    #[inline]
    fn cz(&mut self, qubit0: usize, qubit1: usize) {
        let bits = bits_per_word!(A);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let (wc, sc) = (qubit0 / bits, qubit0 % bits);
        let (wt, st) = (qubit1 / bits, qubit1 % bits);
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_slice();
            let xa = (xp[wc] >> sc) & one;
            let xb = (xp[wt] >> st) & one;
            let zp = pw.zbits.data.as_raw_mut_slice();
            let za = (zp[wc] >> sc) & one;
            let zb = (zp[wt] >> st) & one;
            let phase_flip = (xa & xb) & (za ^ zb);
            pw.phase ^= (((phase_flip & one) != zero) as u8) << 1;
            zp[wc] = zp[wc] ^ (xb << sc);
            zp[wt] = zp[wt] ^ (xa << st);
        });
    }
}

/// The extension gate set, fused per gate.
///
/// | Gate | `X` | `Y` | `Z` |
/// |:---:|:---:|:---:|:---:|
/// | `s` | `Y` | `-X` | `Z` |
/// | `s_dag` | `-Y` | `X` | `Z` |
/// | `sqrt_x` | `X` | `Z` | `-Y` |
/// | `sqrt_x_dag` | `X` | `-Z` | `Y` |
/// | `sqrt_y` | `-Z` | `Y` | `X` |
/// | `sqrt_y_dag` | `Z` | `Y` | `-X` |
impl<A: RowStorage, H> CliffordExtensions for Tableau<A, H> {
    #[inline]
    fn s_dag(&mut self, qubit: usize) {
        // The backwards-prop version of `S` is just `S†`: same bit mapping,
        // phase flips where `x & !z`.
        let bits = bits_per_word!(A);
        let (wi, off) = (qubit / bits, qubit % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xw = pw.xbits.data.as_raw_slice()[wi];
            let zp = pw.zbits.data.as_raw_mut_slice();
            let zw = zp[wi];
            pw.phase ^= ((((xw & !zw) & mask) != zero) as u8) << 1;
            zp[wi] = zw ^ (xw & mask);
        });
    }

    #[inline]
    fn sqrt_x(&mut self, qubit: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (qubit / bits, qubit % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let zw = pw.zbits.data.as_raw_slice()[wi];
            let xp = pw.xbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            pw.phase ^= ((((zw & !xw) & mask) != zero) as u8) << 1;
            xp[wi] = xw ^ (zw & mask);
        });
    }

    #[inline]
    fn sqrt_x_dag(&mut self, qubit: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (qubit / bits, qubit % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let zw = pw.zbits.data.as_raw_slice()[wi];
            let xp = pw.xbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            pw.phase ^= ((((xw & zw) & mask) != zero) as u8) << 1;
            xp[wi] = xw ^ (zw & mask);
        });
    }

    #[inline]
    fn sqrt_y(&mut self, qubit: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (qubit / bits, qubit % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            let xb = xw & mask;
            let zb = zw & mask;
            xp[wi] = (xw & !mask) | zb;
            zp[wi] = (zw & !mask) | xb;
            pw.phase ^= ((((xw & !zw) & mask) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn sqrt_y_dag(&mut self, qubit: usize) {
        let bits = bits_per_word!(A);
        let (wi, off) = (qubit / bits, qubit % bits);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mask = one << off;
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            let xb = xw & mask;
            let zb = zw & mask;
            xp[wi] = (xw & !mask) | zb;
            zp[wi] = (zw & !mask) | xb;
            pw.phase ^= ((((zw & !xw) & mask) != zero) as u8) << 1;
        });
    }

    // control: row, target: col
    // | CY  |  I  |  X  |  Y  |  Z  |
    // |:---:|:---:|:---:|:---:|:---:|
    // |  I  | II  | ZX  | IY  | ZZ  |
    // |  X  | XY  | -YZ | XI  | YX  |
    // |  Y  | YY  | XZ  | YI  | -XX |
    // |  Z  | ZI  | IX  | ZY  | IZ  |
    //
    // Bit transforms: xc'=xc, zc'=zc^xt^zt, xt'=xt^xc, zt'=zt^xc
    // Phase +2 when: xc & (xt ^ zt) & !(zc ^ zt)
    #[inline]
    fn cy(&mut self, control: usize, target: usize) {
        let bits = bits_per_word!(A);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let (wc, sc) = (control / bits, control % bits);
        let (wt, st) = (target / bits, target % bits);
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let xc = (xp[wc] >> sc) & one;
            let zc = (zp[wc] >> sc) & one;
            let xt = (xp[wt] >> st) & one;
            let zt = (zp[wt] >> st) & one;
            let phase_flip = (xc & (xt ^ zt)) & (zc ^ zt ^ one);
            pw.phase ^= (((phase_flip & one) != zero) as u8) << 1;
            zp[wc] = zp[wc] ^ ((xt ^ zt) << sc);
            xp[wt] = xp[wt] ^ (xc << st);
            zp[wt] = zp[wt] ^ (xc << st);
        });
    }
}

// ─── Batched mask sweeps ──────────────────────────────────────────────────

impl<A: RowStorage, H> Tableau<A, H> {
    /// Build per-word bitmasks from a list of qubit indices, returning
    /// `(masks, n_words)`.
    ///
    /// One mask per storage word, stack-allocated up to 8 words. Each row is
    /// then swept word-by-word, skipping zero-mask words and accumulating the
    /// phase with `count_ones() & 1` instead of per-qubit flips — which is what
    /// turns 16 individual `sqrt_y` calls (16 full `2n`-row sweeps) into one.
    ///
    /// Replacing the sequential `ℤ/4` phase updates by one parity is
    /// machine-checked in `lean/PPVM/Tableau/Batch.lean`: `two_mul_natCast`
    /// (a `+2` delta is `ℤ/2`-valued, so `2+2 ≡ 0` and XOR-accumulation is
    /// sound) and `seqApply_eq_batchApply` (the fused sweep equals the per-site
    /// loop — under the `Nodup` hypothesis that the indices are **distinct**).
    /// The per-gate sign predicates are pinned to the audited conjugation tables
    /// by `isSitewise_conjH` / `isSitewise_extSqrtY` / … in the same file.
    #[inline]
    fn build_masks(&self, indices: &[usize]) -> Option<(MaskBuf<A>, usize)> {
        if self.data.is_empty() || indices.is_empty() {
            return None;
        }
        let n_words = self.data[0].xbits.data.as_raw_slice().len();
        let bits_per_word = bits_per_word!(A);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mut masks: MaskBuf<A> = smallvec![zero; n_words];
        for &addr0 in indices {
            masks[addr0 / bits_per_word] =
                masks[addr0 / bits_per_word] | (one << (addr0 % bits_per_word));
        }
        Some((masks, n_words))
    }

    /// Build the Y mask while filtering lost sites, avoiding a separate loss
    /// probe followed by a second traversal of `indices`.
    #[inline]
    fn y_many_skipping(&mut self, indices: &[usize], is_lost: &[bool]) {
        if self.data.is_empty() || indices.is_empty() {
            return;
        }
        let n_words = self.data[0].xbits.data.as_raw_slice().len();
        let bits_per_word = bits_per_word!(A);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        let mut masks: MaskBuf<A> = smallvec![zero; n_words];
        let mut any = false;
        for &addr0 in indices {
            if !is_lost[addr0] {
                masks[addr0 / bits_per_word] =
                    masks[addr0 / bits_per_word] | (one << (addr0 % bits_per_word));
                any = true;
            }
        }
        if any {
            self.y_with_masks(&masks);
        }
    }

    #[inline]
    fn y_with_masks(&mut self, masks: &[<A as BitView>::Store]) {
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_slice();
            let zp = pw.zbits.data.as_raw_slice();
            let mut popcount = 0u32;
            for ((&xw, &zw), &mask) in xp.iter().zip(zp).zip(masks) {
                if mask != zero {
                    popcount += ((xw ^ zw) & mask).count_ones();
                }
            }
            pw.phase ^= ((popcount & 1) as u8) << 1;
        });
    }
}

impl<A: RowStorage, H> CliffordBatch for Tableau<A, H> {
    /// `X` is bit-preserving: the phase flips once per masked qubit with `z = 1`.
    #[inline]
    fn x_many(&mut self, indices: &[usize]) {
        let Some((masks, n_words)) = self.build_masks(indices) else {
            return;
        };
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let zp = pw.zbits.data.as_raw_slice();
            let mut popcount = 0u32;
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                popcount += (zp[wi] & mask).count_ones();
            }
            pw.phase ^= ((popcount & 1) as u8) << 1;
        });
    }

    /// `Y` is bit-preserving: the phase flips where `x ⊕ z = 1`.
    #[inline]
    fn y_many(&mut self, indices: &[usize]) {
        let Some((masks, _)) = self.build_masks(indices) else {
            return;
        };
        self.y_with_masks(&masks);
    }

    /// `Z` is bit-preserving: the phase flips where `x = 1`.
    #[inline]
    fn z_many(&mut self, indices: &[usize]) {
        let Some((masks, n_words)) = self.build_masks(indices) else {
            return;
        };
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_slice();
            let mut popcount = 0u32;
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                popcount += (xp[wi] & mask).count_ones();
            }
            pw.phase ^= ((popcount & 1) as u8) << 1;
        });
    }

    /// Forward `S`: the phase flips where `x & z = 1`, then `z ^= x` on the mask.
    #[inline]
    fn s_many(&mut self, indices: &[usize]) {
        let Some((masks, n_words)) = self.build_masks(indices) else {
            return;
        };
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let mut popcount = 0u32;
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let xw = xp[wi];
                let zw = zp[wi];
                popcount += ((xw & zw) & mask).count_ones();
                zp[wi] = zw ^ (xw & mask);
            }
            pw.phase ^= ((popcount & 1) as u8) << 1;
        });
    }

    /// Apply `H` to many qubits with one combined mask sweep per row. `H` swaps
    /// the x/z bits (as `√Y` does) with a different phase: `+2` where `x & z`.
    #[inline]
    fn h_many(&mut self, indices: &[usize]) {
        let Some((masks, n_words)) = self.build_masks(indices) else {
            return;
        };
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let not_mask = !mask;
                let xw = xp[wi];
                let zw = zp[wi];
                let x_bits = xw & mask;
                let z_bits = zw & mask;
                xp[wi] = (xw & not_mask) | z_bits;
                zp[wi] = (zw & not_mask) | x_bits;
                let phase_bits = x_bits & z_bits;
                pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;
            }
        });
    }

    /// Apply `CNOT` to many pairs on raw storage words. Pairs are applied
    /// sequentially per row, so semantics match the per-pair `cnot` loop exactly.
    fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
        let bits = bits_per_word!(A);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let mut phase_flips = zero;
            for &(control, target) in pairs {
                let (wc, sc) = (control / bits, control % bits);
                let (wt, st) = (target / bits, target % bits);
                let xa = (xp[wc] >> sc) & one;
                let za = (zp[wc] >> sc) & one;
                let xb = (xp[wt] >> st) & one;
                let zb = (zp[wt] >> st) & one;
                // +2 phase when x_a & z_b & !(x_b ^ z_a); 2+2 == 0 mod 4, so
                // XOR-accumulate.
                phase_flips = phase_flips ^ ((xa & zb) & (xb ^ za ^ one));
                zp[wc] = zp[wc] ^ (zb << sc);
                xp[wt] = xp[wt] ^ (xa << st);
            }
            pw.phase ^= (((phase_flips & one) != zero) as u8) << 1;
        });
    }

    /// Apply `CZ` to many pairs on raw storage words. `CZ` is symmetric and
    /// touches only z-bits; pairs are applied sequentially per row.
    fn cz_many(&mut self, pairs: &[(usize, usize)]) {
        let bits = bits_per_word!(A);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let mut phase_flips = zero;
            for &(control, target) in pairs {
                let (wc, sc) = (control / bits, control % bits);
                let (wt, st) = (target / bits, target % bits);
                let xa = (xp[wc] >> sc) & one;
                let za = (zp[wc] >> sc) & one;
                let xb = (xp[wt] >> st) & one;
                let zb = (zp[wt] >> st) & one;
                // +2 phase when x_a & x_b & (z_a ^ z_b)
                phase_flips = phase_flips ^ ((xa & xb) & (za ^ zb));
                zp[wc] = zp[wc] ^ (xb << sc);
                zp[wt] = zp[wt] ^ (xa << st);
            }
            pw.phase ^= (((phase_flips & one) != zero) as u8) << 1;
        });
    }
}

impl<A: RowStorage, H> CliffordExtensionsBatch for Tableau<A, H> {
    /// Backward `S` (i.e. `S†`): the `S` bit mapping, phase flips where `x & !z`.
    #[inline]
    fn s_dag_many(&mut self, indices: &[usize]) {
        let Some((masks, n_words)) = self.build_masks(indices) else {
            return;
        };
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let mut popcount = 0u32;
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let xw = xp[wi];
                let zw = zp[wi];
                popcount += ((xw & !zw) & mask).count_ones();
                zp[wi] = zw ^ (xw & mask);
            }
            pw.phase ^= ((popcount & 1) as u8) << 1;
        });
    }

    /// Apply `CY` to many pairs on raw storage words.
    fn cy_many(&mut self, pairs: &[(usize, usize)]) {
        let bits = bits_per_word!(A);
        let one = <A as BitView>::Store::one();
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            let mut phase_flips = zero;
            for &(control, target) in pairs {
                let (wc, sc) = (control / bits, control % bits);
                let (wt, st) = (target / bits, target % bits);
                let xc = (xp[wc] >> sc) & one;
                let zc = (zp[wc] >> sc) & one;
                let xt = (xp[wt] >> st) & one;
                let zt = (zp[wt] >> st) & one;
                // +2 phase when x_c & (x_t ^ z_t) & !(z_c ^ z_t)
                phase_flips = phase_flips ^ ((xc & (xt ^ zt)) & (zc ^ zt ^ one));
                zp[wc] = zp[wc] ^ ((xt ^ zt) << sc);
                xp[wt] = xp[wt] ^ (xc << st);
                zp[wt] = zp[wt] ^ (xc << st);
            }
            pw.phase ^= (((phase_flips & one) != zero) as u8) << 1;
        });
    }

    /// Apply `√Y` to many qubits with one combined mask sweep per row.
    #[inline]
    fn sqrt_y_many(&mut self, indices: &[usize]) {
        let Some((masks, n_words)) = self.build_masks(indices) else {
            return;
        };
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let not_mask = !mask;
                let xw = xp[wi];
                let zw = zp[wi];
                let x_bits = xw & mask;
                let z_bits = zw & mask;
                xp[wi] = (xw & not_mask) | z_bits;
                zp[wi] = (zw & not_mask) | x_bits;
                let phase_bits = x_bits & !z_bits;
                pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;
            }
        });
    }

    /// Apply `(√Y)†` to many qubits with one combined mask sweep per row.
    #[inline]
    fn sqrt_y_dag_many(&mut self, indices: &[usize]) {
        let Some((masks, n_words)) = self.build_masks(indices) else {
            return;
        };
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let not_mask = !mask;
                let xw = xp[wi];
                let zw = zp[wi];
                let x_bits = xw & mask;
                let z_bits = zw & mask;
                xp[wi] = (xw & not_mask) | z_bits;
                zp[wi] = (zw & not_mask) | x_bits;
                let phase_bits = z_bits & !x_bits;
                pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;
            }
        });
    }

    /// Apply `√X` to many qubits with one combined mask sweep per row.
    #[inline]
    fn sqrt_x_many(&mut self, indices: &[usize]) {
        let Some((masks, n_words)) = self.build_masks(indices) else {
            return;
        };
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let xw = xp[wi];
                let zw = zp[wi];
                let phase_bits = (zw & !xw) & mask;
                pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;
                xp[wi] = xw ^ (zw & mask);
            }
        });
    }

    /// Apply `(√X)†` to many qubits with one combined mask sweep per row.
    #[inline]
    fn sqrt_x_dag_many(&mut self, indices: &[usize]) {
        let Some((masks, n_words)) = self.build_masks(indices) else {
            return;
        };
        let zero = <A as BitView>::Store::zero();
        self.invalidate_hash();
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.xbits.data.as_raw_mut_slice();
            let zp = pw.zbits.data.as_raw_mut_slice();
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let xw = xp[wi];
                let zw = zp[wi];
                let phase_bits = (xw & zw) & mask;
                pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;
                xp[wi] = xw ^ (zw & mask);
            }
        });
    }
}

// ─── GeneralizedTableau: loss-aware forwarding ────────────────────────────

// Single-qubit gate: skip a lost qubit, otherwise delegate to the frame's
// canonical (word-level) method.
macro_rules! forward_single {
    ($name:ident) => {
        #[inline]
        fn $name(&mut self, index: usize) {
            if self.is_lost[index] {
                return;
            }
            self.tableau.$name(index);
        }
    };
}

// Two-qubit gate: skip the pair when either endpoint is lost.
macro_rules! forward_pair {
    ($name:ident) => {
        // Keep one compact monomorphized copy of the packed row kernel. Inlining
        // it through the generalized loss guard duplicates the full loop into
        // every gate/alias caller and makes LLVM spill the four addressed bits
        // at multi-word widths.
        #[inline(never)]
        fn $name(&mut self, control: usize, target: usize) {
            if self.is_lost[control] || self.is_lost[target] {
                return;
            }
            self.tableau.$name(control, target);
        }
    };
}

impl<A: RowStorage, I: Bitstring, H> Clifford for GeneralizedTableau<A, I, H> {
    forward_single!(x);
    forward_single!(y);
    forward_single!(z);
    forward_single!(h);
    forward_single!(s);
    forward_pair!(cnot);
    forward_pair!(cz);
}

impl<A: RowStorage, I: Bitstring, H> CliffordExtensions for GeneralizedTableau<A, I, H> {
    forward_single!(s_dag);
    forward_single!(sqrt_x);
    forward_single!(sqrt_x_dag);
    forward_single!(sqrt_y);
    forward_single!(sqrt_y_dag);
    forward_pair!(cy);
}

impl<A: RowStorage, I: Bitstring, H> GeneralizedTableau<A, I, H> {
    /// Fast path: is any qubit in the slice lost?
    #[inline(always)]
    fn any_lost_single(&self, indices: &[usize]) -> bool {
        indices.iter().any(|&i| self.is_lost[i])
    }

    /// Fast path: does any pair have a lost qubit?
    #[inline(always)]
    fn any_lost_pair(&self, pairs: &[(usize, usize)]) -> bool {
        pairs
            .iter()
            .any(|&(c, t)| self.is_lost[c] || self.is_lost[t])
    }
}

// The batched forms allocate **nothing** in the common case: scan first, then
// forward the untouched slice straight through; only an actual loss builds a
// filtered `Vec`. The surviving indices still get the gate.
macro_rules! forward_batch_single {
    ($name:ident) => {
        #[inline(always)]
        fn $name(&mut self, indices: &[usize]) {
            if !self.any_lost_single(indices) {
                self.tableau.$name(indices);
                return;
            }
            let filtered: Vec<usize> = indices
                .iter()
                .copied()
                .filter(|&i| !self.is_lost[i])
                .collect();
            self.tableau.$name(&filtered);
        }
    };
}

macro_rules! forward_batch_pair {
    ($name:ident) => {
        #[inline(always)]
        fn $name(&mut self, pairs: &[(usize, usize)]) {
            if !self.any_lost_pair(pairs) {
                self.tableau.$name(pairs);
                return;
            }
            let filtered: Vec<(usize, usize)> = pairs
                .iter()
                .copied()
                .filter(|&(c, t)| !self.is_lost[c] && !self.is_lost[t])
                .collect();
            self.tableau.$name(&filtered);
        }
    };
}

impl<A: RowStorage, I: Bitstring, H> CliffordBatch for GeneralizedTableau<A, I, H> {
    forward_batch_single!(x_many);
    #[inline]
    fn y_many(&mut self, indices: &[usize]) {
        self.tableau.y_many_skipping(indices, &self.is_lost);
    }
    forward_batch_single!(z_many);
    forward_batch_single!(h_many);
    forward_batch_single!(s_many);
    forward_batch_pair!(cnot_many);
    forward_batch_pair!(cz_many);
}

impl<A: RowStorage, I: Bitstring, H> CliffordExtensionsBatch for GeneralizedTableau<A, I, H> {
    forward_batch_single!(s_dag_many);
    forward_batch_single!(sqrt_x_many);
    forward_batch_single!(sqrt_x_dag_many);
    forward_batch_single!(sqrt_y_many);
    forward_batch_single!(sqrt_y_dag_many);
    forward_batch_pair!(cy_many);
}
