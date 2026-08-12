// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Clifford operations on the frame: the symplectic/phase primitives, the
//! **fused** [`Clifford`]/[`CliffordExtensions`] impls, the frame primitives,
//! and the batched entry points.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits";
//! `tableau-data-structure.md` §"Column-major orientation" and §"Gate access
//! patterns".
//!
//! # One gate, one contiguous sweep
//!
//! In the canonical column-major orientation ([`crate::storage`]) a qubit's X
//! and Z bits over all `2n` generators are contiguous, so every kernel here is
//! a straight-line loop over `n.div_ceil(64)` words per half — `h(q)` swaps two
//! columns and XORs their `AND` into the sign plane, `cnot(c, t)` touches four
//! columns. The replaced implementation walked all `2n` generators to touch one
//! machine word in each, which at `n = 1889` meant streaming ~2 MB to change
//! `2n` bits.
//!
//! The per-gate *predicates* are unchanged — see [`blocks`], which carries the
//! bit and sign tables and their Lean references. Only the access pattern moved.
//!
//! # Why `Tableau` opts *out* of [`BlanketClifford`](ppvm_traits_2::BlanketClifford)
//!
//! The blanket runs the phase primitive and the column primitives as *separate*
//! steps. Each is still a full sweep of the affected columns, so `h` would cost
//! two passes instead of one and `cnot`/`cz` three instead of one, re-reading
//! bits the fused body already has in registers. `Tableau` therefore does what
//! `Phased<W>` does for the same reason (design §"The `BlanketClifford` marker
//! and the fused phased override"): it implements [`SymplecticColumns`] and
//! [`PhaseTrack`] so the primitive contract is honoured and auditable, stays out
//! of the marker, and supplies its own fused impls below.
//!
//! # Batched entry points
//!
//! The replaced layout amortized its `2n`-generator walk by building a per-word
//! qubit mask and sweeping once for a whole list of targets. With the walk gone
//! there is nothing left to amortize: `x_many(&[q0, q1])` is two column sweeps
//! either way. The batch methods are therefore loops over the single-qubit
//! kernels — but over the **distinct** targets, because that, not the sequential
//! loop, is what the mask sweep computed: a repeated index set its mask bit once
//! and so flipped the phase once. Preserving that is deliberate. It is the
//! mechanism behind `G-061` in `docs/lean-gap.md` — a legal `.stim` file with a
//! repeated single-qubit target — and a layout rewrite is not the place to
//! change a behaviour the ledger is still adjudicating.
//!
//! Pair batches read fresh bits per pair in both layouts, so those are plain
//! sequential loops.
//!
//! Each generator is an `Sp(2n, 2)` isometry (machine-checked per generator in
//! `lean/PPVM/Pauli/Symplectic.lean`: `hAct_isometry`, `sAct_isometry`,
//! `cnotAct_isometry`, `czAct_isometry`) with the sign action of
//! `lean/PPVM/Pauli/Conjugation.lean` (`conjH_sign`, `conjS_sign`, …), and it
//! maps a symplectic frame to a symplectic frame
//! (`lean/PPVM/Tableau/Frame.lean`, `IsSymplecticFrame.map`).

use ppvm_traits_2::{
    Clifford, CliffordBatch, CliffordExtensions, CliffordExtensionsBatch, PhaseTrack,
    StabilizerFrame, SymplecticColumns,
};

use crate::data::{Bitstring, GeneralizedTableau, ScratchRow, Tableau, TransposedTableau};
use crate::storage::{HALVES, Half, Plane, TableauData, blocks};

/// Run a one-qubit column kernel over both halves.
///
/// The optional `inverse: <rule>` clause names the [`crate::inverse`] sign rule
/// for this gate, and it runs **first** — the rules are stated over the pre-gate
/// rows. A sweep without one is not a Clifford (the bare `SymplecticColumns` /
/// `PhaseTrack` primitives), so it abandons the inverse signs instead.
macro_rules! sweep1 {
    ($tab:expr, $q:expr, inverse: $rule:ident, $kernel:expr) => {{
        $tab.$rule($q);
        sweep1!($tab, $q, $kernel)
    }};
    ($tab:expr, $q:expr, $kernel:expr) => {{
        $tab.invalidate_hash();
        for half in HALVES {
            let (x, z, ph) = $tab.data.gate1_mut(half, $q);
            #[allow(clippy::redundant_closure_call)]
            $kernel(x, z, ph);
        }
    }};
}

/// Run a two-qubit column kernel over both halves. `inverse:` as [`sweep1`].
macro_rules! sweep2 {
    ($tab:expr, $a:expr, $b:expr, inverse: $rule:ident, $kernel:expr) => {{
        $tab.$rule($a, $b);
        sweep2!($tab, $a, $b, $kernel)
    }};
    ($tab:expr, $a:expr, $b:expr, $kernel:expr) => {{
        $tab.invalidate_hash();
        for half in HALVES {
            let (xa, za, xb, zb, ph) = $tab.data.gate2_mut(half, $a, $b);
            #[allow(clippy::redundant_closure_call)]
            $kernel(xa, za, xb, zb, ph);
        }
    }};
}

// ─── Sp / phase primitives ────────────────────────────────────────────────

/// The `Sp`-part of conjugation: whole-column bit-plane algebra over all `2n`
/// generators. `Tableau` deliberately does **not** opt into
/// [`BlanketClifford`](ppvm_traits_2::BlanketClifford) (see the module note);
/// the primitives are supplied so the design's audited contract holds and so a
/// mixture algebra keyed on `Tableau` can reach them.
impl<H> SymplecticColumns for Tableau<H> {
    #[inline]
    fn n_qubits(&self) -> usize {
        self.data.n_qubits()
    }

    #[inline]
    fn swap_xz(&mut self, q: usize) {
        self.invalidate_inverse();
        sweep1!(self, q, |x: &mut [u64], z: &mut [u64], _ph: &mut [u64]| {
            x.swap_with_slice(z);
        });
    }

    #[inline]
    fn xor_z_from_x(&mut self, q: usize) {
        self.invalidate_inverse();
        sweep1!(self, q, |x: &mut [u64], z: &mut [u64], _ph: &mut [u64]| {
            for (zw, &xw) in z.iter_mut().zip(x.iter()) {
                *zw ^= xw;
            }
        });
    }

    #[inline]
    fn xor_x_col(&mut self, ctrl: usize, tgt: usize) {
        self.invalidate_inverse();
        sweep2!(
            self,
            ctrl,
            tgt,
            |xc: &mut [u64], _zc: &mut [u64], xt: &mut [u64], _zt: &mut [u64], _ph: &mut [u64]| {
                for (dst, &src) in xt.iter_mut().zip(xc.iter()) {
                    *dst ^= src;
                }
            }
        );
    }

    #[inline]
    fn xor_z_col(&mut self, tgt: usize, ctrl: usize) {
        self.invalidate_inverse();
        sweep2!(
            self,
            ctrl,
            tgt,
            |_xc: &mut [u64], zc: &mut [u64], _xt: &mut [u64], zt: &mut [u64], _ph: &mut [u64]| {
                for (dst, &src) in zc.iter_mut().zip(zt.iter()) {
                    *dst ^= src;
                }
            }
        );
    }

    #[inline]
    fn cz_bits(&mut self, a: usize, b: usize) {
        self.invalidate_inverse();
        sweep2!(self, a, b, |xa: &mut [u64],
                             za: &mut [u64],
                             xb: &mut [u64],
                             zb: &mut [u64],
                             _ph: &mut [u64]| {
            for i in 0..za.len() {
                let (ax, bx) = (xa[i], xb[i]);
                za[i] ^= bx;
                zb[i] ^= ax;
            }
        });
    }
}

/// The extension part: the tableau's phase algebra is a `ℤ₂` sign — the high
/// bit plane of the `ℤ/4` phase — plus the Aaronson–Gottesman `g` rule, which
/// lives behind [`StabilizerFrame::row_multiply`].
impl<H> PhaseTrack for Tableau<H> {
    #[inline]
    fn flip_phase_where_xz(&mut self, q: usize) {
        self.invalidate_inverse();
        sweep1!(self, q, |x: &mut [u64], z: &mut [u64], ph: &mut [u64]| {
            for i in 0..ph.len() {
                ph[i] ^= x[i] & z[i];
            }
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
        self.invalidate_inverse();
        sweep2!(
            self,
            ctrl,
            tgt,
            |xc: &mut [u64], zc: &mut [u64], xt: &mut [u64], zt: &mut [u64], ph: &mut [u64]| {
                for i in 0..ph.len() {
                    ph[i] ^= xc[i] & zt[i] & !(xt[i] ^ zc[i]);
                }
            }
        );
    }

    #[inline]
    fn cz_phase(&mut self, a: usize, b: usize) {
        self.invalidate_inverse();
        sweep2!(self, a, b, |xa: &mut [u64],
                             za: &mut [u64],
                             xb: &mut [u64],
                             zb: &mut [u64],
                             ph: &mut [u64]| {
            for i in 0..ph.len() {
                ph[i] ^= xa[i] & xb[i] & (za[i] ^ zb[i]);
            }
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

// ─── Frame primitives ─────────────────────────────────────────────────────

impl<H> StabilizerFrame for Tableau<H> {
    /// `ω(Z_q, sᵢ) = x_{sᵢ}[q]`, so the pivot is the lowest set bit of one
    /// contiguous column.
    #[inline]
    fn anticommuting_pivot(&self, qubit: usize) -> Option<usize> {
        self.find_z_anticommuting_stabilizer(qubit)
    }

    /// Multiply generator `src` into `dst` (Aaronson–Gottesman `g`-rule).
    ///
    /// A whole-generator fold, so it takes the [`TransposedTableau`] guard: one
    /// transpose in, one out. Callers with a *run* of row multiplies should hold
    /// their own guard rather than paying that per call — which is what the
    /// measurement projection does.
    fn row_multiply(&mut self, src: usize, dst: usize) {
        assert_ne!(src, dst, "row_multiply needs distinct rows");
        let n = SymplecticColumns::n_qubits(self);
        self.invalidate_hash();
        // Two anticommuting generators have a non-Hermitian product, and then no
        // `U` exists for the inverse signs to be the signs *of*.
        self.invalidate_inverse();
        let stride = self.data.stride();
        let mut guard = TransposedTableau::new(self);
        let data = guard.data_mut();

        let (src_half, src_i) = Half::split(src, n);
        let (dst_half, dst_i) = Half::split(dst, n);
        let mut row = ScratchRow::zeroed(stride);
        row.x.copy_from_slice(data.major(src_half, Plane::X, src_i));
        row.z.copy_from_slice(data.major(src_half, Plane::Z, src_i));
        row.phase = data.phase_of(src_half, src_i);
        data.multiply_row_by(dst_half, dst_i, &row.x, &row.z, row.phase);
    }

    /// No-op: this representation never leaves canonical form.
    ///
    /// Every gate is an `Sp(2n, 2)` isometry on the frame
    /// (`isSymplecticFrame_*` in `lean/PPVM/Tableau/Frame.lean`) and the
    /// measurement projection restores the destabilizer/stabilizer pairing in
    /// place (`update_tableau_according_to_outcome`), so the `2n` generators are
    /// a symplectic basis after every public operation and there is nothing left
    /// to restore.
    ///
    /// The projection half is **not** covered by `IsSymplecticFrame.map` (that
    /// lemma needs an `ω`-isometry, which the projection is not — it overwrites
    /// two rows). It is machine-checked separately as
    /// `isSymplecticFrame_projectFrame` in `lean/PPVM/Tableau/Frame.lean`, whose
    /// `projectFrame` is exactly this crate's row sweep.
    ///
    /// Kept as an explicit no-op rather than dropped so a caller written against
    /// the trait is portable to a backend that *does* defer canonicalization.
    #[inline]
    fn canonicalize(&mut self) {}
}

// ─── Fused Clifford ───────────────────────────────────────────────────────

/// Single source of truth for the per-gate Clifford phase/bit logic. Every
/// caller — a bare `Tableau`, a `GeneralizedTableau` (which delegates here) and
/// the batch paths — runs through this one implementation, so there is no
/// parallel copy that can silently diverge.
impl<H> Clifford for Tableau<H> {
    #[inline]
    fn x(&mut self, qubit: usize) {
        sweep1!(self, qubit, inverse: prepend_x, |_x: &mut [u64],
                              z: &mut [u64],
                              ph: &mut [u64]| {
            blocks::pauli_x(z, ph)
        });
    }

    #[inline]
    fn y(&mut self, qubit: usize) {
        sweep1!(self, qubit, inverse: prepend_y, |x: &mut [u64],
                              z: &mut [u64],
                              ph: &mut [u64]| {
            blocks::pauli_y(x, z, ph)
        });
    }

    #[inline]
    fn z(&mut self, qubit: usize) {
        sweep1!(self, qubit, inverse: prepend_z, |x: &mut [u64],
                              _z: &mut [u64],
                              ph: &mut [u64]| {
            blocks::pauli_z(x, ph)
        });
    }

    #[inline]
    fn h(&mut self, qubit: usize) {
        sweep1!(self, qubit, inverse: prepend_h, blocks::h);
    }

    #[inline]
    fn s(&mut self, qubit: usize) {
        // NOTE: S is the only Clifford where forward and backward propagation
        // differ (it is non-Hermitian); only the phase rule differs.
        sweep1!(self, qubit, inverse: prepend_s, |x: &mut [u64],
                              z: &mut [u64],
                              ph: &mut [u64]| {
            blocks::s(x, z, ph)
        });
    }

    #[inline]
    fn cnot(&mut self, control: usize, target: usize) {
        sweep2!(
            self,
            control,
            target,
            inverse: prepend_cnot,
            |xc: &mut [u64], zc: &mut [u64], xt: &mut [u64], zt: &mut [u64], ph: &mut [u64]| {
                blocks::cnot(xc, zc, xt, zt, ph)
            }
        );
    }

    #[inline]
    fn cz(&mut self, qubit0: usize, qubit1: usize) {
        sweep2!(
            self,
            qubit0,
            qubit1,
            inverse: prepend_cz,
            |xa: &mut [u64], za: &mut [u64], xb: &mut [u64], zb: &mut [u64], ph: &mut [u64]| {
                blocks::cz(xa, za, xb, zb, ph)
            }
        );
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
impl<H> CliffordExtensions for Tableau<H> {
    #[inline]
    fn s_dag(&mut self, qubit: usize) {
        sweep1!(self, qubit, inverse: prepend_s_dag, |x: &mut [u64],
                              z: &mut [u64],
                              ph: &mut [u64]| {
            blocks::s_dag(x, z, ph)
        });
    }

    #[inline]
    fn sqrt_x(&mut self, qubit: usize) {
        sweep1!(self, qubit, inverse: prepend_sqrt_x, |x: &mut [u64],
                              z: &mut [u64],
                              ph: &mut [u64]| {
            blocks::sqrt_x(x, z, ph)
        });
    }

    #[inline]
    fn sqrt_x_dag(&mut self, qubit: usize) {
        sweep1!(self, qubit, inverse: prepend_sqrt_x_dag, |x: &mut [u64],
                              z: &mut [u64],
                              ph: &mut [u64]| {
            blocks::sqrt_x_dag(x, z, ph)
        });
    }

    #[inline]
    fn sqrt_y(&mut self, qubit: usize) {
        sweep1!(self, qubit, inverse: prepend_sqrt_y, blocks::sqrt_y);
    }

    #[inline]
    fn sqrt_y_dag(&mut self, qubit: usize) {
        sweep1!(self, qubit, inverse: prepend_sqrt_y_dag, blocks::sqrt_y_dag);
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
        sweep2!(
            self,
            control,
            target,
            inverse: prepend_cy,
            |xc: &mut [u64], zc: &mut [u64], xt: &mut [u64], zt: &mut [u64], ph: &mut [u64]| {
                blocks::cy(xc, zc, xt, zt, ph)
            }
        );
    }
}

// ─── Batched entry points ─────────────────────────────────────────────────

impl<H> Tableau<H> {
    /// Apply a one-qubit gate once per **distinct** target, in order of first
    /// appearance.
    ///
    /// See the module note: this is the semantics the replaced per-word mask
    /// sweep had, and repeated targets are the whole difference. Distinct
    /// one-qubit gates touch disjoint columns and XOR into the sign plane, so
    /// the visit order among them is immaterial.
    #[inline]
    fn for_each_distinct(&mut self, indices: &[usize], gate: impl Fn(&mut Self, usize)) {
        if indices.is_empty() {
            return;
        }
        let mut seen = vec![0u64; self.data.stride()];
        for &q in indices {
            if TableauData::bit(&seen, q) {
                continue;
            }
            TableauData::set_bit(&mut seen, q, true);
            gate(self, q);
        }
    }
}

impl<H> CliffordBatch for Tableau<H> {
    #[inline]
    fn x_many(&mut self, indices: &[usize]) {
        self.for_each_distinct(indices, Clifford::x);
    }

    #[inline]
    fn y_many(&mut self, indices: &[usize]) {
        self.for_each_distinct(indices, Clifford::y);
    }

    #[inline]
    fn z_many(&mut self, indices: &[usize]) {
        self.for_each_distinct(indices, Clifford::z);
    }

    #[inline]
    fn s_many(&mut self, indices: &[usize]) {
        self.for_each_distinct(indices, Clifford::s);
    }

    #[inline]
    fn h_many(&mut self, indices: &[usize]) {
        self.for_each_distinct(indices, Clifford::h);
    }

    /// Pairs are applied in order, each reading fresh bits, so overlapping pairs
    /// behave exactly as the per-pair loop — as they did before.
    fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
        for &(control, target) in pairs {
            Clifford::cnot(self, control, target);
        }
    }

    /// Pairs are applied in order; `CZ` is symmetric and touches only z-bits.
    fn cz_many(&mut self, pairs: &[(usize, usize)]) {
        for &(a, b) in pairs {
            Clifford::cz(self, a, b);
        }
    }
}

impl<H> CliffordExtensionsBatch for Tableau<H> {
    #[inline]
    fn s_dag_many(&mut self, indices: &[usize]) {
        self.for_each_distinct(indices, CliffordExtensions::s_dag);
    }

    fn cy_many(&mut self, pairs: &[(usize, usize)]) {
        for &(control, target) in pairs {
            CliffordExtensions::cy(self, control, target);
        }
    }

    #[inline]
    fn sqrt_y_many(&mut self, indices: &[usize]) {
        self.for_each_distinct(indices, CliffordExtensions::sqrt_y);
    }

    #[inline]
    fn sqrt_y_dag_many(&mut self, indices: &[usize]) {
        self.for_each_distinct(indices, CliffordExtensions::sqrt_y_dag);
    }

    #[inline]
    fn sqrt_x_many(&mut self, indices: &[usize]) {
        self.for_each_distinct(indices, CliffordExtensions::sqrt_x);
    }

    #[inline]
    fn sqrt_x_dag_many(&mut self, indices: &[usize]) {
        self.for_each_distinct(indices, CliffordExtensions::sqrt_x_dag);
    }
}

// ─── GeneralizedTableau: loss-aware forwarding ────────────────────────────

// Single-qubit gate: skip a lost qubit, otherwise delegate to the frame's
// canonical (column-sweep) method.
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
        #[inline]
        fn $name(&mut self, control: usize, target: usize) {
            if self.is_lost[control] || self.is_lost[target] {
                return;
            }
            self.tableau.$name(control, target);
        }
    };
}

impl<I: Bitstring, H> Clifford for GeneralizedTableau<I, H> {
    forward_single!(x);
    forward_single!(y);
    forward_single!(z);
    forward_single!(h);
    forward_single!(s);
    forward_pair!(cnot);
    forward_pair!(cz);
}

impl<I: Bitstring, H> CliffordExtensions for GeneralizedTableau<I, H> {
    forward_single!(s_dag);
    forward_single!(sqrt_x);
    forward_single!(sqrt_x_dag);
    forward_single!(sqrt_y);
    forward_single!(sqrt_y_dag);
    forward_pair!(cy);
}

impl<I: Bitstring, H> GeneralizedTableau<I, H> {
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

impl<I: Bitstring, H> CliffordBatch for GeneralizedTableau<I, H> {
    forward_batch_single!(x_many);
    forward_batch_single!(y_many);
    forward_batch_single!(z_many);
    forward_batch_single!(h_many);
    forward_batch_single!(s_many);
    forward_batch_pair!(cnot_many);
    forward_batch_pair!(cz_many);
}

impl<I: Bitstring, H> CliffordExtensionsBatch for GeneralizedTableau<I, H> {
    forward_batch_single!(s_dag_many);
    forward_batch_single!(sqrt_x_many);
    forward_batch_single!(sqrt_x_dag_many);
    forward_batch_single!(sqrt_y_many);
    forward_batch_single!(sqrt_y_dag_many);
    forward_batch_pair!(cy_many);
}
