// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Maintaining the inverse tableau's signs, and reading the decomposition off
//! them.
//!
//! The storage side — why the inverse's *bits* cost nothing and only its `2n`
//! signs have to be tracked — is [`crate::storage::inverse`]. This module is the
//! algebra: one rule per Clifford generator, and the `O(1)` formula that
//! replaces the fold of `k` whole generators in
//! [`GeneralizedTableau::compute_decomposition`](crate::GeneralizedTableau).
//!
//! # Where the gate rules come from
//!
//! A gate `G` sends the frame's Clifford `U` to `G·U` (every generator is
//! conjugated, `g ↦ G g G†`), so the inverse rows go
//!
//! ```text
//! ix_q = U†X_qU  ↦  U†(G†X_qG)U = Φ(G†X_qG)
//! ```
//!
//! where `Φ = U†·U` is the homomorphism sending `X_j ↦ ix_j`, `Z_j ↦ iz_j`. Each
//! rule is therefore just `G†`'s single-qubit conjugation table, read as a word
//! in `X`/`Z` and pushed through `Φ`. Only the gate's own targets move — Stim's
//! `prepend_*` — so a one-qubit gate is two sign writes when its table is a
//! permutation (`H`, the Paulis, `√Y`), and one `g`-rule pass over the four
//! planes at that qubit when the table produces a `Y` (`S`, `√X`), because
//! `Y = i·X·Z` makes the image a *product* of the qubit's two rows.
//!
//! The tables (`G†PG`, read off the forward kernels in
//! [`crate::storage::blocks`], which implement `g ↦ G g G†`):
//!
//! | `G` | `G†XG` | `G†ZG` |
//! |:--|:--|:--|
//! | `X` | `X` | `−Z` |
//! | `Y` | `−X` | `−Z` |
//! | `Z` | `−X` | `Z` |
//! | `H` | `Z` | `X` |
//! | `S` | `−Y` | `Z` |
//! | `S†` | `Y` | `Z` |
//! | `√X` | `X` | `Y` |
//! | `(√X)†` | `X` | `−Y` |
//! | `√Y` | `Z` | `−X` |
//! | `(√Y)†` | `−Z` | `X` |
//! | `CNOT(c,t)` | `X_cX_t` at `c` | `Z_cZ_t` at `t` |
//! | `CZ(a,b)` | `X_aZ_b` at `a`, `X_bZ_a` at `b` | — |
//! | `CY(c,t)` | `X_cY_t` at `c`, `Z_cX_t` at `t` | `Z_cZ_t` at `t` |
//!
//! Every two-qubit gate here is its own inverse, so one table serves `G` and
//! `G†`; entries not listed are fixed.
//!
//! # Hermitian generators
//!
//! `U` exists only if every generator is Hermitian, which in this crate's
//! convention (a generator's stored `ℤ/4` phase multiplies the *Hermitian*
//! Pauli of its bits) means every phase is even. Clifford gates never touch the
//! low phase plane and the measurement projection only multiplies commuting
//! generators together, so the frame stays Hermitian on its own — but
//! [`StabilizerFrame::row_multiply`](ppvm_traits_2::StabilizerFrame::row_multiply)
//! is public and can multiply two anticommuting generators, which cannot. That
//! call therefore [invalidates](Tableau::invalidate_inverse) the signs, and
//! `inverse_valid ⟹ every phase is even` is an invariant of this module.

use ppvm_traits_2::Pauli;

use crate::data::Tableau;
use crate::storage::{HALVES, Half, InvRow, Orientation, Plane, TableauData, blocks};

impl<H> Tableau<H> {
    /// Whether the inverse-row signs can be read.
    #[inline]
    pub(crate) fn inverse_valid(&self) -> bool {
        self.data.inverse_valid()
    }

    /// Whether the phase of `U†PU` can be read off the signs for `pauli`.
    ///
    /// `X` and `Z` are stored signs, so they need nothing but a current cache —
    /// including under a [`TransposedTableau`](crate::data::TransposedTableau)
    /// guard, which is what keeps a whole `measure_all` off the fold. `Y` is the
    /// exception: it is a *product* of the qubit's two rows, and reading a row's
    /// bits needs the canonical orientation, where a major is an inverse row
    /// rather than a generator ([`TableauData::inv_planes`](crate::storage::TableauData)).
    #[inline]
    pub(crate) fn inverse_readable(&self, pauli: Pauli) -> bool {
        self.data.inverse_valid()
            && (pauli != Pauli::Y || self.data.orientation() == Orientation::ColumnMajor)
    }

    /// Give up on the inverse signs; readers fall back to the row fold.
    ///
    /// Always safe, and the price of a frame mutation with no inverse rule: a
    /// missing rule costs speed, never correctness.
    #[inline]
    pub(crate) fn invalidate_inverse(&mut self) {
        self.data.invalidate_inverse();
    }

    /// The decomposition phase, from one inverse row instead of a fold of `k`
    /// generators.
    ///
    /// `stab_selected` selects stabilizers and `destab_selected` destabilizers —
    /// the crossed naming
    /// [`compute_decomposition`](crate::GeneralizedTableau::compute_decomposition)
    /// uses, because `ω(P, dᵢ)` is what selects `sᵢ`.
    ///
    /// # Derivation
    ///
    /// The fold's statement is `P · T = i^φ · I`, with `T` the ordered product of
    /// the selected generators (stabilizers first). Conjugating by `U` turns `T`
    /// into an *unsigned* word, because `U†sᵢU = Zᵢ` and `U†dᵢU = Xᵢ`:
    ///
    /// ```text
    /// U†(P·T)U = (U†PU) · W,   W = ∏_{i∈S} Zᵢ · ∏_{i∈D} Xᵢ
    /// ```
    ///
    /// `W` and `U†PU` have the same bits (that is what the decomposition
    /// *means*), and a Pauli times itself is the identity, so
    /// `φ = phase(W) + phase(U†PU)`. `phase(W)` is a popcount: the two groups are
    /// ordered `Z`s-then-`X`s, so a qubit in `S ∩ D` contributes `Z·X = i·Y` and
    /// every other site contributes nothing. `phase(U†PU)` enters with a `+`
    /// rather than a `−` because `U†PU` is Hermitian, so its phase is `0` or `2`.
    pub(crate) fn decomposition_phase(
        &self,
        qubit: usize,
        pauli: Pauli,
        stab_selected: &[u64],
        destab_selected: &[u64],
    ) -> u8 {
        debug_assert!(self.inverse_valid());
        let w_phase = (blocks::and_count(stab_selected, destab_selected) % 4) as u8;
        (w_phase + self.data.inv_phase_of(qubit, pauli)) % 4
    }

    /// The deterministic (case-b) outcome of measuring `Z_qubit`: one bit.
    ///
    /// Case b is `ω(Z_q, sᵢ) = 0` for every stabilizer, so
    /// [`Self::decomposition_phase`]'s two masks cannot overlap and the whole
    /// phase *is* the inverse row's sign — `±Z_q` is a stabilizer exactly when
    /// `U†Z_qU = ±Z_q`, and the sign is the outcome. Stim reads the same bit
    /// (`inv_state.zs[q].sign_ref()`).
    #[inline]
    pub(crate) fn inverse_outcome(&self, qubit: usize) -> bool {
        debug_assert!(self.inverse_valid());
        self.data.inv_sign(InvRow::Z, qubit) == 2
    }

    // ─── One-qubit rules ──────────────────────────────────────────────────

    /// Apply a one-qubit rule from `G†`'s table.
    ///
    /// Called **before** the forward sweep: the rules are stated over the
    /// pre-gate rows, whose bits that sweep is about to overwrite.
    #[inline]
    fn prepend1(&mut self, q: usize, rule: Rule1) {
        if !self.inverse_valid() {
            return;
        }
        let px = self.data.inv_sign(InvRow::X, q);
        let pz = self.data.inv_sign(InvRow::Z, q);
        let (nx, nz) = match rule {
            Rule1::Signs(dx, dz) => (px + dx, pz + dz),
            Rule1::Swap(dx, dz) => (pz + dx, px + dz),
            // `Φ(±Y_q) = ±i·ix_q·iz_q`, which is what `inv_phase_of(_, Y)`
            // returns — the one family that costs a pass over the four planes.
            Rule1::YAtX(d) => (self.data.inv_phase_of(q, Pauli::Y) + d, pz),
            Rule1::YAtZ(d) => (px, self.data.inv_phase_of(q, Pauli::Y) + d),
        };
        self.data.set_inv_sign(InvRow::X, q, nx % 4);
        self.data.set_inv_sign(InvRow::Z, q, nz % 4);
    }

    /// `G† : X ↦ Z, Z ↦ X` — the qubit's two inverse rows exchange.
    #[inline]
    pub(crate) fn prepend_h(&mut self, q: usize) {
        self.prepend1(q, Rule1::Swap(0, 0));
    }

    /// `G† : X ↦ −Y, Z ↦ Z`.
    #[inline]
    pub(crate) fn prepend_s(&mut self, q: usize) {
        self.prepend1(q, Rule1::YAtX(2));
    }

    /// `G† : X ↦ Y, Z ↦ Z`.
    #[inline]
    pub(crate) fn prepend_s_dag(&mut self, q: usize) {
        self.prepend1(q, Rule1::YAtX(0));
    }

    /// `G† : X ↦ X, Z ↦ −Z`.
    #[inline]
    pub(crate) fn prepend_x(&mut self, q: usize) {
        self.prepend1(q, Rule1::Signs(0, 2));
    }

    /// `G† : X ↦ −X, Z ↦ −Z`.
    #[inline]
    pub(crate) fn prepend_y(&mut self, q: usize) {
        self.prepend1(q, Rule1::Signs(2, 2));
    }

    /// `G† : X ↦ −X, Z ↦ Z`.
    #[inline]
    pub(crate) fn prepend_z(&mut self, q: usize) {
        self.prepend1(q, Rule1::Signs(2, 0));
    }

    /// `G† : X ↦ X, Z ↦ Y`.
    #[inline]
    pub(crate) fn prepend_sqrt_x(&mut self, q: usize) {
        self.prepend1(q, Rule1::YAtZ(0));
    }

    /// `G† : X ↦ X, Z ↦ −Y`.
    #[inline]
    pub(crate) fn prepend_sqrt_x_dag(&mut self, q: usize) {
        self.prepend1(q, Rule1::YAtZ(2));
    }

    /// `G† : X ↦ Z, Z ↦ −X`.
    #[inline]
    pub(crate) fn prepend_sqrt_y(&mut self, q: usize) {
        self.prepend1(q, Rule1::Swap(0, 2));
    }

    /// `G† : X ↦ −Z, Z ↦ X`.
    #[inline]
    pub(crate) fn prepend_sqrt_y_dag(&mut self, q: usize) {
        self.prepend1(q, Rule1::Swap(2, 0));
    }

    // ─── Two-qubit rules ──────────────────────────────────────────────────

    /// `CNOT†`: `X_c ↦ X_cX_t`, `Z_t ↦ Z_cZ_t`, the other two fixed.
    pub(crate) fn prepend_cnot(&mut self, control: usize, target: usize) {
        if !self.inverse_valid() {
            return;
        }
        let nx = self
            .data
            .inv_pair_phase((InvRow::X, control), (InvRow::X, target));
        let nz = self
            .data
            .inv_pair_phase((InvRow::Z, control), (InvRow::Z, target));
        self.data.set_inv_sign(InvRow::X, control, nx);
        self.data.set_inv_sign(InvRow::Z, target, nz);
    }

    /// `CZ†`: `X_a ↦ X_aZ_b`, `X_b ↦ X_bZ_a`, both `Z`s fixed.
    pub(crate) fn prepend_cz(&mut self, a: usize, b: usize) {
        if !self.inverse_valid() {
            return;
        }
        let na = self.data.inv_pair_phase((InvRow::X, a), (InvRow::Z, b));
        let nb = self.data.inv_pair_phase((InvRow::X, b), (InvRow::Z, a));
        self.data.set_inv_sign(InvRow::X, a, na);
        self.data.set_inv_sign(InvRow::X, b, nb);
    }

    /// `CY†`: `X_c ↦ X_cY_t`, `X_t ↦ Z_cX_t`, `Z_t ↦ Z_cZ_t`, `Z_c` fixed.
    ///
    /// The `X_c` row is the crate's only three-row product,
    /// `Φ(X_cY_t) = ix_c · i·ix_t·iz_t`.
    pub(crate) fn prepend_cy(&mut self, control: usize, target: usize) {
        if !self.inverse_valid() {
            return;
        }
        let nxc = (1 + self.data.inv_triple_phase(
            (InvRow::X, control),
            (InvRow::X, target),
            (InvRow::Z, target),
        )) % 4;
        let nxt = self
            .data
            .inv_pair_phase((InvRow::Z, control), (InvRow::X, target));
        let nzt = self
            .data
            .inv_pair_phase((InvRow::Z, control), (InvRow::Z, target));
        self.data.set_inv_sign(InvRow::X, control, nxc);
        self.data.set_inv_sign(InvRow::X, target, nxt);
        self.data.set_inv_sign(InvRow::Z, target, nzt);
    }
}

// ─── The measurement projection ────────────────────────────────────────────

/// One qubit's four site-planes, as the two inverse families read them.
///
/// A "site-plane" is one bit per inverse row `q`, holding that row's bit at a
/// fixed site — and a site of an inverse row is a *generator index*
/// ([`crate::storage::inverse`]), so a site-plane is one forward generator's bit
/// vector: exactly what [`TableauData::gather_row`] materializes.
///
/// The two families read it differently, which is the whole of the mapping:
///
/// | family | its X plane at site `j` | its Z plane at site `j` |
/// |:--|:--|:--|
/// | `ix` | Z bits of `s_j` | Z bits of `d_j` |
/// | `iz` | X bits of `s_j` | X bits of `d_j` |
struct SitePlanes<'a> {
    /// X and Z bits of the stabilizer at this site.
    stab: (&'a mut [u64], &'a mut [u64]),
    /// X and Z bits of the destabilizer at this site.
    destab: (&'a mut [u64], &'a mut [u64]),
}

impl<'a> SitePlanes<'a> {
    /// Carve four planes off the front of a working buffer.
    fn carve(buf: &'a mut [u64], stride: usize) -> (Self, &'a mut [u64]) {
        let (sx, buf) = buf.split_at_mut(stride);
        let (sz, buf) = buf.split_at_mut(stride);
        let (dx, buf) = buf.split_at_mut(stride);
        let (dz, buf) = buf.split_at_mut(stride);
        (
            Self {
                stab: (sx, sz),
                destab: (dx, dz),
            },
            buf,
        )
    }

    /// Materialize site `j`'s planes out of the frame, in either orientation.
    fn gather(&mut self, data: &TableauData, j: usize) {
        data.gather_row(Half::Stab, j, self.stab.0, self.stab.1);
        data.gather_row(Half::Destab, j, self.destab.0, self.destab.1);
    }
}

impl<H> Tableau<H> {
    /// Carry the inverse signs through a case-a measurement projection.
    ///
    /// Must run **before** the projection mutates the frame: the planes it reads
    /// are the pre-projection generators.
    ///
    /// # The projection is a sequence of appends
    ///
    /// The projection replaces the frame `(dᵢ, sᵢ)` with `sᵢ·s_p` for the selected
    /// stabilizers, `dᵢ·s_p` for the selected destabilizers, `d_p ← s_p` and
    /// `s_p ← (−1)^r Z_a`. Every one of those row operations is a Clifford
    /// **appended** to the frame — `U ↦ U·V` — and appending `V` conjugates each
    /// inverse row by `V†` while leaving the row *indices* alone:
    ///
    /// ```text
    /// ix'_q = U'†X_qU' = V†(U†X_qU)V = V†·ix_q·V
    /// ```
    ///
    /// So this is Stim's `collapse_qubit_z`: the elimination is
    /// `append_CX(p, i)` per selected stabilizer and `append_CZ(i, p)` per
    /// selected destabilizer, then `append_S(p)` if needed and `append_H(p)` to
    /// swap the pivot pair, then the sign of one row is *set* to the outcome.
    ///
    /// Conjugating a Pauli by a Clifford is what the [`blocks`] kernels already
    /// do — they take an X plane, a Z plane and a sign plane and apply the gate
    /// to every row at once. Driven over site-planes instead of qubit columns,
    /// the *same kernels* update the inverse's bits and signs together: one call
    /// per family per append, `O(n/64)` words each.
    ///
    /// # Why no phase is read
    ///
    /// The sign change from conjugating a Pauli by a Clifford depends on that
    /// Pauli's *bits* alone, so the frame's `ℤ/4` phases — which the forward
    /// projection is busy folding `g`-rules into — never enter here. The two
    /// bookkeepings are independent computations of the same frame.
    pub(crate) fn project_inverse(&mut self, addr0: usize, pivot: usize, outcome: bool) {
        debug_assert!(self.data.inverse_valid());
        let n = self.n_qubits();
        let stride = self.data.stride();

        // Ten working planes, reused for the life of the frame: the pivot's four,
        // one site's four, and the two selectors. A case-a measurement is
        // frequent enough that ten allocations per projection showed up in the
        // profile.
        let mut scratch = self.data.take_inv_scratch(10);
        let (mut p, rest) = SitePlanes::carve(&mut scratch, stride);
        let (mut site, rest) = SitePlanes::carve(rest, stride);
        let (destab_sel, rest) = rest.split_at_mut(stride);
        let (stab_sel, _) = rest.split_at_mut(stride);

        // `ω(Z_a, g) = x_g[a]`: one column per half, minus the pivot, which is
        // not multiplied into itself.
        let mut selected = [destab_sel, stab_sel];
        for (half, out) in HALVES.into_iter().zip(selected.iter_mut()) {
            self.data.gather_column(half, Plane::X, addr0, out);
            TableauData::set_bit(out, pivot, false);
        }

        // The pivot's planes persist across every append — it is the one site
        // every gate touches, and `d_p` accumulates the elimination's by-product
        // (Stim's, too) until the final swap turns it into the measured `±Z_a`.
        p.gather(&self.data, pivot);

        // One pass, both appends per site, in the forward sweep's order. Site `i`
        // is gathered exactly once and each of its appends writes only its own
        // planes (`CX` the X plane, `CZ` the Z plane), so the copy is the site's
        // live state throughout — which the second append's sign rule reads. The
        // by-products land on the pivot, which is why the pivot's copy is the one
        // that has to persist.
        //
        // Interleaving is also what makes the order irrelevant: `d_i` and `s_i`
        // anticommute, so a different order gives the accumulated `d_p` a
        // different sign — but that sign is *overwritten* by the outcome below,
        // and every other row is order-independent.
        for i in 0..n {
            let stab = TableauData::bit(selected[Half::Stab as usize], i);
            let destab = TableauData::bit(selected[Half::Destab as usize], i);
            if !stab && !destab {
                continue;
            }
            site.gather(&self.data, i);

            // `append_CX(p, i)`: `s_i ← s_i·s_p`, `d_p ← d_p·d_i`.
            if stab {
                blocks::cnot(
                    p.stab.1,
                    p.destab.1,
                    site.stab.1,
                    site.destab.1,
                    self.data.inv_sign_plane_mut(InvRow::X),
                );
                blocks::cnot(
                    p.stab.0,
                    p.destab.0,
                    site.stab.0,
                    site.destab.0,
                    self.data.inv_sign_plane_mut(InvRow::Z),
                );
            }

            // `append_CZ(i, p)`: `d_i ← d_i·s_p`, `d_p ← d_p·s_i`. `CZ` is
            // symmetric in both its bit map and its sign rule, so the endpoint
            // order is free.
            if destab {
                blocks::cz(
                    site.stab.1,
                    site.destab.1,
                    p.stab.1,
                    p.destab.1,
                    self.data.inv_sign_plane_mut(InvRow::X),
                );
                blocks::cz(
                    site.stab.0,
                    site.destab.0,
                    p.stab.0,
                    p.destab.0,
                    self.data.inv_sign_plane_mut(InvRow::Z),
                );
            }
        }

        // Every generator but the pivot pair now commutes with `Z_a`, so `Z_a`
        // expands over that pair alone: `Z_a = i^φ·d_p·s_p^{ω(Z_a, d_p)}`. When
        // that exponent is set, `append_S(p)` (`d_p ← i·d_p·s_p`) clears it, and
        // then `d_p` **is** `±Z_a`.
        if TableauData::bit(p.destab.0, addr0) {
            blocks::s_dag(
                p.stab.1,
                p.destab.1,
                self.data.inv_sign_plane_mut(InvRow::X),
            );
            blocks::s_dag(
                p.stab.0,
                p.destab.0,
                self.data.inv_sign_plane_mut(InvRow::Z),
            );
        }

        // `append_H(p)` swaps the pivot pair, which is the projection's
        // `d_p ← s_p`, `s_p ← ±Z_a`.
        blocks::h(
            p.stab.1,
            p.destab.1,
            self.data.inv_sign_plane_mut(InvRow::X),
        );
        blocks::h(
            p.stab.0,
            p.destab.0,
            self.data.inv_sign_plane_mut(InvRow::Z),
        );
        debug_assert!(
            p.stab.0.iter().all(|&w| w == 0)
                && blocks::first_set(p.stab.1) == Some(addr0)
                && blocks::count_set(p.stab.1) == 1,
            "the appends must leave the pivot stabilizer at ±Z_{addr0}"
        );

        // The outcome is not implied by the appends — it is the *choice* the
        // projection made. One row pins it: `s_p = (−1)^r Z_a` means
        // `iz_a = U'†Z_aU' = (−1)^r Z_p`, so that row's sign is the outcome. If
        // it disagrees, `append_X(p)` flips the sign of every inverse row with a
        // `Z` at site `p` — `iz_a` among them — which is exactly negating `s_p`.
        if (self.data.inv_sign(InvRow::Z, addr0) == 2) != outcome {
            blocks::pauli_x(p.destab.1, self.data.inv_sign_plane_mut(InvRow::X));
            blocks::pauli_x(p.destab.0, self.data.inv_sign_plane_mut(InvRow::Z));
        }
        self.data.restore_inv_scratch(scratch);
    }
}

/// A one-qubit `G†` table, in the shape the sign update needs.
#[derive(Clone, Copy, Debug)]
enum Rule1 {
    /// `X ↦ ±X`, `Z ↦ ±Z`: the rows keep their identity, only signs move.
    Signs(u8, u8),
    /// `X ↦ ±Z`, `Z ↦ ±X`: the rows exchange.
    Swap(u8, u8),
    /// `X ↦ ±Y`, `Z ↦ Z`.
    YAtX(u8),
    /// `X ↦ X`, `Z ↦ ±Y`.
    YAtZ(u8),
}

/// Test-only oracle: every inverse row, conjugated *forward* through the frame,
/// must give back the basis Pauli it came from.
///
/// This checks the rules against their definition (`U·(U†PU)·U† = P`) rather
/// than against a rearrangement of the formula they feed, so a sign error in any
/// rule fails it. [`Self::apply_frame`] is the independent leg: it substitutes
/// `X_i ↦ dᵢ`, `Z_i ↦ sᵢ` with the same `g`-rule the fold uses.
#[cfg(test)]
impl<H> Tableau<H> {
    /// `U R U†` for `R = (x, z, phase)`, a phased Pauli over the frame's qubits.
    pub(crate) fn apply_frame(&self, x: &[u64], z: &[u64], phase: u8) -> (Vec<u64>, Vec<u64>, u8) {
        use crate::data::ScratchRow;
        use crate::storage::Half;

        let stride = self.data.stride();
        let mut acc = ScratchRow::zeroed(stride);
        acc.add_phase(phase);
        let mut src = ScratchRow::zeroed(stride);
        for i in 0..self.n_qubits() {
            let (xi, zi) = (TableauData::bit(x, i), TableauData::bit(z, i));
            if xi && zi {
                // `Y_i = i·X_i·Z_i`: substitute both factors, keep the `i`.
                acc.add_phase(1);
            }
            if xi {
                acc.mul_generator(&self.data, Half::Destab, i, &mut src);
            }
            if zi {
                acc.mul_generator(&self.data, Half::Stab, i, &mut src);
            }
        }
        (acc.x, acc.z, acc.phase)
    }

    /// Re-derive all `2n` signs from the frame and declare them current.
    ///
    /// `U·(U†PU)·U† = P` pins each sign: conjugating the row's bits forward with
    /// phase `0` gives `i^φ·P`, so the row's own phase must be `φ` (even, hence
    /// its own negation). `O(n³/64)` — a test tool, and the reason a projection
    /// abandons the signs rather than rebuilding them.
    pub(crate) fn rebuild_inverse_signs(&mut self) {
        for q in 0..self.n_qubits() {
            for row in [InvRow::X, InvRow::Z] {
                let (rx, rz) = self.data.inv_planes(row, q);
                let (rx, rz) = (rx.to_vec(), rz.to_vec());
                let phase = self.apply_frame(&rx, &rz, 0).2;
                self.data.set_inv_sign(row, q, phase);
            }
        }
        self.data.revalidate_inverse();
    }

    /// Assert the inverse signs describe this frame.
    pub(crate) fn assert_inverse_consistent(&self) {
        assert!(self.inverse_valid(), "the inverse signs were abandoned");
        let stride = self.data.stride();
        for q in 0..self.n_qubits() {
            for (row, pauli) in [(InvRow::X, Pauli::X), (InvRow::Z, Pauli::Z)] {
                let (rx, rz) = self.data.inv_planes(row, q);
                let (rx, rz) = (rx.to_vec(), rz.to_vec());
                let got = self.apply_frame(&rx, &rz, self.data.inv_sign(row, q));
                let mut want = (vec![0u64; stride], vec![0u64; stride], 0);
                TableauData::set_bit(
                    if pauli == Pauli::X {
                        &mut want.0
                    } else {
                        &mut want.1
                    },
                    q,
                    true,
                );
                assert_eq!(got, want, "U·(U†{pauli:?}_{q}U)·U† is not {pauli:?}_{q}");
            }
        }
    }
}
