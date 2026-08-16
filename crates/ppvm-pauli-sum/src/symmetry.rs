// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lattice translation symmetry groups for operator-space Pauli evolution.
//!
//! A [`TranslationGroup`] represents a finite abelian group `G` acting on
//! qubit positions by permutations. Given such a group, every Pauli word
//! belongs to a translation orbit, and operator dynamics that commute
//! with `G` can be tracked using **one canonical representative per
//! orbit** instead of all `|G|` orbit members — reducing per-step memory
//! and compute by a factor up to `|G|`.
//!
//! Following Teng, Chang, Rudolph, and Holmes (arXiv:2512.12094), this
//! module implements **plain (real-coefficient) merging** of Pauli sums
//! into orbit-representative form — see [`canonicalize_pauli_sum`] and
//! [`symmetry_merge_pauli_sum`]. This handles observables in the trivial
//! (`k=0`) symmetry sector, e.g. sums of single-Z operators over the
//! lattice.
//!
//! **Non-trivial momentum sectors (`k ≠ 0`)** are handled by
//! [`canonicalize_pauli_sum_complex`], which folds with the character
//! phase `χ_k(g)` of each translation. On the Python side, an operator in
//! sector `k` is carried as a *real pair* (real + imaginary components, two
//! real `PauliSum`s) and merged via `PauliSum.momentum_merge`, which reuses
//! this routine — letting gate-based Trotter evolution stay symmetry-
//! compressed in any momentum sector with real coefficients throughout.
//!
//! ## Data model
//!
//! A `TranslationGroup` is specified by a list of generator permutations
//! and their cyclic orders. The group order is the product of the orders.
//! For instance, a 2D `L × L` torus has two generators (translation in
//! x and y) each of order `L`.
//!
//! ## Canonicalization
//!
//! [`TranslationGroup::canonicalize`] returns the **lex-minimum** Pauli
//! word reachable from the input via group action. The ordering is the
//! standard `Ord` impl on `PauliWord` (compare `xbits`, then `zbits`).
//! All orbit members canonicalize to the same representative; orbits are
//! disjoint by construction, so the rep uniquely identifies the orbit.
//!
//! For the chain and ladder layouts — a single generator shifting
//! contiguous blocks of qubits — this is computed in `O(N)` by a staged
//! least-rotation (Booth/Duval) scan rather than by walking all `|G|`
//! group elements; see `canonicalize_block_cyclic`. Other groups fall
//! back to an `O(|G| × N)` odometer walk, which also serves as the test
//! oracle for the fast path. Both return the identical representative.
//!
//! ## Merging
//!
//! [`canonicalize_pauli_sum`] takes parallel `Vec<Word>` / `Vec<f64>`
//! buffers (the representation used by ppvm-lindblad's adaptive
//! evolution) and replaces each Pauli by its canonical rep, summing
//! coefficients for collisions. The output is an orbit-rep basis with
//! coefficients equal to the sum of the input coefficients over each
//! orbit's members. For dynamics that commute with `G` and initial
//! states that are also `G`-invariant, this preserves the expectation
//! value of any `G`-invariant observable (Theorem 1 of arXiv:2512.12094).
//!
//! See the dedicated tests for correctness against full-basis evolution
//! on small systems with no truncation.

use crate::sum::PauliSum;
use fxhash::FxHashMap;
use num::Complex;
use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::Config;
use ppvm_traits::{HashFinalize, PauliStorage, PauliWordTrait};
use std::f64::consts::PI;
use std::hash::BuildHasher;

/// A finite abelian symmetry group acting on qubit positions by
/// permutations.
///
/// Build via the convenience constructors [`Self::chain_1d`],
/// [`Self::torus_2d`], [`Self::torus_3d`], [`Self::ladder`], or
/// [`Self::from_generators`] for an arbitrary list of generator
/// permutations.
///
/// `perms[g]` is the permutation that **generator `g`** applies to qubit
/// indices: a qubit at position `q` moves to position `perms[g][q]`
/// under one application of generator `g`. `orders[g]` is the cyclic
/// order of generator `g` (i.e. applying it `orders[g]` times returns
/// the identity). The full group is the direct product of the cyclic
/// subgroups, with size `Π orders[g]`.
///
/// Only the **generators** are stored; [`Self::canonicalize`] either runs
/// the `O(N)` least-rotation scan (chain/ladder layouts) or walks the
/// group via mixed-radix increments.
#[derive(Debug, Clone)]
pub struct TranslationGroup {
    /// Number of qubits the group acts on.
    n_qubits: usize,
    /// One permutation per generator. `perms[g][q]` is the position
    /// that qubit `q` maps to under one application of generator `g`.
    perms: Vec<Vec<u32>>,
    /// Cyclic order of each generator.
    orders: Vec<u32>,
    /// Set when the group is a *single* generator acting as a cyclic
    /// shift inside contiguous, aligned blocks of qubits — i.e. exactly
    /// the [`Self::chain_1d`] and [`Self::ladder`] layouts. Enables the
    /// `O(N)` least-rotation canonicalizer (see
    /// [`Self::canonicalize_block_cyclic`]).
    block_cyclic: Option<BlockCyclic>,
}

/// Layout of a single-generator group acting as a cyclic shift within
/// `n_blocks` contiguous, aligned blocks of `len` qubits each: qubit
/// `b * len + j` maps to `b * len + (j + 1) % len`.
#[derive(Debug, Clone, Copy)]
struct BlockCyclic {
    n_blocks: usize,
    len: usize,
}

/// Detect the [`BlockCyclic`] layout, if the generators have it.
fn detect_block_cyclic(n_qubits: usize, perms: &[Vec<u32>], orders: &[u32]) -> Option<BlockCyclic> {
    if perms.len() != 1 {
        return None;
    }
    let len = orders[0] as usize;
    if len == 0 || n_qubits == 0 || !n_qubits.is_multiple_of(len) {
        return None;
    }
    let n_blocks = n_qubits / len;
    let perm = &perms[0];
    for b in 0..n_blocks {
        for j in 0..len {
            if perm[b * len + j] as usize != b * len + (j + 1) % len {
                return None;
            }
        }
    }
    Some(BlockCyclic { n_blocks, len })
}

/// Start index of the lexicographically smallest rotation of an abstract
/// `m`-symbol cyclic sequence, via the two-pointer (Booth/Duval) scan.
///
/// `cmp(a, b)` compares the symbols at positions `a` and `b`. `O(m)`
/// comparisons, no allocation.
fn least_rotation<F>(m: usize, cmp: &F) -> usize
where
    F: Fn(usize, usize) -> std::cmp::Ordering,
{
    let (mut i, mut j, mut k) = (0usize, 1usize, 0usize);
    while i < m && j < m && k < m {
        match cmp((i + k) % m, (j + k) % m) {
            std::cmp::Ordering::Equal => {
                k += 1;
                continue;
            }
            std::cmp::Ordering::Greater => i += k + 1,
            std::cmp::Ordering::Less => j += k + 1,
        }
        if i == j {
            j += 1;
        }
        k = 0;
    }
    i.min(j)
}

/// Period of the cyclic sequence `t ↦ start + t (mod m)` — the smallest
/// `p` dividing `m` with `s[t] == s[t + p]` for all `t`.
///
/// Computed as the length of the first Lyndon factor (Duval): the minimal
/// rotation of a sequence is a power `w^{m/|w|}` of a Lyndon word `w`, and
/// `|w|` is the period. `O(m)` comparisons, no allocation. Callers pass the
/// `start` returned by [`least_rotation`]; the count of rotations achieving
/// the minimum is then `m / period`, spaced `period` apart.
fn minimal_rotation_period<F>(m: usize, start: usize, cmp: &F) -> usize
where
    F: Fn(usize, usize) -> std::cmp::Ordering,
{
    let at = |t: usize| (start + t) % m;
    let (mut j, mut k) = (1usize, 0usize);
    while j < m {
        match cmp(at(k), at(j)) {
            std::cmp::Ordering::Less => {
                k = 0;
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                k += 1;
                j += 1;
            }
            std::cmp::Ordering::Greater => break,
        }
    }
    let len = j - k;
    if m.is_multiple_of(len) { len } else { m }
}

impl TranslationGroup {
    /// Construct from explicit generator permutations and orders.
    ///
    /// Each `perm` must be a permutation of `0..n_qubits`. Each `order`
    /// must satisfy `perm^order == identity`.
    pub fn from_generators(n_qubits: usize, perms: Vec<Vec<u32>>, orders: Vec<u32>) -> Self {
        assert_eq!(perms.len(), orders.len(), "perms and orders must match");
        for (g, perm) in perms.iter().enumerate() {
            assert_eq!(
                perm.len(),
                n_qubits,
                "generator {g} permutation has length {} != n_qubits {n_qubits}",
                perm.len()
            );
            let mut seen = vec![false; n_qubits];
            for &p in perm {
                assert!(
                    (p as usize) < n_qubits,
                    "generator {g} maps to out-of-range position {p}"
                );
                assert!(
                    !seen[p as usize],
                    "generator {g} is not a permutation (duplicate target {p})"
                );
                seen[p as usize] = true;
            }
        }
        let block_cyclic = detect_block_cyclic(n_qubits, &perms, &orders);
        Self {
            n_qubits,
            perms,
            orders,
            block_cyclic,
        }
    }

    /// 1D chain of `n` sites with periodic boundary conditions.
    /// Single generator: cyclic shift by one site.
    pub fn chain_1d(n: usize) -> Self {
        let perm: Vec<u32> = (0..n).map(|q| ((q + 1) % n) as u32).collect();
        Self::from_generators(n, vec![perm], vec![n as u32])
    }

    /// 2D `lx × ly` torus, qubit at `(i, j)` indexed as `j*lx + i`.
    /// Two generators: x-shift (i → i+1 mod lx) and y-shift (j → j+1 mod ly).
    pub fn torus_2d(lx: usize, ly: usize) -> Self {
        let n = lx * ly;
        let perm_x: Vec<u32> = (0..n)
            .map(|q| {
                let (i, j) = (q % lx, q / lx);
                (j * lx + (i + 1) % lx) as u32
            })
            .collect();
        let perm_y: Vec<u32> = (0..n)
            .map(|q| {
                let (i, j) = (q % lx, q / lx);
                (((j + 1) % ly) * lx + i) as u32
            })
            .collect();
        Self::from_generators(n, vec![perm_x, perm_y], vec![lx as u32, ly as u32])
    }

    /// 3D `lx × ly × lz` torus, qubit at `(i, j, k)` indexed as
    /// `k*lx*ly + j*lx + i`.
    pub fn torus_3d(lx: usize, ly: usize, lz: usize) -> Self {
        let n = lx * ly * lz;
        let perm_x: Vec<u32> = (0..n)
            .map(|q| {
                let i = q % lx;
                let j = (q / lx) % ly;
                let k = q / (lx * ly);
                (k * lx * ly + j * lx + (i + 1) % lx) as u32
            })
            .collect();
        let perm_y: Vec<u32> = (0..n)
            .map(|q| {
                let i = q % lx;
                let j = (q / lx) % ly;
                let k = q / (lx * ly);
                (k * lx * ly + ((j + 1) % ly) * lx + i) as u32
            })
            .collect();
        let perm_z: Vec<u32> = (0..n)
            .map(|q| {
                let i = q % lx;
                let j = (q / lx) % ly;
                let k = q / (lx * ly);
                (((k + 1) % lz) * lx * ly + j * lx + i) as u32
            })
            .collect();
        Self::from_generators(
            n,
            vec![perm_x, perm_y, perm_z],
            vec![lx as u32, ly as u32, lz as u32],
        )
    }

    /// Multi-leg ladder: `l` sites along the chain × `n_legs` legs.
    /// Single generator: cyclic shift along the chain direction (all
    /// legs simultaneously). Qubit at `(leg, j)` indexed as
    /// `leg * l + j`. No translation along the leg axis (legs are
    /// distinguished).
    pub fn ladder(l: usize, n_legs: usize) -> Self {
        let n = l * n_legs;
        let perm: Vec<u32> = (0..n)
            .map(|q| {
                let leg = q / l;
                let j = q % l;
                (leg * l + (j + 1) % l) as u32
            })
            .collect();
        Self::from_generators(n, vec![perm], vec![l as u32])
    }

    /// Number of qubits the group acts on.
    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Number of generators (rank of the group as an abelian product).
    pub fn n_generators(&self) -> usize {
        self.perms.len()
    }

    /// Total group order: `Π orders[g]`.
    pub fn order(&self) -> usize {
        self.orders.iter().map(|&o| o as usize).product()
    }

    /// Permutation associated with the `g`-th generator (one application).
    pub fn generator_perm(&self, g: usize) -> &[u32] {
        &self.perms[g]
    }

    /// Cyclic order of the `g`-th generator.
    pub fn generator_order(&self, g: usize) -> u32 {
        self.orders[g]
    }

    /// Apply a single generator's permutation to a Pauli word, returning
    /// the resulting word.
    ///
    /// For each qubit `q` of the input, the corresponding `(xbit, zbit)`
    /// pair is placed at position `perm[q]` of the output.
    fn apply_generator<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
        g: usize,
    ) -> PauliWord<A, S, R>
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        let perm = &self.perms[g];
        let mut out: PauliWord<A, S, R> = PauliWord::new(self.n_qubits);
        for (q, &pq) in perm.iter().enumerate().take(self.n_qubits) {
            let xb = w.get_xbit(q);
            let zb = w.get_zbit(q);
            if xb {
                out.set_xbit(pq as usize, true);
            }
            if zb {
                out.set_zbit(pq as usize, true);
            }
        }
        out.rehash();
        out
    }

    /// Odometer step: advance `cur` from the group element with
    /// mixed-radix index `idx - 1` to the one with index `idx`.
    ///
    /// Generator `0` is the fastest-varying digit, so it advances on
    /// every step; digit `g` advances only when all lower digits roll
    /// over, i.e. when `idx` is a multiple of `orders[0..=g-1]`. Applying
    /// generator `g` once always moves digit `g` forward *cyclically*
    /// (the `orders[g]`-th application is the identity), so a roll-over
    /// is just one more application — no rebuild from the identity.
    ///
    /// Cost: `O(1)` generator applications amortised, hence `O(|G| × N)`
    /// for a full walk instead of the `O(|G|² × N)` of rebuilding each
    /// element from scratch.
    #[inline]
    fn advance<A, S, const R: bool>(&self, cur: &mut PauliWord<A, S, R>, idx: usize)
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        let mut p = 1usize;
        for (g, &o) in self.orders.iter().enumerate() {
            *cur = self.apply_generator(cur, g);
            p *= o as usize;
            if !idx.is_multiple_of(p) {
                break;
            }
        }
    }

    /// Lex-min canonical representative of `w`'s translation orbit
    /// under this group.
    ///
    /// For chain/ladder layouts this is `O(N)` via the least-rotation
    /// canonicalizer ([`Self::canonicalize_block_cyclic`]); otherwise it
    /// walks the full group as a mixed-radix odometer, `O(|G| × N)`.
    pub fn canonicalize<A, S, const R: bool>(&self, w: &PauliWord<A, S, R>) -> PauliWord<A, S, R>
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        debug_assert_eq!(
            w.n_qubits(),
            self.n_qubits,
            "word and group must agree on n_qubits"
        );
        if self.perms.is_empty() {
            return *w;
        }
        if let Some(bc) = self.block_cyclic {
            return self.canonicalize_block_cyclic(w, bc).0;
        }
        self.canonicalize_odometer(w).0
    }

    /// Reference canonicalizer: walk the whole group as a mixed-radix
    /// odometer (see [`Self::advance`]), keeping the smallest word seen.
    /// Returns the rep and the index of the group element mapping it back
    /// to `w`. `O(|G| × N)`; used for groups without a
    /// [`BlockCyclic`] layout, and as the test oracle for the fast path.
    fn canonicalize_odometer<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
    ) -> (PauliWord<A, S, R>, usize)
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        let mut best = *w;
        let mut best_idx = 0usize;
        let mut cur = *w;
        for idx in 1..self.order() {
            self.advance(&mut cur, idx);
            if cur < best {
                best = cur;
                best_idx = idx;
            }
        }
        // The walk found `best = g·w` at index `best_idx`, so `w = g⁻¹·best`
        // and the element we must report is the inverse. In an abelian
        // product of cyclic groups that is `(orders[g] − c[g]) mod orders[g]`
        // componentwise.
        (best, self.invert_index(best_idx))
    }

    /// `O(N)` canonicalizer for single-generator cyclic-block groups
    /// (chain, ladder): returns the same rep as the odometer walk — the
    /// `Ord`-lex-min of the orbit — and the index `r` of the group element
    /// with `g^r · rep = w`.
    ///
    /// ## Why this is not one Booth call
    ///
    /// `PauliWord`'s `Ord` compares the whole x-bit plane in qubit order,
    /// *then* the whole z-bit plane. Under a shift by `r`, the comparison
    /// key is therefore the concatenation
    /// `rot_r(x_block0) ‖ … ‖ rot_r(z_block0) ‖ …` — `2 · n_blocks` strings
    /// rotated *together*, not one rotated string, so lex-min over rotations
    /// is not a single least-rotation problem. (Running Booth on an
    /// interleaved per-site symbol would be one call, but it minimises a
    /// different order and so would silently change which orbit member is
    /// canonical.)
    ///
    /// Instead we refine the candidate rotation set plane by plane. After
    /// each plane the surviving rotations form a residue class
    /// `{start + i·step}` of size `m = L / step`, because the rotations
    /// achieving a minimum are exactly those spaced by the *period* of that
    /// minimal rotation. Plane `p + 1` then compares its own string only at
    /// those rotations — which is again a least-rotation problem, over `m`
    /// super-symbols of `step` bits each. Every plane costs `O(L)` symbol
    /// comparisons of `O(step)` bits = `O(L)`, so the whole call is
    /// `O(n_blocks · L) = O(N)`, allocation-free apart from the output word.
    fn canonicalize_block_cyclic<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
        bc: BlockCyclic,
    ) -> (PauliWord<A, S, R>, usize)
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        let l = bc.len;
        // Surviving rotations: { (start + i·step) mod l : i < m }, with
        // step · m == l throughout, and `start < step` (the smallest one).
        let (mut start, mut step, mut m) = (0usize, 1usize, l);
        for plane in 0..2 * bc.n_blocks {
            if m == 1 {
                break;
            }
            let is_x = plane < bc.n_blocks;
            let base = (if is_x { plane } else { plane - bc.n_blocks }) * l;
            // Symbol `j` is the run of `step` bits of this plane starting at
            // rotation offset `start + j·step`.
            let bit = |j: usize, t: usize| -> bool {
                let pos = base + (start + j * step + t) % l;
                if is_x {
                    w.get_xbit(pos)
                } else {
                    w.get_zbit(pos)
                }
            };
            let cmp = |a: usize, b: usize| -> std::cmp::Ordering {
                for t in 0..step {
                    let (x, y) = (bit(a, t), bit(b, t));
                    if x != y {
                        // `false < true`, matching bit-slice lex order.
                        return x.cmp(&y);
                    }
                }
                std::cmp::Ordering::Equal
            };
            let j0 = least_rotation(m, &cmp);
            let period = minimal_rotation_period(m, j0, &cmp);
            start = (start + j0 * step) % l;
            step *= period;
            m /= period;
            start %= step; // smallest member of the surviving residue class
        }
        // Tie-break exactly as the odometer does: it keeps the *first*
        // minimal word it meets, i.e. the smallest number of generator
        // applications `idx = (l − r) mod l`. That is `r = 0` when `r = 0`
        // survives, and otherwise the largest surviving `r`.
        let r = if start == 0 {
            0
        } else {
            start + (m - 1) * step
        };
        // rep = g^{−r}·w, i.e. rep[base + j] = w[base + (j + r) mod l].
        let mut rep: PauliWord<A, S, R> = PauliWord::new(self.n_qubits);
        for b in 0..bc.n_blocks {
            let base = b * l;
            for j in 0..l {
                let src = base + (j + r) % l;
                if w.get_xbit(src) {
                    rep.set_xbit(base + j, true);
                }
                if w.get_zbit(src) {
                    rep.set_zbit(base + j, true);
                }
            }
        }
        rep.rehash();
        (rep, r)
    }

    /// Lex-min canonical representative `r` of `w` together with the
    /// **mixed-radix counter** `c = (c_0, c_1, …)` of the group element
    /// `g` such that `g·r = w`.
    ///
    /// In other words: if `r = self.canonicalize(w)`, this returns
    /// `(r, c)` where applying generator `i` exactly `c[i]` times in
    /// sequence to `r` produces `w`. The counter is used to compute
    /// momentum phases by the phase-aware merge routines.
    ///
    /// Same `O(|G| × n_qubits)` cost as `canonicalize`.
    ///
    /// Allocates the counter `Vec`. Hot paths that only need the momentum
    /// phase should prefer [`Self::canonicalize_with_index`] together with
    /// [`Self::character_table`] — allocation-free, and no transcendental
    /// per call.
    pub fn canonicalize_with_shift<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
    ) -> (PauliWord<A, S, R>, Vec<u32>)
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        if self.perms.is_empty() {
            return (*w, Vec::new());
        }
        let (rep, inv_idx) = self.canonicalize_with_index(w);
        (rep, self.counter_from_index(inv_idx))
    }

    /// Lex-min canonical representative `r` of `w` together with the
    /// **mixed-radix index** of the group element `g` such that `g·r = w`
    /// — i.e. the index of the counter returned by
    /// [`Self::canonicalize_with_shift`].
    ///
    /// The index is directly usable as a subscript into
    /// [`Self::character_table`], which is how the phase-aware evolution
    /// gets `χ_k(g)` without decoding a counter or calling `sin`/`cos`
    /// per term.
    ///
    /// Cost: `O(N)` for chain/ladder layouts, else `O(|G| × N)`.
    /// Allocation-free apart from the returned word.
    pub fn canonicalize_with_index<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
    ) -> (PauliWord<A, S, R>, usize)
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        debug_assert_eq!(w.n_qubits(), self.n_qubits);
        if self.perms.is_empty() {
            return (*w, 0);
        }
        match self.block_cyclic {
            Some(bc) => self.canonicalize_block_cyclic(w, bc),
            None => self.canonicalize_odometer(w),
        }
    }

    /// Decode a group-element index (mixed-radix, generator `0` fastest)
    /// into its per-generator counter.
    pub fn counter_from_index(&self, idx: usize) -> Vec<u32> {
        let mut rem = idx;
        let mut counter: Vec<u32> = Vec::with_capacity(self.perms.len());
        for &o in &self.orders {
            counter.push((rem as u32) % o);
            rem /= o as usize;
        }
        counter
    }

    /// Index of the inverse of the group element with index `idx`.
    fn invert_index(&self, idx: usize) -> usize {
        let mut rem = idx;
        let mut out = 0usize;
        let mut stride = 1usize;
        for &o in &self.orders {
            let c = (rem as u32) % o;
            rem /= o as usize;
            out += (((o - c) % o) as usize) * stride;
            stride *= o as usize;
        }
        out
    }

    /// All `|G|` momentum-sector characters, indexed by group-element
    /// index: `table[idx] == self.character(k_modes,
    /// &self.counter_from_index(idx))`.
    ///
    /// Build this once per evolution step and index it with the value from
    /// [`Self::canonicalize_with_index`]; the alternative — calling
    /// [`Self::character`] per action term — costs a `sin`/`cos` pair and a
    /// counter `Vec` every time.
    pub fn character_table(&self, k_modes: &[i32]) -> Vec<Complex<f64>> {
        assert_eq!(
            k_modes.len(),
            self.perms.len(),
            "k_modes length {} != number of generators {}",
            k_modes.len(),
            self.perms.len()
        );
        (0..self.order())
            .map(|idx| {
                let mut rem = idx;
                let mut phase = 0.0_f64;
                for (&k, &o) in k_modes.iter().zip(self.orders.iter()) {
                    let c = (rem as u32) % o;
                    rem /= o as usize;
                    phase += 2.0 * PI * (k as f64) * (c as f64) / (o as f64);
                }
                Complex::from_polar(1.0, phase)
            })
            .collect()
    }

    /// Momentum-sector character `χ_k(g) = exp(i Σ_g 2π · k[g] · counter[g] / orders[g])`
    /// where `k[g] ∈ ℤ` is the integer momentum mode along generator `g`
    /// (the corresponding wavenumber is `2π · k[g] / orders[g]`).
    ///
    /// `k.len()` must equal `self.n_generators()`. The character of the
    /// identity element (`counter = [0, …]`) is `1`. For the trivial
    /// (`k = [0, …]`) sector all characters are `1` — phase-aware merging
    /// reduces to plain merging.
    pub fn character(&self, k_modes: &[i32], counter: &[u32]) -> Complex<f64> {
        debug_assert_eq!(k_modes.len(), self.perms.len());
        debug_assert_eq!(counter.len(), self.perms.len());
        let mut phase = 0.0_f64;
        for ((&k, &c), &o) in k_modes.iter().zip(counter.iter()).zip(self.orders.iter()) {
            phase += 2.0 * PI * (k as f64) * (c as f64) / (o as f64);
        }
        Complex::from_polar(1.0, phase)
    }

    /// Iterate over all group elements applied to `w`, in mixed-radix
    /// index order. Yields `|G|` Pauli words (`w` itself first, for the
    /// identity element).
    ///
    /// Walks the odometer incrementally (see [`Self::advance`]): `O(|G| × N)`
    /// for the whole orbit.
    pub fn orbit<'a, A, S, const R: bool>(
        &'a self,
        w: &'a PauliWord<A, S, R>,
    ) -> impl Iterator<Item = PauliWord<A, S, R>> + 'a
    where
        A: PauliStorage + 'a,
        S: BuildHasher + Clone + Default + HashFinalize + 'a,
    {
        let mut cur = *w;
        (0..self.order()).map(move |idx| {
            if idx > 0 {
                self.advance(&mut cur, idx);
            }
            cur
        })
    }
}

/// Replace `(basis, coeffs)` in-place with the orbit-representative
/// form: each Pauli word becomes its canonical rep, and coefficients
/// of words that collapse to the same rep are summed.
///
/// Output length ≤ input length. Entries whose summed coefficient
/// equals zero exactly are *not* removed — caller should run a final
/// `drop_tol` prune if desired.
///
/// For dynamics that commute with `group` and initial states that are
/// `group`-invariant (i.e. in the trivial momentum sector), this
/// preserves all `G`-invariant expectation values.
pub fn canonicalize_pauli_sum<A, S, const R: bool>(
    basis: &mut Vec<PauliWord<A, S, R>>,
    coeffs: &mut Vec<f64>,
    group: &TranslationGroup,
) where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    assert_eq!(
        basis.len(),
        coeffs.len(),
        "basis and coeffs length mismatch"
    );
    let mut merged: FxHashMap<PauliWord<A, S, R>, f64> =
        FxHashMap::with_capacity_and_hasher(basis.len(), Default::default());
    for (w, &c) in basis.iter().zip(coeffs.iter()) {
        let rep = group.canonicalize(w);
        *merged.entry(rep).or_insert(0.0) += c;
    }
    basis.clear();
    coeffs.clear();
    basis.reserve(merged.len());
    coeffs.reserve(merged.len());
    for (w, c) in merged {
        basis.push(w);
        coeffs.push(c);
    }
}

/// Replace `(basis, complex_coeffs)` in-place with the orbit-rep form
/// **projected onto momentum sector `k_modes`**.
///
/// Each Pauli `p` is replaced by its canonical rep `r`; the contribution
/// is `(1/|G|) · χ_k(g) · c_p` where `g` is the group element such that
/// `g · r = p` and `χ_k(g) = exp(2πi · Σ_g k_modes[g] · counter[g] / orders[g])`.
///
/// If the input was already a momentum-`k_modes` eigenstate (i.e. the
/// coefficients satisfy `c_{g·p} = χ_k(g)⁻¹ · c_p` for every orbit),
/// the output is the orbit-rep coefficients of that state unchanged.
/// Otherwise the merge discards the components in other sectors —
/// use [`check_momentum_sector`] beforehand to validate.
///
/// For the `k_modes = [0, 0, …]` (trivial) sector this reduces to plain
/// [`canonicalize_pauli_sum`] (real coefficients work, but on complex
/// input the result is complex with vanishing imaginary part).
pub fn canonicalize_pauli_sum_complex<A, S, const R: bool>(
    basis: &mut Vec<PauliWord<A, S, R>>,
    coeffs: &mut Vec<Complex<f64>>,
    group: &TranslationGroup,
    k_modes: &[i32],
) where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    assert_eq!(
        basis.len(),
        coeffs.len(),
        "basis and coeffs length mismatch"
    );
    assert_eq!(
        k_modes.len(),
        group.n_generators(),
        "k_modes length {} != number of generators {}",
        k_modes.len(),
        group.n_generators()
    );
    let inv_g: f64 = 1.0 / (group.order() as f64);
    let chi_table = group.character_table(k_modes);
    let mut merged: FxHashMap<PauliWord<A, S, R>, Complex<f64>> =
        FxHashMap::with_capacity_and_hasher(basis.len(), Default::default());
    for (w, &c) in basis.iter().zip(coeffs.iter()) {
        let (rep, idx) = group.canonicalize_with_index(w);
        let chi = chi_table[idx];
        let contrib = inv_g * chi * c;
        *merged.entry(rep).or_insert(Complex::new(0.0, 0.0)) += contrib;
    }
    basis.clear();
    coeffs.clear();
    basis.reserve(merged.len());
    coeffs.reserve(merged.len());
    for (w, c) in merged {
        basis.push(w);
        coeffs.push(c);
    }
}

/// Verify that a `(basis, complex_coeffs)` Pauli sum lies entirely in
/// the momentum sector `k_modes` under `group`.
///
/// Concretely: for every orbit represented in the basis, all members
/// must satisfy `c_{g·r} = χ_k(g)⁻¹ · c_r` for some choice of orbit-rep
/// coefficient `c_r`.
///
/// Returns `Ok(())` on pass; `Err(SectorCheckError)` on fail with the
/// offending orbit-rep, expected coefficient, and actual coefficient.
///
/// Use this on a user-supplied initial state before feeding it to a
/// phase-aware merging pipeline — silently projecting a wrongly-typed
/// input throws away meaningful physics.
pub fn check_momentum_sector<A, S, const R: bool>(
    basis: &[PauliWord<A, S, R>],
    coeffs: &[Complex<f64>],
    group: &TranslationGroup,
    k_modes: &[i32],
    tol: f64,
) -> Result<(), SectorCheckError<A, S, R>>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    assert_eq!(basis.len(), coeffs.len());
    assert_eq!(k_modes.len(), group.n_generators());

    // Group entries by orbit rep, picking the first-seen member as
    // reference and checking later members against it.
    let chi_table = group.character_table(k_modes);
    let mut reference: FxHashMap<PauliWord<A, S, R>, (Complex<f64>, Vec<u32>)> =
        FxHashMap::default();
    for (p, &c) in basis.iter().zip(coeffs.iter()) {
        let (rep, idx) = group.canonicalize_with_index(p);
        let cnt = group.counter_from_index(idx);
        let chi = chi_table[idx];
        // expected c_p given the rep coefficient c_r:
        //   c_p = χ_k(g)⁻¹ · c_r,  where p = g·r
        // equivalently, c_r = χ_k(g) · c_p (a rearrangement).
        let implied_rep_coeff = chi * c;
        if let Some((rep_coeff, _ref_cnt)) = reference.get(&rep) {
            if (implied_rep_coeff - rep_coeff).norm() > tol * rep_coeff.norm().max(1.0) {
                return Err(SectorCheckError {
                    rep,
                    expected: *rep_coeff,
                    got_implied: implied_rep_coeff,
                    offending_pauli: *p,
                    offending_coeff: c,
                    shift: cnt.clone(),
                });
            }
        } else {
            reference.insert(rep, (implied_rep_coeff, cnt));
        }
    }
    Ok(())
}

/// Detail report for a failed [`check_momentum_sector`].
pub struct SectorCheckError<A: PauliStorage, S, const R: bool> {
    /// Canonical orbit representative for which the check failed.
    pub rep: PauliWord<A, S, R>,
    /// Coefficient that the *first* basis entry implied for `rep`.
    pub expected: Complex<f64>,
    /// Coefficient that `offending_pauli` implies for `rep` under the
    /// purported momentum sector.
    pub got_implied: Complex<f64>,
    /// The basis entry whose coefficient is inconsistent with the
    /// expected `rep` value.
    pub offending_pauli: PauliWord<A, S, R>,
    /// Original coefficient of `offending_pauli` in the input basis.
    pub offending_coeff: Complex<f64>,
    /// Counter encoding the group element `g` such that
    /// `g · rep == offending_pauli`.
    pub shift: Vec<u32>,
}

impl<A: PauliStorage, S, const R: bool> std::fmt::Debug for SectorCheckError<A, S, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SectorCheckError {{ rep: <Word>, expected: {:?}, got_implied: {:?}, \
             offending: <Word>, offending_coeff: {:?}, shift: {:?} }}",
            self.expected, self.got_implied, self.offending_coeff, self.shift,
        )
    }
}

impl<A: PauliStorage, S, const R: bool> std::fmt::Display for SectorCheckError<A, S, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "input not in target momentum sector: orbit rep expected c={:?}, but \
             orbit member (shift {:?}, coeff {:?}) implies c={:?}",
            self.expected, self.shift, self.offending_coeff, self.got_implied,
        )
    }
}

/// Symmetry-merge a [`PauliSum`] in place: each Pauli word becomes its
/// canonical orbit representative, and entries collapsing to the same
/// rep accumulate coefficients.
///
/// This is the Trotter-mode counterpart to [`canonicalize_pauli_sum`]
/// (which operates on the `Vec<Word>, Vec<f64>` representation used by
/// `ppvm-lindblad`'s adaptive evolution). Same semantics: preserves all
/// `G`-invariant expectation values when the dynamics commutes with
/// `group` and the initial state is `group`-invariant.
///
/// Generic over the [`Config`] but constrained to PauliWord-backed
/// representations (i.e. not the loss-aware variant) since
/// canonicalization needs raw `(xbit, zbit)` access.
pub fn symmetry_merge_pauli_sum<T, A, S, const R: bool>(
    psum: &mut PauliSum<T>,
    group: &TranslationGroup,
) where
    T: Config<PauliWordType = PauliWord<A, S, R>>,
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    psum.map_add(|word, coeff| (group.canonicalize(word), coeff.clone()));
}

#[cfg(test)]
mod tests {
    use super::*;

    type W = PauliWord<[u8; 1], fxhash::FxBuildHasher, true>;

    fn word(s: &str) -> W {
        W::from(s)
    }

    #[test]
    fn chain_1d_canonicalizes_via_cyclic_shift() {
        let g = TranslationGroup::chain_1d(4);
        // All cyclic shifts of "IIXY" should canonicalize to the same rep.
        let candidates = ["IIXY", "IXYI", "XYII", "YIIX"];
        let canon: Vec<W> = candidates
            .iter()
            .map(|s| g.canonicalize(&word(s)))
            .collect();
        for c in &canon[1..] {
            assert_eq!(
                *c, canon[0],
                "all cyclic shifts must canonicalize to same rep"
            );
        }
    }

    #[test]
    fn chain_1d_canonicalize_is_lex_min() {
        let g = TranslationGroup::chain_1d(4);
        let canon = g.canonicalize(&word("YIIX"));
        let orbit: Vec<W> = g.orbit(&word("YIIX")).collect();
        let min = orbit.iter().min().unwrap();
        assert_eq!(canon, *min);
    }

    #[test]
    fn orbit_has_correct_size_for_chain() {
        let g = TranslationGroup::chain_1d(4);
        // "XIII" has orbit of size 4 (full chain).
        let orbit: Vec<W> = g.orbit(&word("XIII")).collect();
        assert_eq!(orbit.len(), 4);
        // "XIXI" has orbit of size 2 (period-2 invariant); 4 elements
        // total in the orbit iterator, but only 2 unique.
        let orbit: Vec<W> = g.orbit(&word("XIXI")).collect();
        assert_eq!(orbit.len(), 4); // iterator yields |G|, including duplicates
        let unique: std::collections::HashSet<W> = orbit.into_iter().collect();
        assert_eq!(unique.len(), 2);
    }

    #[test]
    fn torus_2d_canonicalize() {
        // 3x2 torus, 6 qubits.
        let g = TranslationGroup::torus_2d(3, 2);
        assert_eq!(g.n_qubits(), 6);
        assert_eq!(g.order(), 6);
        // X at (0,0) — orbit is all 6 single-X positions.
        let w = word("XIIIII");
        let orbit: Vec<W> = g.orbit(&w).collect();
        let unique: std::collections::HashSet<W> = orbit.into_iter().collect();
        assert_eq!(unique.len(), 6);
        // All canonicalize to the same rep.
        let canon = g.canonicalize(&w);
        for u in &unique {
            assert_eq!(g.canonicalize(u), canon);
        }
    }

    #[test]
    fn ladder_canonicalize() {
        // 2-leg ladder, L=3 → 6 qubits, group order 3 (no swap of legs).
        let g = TranslationGroup::ladder(3, 2);
        assert_eq!(g.n_qubits(), 6);
        assert_eq!(g.order(), 3);
        // X on leg 0 site 0: orbit = {(0,0), (0,1), (0,2)}, NOT including leg 1 sites.
        let w = word("XIIIII"); // qubit 0 = X
        let orbit: Vec<W> = g.orbit(&w).collect();
        assert_eq!(orbit.len(), 3);
        let unique: std::collections::HashSet<W> = orbit.into_iter().collect();
        assert_eq!(unique.len(), 3);
        // The orbit should be {qubit 0=X, qubit 1=X, qubit 2=X} — all leg 0.
        let expected: std::collections::HashSet<W> = ["XIIIII", "IXIIII", "IIXIII"]
            .iter()
            .map(|s| word(s))
            .collect();
        assert_eq!(unique, expected);
    }

    #[test]
    fn canonicalize_pauli_sum_merges_orbit_members() {
        let g = TranslationGroup::chain_1d(4);
        let mut basis: Vec<W> = vec![word("XIII"), word("IXII"), word("IIXI"), word("IIIX")];
        let mut coeffs: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];
        canonicalize_pauli_sum(&mut basis, &mut coeffs, &g);
        // All four collapse to one rep with coeff 1+2+3+4 = 10.
        assert_eq!(basis.len(), 1);
        assert!((coeffs[0] - 10.0).abs() < 1e-12);
    }

    #[test]
    fn canonicalize_pauli_sum_keeps_distinct_orbits() {
        let g = TranslationGroup::chain_1d(4);
        // Two distinct orbits: {XIII, ...} (size 4) and {ZIII, ...} (size 4).
        let mut basis: Vec<W> = vec![word("XIII"), word("IXII"), word("ZIII"), word("IZII")];
        let mut coeffs: Vec<f64> = vec![1.0, 1.0, 2.0, 2.0];
        canonicalize_pauli_sum(&mut basis, &mut coeffs, &g);
        assert_eq!(basis.len(), 2);
        // Coefficients should be {2.0, 4.0} in some order.
        let mut cs = coeffs.clone();
        cs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((cs[0] - 2.0).abs() < 1e-12);
        assert!((cs[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn canonicalize_with_shift_round_trip() {
        // For each cyclic shift of "IIXY" by `a` positions, the shift
        // counter returned should reproduce the original word when
        // applied to the canonical rep.
        let g = TranslationGroup::chain_1d(4);
        for src in ["IIXY", "IXYI", "XYII", "YIIX"] {
            let w = word(src);
            let (rep, cnt) = g.canonicalize_with_shift(&w);
            // Apply gen 0 `cnt[0]` times to rep, should equal w.
            let mut cur = rep;
            for _ in 0..cnt[0] {
                cur = g.apply_generator(&cur, 0);
            }
            assert_eq!(cur, w, "shift {cnt:?} doesn't reproduce {src}");
        }
    }

    #[test]
    fn block_cyclic_canonicalizer_matches_the_odometer() {
        // The O(N) least-rotation path must return *bit-identical* results
        // to the O(|G|·N) group walk — same rep AND same shift index, since
        // the shift index sets the momentum phase. Stabilised words (words
        // with a nontrivial period) are the interesting case: there the
        // minimising rotation is not unique and the two paths must agree on
        // which one to report.
        type W32 = PauliWord<[u8; 4], fxhash::FxBuildHasher, true>;
        let alphabet = ['I', 'X', 'Z', 'Y'];
        let mut rng = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = move || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };
        for (g, n) in [
            (TranslationGroup::chain_1d(6), 6),
            (TranslationGroup::chain_1d(7), 7), // prime order: no proper periods
            (TranslationGroup::chain_1d(8), 8),
            (TranslationGroup::ladder(5, 2), 10),
            (TranslationGroup::ladder(6, 2), 12),
            (TranslationGroup::ladder(6, 3), 18),
        ] {
            assert!(g.block_cyclic.is_some(), "expected the fast path for n={n}");
            let l = g.order();
            let mut cases: Vec<String> = Vec::new();
            // Structured words: empty planes, and every period dividing L.
            cases.push("I".repeat(n));
            cases.push("Z".repeat(n)); // fully stabilised
            cases.push("X".repeat(n));
            for p in 1..=l {
                if l % p == 0 {
                    // period-p pattern, repeated over the whole register
                    let cell: String = (0..p).map(|j| alphabet[(j + 1) % 4]).collect();
                    let mut s = String::new();
                    while s.len() < n {
                        s.push_str(&cell);
                    }
                    s.truncate(n);
                    cases.push(s);
                }
            }
            // Random words, including sparse ones (few non-identity sites).
            for _ in 0..400 {
                let r = next();
                let sparse = r & 1 == 0;
                let s: String = (0..n)
                    .map(|q| {
                        let v = (next() >> (q % 32)) as usize;
                        if sparse && !v.is_multiple_of(4) {
                            'I'
                        } else {
                            alphabet[v % 4]
                        }
                    })
                    .collect();
                cases.push(s);
            }
            for s in cases {
                let w = W32::from(s.as_str());
                // Compare on every orbit member, not just the seed: the two
                // paths must agree pointwise, which also pins canonicality.
                for member in g.orbit(&w) {
                    let fast = g.canonicalize_with_index(&member);
                    let slow = g.canonicalize_odometer(&member);
                    assert_eq!(fast.0, slow.0, "rep mismatch on {s} (n={n})");
                    assert_eq!(fast.1, slow.1, "shift mismatch on {s} (n={n})");
                    // …and the reported shift really maps the rep back.
                    let mut cur = fast.0;
                    for _ in 0..fast.1 {
                        cur = g.apply_generator(&cur, 0);
                    }
                    assert_eq!(cur, member, "shift {} does not reproduce {s}", fast.1);
                }
            }
        }
    }

    #[test]
    fn multi_generator_groups_keep_the_odometer_path() {
        // No cyclic-block layout ⇒ no fast path, and canonicalization must
        // still be correct (covered by torus_2d_canonicalize).
        assert!(TranslationGroup::torus_2d(2, 3).block_cyclic.is_none());
        assert!(TranslationGroup::torus_3d(2, 2, 2).block_cyclic.is_none());
        // A single generator that permutes qubits in a 4-cycle, but not
        // the block-aligned `j → j+1` one: 0→2→1→3→0.
        let g = TranslationGroup::from_generators(4, vec![vec![2u32, 3, 1, 0]], vec![4]);
        assert!(g.block_cyclic.is_none());
        // …and it still canonicalizes correctly through the odometer.
        for member in g.orbit(&word("XZII")) {
            assert_eq!(g.canonicalize(&member), g.canonicalize(&word("XZII")));
        }
    }

    #[test]
    fn canonicalize_with_index_matches_shift_and_table() {
        // The index form must agree with the counter form on both the
        // rep and the momentum phase, for every element of a group with
        // two generators of different orders (mixed-radix odometer).
        let g = TranslationGroup::torus_2d(2, 3);
        assert_eq!(g.order(), 6);
        for k in [[0, 0], [1, 0], [0, 2], [1, 1]] {
            let table = g.character_table(&k);
            assert_eq!(table.len(), g.order());
            for (idx, chi) in table.iter().enumerate() {
                let cnt = g.counter_from_index(idx);
                assert!((chi - g.character(&k, &cnt)).norm() < 1e-12);
            }
            for src in ["XIIIII", "IXZIII", "IIYIXI", "ZIIIIY"] {
                let w = word(src);
                let (rep_s, cnt) = g.canonicalize_with_shift(&w);
                let (rep_i, idx) = g.canonicalize_with_index(&w);
                assert_eq!(rep_s, rep_i, "reps disagree for {src}");
                assert_eq!(
                    cnt,
                    g.counter_from_index(idx),
                    "counters disagree for {src}"
                );
                assert!((table[idx] - g.character(&k, &cnt)).norm() < 1e-12);
            }
        }
    }

    #[test]
    fn odometer_walk_covers_the_whole_group() {
        // `orbit` must still enumerate every group element exactly once,
        // and `canonicalize` must be the lex-min of that enumeration —
        // the property the incremental walk could plausibly break.
        for g in [
            TranslationGroup::chain_1d(5),
            TranslationGroup::ladder(3, 2),
            TranslationGroup::torus_2d(2, 3),
        ] {
            let n = g.n_qubits();
            let mut w: W = PauliWord::new(n);
            w.set_xbit(0, true);
            w.set_zbit(1, true);
            w.rehash();
            let members: Vec<W> = g.orbit(&w).collect();
            assert_eq!(members.len(), g.order());
            assert_eq!(members[0], w, "identity element must come first");
            // Every member canonicalizes to the same rep = lex-min member.
            let lex_min = *members.iter().min().unwrap();
            assert_eq!(g.canonicalize(&w), lex_min);
            for m in &members {
                assert_eq!(g.canonicalize(m), lex_min);
                // …and the reported shift reproduces the member.
                let (rep, cnt) = g.canonicalize_with_shift(m);
                let mut cur = rep;
                for (gi, &c) in cnt.iter().enumerate() {
                    for _ in 0..c {
                        cur = g.apply_generator(&cur, gi);
                    }
                }
                assert_eq!(cur, *m, "shift {cnt:?} does not reproduce the member");
            }
        }
    }

    #[test]
    fn character_trivial_sector_is_one() {
        let g = TranslationGroup::chain_1d(4);
        // k=0 mode → character is always 1.
        for cnt in [vec![0u32], vec![1u32], vec![2u32], vec![3u32]] {
            let chi = g.character(&[0], &cnt);
            assert!((chi - Complex::new(1.0, 0.0)).norm() < 1e-12);
        }
    }

    #[test]
    fn character_obeys_unit_modulus() {
        let g = TranslationGroup::chain_1d(4);
        for k in 0..4 {
            for a in 0..4 {
                let chi = g.character(&[k], &[a as u32]);
                assert!(
                    (chi.norm() - 1.0).abs() < 1e-12,
                    "|χ_{k}(T^{a})| should be 1, got {}",
                    chi.norm()
                );
            }
        }
    }

    #[test]
    fn momentum_zero_complex_merge_matches_real_merge() {
        // k=0 sector: complex merge with all-real input should give
        // real-valued orbit-rep coefficients equal to the plain
        // canonicalize_pauli_sum result.
        let g = TranslationGroup::chain_1d(4);
        let basis: Vec<W> = vec![word("XIII"), word("IXII"), word("IIXI"), word("IIIX")];
        let real_coeffs = vec![1.0, 2.0, 3.0, 4.0];

        let mut basis_real = basis.clone();
        let mut coeffs_real = real_coeffs.clone();
        canonicalize_pauli_sum(&mut basis_real, &mut coeffs_real, &g);

        let mut basis_c = basis.clone();
        let mut coeffs_c: Vec<Complex<f64>> =
            real_coeffs.iter().map(|&v| Complex::new(v, 0.0)).collect();
        canonicalize_pauli_sum_complex(&mut basis_c, &mut coeffs_c, &g, &[0]);

        // Plain merge sums all coefficients onto the single orbit-rep:
        // 1+2+3+4 = 10. Complex merge does the same with a 1/|G|
        // prefactor, so we expect 10/4 = 2.5 on the rep.
        assert_eq!(basis_real.len(), 1);
        assert_eq!(basis_c.len(), 1);
        assert!((coeffs_real[0] - 10.0).abs() < 1e-12);
        assert!((coeffs_c[0].re - 2.5).abs() < 1e-12);
        assert!(coeffs_c[0].im.abs() < 1e-12);
    }

    #[test]
    fn momentum_eigenstate_check_passes() {
        // O = Σ_j e^{ikj} Z_j for k = 2π/4 (mode 1) is a momentum-k
        // eigenstate. check_momentum_sector should accept.
        let g = TranslationGroup::chain_1d(4);
        let basis: Vec<W> = vec![word("ZIII"), word("IZII"), word("IIZI"), word("IIIZ")];
        let k_mode: i32 = 1;
        // Sector condition: c_{T^a p} = e^{-2πi k a / N} c_p.
        // Picking c_{Z_0} = 1: c_{Z_a} = e^{-2πi · 1 · a / 4} = (-i)^a.
        let coeffs: Vec<Complex<f64>> = (0..4_i32)
            .map(|a| Complex::from_polar(1.0, -2.0 * PI * (k_mode as f64) * (a as f64) / 4.0))
            .collect();
        let res = check_momentum_sector(&basis, &coeffs, &g, &[k_mode], 1e-10);
        assert!(
            res.is_ok(),
            "valid k-eigenstate failed sector check: {res:?}"
        );
    }

    #[test]
    fn momentum_eigenstate_check_fails_for_wrong_sector() {
        // Same eigenstate as above, but check against the wrong momentum.
        let g = TranslationGroup::chain_1d(4);
        let basis: Vec<W> = vec![word("ZIII"), word("IZII"), word("IIZI"), word("IIIZ")];
        let coeffs: Vec<Complex<f64>> = (0..4_i32)
            .map(|a| Complex::from_polar(1.0, -2.0 * PI * 1.0 * (a as f64) / 4.0))
            .collect();
        // Check against k=0 (constant) — should fail.
        let res = check_momentum_sector(&basis, &coeffs, &g, &[0], 1e-10);
        assert!(res.is_err(), "k=1 eigenstate wrongly passed as k=0 sector");
    }

    #[test]
    fn momentum_eigenstate_round_trip_merge_preserves_rep_coeff() {
        // Merge a k=1 eigenstate; the orbit-rep coefficient should be
        // unchanged (= 1.0 for our chosen normalization, picking
        // c_{Z_0} = 1).
        let g = TranslationGroup::chain_1d(4);
        let mut basis: Vec<W> = vec![word("ZIII"), word("IZII"), word("IIZI"), word("IIIZ")];
        let mut coeffs: Vec<Complex<f64>> = (0..4_i32)
            .map(|a| Complex::from_polar(1.0, -2.0 * PI * 1.0 * (a as f64) / 4.0))
            .collect();
        canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &g, &[1]);
        assert_eq!(basis.len(), 1);
        // The canonical rep of single-Z orbit is Z_0 (lex-min of
        // {ZIII, IZII, IIZI, IIIZ} is IIIZ since 'I' < 'Z' lex-wise on
        // the (xbits, zbits) tuple; let's just check we got a single
        // entry with norm 1.
        assert!(
            (coeffs[0].norm() - 1.0).abs() < 1e-10,
            "expected |c_rep|=1, got {}",
            coeffs[0].norm()
        );
    }

    /// Trotter-mode end-to-end check that `PauliSum::symmetry_merge`
    /// matches plain Trotter evolution post-canonicalized.
    ///
    /// Setup: n=4 qubit chain, PBC, XY rotations on each bond. Initial
    /// operator `O(0) = Σ_j Z_j` is translation-invariant.
    ///
    /// **dt must be tiny.** First-order Trotter on a chain with PBC is
    /// only translation-equivariant up to `O(dt^2)` (gate-order
    /// commutator errors are NOT themselves T-symmetric). The
    /// "merge-after-each-step" trajectory and the "merge-at-end"
    /// trajectory therefore diverge by an amount proportional to that
    /// Trotter error. We test in the dt → 0 limit where the divergence
    /// is below FP noise.
    #[test]
    fn pauli_sum_symmetry_merge_matches_plain_trotter() {
        use crate::config::indexmap::ByteFxHashF64;
        use crate::prelude::*;

        type Cfg = ByteFxHashF64<1>;

        let n: usize = 4;
        // Tiny dt — Trotter per-step error scales as dt^2 and shows up
        // as a translation-non-equivariant correction; we want it below
        // FP noise at the tolerance we assert below (1e-7).
        let dt = 1e-5_f64;
        let n_steps = 2usize;
        let group = TranslationGroup::chain_1d(n);

        // Total-Z initial: O(0) = Σ_j Z_j (translation-invariant).
        let mut o_u: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
        let mut o_m: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
        for j in 0..n {
            let mut s: Vec<char> = vec!['I'; n];
            s[j] = 'Z';
            let st: String = s.into_iter().collect();
            o_u += (st.as_str(), 1.0);
            o_m += (st.as_str(), 1.0);
        }
        assert_eq!(o_u.len(), n);
        assert_eq!(o_m.len(), n);

        // Apply XY Trotter steps to both copies. With merging, call
        // symmetry_merge_pauli_sum after each step.
        for _ in 0..n_steps {
            for j in 0..n {
                let nxt = (j + 1) % n;
                o_u.rxx(j, nxt, dt);
                o_u.ryy(j, nxt, dt);
                o_m.rxx(j, nxt, dt);
                o_m.ryy(j, nxt, dt);
            }
            symmetry_merge_pauli_sum(&mut o_m, &group);
        }

        // Canonicalize the un-merged result once at the end.
        symmetry_merge_pauli_sum(&mut o_u, &group);

        // Compare as (word → coeff) maps, FP tolerance.
        let u: FxHashMap<_, f64> = o_u.iter().map(|(w, c)| (*w, *c)).collect();
        let m: FxHashMap<_, f64> = o_m.iter().map(|(w, c)| (*w, *c)).collect();
        assert_eq!(
            u.len(),
            m.len(),
            "post-merge basis sizes differ: u={} vs m={}",
            u.len(),
            m.len()
        );
        let mut max_diff = 0.0_f64;
        for (w, &cu) in &u {
            let cm = *m.get(w).unwrap_or_else(|| {
                panic!("rep present in u but not in m: {:?}", w);
            });
            max_diff = max_diff.max((cu - cm).abs());
        }
        // At dt = 1e-5 over 2 steps, accumulated Trotter
        // commutator-induced T-eq error is ~2·dt^2·|H|^2 ≈ 1e-9; we
        // assert 1e-7 to leave safety margin.
        assert!(
            max_diff < 1e-7,
            "Trotter with-merging diverged from without-merging: max |Δc| = {max_diff:e}"
        );
    }
}
