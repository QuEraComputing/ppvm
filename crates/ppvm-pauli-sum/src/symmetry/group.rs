// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::{HashFinalize, PauliStorage, PauliWordTrait};
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
/// Only the **generators** are stored; the algorithm in
/// [`Self::canonicalize`] walks the group via mixed-radix increments.
#[derive(Debug, Clone)]
pub struct TranslationGroup {
    /// Number of qubits the group acts on.
    n_qubits: usize,
    /// One permutation per generator. `perms[g][q]` is the position
    /// that qubit `q` maps to under one application of generator `g`.
    pub(super) perms: Vec<Vec<u32>>,
    /// Cyclic order of each generator.
    pub(super) orders: Vec<u32>,
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
        Self {
            n_qubits,
            perms,
            orders,
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
    pub(super) fn apply_generator<A, S, const R: bool>(
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

    /// Lex-min canonical representative of `w`'s translation orbit
    /// under this group. Walks the full group via mixed-radix counters,
    /// keeping the smallest word seen.
    ///
    /// Total cost: `O(|G| × n_qubits)` per call.
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
        // Mixed-radix counter `(c[0], c[1], …)` ranges over
        // `0..orders[0] × 0..orders[1] × …`. We track the "current"
        // word obtained by applying generator `g` once each time
        // `c[g]` increments; rolling over `c[g]` means we apply
        // generator `g` exactly `orders[g]` times (= identity), so
        // `cur` returns to the orbit member that had `c[g..]` as its
        // tail and `0` in slots 0..g.
        //
        // The simplest correct implementation just enumerates: for each
        // group element index, build the corresponding word from scratch
        // by applying the right number of each generator.
        let mut best = *w;
        let order = self.order();
        let mut idx = 0usize;
        while idx < order {
            // Decode `idx` to mixed-radix counter `c`
            let mut rem = idx;
            let mut counters: Vec<u32> = Vec::with_capacity(self.perms.len());
            for &o in &self.orders {
                counters.push((rem as u32) % o);
                rem /= o as usize;
            }
            // Construct the group element's permutation by composing
            // `generator g` applied `c[g]` times, for each g.
            // We do this lazily by iterating over qubits.
            let mut cur = *w;
            for (g, &c) in counters.iter().enumerate() {
                for _ in 0..c {
                    cur = self.apply_generator(&cur, g);
                }
            }
            if cur < best {
                best = cur;
            }
            idx += 1;
        }
        best
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
    pub fn canonicalize_with_shift<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
    ) -> (PauliWord<A, S, R>, Vec<u32>)
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        debug_assert_eq!(w.n_qubits(), self.n_qubits);
        if self.perms.is_empty() {
            return (*w, Vec::new());
        }
        let mut best = *w;
        let mut best_counter: Vec<u32> = vec![0; self.perms.len()];
        let order = self.order();
        for idx in 0..order {
            // Decode `idx` to mixed-radix counter.
            let mut rem = idx;
            let mut counter: Vec<u32> = Vec::with_capacity(self.perms.len());
            for &o in &self.orders {
                counter.push((rem as u32) % o);
                rem /= o as usize;
            }
            // Build the candidate by applying generator `g` exactly
            // `counter[g]` times.
            let mut cur = *w;
            for (g, &c) in counter.iter().enumerate() {
                for _ in 0..c {
                    cur = self.apply_generator(&cur, g);
                }
            }
            if cur < best {
                best = cur;
                // We need the counter such that g·best = w. The loop
                // above computed cur = g·w with counter, so w = g^{-1}·cur.
                // For abelian cyclic groups, g^{-1} = g^{order-1}, i.e.
                // the counter `(orders[g] - counter[g]) mod orders[g]`.
                best_counter = counter
                    .iter()
                    .zip(self.orders.iter())
                    .map(|(&c, &o)| (o - c) % o)
                    .collect();
            }
        }
        (best, best_counter)
    }

    /// Iterate over all group elements applied to `w`. Yields `|G|`
    /// Pauli words (including `w` itself for the identity element).
    pub fn orbit<'a, A, S, const R: bool>(
        &'a self,
        w: &'a PauliWord<A, S, R>,
    ) -> impl Iterator<Item = PauliWord<A, S, R>> + 'a
    where
        A: PauliStorage + 'a,
        S: BuildHasher + Clone + Default + HashFinalize + 'a,
    {
        let order = self.order();
        (0..order).map(move |idx| {
            let mut rem = idx;
            let mut cur = *w;
            for (g, &o) in self.orders.iter().enumerate() {
                let c = (rem as u32) % o;
                rem /= o as usize;
                for _ in 0..c {
                    cur = self.apply_generator(&cur, g);
                }
            }
            cur
        })
    }
}
