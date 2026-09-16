// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::{HashFinalize, PauliStorage, PauliWordTrait};
use std::hash::BuildHasher;

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn checked_lcm(a: usize, b: usize, context: &str) -> usize {
    a.checked_div(gcd(a, b))
        .and_then(|q| q.checked_mul(b))
        .unwrap_or_else(|| panic!("{context} overflow"))
}

fn permutation_order(perm: &[u32], generator: usize) -> u32 {
    let mut seen = vec![false; perm.len()];
    let mut order = 1usize;
    for start in 0..perm.len() {
        if seen[start] {
            continue;
        }
        let mut length = 0usize;
        let mut q = start;
        loop {
            assert!(!seen[q], "generator {generator} contains a malformed cycle");
            seen[q] = true;
            length += 1;
            q = perm[q] as usize;
            if q == start {
                break;
            }
        }
        order = checked_lcm(order, length, "permutation order");
    }
    u32::try_from(order).unwrap_or_else(|_| {
        panic!("generator {generator} exact permutation order does not fit in u32")
    })
}

fn permutations_commute(left: &[u32], right: &[u32]) -> bool {
    (0..left.len()).all(|q| left[right[q] as usize] == right[left[q] as usize])
}

pub(super) fn checked_group_order(orders: &[u32]) -> usize {
    orders.iter().enumerate().fold(1usize, |acc, (g, &value)| {
        acc.checked_mul(value as usize)
            .unwrap_or_else(|| panic!("group order overflows usize at generator {g}"))
    })
}

pub(super) fn validate_site_count(n: usize, context: &str) {
    let max_index = n
        .checked_sub(1)
        .unwrap_or_else(|| panic!("{context}: site count must be positive"));
    u32::try_from(max_index)
        .unwrap_or_else(|_| panic!("{context}: site count {n} exceeds the u32-addressable range"));
}

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
/// under one application of generator `g`. `orders[g]` is the exact
/// cyclic order of generator `g`. The abstract group is the direct
/// product of these cyclic groups, with order `Π orders[g]`. Its combined
/// permutation action may have a kernel, so distinct group elements can
/// act identically.
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
    order: usize,
    phase_modulus: usize,
}

impl TranslationGroup {
    /// Construct from explicit generator permutations and orders.
    ///
    /// Each `perm` must be a permutation of `0..n_qubits`. Each `order`
    /// must be the permutation's exact cyclic order, not merely a
    /// multiple for which `perm^order == identity`. Generators must
    /// commute, but their combined action may still have a kernel.
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
        for (g, &declared) in orders.iter().enumerate() {
            assert!(declared != 0, "generator {g} order must be nonzero");
            let exact = permutation_order(&perms[g], g);
            assert_eq!(
                declared, exact,
                "generator {g} declared order {declared} != exact permutation order {exact}",
            );
        }
        for left in 0..perms.len() {
            for right in left + 1..perms.len() {
                assert!(
                    permutations_commute(&perms[left], &perms[right]),
                    "generators {left} and {right} do not commute",
                );
            }
        }
        let order = checked_group_order(&orders);
        let phase_modulus = orders.iter().fold(1usize, |acc, &value| {
            checked_lcm(acc, value as usize, "character phase modulus")
        });
        Self {
            n_qubits,
            perms,
            orders,
            order,
            phase_modulus,
        }
    }

    /// 1D chain of `n` sites with periodic boundary conditions.
    /// Single generator: cyclic shift by one site.
    pub fn chain_1d(n: usize) -> Self {
        assert!(n > 0, "chain_1d: n must be positive");
        let order =
            u32::try_from(n).unwrap_or_else(|_| panic!("chain_1d: n={n} does not fit in u32"));
        let perm: Vec<u32> = (0..n)
            .map(|q| {
                u32::try_from((q + 1) % n).expect("chain_1d: target index does not fit in u32")
            })
            .collect();
        Self::from_generators(n, vec![perm], vec![order])
    }

    /// 2D `lx × ly` torus, qubit at `(i, j)` indexed as `j*lx + i`.
    /// Two generators: x-shift (i → i+1 mod lx) and y-shift (j → j+1 mod ly).
    pub fn torus_2d(lx: usize, ly: usize) -> Self {
        assert!(lx > 0, "torus_2d: lx must be positive");
        assert!(ly > 0, "torus_2d: ly must be positive");
        let n = lx
            .checked_mul(ly)
            .unwrap_or_else(|| panic!("torus_2d: lx * ly overflow"));
        validate_site_count(n, "torus_2d");
        let lx_u32 =
            u32::try_from(lx).unwrap_or_else(|_| panic!("torus_2d: lx={lx} does not fit in u32"));
        let ly_u32 =
            u32::try_from(ly).unwrap_or_else(|_| panic!("torus_2d: ly={ly} does not fit in u32"));
        let perm_x: Vec<u32> = (0..n)
            .map(|q| {
                let (i, j) = (q % lx, q / lx);
                u32::try_from(j * lx + (i + 1) % lx)
                    .expect("torus_2d: x-shift target index does not fit in u32")
            })
            .collect();
        let perm_y: Vec<u32> = (0..n)
            .map(|q| {
                let (i, j) = (q % lx, q / lx);
                u32::try_from(((j + 1) % ly) * lx + i)
                    .expect("torus_2d: y-shift target index does not fit in u32")
            })
            .collect();
        Self::from_generators(n, vec![perm_x, perm_y], vec![lx_u32, ly_u32])
    }

    /// 3D `lx × ly × lz` torus, qubit at `(i, j, k)` indexed as
    /// `k*lx*ly + j*lx + i`.
    pub fn torus_3d(lx: usize, ly: usize, lz: usize) -> Self {
        assert!(lx > 0, "torus_3d: lx must be positive");
        assert!(ly > 0, "torus_3d: ly must be positive");
        assert!(lz > 0, "torus_3d: lz must be positive");
        let n = lx
            .checked_mul(ly)
            .and_then(|v| v.checked_mul(lz))
            .unwrap_or_else(|| panic!("torus_3d: lx * ly * lz overflow"));
        validate_site_count(n, "torus_3d");
        let lx_u32 =
            u32::try_from(lx).unwrap_or_else(|_| panic!("torus_3d: lx={lx} does not fit in u32"));
        let ly_u32 =
            u32::try_from(ly).unwrap_or_else(|_| panic!("torus_3d: ly={ly} does not fit in u32"));
        let lz_u32 =
            u32::try_from(lz).unwrap_or_else(|_| panic!("torus_3d: lz={lz} does not fit in u32"));
        let perm_x: Vec<u32> = (0..n)
            .map(|q| {
                let i = q % lx;
                let j = (q / lx) % ly;
                let k = q / (lx * ly);
                u32::try_from(k * lx * ly + j * lx + (i + 1) % lx)
                    .expect("torus_3d: x-shift target index does not fit in u32")
            })
            .collect();
        let perm_y: Vec<u32> = (0..n)
            .map(|q| {
                let i = q % lx;
                let j = (q / lx) % ly;
                let k = q / (lx * ly);
                u32::try_from(k * lx * ly + ((j + 1) % ly) * lx + i)
                    .expect("torus_3d: y-shift target index does not fit in u32")
            })
            .collect();
        let perm_z: Vec<u32> = (0..n)
            .map(|q| {
                let i = q % lx;
                let j = (q / lx) % ly;
                let k = q / (lx * ly);
                u32::try_from(((k + 1) % lz) * lx * ly + j * lx + i)
                    .expect("torus_3d: z-shift target index does not fit in u32")
            })
            .collect();
        Self::from_generators(
            n,
            vec![perm_x, perm_y, perm_z],
            vec![lx_u32, ly_u32, lz_u32],
        )
    }

    /// Multi-leg ladder: `l` sites along the chain × `n_legs` legs.
    /// Single generator: cyclic shift along the chain direction (all
    /// legs simultaneously). Qubit at `(leg, j)` indexed as
    /// `leg * l + j`. No translation along the leg axis (legs are
    /// distinguished).
    pub fn ladder(l: usize, n_legs: usize) -> Self {
        assert!(l > 0, "ladder: l must be positive");
        assert!(n_legs > 0, "ladder: n_legs must be positive");
        let n = l
            .checked_mul(n_legs)
            .unwrap_or_else(|| panic!("ladder: l * n_legs overflow"));
        validate_site_count(n, "ladder");
        let l_u32 =
            u32::try_from(l).unwrap_or_else(|_| panic!("ladder: l={l} does not fit in u32"));
        let perm: Vec<u32> = (0..n)
            .map(|q| {
                let leg = q / l;
                let j = q % l;
                u32::try_from(leg * l + (j + 1) % l)
                    .expect("ladder: shift target index does not fit in u32")
            })
            .collect();
        Self::from_generators(n, vec![perm], vec![l_u32])
    }

    /// Number of qubits the group acts on.
    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Number of generators (rank of the group as an abelian product).
    pub fn n_generators(&self) -> usize {
        self.perms.len()
    }

    /// Abstract product-group order: `Π orders[g]`.
    ///
    /// This can exceed the number of distinct permutations in the action
    /// when the combined action has a kernel.
    pub fn order(&self) -> usize {
        self.order
    }

    /// Permutation associated with the `g`-th generator (one application).
    pub fn generator_perm(&self, g: usize) -> &[u32] {
        &self.perms[g]
    }

    /// Cyclic order of the `g`-th generator.
    pub fn generator_order(&self, g: usize) -> u32 {
        self.orders[g]
    }

    /// Least common multiple of generator orders; denominator for exact
    /// character phase arithmetic.
    pub(super) fn phase_modulus(&self) -> usize {
        self.phase_modulus
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

    pub(super) fn orbit_with_counters<'a, A, S, const R: bool>(
        &'a self,
        word: &'a PauliWord<A, S, R>,
    ) -> GroupOrbit<'a, A, S, R>
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        assert_eq!(
            word.n_qubits(),
            self.n_qubits,
            "word and group must agree on n_qubits"
        );
        GroupOrbit {
            group: self,
            current: *word,
            counter: vec![0; self.orders.len()],
            remaining: self.order,
        }
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
        let mut traversal = self.orbit_with_counters(w);
        let (mut best, _) = traversal
            .next()
            .expect("a finite group contains the identity");
        for (candidate, _) in traversal {
            if candidate < best {
                best = candidate;
            }
        }
        best
    }

    /// Lex-min canonical representative `r` of `w` together with the
    /// **mixed-radix counter** `c = (c_0, c_1, …)` of the group element
    /// `g` such that `g·r = w`.
    ///
    /// In other words: if `r = self.canonicalize(w)`, this returns
    /// `(r, c)` where applying generator `i` exactly `c[i]` times in
    /// sequence to `r` produces `w`. It returns the first valid counter
    /// selected by the deterministic mixed-radix traversal. Counters are
    /// not unique when `r` has a non-trivial stabilizer (or when the
    /// combined action has a kernel). The counter is used to compute
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
        let mut traversal = self.orbit_with_counters(w);
        let (mut best, mut counter_from_word) = traversal
            .next()
            .expect("a finite group contains the identity");
        for (candidate, counter) in traversal {
            if candidate < best {
                best = candidate;
                counter_from_word = counter;
            }
        }
        let counter_to_word = counter_from_word
            .iter()
            .zip(self.orders.iter())
            .map(|(&counter, &order)| (order - counter) % order)
            .collect();
        (best, counter_to_word)
    }

    /// Iterate over all abstract group elements applied to `w`. Yields
    /// [`Self::order`] Pauli words (including `w` itself for the identity
    /// element).
    ///
    /// Words may repeat when `w` has a stabilizer or the combined action
    /// has a kernel; this is not an iterator over distinct orbit members.
    pub fn orbit<'a, A, S, const R: bool>(
        &'a self,
        w: &'a PauliWord<A, S, R>,
    ) -> impl Iterator<Item = PauliWord<A, S, R>> + 'a
    where
        A: PauliStorage + 'a,
        S: BuildHasher + Clone + Default + HashFinalize + 'a,
    {
        self.orbit_with_counters(w).map(|(candidate, _)| candidate)
    }
}

pub(super) struct GroupOrbit<'a, A, S, const R: bool>
where
    A: PauliStorage,
{
    group: &'a TranslationGroup,
    current: PauliWord<A, S, R>,
    counter: Vec<u32>,
    remaining: usize,
}

impl<A, S, const R: bool> Iterator for GroupOrbit<'_, A, S, R>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    type Item = (PauliWord<A, S, R>, Vec<u32>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }
        let item = (self.current, self.counter.clone());
        self.remaining -= 1;
        if self.remaining == 0 {
            return Some(item);
        }
        for g in 0..self.group.orders.len() {
            if self.group.orders[g] == 1 {
                continue;
            }
            self.current = self.group.apply_generator(&self.current, g);
            self.counter[g] += 1;
            if self.counter[g] < self.group.orders[g] {
                break;
            }
            self.counter[g] = 0;
        }
        Some(item)
    }
}
