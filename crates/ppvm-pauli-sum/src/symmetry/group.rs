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

fn checked_lcm(a: usize, b: usize) -> Option<usize> {
    a.checked_div(gcd(a, b)).and_then(|q| q.checked_mul(b))
}

/// Exact cyclic order of `perm`, or `None` if it does not fit in `u32`.
fn permutation_order(perm: &[u32], generator: usize) -> Option<u32> {
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
        order = checked_lcm(order, length)?;
    }
    u32::try_from(order).ok()
}

fn permutations_commute(left: &[u32], right: &[u32]) -> bool {
    (0..left.len()).all(|q| left[right[q] as usize] == right[left[q] as usize])
}

pub(super) fn checked_group_order(orders: &[u32]) -> Option<usize> {
    orders
        .iter()
        .try_fold(1usize, |acc, &value| acc.checked_mul(value as usize))
}

pub(super) fn validate_site_count(n: usize, context: &str) {
    let max_index = n
        .checked_sub(1)
        .unwrap_or_else(|| panic!("{context}: site count must be positive"));
    u32::try_from(max_index)
        .unwrap_or_else(|_| panic!("{context}: site count {n} exceeds the u32-addressable range"));
}

/// A precondition violation in [`TranslationGroup::try_from_generators`].
///
/// Every variant is caller-supplied-input error, and its [`Display`]
/// text is exactly what [`TranslationGroup::from_generators`] panics
/// with.
///
/// [`Display`]: std::fmt::Display
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GroupError {
    /// `perms` and `orders` describe different numbers of generators.
    LengthMismatch { perms: usize, orders: usize },
    /// A generator's permutation is not `n_qubits` long.
    PermutationLength {
        generator: usize,
        len: usize,
        n_qubits: usize,
    },
    /// A generator maps a qubit outside `0..n_qubits`.
    TargetOutOfRange {
        generator: usize,
        target: u32,
        n_qubits: usize,
    },
    /// A generator maps two qubits to the same position.
    DuplicateTarget { generator: usize, target: u32 },
    /// A generator declares cyclic order zero.
    ZeroOrder { generator: usize },
    /// A generator's declared order is not its exact cyclic order.
    OrderMismatch {
        generator: usize,
        declared: u32,
        exact: u32,
    },
    /// Two generators do not commute, so they generate no abelian group.
    NonCommuting { left: usize, right: usize },
    /// A generator's exact cyclic order does not fit in `u32`.
    PermutationOrderOverflow { generator: usize },
    /// The group order `Π orders[g]` does not fit in `usize`.
    GroupOrderOverflow,
}

impl std::fmt::Display for GroupError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LengthMismatch { perms, orders } => write!(
                f,
                "perms ({perms} generators) and orders ({orders}) must have the same length"
            ),
            Self::PermutationLength {
                generator,
                len,
                n_qubits,
            } => write!(
                f,
                "generator {generator}: permutation length {len} != n_qubits {n_qubits}"
            ),
            Self::TargetOutOfRange {
                generator,
                target,
                n_qubits,
            } => write!(
                f,
                "generator {generator}: target {target} out of range [0, {n_qubits})"
            ),
            Self::DuplicateTarget { generator, target } => write!(
                f,
                "generator {generator}: not a permutation (duplicate target {target})"
            ),
            Self::ZeroOrder { generator } => {
                write!(f, "generator {generator} order must be nonzero")
            }
            Self::OrderMismatch {
                generator,
                declared,
                exact,
            } => write!(
                f,
                "generator {generator} declared order {declared} != exact permutation order {exact}"
            ),
            Self::NonCommuting { left, right } => {
                write!(f, "generators {left} and {right} do not commute")
            }
            Self::PermutationOrderOverflow { generator } => write!(
                f,
                "generator {generator} exact permutation order does not fit in u32"
            ),
            Self::GroupOrderOverflow => write!(f, "group order overflows usize"),
        }
    }
}

impl std::error::Error for GroupError {}

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
/// Only the **generators** are stored; [`Self::canonicalize`] either runs
/// the `O(N)` least-rotation scan (chain/ladder layouts) or walks the
/// group as a mixed-radix odometer.
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
    /// Set when the group is a *single* generator acting as a cyclic
    /// shift inside contiguous, aligned blocks of qubits — i.e. exactly
    /// the [`Self::chain_1d`] and [`Self::ladder`] layouts. Enables the
    /// `O(N)` least-rotation canonicalizer (see
    /// [`Self::canonicalize_block_cyclic`]).
    pub(super) block_cyclic: Option<BlockCyclic>,
    /// Per generator, its block-rotation form when it has one (all lattice
    /// translations do). Enables the masked-shift [`Self::apply_generator`].
    pub(super) rotations: Vec<Option<BlockRotation>>,
}

/// Layout of a single-generator group acting as a cyclic shift within
/// `n_blocks` contiguous, aligned blocks of `len` qubits each: qubit
/// `b * len + j` maps to `b * len + (j + 1) % len`.
#[derive(Debug, Clone, Copy)]
pub(super) struct BlockCyclic {
    n_blocks: usize,
    len: usize,
}

/// A generator that acts as a cyclic shift by `stride` positions within
/// aligned blocks of `block` qubits: `b·block + p ↦ b·block + (p + stride) mod block`.
///
/// Every lattice-translation generator has this form: the fastest axis of a
/// torus is `stride = 1` with `block = lx`, the next is `stride = lx` with
/// `block = lx·ly`, and so on. Recognising it lets the whole permutation be
/// applied as a masked shift of the two bit planes rather than a per-qubit
/// gather (see [`TranslationGroup::apply_block_rotation`]).
#[derive(Debug, Clone)]
pub(super) struct BlockRotation {
    stride: usize,
    block: usize,
    /// Destinations that survive the plain left shift: everything except the
    /// low `stride` slots of each block (which receive the previous block's
    /// spill) and everything at or beyond `n_qubits`.
    keep: Vec<u64>,
    /// Sources that wrap: the top `stride` slots of each block.
    high: Vec<u64>,
}

/// Widest storage the masked-shift path handles, in 64-bit words.
const MAX_ROT_WORDS: usize = 16;

/// Recognise a generator permutation as a [`BlockRotation`], and precompute
/// its masks. Returns `None` for permutations that are not block rotations.
fn detect_block_rotation(n_qubits: usize, perm: &[u32]) -> Option<BlockRotation> {
    if n_qubits == 0 {
        return None;
    }
    let stride = perm[0] as usize;
    if stride == 0 {
        return None;
    }
    // Whatever maps to qubit 0 sits `stride` below the top of block 0.
    let block = perm.iter().position(|&t| t == 0)? + stride;
    if block > n_qubits || stride >= block || !n_qubits.is_multiple_of(block) {
        return None;
    }
    for b in 0..n_qubits / block {
        for p in 0..block {
            if perm[b * block + p] as usize != b * block + (p + stride) % block {
                return None;
            }
        }
    }
    let mut keep = vec![0u64; n_qubits.div_ceil(64)];
    let mut high = keep.clone();
    for q in 0..n_qubits {
        let p = q % block;
        if p >= stride {
            keep[q / 64] |= 1u64 << (q % 64);
        }
        if p >= block - stride {
            high[q / 64] |= 1u64 << (q % 64);
        }
    }
    Some(BlockRotation {
        stride,
        block,
        keep,
        high,
    })
}

/// `dst = src << s` over a little-endian multiword bit array.
#[inline]
fn shl_words(src: &[u64], dst: &mut [u64], s: usize) {
    let (ws, bs) = (s / 64, s % 64);
    for i in (0..src.len()).rev() {
        let lo = if i >= ws { src[i - ws] } else { 0 };
        dst[i] = if bs == 0 {
            lo
        } else {
            let hi = if i > ws {
                src[i - ws - 1] >> (64 - bs)
            } else {
                0
            };
            (lo << bs) | hi
        };
    }
}

/// `dst = src >> s` over a little-endian multiword bit array.
#[inline]
fn shr_words(src: &[u64], dst: &mut [u64], s: usize) {
    let n = src.len();
    let (ws, bs) = (s / 64, s % 64);
    for i in 0..n {
        let hi = if i + ws < n { src[i + ws] } else { 0 };
        dst[i] = if bs == 0 {
            hi
        } else {
            let lo = if i + ws + 1 < n {
                src[i + ws + 1] << (64 - bs)
            } else {
                0
            };
            (hi >> bs) | lo
        };
    }
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
    /// Construct from explicit generator permutations and orders,
    /// panicking on any precondition violation.
    ///
    /// Each `perm` must be a permutation of `0..n_qubits`. Each `order`
    /// must be the permutation's exact cyclic order, not merely a
    /// multiple for which `perm^order == identity`. Generators must
    /// commute, but their combined action may still have a kernel.
    ///
    /// Use [`Self::try_from_generators`] when the generators come from
    /// outside the program (an FFI boundary, a config file) and a
    /// precondition violation should be reported rather than abort.
    pub fn from_generators(n_qubits: usize, perms: Vec<Vec<u32>>, orders: Vec<u32>) -> Self {
        Self::try_from_generators(n_qubits, perms, orders)
            .unwrap_or_else(|err| panic!("TranslationGroup::from_generators: {err}"))
    }

    /// Fallible [`Self::from_generators`]: validates every precondition
    /// on the caller-supplied generators and reports the first violation
    /// as a [`GroupError`] instead of panicking.
    pub fn try_from_generators(
        n_qubits: usize,
        perms: Vec<Vec<u32>>,
        orders: Vec<u32>,
    ) -> Result<Self, GroupError> {
        if perms.len() != orders.len() {
            return Err(GroupError::LengthMismatch {
                perms: perms.len(),
                orders: orders.len(),
            });
        }
        for (generator, perm) in perms.iter().enumerate() {
            if perm.len() != n_qubits {
                return Err(GroupError::PermutationLength {
                    generator,
                    len: perm.len(),
                    n_qubits,
                });
            }
            let mut seen = vec![false; n_qubits];
            for &target in perm {
                if target as usize >= n_qubits {
                    return Err(GroupError::TargetOutOfRange {
                        generator,
                        target,
                        n_qubits,
                    });
                }
                if seen[target as usize] {
                    return Err(GroupError::DuplicateTarget { generator, target });
                }
                seen[target as usize] = true;
            }
        }
        for (generator, &declared) in orders.iter().enumerate() {
            if declared == 0 {
                return Err(GroupError::ZeroOrder { generator });
            }
            let exact = permutation_order(&perms[generator], generator)
                .ok_or(GroupError::PermutationOrderOverflow { generator })?;
            if declared != exact {
                return Err(GroupError::OrderMismatch {
                    generator,
                    declared,
                    exact,
                });
            }
        }
        for left in 0..perms.len() {
            for right in left + 1..perms.len() {
                if !permutations_commute(&perms[left], &perms[right]) {
                    return Err(GroupError::NonCommuting { left, right });
                }
            }
        }
        let order = checked_group_order(&orders).ok_or(GroupError::GroupOrderOverflow)?;
        // lcm divides the product, so this cannot overflow once `order` fits.
        let phase_modulus = orders
            .iter()
            .try_fold(1usize, |acc, &value| checked_lcm(acc, value as usize))
            .ok_or(GroupError::GroupOrderOverflow)?;
        let block_cyclic = detect_block_cyclic(n_qubits, &perms, &orders);
        let rotations = perms
            .iter()
            .map(|p| detect_block_rotation(n_qubits, p))
            .collect();
        Ok(Self {
            n_qubits,
            perms,
            orders,
            order,
            phase_modulus,
            block_cyclic,
            rotations,
        })
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

    /// Apply a single generator's permutation to a Pauli word: for each
    /// qubit `q` of the input, the `(xbit, zbit)` pair is placed at position
    /// `perm[q]` of the output.
    ///
    /// Does **not** refresh the cached hash. Equality and ordering compare
    /// the bit planes, so an unhashed word is safe to compare and to keep as
    /// an intermediate; only words that escape into a hash container need
    /// `rehash`. The odometer walk applies a generator per group element and
    /// hashes just the winner.
    pub(super) fn apply_generator<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
        g: usize,
    ) -> PauliWord<A, S, R>
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        if let Some(rot) = &self.rotations[g]
            && let Some(out) = Self::apply_block_rotation(w, rot)
        {
            return out;
        }
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
        out
    }

    /// Apply a block-rotation generator as a masked shift of both bit
    /// planes: `out = ((in << stride) & keep) | ((in & high) >> (block − stride))`.
    ///
    /// This is the same permutation as the per-qubit gather, in `O(N/64)`
    /// word operations instead of `O(N)` bit operations. The cached hash is
    /// *not* refreshed. Returns `None` on big-endian targets (where the byte
    /// view of the bit planes is not in bit order) or if the storage is wider
    /// than [`MAX_ROT_WORDS`], leaving the caller on the general path.
    pub(super) fn apply_block_rotation<A, S, const R: bool>(
        w: &PauliWord<A, S, R>,
        rot: &BlockRotation,
    ) -> Option<PauliWord<A, S, R>>
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        if !cfg!(target_endian = "little") || size_of::<A>() > MAX_ROT_WORDS * 8 {
            return None;
        }
        let mut out = *w;
        for plane in 0..2 {
            let (src_arr, dst_arr) = if plane == 0 {
                (&w.xbits.data, &mut out.xbits.data)
            } else {
                (&w.zbits.data, &mut out.zbits.data)
            };
            let bytes = bytemuck::bytes_of(src_arr);
            let nw = bytes.len().div_ceil(8);
            let (mut src, mut shifted, mut wrapped) = (
                [0u64; MAX_ROT_WORDS],
                [0u64; MAX_ROT_WORDS],
                [0u64; MAX_ROT_WORDS],
            );
            for (i, chunk) in bytes.chunks(8).enumerate() {
                let mut b = [0u8; 8];
                b[..chunk.len()].copy_from_slice(chunk);
                src[i] = u64::from_le_bytes(b);
            }
            shl_words(&src[..nw], &mut shifted[..nw], rot.stride);
            for (i, s) in src[..nw].iter_mut().enumerate() {
                *s &= rot.high.get(i).copied().unwrap_or(0);
            }
            shr_words(&src[..nw], &mut wrapped[..nw], rot.block - rot.stride);
            for i in 0..nw {
                shifted[i] = (shifted[i] & rot.keep.get(i).copied().unwrap_or(0)) | wrapped[i];
            }
            let dst = bytemuck::bytes_of_mut(dst_arr);
            for (i, chunk) in dst.chunks_mut(8).enumerate() {
                let b = shifted[i].to_le_bytes();
                let n = chunk.len();
                chunk.copy_from_slice(&b[..n]);
            }
        }
        Some(out)
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
    /// for a full walk. Leaves the cached hash of `cur` stale.
    #[inline]
    fn advance<A, S, const R: bool>(&self, cur: &mut PauliWord<A, S, R>, idx: usize)
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        let mut p = 1usize;
        for (g, &o) in self.orders.iter().enumerate() {
            if o > 1 {
                *cur = self.apply_generator(cur, g);
            }
            p *= o as usize;
            if !idx.is_multiple_of(p) {
                break;
            }
        }
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
        self.canonicalize_with_index(w).0
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
    /// Same cost as [`Self::canonicalize`], plus the counter `Vec`. Hot
    /// paths should prefer [`Self::canonicalize_with_index`], which is
    /// allocation-free and indexes a precomputed
    /// [`Self::character_table`](crate::symmetry::TranslationGroup::character_table).
    pub fn canonicalize_with_shift<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
    ) -> (PauliWord<A, S, R>, Vec<u32>)
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        let (rep, idx) = self.canonicalize_with_index(w);
        (rep, self.counter_from_index(idx))
    }

    /// Lex-min canonical representative `r` of `w` together with the
    /// **mixed-radix index** (generator `0` fastest) of the group element
    /// `g` such that `g·r = w` — i.e. the index of the counter returned by
    /// [`Self::canonicalize_with_shift`].
    ///
    /// The index is directly usable as a subscript into a
    /// [`CharacterTable`](crate::symmetry::CharacterTable), which is how the
    /// phase-aware evolution gets `χ_k(g)` without decoding a counter or
    /// calling `sin`/`cos` per term.
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
        assert_eq!(
            w.n_qubits(),
            self.n_qubits,
            "word and group must agree on n_qubits"
        );
        match self.block_cyclic {
            Some(bc) => {
                let (rep, r, _) = self.canonicalize_block_cyclic(w, bc);
                (rep, r)
            }
            None => {
                let (rep, idx, _) = self.canonicalize_odometer(w, |_| true);
                (rep, idx)
            }
        }
    }

    /// Canonical rep, the index of the group element mapping it back to `w`,
    /// and the **stabilizer** of `w` checked against `trivial`: returns
    /// `None` as soon as a stabilizer element `s` (`s·w = w`) with
    /// `!trivial(index(s))` is found, else `Some((rep, index, |stabilizer|))`.
    ///
    /// One traversal gives everything the momentum-sector routines need.
    /// `trivial` is only consulted on stabilizer elements, which are rare
    /// (none but the identity for a free orbit).
    pub(super) fn canonicalize_with_stabilizer<A, S, const R: bool, F>(
        &self,
        w: &PauliWord<A, S, R>,
        trivial: F,
    ) -> Option<(PauliWord<A, S, R>, usize, usize)>
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
        F: Fn(usize) -> bool,
    {
        assert_eq!(
            w.n_qubits(),
            self.n_qubits,
            "word and group must agree on n_qubits"
        );
        match self.block_cyclic {
            Some(bc) => {
                let (rep, r, step) = self.canonicalize_block_cyclic(w, bc);
                // The stabilizer of a single-generator group is the cyclic
                // subgroup generated by `g^step`; its characters are all
                // trivial iff that generator's is.
                if step < bc.len && !trivial(step) {
                    return None;
                }
                Some((rep, r, bc.len / step))
            }
            None => {
                let (rep, idx, stabilizer) = self.canonicalize_odometer(w, trivial);
                (stabilizer != 0).then_some((rep, idx, stabilizer))
            }
        }
    }

    /// Reference canonicalizer: walk the whole group as a mixed-radix
    /// odometer (see [`Self::advance`]), keeping the first smallest word
    /// seen. Returns the rep, the index of the group element mapping it
    /// back to `w`, and the stabilizer size — or `0` for the latter if a
    /// stabilizer element fails `trivial` (the walk stops there).
    ///
    /// `O(|G| × N)`; used for groups without a [`BlockCyclic`] layout, and
    /// as the test oracle for the fast path.
    pub(super) fn canonicalize_odometer<A, S, const R: bool, F>(
        &self,
        w: &PauliWord<A, S, R>,
        trivial: F,
    ) -> (PauliWord<A, S, R>, usize, usize)
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
        F: Fn(usize) -> bool,
    {
        let mut best = *w;
        let mut best_idx = 0usize;
        let mut stabilizer = 1usize;
        let mut cur = *w;
        for idx in 1..self.order {
            self.advance(&mut cur, idx);
            if cur == *w {
                if !trivial(idx) {
                    return (best, 0, 0);
                }
                stabilizer += 1;
            }
            if cur < best {
                best = cur;
                best_idx = idx;
            }
        }
        // `advance` leaves the cached hash stale; the winner escapes to the
        // caller (and into hash containers), so refresh it here.
        best.rehash();
        // The walk found `best = g·w` at index `best_idx`, so `w = g⁻¹·best`
        // and the element we must report is the inverse.
        (best, self.invert_index(best_idx), stabilizer)
    }

    /// `O(N)` canonicalizer for single-generator cyclic-block groups
    /// (chain, ladder): returns the same rep as the odometer walk — the
    /// `Ord`-lex-min of the orbit — the index `r` of the group element
    /// with `g^r · rep = w` (the odometer's choice), and the smallest
    /// `step > 0` with `g^step · w = w` (`step == len` for a free orbit).
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
    /// The final survivors are one coset of the stabilizer, so `step` is
    /// its generator.
    pub(super) fn canonicalize_block_cyclic<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
        bc: BlockCyclic,
    ) -> (PauliWord<A, S, R>, usize, usize)
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
        (rep, r, step)
    }

    /// Decode a group-element index (mixed-radix, generator `0` fastest)
    /// into its per-generator counter.
    pub fn counter_from_index(&self, idx: usize) -> Vec<u32> {
        let mut rem = idx;
        let mut counter: Vec<u32> = Vec::with_capacity(self.orders.len());
        for &o in &self.orders {
            counter.push((rem % o as usize) as u32);
            rem /= o as usize;
        }
        counter
    }

    /// Index of the inverse of the group element with index `idx`. In an
    /// abelian product of cyclic groups that is `(orders[g] − c[g]) mod
    /// orders[g]` componentwise.
    pub(super) fn invert_index(&self, idx: usize) -> usize {
        let mut rem = idx;
        let mut out = 0usize;
        let mut stride = 1usize;
        for &o in &self.orders {
            let o = o as usize;
            let c = rem % o;
            rem /= o;
            out += ((o - c) % o) * stride;
            stride *= o;
        }
        out
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
        // `apply_generator` skips hashing; yielded words may become hash
        // keys, so hash each one on the way out.
        let mut word = self.current;
        word.rehash();
        let item = (word, self.counter.clone());
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
