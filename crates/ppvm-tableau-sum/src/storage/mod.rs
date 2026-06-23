// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

pub mod entry_store;
pub mod map;
pub mod vec;

pub use entry_store::{Branch, EntryStore};
use fxhash::FxHashMap;
use ppvm_traits::traits::Clifford;

// Hasher for the structural `word_fingerprint`. gxhash (AES-based) is fastest on
// native and exposes a `gxhash64` bulk free function, but it needs hardware AES
// and does not build on wasm32, so fall back to fxhash there. The fingerprint is
// a transient in-memory dedup key — collisions are resolved by
// `structurally_equal`, and it is never persisted or compared across builds — so
// the hasher may differ per target without affecting results.
#[cfg(target_arch = "wasm32")]
use fxhash::FxHasher as FingerprintHasher;
use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_tableau::{
    data::GeneralizedTableau, sparsevec::SparseVector, tableau_index::TableauIndex,
};
use ppvm_traits::config::Config;
#[cfg(target_arch = "wasm32")]
use std::hash::Hasher;
use std::ops::AddAssign;

// Reusable per-thread scratch buffer for `word_fingerprint`. Gathering every
// row's word bytes into one contiguous slice lets us hash in a single bulk call
// instead of two tiny `Hash::hash` writes per row (high per-call overhead).
// Cleared (capacity retained) per call, so it adapts to any row count / qubit
// width without re-allocating. Bytes (not the storage word type) because the
// storage element width (`[u8; N]` vs `[u64; N]`) is generic at this call site.
thread_local! {
    static WORD_FP_BUF: std::cell::RefCell<Vec<u8>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// View a `Copy` plain-old-data value's bytes. Sound because `A: PauliStorage`
/// implies `bytemuck::Pod`: no padding, every bit pattern valid, so the bytes
/// are fully initialized and `u8`-aligned.
#[inline]
fn pod_bytes<A: Copy>(value: &A) -> &[u8] {
    // SAFETY: `A` is POD (PauliStorage: bytemuck::Pod); reading its
    // `size_of::<A>()` initialized bytes as `[u8]` is sound, and the borrow is
    // tied to `value`.
    unsafe { std::slice::from_raw_parts(value as *const A as *const u8, std::mem::size_of::<A>()) }
}

/// Hash of the `word` (Pauli content) of every row, in order. This is the
/// expensive component (each word is several machine words wide) and is
/// *invariant* under X/Y/Z and `is_lost` flips, so a branch inherits it from
/// its parent unchanged.
/// NOTE: this inheritance is only valid right now (loss + depolarize channels)
/// but may need re-evaluation in the future when more gates are added
pub fn word_fingerprint<T, I, C>(tab: &GeneralizedTableau<T, I, C>) -> u64
where
    T: Config,
    I:,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    WORD_FP_BUF.with(|cell| {
        let mut buf = cell.borrow_mut();
        // Clear retains capacity; refill with every row's bits as raw bytes.
        buf.clear();
        for row in tab.tableau.data.iter() {
            // Gather the Pauli bits directly: the `PauliWord` hash cache is
            // disabled for tableau rows (`REHASH = false`), so hashing
            // `row.word` would feed a stale zero and make every tableau collide.
            buf.extend_from_slice(pod_bytes(&row.word.xbits.data));
            buf.extend_from_slice(pod_bytes(&row.word.zbits.data));
        }

        #[cfg(not(target_arch = "wasm32"))]
        {
            gxhash::gxhash64(&buf, 0)
        }
        #[cfg(target_arch = "wasm32")]
        {
            let mut hasher = FingerprintHasher::default();
            hasher.write(&buf);
            hasher.finish()
        }
    })
}

/// Per-row mask (splitmix64 of `(index, salt)`); a stable pure function used
/// to build the XOR-combinable [`phase_loss_hash`].
#[inline]
fn row_mask(index: usize, salt: u64) -> u64 {
    let mut z = (index as u64)
        .wrapping_mul(0x9E3779B97F4A7C15)
        .wrapping_add(salt);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// Mask XORed into the phase/loss hash when a row's sign bit (phase bit 1) is
/// set. A Pauli flips only this bit, so a branch updates its hash by XORing
/// `sign_mask(row)` for each row it flips.
#[inline]
pub(crate) fn sign_mask(row: usize) -> u64 {
    row_mask(row, 0xA1A1_A1A1_A1A1_A1A1)
}

/// Mask for a row's imaginary bit (phase bit 0). Stabilizer phases are `±1`,
/// so this is rarely set, but it keeps the hash a faithful function of `phase`.
#[inline]
fn imag_mask(row: usize) -> u64 {
    row_mask(row, 0xB2B2_B2B2_B2B2_B2B2)
}

/// Mask XORed in when qubit `q` is lost. Marking a qubit lost is a single XOR
/// of `loss_mask(q)`.
#[inline]
pub(crate) fn loss_mask(q: usize) -> u64 {
    row_mask(q, 0xC3C3_C3C3_C3C3_C3C3)
}

/// XOR contribution of a single row's phase.
#[inline]
fn phase_contrib(row: usize, phase: u8) -> u64 {
    let mut h = 0;
    if phase & 1 != 0 {
        h ^= imag_mask(row);
    }
    if phase & 2 != 0 {
        h ^= sign_mask(row);
    }
    h
}

/// XOR-combinable hash of `is_lost` plus every row's `phase`, formed as the
/// XOR of per-row contributions (the phase/loss half of [`fingerprint`]). Being
/// XOR-combinable lets a branch inherit its parent's value and update only the
/// rows it changed — a sign flip XORs [`sign_mask`], a loss XORs [`loss_mask`].
pub fn phase_loss_hash<T, I, C>(tab: &GeneralizedTableau<T, I, C>) -> u64
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    let mut h = 0u64;
    for (row, ppw) in tab.tableau.data.iter().enumerate() {
        h ^= phase_contrib(row, ppw.phase);
    }
    for (q, lost) in tab.is_lost.iter().enumerate() {
        if *lost {
            h ^= loss_mask(q);
        }
    }
    h
}

/// Phase/loss hash of a Pauli (depolarize) branch: the parent's hash with
/// [`sign_mask`] XORed in for each row whose phase the Pauli flipped. The
/// branch shares the parent's words and `is_lost`, and a Pauli flips only sign
/// bits, so this single walk over the (already-forked) rows reproduces a
/// from-scratch [`phase_loss_hash`] without re-hashing anything.
pub(crate) fn pauli_branch_phase_loss<T, I, C>(
    parent: &GeneralizedTableau<T, I, C>,
    branch: &GeneralizedTableau<T, I, C>,
    parent_phase_loss: u64,
) -> u64
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    let mut h = parent_phase_loss;
    for (row, (pp, bp)) in parent
        .tableau
        .data
        .iter()
        .zip(branch.tableau.data.iter())
        .enumerate()
    {
        if pp.phase != bp.phase {
            h ^= sign_mask(row);
        }
    }
    h
}

/// Full structural fingerprint: the word component combined with the
/// phase/loss component. Equals `word_fingerprint ^ phase_loss_hash`, which
/// lets branch generation reuse a cached `word_fingerprint` and incrementally
/// update the phase/loss component.
pub(crate) fn fingerprint<T, I, C>(tab: &GeneralizedTableau<T, I, C>) -> u64
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    word_fingerprint(tab) ^ phase_loss_hash(tab)
}

pub(crate) fn structurally_equal<T, I, C>(
    tab0: &GeneralizedTableau<T, I, C>,
    tab1: &GeneralizedTableau<T, I, C>,
    scratch: &mut FxHashMap<I, Complex<T::Coeff>>,
) -> bool
where
    T: Config,
    T::Coeff: One + Zero + Clone + num::Num + PartialOrd,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    // NOTE: comparing is_lost and rows is only necessary to avoid hash collisions

    if tab0.is_lost != tab1.is_lost {
        return false;
    }

    if tab0.coefficients.len() != tab1.coefficients.len() {
        return false;
    }

    // Cheaper row comparison first; coefficient compare is O(K) below.
    for (row0, row1) in tab0.tableau.data.iter().zip(tab1.tableau.data.iter()) {
        if row0.phase != row1.phase || row0.word != row1.word {
            return false;
        }
    }

    // Reuse the caller-owned scratch map instead of allocating per call.
    // Clear retains capacity across invocations.
    scratch.clear();
    scratch.reserve(tab1.coefficients.len());
    for (val, idx) in tab1.coefficients.iter() {
        scratch.insert(*idx, *val);
    }

    let threshold_sq = tab0.coefficient_threshold.clone() * tab0.coefficient_threshold.clone();
    let zero = Complex {
        re: T::Coeff::zero(),
        im: T::Coeff::zero(),
    };
    for (val0, idx0) in tab0.coefficients.iter() {
        let val1 = scratch.get(idx0).copied().unwrap_or(zero);
        if (*val0 - val1).norm_sqr() >= threshold_sq {
            return false;
        }
    }

    true
}

/// A lazily-described branch: a mutation applied to a parent entry. Used so the
/// merge can compute the branch fingerprint / structural identity without
/// cloning, materializing the tableau only for surviving new entries.
#[derive(Clone, Copy, Debug)]
pub enum BranchMutation {
    /// Apply Pauli `op` (1=X, 2=Y, 3=Z) at `addr0`: flips per-row sign bits only.
    Pauli { op: u8, addr0: usize },
    /// Mark qubit `q` lost (set is_lost[q] = true).
    Loss { q: usize },
}

/// Materialize a lazily-described branch into a (cloned) tableau in place.
pub(crate) fn apply_branch_mutation<T, I, C>(
    tab: &mut GeneralizedTableau<T, I, C>,
    m: BranchMutation,
) where
    T: Config,
    I: TableauIndex,
    C: SparseVector<Complex<T::Coeff>, I>,
    GeneralizedTableau<T, I, C>: Clifford,
{
    match m {
        BranchMutation::Pauli { op, addr0 } => match op {
            1 => tab.x(addr0),
            2 => tab.y(addr0),
            3 => tab.z(addr0),
            _ => {}
        },
        BranchMutation::Loss { q } => {
            tab.is_lost[q] = true;
        }
    }
}

/// Like [`structurally_equal`], but compares `existing` against the *virtual*
/// tableau `parent + m` without materializing it. Mirrors `structurally_equal`
/// field-by-field, deriving each field of the virtual tableau from `parent`:
/// - `is_lost`: for `Loss { q }`, equals `parent`'s with index `q` forced true;
///   for `Pauli`, equals `parent`'s unchanged.
/// - `coefficients`: unchanged by both mutations.
/// - rows: for `Loss`, unchanged; for `Pauli`, each row's sign bit (phase bit 1)
///   is flipped per the per-column rule (X: z; Y: x^z; Z: x).
pub(crate) fn structurally_equal_mutated<T, I, C>(
    existing: &GeneralizedTableau<T, I, C>,
    parent: &GeneralizedTableau<T, I, C>,
    m: BranchMutation,
    scratch: &mut FxHashMap<I, Complex<T::Coeff>>,
) -> bool
where
    T: Config,
    T::Coeff: One + Zero + Clone + num::Num + PartialOrd,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    // NOTE: comparing is_lost and rows is only necessary to avoid hash collisions

    match m {
        BranchMutation::Loss { q } => {
            // Virtual is_lost == parent's with index q forced true.
            if existing.is_lost.len() != parent.is_lost.len() {
                return false;
            }
            for (i, (&e, &p)) in existing
                .is_lost
                .iter()
                .zip(parent.is_lost.iter())
                .enumerate()
            {
                let virt = if i == q { true } else { p };
                if e != virt {
                    return false;
                }
            }
        }
        BranchMutation::Pauli { .. } => {
            // Virtual is_lost == parent's, unchanged.
            if existing.is_lost != parent.is_lost {
                return false;
            }
        }
    }

    if existing.coefficients.len() != parent.coefficients.len() {
        return false;
    }

    // Cheaper row comparison first; coefficient compare is O(K) below.
    match m {
        BranchMutation::Loss { .. } => {
            for (re, rp) in existing.tableau.data.iter().zip(parent.tableau.data.iter()) {
                if re.phase != rp.phase || re.word != rp.word {
                    return false;
                }
            }
        }
        BranchMutation::Pauli { op, addr0 } => {
            for (re, rp) in existing.tableau.data.iter().zip(parent.tableau.data.iter()) {
                if re.word != rp.word {
                    return false;
                }
                let x: bool = rp.word.xbits[addr0];
                let z: bool = rp.word.zbits[addr0];
                let flip = match op {
                    1 => z,
                    2 => x ^ z,
                    3 => x,
                    _ => false,
                };
                let virt_phase = rp.phase ^ ((flip as u8) << 1);
                if re.phase != virt_phase {
                    return false;
                }
            }
        }
    }

    // Reuse the caller-owned scratch map instead of allocating per call.
    // Clear retains capacity across invocations. Coefficients are unchanged
    // by both mutations, so compare existing vs parent directly.
    scratch.clear();
    scratch.reserve(parent.coefficients.len());
    for (val, idx) in parent.coefficients.iter() {
        scratch.insert(*idx, *val);
    }

    let threshold_sq =
        existing.coefficient_threshold.clone() * existing.coefficient_threshold.clone();
    let zero = Complex {
        re: T::Coeff::zero(),
        im: T::Coeff::zero(),
    };
    for (val0, idx0) in existing.coefficients.iter() {
        let val1 = scratch.get(idx0).copied().unwrap_or(zero);
        if (*val0 - val1).norm_sqr() >= threshold_sq {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod fingerprint_tests {
    use super::{
        fingerprint, loss_mask, pauli_branch_phase_loss, phase_loss_hash, word_fingerprint,
    };
    use ppvm_pauli_sum::config::fxhash::ByteF64;
    use ppvm_tableau::data::GeneralizedTableau;
    use ppvm_traits::traits::Clifford;

    type Cfg = ByteF64<1>;
    type Tab = GeneralizedTableau<Cfg, u128>;

    fn make() -> Tab {
        let mut t: Tab = GeneralizedTableau::new_with_seed(4, 1e-12, 7);
        t.h(0);
        t.cnot(0, 1);
        t.h(2);
        t
    }

    #[test]
    fn fingerprint_is_word_xor_phase_loss() {
        let t = make();
        assert_eq!(fingerprint(&t), word_fingerprint(&t) ^ phase_loss_hash(&t));
    }

    #[test]
    fn phase_loss_hash_changes_on_sign_flip() {
        // Applying X flips the sign of the rows that anticommute, so the
        // phase/loss hash must change.
        let parent = make();
        let mut branch = parent.clone();
        branch.x(0);
        assert_ne!(phase_loss_hash(&parent), phase_loss_hash(&branch));
    }

    #[test]
    fn phase_loss_hash_sign_flip_delta_matches_recompute() {
        // The core invariant: the incremental branch hash (parent's hash with
        // sign_mask XORed in per flipped row) must equal a from-scratch
        // recompute on the branch. Only holds if phase_loss_hash is
        // XOR-combinable with sign-bit contribution == sign_mask(row).
        let parent = make();
        let mut branch = parent.clone();
        branch.x(0);
        let incremental = pauli_branch_phase_loss(&parent, &branch, phase_loss_hash(&parent));
        assert_eq!(incremental, phase_loss_hash(&branch));
    }

    #[test]
    fn phase_loss_hash_loss_delta_matches_recompute() {
        // Marking a qubit lost must equal XORing loss_mask(q) into the hash.
        let parent = make();
        let mut branch = parent.clone();
        branch.is_lost[1] = true;
        assert_eq!(
            phase_loss_hash(&parent) ^ loss_mask(1),
            phase_loss_hash(&branch)
        );
    }

    #[test]
    fn word_fingerprint_distinguishes_different_words() {
        // H on different qubits produces different Pauli words, so their
        // word-fingerprints must differ. The per-row hash cache is disabled for
        // tableau words, so this only holds if word_fingerprint hashes the bits
        // directly instead of the (stale, zero) cache.
        let mut a: Tab = GeneralizedTableau::new_with_seed(4, 1e-12, 7);
        a.h(0);
        let mut b: Tab = GeneralizedTableau::new_with_seed(4, 1e-12, 7);
        b.h(1);
        assert_ne!(word_fingerprint(&a), word_fingerprint(&b));
    }

    #[test]
    fn pauli_and_loss_preserve_word_fingerprint() {
        // X/Y/Z flip only phase bits and loss flips only is_lost; neither
        // touches `word`. So a branch may inherit its parent's word-hash, and
        // inherited-word XOR fresh-phase_lost must equal a full recompute.
        let parent = make();
        let parent_word = word_fingerprint(&parent);

        for op in 0..3u8 {
            let mut b = parent.clone();
            match op {
                0 => b.x(0),
                1 => b.y(1),
                _ => b.z(2),
            };
            assert_eq!(word_fingerprint(&b), parent_word, "Pauli changed word-hash");
            assert_eq!(
                parent_word ^ phase_loss_hash(&b),
                fingerprint(&b),
                "incremental fingerprint != full recompute after Pauli"
            );
        }

        let mut b = parent.clone();
        b.is_lost[0] = true;
        assert_eq!(word_fingerprint(&b), parent_word, "loss changed word-hash");
        assert_eq!(
            parent_word ^ phase_loss_hash(&b),
            fingerprint(&b),
            "incremental fingerprint != full recompute after loss"
        );
    }
}
