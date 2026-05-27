pub mod entry_store;
pub mod map;
pub mod vec;

pub use entry_store::EntryStore;
use fxhash::{FxHashMap, FxHasher};
use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::config::Config;
use ppvm_tableau::{
    data::GeneralizedTableau, sparsevec::SparseVector, tableau_index::TableauIndex,
};
use std::{
    hash::{Hash, Hasher},
    ops::AddAssign,
};

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
    let mut hasher = FxHasher::default();
    for row in tab.tableau.data.iter() {
        // Hash the Pauli bits directly: the `PauliWord` hash cache is disabled
        // for tableau rows (`REHASH = false`), so `row.word.hash()` would feed
        // a stale zero and make every tableau collide.
        row.word.xbits.data.hash(&mut hasher);
        row.word.zbits.data.hash(&mut hasher);
    }
    hasher.finish()
}

/// Hash of `is_lost` plus every row's `phase`. Cheap (one byte per row) and
/// the only part that changes under the noise branch ops, so it is recomputed
/// per branch while [`word_fingerprint`] is reused.
pub(crate) fn phase_lost_fingerprint<T, I, C>(tab: &GeneralizedTableau<T, I, C>) -> u64
where
    T: Config,
    I:,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    let mut hasher = FxHasher::default();
    tab.is_lost.hash(&mut hasher);
    for row in tab.tableau.data.iter() {
        row.phase.hash(&mut hasher);
    }
    hasher.finish()
}

/// Full structural fingerprint: the word component combined with the
/// phase/loss component. Equals `word_fingerprint ^ phase_lost_fingerprint`,
/// which lets branch generation reuse a cached `word_fingerprint`.
pub(crate) fn fingerprint<T, I, C>(tab: &GeneralizedTableau<T, I, C>) -> u64
where
    T: Config,
    I:,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    word_fingerprint(tab) ^ phase_lost_fingerprint(tab)
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

#[cfg(test)]
mod fingerprint_tests {
    use super::{fingerprint, phase_lost_fingerprint, word_fingerprint};
    use ppvm_runtime::config::fxhash::ByteF64;
    use ppvm_runtime::traits::Clifford;
    use ppvm_tableau::data::GeneralizedTableau;

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
    fn fingerprint_is_word_xor_phase_lost() {
        let t = make();
        assert_eq!(
            fingerprint(&t),
            word_fingerprint(&t) ^ phase_lost_fingerprint(&t)
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
                parent_word ^ phase_lost_fingerprint(&b),
                fingerprint(&b),
                "incremental fingerprint != full recompute after Pauli"
            );
        }

        let mut b = parent.clone();
        b.is_lost[0] = true;
        assert_eq!(word_fingerprint(&b), parent_word, "loss changed word-hash");
        assert_eq!(
            parent_word ^ phase_lost_fingerprint(&b),
            fingerprint(&b),
            "incremental fingerprint != full recompute after loss"
        );
    }
}
