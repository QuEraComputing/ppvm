// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The graded map algebra `impl`'d **directly on the containers** — no wrapper
//! types. Two backends:
//!
//! - [`Vec<(K, C)>`] — an unsorted coordinate list with linear-scan
//!   `accumulate`; best for small support (the `GeneralizedTableau` amplitude
//!   vector). Requires only `K: Eq + Clone` — it never hashes.
//! - [`HashMap<K, C, IdentityBuildHasher>`] — the hash-join `accumulate`; best
//!   for large support (`PauliSum`). Requires `K: Indexable` (the direct digest,
//!   consumed pass-through through [`IdentityBuildHasher`]).
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Backends are containers;
//! columnar is expressible from day one" and §"The map is a graded algebra over
//! `C[K]`". The module laws these impls realize are machine-checked in
//! `lean/PPVM/Algebra/GradedMap.lean` (`accumulate_comm`/`accumulate_assoc`,
//! `reduce_structural`, `scale_scale`, `overlap_comm`).
//!
//! # Friction: these impls live in `ppvm-traits-2`, not `ppvm-pauli-sum-2`
//!
//! The implementation plan's crate-map table lists "graded traits `impl`'d on
//! `Vec`/`HashMap`" under `ppvm-pauli-sum-2`. Rust's orphan rule forbids that:
//! both the traits ([`Support`] …) and the containers (`Vec`, `HashMap`) would be
//! foreign to `ppvm-pauli-sum-2`, and `(K, C)` carries no local type to anchor
//! the impl. The impls must therefore live in the crate that *owns* the graded
//! traits — here. `ppvm-pauli-sum-2` still owns everything the orphan rule
//! permits it to (`Sum`, `Policy`, the producers, `Clifford for Sum`, the
//! aliases); only these container impls move.

use std::collections::HashMap;

use crate::algebra::{Conjugate, ImaginaryUnit, KeyProduct};
use crate::batch::{KeyBatch, TermBatch};
use crate::coefficient::Coefficient;
use crate::graded::{Accumulate, Multiply, Pair, Retain, Scale, Support};
use crate::hash::{IdentityBuildHasher, Indexable};

// ===========================================================================
// Vec<(K, C)> — coordinate-list backend (linear scan; K: Eq + Clone).
// ===========================================================================

impl<K, C> Support for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    type Key = K;
    type Coeff = C;

    #[inline]
    fn len(&self) -> usize {
        self.as_slice().len()
    }

    #[inline]
    fn get(&self, key: &K) -> Option<C> {
        self.as_slice()
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, c)| c.clone())
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = (K, C)> {
        self.as_slice().iter().map(|(k, c)| (k.clone(), c.clone()))
    }

    /// The borrowing scan: no clone at all, so a reader that rejects most terms
    /// pays nothing for the ones it rejects.
    #[inline]
    fn for_each_ref(&self, mut f: impl FnMut(&K, &C)) {
        for (k, c) in self.as_slice() {
            f(k, c);
        }
    }
}

impl<K, C> Accumulate for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    /// Linear-scan hash-join: for each produced term, find the matching key and
    /// add onto it, else append. `O(n·m)` in the support size, which is the
    /// right cost model for the small support this backend targets.
    #[inline]
    fn accumulate_batch(&mut self, terms: &TermBatch<K, C>) {
        for (k, c) in terms.iter() {
            if let Some(slot) = self.iter_mut().find(|(ek, _)| ek == k) {
                slot.1 += c.clone();
            } else {
                self.push((k.clone(), c.clone()));
            }
        }
    }

    /// Drop every zero-coefficient term (`reduce_structural`): canonicalize to
    /// reduced finite support. Runs only at finalize, never inline.
    #[inline]
    fn reduce(&mut self) {
        self.retain(|(_, v)| !v.is_zero());
    }
}

impl<K, C> Scale for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    #[inline]
    fn scale(&mut self, s: &C) {
        for (_, v) in self.iter_mut() {
            *v *= s.clone();
        }
    }
}

impl<K, C> Pair for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    #[inline]
    fn probe_batch(&self, keys: &KeyBatch<K>, out: &mut [Option<C>]) {
        debug_assert!(out.len() >= keys.keys().len());
        for (slot, k) in out.iter_mut().zip(keys.keys().iter()) {
            *slot = Support::get(self, k);
        }
    }

    #[inline]
    fn overlap(&self, other: &Self) -> C {
        self.as_slice()
            .iter()
            .filter_map(|(k, a)| Support::get(other, k).map(|b| a.clone() * b))
            .sum()
    }

    #[inline]
    fn hermitian_overlap(&self, other: &Self) -> C
    where
        C: Conjugate,
    {
        self.as_slice()
            .iter()
            .filter_map(|(k, a)| Support::get(other, k).map(|b| a.conj() * b))
            .sum()
    }
}

impl<K, C> Retain<K, C> for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    #[inline]
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        // Inherent `Vec::retain` shadows the trait method (inherent-first
        // resolution), so this does not recurse.
        self.retain(|(k, v)| keep(k, v));
    }
}

impl<K, C> Multiply for Vec<(K, C)>
where
    K: KeyProduct,
    C: ImaginaryUnit,
{
    /// The twisted convolution `(A·B)[k] = Σ_{p·q = k} A[p]·B[q]·i^{β(p,q)}`,
    /// accumulated into `acc` — the coordinate-list spelling of `twistedConv`
    /// (`lean/PPVM/Algebra/Twisted.lean`), whose monomial case is `tmul`.
    ///
    /// Neither `reduce` nor any truncation runs: `acc` keeps an exact-zero
    /// cancellation, exactly as `twistedConv` (a finitely-supported map is
    /// canonicalized only by an explicit [`Accumulate::reduce`]).
    fn multiply_into(&self, other: &Self, acc: &mut Self) {
        for (p, a) in self.as_slice() {
            for (q, b) in other.as_slice() {
                let (k, phase) = p.key_mul(q);
                let c = phase.apply(&(a.clone() * b.clone()));
                if let Some(slot) = acc.iter_mut().find(|(ek, _)| *ek == k) {
                    slot.1 += c;
                } else {
                    acc.push((k, c));
                }
            }
        }
    }
}

// ===========================================================================
// HashMap<K, C, IdentityBuildHasher> — hash-join backend (K: Indexable).
// ===========================================================================

impl<K, C> Support for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    type Key = K;
    type Coeff = C;

    #[inline]
    fn len(&self) -> usize {
        HashMap::len(self)
    }

    #[inline]
    fn get(&self, key: &K) -> Option<C> {
        HashMap::get(self, key).cloned()
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = (K, C)> {
        HashMap::iter(self).map(|(k, v)| (k.clone(), v.clone()))
    }

    /// The borrowing scan: hands out `(&K, &C)` straight from the buckets, so a
    /// filtering reader never clones a coefficient it is about to reject. Same
    /// order as [`Support::iter`] (the map's own bucket order).
    #[inline]
    fn for_each_ref(&self, mut f: impl FnMut(&K, &C)) {
        for (k, v) in HashMap::iter(self) {
            f(k, v);
        }
    }
}

impl<K, C> Accumulate for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    /// Build side of the hash join: probe each produced term, accumulate its
    /// coefficient onto a matching key, insert on a miss.
    #[inline]
    fn accumulate_batch(&mut self, terms: &TermBatch<K, C>) {
        for (k, c) in terms.iter() {
            self.entry(k.clone())
                .and_modify(|e| *e += c.clone())
                .or_insert_with(|| c.clone());
        }
    }

    /// Drop every zero-coefficient key (`reduce_structural`).
    #[inline]
    fn reduce(&mut self) {
        self.retain(|_, v| !v.is_zero());
    }
}

impl<K, C> Scale for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn scale(&mut self, s: &C) {
        for v in self.values_mut() {
            *v *= s.clone();
        }
    }
}

impl<K, C> Pair for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn probe_batch(&self, keys: &KeyBatch<K>, out: &mut [Option<C>]) {
        debug_assert!(out.len() >= keys.keys().len());
        for (slot, k) in out.iter_mut().zip(keys.keys().iter()) {
            *slot = HashMap::get(self, k).cloned();
        }
    }

    /// `Σ_k self[k]·other[k]`, driven from the **smaller** support.
    ///
    /// The shared support is contained in both, so scanning either side and
    /// probing the other is `O(1)` per candidate and yields the same pair set;
    /// walking the smaller one makes the cost `O(min(|a|, |b|))` instead of
    /// `O(|self|)`. Against old's `data().trace(k)`-per-term
    /// (`ppvm-pauli-sum/src/sum/trace.rs`, a full linear scan of `self` per term
    /// of `other` — anti-feature 13) this is the intended asymptotic improvement;
    /// picking the smaller side is the second half of it, and it matters when the
    /// operands are deliberately unequal (`|a| = 10⁵`, `|b| = 10`).
    ///
    /// The *value* is unchanged — the pairing is symmetric and the left factor
    /// stays `self`'s coefficient either way — only the float **summation order**
    /// depends on the direction, which is why the differential bar on `overlap` is
    /// relative (`1e-12`) rather than bit-exact.
    #[inline]
    fn overlap(&self, other: &Self) -> C {
        if self.len() <= other.len() {
            HashMap::iter(self)
                .filter_map(|(k, a)| HashMap::get(other, k).map(|b| a.clone() * b.clone()))
                .sum()
        } else {
            HashMap::iter(other)
                .filter_map(|(k, b)| HashMap::get(self, k).map(|a| a.clone() * b.clone()))
                .sum()
        }
    }

    /// `Σ_k conj(self[k])·other[k]`, driven from the smaller support — see
    /// [`Pair::overlap`] for why the direction is free to differ.
    #[inline]
    fn hermitian_overlap(&self, other: &Self) -> C
    where
        C: Conjugate,
    {
        if self.len() <= other.len() {
            HashMap::iter(self)
                .filter_map(|(k, a)| HashMap::get(other, k).map(|b| a.conj() * b.clone()))
                .sum()
        } else {
            HashMap::iter(other)
                .filter_map(|(k, b)| HashMap::get(self, k).map(|a| a.conj() * b.clone()))
                .sum()
        }
    }
}

impl<K, C> Retain<K, C> for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        // Inherent `HashMap::retain` shadows the trait method; no recursion.
        self.retain(|k, v| keep(k, v));
    }
}

impl<K, C> Multiply for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable + KeyProduct,
    C: ImaginaryUnit,
{
    /// The twisted convolution `(A·B)[k] = Σ_{p·q = k} A[p]·B[q]·i^{β(p,q)}`,
    /// accumulated into `acc` through the hash join — `twistedConv` of
    /// `lean/PPVM/Algebra/Twisted.lean` (monomial case `tmul`; associative by
    /// `tmul_assoc`/`gtmul_assoc`, and `(A·B)[I] = ⟨A, B⟩` by
    /// `twistedConv_apply_id`, tying L4 back to [`Pair::overlap`]).
    ///
    /// Every `(p, q)` pair contributes: the outer product is `O(|A|·|B|)` and is
    /// accumulated **into a distinct `acc`**, never folded back into an operand.
    /// (That is the bilinearity old's `MulAssign<PauliSum>` loses — see
    /// `crate::graded::Multiply` and `ppvm-pauli-sum-2::multiply`.)
    ///
    /// Neither `reduce` nor any truncation runs here: an exact-zero cancellation
    /// stays in `acc`'s support. Canonicalization is the caller's explicit
    /// [`Accumulate::reduce`], per §"`reduce()` is first-class, and runs only at
    /// finalize".
    fn multiply_into(&self, other: &Self, acc: &mut Self) {
        for (p, a) in HashMap::iter(self) {
            for (q, b) in HashMap::iter(other) {
                let (k, phase) = p.key_mul(q);
                let c = phase.apply(&(a.clone() * b.clone()));
                // Warm the fresh product key's structural digest *before* the
                // probe, for the reason `RotateInPlace` does: `key_mul` returns a
                // key with an empty hash cache, and letting the finalize fold fire
                // lazily inside `entry()` puts its mul-chain latency on the
                // bucket-index critical path with nothing to hide it. Semantic
                // no-op (identical digest).
                let _ = k.key_hash();
                acc.entry(k).and_modify(|e| *e += c.clone()).or_insert(c);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::batch::TermSink;

    /// A minimal [`Indexable`] key, so the hash backend is testable here (the
    /// only real one, `PauliWord`, lives downstream). `Hash` is exactly
    /// `write_u64(key_hash())`, as the contract requires.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct Key(u64);

    impl std::hash::Hash for Key {
        fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
            state.write_u64(self.key_hash());
        }
    }

    impl Indexable for Key {
        fn key_hash(&self) -> u64 {
            self.0.wrapping_mul(0x9E37_79B9_7F4A_7C15)
        }
    }

    fn batch(terms: &[(&str, f64)]) -> TermBatch<String, f64> {
        let mut b = TermBatch::with_capacity(terms.len());
        for (k, c) in terms {
            b.push((*k).to_string(), *c);
        }
        b
    }

    #[test]
    fn vec_accumulate_combines_keys() {
        let mut v: Vec<(String, f64)> = Vec::new();
        v.accumulate_batch(&batch(&[("a", 1.0), ("b", 2.0), ("a", 3.0)]));
        assert_eq!(Support::get(&v, &"a".to_string()), Some(4.0));
        assert_eq!(Support::get(&v, &"b".to_string()), Some(2.0));
        assert_eq!(Support::len(&v), 2);
    }

    #[test]
    fn vec_reduce_drops_zero() {
        let mut v: Vec<(String, f64)> = Vec::new();
        v.accumulate_batch(&batch(&[("a", 1.0), ("a", -1.0), ("b", 2.0)]));
        v.reduce();
        assert_eq!(Support::len(&v), 1);
        assert_eq!(Support::get(&v, &"b".to_string()), Some(2.0));
    }

    #[test]
    fn vec_scale_and_overlap() {
        let mut a: Vec<(String, f64)> = Vec::new();
        a.accumulate_batch(&batch(&[("x", 2.0), ("y", 3.0)]));
        let mut b: Vec<(String, f64)> = Vec::new();
        b.accumulate_batch(&batch(&[("x", 5.0), ("z", 7.0)]));
        a.scale(&2.0);
        // overlap = (2*2)*5 = 20; y and z do not match.
        assert_eq!(Pair::overlap(&a, &b), 20.0);
    }

    /// The borrowing scan must see exactly what [`Support::iter`] sees — same
    /// multiset of `(key, coeff)` pairs — on both backends. It is an *override*
    /// of a defaulted method, so nothing but a test keeps the two in step, and a
    /// reader that silently visited a stale or partial view would corrupt every
    /// trace.
    #[test]
    fn for_each_ref_agrees_with_iter_on_both_backends() {
        let terms = batch(&[("a", 1.0), ("b", 2.0), ("c", -3.0), ("a", 0.5)]);

        let mut v: Vec<(String, f64)> = Vec::new();
        v.accumulate_batch(&terms);
        let mut seen: Vec<(String, f64)> = Vec::new();
        v.for_each_ref(|k, c| seen.push((k.clone(), *c)));
        seen.sort_by(|a, b| a.0.cmp(&b.0));
        let mut want: Vec<(String, f64)> = Support::iter(&v).collect();
        want.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(seen, want);
        assert_eq!(seen.len(), 3);

        let mut m: HashMap<Key, f64, IdentityBuildHasher> = HashMap::default();
        for (k, c) in [(1u64, 1.0), (2, 2.0), (3, -3.0)] {
            m.insert(Key(k), c);
        }
        let mut seen: Vec<(Key, f64)> = Vec::new();
        m.for_each_ref(|k, c| seen.push((*k, *c)));
        seen.sort_by_key(|(k, _)| k.0);
        let mut want: Vec<(Key, f64)> = Support::iter(&m).collect();
        want.sort_by_key(|(k, _)| k.0);
        assert_eq!(seen, want);
        assert_eq!(seen.len(), 3);
    }

    #[test]
    fn vec_retain_filters() {
        let mut v: Vec<(String, f64)> = Vec::new();
        v.accumulate_batch(&batch(&[("keep", 2.0), ("drop", 0.5)]));
        Retain::retain(&mut v, |_, c| *c >= 1.0);
        assert_eq!(Support::len(&v), 1);
        assert_eq!(Support::get(&v, &"keep".to_string()), Some(2.0));
    }
}
