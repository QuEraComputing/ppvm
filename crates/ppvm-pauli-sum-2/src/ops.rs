// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The value-level surface of [`Sum`]: equality (exact, support-sensitive,
//! width-sensitive), approximate equality, and the accumulating operators
//! `+= (key, coeff)` / `+= key` / `+= &sum` / `*= scalar`.
//!
//! Ported from `ppvm-pauli-sum/src/sum/{data.rs, approx.rs, ops.rs}`. These are
//! **user-facing contracts**, not conveniences:
//!
//! * `PartialEq` compares `n_sites` **and** the support exactly — zero-coefficient
//!   terms included. `ppvm-pauli-sum/tests/loss.rs::test_reset_channel` asserts
//!   `state == (state2 *= 0.0)`, which only holds because a zeroed term stays in
//!   the support on both sides and equality counts it.
//! * `+= (key, coeff)` **accumulates** onto an existing key (old's
//!   `ACMapAddAssign::add_assign`); it does not replace. `+= key` adds
//!   coefficient one.
//! * `*= s` scales every coefficient, zeros included, and removes nothing.
//!
//! None of these truncate or [`reduce`](Sum::reduce) — the deferred-truncation
//! contract (`ppvm-pauli-sum/tests/truncation_semantics.rs`) covers every
//! insertion path, not just gates.

use std::ops::AddAssign;

use num::One;
use ppvm_pauli_word_2::{PauliStorage, PauliWord};
use ppvm_traits_2::{Accumulate, Indexable, Word};

use crate::policy::Policy;
use crate::store::{AddTerm, InsertTerm};
use crate::sum::Sum;

/// Sealed marker for the types the *bare* `sum += key` form accepts.
///
/// Old's `IntoPauliWord`/`SealedIntoPauliWord` pattern
/// (`ppvm-pauli-sum/src/sum/ops.rs`), and for the same reason: `AddAssign<K>`
/// and `AddAssign<(K, C)>` would otherwise overlap, since coherence cannot rule
/// out `K` being a tuple. Restricting the bare form to a sealed trait that no
/// tuple implements makes the two impls provably disjoint.
mod sealed {
    /// Types accepted by the bare `sum += key` operator. Sealed: no tuple
    /// implements it, which is what keeps the two `AddAssign` impls disjoint.
    pub trait BareKey {}
}

impl<A: PauliStorage, H> sealed::BareKey for PauliWord<A, H> {}
impl sealed::BareKey for &str {}
impl sealed::BareKey for String {}

// --- Single-term accumulation ------------------------------------------------

impl<S, P> Sum<S, P>
where
    S: Accumulate + AddTerm<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Add `coeff` onto the coefficient at `key`, inserting the term if absent.
    ///
    /// Old's `sum += (word, coeff)`: an **accumulate**, never a replace, and a
    /// zero coefficient is inserted rather than dropped. Neither the policy's
    /// truncation nor `reduce` runs.
    ///
    /// The key's width must match the sum's — a `debug_assert!`, as in old
    /// (`debug_assert_eq!(self.n_qubits(), key.n_qubits())`): a mismatch panics
    /// in debug and is unchecked in release.
    #[inline(always)]
    pub fn add_term(&mut self, key: S::Key, coeff: S::Coeff) {
        debug_assert_eq!(
            key.n_sites(),
            self.n_sites(),
            "term key width must match the sum's n_sites"
        );
        self.storage_mut().add_term(key, coeff);
    }
}

impl<K, S, P> AddAssign<(K, S::Coeff)> for Sum<S, P>
where
    K: Into<S::Key>,
    S: Accumulate + AddTerm<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// `sum += (key, coeff)` — see [`Sum::add_term`]. `key` is anything the key
    /// type converts from, so the old crate's string form
    /// `sum += ("IIZI", 1.0)` carries over verbatim.
    #[inline(always)]
    fn add_assign(&mut self, rhs: (K, S::Coeff)) {
        self.add_term(rhs.0.into(), rhs.1);
    }
}

impl<K, S, P> AddAssign<K> for Sum<S, P>
where
    K: Into<S::Key> + sealed::BareKey,
    S: Accumulate + AddTerm<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
    S::Coeff: One,
{
    /// `sum += key` — add the term with coefficient one (old's bare-word `+=`).
    #[inline]
    fn add_assign(&mut self, rhs: K) {
        self.add_term(rhs.into(), S::Coeff::one());
    }
}

// --- Sum + sum: free-module addition ----------------------------------------

impl<S, P> Sum<S, P>
where
    S: Accumulate + AddTerm<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Add another sum onto this one, **pointwise on coefficients**:
    /// `∀ k. self[k] += other[k]`.
    ///
    /// # The one place this diverges from old, and why
    ///
    /// Old ships `impl AddAssign<PauliSum<T>> for PauliSum<T>` (and the
    /// `&PauliSum` variant) as `self.data_mut().extend(rhs)`
    /// (`ppvm-pauli-sum/src/sum/ops.rs:117-140`). `Extend` for `HashMap`/`IndexMap`
    /// **inserts**, i.e. *replaces* the value on a duplicate key — so old's
    /// `A += B` overwrites the coefficient of every shared key with `B`'s instead
    /// of summing: `1.0·ZZ += 2.0·ZZ` yields `2.0·ZZ`, not `3.0·ZZ`. That
    /// contradicts old's own single-term path (`sum += (word, coeff)`, which
    /// routes through the accumulating `ACMapAddAssign::add_assign`), and no test
    /// in the old crate covers sum-plus-sum.
    ///
    /// The Lean oracle adjudicates: free-module addition of two
    /// finitely-supported maps is pointwise coefficient **addition**
    /// (`lean/PPVM/Algebra/GradedMap.lean` — the `C[K]` module laws, with
    /// `accumulateTerms_add` for partitioned batches), so old is wrong and this
    /// accumulates. It is the sibling of the `MulAssign<PauliSum>` correction in
    /// [`crate::multiply`]: a **documented, Lean-backed** divergence, the only
    /// class of behaviour change this port allows.
    ///
    /// Like every other operation here it runs **no** `reduce` and **no**
    /// truncation, so a key that cancels to exactly zero stays in the support.
    ///
    /// The widths must match — a `debug_assert!`, as on every other insertion
    /// path (old's `debug_assert_eq!(self.n_qubits(), key.n_qubits())`).
    pub fn add_sum(&mut self, other: &Self) {
        debug_assert_eq!(
            self.n_sites(),
            other.n_sites(),
            "sum addition requires equal-width operands"
        );
        // Term-by-term through the single-term door rather than a `TermBatch`
        // round-trip: `AddTerm` *is* old's `add_assign`, and staging the rhs into
        // a batch first would allocate two `Vec`s and clone every key twice for a
        // merge that needs neither.
        for (k, c) in other.iter() {
            self.storage_mut().add_term(k, c);
        }
    }
}

/// `a += &b` — pointwise coefficient addition ([`Sum::add_sum`]); **accumulates**
/// on shared keys, where old's `extend`-based operator overwrote (see
/// [`Sum::add_sum`] for the Lean adjudication).
impl<S, P> AddAssign<&Sum<S, P>> for Sum<S, P>
where
    S: Accumulate + AddTerm<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    #[inline]
    fn add_assign(&mut self, rhs: &Sum<S, P>) {
        self.add_sum(rhs);
    }
}

/// `a += b` — the by-value form of [`Sum::add_sum`]. Old has this overload too
/// (`impl AddAssign<PauliSum<T>> for PauliSum<T>`), so a caller porting from old
/// finds the same two spellings.
impl<S, P> AddAssign<Sum<S, P>> for Sum<S, P>
where
    S: Accumulate + AddTerm<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    #[inline]
    fn add_assign(&mut self, rhs: Sum<S, P>) {
        self.add_sum(&rhs);
    }
}

/// `&a + &b` — the fresh sum, i.e. clone-then-[`Sum::add_sum`].
impl<S, P> std::ops::Add<&Sum<S, P>> for &Sum<S, P>
where
    S: Accumulate + AddTerm<S::Key, S::Coeff> + Clone,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    type Output = Sum<S, P>;

    #[inline]
    fn add(self, rhs: &Sum<S, P>) -> Sum<S, P> {
        let mut out = self.clone();
        out.add_sum(rhs);
        out
    }
}

// --- Scalar multiplication ---------------------------------------------------

/// `sum *= s` / `sum * s` for one concrete coefficient ring.
///
/// Concrete rather than generic over `S::Coeff` for a coherence reason: a
/// blanket `MulAssign<S::Coeff>` overlaps the operator *product*
/// `MulAssign<&Sum>` (`crate::multiply`), since coherence cannot rule out
/// `S::Coeff` being `&Sum`. Old has the same shape for the same reason — a macro
/// instantiated per scalar type (`impl_op_mul_assign_coefficient!(f64)`).
///
/// # Why it is exported
///
/// The orphan rule puts the instantiation for a coefficient ring defined
/// *outside* this crate out of this crate's reach and out of that crate's reach
/// alike: `Sum<S, P>` is foreign to `ppvm-sym-2`, and its type parameters precede
/// the local `Term`, so `ppvm-sym-2` cannot write `MulAssign<Term> for Sum<..>`
/// by hand. Exporting the macro — which places the impl in *this* crate's
/// expansion, syntactically inside the downstream crate — is the same escape
/// hatch old used: old `ppvm-sym` took a real (non-dev) dependency on
/// `ppvm-pauli-sum` for the single line `impl_op_mul_assign_coefficient!(Term)`,
/// which is what made `sum *= Term::from(2.0)` compile. `ppvm-sym-2` does the
/// same with this macro, so that spelling survives the port.
///
/// All paths in the expansion are `$crate`- or `::`-absolute, so the macro needs
/// nothing in scope at the call site:
///
/// ```
/// # use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, PauliWord, Sum};
/// # #[derive(Clone, Copy, Debug, PartialEq, Default)]
/// # struct MyScalar(f64);
/// // (in a downstream crate that owns the coefficient ring)
/// // ppvm_pauli_sum_2::impl_scalar_mul!(MyScalar);
/// let mut s: Sum<HashMapStore<PauliWord<[u8; 8]>, f64>, NoPolicy> = Sum::new(2);
/// s += (PauliWord::from("ZZ"), 1.0);
/// s *= 2.0;
/// assert_eq!(s.get(&PauliWord::from("ZZ")), Some(2.0));
/// ```
#[macro_export]
macro_rules! impl_scalar_mul {
    ($ty:ty) => {
        impl<S, P> ::core::ops::MulAssign<$ty> for $crate::Sum<S, P>
        where
            S: $crate::reexport::Accumulate<Coeff = $ty> + $crate::reexport::Scale,
            P: $crate::Policy<S::Key, $ty>,
            S::Key: $crate::reexport::Word + $crate::reexport::Indexable,
        {
            /// Scale every coefficient, zeros included; nothing is removed, so
            /// `sum *= 0.0` keeps the whole key set at `0.0` (old's
            /// `MulAssign<f64>` over `ACMapScale`, which
            /// `ppvm-pauli-sum/tests/loss.rs::test_reset_channel` compares
            /// against).
            #[inline]
            fn mul_assign(&mut self, rhs: $ty) {
                self.scale(&rhs);
            }
        }

        impl<S, P> ::core::ops::Mul<$ty> for $crate::Sum<S, P>
        where
            S: $crate::reexport::Accumulate<Coeff = $ty> + $crate::reexport::Scale,
            P: $crate::Policy<S::Key, $ty>,
            S::Key: $crate::reexport::Word + $crate::reexport::Indexable,
        {
            type Output = Self;

            /// By-value `*=`, i.e. old's clone-then-`*=` without the forced clone.
            #[inline]
            fn mul(mut self, rhs: $ty) -> Self {
                self *= rhs;
                self
            }
        }
    };
}

impl_scalar_mul!(f64);
impl_scalar_mul!(num::Complex<f64>);

// --- Collection compatibility ------------------------------------------------

impl<S, P> IntoIterator for Sum<S, P>
where
    S: Accumulate,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    type Item = (S::Key, S::Coeff);
    type IntoIter = std::vec::IntoIter<Self::Item>;

    /// Consume the sum and yield its terms in the backend's unspecified order.
    ///
    /// The generic storage contract exposes a synthesized owned iterator rather
    /// than its concrete collection's drain type, so the terms are staged in a
    /// `Vec`. Concrete backends remain hidden from callers.
    fn into_iter(self) -> Self::IntoIter {
        self.iter().collect::<Vec<_>>().into_iter()
    }
}

impl<S, P> Extend<(S::Key, S::Coeff)> for Sum<S, P>
where
    S: Accumulate + InsertTerm<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Extend using old's backing-map semantics: a duplicate key is replaced,
    /// not accumulated. Algebraic addition continues to use [`Sum::add_sum`].
    #[inline(always)]
    fn extend<I>(&mut self, iter: I)
    where
        I: IntoIterator<Item = (S::Key, S::Coeff)>,
    {
        for (key, coeff) in iter {
            self.storage_mut().insert_term(key, coeff);
        }
    }
}

// --- Equality ----------------------------------------------------------------

/// Exact equality: same width **and** the same support, key for key and
/// coefficient for coefficient — including zero-coefficient terms.
///
/// Old's `impl PartialEq for PauliSum` (`n_qubits == n_qubits && data() ==
/// data()`). The storage's own `PartialEq` compares only the primary support, so
/// transient `aux`/`scratch` state is invisible here, matching old (which
/// compares `data()` only, never `aux()`). The capacity hint is likewise not part
/// of the value.
impl<S, P> PartialEq for Sum<S, P>
where
    S: Accumulate + PartialEq,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.n_sites() == other.n_sites() && self.storage() == other.storage()
    }
}

/// Approximate equality: equal width, equal support **size**, and a per-key
/// comparison — so a sum carrying an extra zero-coefficient term is *not*
/// approx-equal to one without it. Old's `AbsDiffEq for PauliSum`
/// (`ppvm-pauli-sum/src/sum/approx.rs`).
impl<S, P> approx::AbsDiffEq for Sum<S, P>
where
    S: Accumulate + PartialEq,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
    S::Coeff: approx::AbsDiffEq,
    <S::Coeff as approx::AbsDiffEq>::Epsilon: Copy,
{
    type Epsilon = <S::Coeff as approx::AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        <S::Coeff as approx::AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        if self.n_sites() != other.n_sites() || self.len() != other.len() {
            return false;
        }
        self.iter()
            .all(|(k, v)| other.get(&k).is_some_and(|ov| v.abs_diff_eq(&ov, epsilon)))
    }
}

/// Old's `RelativeEq for PauliSum`; same width/size gates as
/// [`approx::AbsDiffEq`].
impl<S, P> approx::RelativeEq for Sum<S, P>
where
    S: Accumulate + PartialEq,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
    S::Coeff: approx::RelativeEq,
    <S::Coeff as approx::AbsDiffEq>::Epsilon: Copy,
{
    fn default_max_relative() -> Self::Epsilon {
        <S::Coeff as approx::RelativeEq>::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        if self.n_sites() != other.n_sites() || self.len() != other.len() {
            return false;
        }
        self.iter().all(|(k, v)| {
            other
                .get(&k)
                .is_some_and(|ov| v.relative_eq(&ov, epsilon, max_relative))
        })
    }
}
