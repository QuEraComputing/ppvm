// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The symbolic data structures: the canonical monomial [`Prod`], the
//! monomial-table [`Sum`], the four-way [`Inner`] representation, and the
//! user-facing [`Term`].
//!
//! Ported from `ppvm-sym/src/term.rs`. The **algorithm** is old's, verbatim,
//! down to the drop-at-accumulate truncation discipline in [`Sum::add_term`];
//! only the monomial *layout* changed (see [`Prod`]).

use std::cmp::Ordering;

use fxhash::FxHashMap;

/// One variable's contribution to a monomial: `sin(x_var)^sin · cos(x_var)^cos`.
///
/// Packing both exponents for a variable into a single record is what lets a
/// monomial be a **flat sorted vector** instead of old's two `BTreeMap`s.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Factor {
    /// The variable id.
    pub var: u32,
    /// The exponent of `sin(x_var)` (may be `0` if only `cos` is present).
    pub sin: u32,
    /// The exponent of `cos(x_var)` (may be `0` if only `sin` is present).
    pub cos: u32,
}

/// `<phase> sin^m cos^n` — a single monomial over symbolic variables.
///
/// # Layout
///
/// Old (`ppvm-sym::Prod`) stored `sin: BTreeMap<u32, u32>` and
/// `cos: BTreeMap<u32, u32>`. Here the two maps are fused into one **sorted,
/// deduplicated `Vec<Factor>` in ascending variable order**, which is the layout
/// the integration baseline recommends (perf feature 3): `p1.clone() *
/// p2.clone()` runs once per produced monomial in the `Sum × Sum` loop, and old
/// paid *two* `BTreeMap` deep clones (a node allocation per entry) there; this
/// pays one `Vec` allocation.
///
/// Canonicality — the property the whole design rests on, because it is the only
/// thing that stops the symbolic representation from growing as the raw path
/// count — is preserved: the vector is sorted by `var`, holds at most one entry
/// per variable, and never holds an all-zero entry, so **each monomial has a
/// unique representation** and the derived `Hash`/`Eq` are a valid monomial
/// identity (integration baseline, perf feature 3).
///
/// `sin_pow`/`cos_pow` are cached degree totals maintained *incrementally* by
/// [`Prod::mul_sin`]/[`Prod::mul_cos`]/`MulAssign<Prod>` so the `max_sin`
/// truncation test in [`Sum::add_term`] — the innermost operation of the whole
/// crate — stays `O(1)` (perf feature 2).
///
/// The `phase` byte is `k` in `i^k`, `k ∈ ℤ/4`. It is `ppvm-traits-2`'s
/// [`Phase`](ppvm_traits_2::Phase) exponent (`lean/PPVM/Pauli/Phase.lean`,
/// `phaseExp_eq_ref`); unlike old it is composed on multiplication — see
/// `MulAssign<Prod> for Prod` in `crate::mul`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Prod {
    /// Ascending by `var`, one entry per variable, no all-zero entry.
    pub(crate) factors: Vec<Factor>,
    pub(crate) sin_pow: usize,
    pub(crate) cos_pow: usize,
    /// phase factor mod 4, encoded as:
    /// |  | sign | imag |
    /// |--|------|------|
    /// |+1|    0 |    0 |
    /// |+i|    0 |    1 |
    /// |-1|    1 |    0 |
    /// |-i|    1 |    1 |
    pub(crate) phase: u8,
}

impl Default for Prod {
    fn default() -> Self {
        Self::new()
    }
}

impl Prod {
    /// Construct the empty product (the multiplicative identity, value `1`).
    pub fn new() -> Self {
        Self {
            factors: Vec::new(),
            sin_pow: 0,
            cos_pow: 0,
            phase: 0,
        }
    }

    /// Multiply the phase by `i^phase` (modulo `4`).
    pub fn add_phase(&mut self, phase: u8) {
        self.phase = (self.phase + phase) % 4;
    }

    /// The phase exponent `k` of this monomial's `i^k` factor.
    pub fn phase(&self) -> u8 {
        self.phase
    }

    /// Build the singleton product `sin(x_id)`.
    pub fn sin(id: u32) -> Self {
        Self {
            factors: vec![Factor {
                var: id,
                sin: 1,
                cos: 0,
            }],
            sin_pow: 1,
            cos_pow: 0,
            phase: 0,
        }
    }

    /// Build the singleton product `cos(x_id)`.
    pub fn cos(id: u32) -> Self {
        Self {
            factors: vec![Factor {
                var: id,
                sin: 0,
                cos: 1,
            }],
            sin_pow: 0,
            cos_pow: 1,
            phase: 0,
        }
    }

    /// Total power of all sine and cosine factors.
    #[inline]
    pub fn pow(&self) -> usize {
        self.sin_pow + self.cos_pow
    }

    /// Total power of the sine factors (cached; `O(1)`).
    #[inline]
    pub fn sin_pow(&self) -> usize {
        self.sin_pow
    }

    /// Total power of the cosine factors (cached; `O(1)`).
    #[inline]
    pub fn cos_pow(&self) -> usize {
        self.cos_pow
    }

    /// The packed factors, ascending by variable id.
    #[inline]
    pub fn factors(&self) -> &[Factor] {
        &self.factors
    }

    /// The exponent of `sin(x_var)` in this monomial (`0` if absent).
    pub fn sin_exp(&self, var: u32) -> u32 {
        match self.factors.binary_search_by_key(&var, |f| f.var) {
            Ok(i) => self.factors[i].sin,
            Err(_) => 0,
        }
    }

    /// The exponent of `cos(x_var)` in this monomial (`0` if absent).
    pub fn cos_exp(&self, var: u32) -> u32 {
        match self.factors.binary_search_by_key(&var, |f| f.var) {
            Ok(i) => self.factors[i].cos,
            Err(_) => 0,
        }
    }

    /// Merge `f` into the sorted factor vector (exponents add).
    #[inline]
    pub(crate) fn merge_factor(&mut self, f: Factor) {
        match self.factors.binary_search_by_key(&f.var, |g| g.var) {
            Ok(i) => {
                self.factors[i].sin += f.sin;
                self.factors[i].cos += f.cos;
            }
            Err(i) => self.factors.insert(i, f),
        }
    }

    /// Debug-only canonicality/consistency check: ascending unique variables, no
    /// all-zero entry, and the cached degree totals equal the summed exponents.
    ///
    /// The integration baseline (perf feature 2) calls out that an
    /// incrementally-maintained counter which drifts silently changes truncation,
    /// so every mutating operation asserts this in debug builds.
    #[cfg(debug_assertions)]
    pub(crate) fn debug_check(&self) {
        let mut s = 0usize;
        let mut c = 0usize;
        for (i, f) in self.factors.iter().enumerate() {
            assert!(f.sin > 0 || f.cos > 0, "all-zero factor in Prod");
            if i > 0 {
                assert!(self.factors[i - 1].var < f.var, "Prod factors not sorted");
            }
            s += f.sin as usize;
            c += f.cos as usize;
        }
        assert_eq!(s, self.sin_pow, "sin_pow drifted");
        assert_eq!(c, self.cos_pow, "cos_pow drifted");
        assert!(self.phase < 4, "phase out of ℤ/4");
    }

    #[cfg(not(debug_assertions))]
    #[inline]
    pub(crate) fn debug_check(&self) {}
}

/// A formal sum `c₀ + Σᵢ cᵢ · pᵢ`, where each `pᵢ` is a [`Prod`] and
/// `cᵢ` is an `f64` coefficient.
///
/// # The auxiliary double-buffer
///
/// `aux` is the crate's answer to the integration baseline's perf feature 5 (the
/// `ps2.store.aux` failure mode). Old's `Sum::mul_term` did
/// `let mut old_terms = std::mem::take(&mut self.terms)`, which leaves a
/// **zero-capacity** map behind, so the immediately-following re-insertion
/// reallocates and rehashes the whole table from scratch — once per symbolic
/// multiply, i.e. once per coefficient per gate. Here `terms` and `aux`
/// **ping-pong**: the multiply swaps them, drains the old table into the one that
/// already has capacity, and hands the (now empty, still-allocated) buffer back.
/// Neither allocation is ever freed.
///
/// `aux` is **transient**: it is empty between operations, so [`Clone`] starts it
/// fresh and `PartialEq` (representational, per old) observes only `c0`/`terms`.
/// It is kept *inside* the coefficient so `ppvm-pauli-sum-2`'s engine stays
/// generic over `C: Coefficient` and needs no new bound.
#[derive(Debug, Default)]
pub struct Sum {
    pub(crate) c0: f64,
    pub(crate) terms: FxHashMap<Prod, f64>,
    /// Persistent rebuild buffer for [`Sum::mul_term`]; empty between operations.
    pub(crate) aux: FxHashMap<Prod, f64>,
}

impl Clone for Sum {
    /// Clone the value only. `aux` is workspace, empty between operations, so a
    /// clone starts it fresh rather than copying a buffer.
    fn clone(&self) -> Self {
        Self {
            c0: self.c0,
            terms: self.terms.clone(),
            aux: FxHashMap::default(),
        }
    }
}

/// **Representational** equality, matching old's `#[derive(PartialEq)]`: `c0` by
/// exact `f64` comparison and the monomial table by content. A monomial that
/// cancelled to exactly `0.0` stays in the table and therefore *counts*, exactly
/// as in old (behavioural contract 4/5). `aux` is workspace and is not compared.
impl PartialEq for Sum {
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.terms == other.terms
    }
}

impl Sum {
    /// Construct an empty sum (value `0`).
    pub fn new() -> Self {
        Self {
            c0: 0.0,
            terms: FxHashMap::default(),
            aux: FxHashMap::default(),
        }
    }

    /// The constant part `c₀`.
    #[inline]
    pub fn c0(&self) -> f64 {
        self.c0
    }

    /// The number of monomials in the table (excluding `c₀`).
    #[inline]
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Whether the monomial table is empty (`c₀` is not considered).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    /// Iterate the monomial table in the backend's (deterministic, seed-free)
    /// order.
    pub fn iter(&self) -> impl Iterator<Item = (&Prod, f64)> + '_ {
        self.terms.iter().map(|(p, c)| (p, *c))
    }

    /// The coefficient of `p`, or `0.0` if absent.
    pub fn coeff(&self, p: &Prod) -> f64 {
        self.terms.get(p).copied().unwrap_or(0.0)
    }

    /// Add the constant `c` into the sum's `c₀`, dropping it if
    /// `|c| < min_eps`.
    pub fn add_const(&mut self, c: f64, min_eps: f64) {
        if c.abs() < min_eps {
            return;
        }
        self.c0 += c;
    }

    /// Add `coeff · p` into the sum, subject to the same truncation
    /// constraints used elsewhere (`max` caps the sine power, `min_eps`
    /// drops near-zero coefficients).
    ///
    /// Truncation is **drop-at-accumulate**: a rejected monomial is never
    /// materialized in the table (behavioural contract 2, perf feature 7).
    ///
    /// The two axes are *not* equivalent, and `lean/PPVM/Instantiations/Symbolic.lean`
    /// says why. `max` reads only the monomial, and the sine degree is additive
    /// (`sinDeg_add`), so `{p | sinDeg p > max}` is a monomial **ideal**
    /// (`truncIdeal_mul_right`) and dropping at insert is *exact* — equal to
    /// truncating the finished product (`mulMono_drop_at_insert_eq_drop_at_end`,
    /// via `GradedMap.batchMap_filter_key`). `min_eps` reads the *coefficient*,
    /// and `eps_drop_at_insert_ne_drop_at_end` proves that rule is **not**
    /// interchangeable with a post-pass: two sub-threshold contributions to one
    /// monomial are both dropped here but would survive as their sum. That is
    /// why `min_eps` must stay inside this loop rather than move to a
    /// `Sum::truncate()`.
    ///
    /// **Scope.** That exactness claim is about *this* accumulation, not about a
    /// whole propagation: `mulImpl_not_wellDefined` (same file) proves the shipped
    /// four-way product does not factor through the value a [`Term`] denotes,
    /// because the `One × One` / `Const × One` fast arms never reach this function
    /// (`ppvm-sym-2/src/lib.rs` §"Preserved old quirks"). `set_max_sin` bounds
    /// what the map-backed table retains, not the degree of the result.
    ///
    /// # Divergence from old (`oldSuspectedBugs` #3)
    ///
    /// Old short-circuits on `p.pow() == 0` alone, folding the value into `c0`
    /// and **throwing away `p`'s phase**; combined with
    /// [`Term::mul_phase`](crate::Term::mul_phase)'s `Sum` arm that silently
    /// unphases a symbolic sum's constant part. Here the short-circuit
    /// additionally requires `p.phase() == 0`, so a phase-only monomial
    /// (`pow() == 0`, `phase != 0`) is kept in the table.
    ///
    /// That is exactly `phaseFold_const` in
    /// `lean/PPVM/Instantiations/Symbolic.lean` (the key `(0, 0)` must move to
    /// `(0, k)`), and `phaseFold_eq_iSym_pow_mul` is why: the phase relabelling
    /// `mul_phase` performs *is* the ring product `iᵏ · x`, so leaving one summand
    /// behind is not multiplication by `iᵏ` at all —
    /// `phaseFold_drop_const_ne` proves old's spelling computes a different
    /// function. (`lean/PPVM/Algebra/Twisted.lean` `twistedConv_add_left` /
    /// `twistedConv_add_right` and `iPow_add` are the Pauli-key analogues:
    /// additivity in each argument, and no room to drop an `i^k`.)
    pub fn add_term(&mut self, p: Prod, coeff: f64, max: usize, min_eps: f64) {
        if p.sin_pow() > max || coeff.abs() < min_eps {
            return;
        }

        if p.pow() == 0 && p.phase == 0 {
            self.c0 += coeff;
            return;
        }
        *self.terms.entry(p).or_insert(0.0) += coeff;
    }
}

/// Internal representation of a [`Term`].
///
/// `Sum` holds a full formal sum; `One` is the optimisation for a
/// single weighted product; `Var` is a bare symbolic variable used
/// before it has been wrapped in a `sin` or `cos`; `Const` is a
/// numeric scalar.
///
/// # Why all four forms survive the redesign
///
/// Collapsing these onto one map-backed representation "for cleanliness" is the
/// regression the integration baseline's perf feature 1 names: during
/// propagation the overwhelming majority of coefficients are a *single*
/// monomial — the observable starts as `Const(1.0)`, and every `sin`/`cos` a
/// rotation produces is `One(Prod::sin(u), 1.0)`. `One × One → One` and
/// `Const × anything` never touch a hash map at all, saving one `FxHashMap`
/// allocation per coefficient per gate over a deep circuit.
#[derive(Debug, Clone, PartialEq)]
pub enum Inner {
    /// A general sum of products.
    Sum(Sum),
    /// A single weighted product.
    One(Prod, f64),
    /// A bare symbolic variable (only valid as the argument of `sin`
    /// or `cos`).
    Var(u32),
    /// A numeric constant.
    Const(f64),
}

/// A symbolic polynomial in `sin(x_i)` and `cos(x_i)`.
///
/// `Term` is the public-facing wrapper around its [`Inner`] enum. It
/// also carries two truncation parameters, applied during multiplication
/// and addition:
///
/// * `max_sin` — drop terms whose total sine power exceeds this bound.
/// * `min_eps` — drop terms whose coefficient magnitude falls below
///   this threshold.
///
/// Both are stored **inline on every `Term`** and applied *at accumulation
/// time*, and both travel with the coefficient through `clone()`/`mul_sign()`/
/// `sin_cos()`. That placement is load-bearing: `ppvm-pauli-sum-2`'s `Sum` has
/// nowhere to thread a coefficient-level context through
/// [`Coefficient::mul_sign`](ppvm_traits_2::Coefficient::mul_sign) or
/// `Mul<C> for C`, so moving them to a `Policy` would be both a perf regression
/// and a behaviour change (integration baseline, perf feature 7 / behavioural
/// contract 2).
///
/// They are inherited from the **left-hand side only**; the right-hand side's are
/// silently ignored (behavioural contract 1), which is why seeding
/// `set_max_sin`/`set_min_eps` on the initial observable coefficient propagates
/// through a whole circuit: the engine always writes `v.clone() * sin` and
/// `*v *= cos`, with `v` on the left.
///
/// # Examples
///
/// ```
/// use ppvm_sym_2::Term;
///
/// // sin²(x0) at x0 = π/2 equals 1.
/// let expr = Term::var(0).sin() * Term::var(0).sin();
/// let v = expr.eval(&[std::f64::consts::FRAC_PI_2]).unwrap();
/// assert!((v - 1.0).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Term {
    pub(crate) inner: Inner,
    pub(crate) max_sin: usize, // max sin power
    pub(crate) min_eps: f64,   // min coefficient to keep
}

impl Term {
    /// Set the maximum sine power retained during arithmetic.
    pub fn set_max_sin(&mut self, max: usize) {
        self.max_sin = max;
    }

    /// Set the coefficient cutoff used during arithmetic.
    pub fn set_min_eps(&mut self, eps: f64) {
        self.min_eps = eps;
    }

    /// The maximum sine power retained during arithmetic (`usize::MAX` by
    /// default).
    #[inline]
    pub fn max_sin(&self) -> usize {
        self.max_sin
    }

    /// The coefficient cutoff used during arithmetic (`f64::EPSILON` by default).
    #[inline]
    pub fn min_eps(&self) -> f64 {
        self.min_eps
    }

    /// The internal representation. Exposed read-only for differential tests and
    /// diagnostics; the representation itself is old's and is observable through
    /// `PartialEq` (behavioural contract 5), so it is part of the contract.
    #[inline]
    pub fn inner(&self) -> &Inner {
        &self.inner
    }

    /// Construct a bare symbolic variable.
    ///
    /// A bare variable is **only** valid as the argument of [`Term::sin`] /
    /// [`Term::cos`] (or of [`Angle::sin_cos`](ppvm_traits_2::Angle::sin_cos),
    /// which is what a rotation gate calls); every arithmetic operation on it
    /// panics, exactly as in old (behavioural contract 6).
    pub fn var(u: u32) -> Self {
        Self {
            inner: Inner::Var(u),
            max_sin: usize::MAX,
            min_eps: f64::EPSILON,
        }
    }

    /// Apply `sin(·)` to the term. Only valid on variables and
    /// constants (constant-folding on the latter); panics otherwise.
    pub fn sin(mut self) -> Self {
        match &mut self.inner {
            Inner::Var(u) => {
                self.inner = Inner::One(Prod::sin(*u), 1.0);
            }
            Inner::Const(c) => {
                *c = (*c).sin();
            }
            _ => {
                panic!("only variable or constant can be input of sin");
            }
        }
        self
    }

    /// Apply `cos(·)` to the term. Only valid on variables and
    /// constants (constant-folding on the latter); panics otherwise.
    pub fn cos(mut self) -> Self {
        match &mut self.inner {
            Inner::Var(u) => {
                self.inner = Inner::One(Prod::cos(*u), 1.0);
            }
            Inner::Const(c) => {
                *c = (*c).cos();
            }
            _ => {
                panic!("only variable or constant can be input of cos");
            }
        }
        self
    }

    /// Build a constant term from `c`, with the default truncation parameters
    /// (`max_sin = usize::MAX`, `min_eps = f64::EPSILON`).
    pub fn from_f64(c: f64) -> Self {
        Self {
            inner: Inner::Const(c),
            max_sin: usize::MAX,
            min_eps: f64::EPSILON,
        }
    }

    /// The number of monomials this term denotes: `0` for a constant, `1` for a
    /// single weighted product, and the table size for a general sum.
    ///
    /// A diagnostic for the truncation-sweep workloads (monomials produced vs
    /// retained); not part of old's API.
    pub fn n_monomials(&self) -> usize {
        match &self.inner {
            Inner::Sum(s) => s.terms.len(),
            Inner::One(..) => 1,
            Inner::Var(_) | Inner::Const(_) => 0,
        }
    }

    /// The largest `sin_pow` over this term's monomials (`0` if there are none).
    ///
    /// Used by the truncation tests to assert that a seeded `max_sin` really
    /// propagated through a whole circuit.
    pub fn max_monomial_sin_pow(&self) -> usize {
        match &self.inner {
            Inner::Sum(s) => s.terms.keys().map(|p| p.sin_pow()).max().unwrap_or(0),
            Inner::One(p, _) => p.sin_pow(),
            Inner::Var(_) | Inner::Const(_) => 0,
        }
    }
}

/// Merge two sorted factor vectors (exponents add), preserving canonicality.
pub(crate) fn merge_factors(lhs: &[Factor], rhs: &[Factor]) -> Vec<Factor> {
    let mut out = Vec::with_capacity(lhs.len() + rhs.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < lhs.len() && j < rhs.len() {
        let (a, b) = (lhs[i], rhs[j]);
        match a.var.cmp(&b.var) {
            Ordering::Less => {
                out.push(a);
                i += 1;
            }
            Ordering::Greater => {
                out.push(b);
                j += 1;
            }
            Ordering::Equal => {
                out.push(Factor {
                    var: a.var,
                    sin: a.sin + b.sin,
                    cos: a.cos + b.cos,
                });
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&lhs[i..]);
    out.extend_from_slice(&rhs[j..]);
    out
}
