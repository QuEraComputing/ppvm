// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Symbolic and exact coefficient rings over the `-2` trait tower.
//!
//! `ppvm-sym-2` provides two coefficient types, both implementing
//! `ppvm-traits-2`'s [`Coefficient`](ppvm_traits_2::Coefficient), so
//! `ppvm-pauli-sum-2::Sum` instantiates over them exactly as old `PauliSum`
//! could over `ppvm-sym::Term`:
//!
//! * [`Term`] — a polynomial in sines and cosines of symbolic parameters, used
//!   to propagate Pauli operators through *parametric* circuits (variational
//!   ansätze) without committing to angle values until the very end. Ported
//!   from `ppvm-sym`.
//! * [`GaussianInt`] — the Gaussian integers `ℤ[i]`, an **exact** ring with no
//!   `f64` in its representation. This is the runtime witness for the design
//!   claim that splitting `Halvable`/`Angle`/`ImaginaryUnit`/`Conjugate` off
//!   `Coefficient` keeps exact rings expressible (implementation plan §Phase 5).
//!
//! # Quick example
//!
//! Build the symbolic expression `sin(x0) * cos(x1)`, then evaluate it
//! at a concrete `(x0, x1)`:
//!
//! ```
//! use ppvm_sym_2::Term;
//!
//! let expr = Term::var(0).sin() * Term::var(1).cos();
//!
//! let v = expr.eval(&[0.5, 1.0]).unwrap();
//! let expected = 0.5_f64.sin() * 1.0_f64.cos();
//! assert!((v - expected).abs() < 1e-12);
//! ```
//!
//! # Behaviour parity
//!
//! This is a **behaviour-preserving** port of `ppvm-sym`. The truncation
//! discipline (drop-at-accumulate, parameters carried inline on every `Term` and
//! inherited from the left-hand side only), the four-way [`Inner`]
//! representation, the representational `PartialEq`, the partial `sin`/`cos`, the
//! `eval` error shape and the `Display` rendering are all old's. The individual
//! contracts are documented on the items that implement them; the deliberate
//! divergences are listed below.
//!
//! # Divergences from old (all deliberate, all Lean-adjudicated)
//!
//! 1. **`Term += f64` on a `One` receiver is no longer a silent no-op**
//!    (`ppvm-sym/src/add.rs:24-28` built the promoted `Sum` into a local and
//!    never assigned it). A coefficient ring must be an additive monoid
//!    (`lean/PPVM/Algebra/GradedMap.lean` `accumulate_comm`/`accumulate_assoc`),
//!    so `x + c == x` for `c != 0` is indefensible. See
//!    `impl AddAssign<f64> for Term`.
//! 2. **Monomial multiplication composes the phase** (`ppvm-sym/src/mul.rs:63-73`
//!    dropped it): `i^a · i^b = i^{a+b}`, `lean/PPVM/Algebra/Twisted.lean`
//!    `iPow_add`, and `tmul_assoc` needs it. See `impl MulAssign<Prod> for Prod`.
//! 3. **`Sum::add_term` no longer discards the phase of a phase-only monomial**
//!    (`ppvm-sym/src/term.rs:126-129` short-circuited on `pow() == 0` alone,
//!    which silently unphased a symbolic sum's constant part in `mul_phase`).
//!    `lean/PPVM/Instantiations/Symbolic.lean` `phaseFold_eq_iSym_pow_mul` proves
//!    the key relabelling `mul_phase` performs *is* the ring product `iᵏ · x`
//!    (read-out corollary `evalC_phaseFold`), `phaseFold_const` is the constant
//!    summand's arm, and `phaseFold_drop_const_ne` proves old's spelling computes
//!    a *different function* — so the divergence is forced rather than chosen.
//!    (`twistedConv_add_left`/`twistedConv_add_right` and `iPow_add` are the
//!    supporting Pauli-key analogues.)
//! 4. **A phase-aware [`Term::eval_complex`] is added** alongside the real
//!    [`Term::eval`], which still ignores the phase exactly as old did (so every
//!    real-valued golden master is untouched). Old's `f64`-only `eval` made
//!    `Term`'s complex capability unobservable (`lean/PPVM/Pauli/Matrix.lean`
//!    `star_iU`).
//! 5. **[`Halvable`](ppvm_traits_2::Halvable) is not implemented.** Old's
//!    `Coefficient::half` returned the *constant* `0.5` regardless of `self`,
//!    violating the `-2` law `x.half() + x.half() == x`; it was a workaround for
//!    a caller bug in old `ppvm-pauli-sum/src/sum/proj.rs`
//!    (`let half = v.half(); *v *= half;` computes `v²/2` for `f64`). The
//!    implementation plan's Phase 5 prescribes exactly this: the exact-ring
//!    witness does not implement `Halvable`. Nothing on the live path is
//!    affected — `ppvm-pauli-sum-2` has no `Projection`/`p0`/`p1`.
//! 6. **`Display` tie-break order among monomials sharing `(sin_pow, cos_pow)`
//!    differs from old.** The `Sum` sort key is old's, and is a *non-total*
//!    order; ties fall back to the monomial table's iteration order, which is a
//!    function of the hash values and the bucket occupancy. Both of the demanded
//!    perf features change that: the packed-vector monomial layout changes each
//!    key's digest, and the `mul_term` aux double-buffer changes the table's
//!    capacity history. Everything the ordering is *about* — the monomial set,
//!    the coefficients, the `(sin_pow, cos_pow)` grouping, the `[]` rendering of
//!    an empty `Sum` — is byte-identical to old, verified per key in
//!    `tests/engine_symbolic.rs`, and the `examples/symbolic.rs` snapshot (whose
//!    monomials all have distinct sort keys) matches byte for byte, as do both of
//!    old's `Display` unit-test snapshots. Measured extent: **9.8% of printed
//!    coefficients (397 of 4050) across the real workloads**, pinned by
//!    `tests/sym_diff.rs::display_tie_order_divergence_rate_on_the_real_workloads`
//!    so it can never drift further unnoticed. There is no way to keep both the
//!    aux buffer and old's tie order, and imposing a *total* order would be a
//!    different divergence the baseline also says must be flagged; the divergence
//!    is reported for a human ruling rather than silently absorbed.
//!
//! # Preserved old quirks that look like bugs but are not diverged from
//!
//! * **`max_sin` is only enforced on the map-backed form.** `Sum::add_term` /
//!   `Sum::mul_term` are the only places the bound is consulted, so the
//!   `One × One → One` and `Const × One → One` fast arms let a coefficient that
//!   stays a single monomial for a whole circuit exceed `max_sin` without limit.
//!   Verified against old on the Trotter replay in `tests/engine_symbolic.rs`
//!   (both crates leave exactly one such escapee, with `sin_pow = 7` at
//!   `max_sin = 2`).
//!
//!   This is now an *invariant*, not one observation.
//!   `lean/PPVM/Instantiations/Symbolic.lean` models the four-way `Inner` as an
//!   abstraction function `den` plus the implemented product `mulImpl k`, and
//!   proves `mulImpl_not_wellDefined`: `mulImpl k` does **not** factor through
//!   `den` for finite `k` — two representations of the same polynomial multiply to
//!   different polynomials. `mulImpl_one_one_untruncated` is the positive half
//!   (the fast arm computes the *untruncated* ring product for every `k`) and
//!   `fastArm_escapes_bound` the escapee. Two consequences are load-bearing here:
//!   `mulMono_drop_at_insert_eq_drop_at_end` may **not** be cited as an end-to-end
//!   guarantee that the propagated coefficient equals the truncated ring product
//!   (it covers the map-backed accumulation only, and `set_max_sin` is therefore
//!   not a hard degree bound on the result); and unifying the four `Inner` arms
//!   onto one map-backed representation — the regression the integration
//!   baseline's perf feature 1 warns about — would make the product well-defined
//!   and so change numbers, which is a spec violation rather than a tidy-up.
//!
//! # Deferrals
//!
//! Intentionally **not** ported / not done in this component:
//!
//! * **`Coefficient::magnitude`'s absolute-value law is knowingly violated** for
//!   the symbolic forms (`+∞`), because that is the only value reproducing old's
//!   inert `cutoff`. No absolute value exists on `R[sᵢ, cᵢ]` at all (the natural
//!   `ℓ¹` norm is only sub-multiplicative). The adjudication has landed in
//!   `lean/PPVM/Algebra/Truncation.lean` and settles it *against the law rather
//!   than against parity*: `l1_bound_seminorm` weakens `l1_bound_abv`'s
//!   `AbsoluteValue` to a nonnegative, `0`-vanishing, subadditive,
//!   **sub**-multiplicative seminorm and the truncation bound survives — so an
//!   `ℓ¹` `magnitude` *would* be sound — while `l1_bound_seminorm_needs_zero`
//!   shows the one clause that cannot be dropped is `N 0 = 0`, exactly the one
//!   `+∞` breaks. Switching to `ℓ¹` would start dropping terms old kept, which
//!   the prime directive forbids, so parity wins and the consequence is recorded
//!   explicitly: **`CoefficientThreshold` is inert on symbolic coefficients and
//!   carries no `ℓ¹` error bound there.** Whether that law exemption is accepted
//!   permanently is a human ruling. See
//!   [`Coefficient::magnitude`](ppvm_traits_2::Coefficient::magnitude)'s impl for
//!   `Term`, `tests/sym_lean.rs::symbolic_magnitude_deliberately_violates_the_absolute_value_law`
//!   and `tests/sym_diff.rs::magnitude_reproduces_old_cutoff_exactly`.
//! * **`Inner::Var` is kept**, so all of old's `bare variable is not allowed`
//!   panics are reproduced rather than being made unreachable by construction.
//!   `ppvm-traits-2::Angle<C>` already separates the angle domain and would let
//!   `Var` be a distinct type (turning eight runtime panics into compile errors),
//!   but that is a behaviour change and the prime directive rules it out here.
//!   Recorded for adjudication.
//! * **No lazy `OnceLock<u64>` structural hash.** The design's lazy-hash contract
//!   is stated for *keys* (`Indexable::key_hash`); neither `Term` nor `Prod` is a
//!   `Sum` key, and a cache field would need interior mutability, which the
//!   integration baseline's perf feature 10 rules out for this type. The monomial
//!   table hashes the packed factor vector directly, in one pass, with a
//!   seed-free FxHash-class hasher.
//! * **The `ppvm-pauli-sum-2` dependency is real, not dev-only** — for exactly
//!   the one line old needed it for. The engine is generic over
//!   `C: Coefficient`, so no *algorithm* here needs it, but old `ppvm-sym` took a
//!   real dependency on `ppvm-pauli-sum` solely to instantiate the exported
//!   `impl_op_mul_assign_coefficient!(Term)` macro, which is what made
//!   `sum *= Term::from(2.0)` compile on a symbolic Pauli sum. The orphan rule
//!   forbids writing that impl here by hand (`Sum<S, P>` is foreign and its type
//!   parameters precede the local `Term`), so `src/mul.rs` instantiates
//!   `ppvm_pauli_sum_2::impl_scalar_mul!(Term)` — old's shape exactly.
//! * **Differential tests against old `ppvm-sym`, and the perf benches** (the
//!   `sym.*` integration workloads) live in `ppvm-conformance-2`, which is where
//!   the workspace keeps every old-vs-new comparison; this crate ships the
//!   unit-level contract tests and the engine-instantiation tests only.

mod add;
mod coeff;
mod display;
mod eval;
mod exact;
mod mul;
mod term;

pub use exact::GaussianInt;
pub use term::{Factor, Inner, Prod, Sum, Term};
