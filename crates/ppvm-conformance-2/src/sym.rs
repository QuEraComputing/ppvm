// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Phase-5 symbolic-coefficient differential harness: matched OLD (`ppvm-sym`
//! `Term` inside `ppvm-pauli-sum::PauliSum`) and NEW (`ppvm-sym-2` `Term` inside
//! `ppvm-pauli-sum-2::Sum`) engines driven through **one** circuit description,
//! so a differential test and a same-build benchmark can never accidentally run
//! two different gate sequences.
//!
//! # Matched configuration
//!
//! Correctness is storage-independent, but a perf ratio is not, so both sides
//! are pinned to the *same* algebraic configuration:
//!
//! | | old | new |
//! |:--|:--|:--|
//! | key storage | `[u8; 8]` (`config::fxhash::Byte<8, Term>`) | `[u8; 8]` (`PauliWord<[u8; 8]>`) |
//! | coefficient | `ppvm_sym::Term` | `ppvm_sym_2::Term` |
//! | monomial hasher | `fxhash` (seed-free) | `fxhash` (seed-free) |
//! | truncation policy | `NoStrategy` | `NoPolicy` |
//!
//! Truncation on this path is **not** a sum-level policy at all: it is intrinsic
//! to the coefficient (`max_sin`/`min_eps` seeded on the initial observable and
//! inherited left-to-right through every gate), which is exactly the behavioural
//! contract the differential suite pins.
//!
//! # Why the comparison goes through `Display` + `eval`
//!
//! Old `ppvm-sym` exposes `Prod`/`Sum` as types but keeps their fields
//! `pub(crate)`, so a differential test outside the crate **cannot** walk an old
//! `Term`'s monomial table directly. The two observables it does expose are
//! `Display` and `eval`, and both are part of the user-facing contract, so the
//! harness compares:
//!
//! * the **monomial set** and the **representation form** (`[…]` map-backed sum
//!   vs bare `c * p` single monomial vs bare scalar), parsed out of `Display` by
//!   [`TermView`] — exact, no tolerance; and
//! * `eval` at a fixed seeded angle vector *and* at many randomized ones —
//!   within `1e-12`.
//!
//! Printed coefficients carry `{:.3}` precision, so [`TermView`] compares those
//! exactly at 3 d.p. and leans on the randomized `eval` sweep (32+ independent
//! angle vectors against a handful of monomials) for the remaining digits: two
//! different coefficient vectors over a fixed monomial basis cannot agree to
//! `1e-12` at 32 random points.

use std::collections::BTreeMap;

use rand::RngExt;
use rand::rngs::StdRng;

// --- OLD engine ------------------------------------------------------------
use ppvm_pauli_sum::config::fxhash::Byte as OldByte;
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_pauli_word::pattern::PauliPattern as OldPauliPattern;
use ppvm_sym::Term as OldTermTy;
use ppvm_traits::traits::{
    Clifford as OldClifford, RotationOne as OldRotationOne, Trace as OldTrace,
};

// --- NEW engine ------------------------------------------------------------
use ppvm_pauli_sum_2::{
    HashMapStore, NoPolicy, PauliPattern as NewPauliPattern, PauliWord as NewPauliWord, Sum,
};
use ppvm_sym_2::Term as NewTermTy;
use ppvm_traits_2::{Clifford as NewClifford, RotationOne as NewRotationOne, Trace as NewTrace};

/// The OLD symbolic coefficient.
pub type OldTerm = OldTermTy;
/// The NEW symbolic coefficient.
pub type NewTerm = NewTermTy;

/// OLD symbolic sum: `[u8; 8]` storage, `FxHash`, no sum-level strategy.
pub type OldSymSum = OldPauliSum<OldByte<8, OldTermTy>>;

/// NEW key, storage-matched to [`OldSymSum`].
pub type NewSymKey = NewPauliWord<[u8; 8]>;
/// NEW symbolic sum, configuration-matched to [`OldSymSum`].
pub type NewSymSum = Sum<HashMapStore<NewSymKey, NewTermTy>, NoPolicy>;

/// Build an empty OLD symbolic sum on `n` qubits.
pub fn new_old_sum(n: usize) -> OldSymSum {
    OldPauliSum::builder().n_qubits(n).build()
}

/// Build an empty NEW symbolic sum on `n` qubits.
pub fn new_new_sum(n: usize) -> NewSymSum {
    NewSymSum::new(n)
}

/// The seeded initial coefficient: `1.0` carrying the caller's `max_sin` /
/// `min_eps`. This is the ONLY brake on symbolic growth (there is no sum-level
/// policy on this path), and both crates inherit it left-to-right through every
/// gate.
pub fn old_seed_coeff(max_sin: usize, min_eps: f64) -> OldTerm {
    let mut c = OldTerm::from(1.0);
    c.set_max_sin(max_sin);
    c.set_min_eps(min_eps);
    c
}

/// The NEW twin of [`old_seed_coeff`].
pub fn new_seed_coeff(max_sin: usize, min_eps: f64) -> NewTerm {
    let mut c = NewTerm::from(1.0);
    c.set_max_sin(max_sin);
    c.set_min_eps(min_eps);
    c
}

// ===========================================================================
// A parsed view of a `Term`, recovered from `Display`.
// ===========================================================================

/// Which of the four `Inner` forms a `Term` rendered as.
///
/// The form is user-visible (`[…]` for the map-backed sum, `c * p` for a single
/// weighted monomial, a bare number for a constant, `%u` for a bare variable)
/// and is part of the behavioural contract — notably contract 9, where the
/// `mul_term` zero shortcut must leave an **empty `Sum`** printing `[]`, not a
/// `Const(0.0)` printing `0`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Form {
    /// `Inner::Sum` — bracketed.
    MapBacked,
    /// `Inner::One` — `c * p`, unbracketed.
    SingleMonomial,
    /// `Inner::Const` — a bare number.
    Scalar,
    /// `Inner::Var` — `%u`.
    Variable,
}

/// A `Term`'s observable structure, parsed from its `Display` rendering.
///
/// This is the *only* structural view available for the old crate (its `Sum`
/// and `Prod` fields are `pub(crate)`), so both sides are compared through it —
/// which has the side benefit that the comparison is over the user-facing text,
/// i.e. exactly what a user could notice.
#[derive(Debug, Clone, PartialEq)]
pub struct TermView {
    /// The rendered representation form.
    pub form: Form,
    /// The constant part, as printed (`0.0` when not printed at all).
    pub c0: String,
    /// `monomial string -> printed coefficient`, sorted, so the (non-total)
    /// `Display` tie order does not enter the comparison.
    pub monomials: BTreeMap<String, String>,
}

impl TermView {
    /// Parse a `Display` rendering.
    ///
    /// Handles old's unbalanced `"[3.000 "` rendering of a `Sum` with a non-zero
    /// `c0` and an empty table (it returns before writing the closing bracket) —
    /// reproduced by the new crate, and pinned by a parity test.
    pub fn parse(s: &str) -> Self {
        if let Some(rest) = s.strip_prefix('%') {
            return Self {
                form: Form::Variable,
                c0: rest.to_string(),
                monomials: BTreeMap::new(),
            };
        }
        if let Some(body) = s.strip_prefix('[') {
            let body = body.strip_suffix(']').unwrap_or(body).trim();
            let mut view = Self {
                form: Form::MapBacked,
                c0: "0".to_string(),
                monomials: BTreeMap::new(),
            };
            if body.is_empty() {
                return view;
            }
            for chunk in body.split(" + ") {
                match chunk.split_once(" * ") {
                    Some((c, p)) => {
                        view.monomials.insert(p.trim().to_string(), c.trim().into());
                    }
                    None => view.c0 = chunk.trim().to_string(),
                }
            }
            return view;
        }
        match s.split_once(" * ") {
            Some((c, p)) => {
                let mut monomials = BTreeMap::new();
                monomials.insert(p.trim().to_string(), c.trim().to_string());
                Self {
                    form: Form::SingleMonomial,
                    c0: "0".to_string(),
                    monomials,
                }
            }
            None => Self {
                form: Form::Scalar,
                c0: s.trim().to_string(),
                monomials: BTreeMap::new(),
            },
        }
    }

    /// The number of monomials (the map table size, or 1 for the single-monomial
    /// form).
    pub fn n_monomials(&self) -> usize {
        self.monomials.len()
    }

    /// The largest total sine power over this view's monomials.
    pub fn max_sin_pow(&self) -> usize {
        self.monomials
            .keys()
            .map(|p| sin_pow_of(p))
            .max()
            .unwrap_or(0)
    }
}

/// Total sine power of a rendered monomial like `"sin^3(%1) cos^1(%0)"`.
pub fn sin_pow_of(prod: &str) -> usize {
    prod.split_whitespace()
        .filter_map(|f| f.strip_prefix("sin^"))
        .filter_map(|f| f.split_once('(').map(|(m, _)| m))
        .filter_map(|m| m.parse::<usize>().ok())
        .sum()
}

/// The OLD term's view.
pub fn old_view(t: &OldTerm) -> TermView {
    TermView::parse(&t.to_string())
}

/// The NEW term's view.
pub fn new_view(t: &NewTerm) -> TermView {
    TermView::parse(&t.to_string())
}

// ===========================================================================
// Support extraction.
// ===========================================================================

/// The OLD sum's support as a sorted `(pauli_string, coefficient)` vector.
pub fn old_sym_support(s: &OldSymSum) -> Vec<(String, OldTerm)> {
    let mut v: Vec<(String, OldTerm)> = s
        .data()
        .iter()
        .map(|(k, c)| (k.to_string(), c.clone()))
        .collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

/// The NEW sum's support as a sorted `(pauli_string, coefficient)` vector.
pub fn new_sym_support(s: &NewSymSum) -> Vec<(String, NewTerm)> {
    let mut v: Vec<(String, NewTerm)> = s.iter().map(|(k, c)| (k.to_string(), c)).collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

// ===========================================================================
// Workload 1 — `sym.trace.parametric` (a verbatim port of `examples/symbolic.rs`).
// ===========================================================================

/// `examples/symbolic.rs` on the OLD engine, returning the trace `Term`.
pub fn parametric_trace_old() -> OldTerm {
    let pat: OldPauliPattern = "Z?*".into();
    let mut sum = new_old_sum(2);
    sum += ("ZZ", OldTerm::from(1.0));

    OldRotationOne::rz(&mut sum, 0, OldTerm::var(0));
    OldRotationOne::ry(&mut sum, 0, OldTerm::var(1));
    OldRotationOne::rz(&mut sum, 0, OldTerm::var(0));

    OldRotationOne::rz(&mut sum, 1, OldTerm::var(0));
    OldRotationOne::ry(&mut sum, 1, OldTerm::var(1));
    OldRotationOne::rz(&mut sum, 1, OldTerm::var(0));

    OldClifford::cnot(&mut sum, 0, 1);

    OldRotationOne::rx(&mut sum, 0, OldTerm::var(1));
    OldRotationOne::rx(&mut sum, 1, OldTerm::var(1));

    OldTrace::trace(&sum, &pat)
}

/// The same circuit on the NEW engine.
pub fn parametric_trace_new() -> NewTerm {
    let mut sum = new_new_sum(2);
    sum += (NewSymKey::from("ZZ"), NewTerm::from(1.0));

    NewRotationOne::rz(&mut sum, 0, NewTerm::var(0));
    NewRotationOne::ry(&mut sum, 0, NewTerm::var(1));
    NewRotationOne::rz(&mut sum, 0, NewTerm::var(0));

    NewRotationOne::rz(&mut sum, 1, NewTerm::var(0));
    NewRotationOne::ry(&mut sum, 1, NewTerm::var(1));
    NewRotationOne::rz(&mut sum, 1, NewTerm::var(0));

    NewClifford::cnot(&mut sum, 0, 1);

    NewRotationOne::rx(&mut sum, 0, NewTerm::var(1));
    NewRotationOne::rx(&mut sum, 1, NewTerm::var(1));

    NewTrace::trace(&sum, &NewPauliPattern::zero_state())
}

// ===========================================================================
// Workload 2/3 — `sym.tfim.trotter` and `sym.truncation.sweep`.
// ===========================================================================

/// A deep symbolic TFIM-Trotter circuit: `L` layers, each `rx(i, var(2l))` on
/// every site followed by `rzz(i, i+1, var(2l+1))` along the chain.
///
/// A **fresh** symbolic variable per layer per gate family, so the monomial
/// space grows with depth instead of collapsing onto one variable — which is
/// what makes `max_sin` load-bearing.
#[derive(Debug, Clone, Copy)]
pub struct TrotterSpec {
    /// Number of qubits.
    pub n: usize,
    /// Number of Trotter layers.
    pub layers: u32,
    /// Seeded `max_sin` on the initial coefficient.
    pub max_sin: usize,
    /// Seeded `min_eps` on the initial coefficient.
    pub min_eps: f64,
    /// The initial observable, e.g. `"ZIIIII"` or `"ZZIIII"`.
    pub observable: &'static str,
}

impl TrotterSpec {
    /// The headline workload: `n = 6`, `L = 6`, `Z₀`.
    pub fn headline(max_sin: usize) -> Self {
        Self {
            n: 6,
            layers: 6,
            max_sin,
            min_eps: 1e-12,
            observable: "ZIIIII",
        }
    }

    /// The number of distinct symbolic variables the circuit uses.
    pub fn n_vars(&self) -> usize {
        2 * self.layers as usize
    }
}

/// Propagate [`TrotterSpec`] on the OLD engine.
///
/// `rzz(a, b, θ)` is decomposed as `cnot(a, b); rz(b, θ); cnot(a, b)` on BOTH
/// sides (rather than calling old's built-in `rzz`) so the two benchmarked and
/// diffed circuits are literally the same gate sequence.
pub fn trotter_old(spec: &TrotterSpec) -> OldSymSum {
    let mut sum = new_old_sum(spec.n);
    sum += (spec.observable, old_seed_coeff(spec.max_sin, spec.min_eps));
    for l in 0..spec.layers {
        for q in 0..spec.n {
            OldRotationOne::rx(&mut sum, q, OldTerm::var(2 * l));
        }
        for q in 0..spec.n - 1 {
            OldClifford::cnot(&mut sum, q, q + 1);
            OldRotationOne::rz(&mut sum, q + 1, OldTerm::var(2 * l + 1));
            OldClifford::cnot(&mut sum, q, q + 1);
        }
    }
    sum
}

/// Propagate [`TrotterSpec`] on the NEW engine.
pub fn trotter_new(spec: &TrotterSpec) -> NewSymSum {
    let mut sum = new_new_sum(spec.n);
    sum += (
        NewSymKey::from(spec.observable),
        new_seed_coeff(spec.max_sin, spec.min_eps),
    );
    for l in 0..spec.layers {
        for q in 0..spec.n {
            NewRotationOne::rx(&mut sum, q, NewTerm::var(2 * l));
        }
        for q in 0..spec.n - 1 {
            NewClifford::cnot(&mut sum, q, q + 1);
            NewRotationOne::rz(&mut sum, q + 1, NewTerm::var(2 * l + 1));
            NewClifford::cnot(&mut sum, q, q + 1);
        }
    }
    sum
}

// ===========================================================================
// Workload 6 — `sym.random.circuit` (deep heterogeneous replay).
// ===========================================================================

/// A replayable gate with a **symbolic** angle: the `GateOp` families the
/// integration baseline asks for (`h`, `s`, `cnot`, `cz`, `rx/ry/rz(var)`).
///
/// `ppvm-conformance-2`'s existing [`GateOp`](crate::GateOp) carries `f64`
/// angles, which cannot address a symbolic variable and has no `cz`/`ry`; this
/// is its symbolic twin, generated from the same [`seeded_rng`](crate::seeded_rng).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymGate {
    /// Hadamard.
    H(usize),
    /// Phase gate `S`.
    S(usize),
    /// CNOT `(control, target)`.
    Cnot(usize, usize),
    /// CZ `(q0, q1)`.
    Cz(usize, usize),
    /// `rx(qubit, %var)`.
    Rx(usize, u32),
    /// `ry(qubit, %var)`.
    Ry(usize, u32),
    /// `rz(qubit, %var)`.
    Rz(usize, u32),
}

/// A seeded random symbolic circuit, `len` gates on `n_qubits`, with rotation
/// angles cycling over `n_vars` symbolic variables.
pub fn random_sym_circuit(
    rng: &mut StdRng,
    n_qubits: usize,
    len: usize,
    n_vars: u32,
) -> Vec<SymGate> {
    assert!(n_qubits >= 2, "the two-qubit gate families need n >= 2");
    (0..len)
        .map(|i| {
            let q = rng.random_range(0..n_qubits);
            let mut other = rng.random_range(0..n_qubits);
            while other == q {
                other = rng.random_range(0..n_qubits);
            }
            let var = (i as u32) % n_vars;
            match rng.random_range(0..7usize) {
                0 => SymGate::H(q),
                1 => SymGate::S(q),
                2 => SymGate::Cnot(q, other),
                3 => SymGate::Cz(q, other),
                4 => SymGate::Rx(q, var),
                5 => SymGate::Ry(q, var),
                _ => SymGate::Rz(q, var),
            }
        })
        .collect()
}

/// Replay a symbolic circuit on the OLD engine.
pub fn replay_old(sum: &mut OldSymSum, circuit: &[SymGate]) {
    for g in circuit {
        match *g {
            SymGate::H(q) => OldClifford::h(sum, q),
            SymGate::S(q) => OldClifford::s(sum, q),
            SymGate::Cnot(a, b) => OldClifford::cnot(sum, a, b),
            SymGate::Cz(a, b) => OldClifford::cz(sum, a, b),
            SymGate::Rx(q, v) => OldRotationOne::rx(sum, q, OldTerm::var(v)),
            SymGate::Ry(q, v) => OldRotationOne::ry(sum, q, OldTerm::var(v)),
            SymGate::Rz(q, v) => OldRotationOne::rz(sum, q, OldTerm::var(v)),
        }
    }
}

/// Replay the same symbolic circuit on the NEW engine.
pub fn replay_new(sum: &mut NewSymSum, circuit: &[SymGate]) {
    for g in circuit {
        match *g {
            SymGate::H(q) => NewClifford::h(sum, q),
            SymGate::S(q) => NewClifford::s(sum, q),
            SymGate::Cnot(a, b) => NewClifford::cnot(sum, a, b),
            SymGate::Cz(a, b) => NewClifford::cz(sum, a, b),
            SymGate::Rx(q, v) => NewRotationOne::rx(sum, q, NewTerm::var(v)),
            SymGate::Ry(q, v) => NewRotationOne::ry(sum, q, NewTerm::var(v)),
            SymGate::Rz(q, v) => NewRotationOne::rz(sum, q, NewTerm::var(v)),
        }
    }
}

// ===========================================================================
// Angle grids.
// ===========================================================================

/// A deterministic angle vector, used as the "fixed seeded" evaluation point.
pub fn fixed_angles(n_vars: usize) -> Vec<f64> {
    (0..n_vars).map(|i| 0.3 + 0.17 * i as f64).collect()
}

/// `count` seeded random angle vectors of `n_vars` entries each, in `[-π, π)`.
pub fn angle_grid(rng: &mut StdRng, n_vars: usize, count: usize) -> Vec<Vec<f64>> {
    (0..count)
        .map(|_| {
            (0..n_vars)
                .map(|_| rng.random_range(-std::f64::consts::PI..std::f64::consts::PI))
                .collect()
        })
        .collect()
}
