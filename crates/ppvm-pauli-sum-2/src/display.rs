// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`Display`] and [`Debug`] for [`Sum`] — the **mathematical** rendering, ported
//! from `ppvm-pauli-sum/src/sum/display.rs`.
//!
//! Behaviour parity (contract 13). Old renders a sum as its terms sorted by
//! **ascending Pauli weight**, joined with `" + "`, each written as
//! `"{coeff} * {word}"` — with `{:.3}` for [`Display`] and `{:.8}` for [`Debug`]
//! — and pins the result with `insta` snapshots. An empty sum renders as the
//! empty string. A `Debug` that printed the struct's fields instead would hand
//! every caller (a snapshot test, the CLI/TUI at cutover, a panic message from
//! `assert_eq!`) a silently different string, so this is a user-facing contract,
//! not a convenience.
//!
//! # The one implementation difference from old
//!
//! Old sorts through `itertools`' `sorted_by_key` / `sorted_by_cached_key`; this
//! collects into a `Vec` and uses the standard library's **stable** `sort_by_key`,
//! which is the same algorithm class and the same stability guarantee, without the
//! dependency. Neither is deterministic *within* one weight class: the input order
//! is the hash map's, which is unspecified on both sides — old's own snapshots are
//! taken on an `IndexMap`-backed (insertion-ordered) config for exactly that
//! reason.

use std::fmt::{self, Debug, Display};

use ppvm_traits_2::{Accumulate, Indexable, Word};

use crate::policy::Policy;
use crate::sum::Sum;

/// Write the terms of `sum` sorted by ascending [`Word::weight`], joined by
/// `" + "`, as `"{coeff:.precision$} * {key}"`.
///
/// Shared by [`Display`] and [`Debug`], which differ from each other only in the
/// coefficient precision (old: `{:.3}` vs `{:.8}`).
fn render<S, P>(sum: &Sum<S, P>, f: &mut fmt::Formatter<'_>, precision: usize) -> fmt::Result
where
    S: Accumulate,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable + Display,
    S::Coeff: Display,
{
    let mut terms: Vec<(S::Key, S::Coeff)> = sum.iter().collect();
    // Stable, so terms of equal weight keep the backend's iteration order —
    // old's `sorted_by_key` is likewise stable.
    terms.sort_by_key(|(k, _)| k.weight());

    let mut first = true;
    for (k, v) in &terms {
        if !first {
            write!(f, " + ")?;
        }
        write!(f, "{:.*} * {}", precision, v, k)?;
        first = false;
    }
    Ok(())
}

/// Old's `Display for PauliSum`: weight-sorted terms joined by `" + "`, each
/// `"{:.3} * {word}"`.
impl<S, P> Display for Sum<S, P>
where
    S: Accumulate,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable + Display,
    S::Coeff: Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        render(self, f, 3)
    }
}

/// Old's `Debug for PauliSum`: the same rendering as [`Display`] at `{:.8}`.
///
/// Deliberately **not** a `derive`/struct dump: `assert_eq!` failure output, the
/// TUI and old's `insta` snapshots all read this string, and old prints the
/// mathematics.
impl<S, P> Debug for Sum<S, P>
where
    S: Accumulate,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable + Display,
    S::Coeff: Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        render(self, f, 8)
    }
}
