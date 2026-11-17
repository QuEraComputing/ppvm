// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Terms are sorted by weight; ties retain storage iteration order.
//! Display uses three coefficient decimals, Debug uses eight; empty sums print nothing.

use std::fmt::{self, Debug, Display};

use ppvm_traits_2::{Support, Word};

use crate::Sum;

fn format_terms<S: Support, P>(
    sum: &Sum<S, P>,
    f: &mut fmt::Formatter<'_>,
    precision: usize,
) -> fmt::Result
where
    S::Key: Word + Display,
    S::Coeff: Display,
{
    let mut terms: Vec<_> = sum.iter().collect();
    terms.sort_by_key(|(key, _)| key.weight());
    for (index, (key, coeff)) in terms.iter().enumerate() {
        if index != 0 {
            write!(f, " + ")?;
        }
        write!(f, "{coeff:.precision$} * {key}")?;
    }
    Ok(())
}

impl<S: Support, P> Display for Sum<S, P>
where
    S::Key: Word + Display,
    S::Coeff: Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        format_terms(self, f, 3)
    }
}

impl<S: Support, P> Debug for Sum<S, P>
where
    S::Key: Word + Display,
    S::Coeff: Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        format_terms(self, f, 8)
    }
}
