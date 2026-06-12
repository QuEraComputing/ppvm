// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use itertools::Itertools;

use crate::sum::PauliSum;
use ppvm_runtime::config::Config;
use ppvm_runtime::traits::*;
use std::fmt::{Debug, Display};

impl<T: Config> Debug for PauliSum<T>
where
    T::Coeff: Display,
    T::Map: for<'a> ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        for (k, v) in self.data().iter().sorted_by_key(|(k, _)| k.weight()) {
            if !first {
                write!(f, " + ")?;
            }
            write!(f, "{:.8} * {}", v, k)?;
            first = false;
        }
        Ok(())
    }
}

impl<T: Config> Display for PauliSum<T>
where
    T::Coeff: Display,
    T::Map: for<'a> ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        for (k, v) in self.data().iter().sorted_by_cached_key(|(k, _)| k.weight()) {
            if !first {
                write!(f, " + ")?;
            }
            write!(f, "{:.3} * {}", v, k)?;
            first = false;
        }
        Ok(())
    }
}
