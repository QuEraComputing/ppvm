// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::contains::Contains;
use super::data::PauliPattern;

use crate::loss::LossyPauliWord;
use crate::traits::{PauliStorage, PauliWordTrait, Trace};
use crate::word::PauliWord;
use std::hash::BuildHasher;

impl<'a, A, H> Trace<'a, PauliPattern> for PauliWord<A, H>
where
    A: PauliStorage + 'a,
    H: Default + BuildHasher + Clone + 'a,
{
    type Output = bool;
    fn trace(&'a self, value: &'a PauliPattern) -> Self::Output {
        value.contains(&self)
    }
}

impl<'a, A, H> Trace<'a, PauliPattern> for LossyPauliWord<A, H>
where
    A: PauliStorage + 'a,
    H: Default + BuildHasher + Clone + 'a,
{
    type Output = bool;
    fn trace(&'a self, value: &'a PauliPattern) -> Self::Output {
        value.contains(&self)
    }
}

impl<'a, W> Trace<'a, W> for PauliPattern
where
    W: PauliWordTrait + 'a,
{
    type Output = bool;
    fn trace(&'a self, value: &'a W) -> Self::Output {
        self.contains(value)
    }
}
