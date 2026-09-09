// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

// Loss state associated with word sites
pub trait LossState {
    // Whether site i is lost. LossyPauliWord overrides this
    fn is_lost(&self, _i: usize) -> bool {
        false
    }

    // Number of lost sites
    fn loss_weight(&self) -> usize {
        0
    }

    // Mark index `i` lost
    fn set_lost(&mut self, i: usize);

    // Clear the loss flag at index `i`, returning the site to identity.
    fn clear_lost(&mut self, i: usize);

    /// A copy of this word with the loss flag at `i` cleared.
    fn loss_cleared(&self, i: usize) -> Self
    where
        Self: Sized + Clone,
    {
        let mut out = self.clone();
        out.clear_lost(i);
        out
    }
}
