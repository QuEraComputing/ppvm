// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Broadcasting target sets for gate methods, mirroring stim's `*targets`.

/// A set of qubit indices a gate is applied to. `usize` yields a single
/// target; slices, arrays, `Vec`, and ranges yield each element. Two-qubit
/// gates consume [`Targets::pairs`] (consecutive pairs; odd counts panic).
pub trait Targets {
    /// Iterate the individual target indices in order.
    fn each(self) -> impl Iterator<Item = usize>;

    /// Iterate consecutive `(a, b)` pairs. Panics on an odd count.
    fn pairs(self) -> impl Iterator<Item = (usize, usize)>
    where
        Self: Sized,
    {
        let v: Vec<usize> = self.each().collect();
        assert!(
            v.len() % 2 == 0,
            "two-qubit gate requires an even number of targets, got {}",
            v.len()
        );
        (0..v.len() / 2).map(move |i| (v[2 * i], v[2 * i + 1]))
    }
}

impl Targets for usize {
    fn each(self) -> impl Iterator<Item = usize> {
        std::iter::once(self)
    }
}

impl Targets for &[usize] {
    fn each(self) -> impl Iterator<Item = usize> {
        self.to_vec().into_iter()
    }
}

impl<const N: usize> Targets for [usize; N] {
    fn each(self) -> impl Iterator<Item = usize> {
        self.into_iter()
    }
}

impl Targets for Vec<usize> {
    fn each(self) -> impl Iterator<Item = usize> {
        self.into_iter()
    }
}

impl Targets for std::ops::Range<usize> {
    fn each(self) -> impl Iterator<Item = usize> {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::Targets;

    fn collect(t: impl Targets) -> Vec<usize> {
        t.each().collect()
    }

    #[test]
    fn usize_yields_one() {
        assert_eq!(collect(3usize), vec![3]);
    }

    #[test]
    fn array_and_slice_yield_all() {
        assert_eq!(collect([0usize, 1, 2]), vec![0, 1, 2]);
        let v = vec![5usize, 6];
        assert_eq!(collect(v.as_slice()), vec![5, 6]);
    }

    #[test]
    fn pairs_groups_consecutive() {
        let p: Vec<(usize, usize)> = [0usize, 1, 2, 3].pairs().collect();
        assert_eq!(p, vec![(0, 1), (2, 3)]);
    }

    #[test]
    #[should_panic(expected = "even number of targets")]
    fn pairs_rejects_odd() {
        let _: Vec<_> = [0usize, 1, 2].pairs().collect::<Vec<_>>();
    }
}
