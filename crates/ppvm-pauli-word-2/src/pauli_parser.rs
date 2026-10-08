// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Parsing and text formatting for packed Pauli words.

use std::fmt;
use std::hash::BuildHasher;

use bitvec::array::BitArray;
use ppvm_traits_2::{Pauli, Word};

use crate::{HashFinalize, PauliStorage, PauliWord};

impl<A: PauliStorage, H> fmt::Debug for PauliWord<A, H> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PauliWord")
            .field("nqubits", &self.nqubits)
            .field("word", &self.to_string())
            .finish()
    }
}

impl<A: PauliStorage, H> fmt::Display for PauliWord<A, H> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.nqubits {
            let c = match self.get(i) {
                Pauli::I => 'I',
                Pauli::X => 'X',
                Pauli::Y => 'Y',
                Pauli::Z => 'Z',
            };
            f.write_str(c.encode_utf8(&mut [0u8; 1]))?;
        }
        Ok(())
    }
}

impl<A, H> From<&str> for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    /// Parse a Pauli string of `I`/`X`/`Y`/`Z` symbols (underscores are ignored
    /// separators), mirroring `ppvm-pauli-word`'s `From<&str>`. Panics on any
    /// other character or if the width exceeds the backing storage.
    fn from(value: &str) -> Self {
        let mut xbits = BitArray::<A>::ZERO;
        let mut zbits = BitArray::<A>::ZERO;
        let mut i = 0usize;
        for ch in value.chars() {
            match ch {
                'I' => {}
                'X' => xbits.set(i, true),
                'Z' => zbits.set(i, true),
                'Y' => {
                    xbits.set(i, true);
                    zbits.set(i, true);
                }
                '_' => continue,
                other => panic!("invalid Pauli character: {other}"),
            }
            i += 1;
        }
        Self::from_planes(xbits, zbits, i)
    }
}

impl<A, H> From<String> for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn from(value: String) -> Self {
        Self::from(value.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_get_roundtrip() {
        let w: PauliWord = "XZYI".into();
        assert_eq!(w.n_sites(), 4);
        assert_eq!(w.get(0), Pauli::X);
        assert_eq!(w.get(1), Pauli::Z);
        assert_eq!(w.get(2), Pauli::Y);
        assert_eq!(w.get(3), Pauli::I);
        assert_eq!(w.to_string(), "XZYI");
    }

    #[test]
    fn underscore_is_ignored() {
        let a: PauliWord = "X_Y_Z".into();
        let b: PauliWord = "XYZ".into();
        assert_eq!(a, b);
    }
}
