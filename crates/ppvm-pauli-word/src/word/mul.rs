// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use ppvm_traits::traits::{HashFinalize, PauliStorage, PauliWordTrait};

use crate::word::PauliWord;

impl<A, S> std::ops::Mul for PauliWord<A, S>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    type Output = PauliWord<A, S>;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut output = self;
        output *= rhs;
        output
    }
}

impl<A, S> std::ops::MulAssign for PauliWord<A, S>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    fn mul_assign(&mut self, rhs: Self) {
        for i in 0..self.n_qubits() {
            let x_i = self.xbits[i] ^ rhs.xbits[i];
            let z_i = self.zbits[i] ^ rhs.zbits[i];
            self.xbits.set(i, x_i);
            self.zbits.set(i, z_i);
        }
        self.rehash();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pauli_word_mul() {
        let a: PauliWord<[u8; 3]> = "XXXX_YYYY_ZZZZ_IIII".into();
        let b: PauliWord<[u8; 3]> = "YZXI_YZXI_YZXI_YZXI".into();
        let c: PauliWord<[u8; 3]> = "ZYIX_IXZY_XIYZ_YZXI".into();
        assert_eq!(a * b, c);
    }
}
