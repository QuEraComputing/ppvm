// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::BitView;
use num::PrimInt;

use ppvm_traits::traits::{HashFinalize, PauliStorage, PauliWordTrait};

use crate::{phase::PhasedPauliWord, word::PauliWord};

impl<A, S, const REHASH: bool> std::ops::MulAssign
    for PhasedPauliWord<A, S, PauliWord<A, S, REHASH>>
where
    A: PauliStorage,
    <A as BitView>::Store: PrimInt,
    S: std::hash::BuildHasher + Clone + Default + HashFinalize,
{
    fn mul_assign(&mut self, rhs: Self) {
        *self *= &rhs;
    }
}

impl<A, S, const REHASH: bool> std::ops::MulAssign<&Self>
    for PhasedPauliWord<A, S, PauliWord<A, S, REHASH>>
where
    A: PauliStorage,
    <A as BitView>::Store: PrimInt,
    S: std::hash::BuildHasher + Clone + Default + HashFinalize,
{
    fn mul_assign(&mut self, rhs: &Self) {
        let mut sign_count = 0u32;
        let mut imag_count = 0u32;
        let lhs_x = &mut self.word.xbits.data;
        let lhs_z = &mut self.word.zbits.data;
        let rhs_x = &rhs.word.xbits.data;
        let rhs_z = &rhs.word.zbits.data;
        for i in 0..lhs_x.as_raw_slice().len() {
            let a = lhs_x.as_raw_slice()[i];
            let b = lhs_z.as_raw_slice()[i];
            let c = rhs_x.as_raw_slice()[i];
            let d = rhs_z.as_raw_slice()[i];
            let sign = (a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d);
            let imag = (a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d);
            sign_count += sign.count_ones();
            imag_count += imag.count_ones();
            lhs_x.as_raw_mut_slice()[i] = a ^ c;
            lhs_z.as_raw_mut_slice()[i] = b ^ d;
        }
        self.add_phase(((2 * sign_count + imag_count) % 4) as u8);
        self.word.rehash();
        self.add_phase(rhs.phase);
    }
}

impl<A, S, const REHASH: bool> std::ops::Mul for PhasedPauliWord<A, S, PauliWord<A, S, REHASH>>
where
    A: PauliStorage + Clone,
    <A as BitView>::Store: PrimInt,
    S: std::hash::BuildHasher + Clone + Default + HashFinalize,
{
    // xz xz phase
    // 00 00 00
    // 00 01 00
    // 00 10 00
    // 00 11 00
    //
    // 01 00 00
    // 01 01 00
    // 01 10 01
    // 01 11 11
    //
    // 10 00 00
    // 10 01 11
    // 10 10 00
    // 10 11 01
    //
    // 11 00 00
    // 11 01 01
    // 11 10 11
    // 11 11 00
    type Output = PhasedPauliWord<A, S, PauliWord<A, S, REHASH>>;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut output = self.clone();
        output *= rhs;
        output
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mul() {
        for (lhs, rhs, ans) in [("+X", "+X", "+I"), ("+X", "+Y", "+iZ"), ("+X", "+Z", "-iY")] {
            let x: PhasedPauliWord<u64> = lhs.into();
            let y: PhasedPauliWord<u64> = rhs.into();
            assert_eq!((x * y).to_string(), ans);
        }
    }

    #[test]
    fn test_mul_multi_qubit() {
        for (lhs, rhs, ans) in [
            ("+ZI", "-ZI", "-II"),
            ("+II", "-ZI", "-ZI"),
            ("+XI", "+iXI", "+iII"),
            ("-XX", "-XX", "+II"),
        ] {
            let x: PhasedPauliWord<u64> = lhs.into();
            let y: PhasedPauliWord<u64> = rhs.into();
            assert_eq!((x * y).to_string(), ans);
        }
    }
}
