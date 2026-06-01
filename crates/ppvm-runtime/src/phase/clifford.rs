// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use super::data::PhasedPauliWord;
use crate::traits::{Clifford, CliffordExtensions};
use crate::traits::{PauliStorage, PauliWordTrait};

impl<S, H, W> Clifford for PhasedPauliWord<S, H, W>
where
    S: PauliStorage,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait + Clifford,
{
    #[inline]
    fn x(&mut self, index: usize) {
        if self.word.get_lbit(index) {
            // check for loss
            return;
        }
        let phase = (self.word.get_zbit(index)) as u8;
        self.word.x(index);
        self.add_phase(phase << 1);
    }

    #[inline]
    fn y(&mut self, index: usize) {
        if self.word.get_lbit(index) {
            // check for loss
            return;
        }
        let phase = (self.word.get_xbit(index) ^ self.word.get_zbit(index)) as u8;
        self.word.y(index);
        self.add_phase(phase << 1);
    }

    #[inline]
    fn z(&mut self, index: usize) {
        if self.word.get_lbit(index) {
            // check for loss
            return;
        }
        let phase = (self.word.get_xbit(index)) as u8;
        self.word.z(index);
        self.add_phase(phase << 1);
    }

    #[inline]
    fn h(&mut self, index: usize) {
        if self.word.get_lbit(index) {
            // check for loss
            return;
        }
        let phase = (self.word.get_xbit(index) & self.word.get_zbit(index)) as u8;
        self.word.h(index);
        self.add_phase(phase << 1);
    }

    #[inline]
    fn s(&mut self, index: usize) {
        if self.word.get_lbit(index) {
            // check for loss
            return;
        }
        let phase = (self.word.get_xbit(index) & !self.word.get_zbit(index)) as u8;
        self.word.s(index);
        self.add_phase(phase << 1);
    }

    #[inline]
    fn cnot(&mut self, control: usize, target: usize) {
        // phase = 1x y1 where x xor y = 0
        // xx zz    xx zz
        // 11 11 -> 10 01, 2
        // 10 01 -> 11 11, 2
        if self.word.get_lbit(control) || self.word.get_lbit(target) {
            return;
        }
        let phase = ((self.word.get_xbit(control) & self.word.get_zbit(target))
            & (self.word.get_xbit(target) == self.word.get_zbit(control)))
            as u8;
        self.word.cnot(control, target);
        self.add_phase(phase << 1);
    }

    #[inline]
    fn cz(&mut self, control: usize, target: usize) {
        // phase = 11 10, 11 01 = 11 ab where a ^ b = 1
        // 11 01 -> 11 10, 2
        // 11 10 -> 11 01, 2
        if self.word.get_lbit(control) || self.word.get_lbit(target) {
            return;
        }
        let phase = ((self.word.get_xbit(control) & self.word.get_xbit(target))
            & (self.word.get_zbit(control) ^ self.word.get_zbit(target))) as u8;
        self.word.cz(control, target);
        self.add_phase(phase << 1);
    }
}

impl<S, H, W> CliffordExtensions for PhasedPauliWord<S, H, W>
where
    S: PauliStorage,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait + Clifford + CliffordExtensions,
{
    // |    Gate    |  X  |  Y  |  Z  |
    // |:----------:|:---:|:---:|:---:|
    // |     s      | -Y  |  X  |  Z  |
    // |   s_adj    |  Y  | -X  |  Z  |
    // |   sqrt_x   |  X  | -Z  |  Y  |
    // | sqrt*x*adj |  X  |  Z  | -Y  |
    // |   sqrt_y   |  Z  |  Y  | -X  |
    // | sqrt*y*adj | -Z  |  Y  |  X  |

    fn s_adj(&mut self, addr0: usize) {
        if self.word.get_lbit(addr0) {
            return;
        }
        let phase = (self.word.get_xbit(addr0) & self.word.get_zbit(addr0)) as u8;
        self.word.s_adj(addr0);
        self.add_phase(phase << 1);
    }

    fn sqrt_x(&mut self, addr0: usize) {
        if self.word.get_lbit(addr0) {
            return;
        }
        let phase = (self.word.get_xbit(addr0) & self.word.get_zbit(addr0)) as u8;
        self.word.sqrt_x(addr0);
        self.add_phase(phase << 1);
    }

    fn sqrt_y(&mut self, addr0: usize) {
        if self.word.get_lbit(addr0) {
            return;
        }
        let phase = (!self.word.get_xbit(addr0) & self.word.get_zbit(addr0)) as u8;
        self.word.sqrt_y(addr0);
        self.add_phase(phase << 1);
    }

    fn sqrt_x_adj(&mut self, addr0: usize) {
        if self.word.get_lbit(addr0) {
            return;
        }
        let phase = (!self.word.get_xbit(addr0) & self.word.get_zbit(addr0)) as u8;
        self.word.sqrt_x_adj(addr0);
        self.add_phase(phase << 1);
    }

    fn sqrt_y_adj(&mut self, addr0: usize) {
        if self.word.get_lbit(addr0) {
            return;
        }
        let phase = (self.word.get_xbit(addr0) & !self.word.get_zbit(addr0)) as u8;
        self.word.sqrt_y_adj(addr0);
        self.add_phase(phase << 1);
    }

    // | CY  |  I  |  X  |  Y  |  Z  |
    // |:---:|:---:|:---:|:---:|:---:|
    // |  I  | II  | ZX  | IY  | ZZ  |
    // |  X  | XY  | -YZ | XI  | YX  |
    // |  Y  | YY  | XZ  | YI  | -XX |
    // |  Z  | ZI  | IX  | ZY  | IZ  |
    fn cy(&mut self, addr0: usize, addr1: usize) {
        if self.word.get_lbit(addr0) || self.word.get_lbit(addr1) {
            return;
        }
        // phase = -1 for XX -> -YZ and YZ -> -XX
        let xc = self.word.get_xbit(addr0);
        let zc = self.word.get_zbit(addr0);
        let xt = self.word.get_xbit(addr1);
        let zt = self.word.get_zbit(addr1);
        let phase = (xc & (xt ^ zt) & !(zc ^ zt)) as u8;
        self.word.cy(addr0, addr1);
        self.add_phase(phase << 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // CNOT * II * CNOT == II,  00 00 -> 00 00, 0
    // CNOT * IX * CNOT == IX,  01 00 -> 01 00, 0
    // CNOT * IZ * CNOT == ZZ,  00 01 -> 00 11, 0
    // CNOT * IY * CNOT == ZY,  01 01 -> 01 11, 0

    // CNOT * XI * CNOT == XX,  10 00 -> 11 00, 0
    // CNOT * XX * CNOT == XI,  11 00 -> 10 00, 0
    // CNOT * XY * CNOT == YZ,  11 01 -> 10 11, 0
    // CNOT * XZ * CNOT == -YY, 10 01 -> 11 11, 2

    // CNOT * ZI * CNOT == ZI,  00 10 -> 00 10, 0
    // CNOT * ZX * CNOT == ZX,  01 10 -> 01 10, 0
    // CNOT * ZY * CNOT == IY,  01 11 -> 01 01, 0
    // CNOT * ZZ * CNOT == IZ,  00 11 -> 00 01, 0

    // CNOT * YI * CNOT == YX,  10 10 -> 11 10, 0
    // CNOT * YX * CNOT == YI,  11 10 -> 10 10, 0
    // CNOT * YY * CNOT == -XZ, 11 11 -> 10 01, 2
    // CNOT * YZ * CNOT == XY,  10 11 -> 11 01, 0

    #[test]
    fn test_cnot() {
        for (input, target) in [
            ("+II", "+II"),
            ("+IX", "+IX"),
            ("+IZ", "+ZZ"),
            ("+IY", "+ZY"),
            ("+XI", "+XX"),
            ("+XX", "+XI"),
            ("+XY", "+YZ"),
            ("+XZ", "-YY"),
            ("+ZI", "+ZI"),
            ("+ZX", "+ZX"),
            ("+ZY", "+IY"),
            ("+ZZ", "+IZ"),
            ("+YI", "+YX"),
            ("+YX", "+YI"),
            ("+YY", "-XZ"),
            ("+YZ", "+XY"),
        ] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.cnot(0, 1);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_x() {
        for (input, target) in [("+I", "+I"), ("+X", "+X"), ("+Y", "-Y"), ("+Z", "-Z")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.x(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }

    #[test]
    fn test_y() {
        for (input, target) in [("+I", "+I"), ("+X", "-X"), ("+Y", "+Y"), ("+Z", "-Z")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.y(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }

    #[test]
    fn test_z() {
        for (input, target) in [("+I", "+I"), ("+X", "-X"), ("+Y", "-Y"), ("+Z", "+Z")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.z(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }

    #[test]
    fn test_h() {
        for (input, target) in [("+I", "+I"), ("+X", "+Z"), ("+Y", "-Y"), ("+Z", "+X")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.h(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }

    #[test]
    fn test_s() {
        for (input, target) in [("+I", "+I"), ("+X", "-Y"), ("+Y", "+X"), ("+Z", "+Z")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.s(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }

    #[test]
    fn test_s_adj() {
        for (input, target) in [("+I", "+I"), ("+X", "+Y"), ("+Y", "-X"), ("+Z", "+Z")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.s_adj(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }

    #[test]
    fn test_sqrt_x() {
        for (input, target) in [("+I", "+I"), ("+X", "+X"), ("+Y", "-Z"), ("+Z", "+Y")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.sqrt_x(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }

    #[test]
    fn test_sqrt_x_adj() {
        for (input, target) in [("+I", "+I"), ("+X", "+X"), ("+Y", "+Z"), ("+Z", "-Y")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.sqrt_x_adj(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }

    #[test]
    fn test_sqrt_y() {
        for (input, target) in [("+I", "+I"), ("+X", "+Z"), ("+Y", "+Y"), ("+Z", "-X")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.sqrt_y(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }

    #[test]
    fn test_sqrt_y_adj() {
        for (input, target) in [("+I", "+I"), ("+X", "-Z"), ("+Y", "+Y"), ("+Z", "+X")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.sqrt_y_adj(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }

    #[test]
    fn test_cz() {
        for (input, target) in [
            ("+II", "+II"),
            ("+IX", "+ZX"),
            ("+IY", "+ZY"),
            ("+IZ", "+IZ"),
            ("+XI", "+XZ"),
            ("+XX", "+YY"),
            ("+XY", "-YX"),
            ("+XZ", "+XI"),
            ("+ZI", "+ZI"),
            ("+ZX", "+IX"),
            ("+ZY", "+IY"),
            ("+ZZ", "+ZZ"),
            ("+YI", "+YZ"),
            ("+YX", "-XY"),
            ("+YY", "+XX"),
            ("+YZ", "+YI"),
        ] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.cz(0, 1);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_cy() {
        for (input, target) in [
            ("+II", "+II"),
            ("+IX", "+ZX"),
            ("+IY", "+IY"),
            ("+IZ", "+ZZ"),
            ("+XI", "+XY"),
            ("+XX", "-YZ"),
            ("+XY", "+XI"),
            ("+XZ", "+YX"),
            ("+ZI", "+ZI"),
            ("+ZX", "+IX"),
            ("+ZY", "+ZY"),
            ("+ZZ", "+IZ"),
            ("+YI", "+YY"),
            ("+YX", "+XZ"),
            ("+YY", "+YI"),
            ("+YZ", "-XX"),
        ] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.cy(0, 1);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }
}
