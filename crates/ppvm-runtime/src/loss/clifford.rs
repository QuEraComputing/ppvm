use std::hash::BuildHasher;

use super::data::LossyPauliWord;
use crate::traits::{Clifford, CliffordExtensions, PauliStorage, PauliWordTrait};

impl<A, H> Clifford for LossyPauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Clone + Default,
{
    fn x(&mut self, _index: usize) {
        // X * I * X = I    00 -> 00, 0
        // X * X * X = X    10 -> 10, 0
        // X * Z * X = -Z   01 -> 01, 1
        // X * Y * X = -Y.  11 -> 11, 1
    }
    fn y(&mut self, _index: usize) {
        // Y * I * Y = I    00 -> 00, 0
        // Y * X * Y = -X   10 -> 10, 1
        // Y * Z * Y = -Z   01 -> 01, 1
        // Y * Y * Y = Y    11 -> 11, 0
    }
    fn z(&mut self, _index: usize) {
        // Z * I * Z = I    00 -> 00, 0
        // Z * X * Z = -X   10 -> 10, 1
        // Z * Z * Z = Z    01 -> 01, 0
        // Z * Y * Z = -Y   11 -> 11, 1
    }
    fn h(&mut self, index: usize) {
        // H * I * H = I    00 -> 00, 0
        // H * X * H = Z    10 -> 01, 0
        // H * Z * H = X    01 -> 10, 0
        // H * Y * H = -Y   11 -> 11, 1
        if self.lbits[index] {
            // Hadamard on loss gives loss
            return;
        }
        let index_x = self.xbits[index];
        let index_z = self.zbits[index];
        self.xbits.set(index, index_z);
        self.zbits.set(index, index_x);
        self.rehash();
    }
    fn s(&mut self, index: usize) {
        // S * I * S = I    00 -> 00, 0
        // S * X * S = Y    10 -> 11, 0
        // S * Z * S = Z    01 -> 01, 0
        // S * Y * S = -X   11 -> 10, 1
        if self.lbits[index] {
            // S on loss gives loss
            return;
        }
        let z = self.xbits[index] ^ self.zbits[index];
        self.zbits.set(index, z);
        self.rehash();
    }
    fn cnot(&mut self, control: usize, target: usize) {
        //                          xx zz    xx zz  phase
        // CNOT * II * CNOT == II,  00 00 -> 00 00, 0
        // CNOT * IX * CNOT == IX,  01 00 -> 01 00, 0
        // CNOT * IZ * CNOT == ZZ,  00 01 -> 00 11, 0
        // CNOT * IY * CNOT == ZY,  01 01 -> 01 11, 0

        // CNOT * XI * CNOT == XX,  10 00 -> 11 00, 0
        // CNOT * XX * CNOT == XI,  11 00 -> 10 00, 0
        // CNOT * XY * CNOT == YZ,  11 01 -> 10 11, 0
        // CNOT * XZ * CNOT == -YY, 10 01 -> 11 11, 1

        // CNOT * ZI * CNOT == ZI,  00 10 -> 00 10, 0
        // CNOT * ZX * CNOT == ZX,  01 10 -> 01 10, 0
        // CNOT * ZY * CNOT == IY,  01 11 -> 01 01, 0
        // CNOT * ZZ * CNOT == IZ,  00 11 -> 00 01, 0

        // CNOT * YI * CNOT == YX,  10 10 -> 11 10, 0
        // CNOT * YX * CNOT == YI,  11 10 -> 10 10, 0
        // CNOT * YY * CNOT == -XZ, 11 11 -> 10 01, 1
        // CNOT * YZ * CNOT == XY,  10 11 -> 11 01, 0
        if self.lbits[control] || self.lbits[target] {
            return;
        }
        let control_z = self.zbits[target] ^ self.zbits[control];
        let target_x = self.xbits[control] ^ self.xbits[target];
        self.zbits.set(control, control_z);
        self.xbits.set(target, target_x);
        self.rehash();
    }
    fn cz(&mut self, control: usize, target: usize) {
        // CZ = |0><0| I + |1><1| Z
        // CZ * II * CZ = II,   00 00 -> 00 00, 0
        // CZ * IX * CZ = ZX,   01 00 -> 01 10, 0
        // CZ * IY * CZ = ZY,   01 01 -> 01 11, 0
        // CZ * IZ * CZ = ZZ,   00 01 -> 00 01, 0

        // CZ * XI * CZ = XZ,   10 00 -> 10 01, 0
        // CZ * XX * CZ = YY,   11 00 -> 11 11, 0
        // CZ * XY * CZ = -YX,  11 01 -> 11 10, 1
        // CZ * XZ * CZ = XI,   10 01 -> 10 00, 0

        // CZ * ZI * CZ == ZI,  00 10 -> 00 10, 0
        // CZ * ZX * CZ == IX,  01 10 -> 01 00, 0
        // CZ * ZY * CZ == IY,  01 11 -> 01 01, 0
        // CZ * ZZ * CZ == ZZ,  00 11 -> 00 11, 0

        // CZ * YI * CZ == YZ,  10 10 -> 10 11, 0
        // CZ * YX * CZ == -XY, 11 10 -> 11 01, 1
        // CZ * YY * CZ == XX,  11 11 -> 11 00, 0
        // CZ * YZ * CZ == YI,  10 11 -> 10 10, 0

        // xx: identity
        // zz:
        // xx: 00, identity
        // xx: 01, 00 -> 10, 01 -> 11, 10 -> 00, 11 -> 01
        // xx: 10, 00 -> 01, 01 -> 00, 10 -> 11, 11 -> 10
        // xx: 11, 00 -> 11, 01 -> 10, 10 -> 01, 11 -> 00

        if self.lbits[control] || self.lbits[target] {
            return;
        }

        // flip the control z if target x is 1
        let control_z = self.zbits[control] ^ self.xbits[target];
        self.zbits.set(control, control_z);
        // flip the target z if control x is 1
        let target_z = self.zbits[target] ^ self.xbits[control];
        self.zbits.set(target, target_z);
        self.rehash();
    }
}

impl<A, H> CliffordExtensions for LossyPauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Clone + Default,
{
    // |    Gate    |  X  |  Y  |  Z  |
    // |:----------:|:---:|:---:|:---:|
    // |     s      | -Y  |  X  |  Z  |
    // |   s_adj    |  Y  | -X  |  Z  |
    // |   sqrt_x   |  X  | -Z  |  Y  |
    // | sqrt*x*adj |  X  |  Z  | -Y  |
    // |   sqrt_y   |  Z  |  Y  | -X  |
    // | sqrt*y*adj | -Z  |  Y  |  X  |

    #[inline]
    fn s_adj(&mut self, addr0: usize) {
        if self.lbits[addr0] {
            return;
        }
        self.s(addr0);
    }

    #[inline]
    fn sqrt_x(&mut self, addr0: usize) {
        if self.lbits[addr0] {
            return;
        }
        let x = self.xbits[addr0];
        let z = self.zbits[addr0];
        self.set_xbit(addr0, x ^ z);
        self.rehash();
    }

    #[inline]
    fn sqrt_y(&mut self, addr0: usize) {
        if self.lbits[addr0] {
            return;
        }
        let x = self.xbits[addr0];
        let z = self.zbits[addr0];
        self.set_xbit(addr0, z);
        self.set_zbit(addr0, x);
        self.rehash();
    }

    #[inline]
    fn sqrt_x_adj(&mut self, addr0: usize) {
        if self.lbits[addr0] {
            return;
        }
        let x = self.xbits[addr0];
        let z = self.zbits[addr0];
        self.set_xbit(addr0, x ^ z);
        self.rehash();
    }

    #[inline]
    fn sqrt_y_adj(&mut self, addr0: usize) {
        if self.lbits[addr0] {
            return;
        }
        let x = self.xbits[addr0];
        let z = self.zbits[addr0];
        self.set_xbit(addr0, z);
        self.set_zbit(addr0, x);
        self.rehash();
    }

    // | CY  |  I  |  X  |  Y  |  Z  |
    // |:---:|:---:|:---:|:---:|:---:|
    // |  I  | II  | ZX  | IY  | ZZ  |
    // |  X  | XY  | -YZ | XI  | YX  |
    // |  Y  | YY  | XZ  | YI  | -XX |
    // |  Z  | ZI  | IX  | ZY  | IZ  |

    #[inline]
    fn cy(&mut self, addr0: usize, addr1: usize) {
        if self.lbits[addr0] || self.lbits[addr1] {
            return;
        }
        let xc = self.xbits[addr0];
        let zc = self.zbits[addr0];
        let xt = self.xbits[addr1];
        let zt = self.zbits[addr1];
        self.set_zbit(addr0, zc ^ xt ^ zt);
        self.set_xbit(addr1, xt ^ xc);
        self.set_zbit(addr1, zt ^ xc);
        self.rehash();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_x() {
        for (input, target) in [("I", "I"), ("X", "X"), ("Y", "Y"), ("Z", "Z"), ("L", "L")] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.x(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_y() {
        for (input, target) in [("I", "I"), ("X", "X"), ("Y", "Y"), ("Z", "Z"), ("L", "L")] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.y(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_z() {
        for (input, target) in [("I", "I"), ("X", "X"), ("Y", "Y"), ("Z", "Z"), ("L", "L")] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.z(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    // CNOT * II * CNOT == II,  00 00 -> 00 00, 0
    // CNOT * IX * CNOT == IX,  01 00 -> 01 00, 0
    // CNOT * IZ * CNOT == ZZ,  00 01 -> 00 11, 0
    // CNOT * IY * CNOT == ZY,  01 01 -> 01 11, 0

    // CNOT * XI * CNOT == XX,  10 00 -> 11 00, 0
    // CNOT * XX * CNOT == XI,  11 00 -> 10 00, 0
    // CNOT * XY * CNOT == YZ,  11 01 -> 10 11, 0
    // CNOT * XZ * CNOT == -YY, 10 01 -> 11 11, 1

    // CNOT * ZI * CNOT == ZI,  00 10 -> 00 10, 0
    // CNOT * ZX * CNOT == ZX,  01 10 -> 01 10, 0
    // CNOT * ZY * CNOT == IY,  01 11 -> 01 01, 0
    // CNOT * ZZ * CNOT == IZ,  00 11 -> 00 01, 0

    // CNOT * YI * CNOT == YX,  10 10 -> 11 10, 0
    // CNOT * YX * CNOT == YI,  11 10 -> 10 10, 0
    // CNOT * YY * CNOT == -XZ, 11 11 -> 10 01, 1
    // CNOT * YZ * CNOT == XY,  10 11 -> 11 01, 0

    #[test]
    fn test_cnot() {
        for (input, target) in [
            ("II", "II"),
            ("IX", "IX"),
            ("IZ", "ZZ"),
            ("IY", "ZY"),
            ("XI", "XX"),
            ("XX", "XI"),
            ("XY", "YZ"),
            ("XZ", "YY"),
            ("ZI", "ZI"),
            ("ZX", "ZX"),
            ("ZY", "IY"),
            ("ZZ", "IZ"),
            ("YI", "YX"),
            ("YX", "YI"),
            ("YY", "XZ"),
            ("YZ", "XY"),
            ("IL", "IL"),
            ("XL", "XL"),
            ("YL", "YL"),
            ("ZL", "ZL"),
            ("LI", "LI"),
            ("LX", "LX"),
            ("LY", "LY"),
            ("LZ", "LZ"),
            ("LL", "LL"),
        ] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.cnot(0, 1);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_h() {
        // NOTE: phase on "Y" not added in words
        for (input, target) in [("I", "I"), ("X", "Z"), ("Y", "Y"), ("Z", "X"), ("L", "L")] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.h(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_s() {
        for (input, target) in [("I", "I"), ("X", "Y"), ("Z", "Z"), ("Y", "X"), ("L", "L")] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.s(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_s_adj() {
        for (input, target) in [("I", "I"), ("X", "Y"), ("Z", "Z"), ("Y", "X"), ("L", "L")] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.s_adj(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_sqrt_x() {
        for (input, target) in [("I", "I"), ("X", "X"), ("Y", "Z"), ("Z", "Y"), ("L", "L")] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.sqrt_x(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_sqrt_x_adj() {
        for (input, target) in [("I", "I"), ("X", "X"), ("Y", "Z"), ("Z", "Y"), ("L", "L")] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.sqrt_x_adj(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_sqrt_y() {
        for (input, target) in [("I", "I"), ("X", "Z"), ("Y", "Y"), ("Z", "X"), ("L", "L")] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.sqrt_y(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_sqrt_y_adj() {
        for (input, target) in [("I", "I"), ("X", "Z"), ("Y", "Y"), ("Z", "X"), ("L", "L")] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.sqrt_y_adj(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_cz() {
        for (input, target) in [
            ("II", "II"),
            ("IX", "ZX"),
            ("IY", "ZY"),
            ("IZ", "IZ"),
            ("XI", "XZ"),
            ("XX", "YY"),
            ("XY", "YX"),
            ("XZ", "XI"),
            ("ZI", "ZI"),
            ("ZX", "IX"),
            ("ZY", "IY"),
            ("ZZ", "ZZ"),
            ("YI", "YZ"),
            ("YX", "XY"),
            ("YY", "XX"),
            ("YZ", "YI"),
            ("IL", "IL"),
            ("XL", "XL"),
            ("YL", "YL"),
            ("ZL", "ZL"),
            ("LI", "LI"),
            ("LX", "LX"),
            ("LY", "LY"),
            ("LZ", "LZ"),
            ("LL", "LL"),
        ] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.cz(0, 1);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_cy() {
        for (input, target) in [
            ("II", "II"),
            ("IX", "ZX"),
            ("IZ", "ZZ"),
            ("IY", "IY"),
            ("XI", "XY"),
            ("XX", "YZ"),
            ("XY", "XI"),
            ("XZ", "YX"),
            ("ZI", "ZI"),
            ("ZX", "IX"),
            ("ZY", "ZY"),
            ("ZZ", "IZ"),
            ("YI", "YY"),
            ("YX", "XZ"),
            ("YY", "YI"),
            ("YZ", "XX"),
            ("IL", "IL"),
            ("XL", "XL"),
            ("YL", "YL"),
            ("ZL", "ZL"),
            ("LI", "LI"),
            ("LX", "LX"),
            ("LY", "LY"),
            ("LZ", "LZ"),
            ("LL", "LL"),
        ] {
            let mut output: LossyPauliWord<u64> = LossyPauliWord::from(input);
            output.cy(0, 1);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }
}
