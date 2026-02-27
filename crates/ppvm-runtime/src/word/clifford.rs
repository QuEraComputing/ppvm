use std::hash::BuildHasher;

use super::data::PauliWord;
use crate::traits::{Clifford, PauliStorage, PauliWordTrait};

impl<A, H> Clifford for PauliWord<A, H>
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
        let index_x = self.xbits[index];
        let index_z = self.zbits[index];
        self.xbits.set(index, index_z);
        self.zbits.set(index, index_x);
        self.rehash();
    }
    fn s(&mut self, index: usize) {
        // S * I * S = I    00 -> 00, 0
        // S * X * S = -Y    10 -> 11, 0
        // S * Z * S = Z    01 -> 01, 0
        // S * Y * S = X   11 -> 10, 1
        let z = self.xbits[index] ^ self.zbits[index];
        self.zbits.set(index, z);
        self.rehash();
    }
    fn s_dagger(&mut self, index: usize) {
        // only adds different phase
        self.s(index);
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

        // flip the control z if target x is 1
        let control_z = self.zbits[control] ^ self.xbits[target];
        self.zbits.set(control, control_z);
        // flip the target z if control x is 1
        let target_z = self.zbits[target] ^ self.xbits[control];
        self.zbits.set(target, target_z);
        self.rehash();
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
        ] {
            let mut output: PauliWord<u64> = PauliWord::from(input);
            output.cnot(0, 1);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }

    #[test]
    fn test_h() {
        // NOTE: phase on "Y" not added in words
        for (input, target) in [("I", "I"), ("X", "Z"), ("Y", "Y"), ("Z", "X")] {
            let mut output: PauliWord<u64> = PauliWord::from(input);
            output.h(0);
            assert_eq!((input, output.to_string()), (input, target.to_string()));
        }
    }
}
