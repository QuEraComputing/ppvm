use std::hash::BuildHasher;

use super::data::PhasedPauliWord;
use crate::traits::Clifford;
use crate::traits::PauliStorage;

impl<S, H> Clifford for PhasedPauliWord<S, H>
where
    S: PauliStorage,
    H: BuildHasher + Clone + Default,
{
    fn x(&mut self, index: usize) {
        let phase = (self.word.zbits[index]) as u8;
        self.word.x(index);
        self.add_phase(phase << 1);
    }
    fn y(&mut self, index: usize) {
        let phase = (self.word.xbits[index] ^ self.word.zbits[index]) as u8;
        self.word.y(index);
        self.add_phase(phase << 1);
    }
    fn z(&mut self, index: usize) {
        let phase = (self.word.xbits[index]) as u8;
        self.word.z(index);
        self.add_phase(phase << 1);
    }
    fn h(&mut self, index: usize) {
        let phase = (self.word.xbits[index] & self.word.zbits[index]) as u8;
        self.word.h(index);
        self.add_phase(phase << 1);
    }
    fn s(&mut self, index: usize) {
        let phase = (self.word.xbits[index] & !self.word.zbits[index]) as u8;
        self.word.s(index);
        self.add_phase(phase << 1);
    }
    fn s_dagger(&mut self, index: usize) {
        let phase = (self.word.xbits[index] & self.word.zbits[index]) as u8;
        self.word.s_dagger(index);
        self.add_phase(phase << 1);
    }
    fn cnot(&mut self, control: usize, target: usize) {
        // phase = 1x y1 where x xor y = 0
        // xx zz    xx zz
        // 11 11 -> 10 01, 2
        // 10 01 -> 11 11, 2
        let phase = ((self.word.xbits[control] & self.word.zbits[target])
            & (self.word.xbits[target] == self.word.zbits[control])) as u8;
        self.word.cnot(control, target);
        self.add_phase(phase << 1);
    }
    fn cz(&mut self, control: usize, target: usize) {
        // phase = 11 10, 11 01 = 11 ab where a ^ b = 1
        // 11 01 -> 11 10, 2
        // 11 10 -> 11 01, 2
        let phase = ((self.word.xbits[control] & self.word.xbits[target])
            & (self.word.zbits[control] ^ self.word.zbits[target])) as u8;
        self.word.cz(control, target);
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
    fn test_s_dagger() {
        for (input, target) in [("+I", "+I"), ("+X", "+Y"), ("+Y", "-X"), ("+Z", "+Z")] {
            let mut output: PhasedPauliWord<u64> = PhasedPauliWord::from(input);
            output.s_dagger(0);
            assert_eq!(output.to_string(), target.to_string());
        }
    }
}
