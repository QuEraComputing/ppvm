// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

// Clifford behavior for `LossyPauliWord` is provided by the blanket impl
// `impl<T: PauliWordTrait> Clifford for T` in `crate::traits::clifford`.
// `LossyPauliWord::get_lbit` reports actual loss bits, so the gates correctly
// no-op on lost qubits.

#[cfg(test)]
mod tests {
    use super::super::data::LossyPauliWord;
    use crate::traits::{Clifford, CliffordExtensions};

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
