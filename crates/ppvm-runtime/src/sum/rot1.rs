use crate::traits::*;
use crate::word::PauliWord;
use crate::{config::Config, sum::PauliSum};

impl<T, S, H> RotationOne<T> for PauliSum<T>
where
    S: PauliStorage,
    H: std::hash::BuildHasher + Clone + Default,
    T: Config<Storage = S, BuildHasher = H, PauliWordType = PauliWord<S, H>>,
    T::Coeff: std::ops::MulAssign,
    T::Map: ACMapInsert<S, T::Coeff, H, PauliWord<S, H>> + ACMapConsume,
{
    fn rotate_1(&mut self, axis: crate::char::Pauli, addr0: usize, theta: <T as Config>::Coeff) {
        let (sin, cos) = theta.sin_cos();
        self.map_insert(|k, v| {
            let p_g = k.get(addr0);
            let (eps, p_q) = levi_civita(p_g as u8, axis as u8);
            if eps == 0 {
                return None;
            } else {
                let mut coeff = v.clone();
                *v *= cos.clone();
                let mut new_word = k.clone();
                new_word.xbits.set(addr0, p_q & 0b01 != 0);
                new_word.zbits.set(addr0, p_q & 0b10 != 0);
                new_word.rehash();

                coeff *= sin.mul_sign(eps);
                return Some((new_word, coeff));
            }
        });
    }
}

/// 2-bit Pauli code: 00 I, 01 X, 10 Z, 11 Y
/// Returns \(ε, k\) so that  –i \[P_i, P_j\]/2 = ε · P_k.
/// For every commuting pair it yields (0, 0).
#[inline]
pub fn levi_civita(i: u8, j: u8) -> (i8, u8) {
    // --------------------------------------------------- third Pauli by XOR
    let k = i ^ j; // 0 when i == j

    // ----------- commute ⇔ i==0  OR  j==0  OR  k==0  (no false positives)
    let commute = ((i == 0) | (j == 0) | (k == 0)) as u8; // 1 = commute

    // ------------------------------------------------------ sign ε_{ijk}
    #[inline]
    fn rank(p: u8) -> u8 {
        let b1 = p >> 1; // MSB
        (b1 << 1).wrapping_sub(b1 & (p & 1)) // 0,1,2 for X,Y,Z
    }

    let ri = rank(i);
    let rj = rank(j);

    // diff = (rj - ri) mod 3   without an actual modulus
    let mut diff = rj.wrapping_sub(ri).wrapping_add(3); // 0…5
    diff -= 3 & (0u8.wrapping_sub(diff >> 2)); // if ≥3 subtract 3

    // +1 when diff == 1,  –1 when diff == 2
    let eps_raw = 1i8 - 2 * ((diff >> 1) as i8);

    // --------------------------------------------------- zero when commute
    let eps = eps_raw * (1 - commute as i8); // 0 if commute
    let k = k * (1 - commute); // 0 if commute

    (eps, k)
}

#[cfg(test)]
mod tests {
    use crate::config::fxhash::ByteF64;

    use super::*;

    #[test]
    fn test_rx() {
        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("X", 1.0);
        answer.rx(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("X", 1.0);
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Y", 1.0);
        answer.rx(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Y", 2.1_f64.cos());
        expect += ("Z", -2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Z", 1.0);
        answer.rx(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Z", 2.1_f64.cos());
        expect += ("Y", 2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("I", 1.0);
        answer.rx(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("I", 1.0);
        assert_eq!(answer, expect);
    }

    #[test]
    fn test_ry() {
        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("X", 1.0);
        answer.ry(0, 2.1);

        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("X", 2.1_f64.cos());
        expect += ("Z", 2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Y", 1.0);
        answer.ry(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Y", 1.0);
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Z", 1.0);
        answer.ry(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Z", 2.1_f64.cos());
        expect += ("X", -2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("I", 1.0);
        answer.ry(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("I", 1.0);
        assert_eq!(answer, expect);
    }

    #[test]
    fn test_rz() {
        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("X", 1.0);
        answer.rz(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("X", 2.1_f64.cos());
        expect += ("Y", -2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Y", 1.0);
        answer.rz(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Y", 2.1_f64.cos());
        expect += ("X", 2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Z", 1.0);
        answer.rz(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Z", 1.0);
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("I", 1.0);
        answer.rz(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("I", 1.0);
        assert_eq!(answer, expect);
    }
}
