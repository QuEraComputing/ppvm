use std::hash::BuildHasher;

use crate::{traits::PauliStorage, word::PauliWord};

impl<A, S> std::ops::Mul for PauliWord<A, S>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default,
{
    type Output = PauliWord<A, S>;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut output = self.clone();
        output *= rhs;
        output
    }
}

impl<A, S> std::ops::MulAssign for PauliWord<A, S>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default,
{
    fn mul_assign(&mut self, rhs: Self) {
        for i in 0..self.n_qubits() {
            // TODO: this should actually return 0 if there are L bits in either
            if self.lbits[i] || rhs.lbits[i] {
                panic!("Multiplication involving L bits is not defined");
            }
            let l_i = self.lbits[i] & rhs.lbits[i];
            let x_i = (self.xbits[i] ^ rhs.xbits[i]) & !l_i;
            let z_i = (self.zbits[i] ^ rhs.zbits[i]) & !l_i;
            self.xbits.set(i, x_i);
            self.zbits.set(i, z_i);
            self.lbits.set(i, l_i);
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
