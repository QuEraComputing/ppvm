use crate::{phase::PhasedPauliWord, traits::PauliStorage};

impl<A, S> std::ops::MulAssign for PhasedPauliWord<A, S>
where
    A: PauliStorage,
    S: std::hash::BuildHasher + Clone + Default,
{
    fn mul_assign(&mut self, rhs: Self) {
        *self *= &rhs;
    }
}

impl<A, S> std::ops::MulAssign<&Self> for PhasedPauliWord<A, S>
where
    A: PauliStorage,
    S: std::hash::BuildHasher + Clone + Default,
{
    fn mul_assign(&mut self, rhs: &Self) {
        for i in 0..self.n_qubits() {
            let x_i = self.word.xbits[i] ^ rhs.word.xbits[i];
            let z_i = self.word.zbits[i] ^ rhs.word.zbits[i];
            let a = self.word.xbits[i];
            let b = self.word.zbits[i];
            let c = rhs.word.xbits[i];
            let d = rhs.word.zbits[i];
            let sign = (a && b && c && !d) || (a && !b && !c && d) || (!a && b && c && d);
            let imag = (a && !b && d) || (a && !c && d) || (!a && b && c) || (b && c && !d);
            let exp = (sign as u8) << 1 | (imag as u8);
            self.add_phase(exp);
            self.word.xbits.set(i, x_i);
            self.word.zbits.set(i, z_i);
        }
        self.add_phase(rhs.phase);
    }
}

impl<A, S> std::ops::Mul for PhasedPauliWord<A, S>
where
    A: PauliStorage + Clone,
    S: std::hash::BuildHasher + Clone + Default,
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
    type Output = PhasedPauliWord<A, S>;
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
