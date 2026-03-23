use crate::{phase::PhasedPauliWord, traits::{PauliStorage, PauliWordTrait}};

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
        // Rehash after modifying xbits/zbits directly (bypassing PauliWord::MulAssign
        // which calls rehash internally). Without this, hash_cache is stale and the
        // word cannot be used as a correct HashMap key.
        self.word.rehash();
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
    use std::collections::HashMap;

    /// Verify that after MulAssign the word's hash is consistent with its bits.
    ///
    /// If hash_cache were stale, two products that yield the same word would have
    /// different hashes and a HashMap would treat them as distinct keys instead of
    /// accumulating their values — the regression this test catches.
    #[test]
    fn mul_assign_rehashes_word() {
        use crate::word::PauliWord;
        use std::hash::{BuildHasher, Hash, Hasher};
        type W1 = PauliWord<[u8; 1]>;
        type PPW1 = PhasedPauliWord<[u8; 1]>;
        // X*X = I and Z*Z = I. Start from two different words so the stale
        // hashes (hash_cache of the initial lhs) would differ without the fix.
        let mut a = PPW1::build_from_word(W1::from("X"), 0);
        let b = PPW1::build_from_word(W1::from("X"), 0);
        let mut c = PPW1::build_from_word(W1::from("Z"), 0);
        let d = PPW1::build_from_word(W1::from("Z"), 0);
        a *= b; // a.word = I, derived from X
        c *= d; // c.word = I, derived from Z
        // Both products must equal the identity word.
        assert_eq!(a.word, c.word, "both products must equal I");
        // Hash must be consistent with the bits so HashMap accumulates them.
        let bh = fxhash::FxBuildHasher::default();
        let hash_a = { let mut h = bh.build_hasher(); a.word.hash(&mut h); h.finish() };
        let hash_c = { let mut h = bh.build_hasher(); c.word.hash(&mut h); h.finish() };
        assert_eq!(hash_a, hash_c, "equal words must have equal hashes after MulAssign");
        // Functional check: both should map to the same HashMap entry.
        let mut map: HashMap<W1, i32> = HashMap::new();
        *map.entry(a.word).or_insert(0) += 1;
        *map.entry(c.word).or_insert(0) += 1;
        assert_eq!(map.len(), 1, "equal words must map to the same entry");
        assert_eq!(*map.values().next().unwrap(), 2);
    }

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
