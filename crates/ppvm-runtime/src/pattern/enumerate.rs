use super::data::{Decorated, OpPattern, PauliPattern};
use crate::char::Pauli;
use crate::traits::PauliStorage;
use crate::word::PauliWord;

use itertools::{Itertools, MultiProduct};

impl OpPattern {
    pub fn enumerate_matches(&self) -> EnumMatchesOpPattern<'_> {
        EnumMatchesOpPattern {
            pattern: self,
            current: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct EnumMatchesOpPattern<'a> {
    pattern: &'a OpPattern,
    current: usize,
}

impl<'a> Iterator for EnumMatchesOpPattern<'a> {
    type Item = Pauli;

    fn next(&mut self) -> Option<Self::Item> {
        use OpPattern::*;
        match self.pattern {
            Identity if self.current == 0 => {
                self.current += 1;
                Some(Pauli::I)
            }
            Single(op) if self.current == 0 => {
                self.current += 1;
                Some((*op).into())
            }
            Double(left, right) if self.current < 2 => {
                let result = if self.current == 0 {
                    self.current += 1;
                    (*left).into()
                } else {
                    self.current += 1;
                    (*right).into()
                };
                Some(result)
            }
            XYZ if self.current < 3 => {
                self.current += 1;
                Some(unsafe { std::mem::transmute(self.current as u8) })
            }
            SingleOrIdentity(op) if self.current < 2 => {
                let result = if self.current == 0 {
                    self.current += 1;
                    (*op).into()
                } else {
                    self.current += 1;
                    Pauli::I
                };
                Some(result)
            }
            DoubleOrIdentity(a, b) if self.current < 3 => {
                let result = if self.current == 0 {
                    self.current += 1;
                    (*a).into()
                } else if self.current == 1 {
                    self.current += 1;
                    (*b).into()
                } else {
                    self.current += 1;
                    Pauli::I
                };
                Some(result)
            }
            XYZI => {
                let result = unsafe { std::mem::transmute(self.current as u8) };
                self.current += 1;
                Some(result)
            }
            _ => None,
        }
    }
}

impl PauliPattern {
    pub fn enumerate_matches<A: PauliStorage>(
        &self,
        n_qubits: usize,
    ) -> EnumMatchesPauliPattern<'_, A> {
        let mut start: usize = 0;
        let mut iters = Vec::new();
        let mut patterns = self.0.iter().peekable();
        while let Some(pat) = patterns.next() {
            match pat {
                Decorated::Position(op, pos) => {
                    for _ in start..*pos {
                        iters.push(OpPattern::Identity.enumerate_matches());
                    }
                    iters.push(op.enumerate_matches());
                    start = *pos + 1;
                }
                Decorated::Repeat(op, count) => {
                    for _ in 0..*count {
                        iters.push(op.enumerate_matches());
                    }
                }
                Decorated::Star(_) => panic!("Star patterns are not supported"),
            }
        }

        // match iters.len() to n_qubits
        for _ in iters.len()..n_qubits {
            iters.push(OpPattern::Identity.enumerate_matches());
        }

        EnumMatchesPauliPattern {
            n_qubits,
            iter: iters.into_iter().multi_cartesian_product(),
            _phantom: std::marker::PhantomData,
        }
    }
}

#[derive(Debug, Clone)]
pub struct EnumMatchesPauliPattern<'a, A: PauliStorage> {
    n_qubits: usize,
    iter: MultiProduct<EnumMatchesOpPattern<'a>>,
    _phantom: std::marker::PhantomData<A>,
}

impl<'a, A: PauliStorage> Iterator for EnumMatchesPauliPattern<'a, A> {
    type Item = PauliWord<A>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(item) = self.iter.next() {
            let mut result = PauliWord::new(self.n_qubits);
            for (i, bit) in item.iter().enumerate() {
                result.set(i, *bit);
            }
            return Some(result);
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pattern_enumeration() {
        let pat: PauliPattern = "[XY]1Z3".into();
        let items: Vec<_> = pat.enumerate_matches::<u64>(4).collect();
        assert!(items.contains(&"IXIZ".into()));
        assert!(items.contains(&"IYIZ".into()));

        let pat: PauliPattern = "Z?{4}".into();
        let items: Vec<_> = pat.enumerate_matches::<u64>(4).collect();
        assert!(items.contains(&"ZZZZ".into()));
        assert!(items.contains(&"ZZZI".into()));
        assert!(items.contains(&"ZZIZ".into()));
        assert!(items.contains(&"ZZII".into()));
        assert!(items.contains(&"ZIZZ".into()));
        assert!(items.contains(&"ZIZI".into()));
        assert!(items.contains(&"ZIIZ".into()));
        assert!(items.contains(&"ZIII".into()));
        assert!(items.contains(&"IZZZ".into()));
        assert!(items.contains(&"IZZI".into()));
        assert!(items.contains(&"IZIZ".into()));
        assert!(items.contains(&"IZII".into()));
        assert!(items.contains(&"IIZZ".into()));
        assert!(items.contains(&"IIZI".into()));
        assert!(items.contains(&"IIIZ".into()));
        assert!(items.contains(&"IIII".into()));
    }
}
