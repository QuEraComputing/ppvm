use super::data::{Decorated, NotIdentity, OpPattern, PauliPattern};
use crate::char::Pauli;
use crate::traits::PauliStorage;
use crate::word::PauliWord;
use std::hash::BuildHasher;
use std::iter::Peekable;

pub trait Contains<T> {
    fn contains(&self, item: &T) -> bool;
}

impl Contains<Pauli> for NotIdentity {
    fn contains(&self, item: &Pauli) -> bool {
        (*self as u8) == (*item as u8)
    }
}

impl Contains<Pauli> for OpPattern {
    fn contains(&self, item: &Pauli) -> bool {
        match self {
            OpPattern::Identity => *item == Pauli::I,
            OpPattern::Single(op) => op.contains(item),
            OpPattern::Double(left, right) => left.contains(item) || right.contains(item),
            OpPattern::XYZ => *item == Pauli::X || *item == Pauli::Y || *item == Pauli::Z,
            OpPattern::SingleOrIdentity(op) => op.contains(item) || *item == Pauli::I,
            OpPattern::DoubleOrIdentity(left, right) => {
                left.contains(item) || right.contains(item) || *item == Pauli::I
            }
            OpPattern::XYZI => *item == Pauli::I,
        }
    }
}

enum Step {
    NextPattern,
    DontMatch,
}

fn match_position<I: Iterator<Item = (usize, Pauli)>>(
    chars: &mut Peekable<I>,
    pattern: &OpPattern,
    pos: &usize,
) -> Step {
    while let Some((ch_pos, ch)) = chars.next() {
        if ch_pos == *pos {
            if pattern.contains(&ch) {
                return Step::NextPattern;
            } else {
                return Step::DontMatch;
            }
        } else if ch == Pauli::I {
            continue;
        } else {
            // there are non-identity
            // before matching this pattern
            return Step::DontMatch;
        }
    }
    // all identity and length < pos
    Step::DontMatch
}

fn match_star<I: Iterator<Item = (usize, Pauli)>>(
    chars: &mut Peekable<I>,
    pattern: &OpPattern,
) -> Step {
    while let Some((_, ch)) = chars.peek() {
        if pattern.contains(ch) {
            chars.next();
            continue;
        } else {
            // doesn't match but ok
            // match next pattern
            return Step::NextPattern;
        }
    }
    // exhausted all Pauli all matches
    Step::NextPattern
}

fn match_repeat<I: Iterator<Item = (usize, Pauli)>>(
    chars: &mut Peekable<I>,
    pattern: &OpPattern,
    count: usize,
) -> Step {
    let mut matched = 0;
    while let Some((_, ch)) = chars.peek() {
        if pattern.contains(ch) {
            chars.next();
            matched += 1;
            if matched == count {
                return Step::NextPattern;
            }
        } else {
            // doesn't match
            if matched < count {
                return Step::DontMatch;
            } else {
                return Step::NextPattern;
            }
        }
    }
    // exhausted all Pauli all matches
    if matched == count {
        Step::NextPattern
    } else {
        Step::DontMatch
    }
}

impl<A: PauliStorage, H: Default + BuildHasher + Clone> Contains<PauliWord<A, H>> for PauliPattern {
    fn contains(&self, item: &PauliWord<A, H>) -> bool {
        let mut chars = item.iter().enumerate().peekable();
        let mut patterns = self.iter();
        while let Some(current) = patterns.next() {
            let step = match current {
                Decorated::Position(op, pos) => match_position(&mut chars, op, pos),
                Decorated::Star(op) => match_star(&mut chars, op),
                Decorated::Repeat(op, count) => match_repeat(&mut chars, op, *count),
            };

            match step {
                Step::NextPattern => continue,
                Step::DontMatch => return false,
            }
        }

        // check if there are remaining non-ident in `item`
        while let Some((_, ch)) = chars.next() {
            if ch == Pauli::I {
                continue;
            } else {
                // found a non-identity
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_match_not_ident() {
        assert!(NotIdentity::X.contains(&Pauli::X));
        assert!(NotIdentity::Y.contains(&Pauli::Y));
        assert!(NotIdentity::Z.contains(&Pauli::Z));
    }

    #[test]
    fn test_match() {
        let pat = PauliPattern::parse("X0Y1Z2").unwrap();
        let word: PauliWord<u64> = "XYZ".into();
        assert!(pat.contains(&word));

        let pat = PauliPattern::parse("X?0Y1Z2").unwrap();
        let word: PauliWord<u64> = "XYZ".into();
        assert!(pat.contains(&word));
        let word: PauliWord<u64> = "IYZ".into();
        assert!(pat.contains(&word));

        let pat = PauliPattern::parse("[XY]0Y1Z2").unwrap();
        let word: PauliWord<u64> = "XYZ".into();
        assert!(pat.contains(&word));
        let word: PauliWord<u64> = "YYZ".into();
        assert!(pat.contains(&word));

        let pat = PauliPattern::parse("[XY]?0Y1Z2").unwrap();
        let word: PauliWord<u64> = "XYZ".into();
        assert!(pat.contains(&word));
        let word: PauliWord<u64> = "YYZ".into();
        assert!(pat.contains(&word));
        let word: PauliWord<u64> = "IYZ".into();
        assert!(pat.contains(&word));

        let pat = PauliPattern::parse("[XY]?*").unwrap();
        let word: PauliWord<u64> = "XYX".into();
        assert!(pat.contains(&word));
        let word: PauliWord<u64> = "YYX".into();
        assert!(pat.contains(&word));
        let word: PauliWord<u64> = "IYX".into();
        assert!(pat.contains(&word));

        let pat = PauliPattern::parse("[XY]?{2}Z2").unwrap();
        let word: PauliWord<u64> = "XYZ".into();
        assert!(pat.contains(&word));
        let word: PauliWord<u64> = "YYZ".into();
        assert!(pat.contains(&word));
        let word: PauliWord<u64> = "IYZ".into();
        assert!(pat.contains(&word));
    }

    #[test]
    fn test_not_match() {
        let pat: PauliPattern = "Z?*".into();
        let word: PauliWord<u64> = "XYY".into();
        assert!(!pat.contains(&word));
    }
}
