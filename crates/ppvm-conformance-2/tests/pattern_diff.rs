// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::panic::{AssertUnwindSafe, catch_unwind};

use ppvm_conformance_2::tableau::{NewNarrow, OldNarrow};
use ppvm_lossy_pauli_word_2::LossyPauliWord as NewLossyWord;
use ppvm_pauli_sum_2::{PatternParseError, PauliPattern as NewPattern};
use ppvm_pauli_word::loss::LossyPauliWord as OldLossyWord;
use ppvm_pauli_word::pattern::{Contains, PauliPattern as OldPattern};
use ppvm_pauli_word::word::PauliWord as OldWord;
use ppvm_pauli_word_2::PauliWord as NewWord;
use ppvm_tableau::data::GeneralizedTableau as OldTableau;
use ppvm_tableau_2::GeneralizedTableau as NewTableau;

#[test]
fn parser_accepts_and_canonicalizes_the_old_surface() {
    for source in [
        "",
        "X1Y2Z3",
        "[XY]1Y2_3",
        "X?1Y2Z3",
        "[XY]*Z5",
        "[XY]{3}Z5",
        "[XYZ]?*",
        "[ZY]7",
        "[XX]0",
        "X{2",
    ] {
        let old = OldPattern::parse(source).unwrap_or_else(|error| {
            panic!("old rejected {source:?}: {error}");
        });
        let new = NewPattern::parse(source).unwrap_or_else(|error| {
            panic!("new rejected {source:?}: {error}");
        });
        assert_eq!(old.to_string(), new.to_string(), "source {source:?}");
    }
}

#[test]
fn parser_returns_typed_errors() {
    assert_eq!(
        NewPattern::parse("I0"),
        Err(PatternParseError::ExpectedAtom)
    );
    assert_eq!(
        NewPattern::parse("[]0"),
        Err(PatternParseError::EmptyAlternation)
    );
    assert_eq!(
        NewPattern::parse("X"),
        Err(PatternParseError::ExpectedDecoration)
    );
    assert!(OldPattern::parse("I0").is_err());
    assert!(OldPattern::parse("[]0").is_err());
    assert!(OldPattern::parse("X").is_err());
}

#[test]
fn old_contains_examples_match_for_ordinary_words() {
    let cases = [
        ("X0Y1Z2", "XYZ", true),
        ("X?0Y1Z2", "XYZ", true),
        ("X?0Y1Z2", "IYZ", true),
        ("[XY]0Y1Z2", "XYZ", true),
        ("[XY]0Y1Z2", "YYZ", true),
        ("[XY]?0Y1Z2", "XYZ", true),
        ("[XY]?0Y1Z2", "YYZ", true),
        ("[XY]?0Y1Z2", "IYZ", true),
        ("[XY]?*", "XYX", true),
        ("[XY]?*", "YYX", true),
        ("[XY]?*", "IYX", true),
        ("[XY]?{2}Z2", "XYZ", true),
        ("[XY]?{2}Z2", "YYZ", true),
        ("[XY]?{2}Z2", "IYZ", true),
        ("Z?*", "XYY", false),
    ];
    for (pattern, word, expected) in cases {
        assert_match_parity(pattern, word, expected);
    }
}

#[test]
fn sequential_stars_keep_the_old_greedy_no_backtracking_semantics() {
    for (word, expected) in [
        ("XXZIII", true),
        ("XYZIII", true),
        ("ZIIIII", false),
        ("XXYIII", false),
    ] {
        assert_match_parity("[XY]*Z2", word, expected);
    }
    for (word, expected) in [("XXZZ", true), ("XXYY", true), ("XYZY", false)] {
        assert_match_parity("X*Z*Y*", word, expected);
    }
}

#[test]
fn zero_count_repetition_keeps_old_edge_behavior() {
    for (word, expected) in [
        ("", true),
        ("I", true),
        ("Z", false),
        ("XI", true),
        ("X", false),
    ] {
        assert_match_parity("X{0}", word, expected);
    }
}

#[test]
fn lossy_words_match_exactly_as_old() {
    for pattern in ["Z?*", "X0", "_0", "[XY]*Z3"] {
        for word in ["L", "XL", "ILZI", "XXLZ", "IIZI"] {
            let old_pattern = OldPattern::parse(pattern).unwrap();
            let new_pattern = NewPattern::parse(pattern).unwrap();
            let old_word: OldLossyWord<[u8; 1]> = word.into();
            let new_word: NewLossyWord<[u8; 1]> = word.into();
            assert_eq!(
                old_pattern.contains(&old_word),
                new_pattern.matches(&new_word),
                "{pattern:?} on {word:?}"
            );
        }
    }
}

#[test]
fn bounded_enumeration_matches_old_in_order() {
    for (pattern, width) in [
        ("[XY]1Z3", 4),
        ("Z?{4}", 4),
        ("[XZ]?0Y2", 4),
        ("[XYZ]{2}", 3),
    ] {
        let old: Vec<_> = OldPattern::parse(pattern)
            .unwrap()
            .enumerate_matches::<u64>(width)
            .map(|word| word.to_string())
            .collect();
        let new: Vec<_> = NewPattern::parse(pattern)
            .unwrap()
            .enumerate_matches::<u64>(width)
            .map(|word| word.to_string())
            .collect();
        assert_eq!(old, new, "{pattern:?} at width {width}");
    }
}

#[test]
fn tableau_enumeration_rejects_stars_like_old() {
    let old: OldNarrow = OldTableau::new(3, 1e-12);
    let new: NewNarrow = NewTableau::new(3, 1e-12);
    let old_pattern = OldPattern::parse("Z?*").unwrap();
    let new_pattern = NewPattern::parse("Z?*").unwrap();

    assert!(catch_unwind(AssertUnwindSafe(|| old.trace(&old_pattern))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| new.trace(&new_pattern))).is_err());
}

fn assert_match_parity(pattern: &str, word: &str, expected: bool) {
    let old_pattern = OldPattern::parse(pattern).unwrap();
    let new_pattern = NewPattern::parse(pattern).unwrap();
    let old_word: OldWord<[u8; 1]> = word.into();
    let new_word: NewWord<[u8; 1]> = word.into();
    let old = old_pattern.contains(&old_word);
    let new = new_pattern.matches(&new_word);
    assert_eq!(old, expected, "old: {pattern:?} on {word:?}");
    assert_eq!(new, old, "new: {pattern:?} on {word:?}");
}
