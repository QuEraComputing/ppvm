// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Parse the `device circuit.observable` header value into Pauli-sum terms.
//!
//! Lives in `ppvm-vihaco` rather than `ppvm-runtime` because it only parses
//! header syntax that `ppvm-vihaco` consumes; `ppvm-runtime` has no notion of
//! a textual observable. The unrelated `PauliPattern::parse` in
//! `ppvm-runtime` is a *matcher* grammar (alternation, star, positional
//! anchors) and shares only the `I/X/Y/Z` alphabet with this one.
//!
//! Grammar (informal):
//! ```text
//! sum     = WS* term (WS* sign WS* term)* WS*
//! term    = coefficient (WS* '*')? WS* pauli_word
//!         | pauli_word
//! sign    = '+' | '-'
//! coefficient
//!         = digits ('.' digits?)? ([eE] [+-]? digits)?
//!         | '.' digits ([eE] [+-]? digits)?
//! pauli_word
//!         = [IXYZ]{n_qubits}
//! ```
//!
//! - The first term may have a leading `+` or `-` (no sign means `+`).
//! - Absent coefficient defaults to `1.0`.
//! - The `*` is only legal *after* a coefficient — bare `*ZZ` is rejected.
//! - The word must be exactly `n_qubits` characters from `I/X/Y/Z`.
//!
//! Rejected at parse time:
//! - Empty or whitespace-only input.
//! - Bare coefficients with no Pauli word.
//! - Words shorter or longer than `n_qubits`.
//! - Invalid Pauli characters.
//! - Missing `+`/`-` between terms.

use chumsky::error::Simple;
use chumsky::extra;
use chumsky::prelude::*;
use eyre::{Result, eyre};

type Err<'src> = extra::Err<Simple<'src, char>>;

/// Parse a Pauli-sum string like `"1.0*ZZ + 0.5*XX - 0.3*YY"` into a list of
/// `(word_source, coefficient)` pairs. Callers convert the word source to
/// their preferred Pauli-word type via `PauliWord::from` / `LossyPauliWord::from`.
pub fn parse_pauli_sum_terms(input: &str, n_qubits: usize) -> Result<Vec<(String, f64)>> {
    if n_qubits == 0 {
        return Err(eyre!("n_qubits must be at least 1"));
    }
    if input.trim().is_empty() {
        return Err(eyre!("empty Pauli-sum input"));
    }
    pauli_sum_parser(n_qubits)
        .parse(input)
        .into_result()
        .map_err(|errs| {
            let joined = errs
                .into_iter()
                .map(|e| e.to_string())
                .collect::<Vec<_>>()
                .join("; ");
            eyre!("invalid Pauli-sum `{input}`: {joined}")
        })
}

fn pauli_sum_parser<'src>(
    n_qubits: usize,
) -> impl Parser<'src, &'src str, Vec<(String, f64)>, Err<'src>> {
    // Numeric literal: digits[.digits?]?[exp]?  or  .digits[exp]?
    let digits = text::digits(10).at_least(1);
    let exp = one_of("eE")
        .then(one_of("+-").or_not())
        .then(digits)
        .ignored();
    let mantissa = choice((
        digits
            .then(just('.').then(digits.or_not()).or_not())
            .ignored(),
        just('.').then(digits).ignored(),
    ));
    let coefficient = mantissa.then(exp.or_not()).to_slice().map(|s: &str| {
        s.parse::<f64>()
            .expect("mantissa+exponent grammar validates parseability")
    });

    // Pauli word: exactly n_qubits chars from IXYZ.
    let pauli_word = one_of("IXYZ")
        .repeated()
        .exactly(n_qubits)
        .collect::<String>();

    // Term: coefficient [*] word | bare word.
    let term_with_coeff = coefficient
        .then_ignore(just('*').padded().or_not())
        .then(pauli_word)
        .map(|(c, w)| (w, c));
    let term_bare = pauli_word.map(|w| (w, 1.0));
    let term = choice((term_with_coeff, term_bare));

    // Sign factor: +1 / -1.
    let sign = choice((just('+').to(1.0_f64), just('-').to(-1.0_f64)));

    // First term: optional leading sign.
    let first = sign
        .padded()
        .or_not()
        .then(term)
        .map(|(s, (w, c))| (w, s.unwrap_or(1.0) * c));

    // Subsequent terms: required + or - before the term.
    let rest = sign.padded().then(term).map(|(s, (w, c))| (w, s * c));

    let inner = first
        .then(rest.repeated().collect::<Vec<_>>())
        .map(|(first, mut rest)| {
            rest.insert(0, first);
            rest
        });

    inner.padded().then_ignore(end())
}

#[cfg(test)]
mod tests {
    use super::*;

    // ─── Happy paths ────────────────────────────────────────────────────

    #[test]
    fn single_term_no_coefficient() {
        let got = parse_pauli_sum_terms("ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 1.0)]);
    }

    #[test]
    fn single_term_with_explicit_star() {
        let got = parse_pauli_sum_terms("0.5*ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 0.5)]);
    }

    #[test]
    fn single_term_without_star() {
        let got = parse_pauli_sum_terms("0.5ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 0.5)]);
    }

    #[test]
    fn integer_coefficient() {
        let got = parse_pauli_sum_terms("3XY", 2).unwrap();
        assert_eq!(got, vec![("XY".to_string(), 3.0)]);
    }

    #[test]
    fn leading_negative_no_coefficient() {
        let got = parse_pauli_sum_terms("-ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), -1.0)]);
    }

    #[test]
    fn leading_negative_with_coefficient() {
        let got = parse_pauli_sum_terms("-0.5*ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), -0.5)]);
    }

    #[test]
    fn leading_plus_no_coefficient() {
        let got = parse_pauli_sum_terms("+ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 1.0)]);
    }

    #[test]
    fn multiple_terms_with_addition() {
        let got = parse_pauli_sum_terms("1.0*ZZ + 0.5*XX", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 1.0), ("XX".to_string(), 0.5)]);
    }

    #[test]
    fn multiple_terms_with_subtraction() {
        let got = parse_pauli_sum_terms("ZZ - XX", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 1.0), ("XX".to_string(), -1.0)]);
    }

    #[test]
    fn three_term_mixed() {
        let got = parse_pauli_sum_terms("1.0*ZZ + 0.5*XX - 0.3*YY", 2).unwrap();
        assert_eq!(
            got,
            vec![
                ("ZZ".to_string(), 1.0),
                ("XX".to_string(), 0.5),
                ("YY".to_string(), -0.3),
            ]
        );
    }

    #[test]
    fn scientific_notation_coefficient() {
        let got = parse_pauli_sum_terms("1e-3*ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 1e-3)]);
    }

    #[test]
    fn unsigned_positive_exponent() {
        let got = parse_pauli_sum_terms("1e3*ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 1000.0)]);
    }

    #[test]
    fn uppercase_exponent_with_explicit_sign() {
        let got = parse_pauli_sum_terms("2.5E+2*ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 250.0)]);
    }

    #[test]
    fn coefficient_with_trailing_dot() {
        let got = parse_pauli_sum_terms("1.*ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 1.0)]);
    }

    #[test]
    fn whitespace_tolerance_aggressive() {
        let got = parse_pauli_sum_terms("  1.0  *  ZZ  +  0.5  *  XX  ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 1.0), ("XX".to_string(), 0.5)]);
    }

    #[test]
    fn coefficient_starts_with_dot() {
        let got = parse_pauli_sum_terms(".5*ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 0.5)]);
    }

    #[test]
    fn identity_in_word_is_valid() {
        let got = parse_pauli_sum_terms("IZIZ", 4).unwrap();
        assert_eq!(got, vec![("IZIZ".to_string(), 1.0)]);
    }

    #[test]
    fn single_qubit_n_equals_one() {
        let got = parse_pauli_sum_terms("0.5*X + Y - 0.25Z", 1).unwrap();
        assert_eq!(
            got,
            vec![
                ("X".to_string(), 0.5),
                ("Y".to_string(), 1.0),
                ("Z".to_string(), -0.25),
            ]
        );
    }

    #[test]
    fn duplicate_word_emits_two_terms() {
        // Coefficient collapse is the caller's job (via `PauliSum +=`); the
        // parser stays a thin syntactic layer.
        let got = parse_pauli_sum_terms("0.5*ZZ + 0.5*ZZ", 2).unwrap();
        assert_eq!(got, vec![("ZZ".to_string(), 0.5), ("ZZ".to_string(), 0.5)]);
    }

    // ─── Error paths ────────────────────────────────────────────────────

    #[test]
    fn rejects_empty_input() {
        let err = parse_pauli_sum_terms("", 2).unwrap_err();
        assert!(err.to_string().contains("empty"), "{err}");
    }

    #[test]
    fn rejects_whitespace_only() {
        let err = parse_pauli_sum_terms("   \t  \n  ", 2).unwrap_err();
        assert!(err.to_string().contains("empty"), "{err}");
    }

    #[test]
    fn rejects_zero_qubits() {
        let err = parse_pauli_sum_terms("ZZ", 0).unwrap_err();
        assert!(err.to_string().contains("n_qubits"), "{err}");
    }

    #[test]
    fn rejects_bare_coefficient() {
        let err = parse_pauli_sum_terms("0.5", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }

    #[test]
    fn rejects_bare_sign() {
        let err = parse_pauli_sum_terms("+", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }

    #[test]
    fn rejects_short_word() {
        let err = parse_pauli_sum_terms("Z", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }

    #[test]
    fn rejects_long_word() {
        let err = parse_pauli_sum_terms("ZZZ", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }

    #[test]
    fn rejects_invalid_pauli_character() {
        let err = parse_pauli_sum_terms("ZF", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }

    #[test]
    fn rejects_lowercase_pauli_character() {
        let err = parse_pauli_sum_terms("zz", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }

    #[test]
    fn rejects_missing_separator_between_terms() {
        let err = parse_pauli_sum_terms("ZZ XX", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }

    #[test]
    fn rejects_double_sign() {
        let err = parse_pauli_sum_terms("++ZZ", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }

    #[test]
    fn rejects_star_without_coefficient() {
        // Bare `*ZZ` is rejected — `*` is only legal after a coefficient.
        let err = parse_pauli_sum_terms("*ZZ", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }

    #[test]
    fn rejects_trailing_garbage() {
        let err = parse_pauli_sum_terms("ZZ garbage", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }

    #[test]
    fn rejects_trailing_garbage_after_complete_sum() {
        // Distinct from `rejects_missing_separator_between_terms`: here a
        // complete two-term sum parses successfully and the failure is the
        // `end()` check after the last term.
        let err = parse_pauli_sum_terms("ZZ + XX trailing", 2).unwrap_err();
        assert!(err.to_string().contains("invalid Pauli-sum"), "{err}");
    }
}
