// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Robustness fuzz: arbitrary inputs must never panic the parser.
//!
//! Both `parse` and `parse_extended` are public entry points that consume
//! untrusted text (Stim files, network blobs, REPL prompts). They must
//! convert every input into either `Ok(_)` or `Err(_)` — never panic,
//! unwind, infinite-loop, or stack-overflow.
//!
//! Proptest is used as a lightweight in-process fuzzer: ~256 cases per
//! property is enough to keep the suite fast (<1s) while still
//! exercising the grammar against random and Stim-flavoured inputs.

use proptest::prelude::*;
use stim_parser::prelude::parse;
use stim_parser::prelude::parse_extended;

/// Vocabulary the grammar recognises. Biases the mutator toward
/// almost-valid inputs that exercise more of the parse tree than
/// purely random ASCII would.
fn stim_token() -> impl Strategy<Value = &'static str> {
    prop_oneof![
        Just("H"),
        Just("CX"),
        Just("CNOT"),
        Just("M"),
        Just("MZ"),
        Just("MR"),
        Just("MPAD"),
        Just("X"),
        Just("Y"),
        Just("Z"),
        Just("S"),
        Just("S_DAG"),
        Just("I"),
        Just("DEPOLARIZE1"),
        Just("DEPOLARIZE2"),
        Just("PAULI_CHANNEL_1"),
        Just("PAULI_CHANNEL_2"),
        Just("I_ERROR"),
        Just("HERALDED_ERASE"),
        Just("DETECTOR"),
        Just("OBSERVABLE_INCLUDE"),
        Just("TICK"),
        Just("REPEAT"),
        Just("rec[-1]"),
        Just("[T]"),
        Just("[loss]"),
        Just("[correlated_loss]"),
        Just("[R_X(theta=0.5)]"),
        Just("[U3(theta=0.5,phi=1,lambda=0.25)]"),
        Just("(0.1)"),
        Just("(0.1, 0.2, 0.3)"),
        Just("(pi)"),
        Just("(0.5*pi)"),
        Just("{"),
        Just("}"),
        Just("0"),
        Just("1"),
        Just("17"),
        Just(" "),
        Just("\n"),
        Just("\t"),
        Just("# comment\n"),
    ]
}

fn stim_like() -> impl Strategy<Value = String> {
    prop::collection::vec(stim_token(), 0..40).prop_map(|toks| toks.join(""))
}

proptest! {
    /// Pure-random bytes (anything UTF-8-decodable).
    #[test]
    fn parse_never_panics_on_random_bytes(s in "\\PC{0,512}") {
        let _ = parse(&s);
        let _ = parse_extended(&s);
    }

    /// Random ASCII (the surface the grammar mostly lives on).
    #[test]
    fn parse_never_panics_on_random_ascii(bytes in prop::collection::vec(any::<u8>(), 0..512)) {
        if let Ok(s) = std::str::from_utf8(&bytes) {
            let _ = parse(s);
            let _ = parse_extended(s);
        }
    }

    /// Stim-flavoured token soup: 0–40 random vocabulary tokens
    /// concatenated. Exercises the parser deeper into the tree than
    /// purely random bytes (which mostly fail at `ident()`).
    #[test]
    fn parse_never_panics_on_stim_token_soup(s in stim_like()) {
        let _ = parse(&s);
        let _ = parse_extended(&s);
    }

    /// Pathologically long REPEAT counts must not panic / stack overflow.
    #[test]
    fn repeat_count_never_panics(n in any::<u64>()) {
        let src = format!("REPEAT {n} {{\nH 0\n}}\n");
        let _ = parse(&src);
        let _ = parse_extended(&src);
    }

    /// Many nested REPEATs. The grammar's chumsky `recursive(...)`
    /// descends recursively through REPEAT bodies; `parse` /
    /// `parse_extended` run on a 16 MiB dedicated stack thread, which is
    /// large enough to handle a few hundred levels comfortably in debug
    /// builds. 128 is a generous upper bound for any realistic Stim
    /// program.
    #[test]
    fn nested_repeats_never_panic(depth in 0usize..=128) {
        let mut src = String::new();
        for _ in 0..depth {
            src.push_str("REPEAT 1 {\n");
        }
        src.push_str("H 0\n");
        for _ in 0..depth {
            src.push_str("}\n");
        }
        let _ = parse(&src);
        let _ = parse_extended(&src);
    }
}
