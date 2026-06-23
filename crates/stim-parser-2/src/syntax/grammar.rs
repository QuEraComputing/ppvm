// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Chumsky 0.12 grammar for Stim source.
//!
//! Reads top-to-bottom: whitespace/comments -> numbers -> pi-expressions ->
//! identifiers -> tags -> args -> targets -> instruction line -> REPEAT block ->
//! program. Pure syntax; no table lookups.

use chumsky::error::Rich;
use chumsky::extra;
use chumsky::prelude::*;

type Extra<'src> = extra::Err<Rich<'src, char>>;

/// `# ...` comment, stopping before `\n` if present.
fn line_comment<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    just('#')
        .ignore_then(any().filter(|c: &char| *c != '\n').repeated())
        .ignored()
}

/// Pad: zero or more whitespace characters or `# ...` comments.
/// Includes newlines. Used between instructions.
pub(crate) fn pad<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    let ws = any().filter(|c: &char| c.is_whitespace()).ignored();
    choice((line_comment(), ws)).repeated().ignored()
}

/// Inline pad: spaces/tabs/CR only. Excludes comments and `\n`.
pub(crate) fn inline_pad<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    let ws = any()
        .filter(|c: &char| matches!(*c, ' ' | '\t' | '\r'))
        .ignored();
    ws.repeated().ignored()
}

/// At least one inline whitespace character. Used before each target.
pub(crate) fn inline_ws1<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    any()
        .filter(|c: &char| matches!(*c, ' ' | '\t' | '\r'))
        .repeated()
        .at_least(1)
        .ignored()
}

/// Optional trailing spaces/tabs plus an optional line comment.
pub(crate) fn trailing_pad<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    inline_pad().then(line_comment().or_not()).ignored()
}

/// Identifier: `[A-Za-z_][A-Za-z0-9_]*`. Returns owned `String`.
pub(crate) fn ident<'src>() -> impl Parser<'src, &'src str, String, Extra<'src>> + Clone {
    any()
        .filter(|c: &char| c.is_ascii_alphabetic() || *c == '_')
        .then(
            any()
                .filter(|c: &char| c.is_ascii_alphanumeric() || *c == '_')
                .repeated(),
        )
        .to_slice()
        .map(|s: &str| s.to_string())
}

/// Signed float: common decimal/scientific forms accepted by `f64`
/// parsing, excluding non-finite names such as `NaN`/`inf`.
pub(crate) fn signed_float<'src>() -> impl Parser<'src, &'src str, f64, Extra<'src>> + Clone {
    let digits = any()
        .filter(|c: &char| c.is_ascii_digit())
        .repeated()
        .at_least(1);
    let decimal = choice((
        digits.then(just('.')).then(digits.or_not()).ignored(),
        just('.').then(digits).ignored(),
    ));
    let exp_part = one_of("eE").then(one_of("+-").or_not()).then(digits);
    one_of("+-")
        .or_not()
        .then(choice((decimal, digits.ignored())))
        .then(exp_part.or_not())
        .to_slice()
        .map(|s: &str| s.parse::<f64>().expect("validated by combinator shape"))
}

/// Pi-expression: `pi`, `<num>*pi`, or plain number. Evaluates to f64.
pub(crate) fn pi_expr<'src>() -> impl Parser<'src, &'src str, f64, Extra<'src>> + Clone {
    let pi_kw = just("pi").to(std::f64::consts::PI);
    let num_then_pi = signed_float()
        .then(inline_pad().ignore_then(just("*pi")).or_not())
        .map(|(n, suffix)| {
            if suffix.is_some() {
                n * std::f64::consts::PI
            } else {
                n
            }
        });
    choice((pi_kw, num_then_pi))
}

use crate::ast::shared::{Tag, TagParam};

/// `<ident>=<pi_expr>` (Named) or `<pi_expr>` (Positional).
pub(crate) fn tag_param<'src>() -> impl Parser<'src, &'src str, TagParam, Extra<'src>> + Clone {
    let named = ident()
        .then_ignore(inline_pad())
        .then_ignore(just('='))
        .then_ignore(inline_pad())
        .then(pi_expr())
        .map(|(key, value)| TagParam::Named { key, value });
    let positional = pi_expr().map(TagParam::Positional);
    choice((named, positional))
}

/// Tag: `<ident>` or `<ident>(<tag_param>, ...)`.
pub(crate) fn tag<'src>() -> impl Parser<'src, &'src str, Tag, Extra<'src>> + Clone {
    let params = tag_param()
        .separated_by(inline_pad().then(just(',')).then(inline_pad()))
        .allow_trailing()
        .collect::<Vec<_>>()
        .delimited_by(just('(').then(inline_pad()), inline_pad().then(just(')')));
    ident().then(params.or_not()).map(|(name, params)| Tag {
        name,
        params: params.unwrap_or_default(),
    })
}

/// `[tag, tag, ...]`.
pub(crate) fn tags_block<'src>() -> impl Parser<'src, &'src str, Vec<Tag>, Extra<'src>> + Clone {
    tag()
        .separated_by(inline_pad().then(just(',')).then(inline_pad()))
        .allow_trailing()
        .collect::<Vec<_>>()
        .delimited_by(just('[').then(inline_pad()), inline_pad().then(just(']')))
}

/// `(pi_expr, pi_expr, ...)`.
pub(crate) fn args_block<'src>() -> impl Parser<'src, &'src str, Vec<f64>, Extra<'src>> + Clone {
    pi_expr()
        .separated_by(inline_pad().then(just(',')).then(inline_pad()))
        .allow_trailing()
        .collect::<Vec<_>>()
        .delimited_by(just('(').then(inline_pad()), inline_pad().then(just(')')))
}

use crate::syntax::raw::{RawSyntaxNode, RawTarget};
use chumsky::span::SimpleSpan;

/// One non-whitespace, non-`#`, non-`{`/`}` lexeme. Captures the span
/// so the validator can derive a line number for invalid-target errors.
pub(crate) fn target_lexeme<'src>() -> impl Parser<'src, &'src str, RawTarget, Extra<'src>> + Clone
{
    any()
        .filter(|c: &char| !c.is_whitespace() && *c != '#' && *c != '{' && *c != '}')
        .repeated()
        .at_least(1)
        .to_slice()
        .map_with(|s: &str, e| RawTarget {
            text: s.to_string(),
            span: e.span(),
        })
}

/// `<ident> [<tags>]? (<args>)?`. Returns name, tags, args, and the
/// span of the identifier (used for line-number reporting).
pub(crate) fn instruction_head<'src>()
-> impl Parser<'src, &'src str, (String, Vec<Tag>, Vec<f64>, SimpleSpan<usize>), Extra<'src>> + Clone
{
    ident()
        .map_with(|name, e| (name, e.span()))
        .then(tags_block().or_not())
        .then(args_block().or_not())
        .map(|(((name, span), tags), args)| {
            (
                name,
                tags.unwrap_or_default(),
                args.unwrap_or_default(),
                span,
            )
        })
}

/// End of an instruction: optional trailing spaces/comment followed by
/// newline, EOF, or a `}` that belongs to the enclosing REPEAT parser.
pub(crate) fn line_end<'src>() -> impl Parser<'src, &'src str, (), Extra<'src>> + Clone {
    trailing_pad()
        .then(choice((
            just('\n').ignored(),
            just('}').rewind().ignored(),
            end(),
        )))
        .ignored()
}

/// Instruction line: head + space-separated raw targets, terminated by
/// a newline / EOF / `}`. Newline is consumed; `}` is only peeked so
/// the enclosing REPEAT parser can consume it.
pub(crate) fn instruction_line<'src>()
-> impl Parser<'src, &'src str, RawSyntaxNode, Extra<'src>> + Clone {
    instruction_head()
        .then(
            inline_ws1()
                .ignore_then(target_lexeme())
                .repeated()
                .collect::<Vec<RawTarget>>(),
        )
        .map(
            |((name, tags, args, span), targets)| RawSyntaxNode::Instruction {
                name,
                tags,
                args,
                targets,
                span,
            },
        )
        .then_ignore(line_end())
}

/// `REPEAT <count> { <body> }`. `body` is the parser for a list of
/// nodes (instructions and nested REPEATs).
fn repeat_block<'src>(
    body: impl Parser<'src, &'src str, Vec<RawSyntaxNode>, Extra<'src>> + Clone + 'src,
) -> impl Parser<'src, &'src str, RawSyntaxNode, Extra<'src>> + Clone + 'src {
    let digits = any()
        .filter(|c: &char| c.is_ascii_digit())
        .repeated()
        .at_least(1)
        .to_slice()
        .try_map(|s: &str, span| {
            s.parse::<u64>()
                .map_err(|_| Rich::custom(span, format!("invalid REPEAT count {s:?}")))
        });
    just("REPEAT")
        .map_with(|_, e| e.span())
        .then_ignore(inline_ws1())
        .then(digits)
        .then_ignore(inline_pad())
        .then_ignore(just('{'))
        .then_ignore(pad())
        .then(body)
        .then_ignore(pad())
        .then_ignore(just('}'))
        .map(|((span, count), body)| RawSyntaxNode::Repeat { count, body, span })
}

/// Top-level program parser. Recursively defines the body shared by
/// the program and REPEAT blocks.
pub(crate) fn program_parser<'src>() -> impl Parser<'src, &'src str, Vec<RawSyntaxNode>, Extra<'src>>
{
    recursive(|body| {
        let item = choice((repeat_block(body.clone()), instruction_line()));
        pad()
            .ignore_then(item)
            .repeated()
            .collect::<Vec<RawSyntaxNode>>()
            .then_ignore(pad())
    })
    .then_ignore(end())
}

#[cfg(test)]
mod tests {
    use crate::ast::TagParam;
    use crate::syntax::grammar::*;
    use crate::syntax::raw::RawSyntaxNode;

    fn run<'src, T>(p: impl Parser<'src, &'src str, T, Extra<'src>>, src: &'src str) -> T {
        p.parse(src).into_result().expect("parse failed")
    }

    #[test]
    fn ident_matches_alpha_then_alphanumeric() {
        assert_eq!(run(ident(), "H"), "H");
        assert_eq!(run(ident(), "DEPOLARIZE1"), "DEPOLARIZE1");
        assert_eq!(run(ident(), "_x"), "_x");
        assert_eq!(run(ident(), "R_X"), "R_X");
    }

    #[test]
    fn signed_float_parses_common_shapes() {
        assert_eq!(run(signed_float(), "0"), 0.0);
        assert_eq!(run(signed_float(), "0.5"), 0.5);
        assert_eq!(run(signed_float(), "-0.5"), -0.5);
        assert_eq!(run(signed_float(), "+1.0e-3"), 1.0e-3);
        assert_eq!(run(signed_float(), "42"), 42.0);
        assert_eq!(run(signed_float(), "1."), 1.0);
        assert_eq!(run(signed_float(), ".5"), 0.5);
    }

    #[test]
    fn pi_expr_parses_pi_keyword_coeff_and_plain_number() {
        assert_eq!(run(pi_expr(), "pi"), std::f64::consts::PI);
        assert_eq!(run(pi_expr(), "0.5*pi"), 0.5 * std::f64::consts::PI);
        assert_eq!(run(pi_expr(), "0.5 *pi"), 0.5 * std::f64::consts::PI);
        assert_eq!(run(pi_expr(), "-2*pi"), -2.0 * std::f64::consts::PI);
        assert_eq!(run(pi_expr(), "0.5"), 0.5);
    }

    #[test]
    fn pad_eats_comments_and_newlines() {
        // pad's output is (); we just verify it doesn't error.
        let p = pad().then_ignore(end());
        assert!(p.parse("").into_result().is_ok());
        assert!(p.parse("   \n\t# comment\n").into_result().is_ok());
    }

    #[test]
    fn inline_pad_does_not_consume_comment_or_newline() {
        // inline_pad should leave both comment starts and newlines to the
        // line-ending parser.
        let p = inline_pad().then_ignore(just('#'));
        assert!(p.parse(" \t#").into_result().is_ok());
        let p = inline_pad().then_ignore(just('\n'));
        assert!(p.parse(" \t\n").into_result().is_ok());
    }

    #[test]
    fn trailing_pad_consumes_trailing_comment() {
        let p = trailing_pad()
            .then_ignore(just('\n').or_not())
            .then_ignore(end());
        assert!(p.parse(" # comment\n").into_result().is_ok());
        assert!(p.parse("   ").into_result().is_ok());
    }

    #[test]
    fn tag_with_no_params() {
        let t = run(tag(), "T");
        assert_eq!(t.name, "T");
        assert!(t.params.is_empty());
    }

    #[test]
    fn tag_with_positional_params() {
        let t = run(tag(), "R(0.5, 1.0)");
        assert_eq!(t.name, "R");
        assert_eq!(t.params.len(), 2);
        assert!(matches!(&t.params[0], TagParam::Positional(v) if (v - 0.5).abs() < 1e-12));
    }

    #[test]
    fn tag_with_named_param() {
        let t = run(tag(), "R_X(theta=0.5*pi)");
        assert_eq!(t.name, "R_X");
        assert_eq!(t.params.len(), 1);
        match &t.params[0] {
            TagParam::Named { key, value } => {
                assert_eq!(key, "theta");
                assert!((value - 0.5 * std::f64::consts::PI).abs() < 1e-12);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn tags_block_parses_multiple_tags() {
        let ts = run(tags_block(), "[T, R(0.5)]");
        assert_eq!(ts.len(), 2);
        assert_eq!(ts[0].name, "T");
        assert_eq!(ts[1].name, "R");
    }

    #[test]
    fn args_block_parses_csv_floats() {
        let a = run(args_block(), "(0.1, 0.2, 0.3)");
        assert_eq!(a, vec![0.1, 0.2, 0.3]);
    }

    #[test]
    fn args_block_with_pi_exprs() {
        let a = run(args_block(), "(pi, 0.5*pi)");
        assert!((a[0] - std::f64::consts::PI).abs() < 1e-12);
        assert!((a[1] - 0.5 * std::f64::consts::PI).abs() < 1e-12);
    }

    #[test]
    fn target_lexeme_reads_a_non_whitespace_run() {
        let t = run(target_lexeme(), "0");
        assert_eq!(t.text, "0");
        let t = run(target_lexeme(), "rec[-1]");
        assert_eq!(t.text, "rec[-1]");
    }

    #[test]
    fn target_lexeme_stops_at_brace() {
        // target_lexeme should not consume `}`.
        let p = target_lexeme().then_ignore(just('}'));
        let t = p.parse("0}").into_result().expect("parse failed");
        assert_eq!(t.text, "0");
    }

    #[test]
    fn instruction_head_with_tags_and_args() {
        let (name, tags, args, _span) = run(instruction_head(), "S[T](0.5)");
        assert_eq!(name, "S");
        assert_eq!(tags.len(), 1);
        assert_eq!(tags[0].name, "T");
        assert_eq!(args, vec![0.5]);
    }

    #[test]
    fn instruction_head_no_tags_no_args() {
        let (name, tags, args, _span) = run(instruction_head(), "H");
        assert_eq!(name, "H");
        assert!(tags.is_empty());
        assert!(args.is_empty());
    }

    #[test]
    fn instruction_line_with_targets() {
        let n = run(instruction_line(), "CX 0 1 2 3");
        match n {
            RawSyntaxNode::Instruction { name, targets, .. } => {
                assert_eq!(name, "CX");
                let texts: Vec<_> = targets.iter().map(|t| t.text.clone()).collect();
                assert_eq!(texts, vec!["0", "1", "2", "3"]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn instruction_line_no_targets() {
        let n = run(instruction_line(), "TICK");
        match n {
            RawSyntaxNode::Instruction { name, targets, .. } => {
                assert_eq!(name, "TICK");
                assert!(targets.is_empty());
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn program_parses_multiline() {
        let p = program_parser();
        let nodes = p.parse("X 0\nY 1\n").into_result().expect("parse failed");
        assert_eq!(nodes.len(), 2);
    }

    #[test]
    fn program_handles_blank_lines_and_comments() {
        let p = program_parser();
        let nodes = p
            .parse("\n# header\nX 0\n# mid\nY 1\n")
            .into_result()
            .expect("parse failed");
        assert_eq!(nodes.len(), 2);
    }

    #[test]
    fn program_handles_trailing_comments() {
        let p = program_parser();
        let nodes = p
            .parse("X 0 # flip\nY 1 # measure\n")
            .into_result()
            .expect("parse failed");
        assert_eq!(nodes.len(), 2);
    }

    #[test]
    fn program_parses_repeat_block() {
        let p = program_parser();
        let nodes = p
            .parse("REPEAT 3 {\n    X 0\n    M 0\n}\n")
            .into_result()
            .expect("parse failed");
        assert_eq!(nodes.len(), 1);
        match &nodes[0] {
            RawSyntaxNode::Repeat { count, body, .. } => {
                assert_eq!(*count, 3);
                assert_eq!(body.len(), 2);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn program_parses_nested_repeat() {
        let p = program_parser();
        let nodes = p
            .parse("REPEAT 2 {\n  REPEAT 3 {\n    H 0\n  }\n}\n")
            .into_result()
            .expect("parse failed");
        match &nodes[0] {
            RawSyntaxNode::Repeat { body, .. } => match &body[0] {
                RawSyntaxNode::Repeat { count, .. } => assert_eq!(*count, 3),
                other => panic!("{other:?}"),
            },
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn program_parses_one_line_repeat() {
        let p = program_parser();
        let nodes = p
            .parse("REPEAT 5 { H 0 }")
            .into_result()
            .expect("parse failed");
        match &nodes[0] {
            RawSyntaxNode::Repeat { count, body, .. } => {
                assert_eq!(*count, 5);
                assert_eq!(body.len(), 1);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn program_rejects_oversized_repeat_count_without_panicking() {
        let result = std::panic::catch_unwind(|| {
            program_parser()
                .parse("REPEAT 184467440737095516160 { H 0 }")
                .into_result()
        });
        assert!(result.is_ok(), "parser panicked");
        assert!(result.unwrap().is_err());
    }

    #[test]
    fn program_rejects_malformed_repeat_header() {
        let p = program_parser();
        assert!(p.parse("REPEAT3 { H 0 }").into_result().is_err());
    }

    #[test]
    fn program_rejects_unclosed_repeat_block() {
        let p = program_parser();
        assert!(p.parse("REPEAT 2 {\nH 0\n").into_result().is_err());
    }

    #[test]
    fn program_rejects_unmatched_close_brace() {
        let p = program_parser();
        assert!(p.parse("H 0\n}").into_result().is_err());
    }
}
