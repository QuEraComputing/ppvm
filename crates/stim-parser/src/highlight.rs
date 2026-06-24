// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Syntax highlighter for Stim source text — a public, renderer-agnostic
//! companion to the canonical printer.
//!
//! [`highlight_html`] and [`highlight_ansi`] take any Stim string (typically
//! the canonical output of [`StimPrint::to_stim`](crate::print::StimPrint),
//! but raw input works too) and wrap its tokens for display. The lexer is
//! line-aware — the first identifier on a line is the instruction mnemonic
//! (or `REPEAT`) — so it never has to be a full parser.

/// Token category, used to pick a colour. The first identifier on a line is
/// the instruction mnemonic (or `REPEAT`) and is highlighted as a keyword.
#[derive(Clone, Copy)]
enum Tok {
    Keyword,
    Ident,
    Number,
    Punct,
    /// Whitespace / newlines — emitted verbatim, never coloured.
    Plain,
}

fn lex(src: &str) -> Vec<(Tok, &str)> {
    let b = src.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    // The first identifier of each line is the instruction name / `REPEAT`.
    let mut line_has_keyword = false;
    while i < b.len() {
        let c = b[i];
        if c == b'\n' {
            out.push((Tok::Plain, &src[i..i + 1]));
            i += 1;
            line_has_keyword = false;
            continue;
        }
        if matches!(c, b' ' | b'\t' | b'\r') {
            let start = i;
            while i < b.len() && matches!(b[i], b' ' | b'\t' | b'\r') {
                i += 1;
            }
            out.push((Tok::Plain, &src[start..i]));
            continue;
        }
        if c.is_ascii_alphabetic() || c == b'_' {
            let start = i;
            while i < b.len() && (b[i].is_ascii_alphanumeric() || b[i] == b'_') {
                i += 1;
            }
            let tok = if line_has_keyword {
                Tok::Ident
            } else {
                Tok::Keyword
            };
            line_has_keyword = true;
            out.push((tok, &src[start..i]));
            continue;
        }
        // Number: digits, optionally signed / decimal / exponent. A leading
        // `-`/`+`/`.` only starts a number when followed by a digit (so the
        // `-` in `rec[-1]` is part of the number, but a stray `*` is punct).
        let num_start = c.is_ascii_digit()
            || (matches!(c, b'-' | b'+' | b'.') && i + 1 < b.len() && b[i + 1].is_ascii_digit());
        if num_start {
            let start = i;
            if matches!(b[i], b'-' | b'+') {
                i += 1;
            }
            while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'.') {
                i += 1;
            }
            if i < b.len() && matches!(b[i], b'e' | b'E') {
                i += 1;
                if i < b.len() && matches!(b[i], b'-' | b'+') {
                    i += 1;
                }
                while i < b.len() && b[i].is_ascii_digit() {
                    i += 1;
                }
            }
            out.push((Tok::Number, &src[start..i]));
            continue;
        }
        // Anything else (brackets, parens, braces, `*`, `,`, `=`) is punctuation.
        out.push((Tok::Punct, &src[i..i + 1]));
        i += 1;
    }
    out
}

// Colours chosen to read on a light background (Jupyter's default) while
// staying legible on dark themes. Tweak here to restyle.
fn html_color(tok: Tok) -> Option<&'static str> {
    match tok {
        Tok::Keyword => Some("#8250df"), // purple — instruction / REPEAT
        Tok::Ident => Some("#0550ae"),   // blue — tag names, pauli targets
        Tok::Number => Some("#953800"),  // orange — args, indices, counts
        Tok::Punct => Some("#6e7781"),   // grey — [] () {} * , =
        Tok::Plain => None,
    }
}

fn ansi_code(tok: Tok) -> Option<&'static str> {
    match tok {
        Tok::Keyword => Some("\x1b[35m"), // magenta
        Tok::Ident => Some("\x1b[34m"),   // blue
        Tok::Number => Some("\x1b[33m"),  // yellow
        Tok::Punct => Some("\x1b[90m"),   // bright black / grey
        Tok::Plain => None,
    }
}

fn html_escape(s: &str, out: &mut String) {
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            _ => out.push(ch),
        }
    }
}

/// Render Stim source as a syntax-highlighted HTML `<pre>` block — suitable
/// for a Jupyter `_repr_html_`.
pub fn highlight_html(src: &str) -> String {
    let mut s = String::from(
        "<pre style=\"line-height:1.4; font-family:'JetBrains Mono',ui-monospace,monospace;\">",
    );
    for (tok, text) in lex(src) {
        match html_color(tok) {
            Some(color) => {
                s.push_str("<span style=\"color:");
                s.push_str(color);
                s.push_str("\">");
                html_escape(text, &mut s);
                s.push_str("</span>");
            }
            None => html_escape(text, &mut s),
        }
    }
    s.push_str("</pre>");
    s
}

/// Render Stim source with ANSI colour escapes — suitable for the IPython
/// terminal's `_repr_pretty_`.
pub fn highlight_ansi(src: &str) -> String {
    let mut s = String::new();
    for (tok, text) in lex(src) {
        match ansi_code(tok) {
            Some(code) => {
                s.push_str(code);
                s.push_str(text);
                s.push_str("\x1b[0m");
            }
            None => s.push_str(text),
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn html_colors_keyword_number_and_escapes() {
        let html = highlight_html("H 0\nREPEAT 2 {\n    X 0\n}\n");
        assert!(html.starts_with("<pre"));
        assert!(html.contains("#8250df")); // a keyword was coloured (H / REPEAT / X)
        assert!(html.contains("#953800")); // a number was coloured (0 / 2)
        assert!(html.ends_with("</pre>"));
    }

    #[test]
    fn ansi_wraps_tokens_and_resets() {
        let ansi = highlight_ansi("H 0\n");
        assert!(ansi.contains("\x1b[35m")); // keyword colour
        assert!(ansi.contains("\x1b[0m")); // reset
        assert!(ansi.contains('H'));
    }

    #[test]
    fn rec_minus_index_is_one_number_token() {
        // `-1` inside rec[-1] should be a single number token, not punct+number.
        let toks = lex("CX rec[-1] 0");
        let nums: Vec<&str> = toks
            .iter()
            .filter(|(t, _)| matches!(t, Tok::Number))
            .map(|(_, s)| *s)
            .collect();
        assert!(nums.contains(&"-1"), "got {nums:?}");
    }
}
