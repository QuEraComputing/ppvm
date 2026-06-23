// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod sinks;
pub use sinks::{Collect, FailFast};

mod span;
pub use span::{LineMap, Span};

use std::fmt;
use std::sync::Arc;

/// Severity of a diagnostic. Only `Error` aborts a `FailFast` run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    Error,
    Warning,
}

/// A sink's continuation decision: keep processing the current stage, or
/// abort it as soon as possible. The handler returning `Flow` is how the
/// effect model "handles errors as continuations".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flow {
    Continue,
    Abort,
}

/// One diagnostic, carrying a span so every message can render `line:col`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    pub severity: Severity,
    pub span: Span,
    pub message: String,
    /// Stable short kind tag (e.g. "unknown-instruction") for matching
    /// without string-sniffing.
    pub code: Option<&'static str>,
}

impl Diagnostic {
    pub fn error(span: Span, code: &'static str, message: impl Into<String>) -> Self {
        Diagnostic {
            severity: Severity::Error,
            span,
            message: message.into(),
            code: Some(code),
        }
    }

    pub fn warning(span: Span, code: &'static str, message: impl Into<String>) -> Self {
        Diagnostic {
            severity: Severity::Warning,
            span,
            message: message.into(),
            code: Some(code),
        }
    }
}

/// A handler the pipeline emits diagnostics to. The returned `Flow` tells
/// the emitting stage whether to continue (recover) or abort.
pub trait DiagnosticSink {
    fn emit(&mut self, diagnostic: Diagnostic) -> Flow;
}

/// Marker returned by a pipeline transition that was told to `Abort`.
/// The diagnostics themselves live in the caller-owned sink.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Aborted;

/// Aggregate returned by the Tier-1 `parse`/`parse_extended` functions on
/// failure. Owns the source text, so it can render both terse `line:col`
/// summaries ([`Display`](fmt::Display)) and rich, source-pointing reports
/// ([`render`](Diagnostics::render)).
#[derive(Debug, Clone)]
pub struct Diagnostics {
    items: Vec<Diagnostic>,
    source: Arc<str>,
    line_map: Arc<LineMap>,
}

impl Diagnostics {
    pub fn new(items: Vec<Diagnostic>, source: Arc<str>) -> Self {
        let line_map = Arc::new(LineMap::new(&source));
        Diagnostics {
            items,
            source,
            line_map,
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = &Diagnostic> {
        self.items.iter()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// The source text these diagnostics refer to.
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Render a rich, source-pointing report (à la `rustc` / `ariadne`): each
    /// diagnostic shows the offending source line with a caret under its span.
    /// Plain text — no ANSI colour — so it reads cleanly in tracebacks, logs,
    /// and notebooks. Falls back to the terse [`Display`](fmt::Display) form if
    /// rendering ever fails.
    pub fn render(&self) -> String {
        use ariadne::{Config, Label, Report, ReportKind, Source};

        const ID: &str = "<stim>";
        let mut buf = Vec::new();
        for d in &self.items {
            let kind = match d.severity {
                Severity::Error => ReportKind::Error,
                Severity::Warning => ReportKind::Warning,
            };
            let span = (ID, d.span.start..d.span.end);
            let report = Report::build(kind, span.clone())
                .with_config(Config::default().with_color(false))
                .with_message(&d.message)
                .with_label(Label::new(span).with_message(&d.message))
                .finish();
            if report
                .write((ID, Source::from(self.source.as_ref())), &mut buf)
                .is_err()
            {
                return self.to_string();
            }
        }
        String::from_utf8(buf).unwrap_or_else(|_| self.to_string())
    }
}

impl fmt::Display for Diagnostics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, d) in self.items.iter().enumerate() {
            if i > 0 {
                writeln!(f)?;
            }
            let (line, col) = d.span.line_col(&self.line_map);
            let sev = match d.severity {
                Severity::Error => "error",
                Severity::Warning => "warning",
            };
            write!(f, "{sev} at line {line}, col {col}: {}", d.message)?;
        }
        Ok(())
    }
}

impl std::error::Error for Diagnostics {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn diagnostics_display_renders_line_col() {
        let diag = Diagnostic::error(
            Span::new(4, 12),
            "unknown-instruction",
            "unknown instruction 'BADINSTR'",
        );
        let diags = Diagnostics::new(vec![diag], Arc::from("X 0\nBADINSTR 0"));
        assert_eq!(
            diags.to_string(),
            "error at line 2, col 1: unknown instruction 'BADINSTR'"
        );
    }

    #[test]
    fn diagnostics_is_empty_when_no_items() {
        assert!(Diagnostics::new(vec![], Arc::from("")).is_empty());
    }

    #[test]
    fn render_points_at_the_source_span() {
        // `X` (byte 8) is the offending target in "H 0\nM 0 X\n".
        let diag = Diagnostic::error(Span::new(8, 9), "invalid-target", "invalid target \"X\"");
        let diags = Diagnostics::new(vec![diag], Arc::from("H 0\nM 0 X\n"));
        let report = diags.render();
        assert!(
            report.contains("M 0 X"),
            "report should show the offending source line:\n{report}"
        );
        assert!(
            report.contains("invalid target"),
            "report should include the message:\n{report}"
        );
    }
}
