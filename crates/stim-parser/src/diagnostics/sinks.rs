// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Two provided `DiagnosticSink` policies.

use crate::diagnostics::{Diagnostic, DiagnosticSink, Flow, Severity};

/// Aborts the current stage on the first `Error` (warnings pass through).
/// Used by the Tier-1 `parse`/`parse_extended` functions.
#[derive(Debug, Default)]
pub struct FailFast {
    items: Vec<Diagnostic>,
    saw_error: bool,
}

impl FailFast {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn saw_error(&self) -> bool {
        self.saw_error
    }

    pub fn into_items(self) -> Vec<Diagnostic> {
        self.items
    }
}

impl DiagnosticSink for FailFast {
    fn emit(&mut self, diagnostic: Diagnostic) -> Flow {
        let is_error = diagnostic.severity == Severity::Error;
        self.items.push(diagnostic);
        if is_error {
            self.saw_error = true;
            Flow::Abort
        } else {
            Flow::Continue
        }
    }
}

/// Never aborts; accumulates every diagnostic for one-pass reporting.
#[derive(Debug, Default)]
pub struct Collect {
    items: Vec<Diagnostic>,
}

impl Collect {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn into_items(self) -> Vec<Diagnostic> {
        self.items
    }
}

impl DiagnosticSink for Collect {
    fn emit(&mut self, diagnostic: Diagnostic) -> Flow {
        self.items.push(diagnostic);
        Flow::Continue
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::{Diagnostic, DiagnosticSink, Flow, Severity, Span};

    fn err() -> Diagnostic {
        Diagnostic::error(Span::new(0, 1), "x", "boom")
    }
    fn warn() -> Diagnostic {
        Diagnostic::warning(Span::new(0, 1), "x", "heads up")
    }

    #[test]
    fn fail_fast_aborts_on_first_error_and_keeps_it() {
        let mut s = FailFast::new();
        assert_eq!(s.emit(warn()), Flow::Continue); // warnings don't abort
        assert_eq!(s.emit(err()), Flow::Abort);
        let items = s.into_items();
        assert_eq!(items.len(), 2);
        assert_eq!(items[1].severity, Severity::Error);
    }

    #[test]
    fn collect_never_aborts_and_gathers_all() {
        let mut s = Collect::new();
        assert_eq!(s.emit(err()), Flow::Continue);
        assert_eq!(s.emit(err()), Flow::Continue);
        assert_eq!(s.into_items().len(), 2);
    }
}
