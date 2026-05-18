// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

/// Trait for computing `trace(self * RHS)`.
///
/// If a type implements `Trace`, a corresponding `TraceBy` implementation
/// is also provided automatically.
pub trait Trace<'a, RHS: 'a> {
    /// Numeric output of the trace.
    type Output;
    /// Compute `tr(self · value)`.
    fn trace(&'a self, value: &'a RHS) -> Self::Output;
}
