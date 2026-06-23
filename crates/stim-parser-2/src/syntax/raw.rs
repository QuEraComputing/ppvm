// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Untyped syntactic tree produced by the chumsky grammar, before
//! table-driven validation. Crate-internal plumbing between the grammar
//! and the validate pass.

use chumsky::span::SimpleSpan;

use crate::ast::shared::Tag;

pub(crate) type RawSyntaxTree = Vec<RawSyntaxNode>;

#[derive(Debug, Clone)]
pub(crate) enum RawSyntaxNode {
    Instruction {
        name: String,
        tags: Vec<Tag>,
        args: Vec<f64>,
        targets: Vec<RawTarget>,
        span: SimpleSpan<usize>,
    },
    Repeat {
        count: u64,
        body: Vec<RawSyntaxNode>,
        span: SimpleSpan<usize>,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct RawTarget {
    pub text: String,
    pub span: SimpleSpan<usize>,
}
