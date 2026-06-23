// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Stage 1: pure syntax. The chumsky grammar produces a `RawSyntaxTree`;
//! `run_on_parser_stack` runs the recursive grammar on an oversized stack.

mod grammar;
mod raw;

#[expect(unused_imports, reason = "called by pipeline in Task 11")]
pub(crate) use grammar::program_parser;
#[expect(unused_imports, reason = "used by grammar/pipeline in Tasks 10-11")]
pub(crate) use raw::{RawSyntaxNode, RawSyntaxTree, RawTarget};

/// Stack size for the dedicated parsing thread. The chumsky grammar is built
/// around `recursive(...)`, which descends into REPEAT bodies via recursive
/// parser calls; on the default thread stack, deeply nested programs overflow.
/// Running on an oversized dedicated stack supports thousands of nested
/// REPEATs without rewriting the grammar.
#[cfg(not(target_arch = "wasm32"))]
const PARSER_STACK_SIZE: usize = 16 * 1024 * 1024;

/// Run `f` with a large stack for the recursive chumsky grammar. On targets
/// with OS threads this spawns a dedicated [`PARSER_STACK_SIZE`]-byte thread;
/// on `wasm32` (no `std::thread`) `f` runs inline on the caller's stack.
#[expect(dead_code, reason = "used by grammar/pipeline in Tasks 10-11")]
pub(crate) fn run_on_parser_stack<R, F>(f: F) -> R
where
    R: Send,
    F: FnOnce() -> R + Send,
{
    #[cfg(target_arch = "wasm32")]
    {
        f()
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        std::thread::scope(|s| {
            let handle = std::thread::Builder::new()
                .stack_size(PARSER_STACK_SIZE)
                .name("stim-parser".to_string())
                .spawn_scoped(s, f)
                .expect("failed to spawn parser thread");
            match handle.join() {
                Ok(value) => value,
                Err(payload) => std::panic::resume_unwind(payload),
            }
        })
    }
}
