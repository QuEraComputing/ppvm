// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Composable ratatui components + app state for the `ppvm` TUI. Terminal-
//! agnostic: no code here owns a terminal or runs an event loop, so the
//! `Widget` components and `AppState` can be embedded in another ratatui app.

#[cfg(all(feature = "legacy", feature = "traits-2"))]
compile_error!("features `legacy` and `traits-2` are mutually exclusive");
#[cfg(not(any(feature = "legacy", feature = "traits-2")))]
compile_error!("enable exactly one ppvm-tui backend: `legacy` or `traits-2`");

pub mod app;
pub mod codeview;
pub mod command;
pub mod editor;
pub mod widgets;

pub use app::AppState;
