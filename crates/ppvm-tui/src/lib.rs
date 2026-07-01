// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Composable ratatui components + app state for the `ppvm` TUI. Terminal-
//! agnostic: no code here owns a terminal or runs an event loop, so the
//! `Widget` components and `AppState` can be embedded in another ratatui app.

pub mod app;
pub mod codeview;
pub mod command;

pub use app::AppState;
