// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Let plain `cargo build` link this PyO3 extension module on macOS.
//!
//! With the `extension-module` feature, PyO3 deliberately does not link
//! libpython: the CPython symbols (`_PyExc_*`, `_Py_*`, ...) are supplied by
//! the host interpreter when the module is imported. maturin passes
//! `-undefined dynamic_lookup` so the macOS linker tolerates those undefined
//! symbols; plain `cargo build`/`cargo test` does not, and fails to link the
//! cdylib. Emit the same flag here so the Rust-centric workflows
//! (`cargo build --workspace`, etc.) work without going through maturin.
//!
//! Scope notes:
//! - `rustc-link-arg-cdylib` applies only to this crate's cdylib, so it does
//!   not weaken undefined-symbol detection for the rest of the workspace.
//! - maturin still works: it sets the same flag, and a repeated
//!   `-undefined dynamic_lookup` is a no-op for the macOS linker.
//! - Linux/Windows need nothing extra here — a Linux shared object already
//!   permits undefined symbols, and PyO3's own build logic handles Windows.

fn main() {
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos") {
        println!("cargo:rustc-link-arg-cdylib=-undefined");
        println!("cargo:rustc-link-arg-cdylib=dynamic_lookup");
    }
}
