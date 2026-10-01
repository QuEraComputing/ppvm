// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Link against cuPauliProp and the CUDA runtime when the `cuda` feature is on.
//!
//! - `CUPAULIPROP_LIB_DIR`: directory containing `libcupauliprop.so.0`
//!   (e.g. `<cuquantum>/lib` from the SDK or the `cupauliprop-cu12` wheel).
//! - `CUDART_LIB_DIR`: directory containing `libcudart.so`; defaults to
//!   `$CUDA_HOME/lib64`, then `/usr/local/cuda/lib64`.
//!
//! The library directories must also be on `LD_LIBRARY_PATH` at run time.

use std::env;

fn main() {
    println!("cargo:rerun-if-env-changed=CUPAULIPROP_LIB_DIR");
    println!("cargo:rerun-if-env-changed=CUDART_LIB_DIR");
    println!("cargo:rerun-if-env-changed=CUDA_HOME");

    if env::var_os("CARGO_FEATURE_CUDA").is_none() {
        return;
    }

    let cupp_dir = env::var("CUPAULIPROP_LIB_DIR")
        .expect("the `cuda` feature needs CUPAULIPROP_LIB_DIR (directory of libcupauliprop.so.0)");
    println!("cargo:rustc-link-search=native={cupp_dir}");
    // The wheel only ships the versioned soname, so link it verbatim.
    println!("cargo:rustc-link-lib=dylib:+verbatim=libcupauliprop.so.0");

    let cudart_dir = env::var("CUDART_LIB_DIR").unwrap_or_else(|_| {
        let home = env::var("CUDA_HOME").unwrap_or_else(|_| "/usr/local/cuda".to_owned());
        format!("{home}/lib64")
    });
    println!("cargo:rustc-link-search=native={cudart_dir}");
    println!("cargo:rustc-link-lib=dylib=cudart");
}
