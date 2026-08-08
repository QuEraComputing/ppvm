// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[cfg(all(feature = "legacy", feature = "traits-2"))]
compile_error!("features `legacy` and `traits-2` are mutually exclusive");
#[cfg(not(any(feature = "legacy", feature = "traits-2")))]
compile_error!("enable exactly one Python backend: `legacy` or `traits-2`");

#[cfg(feature = "legacy")]
mod legacy;
#[cfg(feature = "traits-2")]
mod traits_2;

#[cfg(feature = "legacy")]
pub use legacy::*;
#[cfg(feature = "traits-2")]
pub use traits_2::*;

#[cfg(feature = "legacy")]
pub const NAME: &str = "legacy";
#[cfg(feature = "traits-2")]
pub const NAME: &str = "traits-2";
