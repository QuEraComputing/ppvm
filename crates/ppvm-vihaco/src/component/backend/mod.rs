// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[cfg(feature = "legacy")]
mod legacy;
#[cfg(feature = "traits-2")]
mod traits_2;

#[cfg(feature = "legacy")]
pub use legacy::*;
#[cfg(feature = "traits-2")]
pub use traits_2::*;
