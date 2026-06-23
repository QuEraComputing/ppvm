// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod hashmap;

#[cfg(all(feature = "dashmap", not(target_arch = "wasm32")))]
mod dashmap;
