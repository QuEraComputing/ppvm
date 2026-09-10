// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Container algebra, batches, hashing, and traversal preferences.
//!
//! The private backend modules implement these traits for `Vec` and `HashMap`.

mod batch;
mod coordinate_list;
mod graded;
mod hash;
mod hash_join;
mod storage;

pub use batch::{Columnar, KeyBatch, KeyColumn, TermBatch, TermProducer, TermSink};
pub use graded::{Accumulate, Multiply, Pair, Retain, Scale, Support, Trace};
pub use hash::{IdentityBuildHasher, IdentityHasher, Indexable};
pub use storage::RekeyStrategy;
