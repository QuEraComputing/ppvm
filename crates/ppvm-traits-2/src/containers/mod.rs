// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The graded map algebra `impl`'d **directly on the containers** — no wrapper
//! types. Two backends, one submodule each:
//!
//! - [`Vec<(K, C)>`] ([`coordinate_list`]) — an unsorted coordinate list with
//!   linear-scan `accumulate`; best for small support (the `GeneralizedTableau`
//!   amplitude vector). Requires only `K: Eq + Clone` — it never hashes.
//! - [`HashMap<K, C, IdentityBuildHasher>`](std::collections::HashMap)
//!   ([`hash_join`]) — the hash-join `accumulate`; best for large support
//!   (`PauliSum`). Requires `K:` [`Indexable`](crate::Indexable) (the direct
//!   digest, consumed pass-through through
//!   [`IdentityBuildHasher`](crate::IdentityBuildHasher)).
//!
//! The two files carry the same six impls ([`Support`](crate::Support),
//! [`Accumulate`](crate::Accumulate), [`Scale`](crate::Scale),
//! [`Pair`](crate::Pair), [`Retain`](crate::Retain),
//! [`Multiply`](crate::Multiply)) and share nothing but those trait
//! definitions — the split is by backend, so a change to the cost model of one
//! never has to be read past the other.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Backends are containers;
//! columnar is expressible from day one" and §"The map is a graded algebra over
//! `C[K]`". The module laws these impls realize are machine-checked in
//! `lean/PPVM/Algebra/GradedMap.lean` (`accumulate_comm`/`accumulate_assoc`,
//! `reduce_structural`, `scale_scale`, `overlap_comm`).
//!
//! # Friction: these impls live in `ppvm-traits-2`, not `ppvm-pauli-sum-2`
//!
//! The implementation plan's crate-map table lists "graded traits `impl`'d on
//! `Vec`/`HashMap`" under `ppvm-pauli-sum-2`. Rust's orphan rule forbids that:
//! both the traits ([`Support`](crate::Support) …) and the containers (`Vec`,
//! `HashMap`) would be foreign to `ppvm-pauli-sum-2`, and `(K, C)` carries no
//! local type to anchor the impl. The impls must therefore live in the crate
//! that *owns* the graded traits — here. `ppvm-pauli-sum-2` still owns
//! everything the orphan rule permits it to (`Sum`, `Policy`, the producers,
//! `Clifford for Sum`, the aliases); only these container impls move.

mod coordinate_list;
mod hash_join;

#[cfg(test)]
mod tests;
