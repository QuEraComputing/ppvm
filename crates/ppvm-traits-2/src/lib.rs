// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `ppvm-traits-2`: the trait foundation for the second trait-system experiment.
//!
//! This crate holds only trait *definitions* and small leaf types — the split
//! `Coefficient`/`Angle`, `Word`/`Indexable`/`PauliBits`, the Pauli algebra
//! traits, the gate/noise traits, the graded map layers, the algebra
//! capabilities (`KeyProduct`/`ImaginaryUnit`/`Conjugate`), and the batch
//! contract. See [`docs/design/traits-2-configuration-and-hashing.md`] and
//! [`docs/design/traits-2-implementation-plan.md`].
//!
//! At Phase 0 (scaffolding) this crate is an intentional stub; the trait modules
//! land in Phase 1.
