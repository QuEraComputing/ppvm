/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Data.Complex.Basic

/-!
# PPVM

Formalization scaffolding for the `ppvm` quantum-circuit simulator. The Rust
workspace is the source of truth; this Lean development states and proves the
mathematical spec that the Rust bit-level implementation is meant to refine.

The current, planned, and open targets are tracked in `lean/README.md`. The
first concrete work lives in `PPVM.Pauli` (the symplectic representation of the
Pauli group and its phase cocycle).

This module intentionally starts nearly empty; the `Mathlib.Data.Complex.Basic`
import above is here only to confirm the Mathlib dependency resolves and builds.
-/

namespace PPVM

end PPVM
