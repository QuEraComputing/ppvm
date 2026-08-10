// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`Clifford`] propagation for a Pauli-keyed [`Sum`] (the `PauliSum` alias),
//! draining each conjugation's `±1` sign to the coefficient. Two fast paths:
//!
//! * The **pure-sign** single-qubit gates `X`/`Y`/`Z` leave the Pauli word fixed
//!   (`XPX = (−1)^z P`, `YPY = (−1)^{x⊕z} P`, `ZPZ = (−1)^x P` — proven as the
//!   group conjugation `G·P·G⁻¹` in `lean/PPVM/Pauli/Conjugation.lean`,
//!   `conjX`/`conjY`/`conjZ`), so they take the in-place
//!   [`Sum::flip_sign_by_key`](crate::Sum) path — walk the existing entries and
//!   scale each coefficient by the `±1` its own bits demand, **no** map rebuild,
//!   no key movement, no reallocation (the old crate's in-place `scale`, restored).
//! * The **word-changing** gates re-key every term via the move-based
//!   [`Sum::rekey_bijective`](crate::Sum) fast path and direct audited bit/sign
//!   formulas.
//!
//! # The Clifford-sign subtlety
//!
//! A Clifford conjugates each Pauli to `±` another Pauli. The **bare**
//! [`PauliWord`] Clifford (the blanket over `SymplecticColumns` + `PhaseTrack`)
//! is *bit-only*: its `PhaseTrack` is a no-op, so it computes the resulting
//! Pauli's bits but **drops the `±` sign**. The design's generic
//! `impl<S> Clifford for Sum<S, P> where S::Key: Clifford` would therefore lose
//! every conjugation sign for a `PauliSum` — a silent correctness bug.
//!
//! So this impl does **not** dispatch to the key's bare `Clifford`. It applies
//! the same pre-mutation bit/sign formulas as the audited fused phased-word
//! kernels, drains the resulting real sign to the coefficient, and builds each
//! replacement key directly. The two-site path batches all toggles into one word
//! rebuild so packed words refresh their structural hash once.
//!
//! A Clifford re-key is a bijection, so colliding re-keyed terms never occur —
//! which is what lets this path take the plain-`insert`
//! [`Sum::rekey_bijective`] fast path rather than [`Sum::apply`]'s batch
//! round-trip. Nothing here runs `reduce` (a `±1` sign cannot zero a
//! coefficient, and a zero term must survive regardless — old has no `reduce`)
//! and nothing here truncates (caller-driven, [`Sum::truncate`]). The bijection is
//! machine-checked at the phase-stripped word level: each generator's `Sp(2n, 2)`
//! bit map is an involution, hence bijective (`lean/PPVM/Pauli/Symplectic.lean`
//! `hAct_involutive`/`sAct_involutive`/`cnotAct_involutive`/`czAct_involutive`,
//! `*_bijective`).
//!
//! Composing that bijectivity with the `±1`-sign reality above gives the
//! semantic guarantee that ties this Clifford path to [`Sum::overlap`]: the
//! Heisenberg re-key **preserves the Hilbert–Schmidt trace pairing**
//! (`overlap(conj_G A, conj_G B) = overlap(A, B)`), because the drained sign
//! squares out (`s_P² = 1`) while the bijection only permutes the summands. This
//! is machine-checked as `clifford_conjugation_preserves_overlap` in
//! `lean/PPVM/Algebra/GradedMap.lean` (over `overlap_eq_fintype_sum`).
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits"
//! (`PauliSum` is deliberately *not* a `BlanketClifford` implementer: "the sum
//! applies the one-row action pointwise and drains each term's phase delta to its
//! coefficient") and §"apply". The conjugation signs are machine-checked in
//! `lean/PPVM/Pauli/Conjugation.lean` (`conjH_Y`: `HYH = −Y`, `conjSdag_sign`,
//! `conjCNOT_sign`, `conjCZ_sign`); the underlying bit maps are the `Sp(2n, 2)`
//! isometries of `lean/PPVM/Pauli/Symplectic.lean`.

use ppvm_traits_2::{
    Accumulate, Clifford, CliffordBatch, CliffordExtensions, CliffordExtensionsBatch, Coefficient,
    Indexable, PauliBits, Retain, Word,
};

use crate::policy::Policy;
use crate::store::{RekeyBijective, SignFlipByKey, StoreAlloc};
use crate::sum::Sum;

/// `CY`'s conjugation rule on one term: `None` when a lost qubit makes the gate
/// a no-op, else `(toggle z_control, toggle both target planes, negate)`.
///
/// It is a named `#[inline(always)]` function rather than the body of the two
/// closures in [`CliffordExtensions::cy`] only so the two arms of that gate
/// share one audited copy of the rule. It is **not** what keeps the re-key loop
/// fast — see the `#[inline(never)]` on [`CliffordExtensions::cy`] for that.
#[inline(always)]
fn cy_toggles<W>(k: &W, control: usize, target: usize) -> Option<(bool, bool, bool)>
where
    W: Word + PauliBits,
{
    if k.is_lost(control) || k.is_lost(target) {
        return None;
    }
    let control_code = k.pauli_code(control);
    let target_code = k.pauli_code(target);
    let xc = control_code & 1 != 0;
    let zc = control_code & 2 != 0;
    let xt = target_code & 1 != 0;
    let zt = target_code & 2 != 0;
    let target_mix = xt ^ zt;
    Some((target_mix, xc, xc & target_mix & !(zc ^ zt)))
}

impl<S, P, W, C> Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + StoreAlloc + Retain<W, C> + RekeyBijective<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient,
    P: Policy<W, C>,
{
    /// Re-key every term by a **direct bit rewrite** of its own word: `f` mutates
    /// the moved key's X/Z bits at the gate's qubit and returns whether the
    /// conjugation sign is `−1`.
    ///
    /// This is the specialization the *single-qubit* word-changing Cliffords take,
    /// and it is a named baseline feature, not an incidental difference: old
    /// implements `h`/`s`/`s_dag`/`sqrt_*` as `map_add` closures
    /// that read `get_xbit`/`get_zbit`, write the swapped/XOR'd bits, `rehash()`,
    /// and pick the sign from a two-bit test, with the source comment
    ///
    /// > "Single-qubit gates are specialised to bit-level updates so we avoid
    /// > round-tripping through `PhasedPauliWord`; the two-qubit gates still go
    /// > through the macro"
    ///
    /// (`ppvm-pauli-sum/src/sum/clifford.rs`). Only `cnot`/`cz`/`cy` — whose sign
    /// rules genuinely want the phased word's audited fused kernel — keep the
    /// wrapper. The sign is applied as a `Neg` on the *moved* coefficient rather
    /// than [`Coefficient::mul_sign`], which takes `&self` and therefore clones on
    /// a heap-owning coefficient ring, and the `+1` case does not touch the
    /// coefficient at all (old's `v.clone()` vs `-v.clone()` branch).
    #[inline(always)]
    fn rekey_bits<F>(&mut self, f: F)
    where
        F: Fn(&mut W) -> bool + Send + Sync,
    {
        self.rekey_bijective(|mut k: W, c: C| {
            let negate = f(&mut k);
            if negate { (k, -c) } else { (k, c) }
        });
    }

    /// Re-key by mutating a moved source word.
    #[inline(always)]
    fn rekey_owned<F>(&mut self, f: F)
    where
        F: Fn(W) -> (W, bool) + Send + Sync,
    {
        self.rekey_bijective(|key: W, coeff: C| {
            let (key, negate) = f(key);
            if negate { (key, -coeff) } else { (key, coeff) }
        });
    }

    /// Re-key from a borrowed source word, allowing packed implementations to
    /// build the replacement planes directly without cloning a digest cache that
    /// the two-site mutation would immediately invalidate.
    #[inline(always)]
    fn rekey_ref<F>(&mut self, f: F)
    where
        F: Fn(&W) -> (W, bool) + Send + Sync,
    {
        self.rekey_bijective_ref(|key: &W, coeff: &C| {
            let (key, negate) = f(key);
            let coeff = coeff.clone();
            if negate { (key, -coeff) } else { (key, coeff) }
        });
    }
}

/// Clifford propagation on a Pauli-keyed `Sum`. Each gate re-keys the whole
/// support pointwise, folding the conjugation sign into the coefficient.
impl<S, P, W, C> Clifford for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C>
        + StoreAlloc
        + Retain<W, C>
        + RekeyBijective<W, C>
        + SignFlipByKey<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient,
    P: Policy<W, C>,
{
    /// `X` conjugation is a **pure sign**: `XPX = (−1)^z P`. The word is fixed, so
    /// this takes the in-place [`Sum::flip_sign_by_key`] fast path — flipping each
    /// term's coefficient iff its `z` bit at `qubit` is set — instead of rebuilding
    /// the map. Sign matches [`PhaseTrack::x_phase`](ppvm_traits_2::PhaseTrack) and
    /// the phased word's fused `Phased::x` (`ppvm-phased-pauli-word-2`).
    #[inline]
    fn x(&mut self, qubit: usize) {
        self.flip_sign_by_key(move |k| {
            if !k.is_lost(qubit) && k.z_bit(qubit) {
                -1
            } else {
                1
            }
        });
    }

    /// `Y` conjugation is a **pure sign**: `YPY = (−1)^{x⊕z} P`. Word fixed → the
    /// in-place fast path, flipping iff `x ⊕ z` at `qubit`. Sign matches
    /// [`PhaseTrack::y_phase`](ppvm_traits_2::PhaseTrack) and the phased word's
    /// fused `Phased::y`.
    #[inline]
    fn y(&mut self, qubit: usize) {
        self.flip_sign_by_key(move |k| {
            if !k.is_lost(qubit) && (k.x_bit(qubit) ^ k.z_bit(qubit)) {
                -1
            } else {
                1
            }
        });
    }

    /// `Z` conjugation is a **pure sign**: `ZPZ = (−1)^x P`. Word fixed → the
    /// in-place fast path, flipping iff the `x` bit at `qubit` is set. Sign matches
    /// [`PhaseTrack::z_phase`](ppvm_traits_2::PhaseTrack) and the phased word's
    /// fused `Phased::z`.
    #[inline]
    fn z(&mut self, qubit: usize) {
        self.flip_sign_by_key(move |k| {
            if !k.is_lost(qubit) && k.x_bit(qubit) {
                -1
            } else {
                1
            }
        });
    }

    /// `H` **swaps** the X and Z bits; the sign flips only for `Y` (both bits
    /// set): `HXH = Z`, `HZH = X`, `HYH = −Y`. A direct bit rewrite
    /// ([`rekey_bits`](Sum::rekey_bits)), as in old — not a [`Phased`]
    /// round-trip.
    #[inline(always)]
    fn h(&mut self, qubit: usize) {
        self.rekey_bits(move |k| {
            if k.is_lost(qubit) {
                return false;
            }
            let x = k.x_bit(qubit);
            let z = k.z_bit(qubit);
            k.set_xz_bits(qubit, z, x);
            x & z
        });
    }

    /// `S`: `z ← x ⊕ z`; the sign flips only for `X` (x set, z clear):
    /// `SXS† = −Y`, `SYS† = X`, `SZS† = Z`.
    #[inline(always)]
    fn s(&mut self, qubit: usize) {
        self.rekey_bits(move |k| {
            if k.is_lost(qubit) {
                return false;
            }
            let x = k.x_bit(qubit);
            let z = k.z_bit(qubit);
            k.set_z_bit(qubit, x ^ z);
            x & !z
        });
    }

    #[inline(always)]
    fn cnot(&mut self, control: usize, target: usize) {
        if W::PREFER_BORROWED_REKEY {
            self.rekey_ref(move |k| {
                if k.is_lost(control) || k.is_lost(target) {
                    return (k.clone(), false);
                }
                let control_code = k.pauli_code(control);
                let target_code = k.pauli_code(target);
                let xc = control_code & 1 != 0;
                let zc = control_code & 2 != 0;
                let xt = target_code & 1 != 0;
                let zt = target_code & 2 != 0;
                let out = k.toggled_bits2(control, false, zt, target, xc, false);
                (out, xc & zt & (xt == zc))
            });
            return;
        }
        self.rekey_owned(move |k| {
            if k.is_lost(control) || k.is_lost(target) {
                return (k, false);
            }
            let control_code = k.pauli_code(control);
            let target_code = k.pauli_code(target);
            let xc = control_code & 1 != 0;
            let zc = control_code & 2 != 0;
            let xt = target_code & 1 != 0;
            let zt = target_code & 2 != 0;
            let out = k.into_toggled_bits2(control, false, zt, target, xc, false);
            (out, xc & zt & (xt == zc))
        });
    }

    #[inline(always)]
    fn cz(&mut self, qubit0: usize, qubit1: usize) {
        if W::PREFER_BORROWED_REKEY {
            self.rekey_ref(move |k| {
                if k.is_lost(qubit0) || k.is_lost(qubit1) {
                    return (k.clone(), false);
                }
                let code0 = k.pauli_code(qubit0);
                let code1 = k.pauli_code(qubit1);
                let x0 = code0 & 1 != 0;
                let z0 = code0 & 2 != 0;
                let x1 = code1 & 1 != 0;
                let z1 = code1 & 2 != 0;
                let out = k.toggled_bits2(qubit0, false, x1, qubit1, false, x0);
                (out, x0 & x1 & (z0 ^ z1))
            });
            return;
        }
        self.rekey_owned(move |k| {
            if k.is_lost(qubit0) || k.is_lost(qubit1) {
                return (k, false);
            }
            let code0 = k.pauli_code(qubit0);
            let code1 = k.pauli_code(qubit1);
            let x0 = code0 & 1 != 0;
            let z0 = code0 & 2 != 0;
            let x1 = code1 & 1 != 0;
            let z1 = code1 & 2 != 0;
            let out = k.into_toggled_bits2(qubit0, false, x1, qubit1, false, x0);
            (out, x0 & x1 & (z0 ^ z1))
        });
    }
}

/// The extended Clifford set, ported from `ppvm-pauli-sum/src/sum/clifford.rs`
/// with old's own split: the five single-qubit gates are direct bit rewrites
/// ([`Sum::rekey_bits`]) and only the two-qubit `cy` goes through the phased
/// word.
impl<S, P, W, C> CliffordExtensions for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C>
        + StoreAlloc
        + Retain<W, C>
        + RekeyBijective<W, C>
        + SignFlipByKey<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient,
    P: Policy<W, C>,
{
    /// `S†`: the same bit map as `S` (`z ← x ⊕ z`); the sign flips for `Y`
    /// (both bits set): `S†XS = Y`, `S†YS = −X`, `S†ZS = Z`.
    #[inline(always)]
    fn s_dag(&mut self, qubit: usize) {
        self.rekey_bits(move |k| {
            if k.is_lost(qubit) {
                return false;
            }
            let x = k.x_bit(qubit);
            let z = k.z_bit(qubit);
            k.set_z_bit(qubit, x ^ z);
            x & z
        });
    }

    /// `√X`: `x ← x ⊕ z`; the sign flips for `Y`: `X ↦ X`, `Y ↦ −Z`, `Z ↦ Y`.
    #[inline(always)]
    fn sqrt_x(&mut self, qubit: usize) {
        self.rekey_bits(move |k| {
            if k.is_lost(qubit) {
                return false;
            }
            let x = k.x_bit(qubit);
            let z = k.z_bit(qubit);
            k.set_x_bit(qubit, x ^ z);
            x & z
        });
    }

    /// `(√X)†`: the same bit map as `√X`; the sign flips for `Z` (z set, x
    /// clear): `X ↦ X`, `Y ↦ Z`, `Z ↦ −Y`.
    #[inline(always)]
    fn sqrt_x_dag(&mut self, qubit: usize) {
        self.rekey_bits(move |k| {
            if k.is_lost(qubit) {
                return false;
            }
            let x = k.x_bit(qubit);
            let z = k.z_bit(qubit);
            k.set_x_bit(qubit, x ^ z);
            !x & z
        });
    }

    /// `√Y`: swap the X and Z bits; the sign flips for `Z`: `X ↦ Z`, `Y ↦ Y`,
    /// `Z ↦ −X`.
    #[inline(always)]
    fn sqrt_y(&mut self, qubit: usize) {
        self.rekey_bits(move |k| {
            if k.is_lost(qubit) {
                return false;
            }
            let x = k.x_bit(qubit);
            let z = k.z_bit(qubit);
            k.set_xz_bits(qubit, z, x);
            !x & z
        });
    }

    /// `(√Y)†`: swap the X and Z bits; the sign flips for `X`: `X ↦ −Z`,
    /// `Y ↦ Y`, `Z ↦ X`.
    #[inline(always)]
    fn sqrt_y_dag(&mut self, qubit: usize) {
        self.rekey_bits(move |k| {
            if k.is_lost(qubit) {
                return false;
            }
            let x = k.x_bit(qubit);
            let z = k.z_bit(qubit);
            k.set_xz_bits(qubit, z, x);
            x & !z
        });
    }

    /// `CY`: direct two-site rewrite using the same pre-mutation sign predicate
    /// as the audited phased-word kernel. The bit/sign rule lives in
    /// [`cy_toggles`] so both arms share one audited copy.
    ///
    /// # Why `#[inline(never)]` on a gate whose siblings are `#[inline(always)]`
    ///
    /// Every other gate here is `#[inline(always)]`, and so are `rekey_owned`
    /// and `rekey_bijective` — so the whole 192-term re-key loop lands in the
    /// caller. What is *not* forced is the per-term closure the loop calls:
    /// closures carry no inline attribute, so joining the loop is left to LLVM's
    /// cost model, and `CY`'s body — one toggle more than `CNOT`/`CZ`, plus the
    /// eager digest rebuild — sat on the wrong side of the threshold. It was the
    /// **only** re-key closure in the whole conformance binary that compiled out
    /// of line, and out of line is expensive out of all proportion to the call:
    /// a 32-byte `PauliWord` exceeds the AArch64 register-return limit, so the
    /// planes are copied to the stack, the bit toggles become byte
    /// read-modify-writes *in memory*, the digest's 16-byte plane load then
    /// overlaps those byte stores (a store-forwarding stall), and the result
    /// goes back out through an indirect `(word, coeff)` struct return — per
    /// term. `pauli_sum_surface/clifford/cy` measured 1.11x against old and
    /// stayed there under both `-Cllvm-args` layout controls, while `cz`/`cnot`
    /// — same machinery, closures LLVM did inline — sat at parity.
    ///
    /// Shrinking the body to sneak back under the threshold is not a fix that
    /// stays fixed; neither is `#[inline(always)]` on a helper the closure calls
    /// (that inlines *into* the closure and makes it bigger — it was tried, and
    /// left the closure out of line). What decides it is the **number of call
    /// sites**: LLVM always inlines the last call to an internal function.
    /// `#[inline(always)]` here cloned the loop into `cy`'s every caller
    /// (`zcy`, `cy_many`, each user call site), giving the closure several call
    /// sites and no such bonus. Pinning `cy` out of line leaves exactly one copy
    /// of the loop, hence exactly one call to the closure, which LLVM then
    /// folds in unconditionally. The gate costs one extra call per *gate* — not
    /// per term — and buys back a call, a frame, two 32-byte plane copies and a
    /// store-forwarding stall on every term. Four launches each, ratios vs old:
    /// `clifford/cy` 1.107x → **0.896x** (0.891–0.900), `clifford/zcy_alias`
    /// 1.107x → **0.887x**, `clifford_batch/cy` 0.848x → **0.785x**; the
    /// `-align-all-functions=6` control agrees (0.901x / 0.898x / 0.866x), and
    /// `cz`/`cnot`/`h`/`s` are unmoved.
    #[inline(never)]
    fn cy(&mut self, control: usize, target: usize) {
        if W::PREFER_BORROWED_REKEY {
            self.rekey_ref(move |k| {
                let Some((toggle_z_c, toggle_t, negate)) = cy_toggles(k, control, target) else {
                    return (k.clone(), false);
                };
                (
                    k.toggled_bits2(control, false, toggle_z_c, target, toggle_t, toggle_t),
                    negate,
                )
            });
            return;
        }
        self.rekey_owned(move |k| {
            let Some((toggle_z_c, toggle_t, negate)) = cy_toggles(&k, control, target) else {
                return (k, false);
            };
            (
                k.into_toggled_bits2(control, false, toggle_z_c, target, toggle_t, toggle_t),
                negate,
            )
        });
    }

    #[inline(always)]
    fn zcy(&mut self, control: usize, target: usize) {
        self.cy(control, target);
    }
}

/// Old ships the batch forms as empty `impl`s over the loop defaults
/// (`impl<T: Config> CliffordBatch for PauliSum<T> {}`); same here.
impl<S, P, W, C> CliffordBatch for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C>
        + StoreAlloc
        + Retain<W, C>
        + RekeyBijective<W, C>
        + SignFlipByKey<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient,
    P: Policy<W, C>,
{
    #[inline(always)]
    fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
        if W::PREFER_BORROWED_REKEY {
            for &(control, target) in pairs {
                self.cnot(control, target);
            }
            return;
        }
        self.rekey_owned(move |mut key| {
            let mut negate = false;
            for &(control, target) in pairs {
                if key.is_lost(control) || key.is_lost(target) {
                    continue;
                }
                let control_code = key.pauli_code(control);
                let target_code = key.pauli_code(target);
                let xc = control_code & 1 != 0;
                let zc = control_code & 2 != 0;
                let xt = target_code & 1 != 0;
                let zt = target_code & 2 != 0;
                negate ^= xc & zt & (xt == zc);
                key = key.into_toggled_bits2(control, false, zt, target, xc, false);
            }
            (key, negate)
        });
    }
}

impl<S, P, W, C> CliffordExtensionsBatch for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C>
        + StoreAlloc
        + Retain<W, C>
        + RekeyBijective<W, C>
        + SignFlipByKey<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient,
    P: Policy<W, C>,
{
}
