// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Paulis on the active axes, and the rotation that collapses onto one axis.

use ppvm_traits_2::{ApplyGateRule, Clifford};

/// A Pauli on the active axes, unsigned: bit `j` of `x` / `z` is the X / Z
/// factor on axis `j`, and both set is Y.
/// `P|α⟩ = i^{|x∧z|} (−1)^{|α∧z|} |α ⊕ x⟩`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AxisPauli {
    pub x: u64,
    pub z: u64,
}

impl AxisPauli {
    pub fn is_identity(&self) -> bool {
        self.x == 0 && self.z == 0
    }

    fn get(&self, j: usize) -> [bool; 2] {
        [self.x >> j & 1 == 1, self.z >> j & 1 == 1]
    }

    fn set(&mut self, j: usize, [x, z]: [bool; 2]) {
        self.x = self.x & !(1 << j) | u64::from(x) << j;
        self.z = self.z & !(1 << j) | u64::from(z) << j;
    }
}

/// An [`AxisPauli`] with a sign, so the gate rules can run on it.
struct Signed {
    p: AxisPauli,
    negative: bool,
}

impl ApplyGateRule for Signed {
    type B = bool;

    fn apply1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        let [x, z] = self.p.get(q);
        let (bits, flip) = rule(x, z);
        self.p.set(q, bits);
        self.negative ^= flip;
    }

    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        let ([xa, za], [xb, zb]) = (self.p.get(a), self.p.get(b));
        let ([xa, za, xb, zb], flip) = rule([xa, za, xb, zb]);
        self.p.set(a, [xa, za]);
        self.p.set(b, [xb, zb]);
        self.negative ^= flip;
    }
}

/// The change of basis `V` of a collapse, shared by the tableau (`C ← C·V`)
/// and the active state (`φ ← V† φ`), as a gate sequence applied in this
/// order: CX from `pivot` to every axis in `fan_out`, CX from every axis in
/// `fan_in` to `pivot`, CZ between `pivot` and every axis in `cz`, then `S`,
/// `H`, `X` on `pivot` where set.
///
/// Built by [`Rotation::collapse`] so that `V† P V = ±Z_pivot` and the
/// measured branch leaves `pivot` in `|0⟩`; one place decides the gates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Rotation {
    pub pivot: usize,
    pub fan_out: u64,
    pub fan_in: u64,
    pub cz: u64,
    pub s: bool,
    pub h: bool,
    pub x: bool,
}

impl Rotation {
    /// The rotation for measuring `p` with outcome `m` (eigenvalue `(−1)^m`):
    /// after `V`, the projection onto that outcome leaves `pivot` in `|0⟩`.
    ///
    /// With an X part, the pivot is its lowest X axis: CX from it clears the
    /// other X bits, CZ clears the remaining Z bits, and `H` (`S`, `H` for a
    /// Y) turns the pivot's X into Z. A Z-only `p` collects its Z bits onto its
    /// lowest Z axis with CX into it.
    pub fn collapse(p: AxisPauli, m: bool) -> Self {
        assert!(!p.is_identity(), "measuring the identity");
        let mut v = if p.x != 0 {
            let pivot = p.x.trailing_zeros() as usize;
            let fan_out = p.x & !(1 << pivot);
            let cx_only = Self {
                pivot,
                fan_out,
                fan_in: 0,
                cz: 0,
                s: false,
                h: false,
                x: false,
            };
            // After the CXs, the Z bits left on other axes are the CZ targets,
            // and a Z left on the pivot makes it a Y.
            let after_cx = cx_only.conjugate(p).0;
            Self {
                cz: after_cx.z & !(1 << pivot),
                s: after_cx.z >> pivot & 1 == 1,
                h: true,
                ..cx_only
            }
        } else {
            let pivot = p.z.trailing_zeros() as usize;
            Self {
                pivot,
                fan_out: 0,
                fan_in: p.z & !(1 << pivot),
                cz: 0,
                s: false,
                h: false,
                x: false,
            }
        };
        let (reduced, negative) = v.conjugate(p);
        debug_assert_eq!(
            reduced,
            AxisPauli {
                x: 0,
                z: 1 << v.pivot
            }
        );
        // The projection leaves the pivot in |b⟩ with (−1)^b = ±(−1)^m; X maps it to |0⟩.
        v.x = m != negative;
        v
    }

    /// Apply the gates of `V` to `c`, in order.
    pub(crate) fn apply_gates(&self, c: &mut impl Clifford) {
        let p = self.pivot;
        for t in bits(self.fan_out) {
            c.cnot(p, t);
        }
        for ctrl in bits(self.fan_in) {
            c.cnot(ctrl, p);
        }
        for t in bits(self.cz) {
            c.cz(p, t);
        }
        if self.s {
            c.s(p);
        }
        if self.h {
            c.h(p);
        }
        if self.x {
            c.x(p);
        }
    }

    /// `V† p V` and whether it picks up a sign, by running the gate rules on
    /// `p` in gate order.
    pub fn conjugate(&self, p: AxisPauli) -> (AxisPauli, bool) {
        let mut signed = Signed { p, negative: false };
        self.apply_gates(&mut signed);
        (signed.p, signed.negative)
    }
}

/// The set bits of `mask`, ascending.
pub(super) fn bits(mut mask: u64) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        (mask != 0).then(|| {
            let j = mask.trailing_zeros() as usize;
            mask &= mask - 1;
            j
        })
    })
}
