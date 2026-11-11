// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Paulis on the active axes, and the rotation that collapses onto one axis.

use ppvm_pauli_kernels::heisenberg;

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

    /// `V† p V` and whether it picks up a sign, by running the gate rules on
    /// `p`'s bits in gate order.
    pub fn conjugate(&self, mut p: AxisPauli) -> (AxisPauli, bool) {
        let mut negative = false;
        let mut two =
            |p: &mut AxisPauli, c: usize, t: usize, rule: fn([bool; 4]) -> ([bool; 4], bool)| {
                let ([xc, zc], [xt, zt]) = (p.get(c), p.get(t));
                let ([xc, zc, xt, zt], flip) = rule([xc, zc, xt, zt]);
                p.set(c, [xc, zc]);
                p.set(t, [xt, zt]);
                negative ^= flip;
            };
        let p_ = self.pivot;
        for t in bits(self.fan_out) {
            two(&mut p, p_, t, heisenberg::cnot);
        }
        for c in bits(self.fan_in) {
            two(&mut p, c, p_, heisenberg::cnot);
        }
        for t in bits(self.cz) {
            two(&mut p, p_, t, heisenberg::cz);
        }
        let one = |p: &mut AxisPauli, rule: fn(bool, bool) -> ([bool; 2], bool)| {
            let [x, z] = p.get(p_);
            let (bits, flip) = rule(x, z);
            p.set(p_, bits);
            flip
        };
        for (on, rule) in [
            (self.s, heisenberg::s as fn(bool, bool) -> ([bool; 2], bool)),
            (self.h, heisenberg::h),
            (self.x, heisenberg::x),
        ] {
            if on {
                negative ^= one(&mut p, rule);
            }
        }
        (p, negative)
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
