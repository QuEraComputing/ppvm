// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`RotationOne`] propagation for a Pauli-keyed [`Sum`] (the `PauliSum` alias):
//! the non-Clifford **branch**.
//!
//! A single-qubit rotation `exp(−i·θ/2·G)` (with axis Pauli `G` = `X`/`Y`/`Z` for
//! `rx`/`ry`/`rz`) conjugates each stored term `(P, c)` to
//! `c·cosθ·P + c·sinθ·(iGP)` when `G` anticommutes with `P` at the qubit, and
//! leaves it unchanged when they commute. This is a genuine **fan-out** (1 or 2
//! output terms per input), but the fan-out is *lopsided*: the `cos·P` diagonal
//! keeps its key, so only the `sin·(iGP)` branch is a genuinely-new term. Rather
//! than pay [`Sum::apply`]'s batch round-trip — which clones every key, resets the
//! map, and re-hashes the *whole* `2N` fan-out including the untouched diagonals —
//! this drives the fused in-place [`Sum::rotate_in_place`] fast path: it scales
//! each diagonal coefficient by `cosθ` **where it sits** (cached hash intact, no
//! bucket move) and hashes/merges only the ≤`N` branch terms, aggregating any
//! collision (a branch's `iGP` may land on another term). An exact cancellation
//! is **kept** — the merge is old's `add_assign`, which leaves a `0.0` entry in
//! the support — and neither `reduce` nor the policy runs here: truncation is
//! caller-driven ([`Sum::truncate`]) and canonicalization is caller-driven
//! ([`Sum::reduce`]). This mirrors the old crate's single-pass `map_insert`
//! (`ppvm-pauli-sum::sum::rot1`), the pure-sign
//! [`Sum::flip_sign_by_key`], and the diagonal [`Sum::scale_by_key`] paths.
//!
//! The `iGP` branch key is a **real** Pauli — the single-qubit anticommuting
//! product `GP = ±iP'` carries one factor of `i`, which the leading `i` of `iGP`
//! cancels — so the branch coefficient stays real and its `±1` sign is drained
//! through [`Coefficient::mul_sign`] exactly as the old crate's
//! `sin.mul_sign(eps)` (`ppvm-pauli-sum::sum::rot1`).
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Behavioral traits"
//! (`RotationOne`) and §"Every gate is a producer feeding `accumulate`". The
//! branch is a norm-preserving, angle-additive 2-D rotation on the coefficient
//! pair whose new key is genuinely new, machine-checked in
//! `lean/PPVM/Instantiations/Rotation.lean` (`anticommute_new_key`, `rot_norm_sq`,
//! `rot_rot`).

use ppvm_traits_2::{
    Accumulate, Angle, Coefficient, Indexable, Pauli, PauliBits, Retain, RotXY, RotationOne,
    RotationTwo, Word,
};

use crate::store::RotateInPlace;
use crate::sum::Sum;

/// The `[x, z]`-indexed axis alphabet: `(z << 1) | x` ↦ `I`, `X`, `Z`, `Y`.
///
/// Old's `const PAULIS: [Pauli; 4] = [I, X, Z, Y]` (`ppvm-pauli-sum/src/sum/rot2.rs`),
/// which doubled as the old `Pauli` enum's discriminant order. The `-2` `Pauli`
/// is declared `I, X, Y, Z`, so the table is written out rather than cast from
/// the discriminant.
const AXIS_PAULIS: [Pauli; 4] = [Pauli::I, Pauli::X, Pauli::Z, Pauli::Y];

/// 2-bit Pauli code — `00` `I`, `01` `X`, `10` `Z`, `11` `Y` (x in bit 0, z in
/// bit 1) — the encoding [`levi_civita`] is written against.
#[inline(always)]
fn pauli_code(p: Pauli) -> u8 {
    match p {
        Pauli::I => 0b00,
        Pauli::X => 0b01,
        Pauli::Z => 0b10,
        Pauli::Y => 0b11,
    }
}

/// `ε, k` with `−i·[P_i, P_j]/2 = ε · P_k`, over the 2-bit Pauli code
/// (`00` `I`, `01` `X`, `10` `Z`, `11` `Y`). A commuting pair yields `(0, 0)`.
///
/// Ported **verbatim** (branch-free bit tricks included) from
/// `ppvm-pauli-sum/src/sum/rot1.rs::levi_civita`; it is the generic
/// [`RotationOne::rotate_1`] kernel's table, kept so the generic path this crate
/// exposes is the same function old exposed.
#[inline]
pub(crate) fn levi_civita(i: u8, j: u8) -> (i8, u8) {
    let k = i ^ j; // third Pauli by XOR; 0 when i == j

    // commute ⇔ i == 0 OR j == 0 OR k == 0 (no false positives)
    let commute = ((i == 0) | (j == 0) | (k == 0)) as u8;

    #[inline]
    fn rank(p: u8) -> u8 {
        let b1 = p >> 1;
        (b1 << 1).wrapping_sub(b1 & (p & 1)) // 0, 1, 2 for X, Y, Z
    }

    let ri = rank(i);
    let rj = rank(j);

    // diff = (rj − ri) mod 3, without a modulus
    let mut diff = rj.wrapping_sub(ri).wrapping_add(3);
    diff -= 3 & (0u8.wrapping_sub(diff >> 2));

    // +1 when diff == 1, −1 when diff == 2
    let eps_raw = 1i8 - 2 * ((diff >> 1) as i8);

    let eps = eps_raw * (1 - commute as i8);
    let k = k * (1 - commute);
    (eps, k)
}

/// Branch-free two-qubit commutator `[Q, P] / 2i`.
///
/// Each qubit is `[x, z]` bits; returns `(ε, x_out0, z_out0, x_out1, z_out1)`
/// with `ε ∈ {−1, 0, +1}` and the output flags masked to zero when the pair
/// commutes. Ported **verbatim** from `ppvm-pauli-sum/src/sum/rot2.rs::comm_2`,
/// including the pre-computed 16-entry `SIGN_NEG` mask.
#[inline(always)]
pub(crate) fn comm_2(q0: [u8; 2], q1: [u8; 2], p0: [u8; 2], p1: [u8; 2]) -> (i8, u8, u8, u8, u8) {
    let [x_a, z_a] = q0;
    let [x_b, z_b] = q1;
    let [x_c, z_c] = p0;
    let [x_d, z_d] = p1;

    // per-qubit anticommutation bits
    let a0 = (x_a & z_c) ^ (z_a & x_c);
    let a1 = (x_b & z_d) ^ (z_b & x_d);

    // the commutator is present when exactly one qubit anticommutes
    let present = a0 ^ a1;

    // 16-entry bit-mask: 1 → negative orientation, 0 → positive
    const SIGN_NEG: u16 = 0x2840;

    let idx0 = (z_a << 3) | (x_a << 2) | (z_c << 1) | x_c;
    let idx1 = (z_b << 3) | (x_b << 2) | (z_d << 1) | x_d;

    let neg0 = (((SIGN_NEG >> idx0) as u8) & 1) & a0;
    let neg1 = (((SIGN_NEG >> idx1) as u8) & 1) & a1;

    let coeff = (((1 - ((neg0 as i8) << 1)) * (a0 as i8))
        + ((1 - ((neg1 as i8) << 1)) * (a1 as i8)))
        * (present as i8);

    let x_out0 = (x_a ^ x_c) & present;
    let z_out0 = (z_a ^ z_c) & present;
    let x_out1 = (x_b ^ x_d) & present;
    let z_out1 = (z_b ^ z_d) & present;

    (coeff, x_out0, z_out0, x_out1, z_out1)
}

/// The generic single-qubit rotation branch for one stored term — old's
/// `rotate_1_map_insert_closure` (`ppvm-pauli-sum/src/sum/rot1.rs`).
///
/// Scales the diagonal `c` by `cos` and returns the anticommuting branch, or
/// `None` when the term commutes (or the site is lost). Shared by
/// [`RotationOne::rotate_1`]'s generic form **and** by [`RotationTwo`]'s
/// lost-qubit fallbacks, exactly as in old.
#[inline(always)]
fn rotate_1_branch<W, C>(
    k: &W,
    c: &mut C,
    axis: Pauli,
    qubit: usize,
    sin: &C,
    cos: &C,
) -> Option<(W, C)>
where
    W: PauliBits + Clone,
    C: Coefficient,
{
    // Check loss first — avoids reading x/z bits on a lost qubit (old's
    // `get_lbit` early-out, a const `false` for the non-lossy word).
    if k.is_lost(qubit) {
        return None;
    }
    let x = k.x_bit(qubit);
    let z = k.z_bit(qubit);
    let p_g = ((z as u8) << 1) | (x as u8);
    let (eps, p_q) = levi_civita(p_g, pauli_code(axis));
    if eps == 0 {
        return None;
    }
    let branch = c.clone() * sin.mul_sign(eps);
    *c *= cos.clone();
    let new_x = p_q & 0b01 != 0;
    let new_z = p_q & 0b10 != 0;
    Some((k.toggled_bits(qubit, new_x != x, new_z != z), branch))
}

/// Single-qubit rotation propagation on a Pauli-keyed `Sum`. Each rotation drives
/// the fused in-place [`Sum::rotate_in_place`] fast path: it scales every
/// diagonal coefficient by `cosθ` where it sits and merges only the anticommuting
/// `iGP` branch terms. The policy does **not** run: truncation is caller-driven
/// ([`Sum::truncate`]), as in old, so two sub-threshold branches on the same key
/// can still merge into a surviving term.
///
/// The per-axis commute test, flipped bits, and `±1` sign `ε` are ported
/// bit-for-bit from `ppvm-pauli-sum::sum::rot1` (`rx`/`ry`/`rz`). The `ε` column
/// is derived from the Pauli phase of `iGP` in
/// `lean/PPVM/Instantiations/Rotation.lean`
/// (`rx_eps_from_product`/`ry_eps_from_product`/`rz_eps_from_product`, real by
/// `branchExp_isRealPhase`).
impl<S, P, W, C, Ang> RotationOne<C, Ang> for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + RotateInPlace<W, C> + Retain<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient,
    Ang: Angle<C>,
    P: crate::policy::Policy<W, C>,
{
    /// Rotate about `axis` on `qubit` by `theta`, dispatching to the per-axis
    /// fast paths below.
    ///
    /// This is the old crate's axis-generic entry point
    /// (`ppvm-pauli-sum::sum::rot1::rotate_1`), which drove a `levi_civita`
    /// table lookup per term; dispatching once on the axis instead is the same
    /// function — old `rx`/`ry`/`rz` were themselves specialized overrides of
    /// `rotate_1` and are diffed against it. `Pauli::I` commutes with every
    /// term, so `levi_civita(p, I) = (0, _)` made the old pass return `None`
    /// everywhere and mutate nothing (`map_insert` does not truncate); the
    /// no-op below is that same behaviour.
    #[inline]
    fn rotate_1(&mut self, axis: Pauli, qubit: usize, theta: Ang) {
        match axis {
            Pauli::X => self.rx(qubit, theta),
            Pauli::Y => self.ry(qubit, theta),
            Pauli::Z => self.rz(qubit, theta),
            Pauli::I => {}
        }
    }

    /// Rotate about `X` on `qubit` by `theta`. A term commutes iff its `z` bit at
    /// `qubit` is clear (`I`/`X`); an anticommuting `Z`/`Y` (flip `x`) branches to
    /// `Y`/`Z` with `ε = −1` if `x` else `+1`.
    #[inline(always)]
    fn rx(&mut self, qubit: usize, theta: Ang) {
        let (sin, cos) = theta.sin_cos();
        self.rotate_x(qubit, sin, cos);
    }

    /// Rotate about `Y` on `qubit` by `theta`. A term commutes iff its `x` and `z`
    /// bits at `qubit` agree (`I`/`Y`); an anticommuting `X`/`Z` (flip `x` and `z`)
    /// branches to `Z`/`X` with `ε = −1` if `z` else `+1`.
    #[inline]
    fn ry(&mut self, qubit: usize, theta: Ang) {
        let (sin, cos) = theta.sin_cos();
        self.rotate_in_place(move |k: &W, c: &mut C| {
            if k.is_lost(qubit) {
                return None;
            }
            let x = k.x_bit(qubit);
            let z = k.z_bit(qubit);
            if x == z {
                return None;
            }
            let branch = c.clone() * sin.mul_sign(if z { -1 } else { 1 });
            *c *= cos.clone();
            // ry anticommuting branch flips both X and Z (X/Z ↦ Z/X).
            Some((k.toggled_bits(qubit, true, true), branch))
        });
    }

    /// Rotate about `Z` on `qubit` by `theta`. A term commutes iff its `x` bit at
    /// `qubit` is clear (`I`/`Z`); an anticommuting `X`/`Y` (flip `z`) branches to
    /// `Y`/`X` with `ε = +1` if `z` else `−1`.
    #[inline]
    fn rz(&mut self, qubit: usize, theta: Ang) {
        let (sin, cos) = theta.sin_cos();
        self.rotate_in_place(move |k: &W, c: &mut C| {
            if k.is_lost(qubit) || !k.x_bit(qubit) {
                return None;
            }
            let z = k.z_bit(qubit);
            let branch = c.clone() * sin.mul_sign(if z { 1 } else { -1 });
            *c *= cos.clone();
            // rz anticommuting branch flips the Z bit (X/Y ↦ Y/X).
            Some((k.toggled_bits(qubit, false, true), branch))
        });
    }
}

/// Two-qubit Pauli-rotation propagation on a Pauli-keyed `Sum`: the generic
/// [`rotate_2`](RotationTwo::rotate_2) over [`comm_2`] plus the three
/// hand-written diagonal fast paths `rzz`/`rxx`/`ryy`.
///
/// Ported from `ppvm-pauli-sum/src/sum/rot2.rs`, which is architecture feature 5's
/// most load-bearing half: each fast path is **one** traversal that computes
/// commutation from two bits, scales the survivor by `cos` and emits exactly one
/// branch, where the `cnot; rz; cnot` decomposition is three traversals plus two
/// full re-keys. Both forms drive the fused in-place
/// [`Sum::rotate_in_place`](crate::Sum) path, so — as everywhere else in this
/// crate — nothing truncates and nothing reduces: an exactly-zero branch is
/// merged and an exact cancellation stays in the support.
///
/// The `ε` conventions are old's verbatim, including the sign asymmetry between
/// the generic path (`sin.mul_sign(-eps)`, `ε` from [`comm_2`]) and the fast
/// paths (`sin.mul_sign(eps)`, `ε` read off the anticommuting qubit's own bits);
/// old pins the two against each other in `rot2.rs`'s `rxx_matches_generic` /
/// `ryy_matches_generic` / `rzz_matches_generic`, which this crate's tests
/// reproduce — but only at the three **diagonal** axes. The off-diagonal axis
/// pairs are covered by `lean/PPVM/Instantiations/Rotation.lean`, which
/// transcribes [`comm_2`] line for line (`comm2Coeff`, `comm2Key`, the
/// `SIGN_NEG = 0x2840` mask as `signNegMask`/`signNegIdx`) and proves over all
/// `2⁸` (axis, key) bit patterns that the early-out is the `anti2` branch
/// predicate (`comm2Coeff_eq_zero_iff`), that the sign this path applies (`−ε`)
/// is the real branch prefactor of `−i·[G_a ⊗ G_b, P]/2`
/// (`comm2_generic_sign_eq_branchExp2`), and that the masked output bits are the
/// two-site product key (`comm2_key_eq_mulBits2`).
impl<S, P, W, C, Ang> RotationTwo<C, Ang> for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + RotateInPlace<W, C> + Retain<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient,
    Ang: Angle<C>,
    P: crate::policy::Policy<W, C>,
{
    /// Rotate about `P_a ⊗ P_b` by `theta`, each axis given as `[x, z]` bits.
    ///
    /// # Panics
    ///
    /// If any axis component exceeds `1`. Old's guard is `> 3` with the message
    /// "Rotation axis cannot be L" — but the components are single **bits**, so
    /// old lets `2` and `3` through and then indexes its 4-element `PAULIS` table
    /// out of bounds (an index panic instead of the intended message). The bound
    /// is `> 1` here: a documented, deliberate divergence on a *validation*
    /// defect where old panics either way (suspected old bug 5).
    ///
    /// A lost qubit degrades to the single-qubit rotation on the survivor, with
    /// the axis read out of [`AXIS_PAULIS`] — old's fallback, and dead code for
    /// the non-lossy word (`is_lost` is a const `false`).
    fn rotate_2(&mut self, axis_a: [u8; 2], axis_b: [u8; 2], a: usize, b: usize, theta: Ang) {
        let [axis_a_x, axis_a_z] = axis_a;
        let [axis_b_x, axis_b_z] = axis_b;
        assert!(
            axis_a_x <= 1 && axis_a_z <= 1 && axis_b_x <= 1 && axis_b_z <= 1,
            "rotation axis components must be 0 or 1 ([x, z] bits)"
        );
        let (sin, cos) = theta.sin_cos();
        let pauli_a = AXIS_PAULIS[((axis_a_z << 1) | axis_a_x) as usize];
        let pauli_b = AXIS_PAULIS[((axis_b_z << 1) | axis_b_x) as usize];
        self.rotate_in_place(move |k: &W, c: &mut C| {
            // Both-lost is handled by the single-qubit logic (which early-outs on
            // the surviving-side loss check), exactly as in old.
            if k.is_lost(a) {
                return rotate_1_branch(k, c, pauli_b, b, &sin, &cos);
            }
            if k.is_lost(b) {
                return rotate_1_branch(k, c, pauli_a, a, &sin, &cos);
            }
            let (x_a, z_a) = (k.x_bit(a), k.z_bit(a));
            let (x_b, z_b) = (k.x_bit(b), k.z_bit(b));
            let (eps, nx_a, nz_a, nx_b, nz_b) = comm_2(
                axis_a,
                axis_b,
                [x_a as u8, z_a as u8],
                [x_b as u8, z_b as u8],
            );
            if eps == 0 {
                return None;
            }
            let branch = c.clone() * sin.mul_sign(-eps);
            *c *= cos.clone();
            // ONE key build for the two-site branch — old wrote its four bits into
            // a single `k.clone()`; chaining `toggled_bits` would copy both planes
            // and rebuild the word twice per produced term.
            let new_key = k.toggled_bits2(
                a,
                (nx_a == 1) != x_a,
                (nz_a == 1) != z_a,
                b,
                (nx_b == 1) != x_b,
                (nz_b == 1) != z_b,
            );
            Some((new_key, branch))
        });
    }

    /// `exp(−i·θ/2·Z_a Z_b)` — the single-pass fast path.
    ///
    /// `ZZ` commutes with a term iff the two qubits **agree** on carrying an
    /// X-component; the anticommuting qubit is the one with its x bit set, and the
    /// sign is `+1` when that qubit is `Y` (z set), `−1` when it is `X`. The branch
    /// toggles both z bits.
    #[inline(always)]
    fn rzz(&mut self, a: usize, b: usize, theta: Ang) {
        let (sin, cos) = theta.sin_cos();
        self.rotate_zz(a, b, sin, cos);
    }

    /// `exp(−i·θ/2·X_a X_b)` — the single-pass fast path.
    ///
    /// `X` anticommutes with a qubit's Pauli iff that Pauli carries a Z-component
    /// (`Z` or `Y`), so `XX` commutes iff the two qubits agree on it; the sign is
    /// `+1` when the anticommuting qubit is `Z`, `−1` when it is `Y`. The branch
    /// toggles both x bits.
    #[inline]
    fn rxx(&mut self, a: usize, b: usize, theta: Ang) {
        let (sin, cos) = theta.sin_cos();
        self.rotate_in_place(move |k: &W, c: &mut C| {
            if k.is_lost(a) {
                return rotate_1_branch(k, c, Pauli::X, b, &sin, &cos);
            }
            if k.is_lost(b) {
                return rotate_1_branch(k, c, Pauli::X, a, &sin, &cos);
            }
            let za = k.z_bit(a);
            let zb = k.z_bit(b);
            if za == zb {
                return None;
            }
            let xa = k.x_bit(a);
            let xb = k.x_bit(b);
            let x_anti = if za { xa } else { xb };
            let eps: i8 = if x_anti { -1 } else { 1 };
            let branch = c.clone() * sin.mul_sign(eps);
            *c *= cos.clone();
            // One plane copy for both toggled x bits (see `rotate_2`).
            Some((k.toggled_bits2(a, true, false, b, true, false), branch))
        });
    }

    /// `exp(−i·θ/2·Y_a Y_b)` — the single-pass fast path.
    ///
    /// `Y` anticommutes with a qubit's Pauli iff it is `X` or `Z` (x ≠ z), so `YY`
    /// commutes iff the two qubits agree on that; the sign is `+1` when the
    /// anticommuting qubit is `X`, `−1` when it is `Z`. The branch toggles all four
    /// bits.
    #[inline]
    fn ryy(&mut self, a: usize, b: usize, theta: Ang) {
        let (sin, cos) = theta.sin_cos();
        self.rotate_in_place(move |k: &W, c: &mut C| {
            if k.is_lost(a) {
                return rotate_1_branch(k, c, Pauli::Y, b, &sin, &cos);
            }
            if k.is_lost(b) {
                return rotate_1_branch(k, c, Pauli::Y, a, &sin, &cos);
            }
            let xa = k.x_bit(a);
            let za = k.z_bit(a);
            let xb = k.x_bit(b);
            let zb = k.z_bit(b);
            let pa = xa ^ za;
            let pb = xb ^ zb;
            if pa == pb {
                return None;
            }
            let x_anti = if pa { xa } else { xb };
            let eps: i8 = if x_anti { 1 } else { -1 };
            let branch = c.clone() * sin.mul_sign(eps);
            *c *= cos.clone();
            // One plane copy for all four toggled bits (see `rotate_2`).
            Some((k.toggled_bits2(a, true, true, b, true, true), branch))
        });
    }
}

/// In-plane rotation `R(axis_angle, θ)` on a Pauli-keyed `Sum`.
///
/// **The sub-rotation order is Heisenberg (backward)**: `rz(axis_angle)`,
/// `rx(θ)`, `rz(−axis_angle)` — the reverse of the tableau's forward order,
/// because a `Sum` propagates observables backward. The observable consequence,
/// which old pins in `rot1.rs::test_r` and which a forward-ordered impl gets
/// wrong, is that `r(q, π/2, θ) == ry(q, θ)` and **not** `ry(q, −θ)`. Proved in
/// `lean/PPVM/Instantiations/Rotation.lean`: `rotXY_heisenberg_order` shows the
/// composite `M_z(−φ) ∘ M_x(θ) ∘ M_z(φ)` is Rodrigues rotation about
/// `cos φ·X + sin φ·Y` by `θ`, with `rotXY_zero_eq_rx` and `rotXY_halfPi_eq_ry`
/// the two behavioural corollaries.
///
/// Ported from `ppvm-pauli-sum/src/sum/rot1.rs`'s `impl RotXY`.
impl<S, P, W, C, Ang> RotXY<C, Ang> for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + RotateInPlace<W, C> + Retain<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient,
    Ang: Angle<C> + Clone + std::ops::Neg<Output = Ang>,
    P: crate::policy::Policy<W, C>,
{
    #[inline]
    fn r(&mut self, qubit: usize, axis_angle: Ang, theta: Ang) {
        self.rz(qubit, axis_angle.clone());
        self.rx(qubit, theta);
        self.rz(qubit, -axis_angle);
    }
}
