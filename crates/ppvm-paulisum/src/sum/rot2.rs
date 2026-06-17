// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::rot1::rotate_1_map_insert_closure;
use crate::sum::PauliSum;
use ppvm_runtime::char::Pauli;
use ppvm_runtime::config::Config;
use ppvm_runtime::traits::*;

const PAULIS: [Pauli; 4] = [Pauli::I, Pauli::X, Pauli::Z, Pauli::Y];

impl<T: Config> RotationTwo<T> for PauliSum<T>
where
    T::Coeff: std::ops::MulAssign,
    T::Map: ACMapInsert<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType> + ACMapConsume,
{
    fn rotate_2(&mut self, axis_a: [u8; 2], axis_b: [u8; 2], a: usize, b: usize, theta: T::Coeff) {
        let [axis_a_x, axis_a_z] = axis_a;
        let [axis_b_x, axis_b_z] = axis_b;
        if axis_a_x > 3 || axis_a_z > 3 || axis_b_x > 3 || axis_b_z > 3 {
            panic!("Rotation axis cannot be L");
        }
        let (sin, cos) = theta.sin_cos();
        let pauli_a = PAULIS[(axis_a_z << 1 | axis_a_x) as usize];
        let pauli_b = PAULIS[(axis_b_z << 1 | axis_b_x) as usize];
        self.map_insert(|k, v| {
            // NOTE: case of both qubits being lost is handled by single-qubit rotation logic
            if k.get_lbit(a) {
                // fall back to single-qubit rotation on qubit b
                return rotate_1_map_insert_closure::<T>(k, v, pauli_b, b, &sin, &cos);
            }
            if k.get_lbit(b) {
                // fall back to single-qubit rotation on qubit a
                return rotate_1_map_insert_closure::<T>(k, v, pauli_a, a, &sin, &cos);
            }
            let (eps, x_a, z_a, x_b, z_b) = comm_2(
                axis_a,
                axis_b,
                [k.get_xbit(a) as u8, k.get_zbit(a) as u8],
                [k.get_xbit(b) as u8, k.get_zbit(b) as u8],
            );

            if eps == 0 {
                None
            } else {
                let mut coeff = v.clone();
                *v *= cos.clone();

                let mut new_word = k.clone();
                new_word.set_xbit(a, x_a == 1);
                new_word.set_xbit(b, x_b == 1);
                new_word.set_zbit(a, z_a == 1);
                new_word.set_zbit(b, z_b == 1);
                new_word.rehash();

                coeff *= sin.mul_sign(-eps);
                Some((new_word, coeff))
            }
        });
    }

    #[inline]
    fn rzz(&mut self, a: usize, b: usize, theta: impl Into<T::Coeff>) {
        let (sin, cos) = theta.into().sin_cos();
        self.map_insert(|k, v| {
            // Loss fallbacks — identical to the generic `rotate_2` path
            // (axis on the surviving qubit is Z for a ZZ rotation). These
            // branches are dead-code-eliminated for non-lossy PauliWord,
            // since `get_lbit` is a const `false`.
            if k.get_lbit(a) {
                return rotate_1_map_insert_closure::<T>(k, v, Pauli::Z, b, &sin, &cos);
            }
            if k.get_lbit(b) {
                return rotate_1_map_insert_closure::<T>(k, v, Pauli::Z, a, &sin, &cos);
            }
            let xa = k.get_xbit(a);
            let xb = k.get_xbit(b);
            // ZZ commutes iff both qubits agree on having an X-component.
            if xa == xb {
                return None;
            }
            let za = k.get_zbit(a);
            let zb = k.get_zbit(b);
            // The anticommuting qubit is the one with xbit == 1.
            // sign = +1 if it is Y (zbit set), -1 if it is X (zbit clear).
            let z_anti = if xa { za } else { zb };
            let eps: i8 = if z_anti { 1 } else { -1 };

            let mut coeff = v.clone();
            *v *= cos.clone();

            let mut new_word = k.clone();
            new_word.set_zbit(a, !za);
            new_word.set_zbit(b, !zb);
            new_word.rehash();

            coeff *= sin.mul_sign(eps);
            Some((new_word, coeff))
        });
    }

    #[inline]
    fn rxx(&mut self, a: usize, b: usize, theta: impl Into<T::Coeff>) {
        let (sin, cos) = theta.into().sin_cos();
        self.map_insert(|k, v| {
            // Loss fallback: a lost qubit leaves a single-qubit X rotation on
            // the surviving partner (axis on the surviving qubit is X for an
            // XX rotation), matching the generic `rotate_2` path.
            if k.get_lbit(a) {
                return rotate_1_map_insert_closure::<T>(k, v, Pauli::X, b, &sin, &cos);
            }
            if k.get_lbit(b) {
                return rotate_1_map_insert_closure::<T>(k, v, Pauli::X, a, &sin, &cos);
            }
            let za = k.get_zbit(a);
            let zb = k.get_zbit(b);
            // X anticommutes with a qubit's Pauli iff that Pauli carries a
            // Z-component (Z or Y). XX commutes iff both qubits agree on that.
            if za == zb {
                return None;
            }
            let xa = k.get_xbit(a);
            let xb = k.get_xbit(b);
            // The anticommuting qubit is the one carrying Z (its z-bit set).
            // sign = +1 if it is Z (x-bit clear), -1 if it is Y (x-bit set).
            let x_anti = if za { xa } else { xb };
            let eps: i8 = if x_anti { -1 } else { 1 };

            let mut coeff = v.clone();
            *v *= cos.clone();

            let mut new_word = k.clone();
            new_word.set_xbit(a, !xa);
            new_word.set_xbit(b, !xb);
            new_word.rehash();

            coeff *= sin.mul_sign(eps);
            Some((new_word, coeff))
        });
    }

    #[inline]
    fn ryy(&mut self, a: usize, b: usize, theta: impl Into<T::Coeff>) {
        let (sin, cos) = theta.into().sin_cos();
        self.map_insert(|k, v| {
            // Loss fallback: a lost qubit leaves a single-qubit Y rotation on
            // the surviving partner (axis on the surviving qubit is Y for a
            // YY rotation), matching the generic `rotate_2` path.
            if k.get_lbit(a) {
                return rotate_1_map_insert_closure::<T>(k, v, Pauli::Y, b, &sin, &cos);
            }
            if k.get_lbit(b) {
                return rotate_1_map_insert_closure::<T>(k, v, Pauli::Y, a, &sin, &cos);
            }
            let xa = k.get_xbit(a);
            let za = k.get_zbit(a);
            let xb = k.get_xbit(b);
            let zb = k.get_zbit(b);
            // Y anticommutes with a qubit's Pauli iff it is X or Z, i.e. iff
            // its x-bit and z-bit differ. YY commutes iff both qubits agree.
            let pa = xa ^ za;
            let pb = xb ^ zb;
            if pa == pb {
                return None;
            }
            // The anticommuting qubit is X or Z (x-bit != z-bit).
            // sign = +1 if it is X (x-bit set), -1 if it is Z (x-bit clear).
            let x_anti = if pa { xa } else { xb };
            let eps: i8 = if x_anti { 1 } else { -1 };

            let mut coeff = v.clone();
            *v *= cos.clone();

            let mut new_word = k.clone();
            new_word.set_xbit(a, !xa);
            new_word.set_zbit(a, !za);
            new_word.set_xbit(b, !xb);
            new_word.set_zbit(b, !zb);
            new_word.rehash();

            coeff *= sin.mul_sign(eps);
            Some((new_word, coeff))
        });
    }
}

// R_{G}[P] = cos(theta) P - sin(theta) [G, P]/2i
// Encoding for one qubit
//   x = 0/1  (X flag)   z = 0/1  (Z flag)
//   (x,z) = 00 → I , 10 → X , 01 → Z , 11 → Y   (little-endian)
//
// Two–qubit word Q = q0 ⊗ q1  is given as four bits
//   (x_a , z_a)  … qubit-0 of Q
//   (x_b , z_b)  … qubit-1 of Q
// P = p0 ⊗ p1   likewise
//   (x_c , z_c)  … qubit-0 of P
//   (x_d , z_d)  … qubit-1 of P
//
// The function returns
//   (coeff , xx_out , zz_out)
// where  coeff ∈ {-1,0,+1}  and  (xx_out, zz_out) are the usual packed masks
// (bit-0 = qubit-0, bit-1 = qubit-1).
/// Branch-free commutator  [ Q , P ] / (2 i )
///
///   Each qubit is encoded as `[x, z]` bits (0 or 1):
///       `q0 = [x_a, z_a]`  Q qubit-0       `q1 = [x_b, z_b]`  Q qubit-1
///       `p0 = [x_c, z_c]`  P qubit-0       `p1 = [x_d, z_d]`  P qubit-1
///
///   Returns (coeff , x_out0 , z_out0 , x_out1 , z_out1)
///           coeff ∈ { -1 , 0 , +1 }
#[inline(always)]
pub fn comm_2(q0: [u8; 2], q1: [u8; 2], p0: [u8; 2], p1: [u8; 2]) -> (i8, u8, u8, u8, u8) {
    let [x_a, z_a] = q0;
    let [x_b, z_b] = q1;
    let [x_c, z_c] = p0;
    let [x_d, z_d] = p1;
    // ── 1.  per-qubit anticommutation bits  a₀ , a₁  ───────────────────────
    let a0 = (x_a & z_c) ^ (z_a & x_c); // qubit-0 anticommutes?
    let a1 = (x_b & z_d) ^ (z_b & x_d); // qubit-1 anticommutes?

    // overall commutator present when exactly one qubit anticommutes
    let present = a0 ^ a1; // 0 ↔ commute, 1 ↔ anticommute

    // ── 2.  sign of the coefficient (+1 / −1)  ─────────────────────────────
    // 16-entry bit-mask: 1 → “negative orientation” ( −1 ), 0 → “positive” (+1)
    const SIGN_NEG: u16 = 0x2840; // pre-computed once

    let idx0 = (z_a << 3) | (x_a << 2) | (z_c << 1) | x_c; // qubit-0 pair index 0-15
    let idx1 = (z_b << 3) | (x_b << 2) | (z_d << 1) | x_d; // qubit-1 pair index

    let neg0 = (((SIGN_NEG >> idx0) as u8) & 1) & a0; // only meaningful if a₀ = 1
    let neg1 = (((SIGN_NEG >> idx1) as u8) & 1) & a1;

    // coeff = ( +1 or −1 ) from the unique acted-on qubit; 0 if they commute
    let coeff = (((1 - ((neg0 as i8) << 1)) * (a0 as i8))
        + ((1 - ((neg1 as i8) << 1)) * (a1 as i8)))
        * (present as i8); // ensures 0 when commuting

    // ── 3.  output flags  (product of Q and P)  masked to zero if commuting ─
    let x_out0 = (x_a ^ x_c) & present;
    let z_out0 = (z_a ^ z_c) & present;
    let x_out1 = (x_b ^ x_d) & present;
    let z_out1 = (z_b ^ z_d) & present;

    (coeff, x_out0, z_out0, x_out1, z_out1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::fxhash::ByteF64;

    type C = ByteF64<2>;
    const PAULI_CHARS: [char; 4] = ['I', 'X', 'Y', 'Z'];

    fn sum_with(word: &str) -> PauliSum<C> {
        let mut s: PauliSum<C> = PauliSum::builder().n_qubits(2).build();
        s += (word, 1.0);
        s
    }

    /// Every diagonal two-qubit rotation has a hand-written fast path in this
    /// file that bypasses `comm_2`. Assert it produces exactly the same
    /// `PauliSum` as the generic `rotate_2` (the pre-existing `comm_2`-based
    /// path) for every two-qubit input Pauli, at a representative angle.
    fn assert_matches_generic(
        axis: [u8; 2],
        special: impl Fn(&mut PauliSum<C>, usize, usize, f64),
    ) {
        let theta = 0.7_f64;
        for &p in &PAULI_CHARS {
            for &q in &PAULI_CHARS {
                let word = format!("{p}{q}");

                let mut got = sum_with(&word);
                special(&mut got, 0, 1, theta);

                let mut want = sum_with(&word);
                want.rotate_2(axis, axis, 0, 1, theta);

                assert_eq!(got, want, "fast path disagrees with rotate_2 for {word}");
            }
        }
    }

    #[test]
    fn rxx_matches_generic() {
        assert_matches_generic([1, 0], |s, a, b, t| s.rxx(a, b, t));
    }

    #[test]
    fn ryy_matches_generic() {
        assert_matches_generic([1, 1], |s, a, b, t| s.ryy(a, b, t));
    }

    #[test]
    fn rzz_matches_generic() {
        assert_matches_generic([0, 1], |s, a, b, t| s.rzz(a, b, t));
    }

    /// Explicit hand-computed values, independent of `rotate_2`/`comm_2`, so a
    /// bug shared by both the fast path and the generic path can't hide.
    #[test]
    fn rzz_explicit_values() {
        let t = 0.7_f64;

        // X_0 I_1 --rzz--> cos·XI − sin·YZ  (anticommutes; X-carrier → −1)
        let mut s = sum_with("XI");
        s.rzz(0, 1, t);
        let mut want: PauliSum<C> = PauliSum::builder().n_qubits(2).build();
        want += ("XI", t.cos());
        want += ("YZ", -t.sin());
        assert_eq!(s, want);

        // Y_0 I_1 --rzz--> cos·YI + sin·XZ  (anticommutes; Y-carrier → +1)
        let mut s = sum_with("YI");
        s.rzz(0, 1, t);
        let mut want: PauliSum<C> = PauliSum::builder().n_qubits(2).build();
        want += ("YI", t.cos());
        want += ("XZ", t.sin());
        assert_eq!(s, want);

        // ZZ commutes with the ZZ generator → unchanged.
        let mut s = sum_with("ZZ");
        s.rzz(0, 1, t);
        assert_eq!(s, sum_with("ZZ"));
    }

    /// A diagonal rotation on non-adjacent qubits must address the right
    /// slots (exercises `a`/`b` other than `0,1`).
    #[test]
    fn rzz_non_adjacent_addressing() {
        let theta = 0.7_f64;
        for &p in &PAULI_CHARS {
            for &q in &PAULI_CHARS {
                // 3-qubit word with the acted-on Paulis on qubits 0 and 2.
                let word = format!("{p}I{q}");

                let mut got: PauliSum<C> = PauliSum::builder().n_qubits(3).build();
                got += (word.as_str(), 1.0);
                got.rzz(0, 2, theta);

                let mut want: PauliSum<C> = PauliSum::builder().n_qubits(3).build();
                want += (word.as_str(), 1.0);
                want.rotate_2([0, 1], [0, 1], 0, 2, theta);

                assert_eq!(got, want, "rzz(0,2) disagrees with rotate_2 for {word}");
            }
        }
    }
}
