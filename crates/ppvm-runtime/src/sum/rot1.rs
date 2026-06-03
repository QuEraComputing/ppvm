// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::traits::*;
use crate::{char::Pauli, config::Config, sum::PauliSum};

pub(crate) fn rotate_1_map_insert_closure<T: Config>(
    k: &T::PauliWordType,
    v: &mut T::Coeff,
    axis: Pauli,
    addr0: usize,
    sin: &T::Coeff,
    cos: &T::Coeff,
) -> Option<(T::PauliWordType, T::Coeff)> {
    if axis == Pauli::L {
        panic!("Rotation axis cannot be L");
    }
    // Check loss first — avoids reading xbit/zbit on lost qubits
    if k.get_lbit(addr0) {
        return None;
    }
    let xbit = k.get_xbit(addr0);
    let zbit = k.get_zbit(addr0);
    // Reconstruct 2-bit Pauli code (I=00, X=01, Z=10, Y=11) from bits
    let p_g = (zbit as u8) << 1 | (xbit as u8);
    let (eps, p_q) = levi_civita(p_g, axis as u8);
    if eps == 0 {
        None
    } else {
        let mut coeff = v.clone();
        *v *= cos.clone();
        let mut new_word = k.clone();
        new_word.set_xbit(addr0, p_q & 0b01 != 0);
        new_word.set_zbit(addr0, p_q & 0b10 != 0);
        new_word.rehash();

        coeff *= sin.mul_sign(eps);
        Some((new_word, coeff))
    }
}

impl<T: Config> RotationOne<T> for PauliSum<T>
where
    T::Coeff: std::ops::MulAssign,
    T::Map: ACMapInsert<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType> + ACMapConsume,
{
    fn rotate_1(&mut self, axis: Pauli, addr0: usize, theta: <T as Config>::Coeff) {
        let (sin, cos) = theta.sin_cos();
        self.map_insert(|k, v| rotate_1_map_insert_closure::<T>(k, v, axis, addr0, &sin, &cos));
    }

    #[inline]
    fn rx(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        // Axis = X (xbit=1, zbit=0). Commutes when zbit==false (I or X).
        // Anticommuting: Z(xbit=0)→new Y, eps=+1; Y(xbit=1)→new Z, eps=-1
        let (sin, cos) = theta.into().sin_cos();
        self.map_insert(|k, v| {
            if k.get_lbit(addr0) {
                return None;
            }
            let zbit = k.get_zbit(addr0);
            if !zbit {
                return None;
            }
            let xbit = k.get_xbit(addr0);
            let mut coeff = v.clone();
            *v *= cos.clone();
            let mut new_word = k.clone();
            new_word.set_xbit(addr0, !xbit);
            new_word.rehash();
            let eps: i8 = if xbit { -1 } else { 1 };
            coeff *= sin.mul_sign(eps);
            Some((new_word, coeff))
        });
    }

    #[inline]
    fn ry(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        // Axis = Y (xbit=1, zbit=1). Commutes when xbit==zbit (I or Y).
        // Anticommuting: X(zbit=0)→new Z, eps=+1; Z(zbit=1)→new X, eps=-1
        let (sin, cos) = theta.into().sin_cos();
        self.map_insert(|k, v| {
            if k.get_lbit(addr0) {
                return None;
            }
            let xbit = k.get_xbit(addr0);
            let zbit = k.get_zbit(addr0);
            if xbit == zbit {
                return None;
            }
            let mut coeff = v.clone();
            *v *= cos.clone();
            let mut new_word = k.clone();
            new_word.set_xbit(addr0, !xbit);
            new_word.set_zbit(addr0, !zbit);
            new_word.rehash();
            let eps: i8 = if zbit { -1 } else { 1 };
            coeff *= sin.mul_sign(eps);
            Some((new_word, coeff))
        });
    }

    #[inline]
    fn rz(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        // Axis = Z (xbit=0, zbit=1). Commutes when xbit==false (I or Z).
        // Anticommuting: X(zbit=0)→new Y, eps=-1; Y(zbit=1)→new X, eps=+1
        let (sin, cos) = theta.into().sin_cos();
        self.map_insert(|k, v| {
            if k.get_lbit(addr0) {
                return None;
            }
            let xbit = k.get_xbit(addr0);
            if !xbit {
                return None;
            }
            let zbit = k.get_zbit(addr0);
            let mut coeff = v.clone();
            *v *= cos.clone();
            let mut new_word = k.clone();
            new_word.set_zbit(addr0, !zbit);
            new_word.rehash();
            let eps: i8 = if zbit { 1 } else { -1 };
            coeff *= sin.mul_sign(eps);
            Some((new_word, coeff))
        });
    }
}

impl<T: Config> RotXY<T> for PauliSum<T>
where
    PauliSum<T>: RotationOne<T>,
{
    fn r(&mut self, addr0: usize, axis_angle: <T as Config>::Coeff, theta: <T as Config>::Coeff) {
        // R(axis_angle, θ) = RZ(axis_angle)·RX(θ)·RZ(−axis_angle). PauliSum runs
        // in the Heisenberg picture (observables propagate backward), so the
        // sub-rotations are emitted in reverse of the tableau's forward order.
        self.rz(addr0, axis_angle.clone());
        self.rx(addr0, theta);
        self.rz(addr0, -axis_angle);
    }
}

/// 2-bit Pauli code: 00 I, 01 X, 10 Z, 11 Y
/// Returns \(ε, k\) so that  –i \[P_i, P_j\]/2 = ε · P_k.
/// For every commuting pair it yields (0, 0).
#[inline]
pub fn levi_civita(i: u8, j: u8) -> (i8, u8) {
    // --------------------------------------------------- third Pauli by XOR
    let k = i ^ j; // 0 when i == j

    // ----------- commute ⇔ i==0  OR  j==0  OR  k==0  (no false positives)
    let commute = ((i == 0) | (j == 0) | (k == 0)) as u8; // 1 = commute

    // ------------------------------------------------------ sign ε_{ijk}
    #[inline]
    fn rank(p: u8) -> u8 {
        let b1 = p >> 1; // MSB
        (b1 << 1).wrapping_sub(b1 & (p & 1)) // 0,1,2 for X,Y,Z
    }

    let ri = rank(i);
    let rj = rank(j);

    // diff = (rj - ri) mod 3   without an actual modulus
    let mut diff = rj.wrapping_sub(ri).wrapping_add(3); // 0…5
    diff -= 3 & (0u8.wrapping_sub(diff >> 2)); // if ≥3 subtract 3

    // +1 when diff == 1,  –1 when diff == 2
    let eps_raw = 1i8 - 2 * ((diff >> 1) as i8);

    // --------------------------------------------------- zero when commute
    let eps = eps_raw * (1 - commute as i8); // 0 if commute
    let k = k * (1 - commute); // 0 if commute

    (eps, k)
}

#[cfg(test)]
mod tests {
    use crate::config::fxhash::ByteF64;

    use super::*;

    #[test]
    fn test_rx() {
        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("X", 1.0);
        answer.rx(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("X", 1.0);
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Y", 1.0);
        answer.rx(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Y", 2.1_f64.cos());
        expect += ("Z", -2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Z", 1.0);
        answer.rx(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Z", 2.1_f64.cos());
        expect += ("Y", 2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("I", 1.0);
        answer.rx(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("I", 1.0);
        assert_eq!(answer, expect);
    }

    #[test]
    fn test_ry() {
        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("X", 1.0);
        answer.ry(0, 2.1);

        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("X", 2.1_f64.cos());
        expect += ("Z", 2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Y", 1.0);
        answer.ry(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Y", 1.0);
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Z", 1.0);
        answer.ry(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Z", 2.1_f64.cos());
        expect += ("X", -2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("I", 1.0);
        answer.ry(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("I", 1.0);
        assert_eq!(answer, expect);
    }

    #[test]
    fn test_r() {
        use std::f64::consts::FRAC_PI_2;
        let theta = 2.1;

        // r(axis_angle=0, θ) == rx(θ).
        let mut via_r: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        via_r += ("Z", 1.0);
        via_r.r(0, 0.0, theta);
        let mut via_rx: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        via_rx += ("Z", 1.0);
        via_rx.rx(0, theta);
        assert!((via_r.overlap(&via_rx) - 1.0).abs() < 1e-9);

        // r(axis_angle=π/2, θ) must equal ry(θ) — NOT ry(−θ). This is the case
        // that distinguishes the Heisenberg (backward) order from the
        // Schrödinger one: a forward-ordered impl would yield ry(−θ) here.
        let mut via_r: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        via_r += ("Z", 1.0);
        via_r.r(0, FRAC_PI_2, theta);
        let mut via_ry: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        via_ry += ("Z", 1.0);
        via_ry.ry(0, theta);
        assert!((via_r.overlap(&via_ry) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_rz() {
        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("X", 1.0);
        answer.rz(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("X", 2.1_f64.cos());
        expect += ("Y", -2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Y", 1.0);
        answer.rz(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Y", 2.1_f64.cos());
        expect += ("X", 2.1_f64.sin());
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("Z", 1.0);
        answer.rz(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("Z", 1.0);
        assert_eq!(answer, expect);

        let mut answer: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        answer += ("I", 1.0);
        answer.rz(0, 2.1);
        let mut expect: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(1).build();
        expect += ("I", 1.0);
        assert_eq!(answer, expect);
    }
}
