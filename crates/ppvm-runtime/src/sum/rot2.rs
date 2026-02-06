use crate::traits::*;
use crate::{config::Config, sum::PauliSum};

impl<T, S, H> RotationTwo<T> for PauliSum<T>
where
    S: PauliStorage,
    H: std::hash::BuildHasher + Clone + Default,
    T: Config<Storage = S, BuildHasher = H>,
    T::Coeff: std::ops::MulAssign,
    T::Map: ACMapInsert<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType> + ACMapConsume,
{
    fn rotate_2(
        &mut self,
        axis_a_x: u8,
        axis_a_z: u8,
        axis_b_x: u8,
        axis_b_z: u8,
        a: usize,
        b: usize,
        theta: T::Coeff,
    ) {
        let (sin, cos) = theta.sin_cos();
        self.map_insert(|k, v| {
            if k.get_lbit(a) || k.get_lbit(b) {
                // Two-qubit rotations are undefined on loss; keep term unchanged.
                return None;
            }
            let (eps, x_a, z_a, x_b, z_b) = comm_2(
                axis_a_x,
                axis_a_z,
                axis_b_x,
                axis_b_z,
                k.get_xbit(a) as u8,
                k.get_zbit(a) as u8,
                k.get_xbit(b) as u8,
                k.get_zbit(b) as u8,
            );

            if eps == 0 {
                return None;
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
                return Some((new_word, coeff));
            }
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
///   Inputs are eight *bool-as-u8* flags (0 or 1):
///       (x_a,z_a)  Q qubit-0       (x_b,z_b)  Q qubit-1
///       (x_c,z_c)  P qubit-0       (x_d,z_d)  P qubit-1
///
///   Returns (coeff , x_out0 , z_out0 , x_out1 , z_out1)
///           coeff ∈ { -1 , 0 , +1 }
#[inline(always)]
pub fn comm_2(
    x_a: u8,
    z_a: u8, // Q₀
    x_b: u8,
    z_b: u8, // Q₁
    x_c: u8,
    z_c: u8, // P₀
    x_d: u8,
    z_d: u8, // P₁
) -> (i8, u8, u8, u8, u8) {
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
