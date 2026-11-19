// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The exact dense backend: `2^k` amplitudes, bit `j` of an index is axis `j`.

use num_complex::Complex64;

use super::axis::{AxisPauli, RotationCollapse, bits};
use super::state::ActiveState;

#[derive(Clone, Debug, PartialEq)]
pub struct DenseState {
    amplitudes: Vec<Complex64>,
}

impl Default for DenseState {
    /// No axes: the single amplitude 1.
    fn default() -> Self {
        Self {
            amplitudes: vec![Complex64::new(1.0, 0.0)],
        }
    }
}

/// `i^e`.
fn i_pow(e: u32) -> Complex64 {
    [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 1.0),
        Complex64::new(-1.0, 0.0),
        Complex64::new(0.0, -1.0),
    ][(e % 4) as usize]
}

/// `c` with `P|α⟩ = c · |α ⊕ x⟩`.
fn pauli_factor(p: AxisPauli, alpha: usize) -> Complex64 {
    let negative = (alpha as u64 & p.z).count_ones() % 2 == 1;
    let c = i_pow((p.x & p.z).count_ones());
    if negative { -c } else { c }
}

impl DenseState {
    pub fn amplitudes(&self) -> &[Complex64] {
        &self.amplitudes
    }

    /// `φ ← G† φ` for each gate of `v` in order.
    fn apply_rotation_dagger(&mut self, v: &RotationCollapse) {
        let p = 1usize << v.pivot;
        for t in bits(v.fan_out) {
            self.cx(p, 1 << t);
        }
        for c in bits(v.fan_in) {
            self.cx(1 << c, p);
        }
        for t in bits(v.cz) {
            let both = p | 1 << t;
            self.for_each_index(|a, amp| {
                if a & both == both {
                    *amp = -*amp;
                }
            });
        }
        if v.s {
            // S† = diag(1, −i).
            self.for_each_index(|a, amp| {
                if a & p != 0 {
                    *amp *= Complex64::new(0.0, -1.0);
                }
            });
        }
        if v.h {
            let r = std::f64::consts::FRAC_1_SQRT_2;
            for a in (0..self.amplitudes.len()).filter(|a| a & p == 0) {
                let (lo, hi) = (self.amplitudes[a], self.amplitudes[a | p]);
                self.amplitudes[a] = (lo + hi) * r;
                self.amplitudes[a | p] = (lo - hi) * r;
            }
        }
        if v.x {
            for a in (0..self.amplitudes.len()).filter(|a| a & p == 0) {
                self.amplitudes.swap(a, a | p);
            }
        }
    }

    /// CX with control mask `c` and target mask `t` (its own inverse).
    fn cx(&mut self, c: usize, t: usize) {
        for a in (0..self.amplitudes.len()).filter(|a| a & c != 0 && a & t == 0) {
            self.amplitudes.swap(a, a | t);
        }
    }

    fn for_each_index(&mut self, mut f: impl FnMut(usize, &mut Complex64)) {
        for (a, amp) in self.amplitudes.iter_mut().enumerate() {
            f(a, amp);
        }
    }
}

impl ActiveState for DenseState {
    fn num_axes(&self) -> usize {
        self.amplitudes.len().trailing_zeros() as usize
    }

    fn activate(&mut self) {
        let len = self.amplitudes.len();
        self.amplitudes.resize(2 * len, Complex64::new(0.0, 0.0));
    }

    fn rotate(&mut self, p: AxisPauli, a: Complex64, b: Complex64) {
        let x = p.x as usize;
        if x == 0 {
            self.for_each_index(|alpha, amp| *amp *= a + b * pauli_factor(p, alpha));
            return;
        }
        // Pair α with β = α ⊕ x: each new amplitude mixes the two old ones.
        for alpha in 0..self.amplitudes.len() {
            let beta = alpha ^ x;
            if alpha < beta {
                let (ua, ub) = (self.amplitudes[alpha], self.amplitudes[beta]);
                self.amplitudes[alpha] = a * ua + b * pauli_factor(p, beta) * ub;
                self.amplitudes[beta] = a * ub + b * pauli_factor(p, alpha) * ua;
            }
        }
    }

    fn branch_weights(&self, p: AxisPauli) -> [f64; 2] {
        let x = p.x as usize;
        // ‖½(1 ± P)φ‖² summed directly, with (Pφ)[α] = c(α ⊕ x) · φ[α ⊕ x]: never
        // negative, and an impossible branch stays O(ε²) instead of ±ε.
        let (mut plus, mut minus) = (0.0, 0.0);
        for (alpha, &amp) in self.amplitudes.iter().enumerate() {
            let p_amp = pauli_factor(p, alpha ^ x) * self.amplitudes[alpha ^ x];
            plus += (amp + p_amp).norm_sqr();
            minus += (amp - p_amp).norm_sqr();
        }
        [0.25 * plus, 0.25 * minus]
    }

    fn collapse(&mut self, p: AxisPauli, m: bool, v: &RotationCollapse) {
        // Project onto eigenvalue (−1)^m: φ ← ½(φ ± P φ).
        let sign = if m { -0.5 } else { 0.5 };
        self.rotate(p, Complex64::new(0.5, 0.0), Complex64::new(sign, 0.0));
        self.apply_rotation_dagger(v);
        // The pivot is |0⟩ now: keep those amplitudes, removing its bit.
        let h = v.pivot;
        let low = (1usize << h) - 1;
        let kept: Vec<Complex64> = (0..self.amplitudes.len() / 2)
            .map(|a| self.amplitudes[(a & low) | (a & !low) << 1])
            .collect();
        let norm = kept.iter().map(|a| a.norm_sqr()).sum::<f64>().sqrt();
        assert!(norm > 0.0, "collapsed onto a branch of zero weight");
        self.amplitudes = kept.into_iter().map(|a| a / norm).collect();
    }
}

#[cfg(test)]
mod tests {
    //! Checked against plain matrices: Paulis as tensor products of 2 × 2
    //! matrices, `V` as a product of gate matrices. Bit `j` of an index is axis
    //! `j`. No code is shared with the backend.
    use super::*;

    type C = Complex64;
    type Mat = Vec<Vec<C>>;

    fn c(re: f64, im: f64) -> C {
        C::new(re, im)
    }

    fn identity(d: usize) -> Mat {
        (0..d)
            .map(|i| {
                (0..d)
                    .map(|j| c(f64::from(u8::from(i == j)), 0.0))
                    .collect()
            })
            .collect()
    }

    fn mul(a: &Mat, b: &Mat) -> Mat {
        let d = a.len();
        (0..d)
            .map(|i| {
                (0..d)
                    .map(|j| (0..d).map(|k| a[i][k] * b[k][j]).sum())
                    .collect()
            })
            .collect()
    }

    fn dagger(a: &Mat) -> Mat {
        let d = a.len();
        (0..d)
            .map(|i| (0..d).map(|j| a[j][i].conj()).collect())
            .collect()
    }

    fn apply(a: &Mat, v: &[C]) -> Vec<C> {
        a.iter()
            .map(|row| row.iter().zip(v).map(|(x, y)| x * y).sum())
            .collect()
    }

    /// The one-qubit matrix `u` on axis `j` of `k`.
    fn on_axis(u: [[C; 2]; 2], j: usize, k: usize) -> Mat {
        let d = 1 << k;
        (0..d)
            .map(|r| {
                (0..d)
                    .map(|s| {
                        if (r ^ s) & !(1 << j) == 0 {
                            u[r >> j & 1][s >> j & 1]
                        } else {
                            c(0.0, 0.0)
                        }
                    })
                    .collect()
            })
            .collect()
    }

    fn pauli(p: AxisPauli, k: usize) -> Mat {
        let (o, l, i) = (c(0.0, 0.0), c(1.0, 0.0), c(0.0, 1.0));
        (0..k).fold(identity(1 << k), |m, j| {
            let u = match (p.x >> j & 1, p.z >> j & 1) {
                (0, 0) => return m,
                (1, 0) => [[o, l], [l, o]],
                (0, _) => [[l, o], [o, -l]],
                _ => [[o, -i], [i, o]],
            };
            mul(&m, &on_axis(u, j, k))
        })
    }

    /// `V` = the product of `v`'s gate matrices in gate order.
    fn rotation_matrix(v: &RotationCollapse, k: usize) -> Mat {
        let d = 1 << k;
        let (o, l, i) = (c(0.0, 0.0), c(1.0, 0.0), c(0.0, 1.0));
        let r = std::f64::consts::FRAC_1_SQRT_2;
        let cx = |ctl: usize, tgt: usize| -> Mat {
            (0..d)
                .map(|row| {
                    (0..d)
                        .map(|col| {
                            let image = if col >> ctl & 1 == 1 {
                                col ^ 1 << tgt
                            } else {
                                col
                            };
                            c(f64::from(u8::from(row == image)), 0.0)
                        })
                        .collect()
                })
                .collect()
        };
        let cz = |a: usize, b: usize| -> Mat {
            (0..d)
                .map(|row| {
                    (0..d)
                        .map(|col| {
                            let s = if col >> a & 1 == 1 && col >> b & 1 == 1 {
                                -1.0
                            } else {
                                1.0
                            };
                            c(if row == col { s } else { 0.0 }, 0.0)
                        })
                        .collect()
                })
                .collect()
        };
        let mut gates = Vec::new();
        bits(v.fan_out).for_each(|t| gates.push(cx(v.pivot, t)));
        bits(v.fan_in).for_each(|s| gates.push(cx(s, v.pivot)));
        bits(v.cz).for_each(|t| gates.push(cz(v.pivot, t)));
        if v.s {
            gates.push(on_axis([[l, o], [o, i]], v.pivot, k));
        }
        if v.h {
            gates.push(on_axis(
                [[c(r, 0.0), c(r, 0.0)], [c(r, 0.0), c(-r, 0.0)]],
                v.pivot,
                k,
            ));
        }
        if v.x {
            gates.push(on_axis([[o, l], [l, o]], v.pivot, k));
        }
        gates.iter().fold(identity(d), |m, g| mul(&m, g))
    }

    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn unit(&mut self) -> f64 {
            (self.next() >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        }

        fn state(&mut self, k: usize) -> DenseState {
            let amps: Vec<C> = (0..1 << k).map(|_| c(self.unit(), self.unit())).collect();
            let norm = amps.iter().map(|a| a.norm_sqr()).sum::<f64>().sqrt();
            DenseState {
                amplitudes: amps.into_iter().map(|a| a / norm).collect(),
            }
        }

        fn pauli(&mut self, k: usize) -> AxisPauli {
            loop {
                let mask = (1u64 << k) - 1;
                let p = AxisPauli {
                    x: self.next() & mask,
                    z: self.next() & mask,
                };
                if !p.is_identity() {
                    return p;
                }
            }
        }
    }

    fn assert_close(a: &[C], b: &[C]) {
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b) {
            assert!((x - y).norm() < 1e-12, "{a:?} != {b:?}");
        }
    }

    #[test]
    fn activate_appends_zero() {
        let mut rng = Rng(3);
        let phi = rng.state(3);
        let mut grown = phi.clone();
        grown.activate();
        assert_eq!(grown.num_axes(), 4);
        let zeros = vec![c(0.0, 0.0); 8];
        assert_close(&grown.amplitudes[..8], &phi.amplitudes);
        assert_close(&grown.amplitudes[8..], &zeros);
    }

    #[test]
    fn rotate_matches_matrices() {
        let mut rng = Rng(5);
        for k in 1..=5 {
            for _ in 0..20 {
                let (phi, p) = (rng.state(k), rng.pauli(k));
                let (a, b) = (c(rng.unit(), rng.unit()), c(rng.unit(), rng.unit()));
                let mut got = phi.clone();
                got.rotate(p, a, b);
                let pm = pauli(p, k);
                let op: Mat = (0..1 << k)
                    .map(|i| {
                        (0..1 << k)
                            .map(|j| a * identity(1 << k)[i][j] + b * pm[i][j])
                            .collect()
                    })
                    .collect();
                assert_close(&got.amplitudes, &apply(&op, &phi.amplitudes));
            }
        }
    }

    #[test]
    fn branch_weights_match_matrices() {
        let mut rng = Rng(7);
        for k in 1..=5 {
            for _ in 0..20 {
                let (phi, p) = (rng.state(k), rng.pauli(k));
                let pv = apply(&pauli(p, k), &phi.amplitudes);
                let e: C = phi
                    .amplitudes
                    .iter()
                    .zip(&pv)
                    .map(|(x, y)| x.conj() * y)
                    .sum();
                let [w0, w1] = phi.branch_weights(p);
                assert!(
                    (w0 - 0.5 * (1.0 + e.re)).abs() < 1e-12
                        && (w1 - 0.5 * (1.0 - e.re)).abs() < 1e-12
                );
            }
        }
    }

    /// On an eigenstate of `P`, up to rounding, the weights are probabilities:
    /// never negative, summing to 1, with the impossible branch at rounding
    /// level.
    #[test]
    fn branch_weights_on_eigenstates() {
        let mut rng = Rng(13);
        for k in 1..=5 {
            for _ in 0..200 {
                let (psi, p, q) = (rng.state(k), rng.pauli(k), rng.pauli(k));
                if p.x == 0 && p.z == 0 {
                    continue;
                }
                let sign = if rng.next() & 1 == 1 { -1.0 } else { 1.0 };
                let pv = apply(&pauli(p, k), &psi.amplitudes);
                let amps: Vec<C> = psi
                    .amplitudes
                    .iter()
                    .zip(&pv)
                    .map(|(a, b)| a + sign * b)
                    .collect();
                let norm = amps.iter().map(|a| a.norm_sqr()).sum::<f64>().sqrt();
                let mut phi = DenseState {
                    amplitudes: amps.into_iter().map(|a| a / norm).collect(),
                };
                // e^{−iθQ} then e^{iθQ}: the same state, with rounding.
                let theta = rng.unit();
                let (cos, sin) = (theta.cos(), theta.sin());
                phi.rotate(q, c(cos, 0.0), c(0.0, -sin));
                phi.rotate(q, c(cos, 0.0), c(0.0, sin));
                let [w0, w1] = phi.branch_weights(p);
                assert!(w0 >= 0.0 && w1 >= 0.0, "{w0} {w1}");
                assert!((w0 + w1 - 1.0).abs() < 1e-12);
                assert!(w0.min(w1) < 1e-12);
            }
        }
    }

    /// `V† ½(1 + (−1)^m P) V` is the projector onto "pivot in |0⟩".
    #[test]
    fn rotation_maps_the_branch_to_pivot_zero() {
        let mut rng = Rng(11);
        for k in 1..=5 {
            for _ in 0..30 {
                let (p, m) = (rng.pauli(k), rng.next() & 1 == 1);
                let v = RotationCollapse::collapse(p, m);
                let vm = rotation_matrix(&v, k);
                let (pm, id) = (pauli(p, k), identity(1 << k));
                let s = if m { -0.5 } else { 0.5 };
                let proj: Mat = (0..1 << k)
                    .map(|i| (0..1 << k).map(|j| 0.5 * id[i][j] + s * pm[i][j]).collect())
                    .collect();
                let got = mul(&mul(&dagger(&vm), &proj), &vm);
                let z = pauli(
                    AxisPauli {
                        x: 0,
                        z: 1 << v.pivot,
                    },
                    k,
                );
                for i in 0..1 << k {
                    for j in 0..1 << k {
                        let want = 0.5 * id[i][j] + 0.5 * z[i][j];
                        assert!(
                            (got[i][j] - want).norm() < 1e-12,
                            "k = {k}, p = {p:?}, m = {m}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn collapse_matches_matrices() {
        let mut rng = Rng(13);
        for k in 1..=5 {
            for _ in 0..30 {
                let (phi, p) = (rng.state(k), rng.pauli(k));
                let m = rng.next() & 1 == 1;
                if phi.branch_weights(p)[usize::from(m)] < 1e-6 {
                    continue;
                }
                let v = RotationCollapse::collapse(p, m);
                let mut got = phi.clone();
                got.collapse(p, m, &v);
                // Project, apply V†, keep the pivot's |0⟩ half, renormalize.
                let s = if m { -0.5 } else { 0.5 };
                let pv = apply(&pauli(p, k), &phi.amplitudes);
                let projected: Vec<C> = phi
                    .amplitudes
                    .iter()
                    .zip(&pv)
                    .map(|(x, y)| 0.5 * x + s * y)
                    .collect();
                let rotated = apply(&dagger(&rotation_matrix(&v, k)), &projected);
                let kept: Vec<C> = (0..1usize << k)
                    .filter(|a| a >> v.pivot & 1 == 0)
                    .map(|a| rotated[a])
                    .collect();
                let norm = kept.iter().map(|a| a.norm_sqr()).sum::<f64>().sqrt();
                let want: Vec<C> = kept.iter().map(|a| a / norm).collect();
                assert_eq!(got.num_axes(), k - 1);
                assert_close(&got.amplitudes, &want);
            }
        }
    }
}
