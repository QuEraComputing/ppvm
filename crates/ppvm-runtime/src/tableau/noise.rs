use std::fmt::Debug;

use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};

use super::data::{GeneralizedTableau, Tableau};
use super::traits::Measure;
use crate::config::Config;
use crate::tableau::sparsevec::SparseVector;
use crate::tableau::traits::TableauIndex;
use crate::traits::*;
use rand::RngExt;

impl<T: Config> Depolarizing<T> for Tableau<T>
where
    T::Coeff: PartialOrd<f64>,
{
    fn depolarize(&mut self, addr0: usize, p: T::Coeff) {
        debug_assert!(p >= 0.0 && p <= 1.0);
        let r = self.rng.random::<f64>();
        if p <= r {
            return;
        }
        if p > r * 3.0 {
            // p / 3 > r >= 0
            self.x(addr0);
        } else if p > r * 1.5 {
            // 2p/3 > r >= p / 3
            self.y(addr0);
        } else {
            // p > r >= 2p/3
            self.z(addr0);
        }
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> Depolarizing<T>
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: PartialOrd<f64>,
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    fn depolarize(&mut self, addr0: usize, p: T::Coeff) {
        debug_assert!(p >= 0.0 && p <= 1.0);
        let r = self.tableau.rng.random::<f64>();
        if p <= r {
            return;
        }
        if p > r * 3.0 {
            // p / 3 > r >= 0
            self.x(addr0);
        } else if p > r * 1.5 {
            // 2p/3 > r >= p / 3
            self.y(addr0);
        } else {
            // p > r >= 2p/3
            self.z(addr0);
        }
    }
}

impl<T: Config> PauliError<T> for Tableau<T>
where
    T::Coeff: PartialOrd<f64> + Zero,
{
    fn pauli_error(&mut self, addr0: usize, p: [<T as Config>::Coeff; 3]) {
        let r = self.rng.random::<f64>();
        let mut cumulative = T::Coeff::zero();
        for (i, p_) in p.iter().enumerate() {
            cumulative += p_.clone();
            if cumulative > r {
                match i {
                    0 => self.x(addr0),
                    1 => self.y(addr0),
                    _ => self.z(addr0),
                }
                return;
            }
        }
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> PauliError<T>
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: PartialOrd<f64> + Zero,
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    fn pauli_error(&mut self, addr0: usize, p: [<T as Config>::Coeff; 3]) {
        debug_assert!(p.iter().all(|p_| *p_ >= 0.0 && *p_ <= 1.0));
        debug_assert!(p[0].clone() + p[1].clone() + p[2].clone() - 1.0 < 1e-7);
        let r = self.tableau.rng.random::<f64>();
        let mut cumulative = T::Coeff::zero();
        for (i, p_) in p.iter().enumerate() {
            cumulative += p_.clone();
            if cumulative > r {
                match i {
                    0 => self.x(addr0),
                    1 => self.y(addr0),
                    _ => self.z(addr0),
                }
                return;
            }
        }
    }
}

fn two_qubit_pauli_error_impl<T: Config>(
    this: &mut impl Clifford,
    addr0: usize,
    addr1: usize,
    p: [T::Coeff; 15],
    r: f64,
) where
    T::Coeff: PartialOrd<f64> + Zero,
{
    debug_assert!(p.iter().all(|p_| *p_ >= 0.0 && *p_ <= 1.0));
    // debug_assert!(p.iter().sum() - 1.0 < 1e-7);
    let sum = T::Coeff::zero();
    let idx = p
        .iter()
        .scan(sum, |acc, p_| {
            *acc += p_.clone();
            Some(acc.clone())
        })
        .position(|cum_prob| cum_prob > r);

    if let Some(i) = idx {
        #[rustfmt::skip]
        const PAULI_PAIRS: [(u8, u8); 16] = [
            (0,0),(0,1),(0,2),(0,3),
            (1,0),(1,1),(1,2),(1,3),
            (2,0),(2,1),(2,2),(2,3),
            (3,0),(3,1),(3,2),(3,3),
        ];
        let cartesian_index = PAULI_PAIRS[i + 1]; // skip II entry

        match cartesian_index.0 {
            0 => {}
            1 => this.x(addr0),
            2 => this.y(addr0),
            _ => this.z(addr0),
        }

        match cartesian_index.1 {
            0 => {}
            1 => this.x(addr1),
            2 => this.y(addr1),
            _ => this.z(addr1),
        }
    }
}

impl<T: Config> TwoQubitPauliError<T> for Tableau<T>
where
    T::Coeff: PartialOrd<f64> + Zero,
{
    fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [<T as Config>::Coeff; 15]) {
        let r = self.rng.random::<f64>();
        two_qubit_pauli_error_impl::<T>(self, addr0, addr1, p, r);
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> TwoQubitPauliError<T>
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: PartialOrd<f64> + Zero,
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [<T as Config>::Coeff; 15]) {
        if self.is_lost[addr0] || self.is_lost[addr1] {
            return;
        }

        let r = self.tableau.rng.random::<f64>();
        two_qubit_pauli_error_impl::<T>(self, addr0, addr1, p, r);
    }
}

impl<T: Config> Depolarizing2<T> for Tableau<T>
where
    T::Coeff: PartialOrd<f64> + Zero,
{
    fn depolarize2(&mut self, addr0: usize, addr1: usize, p: <T as Config>::Coeff) {
        let r = self.rng.random::<f64>();
        let p_arr: [T::Coeff; 15] = core::array::from_fn(|_| p.clone() * (1.0 / 15.0));
        two_qubit_pauli_error_impl::<T>(self, addr0, addr1, p_arr, r);
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> Depolarizing2<T>
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: PartialOrd<f64> + Zero,
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    fn depolarize2(&mut self, addr0: usize, addr1: usize, p: <T as Config>::Coeff) {
        if self.is_lost[addr0] || self.is_lost[addr1] {
            return;
        }

        let r = self.tableau.rng.random::<f64>();
        let p_arr: [T::Coeff; 15] = core::array::from_fn(|_| p.clone() * (1.0 / 15.0));
        two_qubit_pauli_error_impl::<T>(self, addr0, addr1, p_arr, r);
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> LossChannel<T>
    for GeneralizedTableau<T, I, C>
where
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64> + One + Zero + Clone + num::Num + ToPrimitive + std::fmt::Debug,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat,
    I: Debug,
{
    fn loss_channel(&mut self, addr0: usize, p: <T as Config>::Coeff) {
        if p < self.tableau.rng.random::<f64>() {
            return;
        }

        // NOTE: this is O(n^2) but also potentially removes coefficients, which is nice
        let outcome = self.measure(addr0);
        if outcome {
            // flip back to 0
            self.x(addr0);
        }
        self.is_lost[addr0] = true;
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> ResetLossChannel<T>
    for GeneralizedTableau<T, I, C>
{
    fn reset_loss_channel(&mut self, addr0: usize) {
        self.is_lost[addr0] = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::fxhash::ByteF64;

    type TestConfig = ByteF64<1>;
    type TestTab = GeneralizedTableau<TestConfig>;

    fn tab(n: usize) -> TestTab {
        GeneralizedTableau::new(n, 1e-12)
    }

    // === Depolarizing ===

    #[test]
    fn depolarize_p0_no_change() {
        let mut t = tab(1);
        t.depolarize(0, 0.0);
        assert!(!t.measure(0));
    }

    #[test]
    fn depolarize_p1_does_not_mark_lost() {
        // With p=1.0 an error is always applied; verify is_lost is unaffected
        let mut t = tab(1);
        t.depolarize(0, 1.0);
        assert!(!t.is_lost[0]);
    }

    // === PauliError ===

    #[test]
    fn pauli_error_zero_prob_no_change() {
        let mut t = tab(1);
        t.pauli_error(0, [0.0, 0.0, 0.0]);
        assert!(!t.measure(0));
    }

    #[test]
    fn pauli_error_x_flips_qubit() {
        let mut t = tab(1);
        t.pauli_error(0, [1.0, 0.0, 0.0]); // X|0⟩ = |1⟩
        assert!(t.measure(0));
    }

    #[test]
    fn pauli_error_y_flips_qubit() {
        let mut t = tab(1);
        t.pauli_error(0, [0.0, 1.0, 0.0]); // Y|0⟩ = i|1⟩
        assert!(t.measure(0));
    }

    #[test]
    fn pauli_error_z_no_measurement_change() {
        let mut t = tab(1);
        t.pauli_error(0, [0.0, 0.0, 1.0]); // Z|0⟩ = -|0⟩, still measures 0
        assert!(!t.measure(0));
    }

    #[test]
    fn pauli_error_x_on_excited_qubit_flips_back() {
        let mut t = tab(1);
        t.x(0); // |1⟩
        t.pauli_error(0, [1.0, 0.0, 0.0]); // X|1⟩ = |0⟩
        assert!(!t.measure(0));
    }

    // === TwoQubitPauliError ===

    #[test]
    fn two_qubit_pauli_error_zero_prob_no_change() {
        let mut t = tab(2);
        t.two_qubit_pauli_error(0, 1, [0.0; 15]);
        assert!(!t.measure(0));
        assert!(!t.measure(1));
    }

    #[test]
    fn two_qubit_pauli_error_ix_flips_second_only() {
        // p[0] = 1.0 → IX: I on addr0, X on addr1
        let mut t = tab(2);
        let mut p = [0.0f64; 15];
        p[0] = 1.0;
        t.two_qubit_pauli_error(0, 1, p);
        assert!(!t.measure(0));
        assert!(t.measure(1));
    }

    #[test]
    fn two_qubit_pauli_error_xi_flips_first_only() {
        // p[3] = 1.0 → XI: X on addr0, I on addr1
        let mut t = tab(2);
        let mut p = [0.0f64; 15];
        p[3] = 1.0;
        t.two_qubit_pauli_error(0, 1, p);
        assert!(t.measure(0));
        assert!(!t.measure(1));
    }

    #[test]
    fn two_qubit_pauli_error_xx_flips_both() {
        // p[4] = 1.0 → XX
        let mut t = tab(2);
        let mut p = [0.0f64; 15];
        p[4] = 1.0;
        t.two_qubit_pauli_error(0, 1, p);
        assert!(t.measure(0));
        assert!(t.measure(1));
    }

    #[test]
    fn two_qubit_pauli_error_zz_no_measurement_change() {
        // p[14] = 1.0 → ZZ: Z|0⟩ = -|0⟩ on both, still measures 0
        let mut t = tab(2);
        let mut p = [0.0f64; 15];
        p[14] = 1.0;
        t.two_qubit_pauli_error(0, 1, p);
        assert!(!t.measure(0));
        assert!(!t.measure(1));
    }

    #[test]
    fn two_qubit_pauli_error_both_lost_no_change() {
        let mut t = tab(2);
        t.is_lost[0] = true;
        t.is_lost[1] = true;
        let mut p = [0.0f64; 15];
        p[4] = 1.0; // XX — skipped entirely
        t.two_qubit_pauli_error(0, 1, p);
        assert!(t.is_lost[0]);
        assert!(t.is_lost[1]);
    }

    #[test]
    fn two_qubit_pauli_error_first_lost_no_apply() {
        // addr0 lost; p[0] = 1.0 (IX) → marginal p_x for addr1 = 1.0
        let mut t = tab(2);
        t.is_lost[0] = true;
        let mut p = [0.0f64; 15];
        p[0] = 1.0; // IX
        t.two_qubit_pauli_error(0, 1, p);
        assert!(!t.measure(1)); // nothing applied to addr1
    }

    // === Depolarizing2 ===

    #[test]
    fn depolarize2_p0_no_change() {
        let mut t = tab(2);
        t.depolarize2(0, 1, 0.0);
        assert!(!t.measure(0));
        assert!(!t.measure(1));
    }

    #[test]
    fn depolarize2_both_lost_no_change() {
        let mut t = tab(2);
        t.is_lost[0] = true;
        t.is_lost[1] = true;
        t.depolarize2(0, 1, 1.0);
        assert!(t.is_lost[0]);
        assert!(t.is_lost[1]);
    }

    #[test]
    fn depolarize2_first_lost_p0_second_unchanged() {
        let mut t = tab(2);
        t.is_lost[0] = true;
        t.depolarize2(0, 1, 0.0); // effective p on addr1 = 4/5 * 0 = 0
        assert!(!t.measure(1));
    }

    #[test]
    fn depolarize2_second_lost_p0_first_unchanged() {
        let mut t = tab(2);
        t.is_lost[1] = true;
        t.depolarize2(0, 1, 0.0); // effective p on addr0 = 4/5 * 0 = 0
        assert!(!t.measure(0));
    }

    // === LossChannel ===

    #[test]
    fn loss_channel_p0_qubit_not_lost() {
        let mut t = tab(1);
        t.loss_channel(0, 0.0);
        assert!(!t.is_lost[0]);
    }

    #[test]
    fn loss_channel_p1_qubit_marked_lost() {
        let mut t = tab(1);
        t.loss_channel(0, 1.0);
        assert!(t.is_lost[0]);
    }

    #[test]
    fn loss_channel_p1_qubit_reset_to_zero() {
        // Qubit starts in |1⟩; loss_channel should measure, reset to |0⟩, then mark lost
        let mut t = tab(1);
        t.x(0);
        t.loss_channel(0, 1.0);
        assert!(t.is_lost[0]);
        assert!(!t.measure(0)); // Reset to |0⟩ before marking lost
    }

    #[test]
    fn loss_channel_p1_subsequent_gate_is_noop() {
        let mut t = tab(1);
        t.loss_channel(0, 1.0);
        t.x(0); // No-op: qubit is lost
        assert!(!t.measure(0)); // Still |0⟩
    }

    #[test]
    fn loss_channel_p0_second_qubit_unaffected() {
        let mut t = tab(2);
        t.loss_channel(0, 0.0);
        t.loss_channel(1, 0.0);
        assert!(!t.is_lost[0]);
        assert!(!t.is_lost[1]);
    }

    // === ResetLossChannel ===

    #[test]
    fn reset_loss_channel_clears_lost_flag() {
        let mut t = tab(1);
        t.loss_channel(0, 1.0);
        assert!(t.is_lost[0]);
        t.reset_loss_channel(0);
        assert!(!t.is_lost[0]);
    }

    #[test]
    fn reset_loss_channel_qubit_in_ground_state() {
        // loss_channel resets qubit to |0⟩; after reset_loss_channel it should still be |0⟩
        let mut t = tab(1);
        t.x(0); // |1⟩
        t.loss_channel(0, 1.0); // measures, resets to |0⟩, marks lost
        t.reset_loss_channel(0);
        assert!(!t.measure(0)); // back in |0⟩
    }

    #[test]
    fn reset_loss_channel_gates_work_again() {
        let mut t = tab(1);
        t.loss_channel(0, 1.0);
        t.reset_loss_channel(0);
        t.x(0); // should no longer be a no-op
        assert!(t.measure(0));
    }

    // === Statistical tests ===

    #[test]
    fn depolarize_statistics() {
        // Starting from |0⟩, P(measure 1) = P(X) + P(Y) = p/3 + p/3 = 2p/3.
        // Z leaves |0⟩ unchanged; I leaves |0⟩ unchanged.
        let p = 0.6_f64;
        let expected = 2.0 * p / 3.0; // 0.4
        let trials = 500;

        let ones = (0..trials)
            .filter(|_| {
                let mut t = tab(1);
                t.depolarize(0, p);
                t.measure(0)
            })
            .count();

        let fraction = ones as f64 / trials as f64;
        // tolerance ~5σ: σ = sqrt(expected*(1-expected)/trials) ≈ 0.022
        assert!(
            (fraction - expected).abs() < 0.1,
            "Expected fraction {expected:.3}, got {fraction:.3}"
        );
    }

    #[test]
    fn depolarize2_statistics() {
        // Starting from |00⟩, errors that flip qubit 0 to |1⟩ are X and Y on that qubit:
        // XI, XX, XY, XZ, YI, YX, YY, YZ — 8 out of 15, so P(q0=1) = 8p/15.
        let p = 0.6_f64;
        let expected = 8.0 * p / 15.0; // 0.32
        let trials = 500;

        let ones = (0..trials)
            .filter(|_| {
                let mut t = tab(2);
                t.depolarize2(0, 1, p);
                t.measure(0)
            })
            .count();

        let fraction = ones as f64 / trials as f64;
        // tolerance ~5σ: σ = sqrt(expected*(1-expected)/trials) ≈ 0.021
        assert!(
            (fraction - expected).abs() < 0.1,
            "Expected fraction {expected:.3}, got {fraction:.3}"
        );
    }

    #[test]
    fn test_cnot() {
        let mut t = tab(2);
        t.x(0);
        t.cnot(0, 1);
        t.loss_channel(0, 1.0);
        assert!(!t.measure(0));
        assert!(t.measure(1));

        let mut t = tab(2);
        t.loss_channel(0, 1.0);
        t.x(0);
        t.cnot(0, 1);
        assert!(!t.measure(1));
        assert!(!t.measure(0));
    }

    #[test]
    fn test_ghz_statistics() {
        let mut t = tab(2);
        t.h(0);
        t.cnot(0, 1);

        let trials = 100u64;
        let mut z_avg = 0.0;
        let p = 0.1;
        for i in 0..trials {
            let mut t_trial = t.fork(Some(i));
            t_trial.loss_channel(0, p);

            let outcome0 = t_trial.measure(0);
            let outcome1 = t_trial.measure(1);
            if outcome0 == outcome1 {
                z_avg += 1.0 / trials as f64;
            } else {
                z_avg += -1.0 / trials as f64;
            }
        }

        println!("{}", z_avg);
        assert!((z_avg - (1.0 - p)).abs() < 10.0 / trials as f64);
    }
}
