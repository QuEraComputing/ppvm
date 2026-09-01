// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use fxhash::{FxHashMap, FxHashSet};
use num::Complex;
use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::{HashFinalize, PauliStorage};
use std::f64::consts::PI;
use std::hash::BuildHasher;

use super::group::TranslationGroup;

impl TranslationGroup {
    /// Integer numerator of the character phase before division by
    /// [`Self::phase_modulus`]: `Σ_g (k[g] · counter[g] / orders[g]) mod 1`
    /// expressed as an integer in `[0, phase_modulus)`.
    pub(super) fn character_numerator(&self, k_modes: &[i32], counter: &[u32]) -> usize {
        assert_eq!(
            k_modes.len(),
            self.n_generators(),
            "k_modes length mismatch"
        );
        assert_eq!(
            counter.len(),
            self.n_generators(),
            "counter length mismatch"
        );
        let modulus = self.phase_modulus() as u128;
        let mut numerator = 0u128;
        for g in 0..self.n_generators() {
            let order = self.generator_order(g);
            let k = (k_modes[g] as i64).rem_euclid(order as i64) as u128;
            let count = (counter[g] % order) as u128;
            let reduced = (k * count) % order as u128;
            let factor = self.phase_modulus() as u128 / order as u128;
            numerator = (numerator + reduced * factor) % modulus;
        }
        numerator as usize
    }

    /// Momentum-sector character `χ_k(g) = exp(i Σ_g 2π · k[g] · counter[g] / orders[g])`
    /// where `k[g] ∈ ℤ` is the integer momentum mode along generator `g`
    /// (the corresponding wavenumber is `2π · k[g] / orders[g]`).
    ///
    /// `k.len()` must equal `self.n_generators()`. The character of the
    /// identity element (`counter = [0, …]`) is `1`. For the trivial
    /// (`k = [0, …]`) sector all characters are `1`.
    pub fn character(&self, k_modes: &[i32], counter: &[u32]) -> Complex<f64> {
        let numerator = self.character_numerator(k_modes, counter);
        let phase = 2.0 * PI * numerator as f64 / self.phase_modulus() as f64;
        Complex::from_polar(1.0, phase)
    }
}

/// Replace `(basis, complex_coeffs)` in-place with the orbit-rep form
/// **projected onto momentum sector `k_modes`**.
///
/// For each represented orbit, coefficients on its **distinct** orbit
/// members are averaged with the momentum character weight:
/// `(1/|orbit|) · Σ_{p ∈ orbit} χ_k(g_p) · c_p` where `g_p` is the group
/// element such that `g_p · rep = p`.
///
/// Orbits whose stabilizer is incompatible with `k_modes` (the same orbit
/// member is reached with different character numerators) project to zero
/// and are omitted from the output.
///
/// If the input was already a momentum-`k_modes` eigenstate (i.e. the
/// coefficients satisfy `c_{g·p} = χ_k(g)⁻¹ · c_p` for every orbit),
/// the output is the orbit-rep coefficients of that state unchanged.
/// Otherwise the projection discards the components in other sectors —
/// use [`check_momentum_sector`] beforehand to validate.
///
/// For the `k_modes = [0, 0, …]` (trivial) sector all characters are `1`,
/// so projection averages the distinct orbit members onto each rep. This
/// differs from plain [`super::canonicalize_pauli_sum`], whose real-coefficient
/// merging sums collisions without orbit-size normalization.
pub fn canonicalize_pauli_sum_complex<A, S, const R: bool>(
    basis: &mut Vec<PauliWord<A, S, R>>,
    coeffs: &mut Vec<Complex<f64>>,
    group: &TranslationGroup,
    k_modes: &[i32],
) where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    assert_eq!(
        basis.len(),
        coeffs.len(),
        "basis and coeffs length mismatch"
    );
    assert_eq!(
        k_modes.len(),
        group.n_generators(),
        "k_modes length {} != number of generators {}",
        k_modes.len(),
        group.n_generators()
    );
    let mut input: FxHashMap<PauliWord<A, S, R>, Complex<f64>> = FxHashMap::default();
    for (word, &coeff) in basis.iter().zip(coeffs.iter()) {
        *input.entry(*word).or_insert(Complex::new(0.0, 0.0)) += coeff;
    }
    let reps: FxHashSet<_> = input.keys().map(|word| group.canonicalize(word)).collect();
    let mut projected = FxHashMap::default();

    for rep in reps {
        let mut members: FxHashMap<_, (Vec<u32>, usize)> = FxHashMap::default();
        let mut compatible = true;
        for (member, counter) in group.orbit_with_counters(&rep) {
            let numerator = group.character_numerator(k_modes, &counter);
            match members.entry(member) {
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert((counter, numerator));
                }
                std::collections::hash_map::Entry::Occupied(entry) => {
                    if entry.get().1 != numerator {
                        compatible = false;
                        break;
                    }
                }
            }
        }
        if !compatible {
            continue;
        }
        let orbit_size = members.len() as f64;
        let mut rep_coeff = Complex::new(0.0, 0.0);
        for (member, (counter, _)) in members {
            let coeff = input
                .get(&member)
                .copied()
                .unwrap_or(Complex::new(0.0, 0.0));
            rep_coeff += group.character(k_modes, &counter) * coeff / orbit_size;
        }
        projected.insert(rep, rep_coeff);
    }
    basis.clear();
    coeffs.clear();
    basis.reserve(projected.len());
    coeffs.reserve(projected.len());
    for (w, c) in projected {
        basis.push(w);
        coeffs.push(c);
    }
}

/// Verify that a `(basis, complex_coeffs)` Pauli sum lies entirely in
/// the momentum sector `k_modes` under `group`.
///
/// Concretely: for every orbit represented in the basis, all members
/// must satisfy `c_{g·r} = χ_k(g)⁻¹ · c_r` for some choice of orbit-rep
/// coefficient `c_r`. Orbit members absent from `basis` are treated as
/// having coefficient zero, rather than being ignored.
///
/// An orbit with a stabilizer incompatible with `k_modes` cannot carry
/// that sector and fails with [`SectorCheckError::IncompatibleStabilizer`];
/// the corresponding momentum projection would be zero.
///
/// Returns `Ok(())` on pass; `Err(SectorCheckError)` on fail with the
/// offending orbit-rep, expected coefficient, and actual coefficient.
///
/// Use this on a user-supplied initial state before feeding it to a
/// phase-aware merging pipeline — silently projecting a wrongly-typed
/// input throws away meaningful physics.
pub fn check_momentum_sector<A, S, const R: bool>(
    basis: &[PauliWord<A, S, R>],
    coeffs: &[Complex<f64>],
    group: &TranslationGroup,
    k_modes: &[i32],
    tol: f64,
) -> Result<(), SectorCheckError<A, S, R>>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    assert_eq!(basis.len(), coeffs.len());
    assert_eq!(k_modes.len(), group.n_generators());

    if !tol.is_finite() || tol < 0.0 {
        return Err(SectorCheckError::InvalidTolerance { tol });
    }
    let mut input: FxHashMap<PauliWord<A, S, R>, Complex<f64>> = FxHashMap::default();
    for (pauli, &coeff) in basis.iter().zip(coeffs.iter()) {
        if !coeff.re.is_finite() || !coeff.im.is_finite() {
            return Err(SectorCheckError::NonFiniteCoefficient {
                pauli: *pauli,
                coeff,
            });
        }
        *input.entry(*pauli).or_insert(Complex::new(0.0, 0.0)) += coeff;
    }
    for (&pauli, &coeff) in &input {
        if !coeff.re.is_finite() || !coeff.im.is_finite() {
            return Err(SectorCheckError::NonFiniteCoefficient { pauli, coeff });
        }
    }
    input.retain(|_, coeff| *coeff != Complex::new(0.0, 0.0));
    let reps: FxHashSet<_> = input
        .keys()
        .map(|pauli| group.canonicalize(pauli))
        .collect();

    for rep in reps {
        let mut members: FxHashMap<_, (Vec<u32>, usize)> = FxHashMap::default();
        for (member, counter) in group.orbit_with_counters(&rep) {
            let numerator = group.character_numerator(k_modes, &counter);
            match members.entry(member) {
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert((counter, numerator));
                }
                std::collections::hash_map::Entry::Occupied(entry) => {
                    if entry.get().1 != numerator {
                        return Err(SectorCheckError::IncompatibleStabilizer {
                            rep,
                            shift: counter,
                        });
                    }
                }
            }
        }

        let (reference_word, (reference_counter, _)) = members
            .iter()
            .find(|(member, _)| input.contains_key(*member))
            .expect("represented orbit has a nonzero member");
        let rep_coeff = group.character(k_modes, reference_counter) * input[reference_word];

        for (member, (counter, _)) in members {
            let expected = group.character(k_modes, &counter).conj() * rep_coeff;
            let actual = input
                .get(&member)
                .copied()
                .unwrap_or(Complex::new(0.0, 0.0));
            if (actual - expected).norm() > tol * rep_coeff.norm().max(1.0) {
                return Err(SectorCheckError::CoefficientMismatch {
                    rep,
                    offending_pauli: member,
                    expected,
                    actual,
                    shift: counter,
                });
            }
        }
    }
    Ok(())
}

/// Detail report for a failed [`check_momentum_sector`].
pub enum SectorCheckError<A: PauliStorage, S, const R: bool> {
    InvalidTolerance {
        tol: f64,
    },
    NonFiniteCoefficient {
        pauli: PauliWord<A, S, R>,
        coeff: Complex<f64>,
    },
    CoefficientMismatch {
        rep: PauliWord<A, S, R>,
        offending_pauli: PauliWord<A, S, R>,
        expected: Complex<f64>,
        actual: Complex<f64>,
        shift: Vec<u32>,
    },
    IncompatibleStabilizer {
        rep: PauliWord<A, S, R>,
        shift: Vec<u32>,
    },
}

impl<A: PauliStorage, S, const R: bool> std::fmt::Debug for SectorCheckError<A, S, R>
where
    S: BuildHasher + Clone + Default + HashFinalize,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidTolerance { tol } => {
                write!(f, "SectorCheckError::InvalidTolerance {{ tol: {tol:?} }}")
            }
            Self::NonFiniteCoefficient { pauli, coeff } => write!(
                f,
                "SectorCheckError::NonFiniteCoefficient {{ pauli: {pauli}, coeff: {coeff:?} }}"
            ),
            Self::CoefficientMismatch {
                rep,
                offending_pauli,
                expected,
                actual,
                shift,
            } => write!(
                f,
                "SectorCheckError::CoefficientMismatch {{ rep: {rep}, offending_pauli: \
                 {offending_pauli}, expected: {expected:?}, actual: {actual:?}, shift: \
                 {shift:?} }}"
            ),
            Self::IncompatibleStabilizer { rep, shift } => write!(
                f,
                "SectorCheckError::IncompatibleStabilizer {{ rep: {rep}, shift: {shift:?} }}"
            ),
        }
    }
}

impl<A: PauliStorage, S, const R: bool> std::fmt::Display for SectorCheckError<A, S, R>
where
    S: BuildHasher + Clone + Default + HashFinalize,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidTolerance { tol } => write!(
                f,
                "invalid tolerance {tol:?}: must be finite and non-negative"
            ),
            Self::NonFiniteCoefficient { pauli, coeff } => {
                write!(f, "non-finite coefficient {coeff:?} on Pauli word {pauli}")
            }
            Self::CoefficientMismatch {
                rep,
                offending_pauli,
                expected,
                actual,
                shift,
            } => write!(
                f,
                "input not in target momentum sector: orbit rep {rep} expected c={expected:?}, \
                 but orbit member {offending_pauli} (shift {shift:?}) has c={actual:?}"
            ),
            Self::IncompatibleStabilizer { rep, shift } => write!(
                f,
                "stabilizer incompatible with momentum sector: orbit rep {rep} has conflicting \
                 character numerators (shift {shift:?})"
            ),
        }
    }
}
