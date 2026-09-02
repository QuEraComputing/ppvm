// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::sum::PauliSum;
use fxhash::{FxHashMap, FxHashSet};
use num::Complex;
use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::{ACMapAddAssign, ACMapBase, ACMapIter, Config, HashFinalize, PauliStorage};
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

    /// Everything the phase-aware routines need about `w`'s orbit in
    /// momentum sector `k_modes`, from ONE orbit traversal: the lex-min
    /// representative `r`, the mixed-radix counter of the group element
    /// mapping `r` to `w` (as [`Self::canonicalize_with_shift`]), and the
    /// number of **distinct** orbit members `|orbit|`.
    ///
    /// Returns `None` when the orbit's stabilizer is incompatible with
    /// `k_modes` — i.e. some `s` with `s·w = w` has `χ_k(s) ≠ 1`. Such an
    /// orbit cannot carry this sector: its momentum projection is
    /// identically zero, and the rep coefficient a single traversal would
    /// report depends on which counter the traversal happens to pick.
    ///
    /// `|orbit| = |G| / |stabilizer|` (orbit-stabilizer), and equals
    /// `|G|` only for free orbits.
    ///
    /// Same `O(|G| × n_qubits)` cost as [`Self::canonicalize_with_shift`].
    pub fn canonicalize_in_sector<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
        k_modes: &[i32],
    ) -> Option<(PauliWord<A, S, R>, Vec<u32>, usize)>
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        let mut best: Option<(PauliWord<A, S, R>, Vec<u32>)> = None;
        let mut stabilizer = 0usize;
        for (candidate, counter) in self.orbit_with_counters(w) {
            if candidate == *w {
                if self.character_numerator(k_modes, &counter) != 0 {
                    return None;
                }
                stabilizer += 1;
            }
            if best.as_ref().is_none_or(|(b, _)| candidate < *b) {
                best = Some((candidate, counter));
            }
        }
        let (rep, counter_from_word) = best.expect("a finite group contains the identity element");
        let shift = (0..self.n_generators())
            .map(|g| {
                let order = self.generator_order(g);
                (order - counter_from_word[g]) % order
            })
            .collect();
        Some((rep, shift, self.order() / stabilizer))
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
    let projected = project_onto_reps(&input, group, k_modes);
    basis.clear();
    coeffs.clear();
    basis.reserve(projected.len());
    coeffs.reserve(projected.len());
    for (word, (sum, orbit_size)) in projected {
        basis.push(word);
        coeffs.push(sum / orbit_size as f64);
    }
}

/// Character-weighted fold of `input` onto translation-orbit
/// representatives, the shared core of the two momentum-projection
/// conventions.
///
/// Returns `rep → (Σ_{p ∈ orbit} χ_k(g_p) · c_p, |orbit|)`: the
/// **summing** projector, paired with the number of *distinct* orbit
/// members. Callers pick their convention —
/// [`canonicalize_pauli_sum_complex`] divides by `|orbit|` to average,
/// [`momentum_merge_pauli_sum_pair`] takes the sum as-is.
///
/// `|orbit|` is `group.order()` only for free orbits; an orbit with a
/// non-trivial stabilizer has fewer distinct members, which is exactly
/// why the two conventions must not be related by a global `|G|` factor.
///
/// Orbits whose stabilizer is incompatible with `k_modes` (the same
/// orbit member reached with different character numerators) project to
/// zero and are omitted from the output.
fn project_onto_reps<A, S, const R: bool>(
    input: &FxHashMap<PauliWord<A, S, R>, Complex<f64>>,
    group: &TranslationGroup,
    k_modes: &[i32],
) -> FxHashMap<PauliWord<A, S, R>, (Complex<f64>, usize)>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
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
        let orbit_size = members.len();
        let mut sum = Complex::new(0.0, 0.0);
        for (member, (counter, _)) in members {
            let coeff = input
                .get(&member)
                .copied()
                .unwrap_or(Complex::new(0.0, 0.0));
            sum += group.character(k_modes, &counter) * coeff;
        }
        projected.insert(rep, (sum, orbit_size));
    }
    projected
}

/// Momentum-sector merge of a complex operator carried as a **real
/// pair**: `re` and `im` are the real and imaginary parts of
/// `O = re + i·im`. Both are overwritten in place with the
/// orbit-representative form of `O` projected onto momentum sector
/// `k_modes`.
///
/// This is the momentum-sector counterpart of
/// [`super::symmetry_merge_pauli_sum`], and generalizes it to `k ≠ 0`
/// while keeping real coefficients on both sums — the only complex
/// arithmetic is the internal character-weighted fold, which reuses
/// [`canonicalize_pauli_sum_complex`].
///
/// This is the **summing** projector
/// `Σ_{p ∈ orbit} χ_k(g_p) · c_p` over each orbit's *distinct* members, not the
/// orbit-averaged one that [`canonicalize_pauli_sum_complex`] returns.
/// Summing is what makes the merge idempotent — and hence safe to apply
/// after every Trotter step — for *every* orbit, including orbits with a
/// non-trivial stabilizer, and it reduces exactly to
/// [`super::symmetry_merge_pauli_sum`] at `k = 0`.
///
/// Entries whose component is exactly zero are dropped, so a purely real
/// operator leaves `im` empty.
///
/// # Panics
///
/// If `re` and `im` disagree on qubit count, if either disagrees with
/// `group.n_qubits()`, or if `k_modes.len() != group.n_generators()`.
pub fn momentum_merge_pauli_sum_pair<T, A, S, const R: bool>(
    re: &mut PauliSum<T>,
    im: &mut PauliSum<T>,
    group: &TranslationGroup,
    k_modes: &[i32],
) where
    T: Config<PauliWordType = PauliWord<A, S, R>, Coeff = f64>,
    T::Map: ACMapAddAssign<T::Storage, f64, T::BuildHasher, PauliWord<A, S, R>>,
    for<'a> T::Map: ACMapIter<'a, Item = (&'a PauliWord<A, S, R>, &'a f64)>,
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    assert_eq!(
        re.n_qubits(),
        im.n_qubits(),
        "real and imaginary parts disagree on qubit count"
    );
    assert_eq!(
        re.n_qubits(),
        group.n_qubits(),
        "PauliSum qubit count {} != group qubit count {}",
        re.n_qubits(),
        group.n_qubits()
    );
    assert_eq!(
        k_modes.len(),
        group.n_generators(),
        "k_modes length {} != number of generators {}",
        k_modes.len(),
        group.n_generators()
    );

    // Gather both real components into `word -> re + i·im`.
    let mut combined: FxHashMap<PauliWord<A, S, R>, Complex<f64>> = FxHashMap::default();
    for (word, coeff) in re.data().iter() {
        combined.entry(*word).or_insert(Complex::new(0.0, 0.0)).re += *coeff;
    }
    for (word, coeff) in im.data().iter() {
        combined.entry(*word).or_insert(Complex::new(0.0, 0.0)).im += *coeff;
    }
    let projected = project_onto_reps(&combined, group, k_modes);
    re.data_mut().clear();
    im.data_mut().clear();
    for (word, (sum, _orbit_size)) in projected {
        if sum.re != 0.0 {
            *re += (word, sum.re);
        }
        if sum.im != 0.0 {
            *im += (word, sum.im);
        }
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
