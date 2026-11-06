// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Retention under piecewise-constant driving, checked against dense evolution.

use nalgebra::{DMatrix, DVector};
use num::Complex;
use ppvm_lindblad::{
    CandidateStrategy, Error, JumpInput, LindbladSpec, PcStepConfig, PcStepWorkspace, Sector, Word,
    word_from_codes,
};
use ppvm_pauli_sum::symmetry::TranslationGroup;

fn models() -> Vec<LindbladSpec> {
    let mut models = Vec::new();
    for (amplitude, rate, general) in [(0.7, 0.2, false), (-1.1, 0.5, true), (0.3, 0.0, false)] {
        let jump = if general {
            vec![
                ("X".into(), Complex::new(0.5, 0.0)),
                ("Y".into(), Complex::new(0.0, -0.5)),
            ]
        } else {
            vec![("Z".into(), Complex::new(1.0, 0.0))]
        };
        models.push(
            LindbladSpec::new(
                1,
                &[
                    ("X".into(), amplitude),
                    ("Y".into(), 0.4 - amplitude),
                    ("Z".into(), -0.2 * amplitude),
                ],
                &[JumpInput {
                    lincomb: jump,
                    rate,
                }],
            )
            .unwrap(),
        );
    }
    let mut kossakowski = LindbladSpec::new(1, &[("Y".into(), -0.6)], &[]).unwrap();
    kossakowski
        .add_kossakowski(
            &[vec![
                ("X".into(), Complex::new(0.5, 0.0)),
                ("Y".into(), Complex::new(0.0, -0.5)),
            ]],
            &[vec![Complex::new(0.3, 0.0)]],
        )
        .unwrap();
    models.push(kossakowski);
    models
}

fn words() -> Vec<Word> {
    (0..4)
        .map(|code| word_from_codes(&[code]).unwrap())
        .collect()
}

fn dense_expected(
    model: &LindbladSpec,
    basis: &[Word],
    coefficients: &[Complex<f64>],
    dt: f64,
    horizon: f64,
) -> Vec<(Word, Complex<f64>)> {
    let full = words();
    let mut generator = DMatrix::<f64>::zeros(4, 4);
    for (row, col, value) in model.generator(&full) {
        generator[(row, col)] += value;
    }
    let before: Vec<_> = full
        .iter()
        .map(|word| {
            basis
                .iter()
                .position(|w| w == word)
                .map_or(Complex::default(), |i| coefficients[i])
        })
        .collect();
    let exponential = (generator * dt).exp();
    let real = &exponential * DVector::from_iterator(4, before.iter().map(|c| c.re));
    let imag = &exponential * DVector::from_iterator(4, before.iter().map(|c| c.im));
    let after: Vec<_> = (0..4).map(|i| Complex::new(real[i], imag[i])).collect();
    let score = |i: usize| (before[i] + (after[i] - before[i]) * (horizon / dt)).norm();
    let mut ranked: Vec<_> = (1..4).collect();
    ranked.sort_by(|&a, &b| score(b).total_cmp(&score(a)));
    full.into_iter()
        .enumerate()
        .filter(|(i, _)| *i == 0 || ranked[..2].contains(i))
        .map(|(i, word)| (word, after[i]))
        .collect()
}

fn assert_state(basis: &[Word], values: &[Complex<f64>], expected: &[(Word, Complex<f64>)]) {
    assert_eq!(basis.len(), expected.len());
    for (word, coefficient) in expected {
        let index = basis
            .iter()
            .position(|w| w == word)
            .expect("predicted retained word");
        assert!((values[index] - coefficient).norm() < 2e-12);
    }
}

#[test]
fn changing_generators_and_timesteps_use_current_step_with_reused_workspace() {
    let models = models();
    for proposer in [CandidateStrategy::Streaming, CandidateStrategy::Shortlist] {
        let mut basis = words();
        let mut coefficients = vec![0.2, 0.8, -0.5, 0.1];
        let mut workspace = PcStepWorkspace::default();
        let cfg = PcStepConfig {
            proposer,
            forward_looking: Some(0.1),
            max_basis: 3,
            admit_basis: Some(4),
            ..Default::default()
        };
        for (step, dt) in [0.025, 0.1, 0.15, 0.05, 0.1, 0.025].into_iter().enumerate() {
            let model = &models[step % models.len()];
            let before: Vec<_> = coefficients.iter().map(|&c| Complex::new(c, 0.0)).collect();
            let expected = dense_expected(model, &basis, &before, dt, 0.1);
            model
                .pc_step_with_workspace(
                    &mut basis,
                    &mut coefficients,
                    dt,
                    &words()[..1],
                    &cfg,
                    &mut workspace,
                )
                .unwrap();
            let actual: Vec<_> = coefficients.iter().map(|&c| Complex::new(c, 0.0)).collect();
            assert_state(&basis, &actual, &expected);
        }
    }
}

#[test]
fn orbit_retention_uses_complex_prediction_and_preserves_evolved_values() {
    let models = models();
    let group = TranslationGroup::chain_1d(1);
    let sector = Sector::new(&group, &[0]);
    let mut basis = words();
    let mut coefficients = vec![
        Complex::new(0.2, 0.1),
        Complex::new(0.8, -0.2),
        Complex::new(-0.3, 0.2),
        Complex::new(0.42, -0.12),
    ];
    let cfg = PcStepConfig {
        forward_looking: Some(0.1),
        max_basis: 3,
        admit_basis: Some(4),
        ..Default::default()
    };
    let mut changed_selection = false;
    for (step, dt) in [0.025, 0.1, 0.15, 0.05].into_iter().enumerate() {
        let model = &models[step % models.len()];
        let expected = dense_expected(model, &basis, &coefficients, dt, 0.1);
        let magnitude = dense_expected(model, &basis, &coefficients, dt, dt);
        changed_selection |= expected.iter().map(|(w, _)| w).collect::<Vec<_>>()
            != magnitude.iter().map(|(w, _)| w).collect::<Vec<_>>();
        model
            .pc_step_orbit_rep(
                &mut basis,
                &mut coefficients,
                dt,
                &words()[..1],
                &sector,
                &cfg,
            )
            .unwrap();
        assert_state(&basis, &coefficients, &expected);
    }
    assert!(
        changed_selection,
        "prediction must exercise a distinct retained support"
    );
    let before = (basis.clone(), coefficients.clone());
    for horizon in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        let invalid = PcStepConfig {
            forward_looking: Some(horizon),
            ..cfg
        };
        assert!(matches!(
            models[0].pc_step_orbit_rep(
                &mut basis,
                &mut coefficients,
                0.025,
                &[],
                &sector,
                &invalid
            ),
            Err(Error::InvalidRankingHorizon)
        ));
        assert_eq!((basis.clone(), coefficients.clone()), before);
    }
}
