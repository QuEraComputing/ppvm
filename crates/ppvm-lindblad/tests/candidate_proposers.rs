// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Public strategy integration and borrowed-tail lifetime checks.

use ppvm_lindblad::{
    CandidateBuffer, CandidateProposer, CandidateStrategy, Error, LindbladSpec, PcStepConfig,
    PcStepWorkspace, StreamingProposer, Word, word_from_codes,
};

#[derive(Default)]
struct DelegatingProposer {
    strategy: StreamingProposer,
    calls: usize,
    selections: Vec<usize>,
}

impl CandidateProposer for DelegatingProposer {
    fn propose(
        &mut self,
        model: &LindbladSpec,
        buffer: &mut CandidateBuffer<'_>,
        protected: &[Word],
        maximum: usize,
    ) -> Result<(), Error> {
        self.calls += 1;
        assert_eq!(buffer.basis().len(), buffer.coefficients().len());
        assert!(buffer.candidates().is_empty());
        self.strategy.propose(model, buffer, protected, maximum)
    }

    fn select(&mut self, buffer: &mut CandidateBuffer<'_>, threshold: f64) -> usize {
        let count = self.strategy.select(buffer, threshold);
        self.selections.push(count);
        count
    }
}

#[test]
fn explicit_strategy_matches_default_through_predictor_and_corrector() {
    let model = LindbladSpec::new(
        4,
        &[
            ("XXII".into(), 0.3),
            ("IYYI".into(), -0.7),
            ("IIXX".into(), 0.2),
        ],
        &[],
    )
    .unwrap();
    let seed = word_from_codes(&[2, 0, 0, 0]).unwrap();
    let (mut a, mut b) = (vec![seed], vec![seed]);
    let (mut ca, mut cb) = (vec![1.0], vec![1.0]);
    let mut workspace = PcStepWorkspace::default();
    let mut proposer = DelegatingProposer::default();
    let cfg = PcStepConfig {
        max_basis: 8,
        admit_basis: Some(24),
        num_threads: Some(2),
        ..Default::default()
    };
    let mut correctors = 0;
    for _ in 0..12 {
        model.pc_step(&mut a, &mut ca, 0.1, &[seed], &cfg).unwrap();
        let first_selection = proposer.selections.len();
        model
            .pc_step_with_proposer(
                &mut b,
                &mut cb,
                0.1,
                &[seed],
                &cfg,
                &mut workspace,
                &mut proposer,
            )
            .unwrap();
        correctors += usize::from(
            proposer
                .selections
                .get(first_selection + 1)
                .is_some_and(|count| *count > 0),
        );
        assert_eq!(a, b);
        assert_eq!(ca, cb);
    }
    assert!(proposer.calls > 12);
    assert!(correctors > 0);
}

struct FixedProposer {
    rejected: Word,
    accepted: Word,
}

impl CandidateProposer for FixedProposer {
    fn propose(
        &mut self,
        _: &LindbladSpec,
        buffer: &mut CandidateBuffer<'_>,
        _: &[Word],
        _: usize,
    ) -> Result<(), Error> {
        buffer.push(self.rejected);
        buffer.push(self.accepted);
        buffer.numeric_scratch().extend([0.0, 1.0]);
        Ok(())
    }

    fn select(&mut self, buffer: &mut CandidateBuffer<'_>, threshold: f64) -> usize {
        if buffer.numeric_scratch()[1].abs() > threshold {
            buffer.candidates_mut().swap(0, 1);
            1
        } else {
            0
        }
    }
}

#[test]
fn independent_strategy_reuses_main_tail_and_admission_initializes_zero_coefficients() {
    let model = LindbladSpec::new(2, &[], &[]).unwrap();
    let seed = word_from_codes(&[2, 0]).unwrap();
    let accepted = word_from_codes(&[0, 1]).unwrap();
    let mut words = Vec::with_capacity(64);
    words.push(seed);
    let address = words.as_ptr();
    let mut coefficients = vec![0.75];
    let mut workspace = PcStepWorkspace::default();
    let mut proposer = FixedProposer {
        rejected: word_from_codes(&[1, 0]).unwrap(),
        accepted,
    };
    let cfg = PcStepConfig {
        max_basis: 2,
        num_threads: Some(1),
        ..Default::default()
    };
    model
        .pc_step_with_proposer(
            &mut words,
            &mut coefficients,
            0.0,
            &[seed],
            &cfg,
            &mut workspace,
            &mut proposer,
        )
        .unwrap();
    assert_eq!(words, vec![seed, accepted]);
    assert_eq!(coefficients, vec![0.75, 0.0]);
    assert_eq!(words.as_ptr(), address);
}

struct FailingProposer;

impl CandidateProposer for FailingProposer {
    fn propose(
        &mut self,
        _: &LindbladSpec,
        buffer: &mut CandidateBuffer<'_>,
        _: &[Word],
        _: usize,
    ) -> Result<(), Error> {
        buffer.push(word_from_codes(&[1, 0]).unwrap());
        buffer.numeric_scratch().push(1.0);
        Err(Error::Internal("candidate generation failed".into()))
    }
    fn select(&mut self, _: &mut CandidateBuffer<'_>, _: f64) -> usize {
        unreachable!()
    }
}

#[test]
fn proposal_error_discards_temporary_tail_and_allows_default_step_to_resume() {
    let model = LindbladSpec::new(2, &[("XX".into(), 0.3)], &[]).unwrap();
    let seed = word_from_codes(&[2, 0]).unwrap();
    let mut words = vec![seed];
    let mut coefficients = vec![0.75];
    let mut workspace = PcStepWorkspace::default();
    let cfg = PcStepConfig {
        max_basis: 8,
        num_threads: Some(1),
        ..Default::default()
    };
    assert!(
        model
            .pc_step_with_proposer(
                &mut words,
                &mut coefficients,
                0.1,
                &[seed],
                &cfg,
                &mut workspace,
                &mut FailingProposer
            )
            .is_err()
    );
    assert_eq!(words, vec![seed]);
    assert_eq!(coefficients, vec![0.75]);
    let (mut expected_words, mut expected_coefficients) = (words.clone(), coefficients.clone());
    model
        .pc_step(
            &mut expected_words,
            &mut expected_coefficients,
            0.1,
            &[seed],
            &cfg,
        )
        .unwrap();
    model
        .pc_step_with_workspace(
            &mut words,
            &mut coefficients,
            0.1,
            &[seed],
            &cfg,
            &mut workspace,
        )
        .unwrap();
    assert_eq!(words, expected_words);
    assert_eq!(coefficients, expected_coefficients);
}

#[test]
fn forward_retention_preserves_values_and_selects_the_predicted_top_terms() {
    let model: LindbladSpec = LindbladSpec::new(
        4,
        &[
            ("XXII".into(), 0.7),
            ("IYYI".into(), 0.8),
            ("IIXX".into(), 0.9),
        ],
        &[],
    )
    .unwrap();
    let seed = word_from_codes(&[2, 0, 0, 0]).unwrap();
    let cfg = PcStepConfig {
        max_basis: 3,
        admit_basis: Some(12),
        ..Default::default()
    };
    let mut basis = vec![seed];
    let mut coeffs = vec![1.0];
    model
        .pc_step(&mut basis, &mut coeffs, 0.1, &[seed], &cfg)
        .unwrap();
    for proposer in [CandidateStrategy::Streaming, CandidateStrategy::Shortlist] {
        for dt in [0.025, 0.1, 0.15] {
            let working = PcStepConfig {
                proposer,
                max_basis: 12,
                ..cfg
            };
            let forward = PcStepConfig {
                proposer,
                forward_looking: Some(0.1),
                ..cfg
            };
            let (mut full, mut evolved) = (basis.clone(), coeffs.clone());
            let (mut selected, mut values) = (basis.clone(), coeffs.clone());
            model
                .pc_step(&mut full, &mut evolved, dt, &[seed], &working)
                .unwrap();
            model
                .pc_step(&mut selected, &mut values, dt, &[seed], &forward)
                .unwrap();
            let score = |i: usize| {
                let before = basis
                    .iter()
                    .position(|w| *w == full[i])
                    .map_or(0.0, |j| coeffs[j]);
                if dt == 0.1 {
                    evolved[i].abs()
                } else {
                    (before + 0.1 / dt * (evolved[i] - before)).abs()
                }
            };
            let mut top: Vec<_> = (0..full.len()).filter(|&i| full[i] != seed).collect();
            top.sort_by(|&i, &j| score(j).total_cmp(&score(i)));
            top.truncate(2);
            let expected: Vec<_> = full
                .iter()
                .enumerate()
                .filter(|(i, w)| **w == seed || top.contains(i))
                .map(|(i, w)| (*w, evolved[i]))
                .collect();
            assert_eq!(
                selected.iter().copied().zip(values).collect::<Vec<_>>(),
                expected
            );
        }
    }
    let invalid = PcStepConfig {
        forward_looking: Some(f64::NAN),
        ..cfg
    };
    assert!(
        model
            .pc_step(&mut basis, &mut coeffs, 0.025, &[seed], &invalid)
            .is_err()
    );
}
