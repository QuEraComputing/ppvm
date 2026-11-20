// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Caller-supplied candidate proposers through the public API.

use ppvm_lindblad::{
    CandidateBuffer, CandidateProposer, Error, LindbladSpec, PcStepConfig, StreamingProposer, Word,
    word_from_codes,
};

/// Wraps the default proposer and records how the step drives it.
#[derive(Default)]
struct Recording {
    inner: StreamingProposer,
    proposals: usize,
    selections: Vec<usize>,
}

impl CandidateProposer for Recording {
    fn propose(
        &mut self,
        model: &LindbladSpec,
        buffer: &mut CandidateBuffer<'_>,
        protected: &[Word],
        maximum: usize,
    ) -> Result<(), Error> {
        self.proposals += 1;
        assert_eq!(buffer.basis().len(), buffer.coefficients().len());
        assert!(buffer.candidates().is_empty());
        self.inner.propose(model, buffer, protected, maximum)
    }

    fn select(&mut self, buffer: &mut CandidateBuffer<'_>, threshold: f64) -> usize {
        let count = self.inner.select(buffer, threshold);
        self.selections.push(count);
        count
    }
}

/// Proposes the single-qubit X strings missing from the basis and admits
/// them in reverse order, using only the public buffer methods.
struct SingleQubitX {
    n: usize,
}

impl CandidateProposer for SingleQubitX {
    fn propose(
        &mut self,
        _model: &LindbladSpec,
        buffer: &mut CandidateBuffer<'_>,
        protected: &[Word],
        maximum: usize,
    ) -> Result<(), Error> {
        let room = maximum.saturating_sub(buffer.basis().len());
        for q in 0..self.n {
            let mut codes = vec![0; self.n];
            codes[q] = 1;
            let word = word_from_codes(&codes).unwrap();
            if buffer.candidates().len() < room
                && !buffer.basis().contains(&word)
                && !protected.contains(&word)
            {
                buffer.push(word);
            }
        }
        Ok(())
    }

    fn select(&mut self, buffer: &mut CandidateBuffer<'_>, _threshold: f64) -> usize {
        buffer.candidates_mut().reverse();
        buffer.candidates().len()
    }
}

fn model() -> LindbladSpec {
    LindbladSpec::new(
        4,
        &[
            ("XXII".into(), 0.3),
            ("IYYI".into(), -0.7),
            ("IIXX".into(), 0.2),
        ],
        &[],
    )
    .unwrap()
}

#[test]
fn wrapped_default_proposer_matches_pc_step() {
    let model = model();
    let seed = word_from_codes(&[2, 0, 0, 0]).unwrap();
    let (mut a, mut b) = (vec![seed], vec![seed]);
    let (mut ca, mut cb) = (vec![1.0], vec![1.0]);
    let mut proposer = Recording::default();
    let cfg = PcStepConfig {
        max_basis: 8,
        admit_basis: Some(24),
        num_threads: Some(2),
        ..Default::default()
    };
    for _ in 0..6 {
        model.pc_step(&mut a, &mut ca, 0.1, &[seed], &cfg).unwrap();
        model
            .pc_step_with_proposer(&mut b, &mut cb, 0.1, &[seed], &cfg, &mut proposer)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(ca, cb);
    }
    assert_eq!(proposer.proposals, proposer.selections.len());
    assert!(proposer.proposals >= 6);
}

#[test]
fn custom_proposer_admits_its_selection_with_zero_coefficients() {
    let model = model();
    let seed = word_from_codes(&[2, 0, 0, 0]).unwrap();
    let mut basis = vec![seed];
    let mut coeffs = vec![1.0];
    let cfg = PcStepConfig {
        max_basis: 3,
        admit_basis: Some(3),
        ..Default::default()
    };
    model
        .pc_step_with_proposer(
            &mut basis,
            &mut coeffs,
            0.0,
            &[seed],
            &cfg,
            &mut SingleQubitX { n: 4 },
        )
        .unwrap();
    let x = |q: usize| {
        let mut codes = vec![0; 4];
        codes[q] = 1;
        word_from_codes(&codes).unwrap()
    };
    assert_eq!(basis, vec![seed, x(1), x(0)]);
    assert_eq!(coeffs, vec![1.0, 0.0, 0.0]);
}
