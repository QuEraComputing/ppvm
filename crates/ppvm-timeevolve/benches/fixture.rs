use ppvm_runtime::{config::fxhash::ByteF64, prelude::*, strategy::CoefficientThreshold};
use ppvm_timeevolve::{JumpOp, LadderDirection, LadderOp, LindbladOp, RateMatrix};

pub type S = ByteF64<1, CoefficientThreshold>;

const N: usize = 5;

/// Build the benchmark Lindblad operator:
/// - 5 lowering operators `S_−_i` at qubits 0..4 (using the fast Ladder kernel)
/// - Dense 5×5 rate matrix γ_ij = 1 / (1 + |i − j|)
pub fn build_lindblad() -> LindbladOp<S> {
    let ops: Vec<JumpOp<S>> = (0..N)
        .map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower }))
        .collect();

    // Dense 5×5 rate matrix: γ_ij = 1 / (1 + |i − j|)
    let rates: Vec<Vec<f64>> = (0..N)
        .map(|i| {
            (0..N)
                .map(|j| 1.0 / (1.0 + (i as f64 - j as f64).abs()))
                .collect()
        })
        .collect();

    LindbladOp::new(ops, RateMatrix::Dense(rates))
}

/// Build the benchmark initial state: P = Σ_i Z_i, threshold 1e-6.
pub fn build_initial() -> PauliSum<S> {
    let strat = CoefficientThreshold(1e-6);
    let mut p: PauliSum<S> = PauliSum::builder().n_qubits(N).strategy(strat).build();
    let template = vec!['I'; N];
    for i in 0..N {
        let mut zi = template.clone();
        zi[i] = 'Z';
        let sz: String = zi.into_iter().collect();
        p += (sz, 1.0_f64);
    }
    p
}
